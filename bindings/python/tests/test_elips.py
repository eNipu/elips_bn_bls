"""
The Python surface of elips.

Two kinds of test here and they answer different questions.

The first kind checks the BINDING: that the right exception type comes out of
the right failure, that bytes go in and out unmangled, and that verification
cannot be mistaken for a boolean. The C is already tested in the repository's
own suite; what is new in this package is the layer above it, and this is
where that layer can be wrong on its own.

The second kind is the cross-check against py_ecc, both directions, skipped if
py_ecc is not installed. It is the same check the C suite runs, driven from
Python, so the binding cannot quietly corrupt a value on its way through.
"""
import pytest

import elips


IKM = bytes(range(32))
MSG = b"the message"


@pytest.fixture(scope="module")
def kp():
    sk = elips.keygen(IKM)
    return sk, elips.sk_to_pk(sk)


# --------------------------------------------------------------- basics ----

def test_sizes_are_what_the_constants_say(kp):
    sk, pk = kp
    assert len(sk) == elips.SK_BYTES == 32
    assert len(pk) == elips.PK_BYTES == 48
    assert len(elips.sign(sk, MSG)) == elips.SIG_BYTES == 96
    assert len(elips.pop_prove(sk)) == elips.POP_BYTES == 96


def test_keygen_is_deterministic_in_its_input():
    assert elips.keygen(IKM) == elips.keygen(IKM)
    assert elips.keygen(IKM) != elips.keygen(bytes(32))


def test_keygen_without_argument_is_random():
    assert elips.keygen() != elips.keygen()


def test_roundtrip(kp):
    sk, pk = kp
    assert elips.verify(pk, MSG, elips.sign(sk, MSG)) is None


def test_verify_returns_none_so_it_cannot_be_used_as_a_bool(kp):
    """The single most important property of this API.

    `if verify(...)` must never be a silent accept. None is falsy, so a caller
    who writes that gets a branch that never fires -- loudly wrong on the first
    valid signature, rather than quietly wrong on the first forged one."""
    sk, pk = kp
    assert elips.verify(pk, MSG, elips.sign(sk, MSG)) is None
    assert not elips.verify(pk, MSG, elips.sign(sk, MSG))


def test_bytearray_and_memoryview_are_accepted(kp):
    sk, pk = kp
    sig = elips.sign(bytearray(sk), memoryview(MSG))
    elips.verify(memoryview(pk), bytearray(MSG), bytearray(sig))


# ------------------------------------------- rejection, and its exact type --

def test_tampered_message_raises_invalid_signature(kp):
    sk, pk = kp
    sig = elips.sign(sk, MSG)
    with pytest.raises(elips.InvalidSignature):
        elips.verify(pk, MSG + b"!", sig)


def test_wrong_key_raises_invalid_signature(kp):
    sk, _ = kp
    other = elips.sk_to_pk(elips.keygen(bytes([9]) * 32))
    with pytest.raises(elips.InvalidSignature):
        elips.verify(other, MSG, elips.sign(sk, MSG))


def test_mangled_signature_raises(kp):
    sk, pk = kp
    sig = bytearray(elips.sign(sk, MSG))
    sig[-1] ^= 1
    with pytest.raises(elips.ElipsError):
        elips.verify(pk, MSG, bytes(sig))


def test_identity_public_key_raises_invalid_key():
    """A well-formed encoding of a useless key. Accepting it would let a
    signer claim membership with no key at all."""
    inf = bytes([0xC0]) + bytes(47)
    with pytest.raises(elips.InvalidKey):
        elips.pk_validate(inf)


def test_every_rejection_is_catchable_as_elips_error(kp):
    """A caller who writes `except ElipsError: reject()` must not be bypassed
    by some failure mode raising something else."""
    sk, pk = kp
    sig = elips.sign(sk, MSG)
    inf = bytes([0xC0]) + bytes(47)
    for args in [(pk, MSG + b"x", sig), (inf, MSG, sig)]:
        with pytest.raises(elips.ElipsError):
            elips.verify(*args)


# ---- wrong length or type is the CALLER's bug, so a different exception ----

@pytest.mark.parametrize("bad", [b"", bytes(31), bytes(33)])
def test_wrong_length_secret_key_is_a_value_error(bad):
    with pytest.raises(ValueError):
        elips.sk_to_pk(bad)


def test_wrong_length_public_key_is_a_value_error(kp):
    sk, pk = kp
    with pytest.raises(ValueError):
        elips.verify(pk[:-1], MSG, elips.sign(sk, MSG))


def test_wrong_type_is_a_type_error():
    with pytest.raises(TypeError):
        elips.sk_to_pk("not bytes")
    with pytest.raises(TypeError):
        elips.sign(bytes(32), "not bytes")


def test_short_ikm_is_a_value_error():
    with pytest.raises(ValueError):
        elips.keygen(bytes(31))


def test_value_error_is_not_caught_by_an_elips_error_handler():
    """Deliberate: a malformed length is a bug in the calling code, not a bad
    signature, and conflating the two hides the bug."""
    assert not issubclass(ValueError, elips.ElipsError)
    assert not issubclass(elips.ElipsError, ValueError)


# ------------------------------------------------------------ aggregation --

@pytest.fixture(scope="module")
def group():
    sks = [elips.keygen(bytes([i]) * 32) for i in range(1, 5)]
    pks = [elips.sk_to_pk(s) for s in sks]
    pops = [elips.pop_prove(s) for s in sks]
    return sks, pks, pops


def test_aggregate_of_n_is_still_one_signature(group):
    sks, pks, _ = group
    msgs = [b"m%d" % i for i in range(len(sks))]
    sigs = [elips.sign(s, m) for s, m in zip(sks, msgs)]
    agg = elips.aggregate(sigs)
    assert len(agg) == elips.SIG_BYTES
    assert elips.aggregate_verify(pks, msgs, agg) is None


def test_aggregate_verify_rejects_reordered_messages(group):
    sks, pks, _ = group
    msgs = [b"m%d" % i for i in range(len(sks))]
    agg = elips.aggregate([elips.sign(s, m) for s, m in zip(sks, msgs)])
    swapped = [msgs[1], msgs[0]] + msgs[2:]
    with pytest.raises(elips.InvalidSignature):
        elips.aggregate_verify(pks, swapped, agg)


def test_aggregate_verify_refuses_a_repeated_message(group):
    """Without proofs of possession a repeat is forgeable, so it is refused
    rather than documented. The message names the safe alternative."""
    sks, pks, _ = group
    msgs = [b"m%d" % i for i in range(len(sks))]
    agg = elips.aggregate([elips.sign(s, m) for s, m in zip(sks, msgs)])
    with pytest.raises(ValueError, match="distinct"):
        elips.aggregate_verify(pks, [msgs[0]] * len(pks), agg)


def test_empty_inputs_are_refused(group):
    _, pks, _ = group
    with pytest.raises(ValueError):
        elips.aggregate([])
    with pytest.raises(ValueError):
        elips.aggregate_verify([], [], bytes(96))


def test_mismatched_counts_are_refused(group):
    sks, pks, pops = group
    msgs = [b"m%d" % i for i in range(len(sks))]
    agg = elips.aggregate([elips.sign(s, m) for s, m in zip(sks, msgs)])
    with pytest.raises(ValueError):
        elips.aggregate_verify(pks, msgs[:-1], agg)
    with pytest.raises(ValueError):
        elips.fast_aggregate_verify(pks, pops[:-1], b"x", agg)


# -------------------------------------------------- shared message and PoP --

def test_shared_message_needs_the_proofs_and_accepts_them(group):
    sks, pks, pops = group
    shared = b"one message, many signers"
    agg = elips.aggregate([elips.sign(s, shared) for s in sks])
    assert elips.fast_aggregate_verify(pks, pops, shared, agg) is None


def test_one_bad_proof_stops_the_whole_aggregate(group):
    sks, pks, pops = group
    shared = b"one message, many signers"
    agg = elips.aggregate([elips.sign(s, shared) for s in sks])
    bad = list(pops)
    b = bytearray(bad[2]); b[-1] ^= 1; bad[2] = bytes(b)
    with pytest.raises(elips.ElipsError):
        elips.fast_aggregate_verify(pks, bad, shared, agg)


def test_a_proof_does_not_verify_under_another_key(group):
    _, pks, pops = group
    assert elips.pop_verify(pks[0], pops[0]) is None
    with pytest.raises(elips.ElipsError):
        elips.pop_verify(pks[1], pops[0])


# ----------------------------------------- the same reference the C uses ----

py_ecc = pytest.importorskip("py_ecc.bls", reason="py_ecc is not installed")


def test_agrees_with_py_ecc_both_directions():
    """A signature we produced that py_ecc accepts shows our signing follows
    the draft; one py_ecc produced that we accept shows our verification does.
    An implementation can be self-consistently wrong and pass either alone."""
    from py_ecc.bls import G2ProofOfPossession as BLS

    for i in range(4):
        ikm = bytes([i + 1]) * 32
        msg = b"cross-check %d" % i

        sk = elips.keygen(ikm)
        assert sk == (b"%064x" % BLS.KeyGen(ikm)).decode().encode() or \
            int.from_bytes(sk, "big") == BLS.KeyGen(ikm)

        pk = elips.sk_to_pk(sk)
        assert pk == BLS.SkToPk(int.from_bytes(sk, "big"))

        ours = elips.sign(sk, msg)
        theirs = BLS.Sign(int.from_bytes(sk, "big"), msg)
        assert ours == theirs

        assert BLS.Verify(pk, msg, ours)          # they accept ours
        elips.verify(pk, msg, theirs)             # we accept theirs


def test_agrees_with_py_ecc_on_aggregates():
    from py_ecc.bls import G2ProofOfPossession as BLS

    sks = [elips.keygen(bytes([i]) * 32) for i in range(5, 9)]
    pks = [elips.sk_to_pk(s) for s in sks]
    msgs = [b"agg %d" % i for i in range(len(sks))]
    sigs = [elips.sign(s, m) for s, m in zip(sks, msgs)]

    ours = elips.aggregate(sigs)
    assert ours == BLS.Aggregate(sigs)
    assert BLS.AggregateVerify(pks, msgs, ours)
    elips.aggregate_verify(pks, msgs, ours)
