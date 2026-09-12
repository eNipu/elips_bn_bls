"""
BLS signatures on BLS12-381, from the ELiPS pairing library.

    >>> import elips
    >>> sk  = elips.keygen()
    >>> pk  = elips.sk_to_pk(sk)
    >>> sig = elips.sign(sk, b"hello")
    >>> elips.verify(pk, b"hello", sig)          # returns None, or raises

Everything is ``bytes`` in and ``bytes`` out. There is no curve object, no
field element and no domain separation tag to pass, because getting any of
those wrong is how signature code goes quietly wrong.

VERIFICATION RAISES, IT DOES NOT RETURN A BOOLEAN. This is the single most
important thing on this page. A verify that returns a value can be used as::

    if elips.verify(pk, msg, sig):      # WRONG, and silent
        accept()

and a function returning ``None`` makes that branch never taken, while a
function returning ``False`` makes it look like it works until the day it
matters. So verification returns ``None`` on success and raises on failure::

    try:
        elips.verify(pk, msg, sig)
    except elips.ElipsError:
        reject()

Catch ``ElipsError`` to reject. A ``ValueError`` means the calling code passed
something the wrong length or the wrong type, which is a bug in the caller
rather than a bad signature, so it is deliberately a different exception and
is raised before any cryptography happens.

Scheme: draft-irtf-cfrg-bls-signature, proof of possession, minimal-pubkey-size.
Public keys are 48 bytes, signatures and proofs 96, secret keys 32. This is
what Ethereum uses, and every operation is checked against py_ecc on every
build of the underlying library.

NOT AUDITED. It implements the draft and agrees with an independent
implementation; that is a different claim from having been reviewed.
"""

from __future__ import annotations

from typing import Sequence

from ._elips import ffi as _ffi, lib as _lib

__all__ = [
    "SK_BYTES", "PK_BYTES", "SIG_BYTES", "POP_BYTES",
    "ElipsError", "InvalidSignature", "InvalidKey",
    "keygen", "sk_to_pk", "pk_validate",
    "sign", "verify",
    "aggregate", "aggregate_verify", "fast_aggregate_verify",
    "pop_prove", "pop_verify",
]

SK_BYTES = 32
PK_BYTES = 48
SIG_BYTES = 96
POP_BYTES = 96

# Mirrors the codes in include/elips/bls.h.
_OK = 0
_ERR_INVALID = -1
_ERR_BAD_KEY = -2
_ERR_BAD_SIG = -3
_ERR_VERIFY = -4
_ERR_RANDOM = -5
_ERR_DUP_MESSAGE = -6


class ElipsError(Exception):
    """Base for everything this package raises about the data it was given.

    Catch this to reject a signature. It covers a signature that does not
    verify, one that does not decode, and an unusable public key, because from
    a verifier's point of view all three mean the same thing: do not accept."""


class InvalidSignature(ElipsError):
    """The signature does not verify, or does not decode as a G2 point."""


class InvalidKey(ElipsError):
    """The public key is the identity, or is not in the right subgroup."""


def _raise(code: int) -> None:
    msg = _ffi.string(_lib.elips_bls_strerror(code)).decode()
    if code in (_ERR_VERIFY, _ERR_BAD_SIG):
        raise InvalidSignature(msg)
    if code == _ERR_BAD_KEY:
        raise InvalidKey(msg)
    if code == _ERR_RANDOM:
        raise OSError(msg)
    raise ElipsError(msg)


def _check(name: str, value: object, size: int) -> bytes:
    """Length and type are the caller's responsibility, so a mistake there is a
    ValueError raised before any cryptography runs, not a failed signature."""
    if not isinstance(value, (bytes, bytearray, memoryview)):
        raise TypeError(f"{name} must be bytes, not {type(value).__name__}")
    b = bytes(value)
    if len(b) != size:
        raise ValueError(f"{name} must be exactly {size} bytes, got {len(b)}")
    return b


def _msg(value: object) -> bytes:
    if not isinstance(value, (bytes, bytearray, memoryview)):
        raise TypeError("message must be bytes, not " + type(value).__name__)
    return bytes(value)


# --------------------------------------------------------------------- keys

def keygen(ikm: bytes | None = None) -> bytes:
    """A 32-byte secret key.

    With no argument, the system entropy source supplies the input. Pass your
    own ``ikm`` of at least 32 bytes to derive deterministically: the same
    input always gives the same key, which is what makes seed phrases and test
    vectors work. That also means ``ikm`` is as secret as the key itself."""
    out = _ffi.new("uint8_t[]", SK_BYTES)
    if ikm is None:
        rc = _lib.elips_bls_keygen_random(out)
    else:
        b = _msg(ikm)
        if len(b) < 32:
            raise ValueError(f"ikm must be at least 32 bytes, got {len(b)}")
        rc = _lib.elips_bls_keygen(out, b, len(b))
    if rc != _OK:
        _raise(rc)
    return bytes(_ffi.buffer(out, SK_BYTES))


def sk_to_pk(sk: bytes) -> bytes:
    """The 48-byte public key for a secret key."""
    b = _check("sk", sk, SK_BYTES)
    out = _ffi.new("uint8_t[]", PK_BYTES)
    rc = _lib.elips_bls_sk_to_pk(out, b)
    if rc != _OK:
        _raise(rc)
    return bytes(_ffi.buffer(out, PK_BYTES))


def pk_validate(pk: bytes) -> None:
    """Raise unless ``pk`` is a usable public key.

    Decodes it, rejects the identity, and checks subgroup membership. Worth
    calling once on anything that arrives from the network before you store
    it; ``verify`` does it for you every time."""
    b = _check("pk", pk, PK_BYTES)
    rc = _lib.elips_bls_pk_validate(b)
    if rc != _OK:
        _raise(rc)


# ---------------------------------------------------------- sign and verify

def sign(sk: bytes, msg: bytes) -> bytes:
    """A 96-byte signature over ``msg``."""
    s = _check("sk", sk, SK_BYTES)
    m = _msg(msg)
    out = _ffi.new("uint8_t[]", SIG_BYTES)
    rc = _lib.elips_bls_sign(out, s, m, len(m))
    if rc != _OK:
        _raise(rc)
    return bytes(_ffi.buffer(out, SIG_BYTES))


def verify(pk: bytes, msg: bytes, sig: bytes) -> None:
    """Return None if the signature is good. Raise ``ElipsError`` if not.

    There is no return value to test, on purpose. See the module docstring."""
    p = _check("pk", pk, PK_BYTES)
    s = _check("sig", sig, SIG_BYTES)
    m = _msg(msg)
    rc = _lib.elips_bls_verify(p, m, len(m), s)
    if rc != _OK:
        _raise(rc)


# ------------------------------------------------------------- aggregation

def _concat(name: str, items: Sequence[bytes], size: int) -> tuple[bytes, int]:
    n = len(items)
    if n == 0:
        raise ValueError(f"{name} must not be empty")
    return b"".join(_check(f"{name}[{i}]", x, size) for i, x in enumerate(items)), n


def aggregate(sigs: Sequence[bytes]) -> bytes:
    """Combine signatures into one, still 96 bytes. That is the point of BLS."""
    blob, n = _concat("sigs", sigs, SIG_BYTES)
    out = _ffi.new("uint8_t[]", SIG_BYTES)
    rc = _lib.elips_bls_aggregate(out, blob, n)
    if rc != _OK:
        _raise(rc)
    return bytes(_ffi.buffer(out, SIG_BYTES))


def aggregate_verify(pks: Sequence[bytes], msgs: Sequence[bytes],
                     sig: bytes) -> None:
    """Verify one aggregate against n keys and n DISTINCT messages.

    A repeated message is rejected. Without proofs of possession an attacker
    can otherwise move a signature between signers, so the check is refused
    rather than documented. For one shared message, use
    ``fast_aggregate_verify``, which takes the proofs."""
    pk_blob, n = _concat("pks", pks, PK_BYTES)
    if len(msgs) != n:
        raise ValueError(f"got {n} keys and {len(msgs)} messages")
    s = _check("sig", sig, SIG_BYTES)

    bodies = [_msg(m) for m in msgs]
    keep = [_ffi.new("uint8_t[]", b) for b in bodies]       # keep alive
    arr = _ffi.new("uint8_t*[]", keep)
    lens = _ffi.new("size_t[]", [len(b) for b in bodies])
    rc = _lib.elips_bls_aggregate_verify(pk_blob, n, arr, lens, s)
    if rc == _ERR_DUP_MESSAGE:
        raise ValueError("aggregate_verify needs distinct messages; use "
                         "fast_aggregate_verify for one shared message")
    if rc != _OK:
        _raise(rc)


def fast_aggregate_verify(pks: Sequence[bytes], pops: Sequence[bytes],
                          msg: bytes, sig: bytes) -> None:
    """Verify one aggregate against n keys over ONE shared message.

    It takes the proofs of possession and checks them first, and there is no
    way to call it without them. Summing public keys over a shared message is
    forgeable otherwise: an attacker who registers
    ``pk_evil = [t]G1 - sum(pk_honest)`` can produce an aggregate for a group
    whose keys they never held."""
    pk_blob, n = _concat("pks", pks, PK_BYTES)
    pop_blob, m = _concat("pops", pops, POP_BYTES)
    if n != m:
        raise ValueError(f"got {n} keys and {m} proofs")
    s = _check("sig", sig, SIG_BYTES)
    body = _msg(msg)
    rc = _lib.elips_bls_fast_aggregate_verify(pk_blob, pop_blob, n,
                                              body, len(body), s)
    if rc != _OK:
        _raise(rc)


# -------------------------------------------------- proof of possession

def pop_prove(sk: bytes) -> bytes:
    """A 96-byte proof that you hold the secret key for your public key."""
    b = _check("sk", sk, SK_BYTES)
    out = _ffi.new("uint8_t[]", POP_BYTES)
    rc = _lib.elips_bls_pop_prove(out, b)
    if rc != _OK:
        _raise(rc)
    return bytes(_ffi.buffer(out, POP_BYTES))


def pop_verify(pk: bytes, proof: bytes) -> None:
    """Return None if the proof is good for this key. Raise if not."""
    p = _check("pk", pk, PK_BYTES)
    q = _check("proof", proof, POP_BYTES)
    rc = _lib.elips_bls_pop_verify(p, q)
    if rc != _OK:
        _raise(rc)
