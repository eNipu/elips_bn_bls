/**
 * BLS signatures on BLS12-381, in the browser and in node.
 *
 *   import { load } from "@elips/bls";
 *   const elips = await load();
 *
 *   const sk  = elips.keygen();
 *   const pk  = elips.skToPk(sk);
 *   const sig = elips.sign(sk, new TextEncoder().encode("hello"));
 *   elips.verify(pk, msg, sig);        // returns undefined, or throws
 *
 * Everything is Uint8Array in and Uint8Array out. No curve object, no field
 * element, no domain separation tag.
 *
 * VERIFICATION THROWS, IT DOES NOT RETURN A BOOLEAN. A verify that returns a
 * value invites `if (verify(...)) accept()`, which is wrong in the direction
 * that matters. Returning undefined makes that branch never taken, so the
 * mistake shows up on the first VALID signature rather than the first forged
 * one.
 *
 *   try { elips.verify(pk, msg, sig); }
 *   catch (e) { reject(); }
 *
 * Catch ElipsError to reject. A TypeError or RangeError means the calling code
 * passed the wrong type or length, which is a bug in the caller rather than a
 * bad signature, so it is a different error and is thrown before any
 * cryptography runs.
 *
 * Scheme: draft-irtf-cfrg-bls-signature, proof of possession,
 * minimal-pubkey-size. Public keys 48 bytes, signatures and proofs 96, secret
 * keys 32. The variant Ethereum uses.
 *
 * SPEED: this is materially slower than the native library, and the reason is
 * structural rather than fixable here. wasm32 has no 64x64 -> 128 multiply, so
 * every field multiplication goes through a software helper. Fine for a demo
 * or for verifying the occasional signature; not a basis for a throughput
 * claim. The measured numbers are in the README.
 *
 * NOT AUDITED.
 */
import initWasm from "./dist/elips_wasm.mjs";

export const SK_BYTES = 32;
export const PK_BYTES = 48;
export const SIG_BYTES = 96;
export const POP_BYTES = 96;

/** Base for anything this module raises about the data it was given.
 *  Catch this to reject a signature: it covers one that does not verify, one
 *  that does not decode, and an unusable public key. To a verifier all three
 *  mean the same thing. */
export class ElipsError extends Error {
  constructor(msg) { super(msg); this.name = "ElipsError"; }
}
export class InvalidSignature extends ElipsError {
  constructor(msg) { super(msg); this.name = "InvalidSignature"; }
}
export class InvalidKey extends ElipsError {
  constructor(msg) { super(msg); this.name = "InvalidKey"; }
}

const ERR_INVALID = -1, ERR_BAD_KEY = -2, ERR_BAD_SIG = -3;
const ERR_VERIFY = -4, ERR_RANDOM = -5, ERR_DUP_MESSAGE = -6;

let modulePromise = null;

/**
 * Load the WebAssembly module and return the API. Await it once; the result is
 * cached, so calling load() again is free.
 */
export function load() {
  if (!modulePromise) modulePromise = initWasm().then(build);
  return modulePromise;
}

function build(m) {
  const heap = () => m.HEAPU8;

  function fail(code) {
    const msg = m.UTF8ToString(m._elips_bls_strerror(code));
    if (code === ERR_VERIFY || code === ERR_BAD_SIG) throw new InvalidSignature(msg);
    if (code === ERR_BAD_KEY) throw new InvalidKey(msg);
    if (code === ERR_RANDOM) throw new ElipsError("entropy source failed");
    throw new ElipsError(msg);
  }

  function bytes(name, v) {
    if (!(v instanceof Uint8Array)) {
      throw new TypeError(`${name} must be a Uint8Array`);
    }
    return v;
  }

  function sized(name, v, n) {
    bytes(name, v);
    if (v.length !== n) {
      throw new RangeError(`${name} must be exactly ${n} bytes, got ${v.length}`);
    }
    return v;
  }

  /* One arena per call, freed in a finally, so a thrown error cannot leak
   * WebAssembly memory. */
  function withMemory(sizes, fn) {
    const ptrs = sizes.map((s) => m._malloc(Math.max(s, 1)));
    try {
      return fn(ptrs);
    } finally {
      for (const p of ptrs) m._free(p);
    }
  }

  const put = (ptr, src) => heap().set(src, ptr);
  const take = (ptr, n) => heap().slice(ptr, ptr + n);

  function concat(name, items, size) {
    if (!Array.isArray(items) || items.length === 0) {
      throw new RangeError(`${name} must be a non-empty array`);
    }
    const out = new Uint8Array(items.length * size);
    items.forEach((x, i) => out.set(sized(`${name}[${i}]`, x, size), i * size));
    return out;
  }

  return {
    SK_BYTES, PK_BYTES, SIG_BYTES, POP_BYTES,
    ElipsError, InvalidSignature, InvalidKey,

    /** A 32-byte secret key. With no argument the browser's CSPRNG supplies
     *  the input; pass at least 32 bytes to derive deterministically, in which
     *  case that input is as secret as the key. */
    keygen(ikm) {
      return withMemory([SK_BYTES, ikm ? ikm.length : 0], ([sk, ip]) => {
        let rc;
        if (ikm === undefined) {
          rc = m._elips_bls_keygen_random(sk);
        } else {
          bytes("ikm", ikm);
          if (ikm.length < 32) throw new RangeError("ikm must be at least 32 bytes");
          put(ip, ikm);
          rc = m._elips_bls_keygen(sk, ip, ikm.length);
        }
        if (rc !== 0) fail(rc);
        return take(sk, SK_BYTES);
      });
    },

    skToPk(sk) {
      sized("sk", sk, SK_BYTES);
      return withMemory([PK_BYTES, SK_BYTES], ([pk, s]) => {
        put(s, sk);
        const rc = m._elips_bls_sk_to_pk(pk, s);
        if (rc !== 0) fail(rc);
        return take(pk, PK_BYTES);
      });
    },

    /** Throws unless pk is usable: decodes, is not the identity, is in G1. */
    pkValidate(pk) {
      sized("pk", pk, PK_BYTES);
      withMemory([PK_BYTES], ([p]) => {
        put(p, pk);
        const rc = m._elips_bls_pk_validate(p);
        if (rc !== 0) fail(rc);
      });
    },

    sign(sk, msg) {
      sized("sk", sk, SK_BYTES);
      bytes("msg", msg);
      return withMemory([SIG_BYTES, SK_BYTES, msg.length], ([sig, s, mp]) => {
        put(s, sk); put(mp, msg);
        const rc = m._elips_bls_sign(sig, s, mp, msg.length);
        if (rc !== 0) fail(rc);
        return take(sig, SIG_BYTES);
      });
    },

    /** Returns undefined if the signature is good, throws if not. There is no
     *  value to test, on purpose. */
    verify(pk, msg, sig) {
      sized("pk", pk, PK_BYTES);
      sized("sig", sig, SIG_BYTES);
      bytes("msg", msg);
      withMemory([PK_BYTES, msg.length, SIG_BYTES], ([p, mp, sp]) => {
        put(p, pk); put(mp, msg); put(sp, sig);
        const rc = m._elips_bls_verify(p, mp, msg.length, sp);
        if (rc !== 0) fail(rc);
      });
    },

    /** Combine signatures into one, still 96 bytes. That is the point of BLS. */
    aggregate(sigs) {
      const blob = concat("sigs", sigs, SIG_BYTES);
      return withMemory([SIG_BYTES, blob.length], ([out, bp]) => {
        put(bp, blob);
        const rc = m._elips_bls_aggregate(out, bp, sigs.length);
        if (rc !== 0) fail(rc);
        return take(out, SIG_BYTES);
      });
    },

    /** One aggregate against n keys and n DISTINCT messages. A repeated
     *  message is refused: without proofs of possession an attacker can move a
     *  signature between signers. Use fastAggregateVerify for one shared
     *  message. */
    aggregateVerify(pks, msgs, sig) {
      const pkBlob = concat("pks", pks, PK_BYTES);
      if (!Array.isArray(msgs) || msgs.length !== pks.length) {
        throw new RangeError(`got ${pks.length} keys and ${msgs?.length} messages`);
      }
      msgs.forEach((x, i) => bytes(`msgs[${i}]`, x));
      sized("sig", sig, SIG_BYTES);

      const n = pks.length;
      const sizes = [pkBlob.length, SIG_BYTES, n * 4, n * 4,
                     ...msgs.map((x) => x.length)];
      return withMemory(sizes, (ptrs) => {
        const [pb, sp, arr, lens, ...body] = ptrs;
        put(pb, pkBlob); put(sp, sig);
        const ptr32 = new Uint32Array(n), len32 = new Uint32Array(n);
        msgs.forEach((x, i) => { put(body[i], x); ptr32[i] = body[i]; len32[i] = x.length; });
        heap().set(new Uint8Array(ptr32.buffer), arr);
        heap().set(new Uint8Array(len32.buffer), lens);
        const rc = m._elips_bls_aggregate_verify(pb, n, arr, lens, sp);
        if (rc === ERR_DUP_MESSAGE) {
          throw new RangeError("aggregateVerify needs distinct messages; use " +
                               "fastAggregateVerify for one shared message");
        }
        if (rc !== 0) fail(rc);
      });
    },

    /** One aggregate against n keys over ONE shared message.
     *
     *  It takes the proofs of possession and checks them first, and there is
     *  no way to call it without them. Summing public keys over a shared
     *  message is forgeable otherwise: an attacker who registers
     *  pk_evil = [t]G1 - sum(pk_honest) can produce an aggregate for a group
     *  whose keys they never held. */
    fastAggregateVerify(pks, pops, msg, sig) {
      const pkBlob = concat("pks", pks, PK_BYTES);
      const popBlob = concat("pops", pops, POP_BYTES);
      if (pks.length !== pops.length) {
        throw new RangeError(`got ${pks.length} keys and ${pops.length} proofs`);
      }
      bytes("msg", msg);
      sized("sig", sig, SIG_BYTES);
      withMemory([pkBlob.length, popBlob.length, msg.length, SIG_BYTES],
        ([pb, qb, mp, sp]) => {
          put(pb, pkBlob); put(qb, popBlob); put(mp, msg); put(sp, sig);
          const rc = m._elips_bls_fast_aggregate_verify(
            pb, qb, pks.length, mp, msg.length, sp);
          if (rc !== 0) fail(rc);
        });
    },

    popProve(sk) {
      sized("sk", sk, SK_BYTES);
      return withMemory([POP_BYTES, SK_BYTES], ([out, s]) => {
        put(s, sk);
        const rc = m._elips_bls_pop_prove(out, s);
        if (rc !== 0) fail(rc);
        return take(out, POP_BYTES);
      });
    },

    popVerify(pk, proof) {
      sized("pk", pk, PK_BYTES);
      sized("proof", proof, POP_BYTES);
      withMemory([PK_BYTES, POP_BYTES], ([p, q]) => {
        put(p, pk); put(q, proof);
        const rc = m._elips_bls_pop_verify(p, q);
        if (rc !== 0) fail(rc);
      });
    },
  };
}

export default load;
