"""
05_phi_closure.py  --  CM-endomorphism (phi) closure properties of secp256k1.

All statements below are numerically verified against exact group arithmetic.

FINDINGS
 (F1) phi(P) = (zeta*x, y) is an Fp-endomorphism of E: y^2 = x^3 + 7 with
     phi(G) = lambda*G, lambda public and fixed by lambda^2+lambda+1 = 0 mod n.
 (F2) NEGATIVE RESULT (refutes a plausible attack): a single valid BIP-340
     signature cannot be malleated into a second valid signature of the same
     message via the CM endomorphism.  Reason: the response is linear in the
     secret, s' = lambda*s + (e' - lambda*e)*d, and the d-term does not cancel.
     (The classical (R,s) -> (-R,-s) malleability is separately killed by
     BIP-340's even-y rule.)  =>  BIP-340 signature uniqueness survives CM.
 (F3) POSITIVE RESULT: two signatures whose nonces are CM-related
     (k2 = lambda * k1) yield the private key in closed form:
        Schnorr: d = (s2 - lambda*s1) / (e2 - lambda*e1)
        ECDSA:   d = (lambda*s2*e1 - s1*e2) / (s1*r2 - lambda*s2*r1)
     and the relation is INVISIBLE to the standard "duplicate R / duplicate r"
     detector, because r2 = zeta * r1 (x-only closure) rather than r2 = r1.
 (F4) dual-generator hazard: H = phi(G) has PUBLICLY KNOWN discrete log
     relative to G: dlog_H(G) = lambda^{-1} mod n.  Therefore H = phi(G) is a
     valid substitute for a "second independent generator" in any scheme that
     needs an unknown DL relation (dual-base Schnorr, two-base multisig,
     accountable ring signatures, NUMS-generator setups with a proof of
     knowledge of discrete log -- that PoK can be forged because the log is
     public).
 (F5) x-only transcript space is 3-fold closed: {r, zeta*r, zeta^2*r}.
"""
import hashlib
import json
import os
import random

from sclab import (A, B, N, P, add_aff, cube_roots_of_unity,
                   lambda_from_zeta, mul, on_curve, phi_point)

G = (0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798,
     0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8)

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(os.path.dirname(HERE), "results")


def tagged_hash(tag, msg):
    t = hashlib.sha256(tag.encode()).digest()
    h = hashlib.sha256()
    h.update(t + t + msg)
    return h.digest()


def inv_mod(a, m):
    return pow(a, -1, m)


def lift_x(x):
    if x >= P:
        return None
    y2 = (pow(x, 3, P) + B) % P
    y = pow(y2, (P + 1) // 4, P)
    if pow(y, 2, P) != y2:
        return None
    return (x, y if y % 2 == 0 else P - y)


def bip340_sign(msg, d, aux=b"\x00" * 32):
    Pt = mul(d, G)
    Px = Pt[0]
    t = d.to_bytes(32, "big")
    rand = tagged_hash("BIP0340/aux", aux)
    k0 = int.from_bytes(tagged_hash("BIP0340/nonce",
                                    t + bytes([0]) + Px.to_bytes(32, "big") + msg + rand), "big") % N
    R = mul(k0, G)
    k = k0 if R[1] % 2 == 0 else (N - k0) % N
    rx = R[0]
    e = int.from_bytes(tagged_hash("BIP0340/challenge",
                                    rx.to_bytes(32, "big") + Px.to_bytes(32, "big") + msg), "big") % N
    s = (k + e * d) % N
    return rx, s, k0, e, Px


def bip340_verify(msg, rx, s, Px):
    if s >= N:
        return False, "s out of range"
    e = int.from_bytes(tagged_hash("BIP0340/challenge",
                                    rx.to_bytes(32, "big") + Px.to_bytes(32, "big") + msg), "big") % N
    R = add_aff(mul(s, G), mul((N - e) % N, lift_x(Px)))
    if R is None:
        return False, "R = O"
    if R[1] % 2 != 0:
        return False, "R has odd y"
    if R[0] != rx:
        return False, "x mismatch"
    return True, "ok"


def ecdsa_sign_raw(msg32, d, k):
    R = mul(k, G)
    r = R[0] % N
    e = int.from_bytes(hashlib.sha256(msg32).digest(), "big") % N
    s = inv_mod(k, N) * (e + r * d) % N
    return r, s, e


def main():
    random.seed(20261022)
    out = {}
    Z1, Z2 = cube_roots_of_unity(P)
    LAM = lambda_from_zeta(Z1)
    out["zeta"] = Z1
    out["lambda"] = LAM
    assert (LAM * LAM + LAM + 1) % N == 0
    print("zeta  =", Z1, flush=True)
    print("lambda=", LAM, flush=True)

    H = phi_point(G, Z1)
    assert on_curve(H, 0, B, P) and mul(LAM, G) == H
    out["H_phi_G"] = H
    print("F1 ok: phi(G) = lambda*G on curve", flush=True)

    d = 0x1234567890ABCDEF1234567890ABCDEF1234567890ABCDEF1234567890ABCDEF
    msg = b"attack-at-dawn"
    rx, s, k0, e, Px = bip340_sign(msg, d)
    ok, why = bip340_verify(msg, rx, s, Px)
    assert ok
    neg = []
    for zz, tag, lm in ((Z1, "zeta", LAM), (Z2, "zeta^2", LAM * LAM % N)):
        r2 = zz * rx % P
        e2 = int.from_bytes(tagged_hash("BIP0340/challenge",
                                         r2.to_bytes(32, "big") + Px.to_bytes(32, "big") + msg), "big") % N
        s_guess = (lm * s + e2 - lm * e) % N
        s_true = (lm * s + (e2 - lm * e) * d) % N
        ok_g, _ = bip340_verify(msg, r2, s_guess, Px)
        ok_t, _ = bip340_verify(msg, r2, s_true, Px)
        neg.append({"variant": tag, "multiplier": lm, "naive_guess_valid": ok_g,
                    "true_form_valid": ok_t,
                    "d_term_coefficient": (e2 - lm * e) % N})
        print(f"   F2[{tag}] naive d-free guess valid={ok_g} | true d-bearing form valid={ok_t}",
              flush=True)
    assert all((not n["naive_guess_valid"]) and n["true_form_valid"] for n in neg)
    out["f2_single_sig_malleability_blocked"] = neg
    print("F2 confirmed: single-signature CM malleability requires d -> blocked", flush=True)

    d1 = 0x0FEDCBA9876543210FEDCBA9876543210FEDCBA9876543210FEDCBA987654321
    m1, m2 = b"msg-A", b"msg-B"
    k1 = random.randrange(1, N)
    k2 = (LAM * k1) % N
    R1, R2 = mul(k1, G), mul(k2, G)
    r1, r2 = R1[0], R2[0]
    assert r2 == Z1 * r1 % P
    e1 = int.from_bytes(hashlib.sha256(b"chal" + r1.to_bytes(32, "big") + m1).digest(), "big") % N
    e2 = int.from_bytes(hashlib.sha256(b"chal" + r2.to_bytes(32, "big") + m2).digest(), "big") % N
    s1 = (k1 + e1 * d1) % N
    s2 = (k2 + e2 * d1) % N
    d_rec = ((s2 - LAM * s1) * inv_mod((e2 - LAM * e1) % N, N)) % N
    assert d_rec == d1 and mul(d_rec, G) == mul(d1, G)
    out["f3_schnorr_cm_nonce"] = {"r1": r1, "r2": r2, "r2_eq_zeta_r1": True,
                                  "d_recovered": d_rec, "correct": True,
                                  "naive_duplicate_r_detector_fires": r1 == r2}
    print("F3a schnorr: d recovered ==", d_rec == d1, "| naive dup-r detector fires:",
          r1 == r2, flush=True)

    d2 = 0x1111111111111111222222222222222233333333333333334444444444444444
    k1 = random.randrange(1, N)
    k2 = (LAM * k1) % N
    er1, es1, ee1 = ecdsa_sign_raw(b"ec-one", d2, k1)
    er2, es2, ee2 = ecdsa_sign_raw(b"ec-two", d2, k2)
    assert er2 == Z1 * er1 % P
    den = (es1 * er2 - LAM * es2 * er1) % N
    num = (LAM * es2 * ee1 - es1 * ee2) % N
    d_rec2 = (num * inv_mod(den, N)) % N
    assert d_rec2 == d2
    out["f3_ecdsa_cm_nonce"] = {"r1": er1, "r2": er2, "r2_eq_zeta_r1": True,
                                "d_recovered": d_rec2, "correct": True,
                                "naive_duplicate_r_detector_fires": er1 == er2}
    print("F3b ecdsa: d recovered ==", d_rec2 == d2, "| naive dup-r detector fires:",
          er1 == er2, flush=True)

    lam_inv = inv_mod(LAM, N)
    assert mul(lam_inv, H) == G
    out["f4_dual_generator"] = {"H": H, "dlog_H_G": lam_inv,
                                "note": "publicly known -> dual-base/2-base schemes degenerate"}
    print("F4: dlog_H(G) = lambda^{-1} =", lam_inv, "(public)", flush=True)

    out["f5_x_only_closure"] = {"orbit": [rx, Z1 * rx % P, Z2 * rx % P],
                                 "size": len({rx, Z1 * rx % P, Z2 * rx % P})}
    print("F5 x-only orbit size:", out["f5_x_only_closure"]["size"], flush=True)

    os.makedirs(RESULTS, exist_ok=True)
    outpath = os.path.join(RESULTS, "05_out.json")
    with open(outpath, "w") as f:
        json.dump(out, f, indent=2)
    print(f"WROTE {outpath}", flush=True)


if __name__ == "__main__":
    main()
