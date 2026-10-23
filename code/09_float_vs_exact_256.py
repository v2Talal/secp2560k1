"""256-bit contrast: the SAME published listing does NOT crash, but silently loses precision.
Uses the already-certified secp256k1 CM lattice basis (script 03).
Run: python -u 09_float_vs_exact_256.py
"""
import os, random
log = lambda *a: print(*a, flush=True)

N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
lam = 78074008874160198520644763525212887401909906723592317393988542598630163514318
u1 = (303414439467246543595250775667605759171, -64502973549206556628585045361533709077)
u2 = (64502973549206556628585045361533709077, 367917413016453100223835821029139468248)
a1, b1, a2, b2 = u1[0], u1[1], u2[0], u2[1]
det = a1 * b2 - a2 * b1
log("256-bit basis: det == +-N :", abs(det) == N, "| det bits:", abs(det).bit_length())

def rnd(num, den):
    q, r = divmod(num, den)
    if 2 * r > den or (2 * r == den and (q & 1)):
        q += 1
    return q

rnd_ = random.Random(20261022)
bad = 0
worst = 0
for _ in range(200):
    k = rnd_.getrandbits(256) % N
    c1f = round(b2 * k / det)
    c2f = round(-b1 * k / det)
    c1e = rnd(b2 * k, abs(det))
    c2e = rnd(-b1 * k, abs(det))
    if (c1f, c2f) != (c1e, c2e):
        bad += 1
        worst = max(worst, abs(c1f - c1e))
log("256-bit: cases where float rounding != exact rounding :", bad, "/ 200")
log("256-bit: worst coefficient error (in units of 1) :", worst)
log("  -> float mantissa is 53 bits, but the quotient here is ~2^128 bits: the")
log("     low ~75 bits are simply garbage. No exception is raised.")
log("  -> at 2560 bits the same expression needs ~2^1278 => OverflowError instead.")
log("DONE")
