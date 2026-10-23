"""Exact vs float GLV decompose(): why CORE LISTING 4 breaks (and how).
Run: python -u 08_decompose_check.py
NOTE: JSON glv_basis values are DECIMAL strings (signed); omega/lambda are hex.
"""
import json, os

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "data")
log = lambda *a: print(*a, flush=True)

J = json.load(open(os.path.join(DATA, "secp2560k1.json")))
N = int(J["n_hex"], 16)
g = J["glv_basis"]
a1, b1, a2, b2 = int(g["a1"]), int(g["b1"]), int(g["a2"]), int(g["b2"])
det = a1 * b2 - a2 * b1
lam = int(J["glv"]["lambda_hex"], 16)
omega = int(J["glv"]["omega_hex"], 16)
log("=== 2560-bit curve (secp2560k1) ===")
log("N bits:", N.bit_length(), "| basis bits:", a1.bit_length(), b1.bit_length(), a2.bit_length(), b2.bit_length())
log("det == +-N :", abs(det) == N, " sign:", (det > 0) - (det < 0))
log("a1 + b1*lambda == 0 mod N :", (a1 + b1 * lam) % N == 0)
log("a2 + b2*lambda == 0 mod N :", (a2 + b2 * lam) % N == 0)

k = int.from_bytes(os.urandom(320), "big") % N
log("k bits :", k.bit_length(), " (full-size 2560-bit scalar)")

def published_decompose(k, a1=a1, b1=b1, a2=a2, b2=b2, det=det):
    c1 = round(b2 * k / det)
    c2 = round(-b1 * k / det)
    return k - c1 * a1 - c2 * a2, -(c1 * b1 + c2 * b2)

log("")
log("--- published listing, verbatim, on a real 2560-bit scalar ---")
try:
    r = published_decompose(k)
    log("no exception; k1 bits =", abs(r[0]).bit_length(), " k2 bits =", abs(r[1]).bit_length())
    log("congruence holds ?", (r[0] + r[1] * lam - k) % N == 0)
except OverflowError as e:
    log("OverflowError:", e)
    log("why: b2*k ~ 2^%d ; det ~ 2^%d ; quotient ~ 2^%d" %
        ((b2 * k).bit_length(), abs(det).bit_length(), (b2 * k).bit_length() - abs(det).bit_length()))
    log("     Python float max ~ 2^1024  =>  int/int must become float  =>  OverflowError")

def rnd(num, den):
    """round-half-even of rational num/den, den > 0, EXACT integer arithmetic."""
    q, r = divmod(num, den)
    if 2 * r > den or (2 * r == den and (q & 1)):
        q += 1
    return q

def exact_decompose(k):
    d = abs(det); s = 1 if det > 0 else -1
    c1 = rnd(s * b2 * k, d)
    c2 = rnd(-s * b1 * k, d)
    return k - c1 * a1 - c2 * a2, -(c1 * b1 + c2 * b2)

log("")
log("--- exact integer version (no floats) ---")
k1, k2 = exact_decompose(k)
log("k1 bits :", abs(k1).bit_length(), "  k2 bits:", abs(k2).bit_length())
log("k1 + k2*lambda == k (mod N) :", (k1 + k2 * lam - k) % N == 0)
log("paper's claimed decomp_bits [1279,1277] vs ours:",
    [abs(k1).bit_length(), abs(k2).bit_length()])

log("")
log("--- same listing at 256-bit scale (secp256k1 numbers) ---")
N256 = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
lam256 = 78074008874160198520644763525212887401909906723592317393988542598630163514318

def gauss2(u, v):
    while True:
        if abs(u[0]) < abs(v[0]):
            u, v = v, u
        if v[0] == 0:
            return u, v
        q = (2 * u[0] + v[0]) // (2 * v[0])
        u, v = (u[0] - q * v[0], u[1] - q * v[1]), u
        if abs(v[0]) <= abs(u[0]):
            return u, v

u1, u2 = gauss2((N256, 0), ((-lam256) % N256, 1))
a1b, b1b, a2b, b2b = u1[0], u1[1], u2[0], u2[1]
detb = a1b * b2b - a2b * b1b
log("256-bit basis det == +-N :", abs(detb) == N256)
k256 = int.from_bytes(os.urandom(32), "big") % N256
try:
    c1 = round(b2b * k256 / detb)
    c2 = round(-b1b * k256 / detb)
    k1f, k2f = k256 - c1 * a1b - c2 * a2b, -(c1 * b1b + c2 * b2b)
    log("256-bit: NO exception (float max 2^1024 >> 2^256) -> it silently 'works'")
    log("256-bit: float result satisfies k1+k2*lambda == k (mod N) ?",
        (k1f + k2f * lam256 - k256) % N256 == 0)
    c1e = rnd(b2b * k256, abs(detb))
    c2e = rnd(-b1b * k256, abs(detb))
    k1e, k2e = k256 - c1e * a1b - c2e * a2b, -(c1e * b1b + c2e * b2b)
    log("256-bit: EXACT integer result satisfies the congruence ?",
        (k1e + k2e * lam256 - k256) % N256 == 0)
    log("256-bit: float coefficients == exact coefficients ?", (c1, c2) == (c1e, c2e))
    log("   float c1 =", c1, "| exact c1 =", c1e, "| float mantissa = 53 bits")
    log("   float k1 bits:", abs(k1f).bit_length(), "| exact k1 bits:", abs(k1e).bit_length())
except OverflowError as e:
    log("256-bit: OverflowError:", e)

log("")
log("DONE")
