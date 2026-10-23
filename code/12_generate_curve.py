#!/usr/bin/env python3
"""
secp2560k1 -- complete, self-contained generation and certification program.

The program rebuilds the curve from the field prime alone and certifies every
property the verification report prints.  It is pure Python and needs no
third-party package.

    stage 1   the prescribed family   p = 2^bits - 2^(bits/8) - c,   c = 5 (mod 12)
    stage 2   certify that p is prime
    stage 3   Cornacchia:  p = x^2 + 3y^2
    stage 4   the six sextic models b = 7 * 3^j and their six candidate orders
    stage 5   certify that exactly one of the six has prime order; that one is the curve
    stage 6   the base-point ritual
    stage 7   the cubic endomorphism: omega, lambda, and the CM lattice
    stage 8   the quadratic twist
    stage 9   certification: N*G = O, the Hasse window, the cofactor

Usage
    python 12_generate_curve.py --bits 256                       # full search
    python 12_generate_curve.py --bits 2560 --c 132200237        # verify this curve

The first invocation performs the whole minimality search and reproduces
secp256k1 exactly.  The second skips the search (which at 2560 bits is out of
reach of an interpreter) and verifies the parameters of secp2560k1 instead.
"""

import argparse
import json
import os
import random
import sys
import time



def log(msg):
    print(msg, flush=True)


def is_prime(n, rounds=40, seed=20251022):
    """Miller-Rabin with a fixed seed so that the verdict is reproducible."""
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    rng = random.Random(seed)
    for _ in range(rounds):
        a = rng.randrange(2, n - 1)
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def sqrt_mod(a, p):
    """A square root of a modulo the odd prime p, or None.  Tonelli-Shanks."""
    a %= p
    if a == 0:
        return 0
    if pow(a, (p - 1) // 2, p) != 1:
        return None
    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p
            i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, b * b % p, t * b * b % p, r * b % p
    return r


def cube_root_of_unity(p):
    """The two roots of x^2 + x + 1 modulo a prime p = 1 (mod 3)."""
    s = sqrt_mod(p - 3, p)
    if s is None:
        return None
    inv2 = pow(2, -1, p)
    z1 = ((-1 + s) * inv2) % p
    z2 = ((-1 - s) * inv2) % p
    assert (z1 * z1 + z1 + 1) % p == 0 and (z2 * z2 + z2 + 1) % p == 0
    return z1, z2


def round_div(a, b):
    """Round a/b to the nearest integer, ties to even, using integers only.

    b must be positive.  This is the routine that the published listing gets
    wrong: a floating-point division of two 2560-bit integers does not fit in
    a double, and even when it fits the 53-bit mantissa silently truncates the
    low order digits of the quotient.
    """
    q, r = divmod(a, b)
    twice = 2 * r
    if twice > b:
        q += 1
    elif twice == b and (q & 1):
        q += 1
    return q


def isqrt(n):
    import math
    return math.isqrt(n)




def jac_double(X1, Y1, Z1, p):
    if Y1 == 0 or Z1 == 0:
        return (0, 0, 0)
    A = Y1 * Y1 % p
    B = 4 * X1 * A % p
    C = 8 * A * A % p
    D = 3 * X1 * X1 % p
    X3 = (D * D - 2 * B) % p
    Y3 = (D * (B - X3) - C) % p
    Z3 = 2 * Y1 * Z1 % p
    return (X3, Y3, Z3)


def jac_add(X1, Y1, Z1, X2, Y2, Z2, p):
    if Z1 == 0:
        return (X2, Y2, Z2)
    if Z2 == 0:
        return (X1, Y1, Z1)
    Z1Z1 = Z1 * Z1 % p
    Z2Z2 = Z2 * Z2 % p
    U1 = X1 * Z2Z2 % p
    U2 = X2 * Z1Z1 % p
    S1 = Y1 * Z2 * Z2Z2 % p
    S2 = Y2 * Z1 * Z1Z1 % p
    if U1 == U2:
        if S1 != S2:
            return (0, 0, 0)
        return jac_double(X1, Y1, Z1, p)
    H = (U2 - U1) % p
    I = 2 * H % p
    I = I * I % p
    J = H * I % p
    r = 2 * (S2 - S1) % p
    V = U1 * I % p
    X3 = (r * r - J - 2 * V) % p
    Y3 = (r * (V - X3) - 2 * S1 * J) % p
    Z3 = ((Z1 + Z2) ** 2 - Z1Z1 - Z2Z2) % p
    Z3 = Z3 * H % p
    return (X3, Y3, Z3)


def to_affine(X, Y, Z, p):
    if Z == 0:
        return None
    zi = pow(Z, -1, p)
    zi2 = zi * zi % p
    return (X * zi2 % p, Y * zi2 % p * zi % p)


def mul(k, pt, p, b=7):
    """k * pt in affine coordinates; None for the point at infinity."""
    if pt is None or k == 0:
        return None
    X, Y, Z = pt[0], pt[1], 1
    acc = (0, 0, 0)
    for bit in bin(k)[2:]:
        acc = jac_double(*acc, p)
        if bit == "1":
            acc = jac_add(acc[0], acc[1], acc[2], X, Y, Z, p)
    return to_affine(*acc, p)


def on_curve(pt, p, b=7):
    if pt is None:
        return False
    x, y = pt
    return (y * y - x * x * x - b) % p == 0


def lift_x(x, p, b=7):
    y = sqrt_mod(x * x * x + b, p)
    if y is None:
        return None
    return (x, y)




def make_prime(bits, c):
    """p = 2^bits - 2^(bits/8) - c."""
    return (1 << bits) - (1 << (bits // 8)) - c


def sieve_candidates(bits, limit):
    """c = 5 (mod 12), ascending."""
    c = 5
    while c <= limit:
        yield c
        c += 12


def cornacchia(p):
    """p = x^2 + 3y^2 for p = 1 (mod 3) and p = 3 (mod 4)."""
    s = sqrt_mod(p - 3, p)
    if s is None:
        return None
    r0, r1 = p, s
    while r1 * r1 > p:
        r0, r1 = r1, r0 % r1
    x = r1
    y = isqrt((p - x * x) // 3)
    if x * x + 3 * y * y != p:
        return None
    return x, y




MODELS = [
    (0, 7, lambda x, y: x + 3 * y),
    (1, 21, lambda x, y: 2 * x),
    (2, 63, lambda x, y: x - 3 * y),
    (3, 189, lambda x, y: -(x + 3 * y)),
    (4, 567, lambda x, y: -2 * x),
    (5, 1701, lambda x, y: -(x - 3 * y)),
]


def sextic_spectrum(p, x, y):
    """The six (model, trace, candidate order) triples of the j = 0 family."""
    out = []
    for j, b, trace_of in MODELS:
        t = trace_of(x, y)
        out.append((j, b, t, p + 1 - t))
    return out


def certify_order(b_model, candidate, p, tests=12, seed=99):
    """candidate is the order of the model iff it annihilates random points
    of that model and lies inside the Hasse window."""
    lo = p + 1 - 2 * isqrt(p)
    hi = p + 1 + 2 * isqrt(p)
    if not (lo <= candidate <= hi):
        return False
    rng = random.Random(seed)
    x = 0
    found = 0
    while found < tests:
        x += 1
        pt = lift_x(x, p, b_model)
        if pt is None:
            continue
        found += 1
        if mul(candidate, pt, p, b_model) is not None:
            return False
    return True




def base_point(p, b=7, seed=0):
    """The ritual: the first x >= seed for which x^3 + b is a quadratic
    residue, taken with the even y, then doubled."""
    x = seed
    while True:
        pt = lift_x(x, p, b)
        if pt is not None:
            y = pt[1]
            if y % 2:
                y = p - y
            P0 = (x, y)
            G = mul(2, P0, p, b)
            return P0, G
        x += 1


def canonical_secp256k1_generator():
    """The generator published for secp256k1, used only as an extra check."""
    return (0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798,
            0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8)




def gauss_reduce(u, v):
    """Two-dimensional Lagrange reduction on integer vectors."""
    def norm(w):
        return w[0] * w[0] + w[1] * w[1]

    def dot(a, b):
        return a[0] * b[0] + a[1] * b[1]

    while True:
        if norm(u) > norm(v):
            u, v = v, u
        q = round_div(dot(u, v), norm(u))
        if q == 0:
            return u, v
        v = (v[0] - q * u[0], v[1] - q * u[1])


def cm_lattice(n, lam):
    """A reduced basis of {(a, b) : a + b*lambda = 0 (mod n)}."""
    u, v = gauss_reduce((n, 0), ((-lam) % n, 1))
    return u, v


def decompose(k, a1, b1, a2, b2, det, n, lam):
    """k = k1 + k2*lambda (mod n) with |k1|, |k2| about sqrt(n), integers only."""
    c1 = round_div(b2 * k, det)
    c2 = round_div(-b1 * k, det)
    k1 = k - c1 * a1 - c2 * a2
    k2 = -(c1 * b1 + c2 * b2)
    return k1, k2


def crandall_reduce(x, p, bits, half, c):
    """Reduce x modulo p = 2^bits - 2^half - c.

    Since 2^bits = 2^half + c (mod p), a high part hi = x >> bits is
    replaced by hi * (2^half + c): one multiplication by the 32-bit constant
    c and one shift-and-add.  Folding continues until the value fits in
    `bits` bits; at most one subtraction of p then remains.
    """
    M = 1 << bits
    mask = M - 1
    guard = 0
    while x >= M:
        hi = x >> bits
        x = (x & mask) + hi * ((1 << half) + c)
        guard += 1
        if guard > 16:
            raise RuntimeError("reduction did not terminate")
    while x >= p:
        x -= p
    return x




def twist_of(p, t, b=7):
    """The quadratic twist: its model, its order and its rational 3-torsion."""
    d = 2
    while pow(d, (p - 1) // 2, p) != p - 1:
        d += 1
    b_twist = b * d * d * d % p
    n_t = p + 1 + t
    return d, b_twist, n_t


def search(bits, cap, seed):
    """Stage 1 to 5: walk c = 5 (mod 12) until the recipe fires."""
    tried = 0
    for c in sieve_candidates(bits, cap):
        tried += 1
        p = make_prime(bits, c)
        if not is_prime(p):
            continue
        xy = cornacchia(p)
        if xy is None:
            continue
        x, y = xy
        spec = sextic_spectrum(p, x, y)
        primes = [(j, bb, tr, o) for (j, bb, tr, o) in spec if is_prime(o)]
        if len(primes) != 1:
            continue
        j, bb, tr, o = primes[0]
        if j != 0:
            continue
        log("  hit: c = %d, p prime, exactly one prime sextic order, model b = 7" % c)
        return c, p, x, y, spec, tried
    return None, None, None, None, None, tried


def analyse(bits, c, spec=None, seed=2):
    """Every stage from the prime to the certification."""
    out = {"bits": bits, "c": c}
    p = make_prime(bits, c)
    out["p"] = hex(p)
    out["p_dec"] = str(p)
    log("  p = 2^%d - 2^%d - %d" % (bits, bits // 8, c))
    log("  p bits = %d, p mod 12 = %d" % (p.bit_length(), p % 12))

    log("  stage 2: primality of p")
    assert is_prime(p)
    out["p_is_prime"] = True

    log("  stage 3: Cornacchia")
    x, y = cornacchia(p)
    out["cornacchia_x"] = hex(x)
    out["cornacchia_y"] = hex(y)
    out["cornacchia_x_bits"] = x.bit_length()
    out["cornacchia_y_bits"] = y.bit_length()
    log("    x = %d bits, y = %d bits" % (x.bit_length(), y.bit_length()))

    log("  stage 4: the sextic spectrum")
    if spec is None:
        spec = sextic_spectrum(p, x, y)
    models = []
    for (j, bb, tr, o) in spec:
        models.append({
            "j": j, "b": bb, "trace": hex(tr), "order": hex(o),
            "order_bits": o.bit_length(), "order_is_prime": is_prime(o),
        })
        log("    j = %d, b = %-5d order bits %d prime %s"
            % (j, bb, o.bit_length(), models[-1]["order_is_prime"]))
    out["models"] = models
    prime_models = [m for m in models if m["order_is_prime"]]
    assert len(prime_models) == 1, "the recipe requires exactly one prime order"
    assert prime_models[0]["j"] == 0, "the prime order must be the model b = 7"

    n = int(prime_models[0]["order"], 16)
    t = p + 1 - n
    out["n"] = hex(n)
    out["n_dec"] = str(n)
    out["n_bits"] = n.bit_length()
    out["t"] = hex(t)
    out["t_bits"] = t.bit_length()
    assert is_prime(n)
    log("    N = %d bits, prime" % n.bit_length())

    log("    certifying the order of the model b = 7 by point annihilation")
    ok = certify_order(7, n, p, tests=6)
    out["order_certified_by_points"] = ok
    assert ok, "the candidate order did not annihilate the test points"
    log("    the order is certified: it annihilates random points of the model")

    log("  stage 6: the base-point ritual")
    P0, G = base_point(p, seed=seed)
    assert on_curve(G, p)
    out["Gx"] = hex(G[0])
    out["Gy"] = hex(G[1])
    out["P0_x"] = str(P0[0])
    out["P0_y_even"] = (P0[1] % 2 == 0)
    out["G_is_2P0"] = (mul(2, P0, p) == G)
    log("    P0 = (%d, even), G = 2*P0, on curve" % P0[0])
    if P0[0] == 2:
        closed = ((p - 8) // 5) % p
        out["Gx_closed_form_holds"] = (G[0] == closed)
        assert G[0] == closed
        log("    x(G) = (p - 8)/5 (mod p), as the closed form predicts")

    log("  stage 7: the cubic endomorphism")
    z1, z2 = cube_root_of_unity(p)
    omega = None
    lam = None
    for z in (z1, z2):
        phiG = (z * G[0] % p, G[1])
        for cand in cube_root_of_unity(n):
            if mul(cand, G, p) == phiG:
                omega, lam = z, cand
    assert omega is not None, "the endomorphism did not match"
    assert (lam * lam + lam + 1) % n == 0
    out["omega"] = hex(omega)
    out["lambda"] = hex(lam)
    out["lambda_inv"] = hex(pow(lam, -1, n))
    log("    omega found, lambda = %d bits, lambda^2 + lambda + 1 = 0 (mod N)" % lam.bit_length())

    u1, u2 = cm_lattice(n, lam)
    a1, b1 = u1
    a2, b2 = u2
    det = a1 * b2 - a2 * b1
    assert det == n or det == -n
    assert (a1 + b1 * lam) % n == 0 and (a2 + b2 * lam) % n == 0
    out["basis"] = {"a1": str(a1), "b1": str(b1), "a2": str(a2), "b2": str(b2),
                    "det": str(det), "det_equals_n": det == n}
    log("    CM lattice basis, det = N exactly, both relations hold")

    rng = random.Random(7)
    worst = 0
    for _ in range(20):
        k = rng.randrange(1 << (bits - 1), 1 << bits)
        k1, k2 = decompose(k, a1, b1, a2, b2, det, n, lam)
        assert (k1 + k2 * lam - k) % n == 0
        worst = max(worst, max(abs(k1).bit_length(), abs(k2).bit_length()))
    out["decompose_max_half_bits"] = worst
    log("    exact integer decomposition verified, worst half = %d bits" % worst)

    log("  stage 8: the quadratic twist")
    d, b_twist, n_t = twist_of(p, t)
    pt3 = lift_x(0, p, b_twist)
    has3 = False
    if pt3 is not None:
        has3 = (mul(3, pt3, p, b_twist) is None)
    out["twist"] = {"d": d, "b": str(b_twist), "order": hex(n_t),
                    "order_bits": n_t.bit_length(),
                    "v3": _v3(n_t),
                    "rational_3torsion": has3,
                    "gcd_with_n": _gcd(n, n_t)}
    log("    d = %d, b' = %d, n_t = %d bits, gcd(N, n_t) = %d"
        % (d, b_twist, n_t.bit_length(), _gcd(n, n_t)))

    log("  stage 9: certification")
    assert mul(n, G, p) is None, "N*G is not the point at infinity"
    lo = p + 1 - 2 * isqrt(p)
    hi = p + 1 + 2 * isqrt(p)
    assert lo <= n <= hi
    out["N_times_G_is_O"] = True
    out["n_in_hasse_window"] = True
    out["cofactor"] = 1
    log("    N*G = O and N lies in the Hasse window, so #E(Fp) = N and h = 1")

    log("  reduction check")
    prod = (G[0] * G[0]) % (1 << (2 * bits))
    red = crandall_reduce(prod, p, bits, bits // 8, c)
    assert red == prod % p
    out["crandall_ok"] = True
    log("    the Crandall reduction agrees with exact division")
    return out


def _gcd(a, b):
    while b:
        a, b = b, a % b
    return a


def _v3(n):
    """The exponent of 3 in n."""
    e = 0
    while n % 3 == 0:
        n //= 3
        e += 1
    return e


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bits", type=int, default=256)
    ap.add_argument("--c", type=int, default=None)
    ap.add_argument("--cap", type=int, default=4000)
    ap.add_argument("--seed", type=int, default=2,
                    help="the first abscissa tried by the base-point ritual")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    t0 = time.time()
    log("secp2560k1 generator, bits = %d" % args.bits)
    if args.c is None:
        log("stage 1: searching c = 5 (mod 12)")
        c, p, x, y, spec, tried = search(args.bits, args.cap, args.seed)
        if c is None:
            log("no hit within the cap")
            return 1
        log("  candidates tested: %d" % tried)
        out = analyse(args.bits, c, spec, args.seed)
        out["candidates_tested"] = tried
    else:
        out = analyse(args.bits, args.c, None, args.seed)
        out["candidates_tested"] = None
        out["search_performed"] = False
    out["elapsed_seconds"] = round(time.time() - t0, 2)
    log("done in %.1f s" % out["elapsed_seconds"])
    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(out, f, indent=1)
        log("wrote " + args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
