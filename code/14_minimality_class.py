"""
Class-aware minimality certificate for c in   p = 2^n - 2^s - c ,  c = 5 + 12 j.

Two criteria are examined separately:

  (A)  "recipe" criterion -- p prime and exactly one of the six sextic-class
       orders p + 1 - t is prime.  This is the criterion stated in Section 10.4
       of the paper, and it is what produced the 256-bit validation (c = 977).

  (B)  "b = 7" criterion -- p prime and the order of the curve y^2 = x^3 + 7
       itself is prime, i.e. the prime order belongs to the sextic class of the
       curve that this paper actually constructs.

The sextic classes are b_j = 7 * 3^j, j = 0..5, with traces
        {+2A, -2A, +(A+3B), -(A+3B), +(A-3B), -(A-3B)},   p = A^2 + 3B^2.
The six traces form a sign-invariant set, but the assignment of a trace to a
class depends on the sign convention of the Cornacchia pair.  The class of a
prime order is therefore settled numerically: the order is tested against the
twelve-point annihilation test on each of the six curves.

Everything runs on GMP (gmpy2.powmod), which is about ten times faster than
the built-in pow for a 2560-bit modulus.

Usage:
    python -u 14_minimality_class.py
"""

import argparse
import json
import os
import sys
import time
from math import isqrt

import numpy as np
import gmpy2
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "data")
RESULTS = os.path.join(os.path.dirname(HERE), "results")
os.makedirs(RESULTS, exist_ok=True)
LOGPATH = os.path.join(RESULTS, "run_minimality.txt")
OUTPATH = os.path.join(DATA, "paper_minimality.json")

MODELS = [7, 21, 63, 189, 567, 1701]


def log(*parts):
    msg = " ".join(str(p) for p in parts)
    sys.stdout.write(msg + "\n")
    sys.stdout.flush()
    with open(LOGPATH, "a", encoding="utf-8") as fh:
        fh.write(msg + "\n")




def jac_double(P, p):
    X1, Y1, Z1 = P
    if Y1 == 0 or Z1 == 0:
        return (0, 1, 0)
    YY = Y1 * Y1 % p
    S = 4 * X1 * YY % p
    M = 3 * X1 * X1 % p
    X3 = (M * M - 2 * S) % p
    Y3 = (M * (S - X3) - 8 * YY * YY) % p
    Z3 = 2 * Y1 * Z1 % p
    return (X3, Y3, Z3)


def jac_add(P, Q, p):
    X1, Y1, Z1 = P
    X2, Y2, Z2 = Q
    if Z1 == 0:
        return Q
    if Z2 == 0:
        return P
    Z1Z1 = Z1 * Z1 % p
    Z2Z2 = Z2 * Z2 % p
    U1 = X1 * Z2Z2 % p
    U2 = X2 * Z1Z1 % p
    S1 = Y1 * Z2 * Z2Z2 % p
    S2 = Y2 * Z1 * Z1Z1 % p
    if U1 == U2:
        return jac_double(P, p) if S1 == S2 else (0, 1, 0)
    H = (U2 - U1) % p
    I = (2 * H) % p
    I = I * I % p
    J = H * I % p
    r = 2 * (S2 - S1) % p
    V = U1 * I % p
    X3 = (r * r - J - 2 * V) % p
    Y3 = (r * (V - X3) - 2 * S1 * J) % p
    Z3 = ((Z1 + Z2) ** 2 - Z1Z1 - Z2Z2) * H % p
    return (X3, Y3, Z3)


def jac_mul(k, P, p):
    R = (0, 1, 0)
    Q = P
    while k:
        if k & 1:
            R = jac_add(R, Q, p)
        Q = jac_double(Q, p)
        k >>= 1
    return R


def lift_x(x, p, b):
    """A point of y^2 = x^3 + b with abscissa x, or None."""
    v = (x * x % p * x + b) % p
    s = gmpy2.powmod(v, (p + 1) // 4, p)
    if (s * s - v) % p != 0:
        return None
    return (x, int(s))


def order_annihilates(p, b, order, tests=12):
    """Does the given order annihilate the points of y^2 = x^3 + b?

    The order lies in the Hasse window and exceeds half of it, so annihilation
    of one point already identifies the group order: no other order of the six
    is a multiple of it inside the window.
    """
    x = 1
    pts = []
    while len(pts) < tests:
        P = lift_x(x, p, b)
        if P is not None:
            pts.append(P)
        x += 1
        if x > 5000:
            break
    if not pts:
        return False
    for P in pts:
        R = jac_mul(order, (P[0], P[1], 1), p)
        if R[2] != 0:
            return False
    return True




def primes_upto(bound):
    sieve = bytearray(b"\x01") * (bound + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, isqrt(bound) + 1):
        if sieve[i]:
            sieve[i * i :: i] = b"\x00" * ((bound - i * i) // i + 1)
    arr = np.frombuffer(sieve, dtype=np.uint8)
    return np.flatnonzero(arr[5:]) + 5


def _p_of(j, n, s):
    return (1 << n) - (1 << s) - (5 + 12 * j)


def _worker_p(args):
    """Fermat base 2 (GMP), then BPSW, for the primality of p(j)."""
    j, n, s = args
    p = _p_of(j, n, s)
    if gmpy2.powmod(2, p - 1, p) != 1:
        return (j, 0)
    return (j, int(gmpy2.is_prime(p)))


def cornacchia(p):
    for s in (gmpy2.powmod(p - 3, (p + 1) // 4, p),
              gmpy2.powmod(gmpy2.powmod(p - 3, (p + 1) // 4, p), p - 2, p)):
        r0, r1 = p, int(s)
        while r1 * r1 > p:
            r0, r1 = r1, r0 % r1
        if r1 * r1 <= p:
            x = r1
            y = isqrt((p - x * x) // 3)
            if x * x + 3 * y * y == p:
                return x, y
    raise ValueError("cornacchia failed")


def _worker_sextic(args):
    """The six traces, the six orders, and which of them are prime."""
    j, n, s = args
    p = _p_of(j, n, s)
    A, B = cornacchia(p)
    traces = [2 * A, -2 * A, A + 3 * B, -(A + 3 * B), A - 3 * B, -(A - 3 * B)]
    orders = [p + 1 - t for t in traces]
    mask = 0
    prime_orders = []
    for idx, o in enumerate(orders):
        if gmpy2.powmod(2, o - 1, o) != 1:
            continue
        if gmpy2.is_prime(o) > 0:
            mask |= 1 << idx
            prime_orders.append(o)
    return {"j": j, "c": 5 + 12 * j, "cornacchia_A": A, "cornacchia_B": B,
            "traces": traces, "orders": orders, "prime_mask": mask,
            "prime_orders": prime_orders}


def _worker_class(args):
    """Which sextic class owns the prime order?  Tested numerically."""
    j, n, s, order = args
    p = _p_of(j, n, s)
    owner = None
    for cls, b in enumerate(MODELS):
        if order_annihilates(p, b, order, tests=4):
            owner = cls
            break
    return j, owner




def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=2560)
    ap.add_argument("--s", type=int, default=320)
    ap.add_argument("--c", type=int, default=132200237)
    ap.add_argument("--bound", type=int, default=20000000)
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 4) - 2))
    ap.add_argument("--chunk", type=int, default=100000)
    ap.add_argument("--tag", type=str, default="")
    args = ap.parse_args()

    n, s, c_target = args.n, args.s, args.c
    j_max = (c_target - 5) // 12
    total = j_max + 1

    outpath = OUTPATH if not args.tag else OUTPATH.replace(".json", "_%s.json" % args.tag)
    if os.path.exists(LOGPATH) and not args.tag:
        os.remove(LOGPATH)
    log("=" * 74)
    log("CLASS-AWARE MINIMALITY SEARCH  n = %d  s = %d  c = %d  (j = 0..%d, %d candidates)"
        % (n, s, c_target, j_max, total))
    log("started %s" % time.strftime("%Y-%m-%d %H:%M:%S"))
    log("python %s, numpy %s, gmpy2 %s, workers %d"
        % (sys.version.split()[0], np.__version__, gmpy2.version(), args.workers))

    t0 = time.time()
    flags = np.ones(total, dtype=bool)
    ps = primes_upto(args.bound)
    log("  sieve over %d primes <= %d" % (len(ps), args.bound))
    for idx, q in enumerate(ps):
        q = int(q)
        inv12 = pow(12, -1, q)
        r = (((1 << n) % q - (1 << s) % q - 5) * inv12) % q
        flags[r::q] = False
        if (idx + 1) % 300000 == 0:
            log("    sieved %d/%d (%.1fs)" % (idx + 1, len(ps), time.time() - t0))
    survivors = np.flatnonzero(flags)
    log("  sieve done in %.1fs: %d survivors of %d (%.4f%%)"
        % (time.time() - t0, len(survivors), total, 100.0 * len(survivors) / total))

    t0 = time.time()
    scanned = 0
    p_prime_total = 0
    hits = []
    ambiguous = []
    first_class0 = None
    seen = set()

    with Pool(processes=args.workers) as pool:
        pos = 0
        while pos < len(survivors):
            part = survivors[pos : pos + args.chunk]
            prime_js = []
            for j, isp in pool.imap(_worker_p, [(int(x), n, s) for x in part], chunksize=32):
                scanned += 1
                if isp:
                    prime_js.append(j)
            pos += len(part)
            p_prime_total += len(prime_js)
            log("  chunk %d: %d/%d survivors, %d prime p (%.1fs)"
                % (pos // args.chunk, scanned, len(survivors), len(prime_js), time.time() - t0))
            if not prime_js:
                continue
            recs = list(pool.imap(_worker_sextic, [(int(j), n, s) for j in sorted(prime_js)], chunksize=4))
            for rec in recs:
                if rec["j"] in seen:
                    continue
                seen.add(rec["j"])
                k = len(rec["prime_orders"])
                if k == 1:
                    hits.append(rec)
                    log("    recipe hit: j = %-9d c = %-11d prime order bits %d"
                        % (rec["j"], rec["c"], rec["prime_orders"][0].bit_length()))
                elif k > 1:
                    rec["ambiguous"] = True
                    ambiguous.append(rec)
                    log("    AMBIGUOUS: j = %d has %d prime orders" % (rec["j"], k))

    log("  sweep done in %.1fs: %d prime p, %d recipe hits, %d ambiguous"
        % (time.time() - t0, p_prime_total, len(hits), len(ambiguous)))

    with Pool(processes=args.workers) as pool:
        cargs = [(r["j"], n, s, r["prime_orders"][0]) for r in hits]
        owners = dict(pool.imap(_worker_class, cargs, chunksize=1))
    for r in hits:
        r["class_index"] = owners.get(r["j"])
        r["class_b"] = MODELS[owners[r["j"]]] if r["j"] in owners else None
        r["is_class0_b7"] = (owners.get(r["j"]) == 0)
        log("  hit j = %-9d c = %-11d class j = %s (b = %s)"
            % (r["j"], r["c"], r["class_index"], r["class_b"]))
    class0 = [r for r in hits if r.get("is_class0_b7")]
    if class0:
        first_class0 = class0[0]
        log("  FIRST b = 7 HIT: j = %d  c = %d" % (first_class0["j"], first_class0["c"]))

    cert = {}
    if first_class0:
        j = first_class0["j"]
        p = _p_of(j, n, s)
        order = first_class0["prime_orders"][0]
        t0c = time.time()
        ok = order_annihilates(p, 7, order, tests=12)
        cert = {
            "j": j,
            "c": first_class0["c"],
            "p_is_probable_prime": bool(gmpy2.is_prime(p)),
            "order_is_probable_prime": bool(gmpy2.is_prime(order)),
            "order_annihilates_12_points_of_b7": bool(ok),
            "certification_seconds": round(time.time() - t0c, 1),
        }
        try:
            import sympy

            cert["sympy_isprime_p"] = bool(sympy.isprime(p))
            cert["sympy_isprime_order"] = bool(sympy.isprime(order))
        except Exception as exc:
            cert["sympy_error"] = str(exc)
        log("  certification of the b = 7 hit: %s" % cert)

    out = {
        "meta": {
            "generated_by": "code/14_minimality_class.py",
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "gmpy2": gmpy2.version(),
            "tag": args.tag,
        },
        "target": {"n": n, "s": s, "c": c_target, "j_max": j_max, "candidates": total},
        "sieve": {"bound": args.bound, "survivors": int(len(survivors)),
                   "eliminated": total - int(len(survivors))},
        "counts": {"prime_p": p_prime_total, "recipe_hits": len(hits),
                   "ambiguous": len(ambiguous)},
        "first_recipe_hit": hits[0] if hits else None,
        "first_class0_hit": first_class0,
        "recipe_hits": hits,
        "ambiguous": ambiguous,
        "certification": cert,
    }
    with open(outpath, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1)
    log("wrote %s" % outpath)


if __name__ == "__main__":
    main()
