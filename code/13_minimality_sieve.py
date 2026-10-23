"""
Minimality certificate for the constant c in

        p = 2^n - 2^s - c ,      c = 5 + 12 j ,   j = 0, 1, 2, ...

The recipe of the paper (Section 10) admits a candidate c only when

  (i)  p is prime, and
  (ii) exactly one of the six sextic-class orders  p + 1 - t  is prime.

This script performs an exhaustive search over j in increasing order and reports
the first admissible j.  It is instrumented so that the search is reproducible:

  Phase 1  segmented sieve of all j by the primes q <= SIEVE_BOUND:  a candidate is
           discarded as soon as q divides p(j) = 2^n - 2^s - 5 - 12 j.
  Phase 2  parallel Fermat (base 2) plus BPSW test of p on the survivors.
  Phase 3  Cornacchia decomposition p = A^2 + 3 B^2 of every prime p, the six CM
           traces, the six sextic-class orders, and their primality.
  Phase 4  point-multiplication certification of the winning order against the
           curve y^2 = x^3 + 7 (twelve random points).

Usage:
    python -u 13_minimality_sieve.py                # full 2560-bit search
    python -u 13_minimality_sieve.py --n 256 --s 32 --c 977   # pipeline validation
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

STATE = {}


def log(*parts):
    msg = " ".join(str(p) for p in parts)
    sys.stdout.write(msg + "\n")
    sys.stdout.flush()
    with open(LOGPATH, "a", encoding="utf-8") as fh:
        fh.write(msg + "\n")




def ec_add(P, Q, p, b):
    """Affine addition on y^2 = x^3 + b over F_p; None is the point at infinity."""
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2:
        if (y1 + y2) % p == 0:
            return None
        lam = (3 * x1 * x1) % p * pow(2 * y1 % p, p - 2, p) % p
    else:
        lam = (y2 - y1) % p * pow((x2 - x1) % p, p - 2, p) % p
    x3 = (lam * lam - x1 - x2) % p
    y3 = (lam * (x1 - x3) - y1) % p
    return (x3, y3)


def ec_mul(k, P, p, b):
    R = None
    Q = P
    while k:
        if k & 1:
            R = ec_add(R, Q, p, b)
        Q = ec_add(Q, Q, p, b)
        k >>= 1
    return R


def lift_x(x, p, b):
    """The point of the curve with abscissa x, or None if x^3 + b is not a square."""
    v = (x * x % p * x + b) % p
    s = pow(v, (p + 1) // 4, p) if p % 4 == 3 else pow(v, (p + 1) // 4, p)
    if (s * s - v) % p != 0:
        return None
    return (x, s)


def certify_order(p, b, order, tests=12, seed=20261022):
    """The order annihilates twelve random points of the curve."""
    import random

    rng = random.Random(seed)
    xs = []
    x = 1
    while len(xs) < tests:
        if lift_x(x, p, b) is not None:
            xs.append(x)
        x += 1
        if x > 4000:
            break
    for x in xs:
        P = lift_x(x, p, b)
        if ec_mul(order, P, p, b) is not None:
            return False, x
    return True, None




def primes_upto(bound):
    sieve = bytearray(b"\x01") * (bound + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, isqrt(bound) + 1):
        if sieve[i]:
            sieve[i * i :: i] = b"\x00" * ((bound - i * i) // i + 1)
    arr = np.frombuffer(sieve, dtype=np.uint8)
    ps = np.flatnonzero(arr[5:]) + 5
    return ps


def phase1_sieve(n, s, j_max, bound):
    """Mark every j whose p(j) has a prime factor <= bound."""
    t0 = time.time()
    flags = np.ones(j_max + 1, dtype=bool)
    ps = primes_upto(bound)
    log("  sieve: %d primes <= %d" % (len(ps), bound))
    two_n = 1 << n
    two_s = 1 << s
    for idx, q in enumerate(ps):
        q = int(q)
        inv12 = pow(12, -1, q)
        r = ((two_n % q - two_s % q - 5) * inv12) % q
        flags[r::q] = False
        if (idx + 1) % 250000 == 0:
            log("    sieved %d/%d primes (%.1fs)" % (idx + 1, len(ps), time.time() - t0))
    survivors = np.flatnonzero(flags)
    log("  sieve done in %.1fs: %d candidates, %d survivors (%.4f%%)"
        % (time.time() - t0, j_max + 1, len(survivors), 100.0 * len(survivors) / (j_max + 1)))
    return survivors


def _worker_p(args):
    """Fermat base 2, then BPSW, for the primality of p(j)."""
    j, n, s = args
    p = (1 << n) - (1 << s) - (5 + 12 * j)
    if pow(2, p - 1, p) != 1:
        return (j, 0)
    r = gmpy2.is_prime(p)
    return (j, int(r))


def cornacchia(p):
    """p = A^2 + 3 B^2 for a prime p = 1 mod 3."""
    for s in (pow(-3 % p, (p + 1) // 4, p), pow(pow(-3 % p, (p + 1) // 4, p), p - 2, p)):
        r0, r1 = p, s
        while r1 * r1 > p:
            r0, r1 = r1, r0 % r1
        if r1 * r1 <= p:
            x = r1
            y = isqrt((p - x * x) // 3)
            if x * x + 3 * y * y == p:
                return x, y
    raise ValueError("cornacchia failed")


def sextic_spectrum(p):
    """The six CM traces and the six sextic-class orders of a prime p = 1 mod 3."""
    A, B = cornacchia(p)
    traces = [2 * A, -2 * A, A + 3 * B, -(A + 3 * B), A - 3 * B, -(A - 3 * B)]
    orders = [p + 1 - t for t in traces]
    prime_mask = 0
    prime_orders = []
    for idx, o in enumerate(orders):
        if pow(2, o - 1, o) != 1:
            continue
        if gmpy2.is_prime(o) > 0:
            prime_mask |= 1 << idx
            prime_orders.append(o)
    return A, B, traces, orders, prime_mask, prime_orders


def _worker_sextic(args):
    j, n, s = args
    p = (1 << n) - (1 << s) - (5 + 12 * j)
    A, B, traces, orders, mask, prime_orders = sextic_spectrum(p)
    return {
        "j": j,
        "c": 5 + 12 * j,
        "cornacchia_A": A,
        "cornacchia_B": B,
        "traces": traces,
        "orders": orders,
        "prime_mask": mask,
        "prime_orders": prime_orders,
    }


def run(n, s, c_target, bound, workers, chunk):
    j_max = (c_target - 5) // 12
    total = j_max + 1
    log("=" * 72)
    log("SEARCH  n = %d, s = %d, c = %d  (j = 0 .. %d, %d candidates)"
        % (n, s, c_target, j_max, total))
    log("=" * 72)

    survivors = phase1_sieve(n, s, j_max, bound)
    t0 = time.time()
    scanned = 0
    hits = []
    done = False
    seen = set()

    with Pool(processes=workers) as pool:
        pos = 0
        while pos < len(survivors) and not done:
            part = survivors[pos : pos + chunk]
            args = [(int(j), n, s) for j in part]
            prime_js = []
            for j, isp in pool.imap(_worker_p, args, chunksize=64):
                scanned += 1
                if isp:
                    prime_js.append(j)
            pos += len(part)
            log("  chunk %d: %d/%d survivors scanned, %d prime p in chunk (%.1fs)"
                % (pos // chunk, scanned, len(survivors), len(prime_js), time.time() - t0))

            if prime_js:
                sargs = [(int(j), n, s) for j in sorted(prime_js)]
                for rec in pool.imap(_worker_sextic, sargs, chunksize=4):
                    if rec["j"] in seen:
                        continue
                    seen.add(rec["j"])
                    k = len(rec["prime_orders"])
                    if k == 1:
                        rec["recipe_hit"] = True
                        hits.append(rec)
                        log("  *** recipe hit at j = %d  c = %d  prime order = %d"
                            % (rec["j"], rec["c"], rec["prime_orders"][0]))
                        done = True
                        break
                    rec["recipe_hit"] = False
                    hits.append(rec)
                    if k > 1:
                        log("  note: j = %d has %d prime orders" % (rec["j"], k))

    elapsed = time.time() - t0
    log("  search finished in %.1fs" % elapsed)
    return survivors, hits, elapsed, total, j_max


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=2560)
    ap.add_argument("--s", type=int, default=320)
    ap.add_argument("--c", type=int, default=132200237)
    ap.add_argument("--bound", type=int, default=20000000)
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 4) - 2))
    ap.add_argument("--chunk", type=int, default=40000)
    ap.add_argument("--tag", type=str, default="")
    args = ap.parse_args()

    if os.path.exists(LOGPATH) and not args.tag:
        os.remove(LOGPATH)
    outpath = OUTPATH if not args.tag else OUTPATH.replace(".json", "_%s.json" % args.tag)
    log("minimality search started %s" % time.strftime("%Y-%m-%d %H:%M:%S"))
    log("python %s, numpy %s, gmpy2 %s, cpus %d"
        % (sys.version.split()[0], np.__version__, gmpy2.version(), os.cpu_count()))

    survivors, hits, elapsed, total, j_max = run(
        args.n, args.s, args.c, args.bound, args.workers, args.chunk)

    first = hits[0] if hits else None
    out = {
        "meta": {
            "generated_by": "code/13_minimality_sieve.py",
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "gmpy2": gmpy2.version(),
            "tag": args.tag,
        },
        "target": {
            "n": args.n,
            "s": args.s,
            "c": args.c,
            "j_max": j_max,
            "candidates": total,
        },
        "sieve": {
            "bound": args.bound,
            "survivors": int(len(survivors)),
            "eliminated": total - int(len(survivors)),
        },
        "hits": hits,
        "first_hit": first,
        "runtime_s": round(elapsed, 1),
    }

    if first is not None:
        j = first["j"]
        p = (1 << args.n) - (1 << args.s) - (5 + 12 * j)
        log("first admissible candidate: j = %d, c = %d" % (j, 5 + 12 * j))
        log("  p prime            : %s" % gmpy2.is_prime(p))
        log("  prime order        : %s" % first["prime_orders"])
        order = first["prime_orders"][0]
        ok, bad = certify_order(p, 7, order)
        log("  point certification of the order against y^2 = x^3 + 7: %s%s"
            % (ok, "" if ok else " (failed at x = %s)" % bad))
        out["certification"] = {
            "order_prime": bool(gmpy2.is_prime(order)),
            "order_annihilates_12_points": bool(ok),
            "sympy_isprime_order": None,
        }
        try:
            import sympy

            out["certification"]["sympy_isprime_order"] = bool(sympy.isprime(order))
            out["certification"]["sympy_isprime_p"] = bool(sympy.isprime(p))
            log("  sympy.isprime(order) = %s" % sympy.isprime(order))
        except Exception as exc:
            out["certification"]["sympy_error"] = str(exc)
        out["first_hit"]["c"] = 5 + 12 * j
        out["first_hit"]["p_dec"] = str(p)

    with open(outpath, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1)
    log("wrote %s" % outpath)


if __name__ == "__main__":
    main()
