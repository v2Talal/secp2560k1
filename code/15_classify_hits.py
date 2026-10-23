"""
Classify the 52 recipe hits of the 2560-bit minimality sweep.

The upstream search ran its Pool and lost the per-hit classification, so
it is redone here in a single process, with incremental output.

For every hit we
  * rebuild p, the Cornacchia pair (A, B) and the six sextic traces,
  * rebuild the six orders and pick the unique prime one,
  * identify the sextic class that owns that prime order by annihilating points
    of y^2 = x^3 + b_j for j = 0..5 (one point is enough: the order lies in the
    Hasse window and no other of the six orders is a multiple of it),
  * and finally certify the first hit of class 0 (b = 7), which is the class of
    the curve the paper actually constructs.

Usage:  python -u 15_classify_hits.py
"""

import json
import os
import re
import sys
import time
from math import isqrt

import gmpy2

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "data")
RESULTS = os.path.join(os.path.dirname(HERE), "results")
os.makedirs(RESULTS, exist_ok=True)
LOGPATH = os.path.join(RESULTS, "run_minimality.txt")
OUTPATH = os.path.join(DATA, "paper_minimality.json")

MODELS = [7, 21, 63, 189, 567, 1701]
N = 2560
S = 320
C_TARGET = 132200237




def jac_double(P, p):
    X1, Y1, Z1 = P
    if Y1 == 0 or Z1 == 0:
        return (0, 1, 0)
    YY = Y1 * Y1 % p
    Sx = 4 * X1 * YY % p
    M = 3 * X1 * X1 % p
    X3 = (M * M - 2 * Sx) % p
    Y3 = (M * (Sx - X3) - 8 * YY * YY) % p
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
    v = (x * x % p * x + b) % p
    s = gmpy2.powmod(v, (p + 1) // 4, p)
    if (s * s - v) % p != 0:
        return None
    return (x, int(s))


def order_annihilates(p, b, order, tests=2):
    x = 1
    pts = []
    while len(pts) < tests and x < 20000:
        P = lift_x(x, p, b)
        if P is not None:
            pts.append(P)
        x += 1
    if not pts:
        return False
    for P in pts:
        R = jac_mul(order, (P[0], P[1], 1), p)
        if R[2] != 0:
            return False
    return True




def p_of(j):
    return (1 << N) - (1 << S) - (5 + 12 * j)


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


def sextic(p):
    A, B = cornacchia(p)
    traces = [2 * A, -2 * A, A + 3 * B, -(A + 3 * B), A - 3 * B, -(A - 3 * B)]
    orders = [p + 1 - t for t in traces]
    prime = [o for o in orders if gmpy2.powmod(2, o - 1, o) == 1 and gmpy2.is_prime(o) > 0]
    return A, B, traces, orders, prime


def log(*parts):
    msg = " ".join(str(p) for p in parts)
    sys.stdout.write(msg + "\n")
    sys.stdout.flush()


def parse_hits():
    js = []
    with open(LOGPATH, "r", encoding="utf-8") as fh:
        for line in fh:
            m = re.search(r"recipe hit: j = (\d+)\s+c = (\d+)", line)
            if m:
                js.append((int(m.group(1)), int(m.group(2))))
    seen = set()
    out = []
    for j, c in js:
        if j not in seen:
            seen.add(j)
            out.append((j, c))
    return out


def main():
    t0 = time.time()
    hits = parse_hits()
    log("parsed %d recipe hits from %s" % (len(hits), LOGPATH))

    records = []
    for idx, (j, c) in enumerate(hits, 1):
        p = p_of(j)
        A, B, traces, orders, prime = sextic(p)
        if len(prime) != 1:
            log("  j = %d has %d prime orders (skip)" % (j, len(prime)))
            continue
        o = prime[0]
        owner = None
        for cls, b in enumerate(MODELS):
            if order_annihilates(p, b, o, tests=2):
                owner = cls
                break
        rec = {
            "j": j, "c": c,
            "cornacchia_A": A, "cornacchia_B": B,
            "traces": traces, "orders": orders,
            "prime_order": o, "prime_order_bits": o.bit_length(),
            "prime_order_trace": p + 1 - o,
            "class_index": owner,
            "class_b": MODELS[owner] if owner is not None else None,
            "is_class0_b7": owner == 0,
        }
        records.append(rec)
        log("  [%2d/%2d] j = %-9d c = %-11d class j = %s (b = %s)  order bits %d"
            % (idx, len(hits), j, c, owner, rec["class_b"], o.bit_length()))

    class0 = [r for r in records if r["is_class0_b7"]]
    first_class0 = class0[0] if class0 else None
    first_recipe = records[0] if records else None

    log("first recipe hit      : j = %s c = %s class b = %s"
        % (first_recipe["j"], first_recipe["c"], first_recipe["class_b"]) if first_recipe else "none")
    if first_class0:
        log("FIRST b = 7 (class 0) : j = %d c = %d" % (first_class0["j"], first_class0["c"]))
    else:
        log("FIRST b = 7 (class 0) : NONE FOUND")

    cert = {}
    if first_class0:
        j = first_class0["j"]
        p = p_of(j)
        o = first_class0["prime_order"]
        tc = time.time()
        ok12 = order_annihilates(p, 7, o, tests=12)
        cert = {
            "j": j,
            "c": first_class0["c"],
            "p_is_probable_prime": bool(gmpy2.is_prime(p)),
            "order_is_probable_prime": bool(gmpy2.is_prime(o)),
            "order_annihilates_12_points_of_b7": bool(ok12),
            "certification_seconds": round(time.time() - tc, 1),
        }
        try:
            import sympy
            cert["sympy_isprime_p"] = bool(sympy.isprime(p))
            cert["sympy_isprime_order"] = bool(sympy.isprime(o))
        except Exception as exc:
            cert["sympy_error"] = str(exc)
        log("certification of the first b = 7 hit: %s" % cert)

    j_ours = (C_TARGET - 5) // 12
    p_ours = p_of(j_ours)
    A, B, traces, orders, prime = sextic(p_ours)
    ours = {
        "j": j_ours, "c": C_TARGET,
        "cornacchia_A": A, "cornacchia_B": B,
        "traces": traces, "orders": orders,
        "prime_orders": prime,
        "class_index": None, "class_b": None,
    }
    if prime:
        owner = None
        for cls, b in enumerate(MODELS):
            if order_annihilates(p_ours, b, prime[0], tests=4):
                owner = cls
                break
        ours["class_index"] = owner
        ours["class_b"] = MODELS[owner] if owner is not None else None
    log("the paper's own curve: c = %d -> class j = %s (b = %s)"
        % (C_TARGET, ours["class_index"], ours["class_b"]))

    out = {
        "meta": {
            "generated_by": "code/15_classify_hits.py",
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            "python": sys.version.split()[0],
            "gmpy2": gmpy2.version(),
        },
        "target": {"n": N, "s": S, "c": C_TARGET,
                   "j_max": (C_TARGET - 5) // 12,
                   "candidates": (C_TARGET - 5) // 12 + 1},
        "sweep": {"bound": 20000000, "survivors": 1102457,
                  "prime_p": 18646, "recipe_hits": 52, "ambiguous": 0,
                  "sweep_seconds": 1518.8},
        "counts": {"classified": len(records),
                   "class0_hits": len(class0),
                   "distinct_classes": len({r["class_index"] for r in records})},
        "first_recipe_hit": first_recipe,
        "first_class0_hit": first_class0,
        "paper_curve": ours,
        "hits": records,
        "certification": cert,
    }
    with open(OUTPATH, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1)
    log("wrote %s  (%.1fs total)" % (OUTPATH, time.time() - t0))
    return 0


if __name__ == "__main__":
    sys.exit(main())
