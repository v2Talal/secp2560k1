"""research/12_recipe_reproduces_secp256k1.py

Goal: test the user's "minimality recipe" at 256 bits and check whether it
REPRODUCES secp256k1 exactly (c = 977), which would be strong evidence that
the recipe is general and not curve-specific.

Recipe (from the user's paper/JSON):
    p = 2^256 - 2^32 - c,  c = 5, 17, 29, ... (c == 5 mod 12, stepping 12)
    keep c iff:
      (a) p is prime
      (b) exactly ONE of the six sextic-class orders p+1-t is prime
          (the "sextic twist order technique"), and
      (c) that prime order is the group order of E: y^2 = x^3 + 7

Bounded work only: 82 candidates, 256-bit primality tests (fast), Cornacchia
(uses exact integer sqrt).  No unbounded factoring anywhere.
"""
import sys
import sympy

from sclab import cm_repr, candidate_traces, P, N

C_TARGET = 977


def main():
    print("=" * 74)
    print("Does the minimality recipe (c = 5 mod 12, step 12) reproduce")
    print("secp256k1 (c = 977) at 256 bits?")
    print("=" * 74)
    first_hit = None
    hits = []
    tested = 0
    for c in range(5, C_TARGET + 1, 12):
        tested += 1
        p = 2**256 - 2**32 - c
        if not sympy.isprime(p):
            continue
        try:
            x, y = cm_repr(p)
        except Exception as e:
            print(f"  c={c:4d}  p prime but Cornacchia failed: {e}")
            continue
        traces = candidate_traces(p)
        orders = [p + 1 - t for t in traces]
        prime_orders = [o for o in orders if sympy.isprime(o)]
        if len(prime_orders) == 1:
            o = prime_orders[0]
            t = p + 1 - o
            hits.append((c, p, o, t))
            print(f"  c={c:4d}  p prime  |  exactly one prime sextic order: {o}")
            print(f"            == secp256k1 n ? {o == N}   trace t = {t}")
            if first_hit is None:
                first_hit = c
    print("-" * 74)
    print(f"candidates tested (c = 5 mod 12, c <= 977): {tested}")
    print(f"first c giving a prime p with a unique prime-order sextic member: "
          f"{first_hit}")
    print(f"secp256k1's c: {C_TARGET}   ->  REPRODUCES EXACTLY: "
          f"{first_hit == C_TARGET}")
    if hits:
        print(f"all such c <= 977: {[h[0] for h in hits]}")
    print()
    if first_hit == C_TARGET:
        print("VERDICT: the recipe is GENERAL -- at 256 bits it lands exactly on")
        print("         secp256k1's c = 977.  This is a strong, cheap")
        print("         reproducibility argument to put in the paper.")
    else:
        print("VERDICT: the recipe does NOT land on 977 with criteria (a)+(b);")
        print("         the paper's minimality rule must include further")
        print("         conditions (twist security, embedding degree, GLV, ...).")
    sys.stdout.flush()


if __name__ == "__main__":
    main()
