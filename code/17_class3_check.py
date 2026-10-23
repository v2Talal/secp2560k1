"""Why 13 of the 52 recipe hits have no owner among the six models 7*3^j.

The six models 7*3^j (j = 0..5) represent the six F_p-isomorphism classes of the
j = 0 family only if the class of 3 generates the sextic quotient F_p*/(F_p*)^6.
In Z/6 the class of 3 has order 6/gcd(e,6) where 3 = g^e for a generator g of F_p*.
So it generates iff e is odd (3 is a non-square) and e is not divisible by 3
(3 is not a cube).  This script tests that for all 52 hits and cross-checks the
result against the annihilation test in 15_classify_hits.py.
"""
import json, os, sys, gmpy2

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "data")
d = json.load(open(os.path.join(DATA, "paper_minimality.json")))
MODELS = [7, 21, 63, 189, 567, 1701]
N, S = 2560, 320


def p_of(j):
    return (1 << N) - (1 << S) - (5 + 12 * j)


rows = []
for h in d["hits"]:
    j = h["j"]
    p = p_of(j)
    nonsq = gmpy2.powmod(3, (p - 1) // 2, p) == p - 1
    notcube = gmpy2.powmod(3, (p - 1) // 3, p) != 1
    rows.append({
        "j": j, "c": h["c"],
        "class_index": h["class_index"],
        "is_class0_b7": h["is_class0_b7"],
        "three_is_non_square": bool(nonsq),
        "three_is_not_a_cube": bool(notcube),
        "class3_generates": bool(nonsq and notcube),
    })

agree = 0
for r in rows:
    if (r["class3_generates"]) == (r["class_index"] is not None):
        agree += 1
    else:
        print("MISMATCH", r, flush=True)

from collections import Counter
cnt = Counter((r["class3_generates"], r["class_index"] is not None) for r in rows)
print("cross-tab (generates, has_owner):", dict(cnt), flush=True)
print("agreement rows: %d / %d" % (agree, len(rows)), flush=True)

for r in rows:
    if r["class_index"] is None:
        print("  no owner: j=%-9d c=%-11d nonsq=%s notcube=%s" %
              (r["j"], r["c"], r["three_is_non_square"], r["three_is_not_a_cube"]), flush=True)

print()
print("class-0 rows all generate:", all(r["class3_generates"] for r in rows if r["is_class0_b7"]), flush=True)

json.dump({"rows": rows}, open(os.path.join(DATA, "paper_minimality_class3.json"), "w"), indent=1)
print("saved paper_minimality_class3.json", flush=True)
