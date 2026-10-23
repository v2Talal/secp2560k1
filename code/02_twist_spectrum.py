"""
02_twist_spectrum.py
  Clean, certified 12-curve twist spectrum of secp256k1 + the 3-isogeny class.
"""
import json, os, time
from sclab import *

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(os.path.dirname(HERE), "results")

def say(*a):
    print(*a, flush=True)

t0 = time.time()
out = {}

x, y = cm_repr(P)
cands = candidate_traces(P)
lo, hi = hasse_window(P)

def order_of(bb):
    for tr in cands:
        cand = P + 1 - tr
        if not (lo <= cand <= hi):
            continue
        if certify_order(bb, cand, tests=1)[0] and certify_order(bb, cand, tests=2)[0]:
            return tr, cand
    return None, None

say("=" * 84)
say("THE 12 Fp-DISTINCT j=0 CURVES OVER F_p   (b = 7*3^j)")
say("=" * 84)
say(f"{'j':>2} {'b = 7*3^j (hex)':<22} {'trace':>42}  order-factorstructure")
rows = []
for j in range(12):
    b = (B * pow(3, j, P)) % P
    tr, o = order_of(b)
    fac, cof, cofprime = small_factor_report(o)
    say(f"{j:>2} {hex(b)[:20]:<22} {tr:>42}  {fac if fac else 'no factor <= 1e6'}  cofactor prime? {cofprime}")
    rows.append({"j": j, "b_hex": hex(b), "trace": tr, "order": o,
                 "small_factors": {str(k): v for k, v in fac.items()},
                 "cofactor": cof, "cofactor_prime": bool(cofprime)})
    out.setdefault("curves", []).append(rows[-1])

say("")
say("=" * 84)
say("ISOMORPHISM CLASSES  (b'/b must be a 6th power in Fp*)")
say("=" * 84)
exp6 = (P - 1) // 6
classes = {}
for r in rows:
    b = int(r["b_hex"], 16)
    key = None
    for r2 in rows:
        b2 = int(r2["b_hex"], 16)
        if pow(b * inv(b2, P) % P, exp6, P) == 1:
            key = r2["j"]
            break
    classes.setdefault(key, []).append(r["j"])
for k, v in classes.items():
    say(f"  class rep j={k}: members j = {v}   (order {rows[k]['order']})")

say("")
say("=" * 84)
say("QUADRATIC TWIST PAIRING   (b -> b*d^3, d = 3 nonsquare)")
say("=" * 84)
for r in rows:
    b = int(r["b_hex"], 16)
    bq, nu = quadratic_twist_b(b)
    j2 = next((r2["j"] for r2 in rows if int(r2["b_hex"], 16) == bq), None)
    tr2 = rows[j2]["trace"] if j2 is not None else None
    ok = (tr2 is not None and tr2 == -r["trace"])
    say(f"  j={r['j']:>2}  twist -> j={j2}   trace {r['trace']} -> {tr2}   opposite? {ok}")
    r["twist_j"] = j2
out["classes"] = {str(k): v for k, v in classes.items()}
out["curves"] = rows

say("")
say("=" * 84)
say("THE 3-ISOGENY PARTNER  (E_b -> E_{-27b})")
say("=" * 84)
b27 = (-27 * B) % P
tr27, o27 = order_of(b27)
say(f"  E(-27b): b' = {hex(b27)[:20]}  trace {tr27}  order {o27}")
say(f"  == n ? {o27 == N}")
j27 = next((r["j"] for r in rows if int(r["b_hex"], 16) == b27), None)
say(f"  b' equals 7*3^j for j = {j27}  -> same Fp-isomorphism class as that curve")
say(f"  is -27 a 6th power mod p ? {pow((-27) % P, exp6, P) == 1}")
b4 = (4 * B) % P
tr4, o4 = order_of(b4)
say(f"  other partner E(4b): b' = {hex(b4)[:20]}  trace {tr4}  order {o4}   == n? {o4 == N}")
j4 = next((r["j"] for r in rows if int(r["b_hex"], 16) == b4), None)
say(f"  b' equals 7*3^j for j = {j4}")
out["isogeny"] = {"b_minus27": hex(b27), "order": o27, "trace": tr27, "j": j27,
                  "b_4": hex(b4), "order4": o4, "trace4": tr4, "j4": j4,
                  "minus27_is_6th_power": pow((-27) % P, exp6, P) == 1}

G = find_point(B, A, P, seed=1)
img = isogeny3_point(G, B, P)
say(f"  isogeny3(G) on E(-27b)? {on_curve(img, 0, b27)}")
say(f"  3*G + phi(3G)?  checking n*phi(G) == O : {mul(N, img) is None}")
out["isogeny"]["map_verified"] = bool(on_curve(img, 0, b27))
out["isogeny"]["n_times_image_is_O"] = bool(mul(N, img) is None)

say(f"\ndone in {time.time()-t0:.1f}s")
os.makedirs(RESULTS, exist_ok=True)
outpath = os.path.join(RESULTS, "02_out.json")
with open(outpath, "w") as f:
    json.dump(out, f, indent=1, default=str)
say(f"written -> {outpath}")
