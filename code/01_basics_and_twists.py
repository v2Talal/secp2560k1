"""
01_basics_and_twists.py
  (a) certify the basic invariants of secp256k1
  (b) certify the FULL twist spectrum: all six sextic twists + their quadratic
      twists, with certified group orders and (time-bounded) factorisations.

Everything is time-bounded and prints progress with flush=True so the run
never appears to hang.
"""
import json, os, time
from sclab import *

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(os.path.dirname(HERE), "results")

def say(*a):
    print(*a, flush=True)

t_start = time.time()
out = {}

say("=" * 78)
say("A. FIELD / CURVE INVARIANTS")
say("=" * 78)
say(f"p  = 2^256 - 2^32 - 977   = {P}")
say(f"p  bits {bits(P)} | p mod 3 = {P%3} | mod 4 = {P%4} | mod 8 = {P%8}")
say(f"9 | p-1 ? {(P-1)%9==0}      36 | p-1 ? {(P-1)%36==0}   3 | p-1 ? {(P-1)%3==0}")
say(f"n  bits {bits(N)} | n mod 3 = {N%3} | mod 4 = {N%4} | mod 8 = {N%8} | prime? {isprime(N)}")
say(f"p+1 mod n = {(P+1)%N}")
out["p_mod"] = {"mod3": P%3, "mod4": P%4, "mod8": P%8, "p1_div9": (P-1)%9==0,
                "p1_div36": (P-1)%36==0}
out["n_mod"] = {"mod3": N%3, "mod4": N%4, "mod8": N%8, "prime": bool(isprime(N))}
out["p1_mod_n"] = (P+1) % N

ok, why = certify_order(B, N)
say(f"\ncertify #E(Fp) == n ?  {ok}  ({why})   [{time.time()-t_start:.1f}s]")
assert ok
t = trace_of(N)
say(f"Frobenius trace t = p+1-n = {t}   (bits {bits(abs(t))})")
out["trace"] = t

x, y = cm_repr(P)
say(f"\nCornacchia(d=3): p = x^2 + 3y^2")
say(f"  x = {x}  (bits {bits(abs(x))})")
say(f"  y = {y}  (bits {bits(abs(y))})")
say(f"  check: {x*x + 3*y*y == P}")
assert x*x + 3*y*y == P
cands = candidate_traces(P)
say(f"  six CM trace candidates: {cands}")
say(f"  known trace t in set: {t in cands}")
assert t in cands
out["cm_x"], out["cm_y"], out["cm_candidates"] = x, y, cands

Lc = 2*x
rem = 4*P - Lc*Lc
Mc = isqrt(rem//27) if rem % 27 == 0 else None
if Mc is not None and Mc*Mc*27 == rem and (2*y) % 3 == 0:
    say(f"  classical: 4p = L^2 + 27M^2 with L = {Lc}, M = {(2*y)//3}")
    out["cm_L"], out["cm_M"] = Lc, (2*y)//3
else:
    say("  (classical 27M^2 form needs 3 | 2y; not applicable here)")
    out["cm_L"], out["cm_M"] = None, None

z1, z2 = cube_roots_of_unity(P)
say(f"\ncube roots of unity:  zeta1 = {z1}\n                     zeta2 = {z2}")
out["zeta1"], out["zeta2"] = z1, z2
G = find_point(B, A, P, seed=1)
say(f"base point G on curve: {on_curve(G)}")
lam1 = next((l for l in n_cube_roots(N) if mul(l, G) == phi_point(G, z1)), None)
lam2 = next((l for l in n_cube_roots(N) if mul(l, G) == phi_point(G, z2)), None)
say(f"lambda (zeta1) = {lam1}\nlambda (zeta2) = {lam2}")
assert lam1 is not None and lam2 is not None
assert (lam1*lam1 + lam1 + 1) % N == 0 and lam2 == lam1*lam1 % N
out["lambda_zeta1"], out["lambda_zeta2"] = lam1, lam2
say(f"lambda hex = {hex(lam1)}")

say("")
say("=" * 78)
say("B. TWIST SPECTRUM (the 12 curves Fp-isomorphic to j = 0)")
say("=" * 78)
c0, reps, exp = sixth_power_cosets(P)
say(f"generator c0 of F_p*/(F_p*)^6 : {c0}")

def order_of_curve(bb):
    """Certified order of y^2 = x^3 + bb using the CM candidate traces."""
    lo, hi = hasse_window(P)
    for tr in cands:
        cand = P + 1 - tr
        if not (lo <= cand <= hi):
            continue
        ok, why = certify_order(bb, cand, tests=1)
        if ok:
            ok2, _ = certify_order(bb, cand, tests=2)
            if ok2:
                return tr, cand
    return None, None

twists = []
for idx, c in enumerate(reps):
    bb = (B * c) % P
    tr, ordr = order_of_curve(bb)
    fac, cof, cofprime = small_factor_report(ordr)
    nu_b, nu = quadratic_twist_b(bb)
    tr2, ordr2 = order_of_curve(nu_b)
    fac2, cof2, cof2prime = small_factor_report(ordr2)
    say(f"\nsextic twist #{idx}: b' = {hex(bb)[:22]}...")
    say(f"    trace = {tr}")
    say(f"    order = {ordr}")
    say(f"    small factors <= 1e6 : {fac if fac else 'none'}")
    say(f"    cofactor {cof} is prime? {cofprime}")
    say(f"  quadratic twist (b'*{nu} mod p):")
    say(f"    trace = {tr2}")
    say(f"    order = {ordr2}")
    say(f"    small factors <= 1e6 : {fac2 if fac2 else 'none'}")
    say(f"    cofactor {cof2} is prime? {cof2prime}")
    twists.append({"idx": idx, "b_sextic_hex": hex(bb), "trace_sextic": tr,
                    "order_sextic": ordr, "small_factors_sextic": {str(k): v for k, v in fac.items()},
                    "cofactor_sextic": cof, "cofactor_sextic_prime": bool(cofprime),
                    "b_quad_hex": hex(nu_b), "trace_quad": tr2, "order_quad": ordr2,
                    "small_factors_quad": {str(k): v for k, v in fac2.items()},
                    "cofactor_quad": cof2, "cofactor_quad_prime": bool(cof2prime)})
    say(f"    [{time.time()-t_start:.1f}s]")
out["twists"] = twists

say("\nSUMMARY -- prime-order members of the 12-curve j=0 family:")
for t_ in twists:
    if t_["cofactor_sextic_prime"] and not t_["small_factors_sextic"]:
        say(f"  sextic #{t_['idx']}: order {t_['order_sextic']}  PRIME")
    if t_["cofactor_quad_prime"] and not t_["small_factors_quad"]:
        say(f"  quad of sextic #{t_['idx']}: order {t_['order_quad']}  PRIME")

os.makedirs(RESULTS, exist_ok=True)
outpath = os.path.join(RESULTS, "01_out.json")
with open(outpath, "w") as f:
    json.dump(out, f, indent=1, default=str)
say(f"\ndone in {time.time()-t_start:.1f}s -> {outpath}")
