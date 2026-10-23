"""Why the small-factor structure of secp2560k1 is what it is.

Q1  Does the twist cofactor really have no prime factor below 2e6 (the
    originating paper's claim), or did it only look that way?
Q2  Every one of the six sextic-class orders is a norm in the Eisenstein
    ring Z[omega].  Which primes can divide such a norm, and does that
    predict the observed small factors?
Q3  How surprising is "composite but with no factor below 2e6"?
    (exact Mertens products, and the arithmetic-progression-restricted
    version that the norm structure actually permits)

Integer-only; every loop is capped.
"""
import json, os
from math import gcd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "data")
OUT = {}
def log(*a, **k):
    print(*a, flush=True)
D  = json.load(open(os.path.join(DATA, "paper_data.json")))
X2 = json.load(open(os.path.join(DATA, "paper_extra2.json")))

P  = D["params"]; TW = D["twist"]
p  = int(P["p_dec"]);  n  = int(P["n_dec"])
nt = int(TW["nt"])
cx = int(X2["cornacchia_stats"]["x_hex"], 16)
cy = int(X2["cornacchia_stats"]["y_hex"], 16)

log("p bits", p.bit_length(), " n bits", n.bit_length())
log("cornacchia x bits", cx.bit_length(), " y bits", cy.bit_length())

def small_factors_int(n, bound=10**6):
    f = {}
    d = 2
    while d * d <= n and d <= bound:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f

CANDS = [
    ("trace=+2x",     cx - 1,        cy,        1),
    ("trace=-2x",     cx + 1,        cy,        1),
    ("trace=+(x+3y)", 2*cx - 1,  2*cy - 1,     4),
    ("trace=-(x+3y)", 2*cx + 1,  2*cy + 1,     4),
    ("trace=+(x-3y)", 2*cx - 1,  2*cy + 1,     4),
    ("trace=-(x-3y)", 2*cx + 1,  2*cy - 1,     4),
]
models = X2["sextic_models"]
rows = []
for m in models:
    ordv = int(m["order"])
    hit = None
    for (name, A, B, d) in CANDS:
        if (A*A + 3*B*B)//d == ordv:
            hit = (name, A, B, d); break
    assert hit, "no Eisenstein representation for j=%s" % m["j"]
    name, A, B, d = hit
    g = gcd(A, B)
    gf = small_factors_int(g)
    sf = {int(k): v for k, v in (m.get("small_factors", {}) or {}).items()}
    chk = []
    for q in sorted(sf):
        if q == 3:
            chk.append({"q": q, "kind": "ramified", "divides_both": None})
        elif q % 3 == 1:
            chk.append({"q": q, "kind": "split", "divides_both": None})
        else:
            chk.append({"q": q, "kind": "inert", "divides_both": bool(A % q == 0 and B % q == 0)})
    rows.append(dict(
        j=m["j"], b=m["b_dec"], trace=m["trace"], order_bits=ordv.bit_length(),
        is_prime=bool(m.get("is_prime_model")),
        form=name, d=d, A_bits=A.bit_length(), B_bits=B.bit_length(),
        gcd_AB=g, gcd_AB_bits=g.bit_length(),
        gcd_AB_factors={str(k): v for k, v in gf.items()},
        inert_candidates=sorted(q for q in gf if q % 3 == 2),
        small_factors={str(k): v for k, v in sf.items()},
        factor_check=chk,
        A_mod_q={str(q): (A % q) for q in sorted(sf)},
        _A=A, _B=B,
    ))
    log("j=%s b=%s form=%s gcd(A,B)=%s (bits %d) small=%s" %
        (m["j"], m["b_dec"], name, g, g.bit_length(), sf), flush=True)

def sieve_primes(limit):
    s = bytearray([1]) * (limit + 1)
    s[0] = 0
    if limit >= 1: s[1] = 0
    for i in range(2, int(limit ** 0.5) + 1):
        if s[i]:
            s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(3, limit + 1) if s[i]]

log("sieving primes to 2e6 ...", flush=True)
PR = sieve_primes(2_000_000)
log("primes(3..2e6) =", len(PR), flush=True)

def divisors(M, primes):
    out = {}
    for q in primes:
        if M % q == 0:
            e = 0
            while M % q == 0:
                M //= q; e += 1
            out[q] = e
    return out

rem = nt
v3 = 0
while rem % 3 == 0:
    rem //= 3; v3 += 1
log("v3(nt) =", v3, " remainder bits =", rem.bit_length(), flush=True)
sf_tw = divisors(rem, PR)
log("small factors of twist cofactor below 2e6:", sf_tw, flush=True)
sf_nt = divisors(nt, PR)
log("small factors of nt below 2e6:", sf_nt, flush=True)

for r, m in zip(rows, models):
    sf = divisors(int(m["order"]), PR)
    r["sieve_2e6"] = {str(k): v for k, v in sf.items()}
    log("j=%s b=%s sieve<=2e6: %s" % (r["j"], r["b"], r["sieve_2e6"]), flush=True)

prod_all = 1.0; prod_a1 = 1.0; prod_a2 = 1.0
cnt_a1 = 0; cnt_a2 = 0
for q in PR:
    prod_all *= (1.0 - 1.0/q)
    if q % 3 == 1:
        prod_a1 *= (1.0 - 1.0/q); cnt_a1 += 1
    else:
        prod_a2 *= (1.0 - 1.0/q); cnt_a2 += 1
log("P(no factor<=2e6, any odd prime) =", prod_all, flush=True)
log("P(no factor<=2e6, q=1 mod 3)      =", prod_a1, flush=True)
log("P(no factor<=2e6, q=2 mod 3)      =", prod_a2, flush=True)
log("counts q=1mod3:", cnt_a1, " q=2mod3:", cnt_a2, flush=True)

rules = []
for r, m in zip(rows, models):
    ordv = int(m["order"])
    A = r["_A"]; B = r["_B"]
    three = (ordv % 3 == 0)
    two   = (ordv % 2 == 0)
    r3 = (A % 3 == 0)
    r2 = (A % 2 == B % 2) if r["d"] == 1 else False
    r["rule_3"] = {"three_divides_order": three, "three_divides_A": r3, "consistent": three == r3}
    r["rule_2"] = {"two_divides_order": two, "A_parity_eq_B": r2, "consistent": two == r2}
    r["rule_inert"] = ("no odd prime = 2 mod 3 can divide this order (gcd(A,B) = %d)" % r["gcd_AB"])
    rules.append({"j": r["j"], "b": r["b"], "form": r["form"], **{k: r[k] for k in
                 ("rule_3", "rule_2", "rule_inert", "gcd_AB", "gcd_AB_factors")}})
    log("j=%s b=%s %s | 3|N=%s 3|A=%s | 2|N=%s A=B(2)=%s | gcd=%d" %
        (r["j"], r["b"], r["form"], three, r3, two, r2, r["gcd_AB"]), flush=True)

LIMIT2 = 20_000_000_000
log("sieving primes to 2e7 ...", flush=True)
PR2 = sieve_primes(LIMIT2)
log("primes(3..2e7) =", len(PR2), flush=True)
deep_tw = divisors(rem, PR2)
log("twist cofactor, factors <= 2e7:", deep_tw, flush=True)
for r, m in zip(rows, models):
    sf = divisors(int(m["order"]), PR2)
    r["sieve_2e7"] = {str(k): v for k, v in sf.items()}
    log("j=%s sieve<=2e7: %s" % (r["j"], r["sieve_2e7"]), flush=True)
OUT["rules"] = rules
OUT["deep_sieve_2e7"] = {"limit": LIMIT2,
                         "twist_cofactor": {str(k): v for k, v in deep_tw.items()},
                         "classes": {str(r["j"]): r["sieve_2e7"] for r in rows}}

p_all2 = 1.0; p_a1_2 = 1.0
for q in PR2:
    p_all2 *= (1.0 - 1.0/q)
    if q % 3 == 1:
        p_a1_2 *= (1.0 - 1.0/q)
OUT["probabilities_2e7"] = {"at_2e7_any_odd_prime": p_all2, "at_2e7_q_1_mod_3": p_a1_2}
log("P(no factor<=2e7, any odd prime) =", p_all2, flush=True)
log("P(no factor<=2e7, q=1 mod 3)      =", p_a1_2, flush=True)

OUT["cornacchia"] = {"x_bits": cx.bit_length(), "y_bits": cy.bit_length()}
OUT["classes"] = rows
OUT["twist"] = {"v3": v3, "remainder_bits": rem.bit_length(),
                "cofactor_small_factors_below_2e6": {str(k): v for k, v in sf_tw.items()},
                "nt_small_factors_below_2e6": {str(k): v for k, v in sf_nt.items()}}
OUT["probabilities"] = {"any_odd_prime": prod_all, "q_1_mod_3": prod_a1, "q_2_mod_3": prod_a2,
                        "primes_1_mod_3": cnt_a1, "primes_2_mod_3": cnt_a2}
OUT.setdefault("probabilities", {}).update(OUT.pop("probabilities_2e7", {}))
json.dump(OUT, open(os.path.join(DATA, "paper_smallfac.json"), "w"), indent=1)
log("wrote paper_smallfac.json")
