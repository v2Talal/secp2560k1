"""Independent verification of the user's own curve secp2560k1 (paper_secp2560k1.pdf).
Ground truth oracle: point-order test + exact Hasse window (isqrt, never float sqrt).
NO unbounded factoring. NO float sqrt. All prints flushed.
"""
import json, os, sys, time, hashlib
from math import isqrt
from sympy import isprime

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "data")
RESULTS = os.path.join(os.path.dirname(HERE), "results")

T0 = time.time()
def log(*a):
    print(f"[{time.time()-T0:7.1f}s]", *a, flush=True)

J = json.load(open(os.path.join(DATA, "secp2560k1.json"), encoding="utf-8"))
P = int(J["p_hex"], 16)
N = int(J["n_hex"], 16)
A = int(J["a_hex"], 16)
B = int(J["b_hex"], 16)
Gx = int(J["Gx_hex"], 16); Gy = int(J["Gy_hex"], 16)
G2x = int(J["G_even"]["Gx_hex"], 16); G2y = int(J["G_even"]["Gy_hex"], 16)
omega = int(J["glv"]["omega_hex"], 16)
lam = int(J["glv"]["lambda_hex"], 16)
a1 = J["glv_basis"]["a1"]; b1 = J["glv_basis"]["b1"]
a2 = J["glv_basis"]["a2"]; b2 = J["glv_basis"]["b2"]
R = {}

def jdbl(Pt, p, b):
    X, Y, Z = Pt
    if Y == 0 or Z == 0: return (0, 0, 0)
    A_ = (Y*Y) % p
    B_ = (4*X*A_) % p
    C_ = (8*A_*A_) % p
    D_ = (3*X*X) % p
    X3 = (D_*D_ - 2*B_) % p
    Y3 = (D_*(B_-X3) - C_) % p
    Z3 = (2*Y*Z) % p
    return (X3, Y3, Z3)
def jadd(Pt, Q, p, b):
    X1,Y1,Z1 = Pt; X2,Y2,Z2 = Q
    if Z1 == 0: return Q
    if Z2 == 0: return Pt
    Z1Z1 = (Z1*Z1)%p; Z2Z2 = (Z2*Z2)%p
    U1 = (X1*Z2Z2)%p; U2 = (X2*Z1Z1)%p
    S1 = (Y1*Z2*Z2Z2)%p; S2 = (Y2*Z1*Z1Z1)%p
    if U1 == U2:
        if S1 != S2: return (0,0,0)
        return jdbl(Pt, p, b)
    H = (U2-U1)%p; I = (2*H)%p; I = (I*I)%p
    Jx = (H*I)%p; V = (U1*I)%p
    Y3 = (S2-S1)%p; Y3 = (2*Y3*J)%p
    X3 = (Y3*Y3 - Jx*Jx - 2*V)%p
    Z3 = (((Z1+Z2)*(Z1+Z2) - Z1Z1 - Z2Z2)*H)%p
    return (X3, Y3, Z3)
def to_affine(Pt, p):
    X,Y,Z = Pt
    if Z == 0: return None
    zi = pow(Z, -1, p); zi2 = (zi*zi)%p
    return ((X*zi2)%p, (Y*zi2*zi)%p)
def mul(k, aff, p, b):
    Pt = (aff[0], aff[1], 1)
    if k == 0: return None
    if k < 0: k = (-k) % N
    R_ = (0,0,0)
    for bit in bin(k)[2:]:
        R_ = jdbl(R_, p, b)
        if bit == '1': R_ = jadd(R_, Pt, p, b)
    return to_affine(R_, p)
def on_curve(aff, p, b):
    if aff is None: return False
    x, y = aff
    return (y*y - (x*x*x + b)) % p == 0
def legendre(v, p): return pow(v, (p-1)//2, p)

log("=== 1. FIELD / CURVE ===")
c = 132200237
R["p_form_ok"] = (P == 2**2560 - 2**320 - c)
R["p_bits"] = P.bit_length(); R["n_bits"] = N.bit_length()
R["p_mod3"] = P % 3; R["p_mod4"] = P % 4; R["p_mod12"] = P % 12
R["p_mod8"] = P % 8
R["n_mod3"] = N % 3
R["p_primality_probable"] = bool(isprime(P))
R["n_primality_probable"] = bool(isprime(N))
R["G_on_curve"] = on_curve((Gx, Gy), P, B)
R["G2_on_curve"] = on_curve((G2x, G2y), P, B)
R["G2x_eq_Gx"] = (G2x == Gx)
R["G2_is_negG"] = (G2y == (P - Gy) % P)
R["Gy_parity"] = Gy % 2; R["G2y_parity"] = G2y % 2
log("p_form_ok", R["p_form_ok"], "p_bits", R["p_bits"], "n_bits", R["n_bits"])
log("p%3", R["p_mod3"], "p%4", R["p_mod4"], "p%12", R["p_mod12"], "n%3", R["n_mod3"])
log("p prime?", R["p_primality_probable"], "n prime?", R["n_primality_probable"])
log("G on curve?", R["G_on_curve"], "G2 on curve?", R["G2_on_curve"])
log("G2x==Gx", R["G2x_eq_Gx"], "G2==-G ?", R["G2_is_negG"], "Gy parity", R["Gy_parity"], "G2y parity", R["G2y_parity"])
log("Gy tail hex", hex(Gy)[-8:], " G2y tail hex", hex(G2y)[-8:])

def sqrt_mod(a, p):
    a %= p
    if a == 0: return 0
    if legendre(a, p) != 1: return None
    if p % 4 == 3: return pow(a, (p+1)//4, p)
    return None
y15 = sqrt_mod((2**3 + B) % P, P)
R["x0_2_is_Q"] = (y15 is not None)
R["2P0_x"] = 2 if y15 is not None else None
log("x0=2 on curve (x^3+7 QR)?", R["x0_2_is_Q"], "  -> G.x should be 2 if G=2*P0")
log("G.x == 2 ?", Gx == 2, " G.x hex head", hex(Gx)[:14], " G.x hex tail", hex(Gx)[-10:])

log("=== 2. ORDER / HASSE ===")
t = P + 1 - N
R["trace_t"] = t
R["hasse_lo"] = P + 1 - 2*isqrt(P); R["hasse_hi"] = P + 1 + 2*isqrt(P)
R["n_in_hasse"] = (R["hasse_lo"] <= N <= R["hasse_hi"])
og = mul(N, (Gx, Gy), P, B)
R["N_G_is_O"] = (og is None)
R["order_exactly_N"] = R["N_G_is_O"] and R["n_primality_probable"] and R["n_in_hasse"]
og2 = mul(N, (G2x, G2y), P, B)
R["N_G2_is_O"] = (og2 is None)
log("t =", t, "bits", t.bit_length())
log("N in Hasse window?", R["n_in_hasse"], " N*G==O ?", R["N_G_is_O"], " N*G2==O ?", R["N_G2_is_O"])
log("=> #E(Fp) = N, h=1 :", R["order_exactly_N"])

log("=== 3. CM / SEXTIC SPECTRUM ===")
s = pow((-3) % P, (P+1)//4, P)
assert (s*s) % P == (P-3) % P, "sqrt(-3) failed"
r0, r1 = P, s
while r1*r1 > P:
    r0, r1 = r1, r0 % r1
Acm = r1; rem = P - Acm*Acm
assert rem % 3 == 0
Bcm = isqrt(rem//3)
assert Acm*Acm + 3*Bcm*Bcm == P
R["cornacchia_A"] = Acm; R["cornacchia_B"] = Bcm
traces = sorted({2*Acm, -2*Acm, Acm+3*Bcm, -(Acm+3*Bcm), Acm-3*Bcm, -(Acm-3*Bcm)})
R["six_traces"] = traces
R["t_in_six_traces"] = (t in traces)
log("Cornacchia A bits", Acm.bit_length(), "B bits", Bcm.bit_length())
log("six CM traces:", [hex(x) for x in traces])
log("t in six traces (sextic technique valid)?", R["t_in_six_traces"])
cand = []
for tr in traces:
    o = P + 1 - tr
    cand.append({"trace": tr, "order_bits": o.bit_length(), "prime": bool(isprime(o))})
R["six_orders"] = cand
for c_ in cand:
    log("  trace", hex(c_["trace"])[:12], "order bits", c_["order_bits"], "prime", c_["prime"])
R["n_is_prime_candidate"] = any(c_["prime"] for c_ in cand)
R["paper_table_bitlengths_match"] = sorted(c_["order_bits"] for c_ in cand) == [2560,2560,2560,2561,2561,2561]
log("paper Table-3 bitlength multiset matches?", R["paper_table_bitlengths_match"])

log("=== 4. TWIST ===")
nt = P + 1 + t
R["nt_bits"] = nt.bit_length()
R["nt_equals_p1_plus_t"] = True
R["nt_mod3"] = nt % 3
m = nt; e3 = 0
while m % 3 == 0: m //= 3; e3 += 1
R["twist_3_adic_val"] = e3
R["twist_rest_bits"] = m.bit_length()
R["twist_rest_is_prime"] = bool(isprime(m))
log("nt bits", R["nt_bits"], " nt%3", R["nt_mod3"], " 3-adic val", e3, " rest bits", R["twist_rest_bits"])
log(">>> twist rest (nt/3) PRIME?", R["twist_rest_is_prime"], " <-- paper claims 'composite'")
R["gcd_n_nt"] = __import__("math").gcd(N, nt)
log("gcd(N, nt) =", R["gcd_n_nt"], "(1 => no additive/multiplicative transfer via twist)")
def small_factors(n, bound=2_000_000):
    facs = {}
    for pr in sieve_primes(bound):
        if pr*pr > n: break
        c2 = 0
        while n % pr == 0: n //= pr; c2 += 1
        if c2: facs[pr] = c2
    return facs, n
_sieve = {}
def sieve_primes(bound):
    if bound in _sieve: return _sieve[bound]
    import sympy
    lst = list(sympy.primerange(2, bound))
    _sieve[bound] = lst
    return lst
facs, rest = small_factors(m)
R["twist_rest_small_factors"] = facs
R["twist_rest_after_trial"] = rest
log("twist rest factors < 2e6:", facs, " remaining composite?", not isprime(rest) if rest > 1 else False)

d = 2
while legendre(d, P) == 1: d += 1
R["twist_d"] = d
bt = (7 * pow(d, 3, P)) % P
R["twist_b"] = bt
def find_point(b, p, seed=12345):
    x = seed
    for _ in range(4000):
        rhs = (x*x*x + b) % p
        y = sqrt_mod(rhs, p)
        if y is not None: return (x, y)
        x += 1
    return None
Qtw = find_point(bt, P)
R["twist_point_found"] = Qtw is not None
R["twist_nt_Q_is_O"] = (mul(nt, Qtw, P, bt) is None) if Qtw else None
R["twist_nt2_Q_not_O"] = (mul(nt*2, Qtw, P, bt) is None) if Qtw else None
log("twist b' = 7*d^3 ; nt*Q == O ?", R["twist_nt_Q_is_O"], " (2nt)*Q == O ?", R["twist_nt2_Q_not_O"])
def cube_root_if_exists(a, p):
    a %= p
    if a == 0: return 0
    k = pow(3, -1, (p-1)//3)
    r = pow(a, k, p)
    if (r*r*r) % p == a: return r
    return None
r3 = cube_root_if_exists((-4*bt) % P, P)
R["twist_rational_3torsion_x"] = r3
if r3 is not None:
    y3 = sqrt_mod((r3*r3*r3 + bt) % P, P)
    R["twist_3torsion"] = (r3, y3) if y3 else None
    if y3:
        P3 = (r3, y3)
        R["twist_3P_is_O"] = (mul(3, P3, P, bt) is None)
        log("twist rational 3-torsion found, 3P==O ?", R["twist_3P_is_O"])
R["rational_3torsion_on_E"] = (mul(3, (Gx, Gy), P, B) is not None and 3 % N != 0)
log("E has rational 3-torsion?", False, "(3 does not divide N:", N % 3, ")")

log("=== 5. GLV ===")
R["omega_cubed_1"] = (pow(omega, 3, P) == 1)
R["omega_ne_1"] = (omega != 1)
R["lambda_poly"] = ((lam*lam + lam + 1) % N == 0)
lamG = mul(lam, (Gx, Gy), P, B)
R["phi_G_ok"] = (lamG is not None and lamG[0] == (omega*Gx) % P and lamG[1] == Gy)
log("omega^3==1", R["omega_cubed_1"], " omega!=1", R["omega_ne_1"])
log("lambda^2+lambda+1==0 mod N", R["lambda_poly"], " lambda*G == (omega*Gx, Gy)?", R["phi_G_ok"])
R["a1_rel"] = ((a1 + b1*lam) % N == 0)
R["a2_rel"] = ((a2 + b2*lam) % N == 0)
det = a1*b2 - a2*b1
R["basis_det"] = det; R["basis_det_ok"] = (abs(det) == N)
log("a1+b1*lam==0", R["a1_rel"], " a2+b2*lam==0", R["a2_rel"], " det", det, " == +-N ?", R["basis_det_ok"])
log("basis vector bit sizes: a1", a1.bit_length(), "b1", abs(b1).bit_length(), "a2", a2.bit_length(), "b2", b2.bit_length())

def decompose_paper_float(k):
    return (round(b2*k/det), round(-b1*k/det))
def round_div(num, den):
    if den < 0: num, den = -num, -den
    q, r = divmod(num, den)
    if 2*r > den or (2*r == den and q % 2 == 1): q += 1
    return q
def decompose_exact(k):
    c1 = round_div(b2*k, det); c2 = round_div(-b1*k, det)
    k1 = k - c1*a1 - c2*a2
    k2 = -(c1*b1 + c2*b2)
    return k1, k2
float_bug = None
try:
    decompose_paper_float(12345)
except Exception as e:
    float_bug = f"{type(e).__name__}: {e}"
R["paper_float_decompose_error"] = float_bug
log("paper's printed decompose() on a 2560-bit scalar ->", float_bug)
import random
random.seed(7)
ok = True; bits1 = []; bits2 = []
for _ in range(20):
    k = random.randrange(1, N)
    k1, k2 = decompose_exact(k)
    if (k1 + k2*lam) % N != k: ok = False
    bits1.append(abs(k1).bit_length()); bits2.append(abs(k2).bit_length())
R["decompose_exact_ok"] = ok
R["decompose_bits"] = [max(bits1), max(bits2)]
log("exact-integer decomposition identity holds on 20 random scalars:", ok, " max bits", R["decompose_bits"])

log("=== 6. EMBEDDING DEGREE ===")
v = P % N
found = None
lim = 300000
for k in range(1, lim+1):
    v = (v * (P % N)) % N
    if v == 1:
        found = k; break
R["ord_n_p_lower"] = (found if found else lim+1)
R["ord_n_p_found"] = found
log("ord_n(p) :", ("NOT found up to %d" % lim) if found is None else found, flush=True)
if found is not None:
    log(">>> paper's 'embedding > 20000' claim:", "WRONG (exact value found)" if found <= 20000 else "consistent")
    R["embedding_claim_ok"] = found > 20000
else:
    R["embedding_claim_ok"] = True
log("done, elapsed", round(time.time()-T0,1), "s", flush=True)
os.makedirs(RESULTS, exist_ok=True)
json.dump(R, open(os.path.join(RESULTS, "06_out.json"),"w"), indent=1, default=str)
log("wrote 06_out.json")
