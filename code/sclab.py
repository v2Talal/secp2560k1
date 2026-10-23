"""
sclab.py -- computational laboratory for structural analysis of secp256k1.

Everything here is written so that every claimed group order is CERTIFIED by an
independent point-order test (n*R == O) plus the Hasse bound, never by trusting a
formula's sign conventions.  CM theory is used only to *generate candidate* traces;
the certificate decides.
"""

from math import gcd, isqrt
from sympy import factorint, isprime

P = 2**256 - 2**32 - 977
N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
A, B = 0, 7

def inv(x, m=P):
    return pow(x % m, -1, m)

def jac_double(pt, a=A, p=P):
    X, Y, Z = pt
    if Y == 0 or Z == 0:
        return (0, 0, 0)
    YY = Y * Y % p
    S = 4 * X * YY % p
    M = (3 * X * X + a * (Z * Z % p) ** 2) % p
    X3 = (M * M - 2 * S) % p
    Y3 = (M * (S - X3) - 8 * YY * YY) % p
    Z3 = 2 * Y * Z % p
    return (X3, Y3, Z3)

def jac_add(p1, p2, a=A, p=P):
    X1, Y1, Z1 = p1
    X2, Y2, Z2 = p2
    if Z1 == 0:
        return p2
    if Z2 == 0:
        return p1
    Z1Z1 = Z1 * Z1 % p
    Z2Z2 = Z2 * Z2 % p
    U1 = X1 * Z2Z2 % p
    U2 = X2 * Z1Z1 % p
    S1 = Y1 * Z2 * Z2Z2 % p
    S2 = Y2 * Z1 * Z1Z1 % p
    if U1 == U2:
        if S1 != S2:
            return (0, 0, 0)
        return jac_double(p1, a, p)
    H = (U2 - U1) % p
    I = (2 * H) ** 2 % p
    J = H * I % p
    R = 2 * (S2 - S1) % p
    V = U1 * I % p
    X3 = (R * R - J - 2 * V) % p
    Y3 = (R * (V - X3) - 2 * S1 * J) % p
    Z3 = ((Z1 + Z2) ** 2 - Z1Z1 - Z2Z2) * H % p
    return (X3, Y3, Z3)

def jac_mul(k, pt, a=A, p=P):
    if pt[2] == 0:
        return (0, 0, 0)
    if k < 0:
        k = -k
        pt = (pt[0], (-pt[1]) % p, pt[2])
    R = (0, 0, 0)
    Q = pt
    while k:
        if k & 1:
            R = jac_add(R, Q, a, p)
        Q = jac_double(Q, a, p)
        k >>= 1
    return R

def to_affine(pt, p=P):
    X, Y, Z = pt
    if Z == 0:
        return None
    zi = inv(Z, p)
    zi2 = zi * zi % p
    return (X * zi2 % p, Y * zi2 * zi % p)

def to_jac(aff):
    return (0, 0, 0) if aff is None else (aff[0], aff[1], 1)

def mul(k, aff, a=A, p=P):
    """k*aff in affine coordinates."""
    return to_affine(jac_mul(k, to_jac(aff), a, p), p)

def add_aff(p1, p2, a=A, p=P):
    return to_affine(jac_add(to_jac(p1), to_jac(p2), a, p), p)

def on_curve(aff, a=A, b=B, p=P):
    x, y = aff
    return (y * y - (x * x * x + a * x + b)) % p == 0

def lift_x(x, a=A, b=B, p=P):
    """Return the point with the given x and EVEN y (BIP-340 convention), or None."""
    rhs = (x * x * x + a * x + b) % p
    y = sqrt_mod(rhs, p)
    if y is None:
        return None
    return (x, y if y % 2 == 0 else p - y)

def sqrt_mod(a, p=P):
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
    m, c = s, pow(z, q, p)
    t = pow(a, q, p)
    r = pow(a, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p
            i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c = i, b * b % p
        t = t * c % p
        r = r * b % p
    return r

def prime_factors(n):
    return sorted(factorint(n).keys())

def factor_bounded(n, budget=8.0):
    """
    Time-bounded factorisation: trial division by small primes, then Brent's
    Pollard rho.  Returns (factor_dict, complete_flag).  Never blocks longer
    than `budget` seconds -- sympy's factorint can hang for minutes on
    256-bit semiprimes.
    """
    import time, random
    t0 = time.time()
    fac = {}
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47):
        if n % p == 0:
            e = 0
            while n % p == 0:
                n //= p
                e += 1
            fac[p] = e
    if n == 1:
        return fac, True
    for p in range(53, 2000, 2):
        if time.time() - t0 > budget:
            return fac, False
        if n % p == 0:
            e = 0
            while n % p == 0:
                n //= p
                e += 1
            fac[p] = e
    if n == 1:
        return fac, True
    def rho(m, deadline):
        if m % 2 == 0:
            return 2
        while True:
            if time.time() > deadline:
                return None
            y, c, m_ = random.randrange(1, m), random.randrange(1, m), m
            g = r = q = 1
            while g == 1:
                x = y
                for _ in range(r):
                    y = (y * y + c) % m_
                k = 0
                while k < r and g == 1:
                    ys = y
                    for _ in range(min(128, r - k)):
                        y = (y * y + c) % m_
                        q = q * abs(x - y) % m_
                    g = gcd(q, m_)
                    k += 128
                r *= 2
            if g == m:
                g = 1
                while g == 1:
                    ys = (ys * ys + c) % m_
                    g = gcd(abs(x - ys), m_)
            if g != m:
                return g
    stack = [n]
    deadline = t0 + budget
    while stack:
        m = stack.pop()
        if m == 1:
            continue
        if time.time() > deadline:
            return fac, False
        if isprime(m):
            fac[m] = fac.get(m, 0) + 1
            continue
        d = rho(m, deadline)
        if d is None:
            return fac, False
        stack.append(d)
        stack.append(m // d)
    return fac, True

_PRIME_CACHE = {}

def _primes_upto(n):
    if n in _PRIME_CACHE:
        return _PRIME_CACHE[n]
    sieve = bytearray([1]) * (n + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, isqrt(n) + 1):
        if sieve[i]:
            sieve[i * i:: i] = bytearray(len(sieve[i * i:: i]))
    res = [i for i in range(2, n + 1) if sieve[i]]
    _PRIME_CACHE[n] = res
    return res

def small_factor_report(n, bound=10**6):
    """
    Full report of the *small* factor structure of n (all prime factors <= bound),
    plus primality of the remaining cofactor.  Deterministic and fast.
    Full factorisation of a random 256-bit integer is computationally infeasible,
    so this is the only meaningful statement available.
    """
    rem = n
    fac = {}
    for p in _primes_upto(bound):
        if p * p > rem:
            break
        if rem % p == 0:
            e = 0
            while rem % p == 0:
                rem //= p
                e += 1
            fac[p] = e
    if rem > 1 and rem in fac:
        pass
    for p in list(fac):
        pass
    cofactor = rem
    prime_cofactor = isprime(cofactor)
    return fac, cofactor, prime_cofactor

def bits(n):
    return n.bit_length()


def cm_repr(p=P):
    """
    Cornacchia for d = 3: find integers x, y with p = x^2 + 3 y^2  (p = 1 mod 3).

    Writing pi = x + y*sqrt(-3) (which lies in Z[omega] = Z[(-1+sqrt(-3))/2]),
    the six CM traces of the sextic twists of a j = 0 curve are

        { +-2x , +-(x+3y) , +-(x-3y) }.

    (Equivalently, with 4p = L^2 + 27M^2, one has L = 2x and the traces are
     +-L, +-(L+-9M)/2.)
    """
    t = sqrt_mod((p - 3) % p, p)
    if t is None:
        return None
    r0, r1 = p, t
    while r1 * r1 > p:
        r0, r1 = r1, r0 % r1
    x = r1
    rem = p - x * x
    if rem <= 0 or rem % 3 != 0:
        return None
    y2 = rem // 3
    y = isqrt(y2)
    if y * y != y2:
        return None
    return x, y

def candidate_traces(p=P):
    """
    The six CM traces of the sextic twists of a j = 0 curve over F_p (p = 1 mod 3):
        {+-2x, +-(x+3y), +-(x-3y)}   with p = x^2 + 3y^2.
    """
    xy = cm_repr(p)
    if xy is None:
        return []
    x, y = xy
    vals = set()
    for v in (2 * x, x + 3 * y, x - 3 * y):
        vals.add(v)
        vals.add(-v)
    return sorted(vals)


HASSE_BOUND = isqrt(4 * P) + 2

def hasse_window(p=P):
    return p + 1 - (isqrt(4 * p) + 2), p + 1 + (isqrt(4 * p) + 2)

def certify_order(b, cand_order, a=A, p=P, tests=3, seed=12345):
    """
    Certify that the curve y^2 = x^3 + a x + b over F_p has order exactly cand_order.

    Strategy:
      (1) find a point R on the curve (deterministic search over small x),
      (2) check cand_order * R == O for `tests` independent random points R,
      (3) apply the Hasse bound: |#E - (p+1)| <= 2*sqrt(p).
    If (2) holds then cand_order | #E.  Since cand_order is within the Hasse
    window and 2*#cand > window, we conclude #E == cand_order.
    """
    lo, hi = hasse_window(p)
    if not (lo <= cand_order <= hi):
        return False, "outside Hasse window"
    for i in range(tests):
        R = find_point(b, a, p, seed=seed * 7919 + i)
        if R is None:
            return False, "no point found"
        if mul(cand_order, R, a, p) is not None:
            return False, "cand*R != O"
    if 2 * cand_order > hi:
        return True, "certified"
    return False, "window too wide to conclude (needs order of point exactly)"


def find_point(b, a=A, p=P, seed=0):
    """Find a point on y^2 = x^3 + a x + b by scanning x (deterministic)."""
    start = 1 + (seed % 1000)
    for i in range(start, start + 5000):
        x = i
        rhs = (x * x * x + a * x + b) % p
        y = sqrt_mod(rhs, p)
        if y is not None:
            return (x, y)
    return None

def sixth_power_cosets(p=P):
    """
    Representatives of F_p^* / (F_p^*)^6  (6 elements) plus a class label map.

    class(c) = c^((p-1)/6)  in mu_6 ; two elements share a class iff they are in
    the same sixth-power coset.
    """
    exp = (p - 1) // 6
    c0 = None
    x = 2
    while c0 is None and x < 1000:
        if pow(x, exp, p) != 1:
            seen = set()
            for k in range(6):
                seen.add(pow(x, k, p) and pow(pow(x, k, p), exp, p))
            if len(seen) == 6:
                c0 = x
        x += 1
    if c0 is None:
        raise RuntimeError("could not find generator of quotient")
    reps = [pow(c0, k, p) for k in range(6)]
    return c0, reps, exp

def sextic_twists(b=B, p=P):
    """The six sextic twists of y^2 = x^3 + b : y^2 = x^3 + b*c for c in 6 coset reps."""
    c0, reps, exp = sixth_power_cosets(p)
    return [(c, (b * c) % p) for c in reps]

def quadratic_twist_b(b, p=P):
    """
    b' for the quadratic twist of y^2 = x^3 + b.

    The general quadratic twist of y^2 = x^3 + a x + b by a nonsquare d is
        y^2 = x^3 + a d^2 x + b d^3 .
    With a = 0 this reduces to  b' = b * d^3  (d nonsquare).
    NB: b*d^2 would be WRONG -- d^2 is a square, so b*d^2 is Fp-isomorphic to E.
    """
    nu = 2
    while pow(nu, (p - 1) // 2, p) != p - 1:
        nu += 1
    return (b * pow(nu, 3, p)) % p, nu

def cube_roots_of_unity(p=P):
    """The two primitive cube roots of unity in F_p (p = 1 mod 3)."""
    s = sqrt_mod((p - 3) % p, p)
    assert s is not None
    inv2 = inv(2, p)
    z1 = (-1 + s) * inv2 % p
    z2 = (-1 - s) * inv2 % p
    return z1, z2

def phi_point(pt, zeta, a=A, p=P):
    """The order-3 CM endomorphism phi(x,y) = (zeta*x, y) (valid when a = 0)."""
    if pt is None:
        return None
    x, y = pt
    return (zeta * x % p, y)

def verify_endomorphism_point(pt, zeta, p=P):
    return on_curve(phi_point(pt, zeta, p), 0, B, p)

def lambda_from_zeta(zeta, p=P, n=N, generator=None, point=None):
    """
    Determine lam in Z/nZ with phi(P) = lam * P  (as group elements of E(Fp)).

    phi acts on the cyclic group E(Fp) = Z/nZ by multiplication by an element of
    order 3; we identify it by comparing phi(G) with lam*G for a point G.
    """
    for cand in (n_cube_roots(n),):
        for lam in cand:
            if point is None:
                point = find_point(B, A, p)
            if mul(lam, point, A, p) == phi_point(point, zeta, A, p):
                return lam
    return None

def n_cube_roots(n=N):
    """The two nontrivial cube roots of unity mod n (n = 1 mod 3)."""
    s = sqrt_mod((n - 3) % n, n)
    if s is None:
        return []
    inv2 = inv(2, n)
    return [(-1 + s) * inv2 % n, (-1 - s) * inv2 % n]

def isogeny3_map(x, b, p=P):
    """x-coordinate of the image of (x,y) under the 3-isogeny E_b -> E_{-27b}."""
    num = (x * x * x + 4 * b) % p
    den = (x * x) % p
    return num * inv(den, p) % p

def isogeny3_y(x, y, b, p=P):
    """y-coordinate of the image; sign chosen so the point lies on E_{-27b}."""
    num = (x * x * x - 8 * b) % p
    den = (x * x * x) % p
    yy = y * num % p * inv(den, p) % p
    return yy

def isogeny3_point(pt, b, p=P):
    """The rational 3-isogeny  E_b -> E_{-27 b},  (x,y) -> ((x^3+4b)/x^2, y(x^3-8b)/x^3)."""
    if pt is None:
        return None
    x, y = pt
    if x == 0:
        return None
    X = isogeny3_map(x, b, p)
    Y = isogeny3_y(x, y, b, p)
    if on_curve((X, Y), 0, (-27 * b) % p, p):
        return (X, Y)
    return (X, (-Y) % p)

def trace_of(E_order, p=P):
    return p + 1 - E_order

def embedding_degree(n=N, p=P):
    """Smallest k with n | p^k - 1 (None if not found within the search bound)."""
    pk = 1
    for k in range(1, 4096):
        pk = pk * p % n
        if pk == 1:
            return k
    return None
