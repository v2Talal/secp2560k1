import time
from sclab import *

def tm(label, fn):
    t0 = time.time()
    r = fn()
    print(f"{label}: {time.time()-t0:.3f}s -> {r}", flush=True)
    return r

G = find_point(B, A, P, seed=1)
tm("find_point", lambda: find_point(B, A, P, seed=1))
tm("one 256-bit scalar mult", lambda: mul(N, G) is None)
tm("isprime(N) [known prime]", lambda: isprime(N))
cands = candidate_traces(P)
m = P + 1 - cands[0]
tm("isprime(composite 256-bit)", lambda: isprime(m))
tm("cm_repr", lambda: cm_repr(P) is not None)
tm("sqrt_mod", lambda: sqrt_mod(P-3, P) is not None)
