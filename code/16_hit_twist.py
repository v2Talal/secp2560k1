import json, gmpy2, os, sys, time
DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
d = json.load(open(os.path.join(DATA, "paper_minimality.json")))
MODELS = [7,21,63,189,567,1701]

def primes_upto(bound):
    sieve = bytearray(b"\x01") * (bound + 1)
    sieve[0:2] = b"\x00\x00"
    from math import isqrt
    for i in range(2, isqrt(bound) + 1):
        if sieve[i]:
            sieve[i*i::i] = b"\x00" * ((bound - i*i)//i + 1)
    import numpy as np
    arr = np.frombuffer(sieve, dtype=np.uint8)
    return (np.flatnonzero(arr[5:]) + 5)

PS = [int(x) for x in primes_upto(20000000)]
print("primes:", len(PS))

rows = []
for h in d["hits"]:
    if not h["is_class0_b7"]:
        continue
    j = h["j"]; c = h["c"]
    p = (1<<2560) - (1<<320) - c
    N = h["prime_order"]
    t = p + 1 - N
    nt = p + 1 + t
    g = gmpy2.gcd(N, nt)
    spf = None
    for q in PS:
        if nt % q == 0:
            spf = q; break
    rows.append((c, j, N.bit_length(), t.bit_length(), int(g), spf))
    print(f"c={c:>10} j={j:>9} Nbits={N.bit_length()} tbits={t.bit_length()} gcd(N,nt)={int(g)} twist smallest factor={spf}")
