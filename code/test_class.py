import os
import sys
import importlib.util

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

spec = importlib.util.spec_from_file_location(
    "c32", os.path.join(HERE, "15_classify_hits.py")
)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

p = 2**256 - 2**32 - 977
n = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141

print(m.order_annihilates(p, 7, n, tests=3))

Gx = 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798
Gy = 0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8

R = m.jac_mul(n, (Gx, Gy, 1), p)
print("n*G is infinity:", R[2] == 0)

P = m.lift_x(3, p, 7)
print("point on curve:", P is not None)
