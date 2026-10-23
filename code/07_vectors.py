"""Final independent checks v2: fix int-from-bytes, clean Schnorr verification, add diagnostics."""
import json, os, time, hashlib, math
from sympy import isprime
T0 = time.time()
def log(*a): print(f"[{time.time()-T0:6.1f}s]", *a, flush=True)
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "data")
J = json.load(open(os.path.join(DATA, "secp2560k1.json"), encoding="utf-8"))
P = int(J["p_hex"],16); N = int(J["n_hex"],16)
Gx = int(J["Gx_hex"],16); Gy = int(J["Gy_hex"],16)
G2x = int(J["G_even"]["Gx_hex"],16); G2y = int(J["G_even"]["Gy_hex"],16)

def jdbl(T,p):
    X,Y,Z=T
    if Y==0 or Z==0: return (0,0,0)
    YY=(Y*Y)%p; YYYY=(YY*YY)%p; S=(4*X*YY)%p; M=(3*X*X)%p; Tt=(M*M-2*S)%p
    return (Tt,(M*(S-Tt)-8*YYYY)%p,(2*Y*Z)%p)
def jadd(T,Q,p):
    X1,Y1,Z1=T; X2,Y2,Z2=Q
    if Z1==0: return Q
    if Z2==0: return T
    z1z1=(Z1*Z1)%p; z2z2=(Z2*Z2)%p
    U1=(X1*z2z2)%p; U2=(X2*z1z1)%p
    S1=(Y1*Z2*z2z2)%p; S2=(Y2*Z1*z1z1)%p
    if U1==U2:
        if S1!=S2: return (0,0,0)
        return jdbl(T,p)
    H=(U2-U1)%p; I=(2*H)%p; I=(I*I)%p; Jx=(H*I)%p
    r=(2*((S2-S1)%p))%p; V=(U1*I)%p
    X3=(r*r-Jx-2*V)%p; Y3=(r*(V-X3)-2*S1*Jx)%p
    Z3=(((Z1+Z2)*(Z1+Z2)-z1z1-z2z2)*H)%p
    return (X3,Y3,Z3)
def to_aff(T,p):
    X,Y,Z=T
    if Z==0: return None
    zi=pow(Z,-1,p); zi2=(zi*zi)%p
    return ((X*zi2)%p,(Y*zi2*zi)%p)
def mul(k,aff,p):
    if k==0 or aff is None: return None
    T=(aff[0],aff[1],1); R_=(0,0,0)
    for bit in bin(k)[2:]:
        R_=jdbl(R_,p)
        if bit=='1': R_=jadd(R_,T,p)
    return to_aff(R_,p)
def add_aff(A,B,p):
    if A is None: return B
    if B is None: return A
    return to_aff(jadd((A[0],A[1],1),(B[0],B[1],1),p), p)
def legendre(v,p): return pow(v,(p-1)//2,p)
def sqrt_mod(a,p):
    a%=p
    if a==0: return 0
    if legendre(a,p)!=1: return None
    return pow(a,(p+1)//4,p) if p%4==3 else None
def lift_x(x, even=True):
    y = sqrt_mod((x*x*x+7)%P,P)
    if y is None: return None
    if (y%2==0)!=even: y=(P-y)%P
    return (x,y)
def shake(b,n=320): return hashlib.shake_256(b).digest(n)
def XOF(data,out=320,tag=None):
    if tag is not None:
        th=hashlib.shake_256(tag).digest(64); data=th+th+data
    return hashlib.shake_256(data).digest(out)
def Xb(i): return i.to_bytes(320,"big")

log("=== ECDSA vector #1 ===")
vec = J["audit_evidence"]["6_ecdsa"]["vector"]
r=int(vec["r"],16); s=int(vec["s"],16); pub=vec["pubkey"]; msg=vec["msg"].encode()
Qx=int(pub[2:],16); Qy=sqrt_mod((Qx*Qx*Qx+7)%P,P)
if Qy is not None and (Qy%2)!=(int(pub[:2],16)%2): Qy=(P-Qy)%P
Qp=(Qx,Qy)
Qd = mul(424242,(Gx,Gy),P)
log("  424242*G x  =", hex(Qd[0])[:20], "...")
log("  published x =", hex(Qx)[:20], "...")
log("  equal x ?", Qd[0]==Qx, "| equal y ?", Qd[1]==Qy, "| equal as points ?", Qd==Qp)
log("  published prefix", pub[:2], " y parity of decompressed", Qy%2)
e = int.from_bytes(shake(msg),'big') % N
si=pow(s,-1,N); u1=(e*si)%N; u2=(r*si)%N
RR = add_aff(mul(u1,(Gx,Gy),P), mul(u2,Qp,P), P)
log("  ECDSA#1 verifies (e=SHAKE256(msg,320) mod N) ?", RR is not None and RR[0]%P==r%P)
log("  low-S ?", s <= N//2)
RR2 = add_aff(mul(u1,(G2x,G2y),P), mul(u2,Qp,P), P)
log("  ECDSA#1 verifies with base -G ?", RR2 is not None and RR2[0]%P==r%P)
v2=J["audit_evidence"]["second_vectors"]["ecdsa"]
r2=int(v2["r"],16); s2=int(v2["s"],16); pub2=v2["pubkey"]; msg2=v2["msg"].encode()
Qx2=int(pub2[2:],16); Qy2=sqrt_mod((Qx2*Qx2*Qx2+7)%P,P)
if Qy2 is not None and (Qy2%2)!=(int(pub2[:2],16)%2): Qy2=(P-Qy2)%P
e2=int.from_bytes(shake(msg2),'big')%N; si2=pow(s2,-1,N)
RRb=add_aff(mul((e2*si2)%N,(Gx,Gy),P), mul((r2*si2)%N,(Qx2,Qy2),P), P)
log("  ECDSA#2 verifies ?", RRb is not None and RRb[0]%P==r2%P)

log("=== SCHNORR vectors ===")
for label, sv in (("#1", J["audit_evidence"]["4_schnorr"]["vector"]), ("#2", J["audit_evidence"]["second_vectors"]["schnorr"])):
    px=int(sv["pubkey_x"],16); rr=int(sv["sig_r"],16); ss=int(sv["sig_s"],16); mm=sv["msg"].encode()
    Pev=lift_x(px,True); Rev=lift_x(rr,True)
    if Pev is None or Rev is None:
        log("  ",label,"lift_x failed"); continue
    ev=int.from_bytes(XOF(Xb(Rev[0])+Xb(Pev[0])+mm,tag=b"K2560SCHNORR/challenge"),'big')%N
    for nm, base in (("G2", (G2x,G2y)), ("G", (Gx,Gy))):
        lhs=mul(ss, base, P)
        rhs=add_aff(Rev, mul(ev,Pev,P), P)
        ok = (lhs==rhs)
        log("  ",label,"s*%s == R + e*P ? %s" % (nm, ok))
        if ok: break
    log("  ",label,"pubkey x-even lift ok, R even y:", Rev[1]%2==0)
