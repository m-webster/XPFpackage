import numpy as np
from common import *
from XPCodes import *

## Take a hypercube S, copy and translate, then join corresponding vertices
def TranslateCube(S):
    temp = []
    D = len(S)+1
    ## i represents the dimension of the faces
    for d in range(D):
        r = None
        if d+1 < D:
            Sd = S[d]
            nzeros = ZMatZeros(np.shape(Sd))
            ## make 2 copies of the faces
            a = np.hstack([Sd,nzeros])
            b = np.hstack([nzeros,Sd])
            r = np.vstack([a,b])
        if d > 0:
            ## join the corresponding faces
            r2 = np.hstack([S[d-1],S[d-1]])
            r = r2 if r is None else np.vstack([r,r2]) 
        temp.append(r)
    return temp

## create a hypercube of dimension n
def HyperCube(D):
    S = ZMat([[1]])
    for i in range(D):
        S = TranslateCube(S)
    return S


# setVerbose(True)
## dimension of Hypercube
D = 4
N = 2
n = 1 << D
C = HyperCube(D)
# print(C)
zDim = 2
xDim = D

## X stabilizers are on faces of dimension xDim
SX = makeXP(0,C[xDim],0) if xDim >= 0 and xDim <=D else ZMat([],2*n+1)
## Z stabilizers are on face of dimension zDim
SZ = makeXP(0,0,C[zDim]) if zDim >= 0 and zDim <=D else ZMat([],2*n+1)

G = np.vstack([SX,SZ])
# print(XP2Str(G,N))
P = max(2,n)
G = XPSetN(G,N,P)
N = P
print('Generators G:')
print(XP2Str(G,P))
C = Code(G,P)
# S = getVal(C,'S')
# print('Canonical Generators S:')
# print(XP2Str(S,P))

Em = getVal(C,'Em')
# print('Orbit Representatives Em:')
# print(ZmatPrint(Em,2))

print('Codespace dimension',len(Em))


LD,FD = getVals(C,['LD','FD'])
print('Phase vector generators:')
print(ZmatPrint(FD,2*P))
print('Corresponding logical operators')
print(XP2Str(LD,P))


# LO = getVal(C,'LO')
# print('Logical Operators')
# for L in LO:
#     print(XP2Str(L,P),C.Fvector(L))

