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

