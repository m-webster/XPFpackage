from queue import PriorityQueue
import numpy as np
import itertools
from common import *
from NSpace import *
from XPCodes import *
from XPAlgebra import *


## binary inclusion uv = v
def binInclusion(u,v):
    u,v = ZMat(u),ZMat(v)
    return np.sum(u * v - v,axis=-1)==0

def CPr(CPList):
    if len(CPList) == 0:
        return 0
    return len(CPList[0][2])

def CPRandom(r,m,P):
    vList = np.random.randint(2,size=(m,r))
    pList = np.random.randint(1,P,size=m)
    qList = np.random.randint(2,P,size=m)
    return CPCanonical([(pList[i],qList[i],vList[i]) for i in range(m)])

def CPEmbed(CPList, Sx):
    L = leadingIndices(Sx)
    m,n = np.shape(Sx)
    L = ZMat([set2Bin(n,[l]) for l in L])
    return [(p,q,ZMat(v) @ L) for (p,q,v) in CPList]

## apply CPList operators to state S of precision N
def CPApply(CPList,S,N):
    Sp,Sx,Sz = XPcomponents(S)
    r = CPr(CPList)
    for p,q,v in CPList:
        ix = binInclusion(Sx, v)
        Sp[ix] += 2*N*p//q
    Sp = np.mod(Sp,2*N)
    return XPmergeComponents([Sp,Sx,Sz])

## apply CPList to Sx or |+>^r if Sx not specified
def CPState(CPList,N,Sx=None):
    r = CPr(CPList)
    Sx = np.eye(r,dtype=int) if Sx is None else Sx
    SX = makeXP(0,Sx,0)
    S,ix = OrbitOperator(SX,N)
    return CPApply(CPList,S,N)

## order CPList by increasing weight then lex order
def CPCanonical(CPList):
    temp = []
    for i in range(len(CPList)):
        p,q,v = CPList[i]
        ## ensure q > 0
        if q < 0:
            p,q = -p,-q
        ## p := p modulo q
        p = p % q
        ## ignore trivial phases
        if p > 0 and sum(v) > 0:
            g = np.gcd(p,q)
            p,q = p//g,q//g
            temp.append((p,q,tuple(v))) 
    return temp

def CPEqual(CP1, CP2):
    return set(CPCanonical(CP1)) == set(CPCanonical(CP1))

## return a set of CP operators corresponding to SX of precision N
def phaseFunction(SX,m,N):
    t = logCeil(N,2)
    r = len(SX)+1 if N != 2**t else t+1
    S, uList = OrbitOperator(SX,N,t=r)
    XPm = makeXP(0,m,0)
    S = XPMul(S,XPm,N)
    temp = []
    pList = XPp(S)
    wList = np.sum(uList,axis=-1)
    for i in np.argsort(wList):
        u,p = uList[i],pList[i]
        if p > 0:
            ix = binInclusion(uList, u)
            pList[ix] = np.mod(pList[ix]-p,2*N)
            temp.append((p,2*N,u))
    return CPCanonical(temp)

## columns are binary reps of subsets of n of size 1 to m
def Mrm(r,m):
    A = [set2Bin(r,s) for k in range(m) for s in itertools.combinations(range(r),k+1)]
    return np.transpose(A)

## optimised Mrm include only binary strings required for optimised CP2SX embedding
def MrmOpt(CPList,CZOpt):
    r = CPr(CPList)
    M = set()
    Q = list()
    V = set()
    for i in range(len(CPList)):
        p,q,v = CPList[i]
        m = weight(v)
        s = tuple(bin2Set(v))
        if s not in V:
            Q.append((m,s))
            V.add(s)
            if m > 1 and not CZOpt[i]:
                M.add((m,s))
    while len(Q) > 0:
        m, s = Q.pop()
        if m > 2:
            for s in itertools.combinations(s,m-1):
                if s not in V:
                    Q.append((m-1,s))
                    V.add(s)
                    M.add((m-1,s))
    M = ZMat([set2Bin(r,s) for m,s in sorted(M)],r)
    # print(func_name(),np.shape(M))
    return np.hstack([np.eye(r,dtype=int),np.transpose(M)])
    


## determine if we can use an optimised embedding for CP operator
## CP(pi/qi,vi) where mi = wt(vi)
# def CPopt(mi,pi,qi):
#     return (pi,qi) == (1,2)

## check if the operator (p,q,v) is a controlled Z operator
def isCZ(p,q,v):
    return 2*p - q == 0

## check if we can use the optimisation for generalised CZ operators
def CZOptimisible(CPList):
    CPList = [(p,q,ZMat(v)) for (p,q,v) in CPList]
    temp = ZMatZeros(len(CPList))
    for i in range(len(CPList)):
        pi,qi,vi = CPList[i]
        ## check if we have a CZ operator
        if isCZ(pi,qi,vi):
            temp[i] = 1
            for j in range(len(CPList)):
                pj,qj,vj = CPList[j]
                ## check if supp(vi) \subset supp(vj)
                if i != j and (binInclusion(vj,vi)):
                    temp[i] =0
                    break
    return temp

## convert set of edges, weights and phase base P
## to a set of CP operators
def EW2CP(edges,weights,P):
    ## vertices
    vertices = set()
    ## process edges - update emax and vertices
    for e in edges:
        vertices.update(set(e))
    ## r is the number of distinct vertices
    r = len(vertices)
    ## make dictionary of vertices, indexed by sort order 
    vertices = sorted(vertices)
    vertices = {vertices[i]:i for i in range(len(vertices))}
    ## turn edges into binary representation of subsets of vertices
    edges = [[vertices[v] for v in e] for e in edges]
    edges = [set2Bin(r,e) for e in edges]
    ## update weights - default to 1 for hypergraph states
    if weights is None:
        weights = [1]*len(edges)
    if P is None:
        wMax = max(weights) + 1
        t = logCeil(wMax)
        P = 1 << t
    temp = [(weights[i],P,edges[i]) for i in range(len(edges))]
    return CPCanonical(temp)

def alternating_vector(M):
    return (-1)**np.mod(np.sum(M,axis=0),2)

def CP2Str(CPList):
    temp = [(f'CP({p}/{q},{ZMat2str(v)})') for (p,q,v) in CPList]
    return "\n".join(temp)

def CP2SX(CPList,optimised):
    CPList = CPCanonical(CPList)
    r = CPr(CPList)
    m = 1
    N = 1
    
    CZOpt = CZOptimisible(CPList) if optimised else [0] * len(CPList)
    for i in range(len(CPList)):
        pi,qi,vi = CPList[i]
        mi = weight(vi)
        if mi == 1 :
            Ni = qi // 2 if qi > 2 and qi % 2 == 0 else qi 
        else:
            Ni = qi * (1 << (mi-2))
        if CZOpt[i]:
            mi = mi - 1
        m = max(m,mi)
        N = np.lcm(N,Ni)
    if m > 1 and N % 2 == 1:
        N = 2 * N
    ## M^r_m for defining Embedding Operator
    M = MrmOpt(CPList,CZOpt) if optimised else Mrm(r,m)
    ## n - number of qubits in embedding
    r,n = np.shape(M)
    ## alternating vector
    a = alternating_vector(M)

    SXx = M 
    SXz = ZMatZeros((r,n))
    SXp = ZMatZeros(r)
    for i in range(len(CPList)):
        (pi,qi,vi) = CPList[i]
        wi = ZMat([1 if np.sum(u*vi - u) == 0 else 0 for u in np.transpose(M)])
        mi = weight(vi)
        for j in range(r):
            xj = M[j]
            if mi == 1:
                zj = pi*wi*xj*2*N//qi
                SXz[j] = SXz[j] - zj 
                SXp[j] += np.sum(zj)
            else:
                axj = vi[j]*a * (xj - 1) if CZOpt[i] else a * xj
                zj = np.mod(pi*2*N//(qi * 2**(mi-1))*axj*wi,N)
                SXz[j] = SXz[j] + zj 
    ## Eliminate any qubits where no phase is applied by SX
    if optimised:
        ix = [i for i in range(n) if i < r or np.sum(SXz,axis=0)[i]> 0]
        SXz = SXz[:,ix]
        SXx = SXx[:,ix]
        r, n = np.shape(SXx)
    
    SX = XPmergeComponents([SXp,SXx,SXz])
    SX = XPRound(SX,N)
    ## construct K^r_m = Ker(M^r_m)
    SZz = kerIA(SXx)
    SZ = makeXP(0,0,SZz * N//2)
    return SX,SZ,N
