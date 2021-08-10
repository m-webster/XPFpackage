import numpy as np
import NSpace as ns
import itertools
from common import *
from codes import *

def W(n,m):
    ix = [s for s in itertools.combinations(range(n),m)]
    temp = np.zeros((n,len(ix)),dtype=int)
    for i in range(len(ix)):
        for j in ix[i]:
            temp[j,i] = 1
    return temp

def ProdComb(A,r):
    if r ==0:
        return np.ones((1,len(A[0])),dtype=int)
    if r == 1:
        return A
    m = len(A)
    AList = []
    for s in itertools.combinations(range(m),r):
        # print(A[s,:])
        AList.append(np.prod(A[s,:],axis=0))
    # print(r,AList)
    return np.vstack(AList)

def RM(r,m,p=0):
    A = np.hstack([W(m,s) for s in range(m+1)])
    n = np.shape(A)[-1]
    E = np.vstack([ProdComb(A,i) for i in range(r+1)]) if r+1 > 0 else ZMatZeros((1,n))
    return E[p:,p:]

def self_dual(A):
    if isZero(A):
        return 1
    SD = A @ np.transpose(A)
    N = 1
    while np.sum(np.mod(SD,N*2)) == 0 and N < 64:
        N = N*2
    return N

def binom(n,k):
    k = min(k,n-k)
    a = 1
    for i in range(k):
        a = a*(n-i)//(i+1)
    return a

def maxp(m,r):
    rMin = min(r,m-r-1)
    return np.sum([binom(m,i) for i in range(rMin+1)])

def checkParameters(m,r,p):
    C = [m,r,p]
    if min(C) < 0:
        print('Parameter error: all parameters must be positive')
        return False
    if r > m:
        print('Parameter error: r cannot be greater than m')
        return False
    pmax = maxp(m,r)
    print('Maximum value of p =',pmax)
    if p > pmax:
        print('Parameter error: p too large')
        return False
    return True


# ## Reed Muller parameters
# ## m determines the number of qubits - without punctures, the code is on 2^m qubits
# m = 4
# ## r determines the number of X genreators - |SX| = \sum_{0 \le i le r} \binom{m}{i}
# r = 1
# # p is the number of punctures - remove the first p columns and rows
# p = 1
# ## set debugging level
# setVerbose(False)


# def codeSummary(G,N):
#     ## make an XP code
#     C = Code(G,N)
#     M = getVal(C,'LI')
#     LO = getVal(C,'LO')
#     print(f'Logical Operators: Precision N={N}')
#     if len(LO) > 0:
#         for A in LO:
#             print(XP2Str(A,N),isLO(A,M,N),'f vector',ZMat2str(C.Fvector(A),max(11,2*N)))
#     else:
#         print('None')

# if checkParameters(m,r,p):
#     print(f'Reed Muller code with parameters m={m}, r={r}, p={p}')
#     ## We are building a CSS code:
#     ## X components of non-diagonal generators
#     SXx = RM(r,m,p)
#     ## Z components of diag generators - this is dual to SXx
#     SZz = RM(m-r-1,m,p)
#     ## set precision
#     N = 2
#     ## make generators
#     SX = makeXP(0,SXx,0)
#     SZ = makeXP(0,0,SZz)
#     G = np.vstack([SX,SZ])

#     print('Stabilizer Generators')
#     print(XP2Str(G,N),"\n")

#     codeSummary(G,N)
#     P = self_dual(SXx)
#     if P > 1:
#         N2 = P*2
#         print('Natural precision N =',N2)
#         G = XPSetN(G,N,N2)
#         codeSummary(G,N2)