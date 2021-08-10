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

