import numpy as np
import random as random
import math as math
# import queue
from common import *

#######################################
#### General GCD/Modulo Arithmetic ####
#######################################

def Quo(a,b,N,check=False):
    if check is None:
        return True
    if check is not False:
        r = (a - b * check)
        return r >= 0 and r < b
    a = a % N
    b = b % N
    if b == 0:
        return None
    return (a // b) % N

def Div(a,b,N,check=False):
    a = a % N
    b = b % N
    if check is None:
        return True
    if check is not False:
        g = np.gcd(b,N)
        if a % g > 0:
            return True
        return (b * check) % N == a
    if b == 0:
        return None
    g = np.gcd(b,N)
    # print('g',g)
    if a % g == 0:
        # print(b / g, a /g)
        r = a % b
        while r > 0:
            a += N
            r = a % b
        return a // b % N
    return None

def Gcdex(a,b,N,check=False):
    if check is None:
        return True
    if check is not False:
        (g,s,t,u,v) = check
        result = (s * a + t * b) % N == g 
        result &= (u * a + v * b) % N == 0 
        result &= (s * v - t * u) % N == 1
        return result
    s = 0
    old_s = 1
    t = 1
    old_t = 0
    r = b
    old_r = a
    while r != 0:
        quotient = old_r // r
        (old_r, r) = (r, old_r - quotient * r)
        (old_s, s) = (s, old_s - quotient * s)
        (old_t, t) = (t, old_t - quotient * t)
    p = np.sign(t * old_s - s * old_t)
    u = p * s
    v = p * t
    g = old_r
    s = old_s
    t = old_t

    return(g,s,t,u,v)

def inverse(a,N,check=False):
    if check is None:
        return True
    if check is not False:
        return a * check % N == 1
    a = a % N
    if a == 0:
        return None
    t = 0
    newt = 1
    r = N
    newr = a
    while newr > 0:
        quotient = r // newr
        (t, newt) = (newt, t - quotient * newt) 
        (r, newr) = (newr, r - quotient * newr)

    if r > 1:
        return None
    if t < 0:
        t = t + N

    return t

def Ann(a,N,check=False):
    if check is None:
        return True
    if check is not False:
        return (N // np.gcd(a,N) - check) % N == 0
    a = a % N
    if a == 0:
        return 1
    u = N // np.gcd(a,N)
    return u % N

def Split(a,N,check=False):
    if check is None:
        return True
    if check is not False:
        primes = PrimeNumbers(N)
        for p in primes:
            if check % p == 0 and a % p ==0:
                return False
            if N // check % p == 0 and a % p !=0:
                return False
        return True
    if N ==0:
        return 0
    a = a % N
    if a == 0:
        return 1
    r = math.ceil(np.log2(np.log2(N))) if N > 1 else 1
# print('N,r',N,r)
    for i in range(r):
        a = a*a % N
    return N // np.gcd(a,N)

def Factorize(a):
    f = []
    m = []
    for x in range(2,math.ceil(a ** 0.5)+1):
        if a % x == 0:
            i = 0
            f.append(x)
            while a % x == 0:
                i += 1
                a = a // x
            m.append(i)
    return f,m

# def Factors(a):
#     f = set()
#     for x in range(2,math.ceil(a ** 0.5)+1):
#         if a % x == 0:
#             f.add(x)
#             f.add(a //x)
#     f.add(a)
#     return sorted(f)

def PrimeNumbers(a):
    if a < 2:
        return []
    f = [2]
    for x in range(1,a // 2):
        isPrime = True
        xi = 2 * x + 1
        i = 0
        maxi = math.ceil(xi ** 0.5)
        while isPrime and i < len(f):
            isPrime = xi % f[i] > 0
            i += 1
        if isPrime:
            f.append(xi)
    return f

def Stab(a,b,N,check=False):
    if check is not False:
        return (np.gcd(a + check * b,N) - np.gcd.reduce([a,b,N])) % N == 0
    a = a % N
    b = b % N
    g = np.gcd.reduce([a,b,N])
    c = Split(a//g,N//g)
    return c % N

# def Stabex(a,N,check=False):
#     a = makeMat(a)
#     if check is not False:
#         check = makeMat(check)
#         g = np.gcd.reduce(list(a) + [N])
#         temp = np.sum(a * check) % N
#         return g == np.gcd(temp,N)
#     if len(a) == 0:
#         return []
#     c = [1]
#     temp = a[0]
#     for i in range(1,len(a)):
#         c.append(Stab(temp,a[i],N))
#         temp = (temp + c[i] * a[i]) % N
#     return c

def Unit(a,N,check=False):
    if check is None:
        return True
    if check is not False:
        return np.gcd(check,N) == 1 and (check * a - np.gcd(a,N)) % N == 0
    a = a % N
    if a == 0:
        return 1
    g = np.gcd(a,N)
    s = Div(g,a,N)
# print('g,s',g,s)
    if g == 1:
        return s
    d = Stab(s,N//g,N)
# print('d,N //g',d,N //g)
    c = (s + d * N // g) % N
    return c

#######################################
####          Testing Functions    ####
#######################################

def BinaryTest(f,N=False):
    if N is False:
        N = random.randint(2,20)
    print('Testing',f,'N=',N)
    temp = [['a','c','Check']]    
    for a in range(N):
        for b in range(N):
            check = f(a,b,N)
            r = f(a,b,N,check)
            print(ZMat2str([a,b,check,r]))
            if not r:
                temp.append([a,b,check,r])
    if len(temp) > 1:
        for r in temp:
            print(ZMat2str(r))
    else:
        print("All OK")
        
def UnaryTest(f,N=False):
    if N is False:
        N = random.randint(2,20)
    print('Testing',f,'N=',N)
    temp = [['a','c','Check']]    
    for a in range(N):
        check = f(a,N)
        r = f(a,N,check)
        print(ZMat2str([a,check,r]))
        if not r:
            temp.append([a,check,f])
    if len(temp) > 1:
        for r in temp:
            print(ZMat2str(r))
    else:
        print("All OK")

def TestRingOps(N=False):
    if N is False:
        N = random.randint(2,20)
    Binaries = [Quo,Div,Gcdex,Stab]
    for f in Binaries:
        BinaryTest(f,N)
    Unaries = [inverse,Ann,Split,Unit]
    for f in Unaries:
        UnaryTest(f,N)    

#######################################
####   Howellization of Matrixes   ####
#######################################
        
def doOperation(A,opData,N=1):
    op, data = opData
    ## eliminate rows j,m
    if op == 'u':
        (j,m,s,t,u,v) = data
        B = ZMat2D([[s,t],[u,v]])
        R = ZMat2D([A[j],A[m]])
        C = matMul(B, R, N)
        A[j] = C[0]
        A[m] = C[1]
    ## replace row[m] with c * row[j]
    if op == 'r':
        (j,m,c) = data
        A[m] = np.mod(c * A[j],N)
    ## replace row[i] with row[i] + c * row[j]
    if op == 'a':
        (j,i,c) = data
        A[i] = np.mod(A[i] + c * A[j],N)
    ## replace row[j] with c * row[j] 
    if op == 'm':
        (j,c) = data
        A[j] = np.mod(c * A[j],N)
    return A


def ReducedEchelon(A,N,check=False):
    ## put matrix A in reduced echelon assuming it is already in echelon form
    ## ie all entries above leading element aj in row j are < aj
    if check is not False:
        T,E = check
        return isZero(E @ A - T,N)
    m,n = getmn(A)
    T = A.copy()
    E = ZMatI(n)
    r = 0
    operations = []
    for k in range(m):
        if T[r][k] != 0:
            c = Unit(T[r][k],N)
            if c != 1:
                opData = ('r',(r,r,c))
                operations.append(opData)
                doOperation(T,opData,N)
                doOperation(E,opData,N)
            for i in range(r):
                c = (-Quo(T[i][k],T[r][k],N)) % N
                if c != 0:
                    opData = ('a',(r,i,c))
                    operations.append(opData)
                    doOperation(T,opData,N)
                    doOperation(E,opData,N)
            r = r+1
    return T,E

def Echelon(A,N,check=False):
    ## put matrix A in echelon (upper triangular) form
    if check is not False:
        T,E = check
        return isZero(E @ A - T,N)
    m,n = getmn(A)
    T = A.copy()
    E = ZMatI(n)
    r = 0
    operations = []
    for k in range(m):
        for i in range(r+1,n):
            g,s,t,u,v = Gcdex(T[r][k],T[i][k],N)
            opData = ('u',(r,i,s,t,u,v))
            operations.append(opData)
            doOperation(T,opData,N)
            doOperation(E,opData,N)            
        if T[r][k] != 0:
            r = r+1
    E = E[:r]
    T = T[:r]
    return T,E                            

def HowellRecursive(A,N,check=False):
    ## recursive version of Howellization
    ## H = P @ A
    ## K = Ker(At)
    m,n = getmn(A)
    if check is not False:
        P,H,K = check
        return isZero(P@A - H,N) and isZero(K @ A,N)
    A1 = A.copy()
    if n < m:
        A1 = matResize(A,m,m)
    Q,U,C,W,r = WeakHowell(A1,0,N)
    P = np.mod(Q @ U @ C,N)
    # print('W',W)
    K = matMul(W,P,N)
    K = K[:,:n]
    H = matMul(P,A,N)
    ## H is in Echelon form - turn into reduced echelon form
    H,E = ReducedEchelon(H,N)
    ## Update P
    P = matMul(E,P,N)
    P = P[:,:n]
    result = (P,H,K)
    return result


## part of Recursive Howell Algorithm
def WeakHowell(A,k,N,check=False):
    m,n = getmn(A)
    if check is not False:
        Q,U,C,W,r = check
        T = np.mod(Q @ U @ C @ A,N)
        t = k+r
        H = T[k:t,:]
        K = W[k:t,k:t]
        S = -W[:k,k:t]
        A_ = A[:k,:]
        return isZero(K@H,N) and isZero(W@T,N) and isZero(A_ - S@H,N)
        
    if isZero(A,N):
        # print('WeakHowell - Zero')
        return ZMatI(n),ZMatI(n),ZMatI(n),ZMatI(n),0
    if m == 1:
        Q,U,C,W,r,a = ZMatI(n),ZMatI(n),ZMatI(n),ZMatI(n),1,A[k][0]
        if A[k][0] == 0:
            T,E = Echelon(A,N)
        for i in range(n):
            if i != k:
                if A[k][0] == 0:
                    c = E[0][i]
                else:
                    c = Stab(a,A[i][0],N)
                C[k][i] = c
                a = (a + c * A[i][0]) % N
        W[k][k] = Ann(a,N)
        if a > 0:
            for i in range(k):
                W[i][k] = -Div(A[i][0],a,N) % N
            for i in range(k+1,n):
                Q[i][k] = -Div(A[i][0],a,N) % N
    else:
        m1 = m // 2
        A1 = A[:,:m1]
        B = A[:,m1:]
        Q1,U1,C1,W1,r1 = WeakHowell(A1,k,N)
        A2 = np.mod(W1 @ Q1 @ U1 @ C1 @ B,N)
        Q2,U2,C2,W2,r2 = WeakHowell(A2,k+r1,N)
        Q,U,C = Combine(Q1,U1,C1,Q2,U2,C2,W1,k,r1,r2,N)
        W = np.mod(W1 + W2 - ZMatI(n), N)
        r = r1 + r2
# print('Checking WeakHowell',WeakHowell(A,k,N,check=(Q,U,C,W,r)))
    return Q,U,C,W,r
                
## combine matrices for the purpose of calculating weak Howell form
def Combine(Q1,U1,C1,Q2,U2,C2,W1,k,r1,r2,N,check=False):
    n = len(Q1)

    k1 = n - r1 - r2 - k

    if check is not False:
        Q,U,C = check
        return isZero(np.mod(Q2 @ U2 @ ((C2-ZMatI(n))@W1 + ZMatI(n)) @ Q1 @ U1 @ C1 - Q @ U @ C,N))

    q_1 = Q1[k+r1:k+r1+r2,k:k+r1]
    q1 = Q1[k+r1+r2:,k:k+r1]
    u1 = U1[k:k+r1,k:k+r1]
    c1 = C1[k:k+r1,:k]
    d_1 = C1[k:k+r1,k+r1:k+r1+r2]
    d1 = C1[k:k+r1,k+r1+r2:]
    q2 = Q2[k+r1+r2:,k+r1:k+r1+r2]
    u2 = U2[k+r1:k+r1+r2,k+r1:k+r1+r2]
    c2 = C2[k+r1:k+r1+r2,:k]
    c_2 = C2[k+r1:k+r1+r2,k:k+r1]
    d2 = C2[k+r1:k+r1+r2,k+r1+r2:]
    w1 = W1[:k,k:k+r1]
    w_1 = W1[k:k+r1,k:k+r1]

    c11 = np.mod(c1 - d_1 @ c2,N) 
    c12 = np.mod(d1 - d_1 @ d2,N)
    u21 = np.mod(u2 @ (q_1 + d2 @ q1 + c2 @ w1 + c_2@w_1) @ u1,N)
    u12 = np.mod(u1 @ d_1,N)
    u22 = np.mod(u2 + u21 @ d_1,N)

    Q = BlockCol([
        BlockRow([ZMatI(n-k1),Zero(n-k1,k1)]),
        BlockRow([Zero(k1,k),q1,q2,ZMatI(k1)])
        ])

    U = BlockCol([
        BlockRow([ZMatI(k),Zero(k,n-k)]),
        BlockRow([Zero(r1,k),u1,u12,Zero(r1,k1)]),
        BlockRow([Zero(r2,k),u21,u22,Zero(r2,k1)]),
        BlockRow([Zero(k1,n-k1),ZMatI(k1)])
        ])

    C = BlockCol([
        BlockRow([ZMatI(k),Zero(k,n-k)]),
        BlockRow([c11,ZMatI(r1),Zero(r1,r2),c12]),
        BlockRow([c2,Zero(r2,r1),ZMatI(r2),d2]),
        BlockRow([Zero(k1,n-k1),ZMatI(k1)])
        ])

# print('Combine Check',Combine(Q1,U1,C1,Q2,U2,C2,W1,k,r1,r2,N,check=(Q,U,C)))
    return Q,U,C


def Zero(m,n=-1):
    ## all zeros n x n matrix
    if n == -1:
        n = m
    return np.zeros((m,n),dtype=int)

def BlockRow(A):
    return np.hstack(A)

def BlockCol(A):
    return np.vstack(A)

## return an index of the independent rows of A
def HowellIncremental(A,N):
    # A = RemoveZeroRows(A)
    H = None
    n = len(A[0])
    for i in range(len(A)):
        a = A[i]
        if H is None:
            B = [a]
            C = []
            H = np.zeros((0,n),dtype=int)
        else:
            B = np.vstack([H,a])

        H1,ops = HowellForm(B,N)
        # print('H1')
        # print(H1)
        H1 = RemoveZeroRows(H1)
        if not matEqual(H,H1,N):
            H = H1
            C.append(i)
            # if len(C) == n:
            #     break
    return C

## Non-recursive version of Howell Algorithm
def HowellForm(A,N):
    operations = []
    A = A.copy()
    A = np.mod(A,N)
    if len(A) == 0:
        return A,operations
    n = len(A[0])
    m = len(A)
    ##Augment A with zero rows to make it square;
    for i in range(m,n):
        A = ZMatAddZeroRow(A)
    m = max(m,n)

    ## Put A in upper triangular form
    for j in range(n):
        for i in range(j+1,m):
            (g, s, t, u, v) = Gcdex(A[j][j],A[i][j],N)
            if (s,t,u,v) != (1,0,0,1):
                opData = ('u',(j,i,s,t,u,v))
                operations.append(opData)
                doOperation(A,opData,N)

    ## Put A in Howell form
    ##Augment A with one zero row;
    A = ZMatAddZeroRow(A)
    
    for j in range(n):
        if A[j][j] > 0:
            c = Unit(A[j][j],N)
            if c!=1:
                opData = ('r',(j,j,c))
                operations.append(opData)
                doOperation(A,opData,N)

            for i in range(j):
                c = (-Quo(A[i][j],A[j][j],N)) % N
                if c != 0:
                    opData = ('a',(j,i,c))
                    operations.append(opData)
                    doOperation(A,opData,N)
            c = Ann(A[j][j],N)

            opData = ('r',(j,m,c))
            operations.append(opData)
            doOperation(A,opData,N)
        else:
            opData = ('r',(j,m,1))
            operations.append(opData)
            doOperation(A,opData,N)
        if not isZero(A[m]):
            for i in range(j+1,n):
                (g, s, t, u, v) = Gcdex(A[i][i],A[m][i],N)
                if (s,t,u,v) != (1,0,0,1): 
                    opData = ('u',(i,m,s,t,u,v))
                    operations.append(opData)
                    doOperation(A,opData,N)

    A = A[:-1]
    return np.mod(A,N),operations

def Howellize(A,N,check=False):
    if check is not False:
        H,P = check
        return matEqual(matMul(P, A,N),H,N)
    m = len(A)
    n = len(A[0]) if m else 0
    ## H is the Howell Basis
    H, rowOps = HowellForm(A,N)
    ## P is the matrix transforming A to H
    P = ZMatI(m)
    P = matResize(P,max(n+1,m+1),m)
    for opData in rowOps:

        doOperation(P,opData,N)
    return (H,P,rowOps)
        

#######################################
####    Linear Algebra Modulo N    ####
#######################################

## Calculate intersection of two affine spaces
## U1 = o1 + <A1> mod N
## U2 = o2 + <A2> mod N
def affineIntersection(A1,o1,A2,o2,N,C=False):
    if C is not False:
        A,o = C
        tocheck = np.mod(o + A,N)
        tocheck = np.vstack([[o],tocheck])
        for v in tocheck:
            v1 = np.mod(v - o1,N)
            b,u = matResidue(A1,v1,N)
            if not isZero(b):
                print(v1, 'Not in span of A1')
                return False
            v2 = np.mod(v - o2,N)
            b,u = matResidue(A2,v2,N)
            if not isZero(b):
                print(v2, 'Not in span of A2')
                return False
        return True

    nsp = NSpace(np.vstack([A1,A2]),N)
    nsp.simplifyH()
    v = np.mod(o1-o2,N)
    b,u = matResidue(nsp.H,v,N)
    if not isZero(b):
        ## there is no solution o1-o2 = - v1@A1 + v2@A2 <=> o1 + v1@A1 = o2 + v2@A2
        return False
    ## residue of o1-o2 wrt A1 = v2@A2
    b,u = matResidue(A1,v,N)
    ## new offset is o2 + v2@A2
    o = np.mod(b + o2,N)
    ## new affine space is intersection
    A = nsIntersection([A1,A2],N)
    return A,o

## Calculate intersection of two affine spaces
## U1 = o1 + <A1> mod N
## U2 = o2 + <A2> mod N
## return o which is in both U1 and U2 or False
def affineIntercept(A1,o1,A2,o2,N,C=False):
    nsp = NSpace(np.vstack([A1,A2]),N)
    nsp.simplifyH()
    v = np.mod(o1-o2,N)
    b,u = matResidue(nsp.H,v,N)
    if not isZero(b):
        ## there is no solution o1-o2 = - v1@A1 + v2@A2 <=> o1 + v1@A1 = o2 + v2@A2
        return False
    ## residue of o1-o2 wrt A1 = v2@A2
    b,u = matResidue(A1,v,N)
    ## new offset is o2 + v2@A2
    o = np.mod(b + o2,N)
    return o

## intersection of multiple N spaces
def nsIntersection(Alist,N):
    if len(Alist) == 0:
        return False
    A = Alist[0]
    for B in Alist[1:]:
        nsp = NSpace(np.vstack([A,B]),N) 
        nsp.getVal('Kt')
        C = matMul(nsp.Kt,np.vstack([A,np.zeros(B.shape,dtype=int)]),N)
        nsp = NSpace(C,N)
        nsp.simplifyH()
        A = nsp.H
    return A

## union of multiple N spaces
def nsUnion(Alist,N):
    if len(Alist) == 0:
        return False
    nsp = NSpace(np.vstack(Alist),N) 
    nsp.simplifyH()
    return nsp.H

## main class of this module for linear algebra modulo N
class NSpace:

    def __init__(self,A,N):
        self.A = ZMat2D(A)
        self.At = np.transpose(A)
        self.N = N
        self.m = len(A)
        self.n = len(A[0]) if self.m > 0 else 0  
        self.getVal('H')   

    # def checkBasis(self,check):
    #     At,H,P,K,T,S,Kt = check
    #     temp = []
    #     if not matEqual(H, matMul(P,self.A, self.N),self.N):
    #         temp.append(['H=P@A', mat2str(H),mat2str(P), mat2str(self.A)]) 
    #     if not isZero(matMul(self.A,np.transpose(K), self.N)):
    #         temp.append(['A@K', mat2str(K), mat2str(self.A)])
    #     if not matEqual(T, matMul(S,At, self.N),self.N):
    #         temp.append(['T=S@At', mat2str(T),mat2str(S), mat2str(At)])    
    #     return temp

    # def makeBasis(self,check=False):
    #     if check is not False:
    #         return self.checkBasis(check)
    #     (P,H,Kt) = HowellRecursive(self.A,self.N)
    #     At = np.transpose(self.A)
    #     (S,T,K) = HowellRecursive(At,self.N)
    #     check = (At,H,P,K,T,S,Kt)
    #     return check

    def getVal(self,a):
        if hasattr(self, a):
            return getattr(self,a)
        if a == 'H':
            self.H,ops = HowellForm(self.A,self.N)
        if a == 'P':
            self.H,self.P,ops = Howellize(self.A,self.N)
        if a in ['T','S']:
            self.T,self.S,ops = Howellize(self.At,self.N)
        if a == 'K':
            Ht = RemoveZeroRows(self.getVal('H'))
            Ht = matResize(Ht,self.n,self.n)
            Ht = np.transpose(Ht)
            P1,H1,self.K = HowellRecursive(Ht,self.N)
        if a == 'Kt':
            Ht = RemoveZeroRows(self.getVal('T'))
            Ht = np.transpose(Ht)
            P1,H1,self.Kt = HowellRecursive(Ht,self.N)
        return getattr(self,a)
            
    # def makeBasisAlt(self,check=False):
    #     if check is not False:
    #         return self.checkBasis(check)
    #     At = np.transpose(self.A)
    #     H,P,ops = Howellize(self.A,self.N)

    #     T,S,ops = Howellize(At,self.N)
    #     Ht = RemoveZeroRows(H)
    #     Ht = matResize(Ht,self.n,self.n)
    #     Ht = np.transpose(Ht)
    #     P1,H1,K = HowellRecursive(Ht,self.N)
    #     Kt=K
    #     check = (At,H,P,K,T,S,Kt)
    #     return check

    def simplifyH(self):
        ix = [i for i in range(len(self.H)) if np.sum(self.H[i]) > 0]
        self.H = self.H[ix]
        if hasattr(self,'P'):
            self.P = self.P[ix]

    def simplifyKer(self):
        self.getVal('K')
        # P,H,K = HowellRecursive(self.K,self.N)
        H,ops = HowellForm(self.K,self.N)
        self.K = RemoveZeroRows(H)
        if isZero(self.K):
            self.K = np.array([H[0]])                  

    ## solve Ax = b modulo N
    def makeOffset(self,b,check=False):
        self.getVal('S')
        if check is not False:
            v,d = check
            return matEqual(d,self.checkOffset(b,v),self.N)
        ## check if b is in the span of T
        o,a = matResidual(self.T,b,self.N)
        if not isZero(o):
            print("makeOffset: not in span residual",o)
        c = matLinearComb(self.S,a,self.N)
        c = matResize(c,1,self.n)
        ## o and d should be all zero if exact solution exists
        ## otherwise, difference between estimate and solution
        d = self.checkOffset(b,c)
        return c[0],d[0]
    
    def checkOffset(self,b,c):
        ## result is all zero if c is an exact solution
        ## otherwise, difference between estimate and solution
        b = colVector(b)
        c = colVector(c)
        v = matMul(self.A,c,self.N)
        return matAdd(b, -v,self.N)


def matResidue(A,v,N):
    v = mat2list(v)
    v = np.hstack([[1],v])
    A = np.hstack([np.zeros((len(A),1),dtype=int),A])
    A = np.vstack([v,A])
    A,P,operations = Howellize(A,N)
    b = A[0,1:]
    # print(operations)
    a = P[0,1:]
    return b,a

def matResidual(A,b,N,check=False):
    ## take vector b
    ## return resdidual after subtracting Howell basis vectors
    ## and vector with multiples of each basis vector
    # b1,a1 = matResidue(A,b,N)
    # print(func_name(),A)
    b = rowVector(b)
    A = ZMat2D(A)
    if check is not False:
        x, o = check
        return matEqual(b,matLinearComb(A,x,N,o))
    b = makeList(np.mod(b,N))
    if len(A) == 0 or len(A[0]) == 0:
        return b,ZMat2D([])
    temp = []
    for i in range(len(A)):
        ix = Highbit(A[i])
        c = 0
        if ix >= 0:
            d = Quo(b[ix],A[i][ix],N)
            c = d if d else 0
        temp.append(c)
        b = np.mod(b - c * A[i],N)
    a = ZMat2D(temp)
    
    # print('a,a1',a,a1)
    # print('b,b1',b,b1)
    return b,a

def matLinearComb(A,b,N,o=[0],check=False):
    ## take vector b representing pattern of basis vectors
    ## return linear combination of basis + offset o
    if check is not False:
        b1,o1 = matResidual(A,check,N)
        return matEqual(b,b1) and matEqual(o,o1)
    o = rowVector(o)
    b =  rowVector(matMul(b, A, N))
    return np.mod(b + o, N)

def testNSP(A,bList,N):
    A = ZMat2D(A)
    nsp = NSpace(A,N)
    nsp.getVal('K')
    nsp.getVal('P')
    nsp.getVal('Kt')
    
    print('A\n',ZMat2str(nsp.A))
    print('H\n',ZMat2str(nsp.H))
    print('P\n',ZMat2str(nsp.P))
    print('Kt\n',ZMat2str(nsp.Kt))
    print('At\n',ZMat2str(nsp.At))
    print('T\n',ZMat2str(nsp.T))
    print('S\n',ZMat2str(nsp.S))
    print('K\n',ZMat2str(nsp.K) )
    # print('H Check\n',ZMat2str(matMul(nsp.P, nsp.A,N)))
    
    for b in bList:
        print("Solving for",b)
        c,d = nsp.makeOffset(b)
        if isZero(d):
            print("Solution Found",c)
        else:
            print("Failed")

    n = nsp.n
    m = nsp.m



