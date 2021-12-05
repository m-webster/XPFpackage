import numpy as np
import itertools
import sys

from numpy.core.fromnumeric import repeat
from common import *
# from NSpace import *

###############################################
##            XP Operator Algebra            ##
###############################################

## return number of qubits XP operator A acts on
def XPn(A):
    s = np.size(A,-1)
    if s % 2 == 1:
        return (s-1)//2
    else:
        return None

## return components of A
## note that A can be a single operator or array of operators
def XPcomponents(A,n=None):
    A = ZMat(A)
    if n is None:
        n = XPn(A)
    slices = [-1,slice(n),slice(n,2*n)]
    if np.ndim(A) > 1:
        ## multi-dimensional array
        return [A[:,si] for si in slices]
    ## 1-dimensional array
    return [A[si] for si in slices]

# merge components strict
# make an XP operator, dimensions of components must match
def XPmergeComponents(C):
    n = np.shape(C[1])[-1]
    C[0] = np.transpose([C[0]])
    return ZMat(np.hstack([C[1],C[2],C[0]]),2*n+1)

# merge components non-strict
# tries to make an XP operator, even if dimensions of components vary
def makeXP(p,x,z):
    C = [p,x,z]
    ## if p is a vector, turn it into a 2D column vector
    if np.ndim(C[0]) == 1:
        C[0] = colVector(C[0])
    ## find max dimension of components
    maxDim = np.amax([np.ndim(c) for c in C])
    ## convert x,z into 2D arrays
    C[1],C[2] = ZMat2D(C[1]),ZMat2D(C[2])
    ## broadcast to same shape
    C = [np.array(a) for a in np.broadcast_arrays(*C)]
    ## get first column of p
    C[0] = C[0][:,0]
    C = XPmergeComponents(C)
    ## If we were given 0D or 1D inputs, return single row
    if maxDim < 2:
        C = C[0]
    return C

## force list of XP operators A into 2D form
## oneD records whether original arrays were both 1D
def XP2D(A):
    A1D = typeName(A) not in {'tuple','list'}
    if A1D:
        A = [A]
    dims = [np.ndim(a) for a in A]
    maxdim = np.amax(dims)
    mindim = np.amin(dims)
    oneD = maxdim == 1 and mindim == 1
    A = [ZMat2D(a) for a in A]
    if A1D:
        A = A[0]
    return A, oneD

# def XPsetComponent(A,i,val,n=None):
#     C = XPcomponents(A,n)
#     C[i] = val
#     return XPmergeComponents(C)

def XPp(A,n=None):
    return XPcomponents(A,n)[0]

def XPx(A,n=None):
    return XPcomponents(A,n)[1]

def XPz(A,n=None):
    return XPcomponents(A,n)[2]

## put XP operator A into unique vector form - ie p in Z_2N, x in Z_2, z in Z_N
def XPRound(A,N):
    p,x,z = XPcomponents(A)
    return XPmergeComponents([np.mod(p,2*N),np.mod(x,2),np.mod(z,N)])

## rescale XP operator g with precision N to precision P
## return False if this is not possible
def XPSetNsingle(A,N,P):
    A = XPRound(A,N)
    p,x,z = XPcomponents(A)
    if P > N:
        if P % N > 0:
            return False
        F = P // N
        p = F*p
        z = F*z      
    if N > P:
        F = N // P
        if N % P > 0 or np.sum(np.mod(p,F))> 0 or np.sum(np.mod(z,F)) > 0:
            return False
        p = p//F
        z = z//F
    return XPmergeComponents([p,x,z])

def XPSetN(A,N,P):
    n = XPn(A)
    ## convert A to 2D
    A,oneD = XP2D(A)
    m = len(A)
    ## convert N, P to 1D
    N,P = ZMat1D(N,m),ZMat1D(P,m)
    temp = []
    for i in range(m):
        res = XPSetNsingle(A[i],N[i],P[i])
        if res is False:
            return res
        temp.append(res)
    if oneD:
        temp = temp[0]
    return ZMat(temp,2*n +1)

## Return string representation of XP operator A with precision N
def XP2Str(A,N):
    A = ZMat2D(A)
    if len(A) == 0:
        return 'None'
    C = XPcomponents(A)
    C[0] = np.transpose([C[0]])
    cMod = [2*N,2,N]
    p,x,z = [ZMat2str(C[i],cMod[i]) for i in range(3)]
    temp = [f'XP_{str(N)}({p[i]}|{x[i]}|{z[i]})' for i in range(len(p))]
    return "\n".join(temp)

## x is string of digits, c is the operator label (eg X,Z,S,T etc) - for insertion into TeX
def row2strAlt(x,c):
    if np.ndim(x) == 1:
        x = ZMat([x])
    p = vp2strAlt(x)
    return ["".join(["" if x[j][i] ==0 else f'{c}_{{{i+1}}}{p[j][i]}' for i in range(len(x[j]))]) for j in range(len(x))]

def p2strAlt(p):
    return  "" if p < 2 else f'^{{{p}}}'

def phase2strAlt(p):
    w = '' if p ==0 else '\omega' 
    return  w if p < 2 else f'{w}^{{{p}}}'    

vp2strAlt = np.vectorize(p2strAlt)
vphase2strAlt = np.vectorize(phase2strAlt)

## alternative way of printing XP operators - using old-school X,Z,S,T etc - for insertion into TeX
def XP2StrAlt(A,N):
    if np.ndim(A) == 1:
        A = ZMat([A])
    # print('A',A)
    ## for diag operator use Z,S,T or P otherwise 
    symdict = {2:'Z',4:'S',8:'T'}
    c = symdict[N] if N in symdict else 'P'
    p,x,z = XPcomponents(A)
    x,z = row2strAlt(x,'X'),row2strAlt(z,c)
    p = vphase2strAlt(p)
    temp = [f'{p[i]}{x[i]}{z[i]}' for i in range(len(p))]
    return "\\\\\n".join(temp)    

# Input: string representing XP operator in format 'XP_N(p|x|z)'
# Output: XP operator in vector format
def str2XP(mystr):
    mystr = mystr.replace("(","|")
    mystr = mystr.replace("\n",",")
    mystr = mystr.replace(";",",")
    punc = ["|"," ",","]
    ## strip everything except allowed punctuation and numbers
    mystr = "".join([s for s in mystr if s in punc or s.isnumeric()])
    ## split into lines each with an XP operator
    mystr = mystr.split(",")
    p,x,z,N = [],[],[],[]
    for s in mystr:
        s = s.strip()
        Ni,pi,xi,zi = [a.strip() for a in s.split("|")]
        N.append(int(Ni))
        p.append(int(pi))
        x.append(str2ZMat(xi))
        z.append(str2ZMat(zi))
    A = XPmergeComponents([p,x,z])
    if len(p) == 1:
        A,P = A[0],N[0]
    else:
        P = np.lcm.reduce(N)
        A = XPSetN(A,N,P)
    return A,P

## return identity XP operator on n qubits
def XPI(n):
    return ZMatZeros(2*n+1)

## check if XP operator A is diagonal
def XPisDiag(A):
    return np.where(np.sum(XPx(A),axis=-1)==0,1,0)

## check if XP operator A is identity
def XPisI(A):
    return np.where(np.sum(A,axis=-1)==0,1,0)

## return D(z) where D is the antisymmetric operator
def XPD(z):
    z = ZMat(z)
    x = ZMatZeros(np.shape(z))
    return XPmergeComponents([np.sum(z,axis=-1),x,-z])

## test binary function on XP operators
def XPTestBinary(A,B,N,C,f):
    # print(func_name(),A,B,N,C,f)
    if np.ndim(A) > 1:
        return [f(A[i],B,N,C[i]) for i in range(len(A))]
    if np.ndim(B) > 1:
        return [f(A,B[i],N,C[i]) for i in range(len(B))]
    return f(A,B,N,C)

## test function for XPMul
def XPMulTest(A,B,N,C):
    # print(func_name(),A,B,N,C)
    AM = XP2Mat(A,N)
    BM = XP2Mat(B,N)
    CM = XP2Mat(C,N)
    # print(AM,BM,CM)
    return np.isclose(AM @ BM, CM).all()

## Multiply two XP operators A, B with precision N
def XPMul(A,B,N,C=False):
    if C is not False:
        check = XPTestBinary(A,B,N,C,XPMulTest)
        return np.all(check)
    [A,B],oneD = XP2D([A,B])
    z1 = XPz(A)
    x2 = XPx(B)
    C = XPRound(A+B+XPD(2*x2*z1),N)
    return C[0] if oneD else C

## test function for XPPower
def XPPowerTest(A,N,d,C):
    return np.isclose(XP2Mat(C,N), np.linalg.matrix_power(XP2Mat(A,N),d)).all()

## Raise XP operator A of precision N to power d
def XPPower(A,N,d,C=False):
    if C is not False:
        if np.ndim(A) == 1:
            A,C,d = ZMat([A]),ZMat([C]),ZMat([d])
        if np.size(d) == 1:
            d = np.broadcast_to(d,len(A))
            # print(func_name(),d)
        check = [XPPowerTest(A[i],N,d[i],C[i]) for i in range(len(A))]
        return np.all(check)
    p,x,z = XPcomponents(A)
    d = np.transpose([d])
    a = d % 2
    return XPRound(d*A + XPD((d-a)*x*z),N) 

## return A^-1
def XPInverse(A,N,C=False):
    if C is not False:
        AC = XPMul(A,C,N)
        check = XPisI(AC)
        return np.all(check)
    p,x,z = XPcomponents(A)
    return XPRound(-A-XPD(2*x*z),N)

## square of XP operator
def XPSquare(A,N,C=False):
    if C is not False:
        check = ZMatEqual(XPMul(A,A,N),C)
        return np.all(check)
    p,x,z = XPcomponents(A)
    return XPRound(2*A+XPD(2*x*z),N)

## commmutator of XP operators A and B
def XPCommutator(A,B,N,C=False):
    if C is not False:
        AB = XPMul(A,B,N)
        AinvBinv = XPMul(XPInverse(A,N),XPInverse(B,N),N)
        check= ZMatEqual(C,XPMul(AB,AinvBinv,N))
        return np.all(check)
    [A,B],oneD = XP2D([A,B])
    p1,x1,z1 = XPcomponents(A)
    p2,x2,z2 = XPcomponents(B)
    z = x1*z2-x2*z1+2*x1*x2*(z1-z2)
    C= XPRound(XPD(2*z),N)
    return C[0] if oneD else C

def XPConjugate(A,B,N,C=False):
    if C is not False:
        check = ZMatEqual(C,XPMul(XPMul(A,B,N),XPInverse(A,N),N))
        return np.all(check)
    [A,B],oneD = XP2D([A,B])
    p1,x1,z1 = XPcomponents(A)
    p2,x2,z2 = XPcomponents(B)
    z=x1*z2+x2*z1-2*x1*x2*z1
    C= XPRound(B+XPD(2*z),N)
    return C[0] if oneD else C

def CW2Dict(CW):
    CWphases = dict()
    for i in range(len(CW)):
        p0,x0,z0 = XPcomponents(CW[i])
        for j in range(len(x0)):
            CWphases[tuple(x0[j])] = p0[j]
    return CWphases

## apply XP operator A to S and return change in phase component
def Fvector(A,S,N):
    # report('S')
    # report(State2Str(S,N))
    phaseDict = CW2Dict([S])
    # report(func_name(),'phaseDict')
    # report(phaseDict)
    F = XPMul(A,S,N)
    Fx = XPx(F)
    p0 = ZMat([phaseDict[tuple(x)] for x in Fx])
    p1 = XPp(F)
    # report('p0',p0)
    # report('p1',p1)
    return np.mod(p1-p0,2*N)

## distance of XP operator A - ie count of qubits where either z or x component is nonzero
def XPdistance(A):
    p,x,z = XPcomponents(A)
    s = np.where(x+z > 0,1,0)
    return np.sum(s,axis=-1)

## helper function for XPDegree
def XPdeg(A,N):
    g = np.gcd(XPz(A),N)
    return np.lcm.reduce(N//g,axis=-1)

## calculate degre of XP operator with precision N
def XPDegree(A,N,C=False):
    if C is not False:
        AC = XPPower(A,N,C)
        if not (np.sum(XPx(AC)) + np.sum(XPz(AC))) == 0:
            return False
        if C % 2 == 0:
            AC = XPPower(A,N,C//2)
            if (np.sum(XPx(AC)) + np.sum(XPz(AC))) == 0:
                return False
        return True
    M = XPdeg(A,N)
    M2 = 2* XPdeg(XPSquare(A,N),N)
    return np.where(XPisDiag(A),M,M2)

## return fundamental phase of XP operator A with precision N
def XPFundamentalPhase(A,N):
    d = XPDegree(A,N)
    Ad = XPPower(A,N,d)
    return XPp(Ad), d

## Return possible eigenvalues of XP operator A of precision N
## Single value A only
def XPEigenvalues(A,N,C=False):
    Np = 2 * N
    if C is not False:
        AM = XP2Mat(A,N)
        evals = np.linalg.eig(AM)[0]
        evals = {C2phase(z,Np) for z in evals}
        # report('evals',evals)
        C = set(C)
        # report('C',C)
        return evals <= C
    p,d = XPFundamentalPhase(A,N)
    q = np.mod(-p,2*N)
    return [np.mod((Np * j + q) // d,Np) for j in range(d)]

## return m random XP operators on n qubits, of precision N
## N: z component is in Z_N
def XPRandom(n,N,m=None):
    Np = 2 * N
    s = n if m is None else (m,n)
    x = np.random.randint(2, size=s)
    z = np.random.randint(N, size=s)
    p = np.random.randint(Np) if m is None else np.random.randint(Np,size=m)
    return XPmergeComponents([p,x,z])

## Adjust the phase of XP operator A so 1 is an Eigenvalue
def XPSetEval(A,N):
    # report(func_name(),'A',XP2Str(A,N))
    p,x,z = XPcomponents(A)
    f,d = XPFundamentalPhase(A,N)
    # report('f,d',f,d)
    Np = 2*N
    p = np.mod(p - f // d,Np)
    return XPmergeComponents([p,x,z])

## convert XP operator A with precision N into a complex matrix - for testing
def XP2Mat(A,N):
    p,x,z = XPcomponents(A)
    # print(func_name(),p,x,z)
    p = Phase2C(p,N)
    T = None
    for i in range(XPn(A)):
        mi = XPMat(x[i],z[i],N)
        if T is None:
            ## temp <- operator on 1st qubit 
            T = p * mi
        else:
            ## tensor product
            T = np.kron(T,mi)
    return T    

## convert a cos phase into a complex number 
def Cos2C(c,N,C=False):
    c = ZMat(c)
    ## return complex equivalent of cos(c*pi/2N)
    if C is not False:
        return c == C2Cos(C,N)
    ## 2*pi is a full rotation
    return np.where(c ==-1,0.5,np.cos(np.pi * c / N))

## convert a complex number into a cos phase
def C2Cos(a,N,C=False):
    if C is not False:
        return np.isclose(a,Cos2C(C,N))
    a = np.abs(a)
    p = N * np.arccos(a) / np.pi
    ## round phase to nearest int
    pRound = ZMat(np.round(p))
    valid = np.logical_not(np.isclose(a, 0.5))
    return np.where(valid,np.mod(pRound,2*N),-1)

## convert a phase into a complex number 
def Phase2C(p,N,C=False):
    ## return complex equivalent of phase p/2N
    if C is not False:
        return p == C2phase(C,N)
    ## 2*pi is a full rotation
    return np.exp(np.pi * 1j * p / N)

## convert a complex number into a phase
def C2phase(a,N,C=False):
    if C is not False:
        return np.isclose(a,Phase2C(C,N))
    p = N * np.angle(a) / np.pi
    ## round phase to nearest int
    pRound = ZMat(np.round(p))
    lengthValid = np.logical_not(np.isclose(np.abs(a), 0))
    angleValid = np.isclose(np.abs(p-pRound), 0)
    return np.where(np.logical_and(lengthValid,angleValid),np.mod(pRound,2*N),None)

## convert single qubit Z(p) to matrix representation
def PMat(z,N):
    ## P^z is diag(1,w^2z)
    return np.diag([1,Phase2C(2*z,N)]) 

## convert single qubit p X(x) Z(z) to matrix representation
def XPMat(x,z,N):
    temp = PMat(z,N)
    if x == 1:
        ## if x is set, swap rows
        temp = np.flip(temp,0)
    return temp

## action of projector onto +1 Eigenspace of A on state S
def XPProj(A,S,N,C=False):
    # print(func_name(),'S',S)
    if np.ndim(S) == 1:
        S = ZMat([S])
    if C is not False:
        S1,c = C
        # print(func_name(),'S1,c ',S1,c  )
        # print(func_name(),'S1',State2Str(S1,2*N,c))
        S1 = State2C(S1,2*N,c)
        S2 = State2C(S,N)
        S3 = matProjector(A,N) @ S2
        # print('S3',S3)
        S4,c4 = C2state(S3,2*N)
        # print('S4',State2Str(S4,2*N,c4))
        return np.isclose(S1,S3).all()
    if XPisDiag(A):
        AS = XPMul(A,S,N)
        # report('S',S,XPp(S))
        # report('AS',State2Str(AS,N),XPp(AS))
        ## dp = change in phase between S and AS
        dp = np.mod(XPp(AS) - XPp(S),2*N)
        # print(func_name(),'dp',dp,dp==0)
        ## return only elements |e> of S where A|e> = |e>
        ## Double precision
        S = XPSetN(S,N,2*N)
        return S[dp == 0,:],None
    ## filter S by elements where A^2|e> = |e>
    A2 = XPSquare(A,N)
    # print(func_name(),'A2',XP2Str(A2,N))
    S1,c = XPProj(A2,S,N)
    # print(func_name(),'S1',S1)
    S2 = XPMul(A,S1,N)
    # print(func_name(),'A',A)
    # print(func_name(),'S2',S2)
    return StateAdd(S1,S2,N)

## complex version of projector of XP operator A
## return +1 projector 
def matProjector(A,N):
    I = np.eye(1 << XPn(A),dtype=complex)
    evals = XPEigenvalues(A,N)
    vList = [Phase2C(p,N) for p in evals]
    Amat = XP2Mat(A,N)
    temp = dict()
    d = len(vList)
    ## calculate projectors using Schwinger equation
    vi = 1
    B = I
    for vj in vList:
        if not np.isclose(vj,vi).all():
            B = (B @ (vj*I - Amat))/(vj - vi)
    return B


############ States ###############

## generate a random state
def StateRandom(N,n,m=1):
    ## ensure that each x component is unique
    x = np.random.choice(1<<n, size=m, replace=False)
    x = int2ZMat(x,2,n)
    ## phase can be in range [0,2N]
    p = np.random.randint(2*N,size=m)
    ## coefficients are in range [0,N//2] as we require cos(pi*c/N) to be positive
    c = np.random.randint(-1,N//2,size=m)
    S = XPmergeComponents([p,x,ZMatZeros(np.shape(x))])
    return S,c

## Product state (|0> + w^q|1>)^n
def StatePlus(N,n,q=0):
    x = ZMat(list(itertools.product([0,1],repeat=n)))
    S = makeXP(0,x,0)
    if q != 0:
        A = makeXP(0,0,ZMat([q//2]*n))
        S = XPMul(A,S,N)
    return S

def State2Dict(S):
    if np.ndim(S) == 1:
        S = ZMat([S])
    p,x,z = XPcomponents(S)
    # print(func_name(),'x',x)
    x = ZMat2int(x)
    # print(func_name(),'x',x)
    return {x[i]:p[i] for i in range(len(p))}

## calculate (w^p + w^q)/2 - used for projectors of nondiagonal operators
## note that the resulting state is precision 2*N
def PhaseAdd(p,q,N,C=False):
    if C is not False:
        p1,c = C
        z = Phase2C(p,N)+Phase2C(q,N)
        return np.isclose(z,Phase2C(p1,N,c=c))
    ## c is the coefficient cos(c pi/2N)
    c = np.abs(p - q)
    d= np.mod(c,2*N)
    ## resulting phase
    p = p + q
    ## where d > N, we multiply both components by -1 = w^2N
    Nadd = np.where(np.greater(d,N),2*N,0)
    return np.mod(p+Nadd,4*N), np.mod(c+Nadd,4*N) 

## calculate (S1 + S2)/2
def StateAdd(S1,S2,N,C=False):
    if C is not False:
        S3,c = C
        # print('S1+S2',(State2C(S1,N)+State2C(S2,N)))
        # print('S3',State2C(S3,2*N,c))
        return np.isclose((State2C(S1,N)+State2C(S2,N))/2,State2C(S3,2*N,c)).all()
    S1,S2 = ZMat2D(S1), ZMat2D(S2)
    ## check if we have empty states
    S = None
    if len(S1) == 0:
        S = S2
    if len(S2) == 0:
        S = S1
    if S is not None:
        ## double precision
        S = XPSetN(S,N,2*N)
        return S,[-1]*len(S)
    ## make state dictionaries
    D1,D2 = State2Dict(S1),State2Dict(S2)
    ## make sets of the dictionaries
    SD1,SD2 = set(D1.keys()),set(D2.keys())
    ## A - where x is in both S1 and S2
    A = ZMat([[D1[x],D2[x],x] for x in SD1.intersection(SD2)],3)
    ## Add the phases in A to find p, c
    p,c = PhaseAdd(A[:,0],A[:,1],N)
    x = A[:,2]
    ## Otherwise, just double the phases and c is zero 
    B = ZMat([[2*D1[x],-1,x] for x in SD1-SD2],3)
    C = ZMat([[2*D2[x],-1,x] for x in SD2-SD1],3)
    mydata = np.vstack([np.transpose([p,c,x]),B,C])
    ## Exclude x with zero coefficient - cos N pi/2N = 0
    ix = mydata[:,1] != N
    # ix = [a is not None for a in mydata[:,0]]
    mydata = mydata[ix,:]
    p,c,x = [ZMat(mydata[:,i]) for i in range(3)]
    ## convert x from int to array of 0/1
    x = int2ZMat(x,2,n=XPn(S1))
    ## merge components back into state format
    z = ZMatZeros(np.shape(x))
    S = XPmergeComponents([p,x,z])
    return S,c
     

## convert a state into a  complex vector
## convert the state to a tuple representing the state in Cartesian coordinates
def State2C(S,N,c=None,C=False):
    if C is not False:
        S1,c1 = C2state(C,N)
        C1 = State2C(S1,N,c1)
        return np.isclose(C, C1).all()
    S = ZMat2D(S)
    # print(func_name(),'S',np.shape(S))
    # if c is not None:
    #     print(func_name(),'c',np.shape(c))
    ## there are 2**n basis vector, so need a tuple with this number of elements
    temp = np.zeros(1 << XPn(S),dtype=complex)
    if len(S) == 0:
        return temp
    p,x,z = XPcomponents(S)
    x= ZMat2int(x,2)
    c = 1 if c is None else Cos2C(c,N) 
    temp[x] = c * Phase2C(p,N)
    return temp

## convert a complex vector into a state
def C2state(val,N,C=False):
    if C is not False:
        S,c = C
        S1,c1 = C2state(State2C(S,N,c),N)
        return StateEqual(S,S1) and np.array_equal(c,c1)
    n = logCeil(len(val)-1)
    ## absolute value of val vector
    vabs = np.abs(val)
    ## exclude abs values outside err from vmax
    err = 1e-6
    ix = [i for i in range(len(val)) if vabs[i] > err]
    if len(ix) == 0:
        return ZMat([],2*n+1),None
    x = int2ZMat(ix,2,n)
    p = C2phase(val[ix],N)
    c = C2Cos(val[ix],N)
    z = ZMatZeros(np.shape(x))
    S = XPmergeComponents([p,x,z])
    return S,c

def stateAmplitude(S,N,c=None,C=False):
    if C is not False:
        phi = State2C(S,N,c)
        A = np.sum(phi * np.conj(phi))
        return np.isclose(A,C)
    S = ZMat2D(S)
    if c is None:
        ## in this case, it's a state of form \sum_i w^p_i|e_i>
        ## amplitude is the size of the Z-support
        return len(S) 
    ## otherwise, there is a scaling factor which needs to be taken into account
    return np.sum(Cos2C(c,N)**2)

def StateSort(s):
    p,x,z = XPcomponents(s)
    z = ZMatZeros(np.shape(x))
    return ZMatSort(XPmergeComponents([p,x,z])) 

def StateEqual(s1,s2):
    return np.array_equal(StateSort(s1), StateSort(s2))

def State2Str(S,N,c=None,C=False):
    S = ZMat2D(S)
    if len(S) == 0:
        return "None"
    if c is None:
        coeff = [""] * len(S)
    else:
        # coeff = ["" if ci == -1 else f'c{ci}' for ci in c]
        ## display coeff as float to 3 decimal points
        coeff = Cos2C(c,N)
        coeff = np.char.mod("%0.3f",coeff)
        ix = [ci ==0 for ci in c]
        coeff[ix] = ""
    coeff = np.array(coeff)
    p,x,z = XPcomponents(S)
    x = ZMat2str(x,2)
    ix = np.argsort(x)
    x, p,coeff = x[ix],p[ix],coeff[ix]
    omega = ['' if p_i ==0 else f'w{p_i}/{2*N}' for p_i in p]
    eStr =  [f'{coeff[i]}{omega[i]}|{x[i]}>' for i in range(len(S))]
    return "+".join(eStr) 

## transform basis vector e to XP operator with zero phase and z components
def vec2state(e):
    s = np.shape(e)
    p = 0 if len(s)<2 else ZMatZeros(s[0])
    z = ZMatZeros(s)
    return XPmergeComponents([p,e,z])
