import numpy as np
import sys
import itertools

## object helper functions
## simplify display of complex np.arrays
np.set_printoptions(precision=3,suppress=True)

# def argsort(seq):
#     # http://stackoverflow.com/questions/3071415/efficient-method-to-calculate-the-rank-vector-of-a-list-in-python
#     return sorted(range(len(seq)), key=seq.__getitem__)

def func_name():
    """
    :return: name of caller
    """
    return sys._getframe(1).f_code.co_name

def typeName(val):
    return type(val).__name__


#######################################
####            ZMat               ####
#######################################
## Integer Matrices

## make an integer numpy array
## if n is set, check that rows have length n
def ZMat(A,n=None):
    if typeName(A) != 'ndarray' or A.dtype != int:
        A = np.array(A,dtype=int)
    if n is not None:
        s = list(A.shape)
        if s[-1] == 0:
            A= np.empty((0,n),dtype=int)
    return A

## ensure A is a 2-dimensional integer numpy array
def ZMat2D(A):
    A = ZMat(A)
    if np.ndim(A) == 2:
        return A
    if np.ndim(A) == 0:
        return ZMat([[A]])
    if np.ndim(A) == 1:
        return ZMat([A])    
    d = np.shape(A)[-1]
    return np.reshape(A,(-1,d))

## ensure A is a 1-dimensional integer numpy array of length n
def ZMat1D(A,n):
    A = ZMat(A)
    if np.ndim(A) != 1:
        A = np.reshape(A,(-1,))
    if len(A) < n:
        A = np.tile(A,n//len(A))
    if len(A) > n:
        A = A[:n]
    return A



# Input: string of single digit numbers, split by spaces
# Output: integer array
def str2ZMat(mystr):
    if mystr.find(" ") > 0:
        mystr = mystr.split()
    return ZMat([int(s) for s in mystr])

def int2ZMat(A,N=2,n=None):
    ## return an array representation of integer x
    ## x has n bits in base N
    if n is None:
        n = logCeil(np.amax(A),N)
    d = np.ndim(A)
    B = np.expand_dims(A, axis=d)
    B = np.repeat(B,n,axis=d)
    Ni = N ** np.arange(n-1,-1,-1)
    return np.apply_along_axis(func1d=modDiv,axis=d,arr=B,b=Ni,N=N)

def modDiv(a,b,N):
    return np.mod(a//b,N)

def ZMat2int(A,N=2):
    A = ZMat(A)
    n = np.shape(A)[-1]
    Ni = N ** np.arange(n-1,-1,-1)
    return np.apply_along_axis(func1d=np.dot,axis=-1,arr=A,b=Ni)

## print table, inserting dividers at rowdiv, coldiv
def ZmatTable(data,rowdiv=None,coldiv=None):
    # get max lengths of columns
    colLen = np.amax(np.char.str_len(data),axis=0)
    m,n = np.shape(data)
    for i in range(n):
        data[:,i] = np.char.rjust(data[:,i],colLen[i])
    # insert column dividers
    if coldiv is not None:
        mycol = "|"
        data = np.insert(data,coldiv,mycol,axis=-1)
    # merge rows with " " separators
    data = np.array([" ".join(data[i]) for i in range(m)])
    # insert row dividers
    if rowdiv is not None:
        myrow = "-" * len(data[0])
        data = np.insert(data,rowdiv,myrow,axis=0)
    # join rows 
    return "\n".join(data)

def ZMat2str(A,N=None):
    if np.size(A) == 0:
        return ""
    S = np.char.mod('%d', A)
    sep = ""
    if N is None:
        N = np.amax(A) + 1
    if N > 10:
        Nw= len(str(N-1))
        S = np.char.rjust(S,Nw)
        sep = " "
    return np.apply_along_axis(func1d=sepjoin,axis=-1,arr=S,sep=sep)

def sepjoin(a,sep):
    return sep.join(a)

def ZmatPrint(A,N=None):
    return "\n".join(ZMat2str(A,N))

def isiter(vals):
    # can we turn vals into an array?
    return hasattr(vals,'count') or hasattr(vals,'__array__')

def isint(val):
    return val == int(val)

## is x a power of N?
def isPower(x,N=2):
    t = logCeil(x,N)
    return N ** (t-1) == x

## return max(t) where x = N^t
def logCeil(x,N=2):
    i = 0
    while x > 0:
        x = x // N
        i = i+1
    return i

## argsort but allowing for sorting of tuples
def argsort(seq,reverse=False):
    # http://stackoverflow.com/questions/3071415/efficient-method-to-calculate-the-rank-vector-of-a-list-in-python
    return sorted(range(len(seq)), key=seq.__getitem__,reverse=reverse)

def weight(A):
    return sum([1 if a > 0 else 0 for a in A])

def getmn(A):
    n = len(A)
    m = len(A[0]) if n else 0
    return m,n

def makeList(A):
    return rowVector(A)[0]

def ZMatAddZeroRow(A):
    ## append a row of all zeros at the end of A
    A = ZMat2D(A)
    d = np.shape(A)[-1]
    return np.vstack([A,ZMatZeros(d)])

def RemoveZeroRows(A,N=False):
    ## remove any zero rows
    A = ZMat2D(A)
    ix = np.logical_not([isZero(a,N) for a in A])
    return A[ix]

def isConstant(A):
    return np.amax(A) == np.amin(A)

def isZero(A,N=False):
    ## check if the row is all zeros % N
    if N:
        A = np.mod(A,N)
    return np.all(A == 0)

def ZMatI(n):
    ## identity n x n matrix
    return np.eye(n,dtype=int)

## return np array of zeros of length/shape s
def ZMatZeros(s):
    return np.zeros(s,dtype=int)

## sort list of XP operators 
def ZMatSort(A,reverse=False):
    if np.ndim(A) == 1:
        return A
    A = ZMat(A)
    ix = argsort(ZMat2tuple(A),reverse=reverse)
    return A[ix,:]

def ZMat2tuple(A):
    # print(func_name(),A)
    n = np.shape(A)[-1]
    A = np.reshape(A,(-1,n))
    return [tuple(a) for a in A]

## check if two sets of XP operators A and B are equal
def ZMatEqual(A,B):
    if np.shape(A) != np.shape(B):
        return False
    A,B = ZMatSort(A),ZMatSort(B)

    return np.array_equal(A,B)

def leadingIndex(a):
    # report('a',a)
    i = 0
    n = len(a)
    while i < n and a[i]==0:
        i+=1
    return i

def set2Bin(n,t):
    temp = ZMatZeros(n)
    temp[list(t)] = 1
    return temp

## iterator - rows correspond to subsets of [0..n-1] of size between w1 and w2
def BinPowerset(n,w1=None,w2=None):
    w1 = n if w1 is None else min(w1,n)
    wrange = range(w1) if w2 is None else range(w1,w2)
    for w in wrange:
        for t in itertools.combinations(range(n),w):
            yield set2Bin(n,t) 


def RowSpan(A,N=2):
    A = ZMat(A)
    g = ZMat([np.lcm.reduce(N // np.gcd(a,N)) for a in A])
    # report('g',g)
    G = [range(a) for a in g]
    for ix in itertools.product(*G):
        # ix = list(ix)
        yield np.mod(ix @ A,N)

# def SS2row(n,a):
#     b = np.zeros(n,dtype=int)
#     b[list(a)] = 1
#     return b

# def RowSpanComb(A,N=2,wMin=0,wMax=None,returnIx=False):
#     A = ZMat(A)
#     n = len(A)
#     if wMax is None:
#         wMax = n
#     for m in range(wMin,wMax+1):
#         for a in itertools.combinations(range(n),m):
#             ix = SS2row(n,a)
#             x= np.mod(ix @ A,N)
#             if returnIx:
#                 yield x,ix
#             else:
#                 yield x

def colVector(b):
    return  np.reshape(ZMat2D(b),(-1,1))

def rowVector(b):
    return  np.reshape(ZMat2D(b),(1,-1))

def matEqual(A,B,N):
    temp = matAdd(A,-B,N)
    return isZero(temp,N)

def matAdd(A,B,N):
## resize to allow valid addition
    ma,na = A.shape
    mb,nb = B.shape
    m = max(ma,mb)
    n = max(na,nb)
    A = matResize(A,m,n)
    B = matResize(B,m,n)
    return np.mod(A+B,N)

def matMul(A,B,N=False):
    A = ZMat2D(A)
    B = ZMat2D(B)
    ma,na = A.shape
    mb,nb = B.shape
    n = max(na,mb)
    if na < n:
        A = matResize(A,ma,n)
    if mb < n:
        B = matResize(B,n,nb)
    if N is False:
        return A @ B
    else:
        return np.mod(A @ B, N)
    
# def ZMat2D(A,square=True):
#     ## convert general object into np.array over int ##
#     A = np.array(A,dtype=int)
#     if len(A) == 0:
#         return np.array([[]])
#     s = A.shape
#     if square and len(s) == 1:
#         A = np.array([A],dtype=int)
#     return A

def matResize(A,m,n,check=False):
    if check is not False:
        ma,na = check.shape
        return ma == m and na == n
    ma,na = A.shape
    dn = n-na
    dm = m-ma
    if dn == 0 and dm == 0:
        return A
    temp = []
    for i in range(min(m,ma)):
        r = mat2list(A[i])
        if dn > 0:
            r = r + [0] * dn
        if dn < 0:
            r = r[:n]
        temp.append(r)
    for i in range(dm):
        temp.append([0] * n)
    return ZMat2D(temp)

def mat2list(b):
    b = ZMat2D(b)
    return b.flatten().tolist()

# def arr2str(r):
#     return "".join([str(a) for a in r])

# def ZMat2str(A,sep='\n'):
#     temp = [arr2str(r)  for r in A]
#     return sep.join(temp)

def Highbit(r):
# print('r',r)
    i = 0
    while i < len(r) and r[i] == 0:
        i += 1
    return i if i < len(r) else -1

# def matDiag(A, N):
#     m,n = A.shape
#     if not m:
#         return A
#     A = np.mod(A,N)
#     temp = rowVector([0])
#     temp = matResize(temp,n)
#     for i in range(len(A)):
#         r = A[i]
#         ix = Highbit(r)
#         if ix > -1:
#             temp[ix] = r
#     return temp

# def indepSet(A,N):
#     temp = []
#     for b in A:
#         r,a = matResidual(A,b,N)
#         if not isZero(r,N):
#             temp.append(r)
#     return temp

###### Debugging #############

verbose = False

def report(*args ):
    ## print, but only if global verbose setting is True
    global verbose
    if verbose:
        print(*args)

verbose_old = []

def setVerbose(val):
    ## set verbose variable, keep previous value
    global verbose, verbose_old
    verbose_old.append(verbose)
    verbose = val

def unsetVerbose():
    ## return verbose to previous setting
    global verbose, verbose_old
    if len(verbose_old):
        verbose = verbose_old.pop() 

### variable storage ####

def getVal(obj,label):
    # report('getVal',label)
    
    if checkVal(obj,label):
        return getattr(obj,label)
    getlabel = 'get'+label
    if getlabel not in dir(obj):
        return False
    f = obj.__getattribute__(getlabel)
    if type(f).__name__ != 'method':
        return False
    return setVal(obj,label,f())

def checkVal(obj,label):
    if not hasattr(obj,label):
        return False
    return getattr(obj,label) is not None

def setVal(obj,label,val=None):
    setattr(obj,label,val)
    return val

def setVals(obj,labels,vals):
    return [setVal(obj,labels[i],vals[i]) for i in range(len(labels))]

def getVals(obj,labels):
    return [getVal(obj,label) for label in labels]

# N = 8
# n = 3
# A = np.random.randint(N,size=(n,n))
# B = int2ZMat(A)
# print(A)
# print(B)
# print(ZMat2int(B,2))
# print(ZMat2str(A))

# print(ZMat2D(B))

# D = ZMat([],n)
# D = ZMatAddZeroRow(D)
# print(D)
# print(ZMatRemoveZeroRows(D))

# A = ZMatAddZeroRow(A)
# print(A)
# A = ZMatRemoveZeroRows(A)
# print(A)