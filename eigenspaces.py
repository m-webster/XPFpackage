from common import *
from XPAlgebra import *
from XPCodes import *
import itertools as iter

## turn Z component into partition
def z2part(z,N):
    temp = ZMatZeros(N)
    for a in z:
        temp[a]+=1
    return temp

## eigenspace of N-z is the same as z
## only include one such operator
def minPartition(z,N):
    z = ZMat(z)
    z0 = np.mod(N - z,N)
    p = z2part(z,N)
    p0 = z2part(z0,N)
    return ZMat2str(p0) <= ZMat2str(p)

def diagOps(N,n):
    r = ZMat(list(range(N))) 
    m = [x if x < N//2 else N-x for x in r]
    ix = argsort(m)
    r = r[ix]
    m = np.gcd(N,r)
    # print('m',m)
    ix = argsort(m,reverse=True)
    r = r[ix]
    for z in iter.combinations_with_replacement(r, n):
        if minPartition(z,N):
        # if True:
            yield XPmergeComponents([0,ZMatZeros(np.shape(z)),z])

def diagOpReport(N,n):
    Eclass = dict()
    Eqclass = dict()
    LXxclass = dict()
    data = []
    data.append(['Operator','degree','|Em|','|Eq|','|LXx|','E Class','Eq Class','LXx Class','Eq','LXx'])
    for A in diagOps(N,n):
        Em = getEm(ZMat2D([A]),N,display=False)
        Eq,LXx  = cosetDecomposition(Em)
        if isZero(LXx):
            LXx = ZMat([])
        Eqstr = ",".join(ZMat2str(Eq))
        LXxstr = ",".join(ZMat2str(LXx))
        if Eqstr not in Eqclass:
            Eqclass[Eqstr] = len(Eqclass)
        if LXxstr not in LXxclass:
            LXxclass[LXxstr] = len(LXxclass)
        c = (Eqclass[Eqstr],LXxclass[LXxstr])
        if c not in Eclass:
            Eclass[c] = len(Eclass)
        data.append([XP2Str(A,N),XPDegree(A,N),len(Em),len(Eq),len(LXx),Eclass[c],Eqclass[Eqstr],LXxclass[LXxstr]," "+Eqstr," "+LXxstr])
    print(ZmatTable(np.array(data),rowdiv=1,coldiv=1))

def espaceDims(N,n):
    dims = dict()
    for A in diagOps(N,n):
        Em = getEm(ZMat2D([A]),N)
        d = len(Em)
        dims[d] = XP2Str(A,N)
        # if d not in dims:
        #     dims[d] = XP2Str(A,N)
    data = []
    data.append(['Operator','|Em|'])
    for d in sorted(dims.keys()):
        data.append([dims[d],str(d)])
    print(ZmatTable(np.array(data),rowdiv=1,coldiv=1))


# N,n = 4,5
# diagOpReport(N,n)