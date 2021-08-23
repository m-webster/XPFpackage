import numpy as np
import itertools
from common import *
from NSpace import *
from XPCodes import *

## columns are binary reps of subsets of n of size 1 to m
def WC(n,m):
    A = [set2Bin(n,s) for k in range(m) for s in itertools.combinations(range(n),k+1)]
    return np.transpose(A)

## return span of WC(n,t) modulo 2
def basisElts(n,s=1,t=None):
    if t is None:
        t = n
    return RowSpan(WC(n,t))

class graphState:
    def __init__(self,edges,weights=None,P=None):
        ## vertices
        self.vertices = set()
        ## maximum edge size
        self.emax = 0
        ## process edges - update emax and vertices
        for e in edges:
            self.emax = max(self.emax,len(e))
            self.vertices.update(set(e))
        ## hypergraph if no weights or edge weight > 2
        self.hyper = self.emax > 2 or weights is None
        ## number of qubits is the number of distinct vertices
        self.n = len(self.vertices)
        ## make dictionary of vertices, indexed by sort order 
        self.vertices = sorted(self.vertices)
        self.vertices = {self.vertices[i]:i for i in range(len(self.vertices))}
        ## turn edges into binary representation of subsets of vertices
        self.edges = [self.binEdge(e) for e in edges]
        ## update weights - default to 1 for hypergraph states
        self.weights = [1]*len(edges) if self.hyper else weights
        ## P-1 is the max weight of an edge
        self.P = 1 << max(self.weights).bit_length() if (P is None or self.hyper) else P
        ## t+1: max subset size
        ## N: precision of resulting code
        if self.hyper:
            self.t = self.emax - 1
            self.N = 1 << self.t
        else:
            self.t = 2
            self.N = self.P
    
    ## turn edges into binary representation of subsets of vertices
    def binEdge(self,e):
        b = [0]*self.n
        for i in e:
            b[self.vertices[i]] = 1
        return np.array(b,dtype=int)

    ## apply controlled phase operator represented by u, w
    ## to state represented by pList,eList
    def applyPhase(self,pList,eList,u,w):
        for i in range(len(eList)):
            e0 = eList[i][:self.n]
            if tuple(e0 * u) == tuple(u):
                pList[i] += 2*self.N * w // self.P
        return pList
    
    ## calculate state specified by the graphState object
    def State(self):
        ## state with no phases applied
        eList = list(basisElts(self.n,1,self.t))
        pList = np.array([0]*len(eList),dtype=int)
        ## each edge represents a controlled phase operator 
        ## apply to the state successively
        for i in range(len(self.edges)):
            e = self.edges[i]
            w = self.weights[i]
            pList = self.applyPhase(pList,eList,e,w)
        ## round the phase list
        pList = np.mod(pList,2*self.N)
        ## return the final state
        return XPmergeComponents([pList,eList,ZMatZeros(np.shape(eList))])

    ## naive method for finding the stabilizers of the state
    def XPCodeNaive(self):
        ## find the state represented
        s = self.State()
        ## CW2LI calculates the stabilizers of the state
        return CW2LI(ZMat([s]),self.N)

    ## direct method for finding the stabilizers of the state
    def XPCode(self):
        ## X components of non-diagonal stabilizer generators
        SXx = WC(self.n,self.t)
        n = len(SXx[0])

        ## Initialise SX
        SX = makeXP(0,SXx,0)

        ## Initialise SZ
        ## SZz corresponds to a basis Kernel of SXx
        nsp = NSpace(SXx,2)
        nsp.simplifyKer()
        K = nsp.K if len(nsp.K) > 0 else ZMatAddZeroRow(nsp.K)
        SZ = makeXP(0,0,K)
        ## update Sz to be precision N
        SZ = XPSetN(SZ,2,self.N)

        ## Conjugation of Sx by controlled phase operators corresponding to each edge
        for k in range(len(self.edges)):
            u = self.edges[k]
            w = self.weights[k]
            report(f'Applying Controlled Phase Operator CP(u={u},w={w}/{self.P})')
            I = np.eye(self.n,dtype=int)
            ## iterate through each non-diagional stabilizer generator
            for i in range(len(SX)):
                A = SX[i]
                ## check if u[i] is set
                if u[i] == 0:
                    continue
                ## u_i is the same as u, but with u_i[i] = 0
                u_i = u - I[i]
                report(f'\nConjugation of CP(u={u},w={w}/{self.P}) with generator SX[{i}]={XP2Str(A,self.N)}; u_i = {u_i}')
                temp = [['j','c_j',' uinc','z[j]']]
                ## adjustment to z component of SX[i]
                zAdj = np.zeros(n,dtype=int)
                ## j iterates over all qubits
                for j in range(n):      
                    ## c_j is column j of SXx             
                    c_j = SXx[:,j]
                    
                    if self.hyper:
                        uinc = sum(c_j - u_i * c_j)==0
                        w = 1 << (self.emax - sum(u))
                        zAdj[j] = w*(-1)**(sum(c_j)-1) if uinc else 0
                    else:
                        uinc = np.array_equal(c_j, I[i]) or np.array_equal(c_j, u)
                        zAdj[j] = w * (-1)**sum(c_j) if uinc else 0
                    temp.append([str(j),str(c_j),str(uinc),str(zAdj[j])])
                temp = np.array(temp)
                report(ZmatTable(temp,1,[1,2,3]))
                zAdj = XPRound(makeXP(0,0,zAdj),self.N)
                SX[i] = SX[i]+zAdj
        G = np.vstack([SX,SZ])
        return XPRound(G,self.N)              
        