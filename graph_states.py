import numpy as np
import itertools
from common import *
from NSpace import *
from XPCodes import *

## columns are binary reps of subsets of n of size 1 to m
def WC(n,m):
    A = [set2Bin(n,s) for k in range(m) for s in itertools.combinations(range(n),k+1)]
    return np.transpose(A)

class graphState:
    def __init__(self,edges,weights=None,P=None,optimised=True):
        ## vertices
        self.vertices = set()
        self.edgeList = edges
        ## optimise embedding operator
        self.optimised = optimised
        ## maximum edge size
        self.emax = 0
        ## process edges - update emax and vertices
        for e in edges:
            self.emax = max(self.emax,len(e))
            self.vertices.update(set(e))
        ## hypergraph if no weights or edge weight > 2
        self.hyper = self.emax > 2 or weights is None or len(set(weights)) < 2
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
    
    def WCOptimised(self):
        if not self.hyper:
            ## for weighted graph states, we just want the vertices and bin reps of edges
            W = np.vstack([np.eye(self.n,dtype=int),self.edges])
            return np.transpose(W)

        ## subsets of vertices for each edge
        vertexSubsets = dict()
        ## for hypergraphs, we want subsets of size m-1 where m is the edge size
        ## go through the edges
        for e in self.edgeList:
            smax = len(e)
            for k in range(2,smax):
                for s in itertools.combinations(e,k):
                    if k not in vertexSubsets: 
                        vertexSubsets[k] = set()
                    vertexSubsets[k].add(s)
        ## identity represents size 1 subsets - ie vertices
        W = [np.eye(self.n,dtype=int)]
        ## go through vertex subsets of edges of size k
        for k,Sk in vertexSubsets.items():
            ## add binEdge for each subset of vertices
            W.append([self.binEdge(s) for s in sorted(Sk)])
        ## combine and return transpose
        W = np.vstack(W)
        return np.transpose(W)

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
    
    ## Embedded state |\phi> specified by the graphState object
    def EmbeddedState(self,embed=True):
        ## state with no phases applied
        SXx = getVal(self,'SXx') if embed else np.eye(self.n,dtype=int)
        eList = list(RowSpan(SXx))
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

    ## State |\phi> specified by the graphState object on n qubits
    def State(self):
        return self.EmbeddedState(embed=False)

    ## naive method for finding the stabilizers of the state
    def XPCodeNaive(self):
        ## find the state represented
        s = self.EmbeddedState()
        ## CW2LI calculates the stabilizers of the state
        return CW2LI(ZMat([s]),self.N)

    def getSXx(self):
        return self.WCOptimised() if self.optimised else WC(self.n,self.t)

    ## direct method for finding the stabilizers of the state
    def XPCode(self):
        verbose = getVerbose()
        ## X components of non-diagonal stabilizer generators
        SXx = getVal(self,'SXx')
        n = len(SXx[0])
        if verbose:
            report("\nStep 1: Embedding Operator")
            report(f'The embedding operator is E^{self.n}_{self.t} which is based on the matrix W^{self.n}_{self.t} =')
            report(ZmatPrint(SXx,2))

        ## Initialise SX
        SX = makeXP(0,SXx,0)

        ## Initialise SZ
        ## SZz corresponds to a basis Kernel of SXx
        nsp = NSpace(SXx,2)
        nsp.simplifyKer()
        K = nsp.K if len(nsp.K) > 0 else ZMatAddZeroRow(nsp.K)
        SZ = makeXP(0,0,K)

        if verbose:
            report(f'\nStep 2: Generators for CSS Code of precision {self.N} Stabilizing |\psi_0>:')
            report(f'We determine stabilizer generators for the initial embedded state with no relative phases - ie |\psi_0>=E^{self.n}_{self.t}{"|+>"*self.n}.')
            report(f'The X components of the RX are the rows of W^{self.n}_{self.t}.')
            report(f'The Z components of the RZ are a basis of Ker_Z2(W^{self.n}_{self.t}).')
            report('RX =')
            report(XP2Str(SX,2))
            report('RZ =')
            report(XP2Str(SZ,2))

        ## update Sz to be precision N
        SZ = XPSetN(SZ,2,self.N)
        if verbose:
            report('\nStep 3: Transformation of CSS Stabilizers')
            report('We now transform RX, RZ by the controlled Z or controlled Phase operators U_i defining the graph state.')
            report(f'RZ commutes with the U_i so is not changed. We rescale RZ to precision N={self.N} by multiplying Z components by {self.N//2}.')
            report('We conjugate each of the non-diagonal stabilizer generators RX by each of the operators U_i.')
        uincText = 'c_j u_i = c_j' if self.hyper else 'c_j = u or i'
        ## Conjugation of Sx by controlled phase operators corresponding to each edge
        for k in range(len(self.edges)):
            u = self.edges[k]
            w = self.weights[k]
            if verbose:
                report(f'\nApplying Controlled Phase Operator CP(u={ZMat2str(u)},q={w}/{self.P})')
            I = np.eye(self.n,dtype=int)
            ## iterate through each non-diagional stabilizer generator
            for i in range(len(SX)):
                A = SX[i]
                
                if verbose:
                    report(f'\nConjugation of CP(u={ZMat2str(u)},q={w}/{self.P}) with generator SX[{i}]={XP2Str(A,self.N)}')
                ## check if u[i] is set
                if u[i] == 0:
                    if verbose:
                        report(f'u[{i}] = 0 - No update required for this generator.')
                    continue
                ## u_i is the same as u, but with u_i[i] = 0
                if self.hyper:
                    u_i = u - I[i]
                    if verbose:
                        report(f'u_i = {ZMat2str(u_i)}')   
                temp = [['j','c_j',f' {uincText}','z[j]']]
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
                    ## only update temp in verbose mode
                    if verbose:
                        temp.append([str(j),ZMat2str(c_j,2),str(uinc),str(zAdj[j])])
                temp = np.array(temp)
                report(ZmatTable(temp,rowdiv = [1],coldiv=[1,2,3]))
                zAdj = makeXP(0,0,zAdj)
                SX[i] = XPRound(SX[i]+zAdj,self.N)
 
        G = np.vstack([SX,SZ])
        return XPRound(G,self.N)              
        