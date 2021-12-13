from common import *
from XPAlgebra import *
from XPCodes import *
from cellulation import *
import itertools as iter

## control qubits are on even coordinates - origin at (1,1), cell size 2
def isControl(v):
    return v[0] % 2 == 0

## SW qubits - origin at (1,1), cell size 2
def isSW(v):
    if isControl(v):
        return False
    return (np.sum(v) // 2) % 2 == 1

## NW qubits - origin at (1,1), cell size 2 
def isNW(v):
    if isControl(v):
        return False
    return (np.sum(v) // 2) % 2 == 0

## return qubit index of all qubits of class qtype
def UJQubitClass(g,qtype='SW'):
    if qtype=='Control':
        test_func = isControl 
    elif qtype == 'NW':
        test_func = isNW
    else: 
        test_func = isSW
    return [g.cell2index[v] for v in g.cells[0] if test_func(v)]

## non-diagonal stabilizers for Union Jack State on cellulation g
## default version based on hypergraph
def UJNonDiag(g):
    n = len(g.cells[0]) + len(g.cells[1])
    G = [0]*len(g.cells[0])
    X = makeXP(0,1,0)
    S = makeXP(0,0,1)
    S3 = makeXP(0,0,3)
    ## vertex operators - X, Z and S3 operators
    for vlabel,vCell in g.cells[0].items():
        Six = []
        Xix = [vCell.index] + [g.cell2index[e] for e in vCell.adj[1]]
        S3ix = []
        for f in vCell.adj[2]:
            for v in g.cells[2][f].adj[0]:
                vix = g.cell2index[v]
                if vix not in Xix:
                    Six.append(vix)
            for e in g.cells[2][f].adj[1]:
                eix = g.cell2index[e]
                if eix not in Xix:
                    S3ix.append(eix)
        gX = applyOperator(n,Xix,X)
        gS = applyOperator(n,Six,S)
        gS3 = applyOperator(n,S3ix,S3)
        G[vCell.index] = (gX + gS + gS3)
    return G

## apply displacement d to vertex v and return a vertex in cellulation g if possible
def findVertex(g,v,d,m=None):
    vi = v + d
    ## cater for torus - origin at (1,1), cell size 2
    if m is not None:
        vi = np.mod(vi-1,2*(m))+1
    return g.get_alias(tuple(vi))

## Modified non-diagonal stabilizers for Union Jack State on cellulation g
## Apply alternating edge diagional stabilizers to eliminate Z operators
def UJNonDiagMod(g,m=None):
    n = len(g.cells[0]) + len(g.cells[1])
    G = [0]*len(g.cells[0])
    X = makeXP(0,1,0)
    S = makeXP(0,0,1)
    S3 = makeXP(0,0,3)
    ## vertex operators - X, Z and S3 operators
    for vlabel,vCell in g.cells[0].items():
        vlabel = ZMat(vlabel)
        Xix = [vCell.index] + [g.cell2index[e] for e in vCell.adj[1]]
        Six = [[],[]]
        CEdges = ZMat([[1,1],[-1,1],[-1,-1],[1,-1]])
        WEdges = ZMat([[2,0],[1,1],[0,2],[-1,1],[-2,0],[-1,-1],[0,-2],[1,-1]])

        ## Control qubits are on even coordinates
        if vlabel[0] % 2 == 0:
            delta = CEdges
            ## create checkerboard pattern of control qubit stabilizers
            SPower = 0 if sum(vlabel)//2 % 2 == 0 else 1
        else:
            ## SW/NW stabilizers
            delta = WEdges
            SPower = 0
        vi = findVertex(g,vlabel,delta[-1],m)

        for j in range(len(delta)):
            vj = findVertex(g,vlabel,delta[j],m)
            e = tuple(sorted([vi,vj]))
            if e in g.cells[1]:
                Six[SPower].append(g.cell2index[e])
            else:
                # print(vi,vj)
                if vi in g.cells[0]:
                    Six[SPower].append(g.cell2index[vi])
                if vj in g.cells[0]:
                    Six[SPower].append(g.cell2index[vj])
            vi = vj
            SPower = 1-SPower

        gX = applyOperator(n,Xix,X)
        gS = applyOperator(n,Six[0],S)
        gS3 = applyOperator(n,Six[1],S3)
        G[vCell.index] = XPRound(gX + gS + gS3,4)
    return G

## Return Stabilizers for Union Jack code on cellulation g
def UnionJackCode(g,mod=False,m=None):
    g.makeCellIndex([0,1])
    N = 4
    Z = makeXP(0,0,2)
    n = len(g.cells[0]) + len(g.cells[1])
    G = UJNonDiagMod(g,m) if mod else UJNonDiag(g)

    ## apply Z operator to edges around each face
    for flabel,fCell in g.cells[2].items():
        ix = [g.cell2index[e] for e in fCell.adj[1]]
        G.append(applyOperator(n,ix,Z))
    
    ## apply Z operator on each edge
    for elabel,eCell in g.cells[1].items():
        ix = [g.cell2index[v] for v in eCell.adj[0]]
        ix.append(g.cell2index[elabel])
        G.append(applyOperator(n,ix,Z))
    return G,N    

def UJAnalysis(m,torus=False):
    save = True
    ## first make a celluation of the m x m square plane
    o,r,c = (1,1),m,m
    faces = tri_square_plane(o,r,c)
    g = Cellulation(faces)
    if torus:
        g = Torus(g)
    
    ## Union Jack State
    Edges = [(0,1,2)]
    Weights = None 
    ## generate stabilizers for Union Jack State using graphState algorithm
    E,W = applyWeightedGraph(g,Edges,Weights)
    GState = graphState(E,W)
    G,N,n = GState.XPCode(),GState.N,GState.n
    g.SX2cell(GState.SXx)
    g.display(save=save,show_qubits=False,show_operator=[G[0],G[28],G[52],G[11],G[38],G[270],G[148]],title='Examples of Stabilizers')
    # Save images of all generators to file
    # for A in G:
    #     g.display(show_operator=A,save=True)

    g = Cellulation(faces)
    if torus:
        g = Torus(g)
    ## canonical generators
    S1 = CanonicalGenerators(G,N)
    print('Canonical Generators - S1')
    print(len(S1))
    ## Modified version of stabilizers make it easier to see the symmetries for open boundaries
    G,N = UnionJackCode(g,mod=True,m=m) if torus else UnionJackCode(g,mod=True)
    ## canonical generators
    S2 = CanonicalGenerators(G,N)
    print('Canonical Generators - S2')
    print(len(S2))
    print('Check that Canonical Genrators are Equal for both methods: S1 == S2')
    print(np.all(np.isclose(S1,S2)))
    # g.display(show_qubits=False,show_operator=[G[0],G[4],G[28],G[59],G[74],G[399],G[412],G[203]],title='Modified Stabilizers')
    g.display(save=save,show_qubits=False,show_operator=[G[0],G[28],G[52],G[11],G[38],G[399],G[412],G[203]],title='Modified Stabilizers')

    # Save images of all generators to file
    # for A in G:
    #     g.display(show_operator=A,save=True)

    ## Here we display the Z_2 symmetries of the Union Jack State
    ## There are 3 classes of qubits - SW, NW and Control
    qClasses = ['SW','NW','Control']
    ## we make Symmetry operators by multiplying together the non-diagonal stabilizer generators
    ## centred on all qubits of a given class
    SymmOps = []
    SX,SZ = splitDiag(G)
    n = XPn(SX)
    # SXx = XPx(SX)
    # n = len(SXx)
    for qcl in qClasses:
        QList = UJQubitClass(g,qcl)
        u = set2Bin(n,QList)
        A = GeneratorProduct(SX,u,N)
        SymmOps.append(A)
        g.display(save=save,show_operator=A,title=qcl + ' Qubit Symmetry')

    ## Now display the product of all 3 symmetries
    A = XPMul(SymmOps[0],SymmOps[1],N)
    B = XPMul(A,SymmOps[2],N)
    g.display(save=save,show_operator=B,title='Product of all 3 Symmetries')

    ## We now consider the Symmetries restricted to the bulk
    ## The easiest way to see these is to calculate the product of the non-diagonal stabilizers 
    ## but exclude those centred on qubits on the boundary
    ## Boundary Qubits
    Bdy = g.boundary()
    if len(Bdy) > 0:
        boundaryQubits = {g.cells[0][g.get_alias(a)].index for a in Bdy[0] }
        B = None
        ## there are no Control qubits on the boundary, so just look at each of the SW and NW classes
        for qcl in ['SW','NW','Control']:
            QList = set(UJQubitClass(g,qcl))
            # QList.intersection_update(boundaryQubits)
            QList.difference_update(boundaryQubits)
            u = set2Bin(n,QList)
            H = GeneratorProduct(SX,u,N)
            # H = None
            # for A in G:
            #     x = XPx(A)
            #     if leadingIndex(x) in QList:
            #         H = A if H is None else XPMul(H,A,N)
            
            g.display(save=save,show_operator=H,title=qcl + ' Symmetry on Bulk')
            B = H if B is None else XPMul(B,H,N)
        ## Now we display the product of all non-diagonal stabilziers centred on qubits on the boundary
        g.display(save=save,show_operator=B,title='Product of all Symmetries on Bulk')


# UJAnalysis(m)