import matplotlib
import numpy as np
from XPCodes import *
from common import *
from NSpace import *
from XPAlgebra import *
from graph_states import *
import pylab as pl
from matplotlib import collections  as mc
import matplotlib.pyplot as plt

plotNo = 0

def midpoint(label):
    m = np.reshape(ZMat(label),(-1,2))
    return tuple(np.average(m,axis=0))

class HyperPlane:
    def __init__(self,label,d):
        self.label = label
        self.d = d
        self.midpoint =  midpoint(label)
        self.adj = [set(), set(), set()]

# def setOrder(HList):
#     mykeys = list(HList.keys())
#     # midpoints = [(HList[k].midpoint[1],HList[k].midpoint[0]) for k in mykeys]
#     midpoints = [(HList[k].midpoint) for k in mykeys]
#     indices = argsort(midpoints)

#     for i in range(len(indices)):
#         HList[mykeys[indices[i]]].index = i

class Cellulation:
    def __init__(self,faces):
        self.clock = 0
        self.aliases = dict()
        self.faceset = faces
        self.generateCells()

    def boundary(self):
        eList = [e for e in self.cells[1] if len(self.cells[1][e].adj[2]) < 2]
        vDict = dict()
        for (a,b) in eList:
            if a not in vDict:
                vDict[a] = set()
            if b not in vDict:
                vDict[b] = set()   
            vDict[a].add(b) 
            vDict[b].add(a)     
        faces = []
        while(len(vDict)> 0):
            a = min(vDict.keys())
            f = []
            while a in vDict:
                b = min(vDict[a])
                vDict[a].discard(b)
                if len(vDict[a]) ==0:
                    del(vDict[a])
                vDict[b].discard(a)
                if len(vDict[b]) ==0:
                    del(vDict[b])
                f.append(b)
                a = b
            faces.append(f)
        return faces

    ## legs of a face ie edged adjcent to vertices of the face, but not in the face itself
    def legs(self,a):
        f = self.cells[2][a]
        # a = {self.get_alias(x) for x in a}
        # print(func_name(),'a',a)
        temp = set()
        for p in a:
            for (x,y) in self.cells[0][p].adj[1]:
                # x = self.get_alias(x)
                # y = self.get_alias(y)
                x, y = min([x,y]),max([x,y])
                xin = 1 if x in a else 0
                yin = 1 if y in a else 0
                if xin + yin == 1:
                    temp.add((x,y))
        # print(func_name(),'temp',temp)
        return temp

    def generateCells(self):
        self.cells = [dict() for i in range(3)]
        facesDone = set()
        facesetNew = []
        for f in self.faceset:
            f = tuple([self.get_alias(v) for v in f])
            fset = tuple(sorted(f))
            if fset not in facesDone:
                facesDone.add(fset)
                facesetNew.append(f)
                self.cells[2][f] = HyperPlane(f,2)
                for i in range(len(f)):
                    a = f[i]
                    j = (i + 1) % len(f)
                    b = f[j]
                    if a > b:
                        (a,b) = (b,a)
                    e = (a,b)
                    if a not in self.cells[0]:
                        self.cells[0][a] = HyperPlane(a,0)
                    if b not in self.cells[0]:
                        self.cells[0][b] = HyperPlane(b,0)
                    if e not in self.cells[1]:
                        self.cells[1][e] = HyperPlane(e,1)
                        
                    # update face
                    self.cells[2][f].adj[0].add(a)
                    self.cells[2][f].adj[0].add(b)
                    self.cells[2][f].adj[1].add(e)

                    # update edge
                    self.cells[1][e].adj[0].add(a)
                    self.cells[1][e].adj[0].add(b)
                    self.cells[1][e].adj[2].add(f)

                    # update vertices
                    self.cells[0][a].adj[1].add(e)
                    self.cells[0][a].adj[2].add(f)

                    self.cells[0][b].adj[1].add(e)
                    self.cells[0][b].adj[2].add(f)     
        for i in range(3):
            # setOrder(self.cells[i])
            self.setIndex(i)
        # print(func_name(), facesDone)
        self.faceset = facesetNew

    ## for qubit index for cells of dimension i
    def setIndex(self,i):
        HList = self.cells[i]
        mykeys = list(HList.keys())
        if i == 0:
            ## When i==0, we look at location of vertices and order by this
            midpoints = [HList[k].midpoint for k in mykeys]
        else:
            ## When i > 0, we form tuples of the indices of the vertices making up the face
            midpoints = [tuple(sorted([self.cells[0][v].index for v in k])) for k in mykeys]
        indices = argsort(midpoints)
        
        for i in range(len(indices)):
            HList[mykeys[indices[i]]].index = i

    # def addTail(self):
    #     faces = self.boundary()
    #     i = 0
    #     for f in faces:
    #         i += 1
    #         p = (0,-i)
    #         q = f[0]
    #         e = (p,q)
    #         ## add e as a neighbour to existing point q
    #         self.cells[0][q].adj[1].add(e)
    #         ## add new point
    #         c = HyperPlane(p,0)
    #         c.adj[1].add(e)
    #         self.cells[0][p] = c
    #         ## add new edge
    #         c = HyperPlane(e,1)
    #         c.adj[0].add(p)
    #         c.adj[0].add(q)
    #         self.cells[1][e] = c
            
    #         print(func_name(),f[0])
    #     for i in range(2):
    #         setOrder(self.cells[i])
    #     return True

    def genus(self):
        F = len(self.cells[2])
        E = len(self.cells[1])
        V = len(self.cells[0])
        return("F = {}, E= {},V = {},F - E + V = {}".format(F,E,V,F - E + V))

    def get_alias(self,a):
        if a in self.aliases:
            if not isinstance(self.aliases[a],set):
                a = self.aliases[a]
        return a

    def identify(self,a,b):
        if a in self.aliases:
            if not isinstance(self.aliases[a],set):
                a = self.aliases[a]
        else:
            self.aliases[a] = set()
        
        if b in self.aliases:
            if not isinstance(self.aliases[b], set):
                b = self.aliases[b]
            self.aliases[a].update(self.aliases[b])
            for c in self.aliases[b]:
                self.aliases[c] = a
                
        self.aliases[a].add(b)
        self.aliases[b] = a


    def displayOp(self,z,Xcomponent=False):
        colList = ['gray','red','orange','blue']
        n = len(z)
        offset = 0.05
        xCoord,yCoord,cList = [],[],[]
        for i in range(n):
            zi = z[i]
            if zi > 0:
                (a,b) = midpoint(self.index2cell[i])
                if Xcomponent:
                    a,b = a+offset,b-offset
                xCoord.append(a)
                yCoord.append(b)
                c = colList[0] if Xcomponent else colList[zi] 
                cList.append(c)
        return xCoord,yCoord,cList

    def display(self,show_qubits=False,show_operator=None,save=False,title=None,showaxis=False):
        eList = list(self.cells[1].keys())
        lc = mc.LineCollection(eList, linewidths=0.5,colors='gray', linestyle='dashed')
        fig, ax = pl.subplots()
        ax.add_collection(lc)
        ax.autoscale()
        ax.margins(0.1)
        if show_qubits:
            xCoord,yCoord,txt = [],[],[]
            for i,l in self.index2cell.items():
                (a,b) = midpoint(l)
                xCoord.append(a)
                yCoord.append(b)
                txt.append(i)
            ax.scatter(xCoord, yCoord,c='dimgray',s=4)
            for i in range(len(txt)):
                ax.annotate(txt[i], (xCoord[i], yCoord[i]))
        if show_operator is not None:
            for A in ZMat2D(show_operator) :
                p,x,z = XPcomponents(A)
                # print(func_name(),x)
                xCoord,yCoord,cList = self.displayOp(x,Xcomponent=True)
                ax.scatter(xCoord, yCoord,c=cList,s=36)
                xCoord,yCoord,cList = self.displayOp(z,Xcomponent=False)
                ax.scatter(xCoord, yCoord,c=cList,s=36)
        global plotNo
        if title is not None:
            plt.title(title)
        if showaxis is False:
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
        if save:
            fig.savefig('./images/' + "{0:03d}".format(plotNo) + '.png', transparent=False, dpi=100, bbox_inches="tight")
            plotNo += 1
            matplotlib.pyplot.close(fig)       
        else: 
            plt.show()

    def makeCellIndex(self,dims=None):
        ## default to edges
        if dims is None:
            dims = [1]
        self.cell2index = dict()
        self.index2cell = dict()
        iMax = 0
        for d in dims:
            for l,c in self.cells[d].items():
                ix = c.index + iMax
                self.cell2index[l] = ix
                self.index2cell[ix] = l
            iMax += len(self.cells[d])

    ## after making a hypergraph, relate index back to original cells
    def SX2cell(self,SXx):
        self.index2cell = sorted(self.cells[0].keys())
        m,n = np.shape(SXx)
        ix2cell = dict()
        cell2ix1 = dict()
        for i in range(n):
            s = SXx[:,i]
            r = [self.index2cell[j] for j in range(m) if s[j]==1]
            r = tuple(sorted(r)) if len(r) > 1 else r[0]
            ix2cell[i] = r
            cell2ix1[r] = i
        cell2ix2 = dict()
        for d in range(3):
            for l,c in self.cells[d].items():
                # r = [l] if d == 0 else l
                r = tuple(sorted(l)) if d > 0 else l
                if r in cell2ix1:
                    ix = cell2ix1[r]
                    cell2ix2[l] = ix
                    c.index = ix
                else:
                    c.index = None
        self.index2cell = ix2cell
        self.cell2index = cell2ix2

def square_unit(o):
    offsets = ZMat([[0,0,2,2],[0,2,2,0]])
    return tuple(tuple(o + x) for x in np.transpose(offsets))

def square_plane(o,r,c):
    (x,y) = o
    faces = []
    for i in range(r):
        for j in range(c):
            faces.append(square_unit((x+2*j,y +2*i)))
    return faces

# hexagonal tiling
def hex_unit(o):
    offsets = ZMat([[2,1,0,1,2,3],[-1,-1,0,1,1,0]])
    return tuple(tuple(o + x) for x in np.transpose(offsets))

def hex_plane(o,r,c):
    faces = []
    loc = ZMat(o)
    for i in range(r):
        loc[0] = o[0]
        loc[1] = o[1] + i * 2
        flip = False
        for j in range(c):
            faces.append(hex_unit(loc))
            loc[0] +=2
            loc[1] += (-1 if flip else 1)
            flip = not flip
    return faces 

def tri_square_plane(o,r,c):
    faces = []
    loc = ZMat(o)
    for i in range(r):
        for j in range(c):
            loc = o + 2 * ZMat([i,j])
            faces.extend(tri_square_unit(loc))
    return faces 

# triangualar tiling within squares
def tri_square_unit(o):
    faces = []
    a = o
    b = o + ZMat([1,1])
    z = ZMat([[0,2],[2,2],[2,0],[0,0]])
    for i in range(4):
        c = o + z[i]
        faces.append((tuple(a),tuple(b),tuple(c)))
        a = c
    return faces

# triangualar tiling
def tri_unit(o,flip=False):
    offsets = ZMat([[1,2,0],[0,2,2]])  if flip else ZMat([[0,1,2],[0,2,0]]) 
    return tuple(tuple(o + x) for x in np.transpose(offsets))

def tri_plane(o,r,c):
    faces = []
    loc = ZMat(o)
    for i in range(r):
        loc[0] = o[0]
        flip = i % 2 == 1
        for j in range(c):
            faces.append(tri_unit(loc,flip))
            loc[0] += 1
            flip = not flip
        loc[1] += 2
    return faces    

def squareTorus(o,r,c):
    faces = square_plane(o,r,c)
    g = Cellulation(faces)
    return Torus(g)

def Torus(g):
    Bdy = g.boundary()[0]
    segments = UDBoundary(Bdy) + LRBoundary(Bdy)
    print('Bdy Segments')
    print(segments)
    print(g.genus())

    S = [(1,0),(3,2)]
    for s1, s2 in S:
        for i in range(len(segments[s1])):
            a, b = segments[s1][i], segments[s2][i]
            g.identify(a,b)
    g.generateCells()
    print(g.genus())
    return g

def Sphere(g):
    Bdy = g.boundary()[0]
    segments = UDBoundary(Bdy) + LRBoundary(Bdy)
    print('Bdy Segments')
    print(segments)
    print(g.genus())

    S = [(0,2),(1,3)]
    for s1, s2 in S:
        for i in range(len(segments[s1])):
            a, b = segments[s1][i], segments[s2][i]
            g.identify(a,b)
    g.generateCells()
    print(g.genus())
    return g

def swapTuple(A):
    return [(b,a) for (a,b) in A]

def LRBoundary(vertices,lr=False):
    top, bottom = UDBoundary(swapTuple(vertices))
    return swapTuple(top),swapTuple(bottom)

def UDBoundary(vertices):
    Xmax = dict()
    Xmin = dict()
    for (x,y) in vertices:
        if x not in Xmax:
            Xmin[x],Xmax[x] = y,y
        if y < Xmin[x]:
            Xmin[x] = y
        if y > Xmax[x]:
            Xmax[x] = y
    XCoord =  sorted(Xmax.keys())
    Ymin = min(Xmin.values())
    Ymax = max(Xmax.values())
    top,bottom = [],[]
    for x in XCoord:
        # if Xmax[x] in Ymax and Xmin[x] in Ymin:
        if Xmax[x] > Ymax - 2 and Xmin[x] < Ymin + 2:
            top.append((x,Xmax[x]))
            bottom.append((x,Xmin[x]))
    return top,bottom

def Corners(vertices):
    s = 1e6
    c = 1
    Cmax = [-s,-s]
    Cmin = [s,s]
    Tmax = [0,0]
    Tmin = [0,0]
    for (x,y) in vertices:
        for i in range(2):
            p = (-1)^i
            a = x + y*c
            if a < Cmin[i]:
                Cmin[i] = a 
                Tmin[i] = (x,y)
            if a > Cmax[i]:
                Cmax[i] = a
                Tmax[i] = (x,y)
            c = - 1/c
    return Tmin + Tmax

def applyOperator(n,ix,A):
    p,x,z = XPcomponents(A)
    bp = len(ix)*p
    bx = ZMatZeros(n)
    bz = ZMatZeros(n) 
    # print(func_name(),A)
    # print(func_name(),ix)
    for i in ix:
        bx[i] = bx[i] + x
        bz[i] = bz[i] + z
    return XPmergeComponents([bp,bx,bz])

def vertex_operator(g,A,B=None):
    G = []
    n = len(g.cells[1])
    for a in g.cells[0].keys():
        ix = [g.cells[1][b].index for b in g.cells[0][a].adj[1]]
        op = applyOperator(n,ix,A)
        if B is not None:
            ## edges on neighbouring plaquettes
            ix2 = {g.cells[1][e].index for b in g.cells[0][a].adj[2] for e in g.cells[2][b].adj[1] }
            ## exclude those already in ix
            ix2.difference_update(set(ix))
            ix2 = list(ix2)
            op += applyOperator(n,ix2,B)
        G.append(op) 
    return G

def plaquette_operator(g,A,B=None):
    G = []
    n = len(g.cells[1])
    for a in g.cells[2].keys():
        ## A is the operator applied to the edges of the face
        ix = [g.cells[1][b].index for b in g.cells[2][a].adj[1]]
        op = applyOperator(n,ix,A)
        ## B is the operator applied to the legs of the face
        if B is not None:
            ix = [g.cells[1][b].index  for b in g.legs(a)]
            op += applyOperator(n,ix,B)
        G.append(op) 
    return G

def SurfaceCode(g,diagOnly=False):
    N = 2
    ## Z operator - apply to edges of vertices
    Z = makeXP(0,0,1)
    G = vertex_operator(g,Z)
    if not diagOnly:
        ## X operator - apply to edges of plaquettes
        X = makeXP(0,1,0)
        G += plaquette_operator(g,X)
    return G,N    

def SemionCode(g,diagOnly=False):
    N = 4
    ## Z operator - apply to edges of vertices
    Z = makeXP(0,0,2)
    G = vertex_operator(g,Z)
    if not diagOnly:
        ## X operator - apply to edges of plaquettes
        X = makeXP(0,1,0)
        ## S operator - apply to legs of plaquettes
        S = makeXP(0,0,1)
        G += plaquette_operator(g,X,S)
    return G,N    

def TQDZ2Code(g):
    N = 4
    ## Z operator - apply to edges of plaquettes
    Z = makeXP(0,0,2)
    G = plaquette_operator(g,Z)
    ## XZ operator - apply to neigbouring edges of vertices
    XZ = makeXP(0,1,2)
    ## S operator - apply to edges of neighbouring plaquettes, which aren't neigbours
    S = makeXP(0,0,1)
    G += vertex_operator(g,XZ,S)
    return G,N

## apply weighted graph to a cellulation
## return stabilizer generators and precision
def applyWeightedGraph(g,Edges,Weights=None):
    ## default weights if none specified
    Weights = [1]*len(Edges) if Weights is None else Weights
    ## dictionary of edges and weights
    temp = dict()
    for f in g.faceset:
        for i in range(len(Edges)):
            e = tuple(sorted([f[j] for j in Edges[i]]))
            temp[e] = Weights[i]
    E,W = [],[]
    for e,w in temp.items():
        E.append(e)
        W.append(w)
    return E,W


def main():
    m = 4
    torus = False
    o,r,c = (1,1),m,m

    # faces = tri_plane(o,r,c)
    # faces = square_plane(o,r,c)
    # faces = hex_plane(o,r,c)

    ## Union Jack State
    faces = tri_square_plane(o,r,c)
    Edges = [(0,1,2)]
    

    ## 2D Cluster State
    # faces = square_plane(o,r,c)
    # Edges = [(0,1),(1,2),(2,3),(3,0)]    

    g = Cellulation(faces)
    
    if torus:
        g = Torus(g)
    # g.addTail()

    ## Apply weighted graph to g
    # E,W = applyWeightedGraph(g,Edges)
    # GState = graphState(E,W)
    # G,N,n = GState.XPCode(),GState.N,GState.n
    # g.SX2cell(GState.SXx)

    G,N = UnionJackCode(g)
    # G,N = SurfaceCode(g)
    # G,N = SemionCode(g)
    # G,N = TQDZ2Code(g)
    # G,N = SurfaceCode(g,True)


    # print('State')
    # print(State2Str(GState.State(),N))

    print('Generators',len(G))
    # print(XP2Str(G,N))

    C = Code(G,N)

    # print("Symmetry Check")
    # LI = getVal(C,'LI')
    # Symm = C.ZSymmetries()
    
    # for x in Symm:
    #     A = makeXP(0,x,0)
    #     r,u = XPResidual(LI,A,N)
    #     s = str(ZMat2str(x[:n]))
    #     print(s,"is logical identity",isZero(r))
    #     s = s.replace('0',' ')
    #     M = m + (0 if torus else 1)
    #     for i in range(len(s)//M):
    #         print(s[M*i:M*i + M])
    

    # # setVerbose(True)
    S = getVal(C,'S')
    print('Em')
    Em = getVal(C,'Em')
    print(len(Em))
    print('Eq')
    Eq,LXx = cosetDecomposition(Em)
    print('len(Eq)',len(Eq))
    # print(ZmatPrint(Eq))
    print('LXx')
    print(ZmatPrint(LXx))
    # # Eq,LXx = getVals(C,['Eq','LXx'])
    # # print('Eq',len(Eq))
    # # print(ZmatPrint(Eq))
    # # print('LXx',len(LXx))
    
    # # # return False
    # # print('Now for the slow stuff')
    # S,LI,LO,LD,FD = getVals(C,['S','LI','LO','LD','FD'])
    # print('CanonicalGenerators',CanonicalGenerators(G,N,S))
    # print('LI')
    # print(XP2Str(LI,N))

    # print('FD')
    # print(FD)
    # # print('LD')
    # # print(XP2Str(LD,N))

    # print('LI Precision',N,len(LI))
    # # print(XP2Str(LI,N))
    # P = 2*N
    # LI_P = C.getLI(P)
    # print('LI Precision',P,len(LI_P))
    # # print(XP2Str(LI_P,P))

    # print('LO Precision',N,len(LO))
    # for A in LO:
    #     print(XP2Str(A,N),isLO(A,LI,N))

    # LO_P = C.getLO(P)
    # print('LO Precision',P,len(LO_P))
    # for A in LO_P:
    #     print(XP2Str(A,P),isLO(A,LI_P,P))

    #t = squareTorus((0,0),3,3)
    
    #print(semionZ2(3))
    ##s = colour_2d_triangular_SX()
    ##print(s)

    # g.display(show_qubits=False,show_operator=[G[0],G[20],G[13],G[32]])
    for A in G:
        g.display(show_qubits=False,show_operator=A)

# import cProfile
# cProfile.run('main()')
# main()