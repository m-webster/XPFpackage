import matplotlib
import numpy as np
from XPCodes import *
from common import *
from NSpace import *
from XPAlgebra import *
import pylab as pl
from matplotlib import collections  as mc
import matplotlib.pyplot as plt

plotNo = 0

def midpoint(label):
    m = np.reshape(ZMat(label),(-1,2))
    return tuple(np.average(m,axis=0))

def pointDist2(v1,v2):
    return np.sum([(v1[i]-v2[i])**2 for i in range(len(v1))])

def pointAngle(v1,v2):
    cellSize = 3
    x, y = v2[0]-v1[0],v2[1]-v1[1]
    a = np.angle(complex(x,y))
    if abs(x) > cellSize:
        a = np.pi - a 
    if abs(y) > cellSize:
        a = 2*np.pi - a 
    if np.isclose(a,0):
        a = 0
    if a < 0:
        a += 2 * np.pi
    if a > 2 * np.pi:
        a = a - 2 * np.pi
    return a

## for vertex v and cellulation g
## return points at the end of edges incident on v
## sorted by the angle made with respect to v
## takes into account cellulations on torus
def pointAngleSort(g,v):
    vCell = g.cells[0][v]
    vList = set()
    for e in vCell.adj[1]:
        vList.add(e[0])
        vList.add(e[1])
    vList.remove(v)
    vList = list(vList)
    aList = [pointAngle(v,v2) for v2 in vList]
    ix = argsort(aList)
    aList = [aList[i] for i in ix]
    vList = [vList[i] for i in ix]
    return vList,aList

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

## return an oriented faces from list of edges
def edges2face(eList):
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

class Cellulation:
    def __init__(self,faces):
        self.clock = 0
        self.aliases = dict()
        self.faceset = faces
        self.generateCells()
        

    def boundary(self):
        eList = [e for e in self.cells[1] if len(self.cells[1][e].adj[2]) < 2]
        return edges2face(eList)


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
                # x, y = min([x,y]),max([x,y])
                xin = 1 if x in a else 0
                yin = 1 if y in a else 0
                if xin + yin == 1:
                    temp.add((x,y))
        # print(func_name(),'temp',temp)
        return temp

    def generateCells(self):
        ## to determine if two points are nonlocal
        longRange = 3
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
                    j = (i + 1) % len(f)
                    a, b = f[i],f[j]
                    d2 = pointDist2(a,b)
                    ## don't reverse for identified points...not currently used
                    if (d2 <= longRange and a > b) or (d2 > longRange and a < b):
                        (a,b) = (b,a)
                    ## default - sort alphanumerically
                    # if a > b:
                    #     (a,b) = (b,a)
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

    ## remove cell c of dimension d
    def remCell(self,c,d):
        if c not in self.cells[d]:
            return False
        # print(func_name(),'before',self.genus())
        mycell = self.remAdj(c,d)
        for i in range(d,3):
            if i > d:
                for k in mycell.adj[i]:
                    self.remAdj(k,i)
            self.setIndex(i)
        # print(func_name(),'after',self.genus())

    ## remove all adj references to cell c of dimension d
    def remAdj(self,c,d):
        if c not in self.cells[d]:
            return False
        mycell = self.cells[d].pop(c)
        for i in range(3):
            for k in mycell.adj[i]:
                if k in self.cells[i] and c in self.cells[i][k].adj[d]:
                    self.cells[i][k].adj[d].remove(c)
        if d == 2:
            self.faceset.remove(c)
        return mycell


    ## for qubit index for cells of dimension i
    def setIndex(self,i):
        HList = self.cells[i]
        mykeys = list(HList.keys())
        if i == 0:
            ## When i==0, we look at location of vertices and order by this
            midpoints = [HList[k].midpoint for k in mykeys]
        else:
            ## When i > 0, we form tuples of the indices of the vertices making up the face
            midpoints = [tuple(sorted([self.cells[0][v].index for v in self.cells[i][k].adj[0]])) for k in mykeys]
        indices = argsort(midpoints)
        
        for i in range(len(indices)):
            HList[mykeys[indices[i]]].index = i

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
        colList = ['#999999','#e41a1c','#ff7f00','#377eb8','#4daf4a','#984ea3','#ffff33','#a65628','#f781bf']
        n = len(z)
        offset = 0.1
        xCoord,yCoord,cList = [],[],[]
        for i in range(n):
            zi = z[i]
            if zi > 0:
                # (a,b) = midpoint(self.index2cell[i])
                (a,b) = self.qubitLocation[i]
                if Xcomponent:
                    # a,b = a+offset,b-offset
                    a,b = a+offset,b
                xCoord.append(a)
                yCoord.append(b)
                c = colList[0] if Xcomponent else colList[zi] 
                cList.append(c)
        return xCoord,yCoord,cList

    def display(self,show_qubits=False,show_operator=None,save=False,title=None,showaxis=False,dpi=100):
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
            fig.savefig('./images/' + "{0:03d}".format(plotNo) + '.png', transparent=False, dpi=dpi, bbox_inches="tight")
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
        self.qubitLocation = dict()
        iMax = 0
        for d in dims:
            for l,c in self.cells[d].items():
                ix = c.index + iMax
                self.cell2index[l] = ix
                self.index2cell[ix] = l
                self.qubitLocation[ix] = midpoint(l)
            iMax += len(self.cells[d])

    ## after making a hypergraph, relate index back to original cells
    def SX2cell(self,SXx):
        self.index2cell = sorted(self.cells[0].keys())
        m,n = np.shape(SXx)
        ix2cell = dict()
        ix2loc = dict()
        cell2ix1 = dict()
        for i in range(n):
            s = SXx[:,i]
            r = [self.index2cell[j] for j in range(m) if s[j]==1]
            r = tuple(sorted(r)) if len(r) > 1 else r[0]
            ix2cell[i] = r
            ix2loc[i] = midpoint(r)
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
        self.qubitLocation = ix2loc
    
    
    ## Has the effect of slicing off the boundary along dual lattice
    ## Removes vertices and edges on the boundary
    ## Keeps faces and edges sliced in half, but updates adjacency matrix
    def removeBoundary(self,BList):
        BLabels = 'UDLR'
        Bdy = self.boundary()[0]
        segments = UDBoundary(Bdy) + LRBoundary(Bdy)
        for s in BList:
            seg = segments[BLabels.index(s)]
            for i in range(len(seg)):
                if i + 1 < len(seg):
                    e = (seg[i],seg[i+1])
                    if e not in self.cells[1]:
                        (x,y) = e
                        e = (y,x)
                    self.remAdj(e,1)
                self.remAdj(seg[i],0)
        for i in range(2):
            self.setIndex(i)


    ## remove pCount points from iterior
    def removePoints(self,pCount):
        if pCount < 1:
            return True
        bdy = set(self.boundary()[0])
        toRemove = []
        for v,c in self.cells[0].items():
            vBdy = set()
            for f in c.adj[2]:
                vBdy.update(set(f))
            vBdy.remove(v)
            if len(vBdy.intersection(bdy)) == 0:
                toRemove.append(v)
                bdy.update(vBdy)
        i = 0
        for v in toRemove:
            self.remCell(v,0)
            i+=1
            if i >= pCount:
                break
        for i in range(3):
            self.setIndex(i)

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
        for j in range(2*c):
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
    # print('Bdy Segments')
    # print(segments)
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


def XZZX_code(g):
    g.makeCellIndex([0])
    N = 2
    G = []
    n = len(g.cells[0])
    for p in g.faceset:
        A = ZMatZeros(n),ZMatZeros(n)
        for i in range(4):
            # v = g.cells[0][p[i]]
            ix = g.cell2index[p[i]]
            j = i % 2
            A[j][ix] = 1
        G.append(makeXP(0,*A))
    return G,N

def surface_code(g,diagOnly=False):
    g.makeCellIndex([1])
    N = 2
    ## Z operator - apply to edges of vertices
    Z = makeXP(0,0,1)
    G = vertex_operator(g,Z)
    if not diagOnly:
        ## X operator - apply to edges of plaquettes
        X = makeXP(0,1,0)
        G += plaquette_operator(g,X)
    return G,N    

def semion_code(g,diagOnly=False):
    g.makeCellIndex([1])
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

#### UNION JACK CODE ####

## return class of vertex in tri_square cellulation
## 0: SW
## 1: Control0
## 2: NW
## 3: Control1
def UJVertexClasses(g):
    vClasses = [[],[],[],[]]
    r,c = 0,0
    vLast = None
    for i in range(len(g.cells[0])):
        v = g.index2cell[i]
        ## check x coordinate to check if we are in a new column
        if vLast is not None:
            if v[0] > vLast[0]:
                c = (c+1) % 4
                r = 0
        k = (c + 2 * r) % 4
        vClasses[k].append(v)
        vLast = v
        r = 1-r
    return vClasses

## Modified non-diagonal stabilizers for Union Jack State on cellulation g
## Apply alternating edge diagional stabilizers to eliminate Z operators
def UJNonDiag(g):
    n = len(g.cells[0]) + len(g.cells[1])
    G = [0]*len(g.cells[0])
    X = makeXP(0,1,0)
    S = makeXP(0,0,1)
    S3 = makeXP(0,0,3)
    vClasses = UJVertexClasses(g)
    for qClass in range(4):
        for v in vClasses[qClass]:
            vCell = g.cells[0][v]
            Xix = [vCell.index] + [g.cell2index[e] for e in vCell.adj[1]]
            Six = [[],[]]
            ## get points at the end of edges centred on v, ordered by angle they make with v
            vList,angleList = pointAngleSort(g,v)
            if qClass % 2 == 1:
                ## Control qubits - create checkerboard pattern
                SPower = 0 if qClass == 1 else 1
            else:
                ## SW/NW stabilizers - 0 in all cases
                SPower = 0
            vLen = len(vList)
            for i in range(vLen):
                j = (i+1) % vLen
                vi,vj = vList[i],vList[j]
                e = (vi,vj)
                if e not in g.cell2index:
                    e = (vj,vi)
                if e in g.cell2index:
                    Six[SPower].append(g.cell2index[e])
                    SPower = 1-SPower
                else: 
                    Six[SPower].append(g.cell2index[vi])
                    # SPower = 0 if np.isclose(np.cos(angleList[j]),0)  else 1
                    Six[1-SPower].append(g.cell2index[vj])

            gX = applyOperator(n,Xix,X)
            gS = applyOperator(n,Six[0],S)
            gS3 = applyOperator(n,Six[1],S3)
            G[vCell.index] = XPRound(gX + gS + gS3,4)
    return G

## Return Stabilizers for Union Jack code on cellulation g
def union_jack_code(g):
    g.makeCellIndex([0,1])
    N = 4
    Z = makeXP(0,0,2)
    n = len(g.cells[0]) + len(g.cells[1])
    G = UJNonDiag(g)

    Sz = kerIA(XPx(G))
    SZ = makeXP(0,0,Sz*N//2)
    G = np.vstack([G,SZ])

    return G,N    


#### TWISTED QUANTUM DOUBLE CODES ####

## SINGLE LEVEL TWIST ##
def TQDZ2Code(g):
    g.makeCellIndex([1])
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

def omegaPart(w,V,t):
    A = []
    for i in range(3):
        if V[i] == '*':
            if t[w[i]] == 0:
                return [] 
        else:
            A.append((V[i],w[i]))
    return [tuple(A)]

def omega(w,V,t):
    temp = []
    temp.extend(omegaPart(w,V,t))
    tPos = V.index('*')
    if len(temp) > 0 and tPos > 0:
        V[tPos-1] = "*"
        temp.extend(omegaPart(w,V,t))
    return temp

def APhase(w,t):
    VPaths = [[7,2,'*',3,9], [6,0,'*',5,10]]
    temp = []
    for VList in VPaths:
        for i in range(3):
            temp.extend(omega(w,VList[i:i+3],t))
    return temp

def allocatePos(vList,angleList):
    temp = ['*' for i in range(6)]
    pi = np.pi
    for i in range(len(vList)):
        v,a = vList[i],angleList[i]
        if np.isclose(0,a):
            ix = 0
        elif np.isclose(pi,a):
            ix = 3
        elif a < pi/2:
            ix = 1
        elif a < pi:
            ix = 2
        elif a < 3*pi/2:
            ix = 4
        else:
            ix = 5
        temp[ix] = v
    return temp

def applyQTD(g,w):
    toRemove = [k for k,c in g.cells[0].items() if len(c.adj[2]) == 1]
    ## remove points shared by only one face
    # print(func_name(),'toRemove',toRemove)
    for k in toRemove:
        g.remCell(k,0)
    D = max(w) +1
    APhases = []
    offset = np.array([0,0.1])
    eList = sorted(g.cells[1].keys())
    index2cell = [(e,i) for i in range(D) for e in eList]
    cell2index = {index2cell[i]:i for i in range(len(index2cell))}
    ## allow for any missing lattice edges
    for i in range(D):
        cell2index[('*',i)] = -1
    qubitLocation = [midpoint(e) + i * offset for i in range(D) for e in eList]
    n = len(index2cell)
    # print('n',n)
    ## ancilla qubits
    aList = []
    ## phases for each non-diagonal generator
    for i in range(D):
        t = set2Bin(D,[i])
        AP = APhase(w,t)
        for a in AP:
            ## phase involves edges from 2 different levels
            if len(a) == 2 and a[0][1] != a[1][1]:
                aList.append(a)
        APhases.append(AP)
    # print('Ancilla Qubits',aList)
    # print('APhases',APhases)
    ## set of ancilla qubits to be added
    aQubits = set()
    ## iterate through each vertex
    ## list of edges for each vertex
    vLists = []
    for v in g.cells[0].keys():
        ## generate list of vertices around v ordered by angle
        vList,angleList = pointAngleSort(g,v)
        if len(vList) < 6:
            print(func_name(),vList,angleList)
            vList = allocatePos(vList,angleList)
            print(vList)
        
        ## qubits on edges adjacent to the vertex v
        inner = []
        ## qubits on edges on boundary of hexagon centred on v
        outer = []
        vLen = len(vList)
        ## generate list of edges around the vertex in order [0..11] required for the TQD
        for i in range(vLen):
            j = (i+1) % vLen
            vi,vj = vList[i],vList[j]
            e = (v,vi)
            if e not in g.cells[1]:
                e = (vi,v)
            if vi == '*':
                e = '*'
            inner.append(e)
            e = (vj,vi)
            if e not in g.cells[1]:
                e = (vi,vj) 
            if vi == '*' or vj == '*':
                e = '*'
            outer.append(e)  
        eList = inner + outer 
        vLists.append(eList)
        ## generate labels for ancilla qubits to be added
        for (a,b) in aList:
            ## aix, bix are the qubit numbers for the respective edges
            aix = cell2index[(eList[a[0]],a[1])]
            bix = cell2index[(eList[b[0]],b[1])]
            # aix = (eList[a[0]],a[1])
            # bix = (eList[b[0]],b[1])
            if aix > bix:
                aix,bix = bix,aix
            ## only add ancilla qubit if both vertices are in the lattice
            if aix >= 0 and bix >= 0:  
                aQubits.add((aix,bix))
    ## sort the ancilla qubits by label and update cell2index
    aQubits = sorted(aQubits)
    index2cell.extend(aQubits)
    ## for each primary qubit a, Qubit index of any ancilla qubits
    aDict = {a:[] for a in range(n)}
    for i in range(len(aQubits)):
        a,b = aQubits[i]
        ix = i + n
        cell2index[(a,b)] = ix
        for c in [a,b]:
            aDict[c].append(ix)
        loc = (qubitLocation[a] + qubitLocation[b])/2
        qubitLocation.append(loc)
    ## n is the total number of qubits
    # n = len(cell2index)
    n += len(aQubits)
    ## G - list of generators to be returned
    G = []
    N = 4
    X = makeXP(0,1,0)
    S = makeXP(0,0,1)
    Z = makeXP(0,0,2)
    S3 = makeXP(0,0,3)
    ## Non-Diagonal Stabilizers
    for vList in vLists:
        for z in range(D):
            Six = []
            Zix = []
            S3ix = []
            ## apply X operators around the internal edges
            Xix = [cell2index[(e,z)] for e in vList[:6]]
            ## remove any edges not in lattice
            Xix = [x for x in Xix if x >= 0]
            ## associated ancilla Qubits
            Aix = []
            for x in Xix:
                Aix.extend(aDict[x])
            Xix = Xix + Aix
            # print('Xix',Xix)
            for A in APhases[z]:
                ## convert to qubit indices for each edge
                B = sorted([cell2index[(vList[a[0]],a[1])] for a in A])
                if B[0] >= 0:
                    if A in aList:
                        ## 2 edges, different levels
                        ## qubit index of ancilla
                        Aix = cell2index[tuple(B)]
                        ## apply S and X to ancilla
                        # Xix.append(Aix)
                        Six.append(Aix)
                        ## apply S^3 to the other 2 edges in the phase
                        S3ix.extend(B)
                    elif len(A) == 1:
                        ## 1 edge, apply Z
                        Zix.append(B[0])
                    else:
                        ## 2 edges, same level
                        ## find the third edge c of the triangle 
                        a,b = A
                        az = a[1]
                        a,b = vList[a[0]],vList[b[0]]
                        c = tuple(list(set(a).symmetric_difference(set(b))))
                        if c not in g.cells[1]:
                            c = c[1],c[0]
                        # print('a,b,c',a,b,c)
                        ## apply S^3 to the 2 edges in the phase
                        S3ix.extend(B)
                        ## and S to the third edge
                        Six.append(cell2index[(c,az)])
            ## make the generator A
            A = applyOperator(n,Xix,X)
            A += applyOperator(n,Six,S)
            A += applyOperator(n,Zix,Z)
            A += applyOperator(n,S3ix,S3)
            ## make unique vector rep
            A = XPRound(A,N)
            # print(XPdistance(A))
            # print(XP2Str(A,N))
            G.append(A)

    ## Diagonal stabilizers
    ## apply Z to edges around the faces
    for f,c in g.cells[2].items():
        for z in range(D):
            ix = [cell2index[(e,z)] for e in c.adj[1]]
            G.append(applyOperator(n,ix,Z))

    ## apply Z on ancilla qubits, and associated edges
    for e in aQubits:
        ix = [cell2index[e]]
        ix.extend(list(e))
        G.append(applyOperator(n,ix,Z))
    g.qubitLocation = qubitLocation
    g.cell2index = cell2index
    g.index2cell = index2cell
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

