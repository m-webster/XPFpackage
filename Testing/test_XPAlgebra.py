import add_parent_dir
from XPAlgebra import *
from common import *
n = 7
m = 10
N = 8
P = 16
Q = 2
A = XPRandom(n,N,m)
print('List of Random XP operators A')
print(XP2Str(A,N))
print('Ensure all have +1 as Eigenvalue')
A = XPSetEval(A,N)
print(XP2Str(A,N))
print('Rescale to precision P=',P)
C = XPSetN(A,N,P)
print(XP2Str(C,P))
print('Rescaled to precision Q=',Q)
C = XPSetN(A,N,Q)
print(C)
print('Operator Distance')
print(XPdistance(A))
print('Operator Diagonal')
print(XPisDiag(A))
print('Operator Degree')
print(XPDegree(A,N))
print('Operator Fundamental Phase')
print(XPFundamentalPhase(A,N))

print('Test Inverse: A^-1')
res = XPInverse(A,N)
print(XP2Str(res,N))
print(XPInverse(A,N,res))
print('Test Square: A^2')
res= XPSquare(A,N)
print(XP2Str(res,N))
print(XPSquare(A,N,res))
d = 5
print('Test Power: A^',d)
res = XPPower(A,N,d)
print(XP2Str(res,N))
print(XPPower(A,N,d,res))

B = XPRandom(n,N)
print('B:',XP2Str(B,N))
print('Test Multiplication: A*B')
res = XPMul(A,B,N)
print(XP2Str(res,N))
print(XPMul(A,B,N,res))

print('Test Conjugation: ABA^-1')
res = XPConjugate(A,B,N)
print(XP2Str(res,N))
print(XPConjugate(A,B,N,res))

print('Test Commutator: [A,B]')
res = XPCommutator(A,B,N)
print(XP2Str(res,N))
print(XPCommutator(A,B,N,res))


## 0D,0D,0D
C = [1,2,3]
## 0D,0D,1D
C = [1,2,[3,4]]
# ## 0D,1D,0D
C = [1,[2,3],4]
# ## 1D,0D,0D
C = [[1,2],3,4]
# ## 0D,1D,1D
C = [1,[2,3,4],[5,6]]
# ## 1D,1D,0D
C = [[1,2],[3,4,5],[6]]
# ## 1D,2D,1D
C = [[1,2],[[3,4],[5,6]],[7,8]]

print(makeXP(*C))