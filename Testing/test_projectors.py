import add_parent_dir
from XPAlgebra import *
from common import *
print('log2ceil(128)',logCeil(128))
n = 5
m = 1 << (n-2)
N = 4
P = 16
Q = 2

S1,c = StateRandom(N,n,m)
print('Random State S1')
print(State2Str(S1,N))
res = State2C(S1,N)
# print(res)
print('Test conversion to complex nparray',State2C(S1,N,C=res))
a = stateAmplitude(S1,N)
print('Test amplitude',stateAmplitude(S1,N,C=a))

res = State2C(S1,N,c)
print('Test conversion to complex nparray with weights',State2C(S1,N,c,C=res))
a = stateAmplitude(S1,N,c)
print('Test amplitude with weights',stateAmplitude(S1,N,c,C=a))

S2,c = StateRandom(N,n,m)
print('Random State S2')
print(State2Str(S2,N))
res = State2C(S2,N)
print('Test conversion to complex nparray',State2C(S2,N,C=res))

S12,c12 = StateAdd(S1,S2,N)
print('(S1 + S2)')
print(State2Str(S12,2*N,c12))
print('Checking (S1 + S2)/2',StateAdd(S1,S2,N,(S12,c12)))

A = XPRandom(n,N)
print('Measure A')
print(XP2Str(A,N))
print(XP2StrAlt(A,N))
Evals = XPEigenvalues(A,N)
print('Eigenvalues of A',Evals)
l = np.random.choice(Evals)
# l = 0
print('Measure Eigenvalue',l)
A[-1] = np.mod(A[-1] - l,2*N)
print('Adjusted operator A',XP2Str(A,N))
Evals = XPEigenvalues(A,N)
print('Eigenvalues of A',Evals)
res = XPProj(A,S1,N)
# print('Result',res)
Sr,cr = res
print('Projector Check',XPProj(A,S1,N,res))
# print(State2Str(Sr,2*N,cr))

# n=3
# T = StatePlus(N,n,0)
# print(State2Str(T,N))
# A = intMat([2,1,1,1,1,2,3])
# print('Measure A')
# print(XP2Str(A,N))
# print(XP2StrAlt(A,N))
# Evals = XPEigenvalues(A,N)
# print('Eigenvalues of A',Evals)
# l = np.random.choice(Evals)
# # l = 0
# print('Measure Eigenvalue',l)
# A = XPSetEval(A,N,l)
# print('Adjusted operator A',XP2Str(A,N))
# res = XPProj(A,T,N)
# # print('Result',res)
# Sr,cr = res
# print('Projector Check',XPProj(A,T,N,res))
# print(State2Str(Sr,2*N,cr))
# print('Projected State',State2Str(S,N,c))

# A = intMat([[8,0,0,0,0,0,0,0,6,5,5,4,4,4,4],[7,1,1,1,1,1,1,1,1,2,4,1,2,3,4],[1,1,1,1,0,0,0,0,3,1,3,4,4,4,4]])
# print(XP2StrAlt(A,8))