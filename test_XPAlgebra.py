{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 1,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "List of XP operators A\n",
      "[[6 3 5 0 1 4 1 0 0 0 3 0 4 6 2]\n",
      " [2 3 2 2 5 6 4 6 2 0 4 3 0 0 3]\n",
      " [0 0 6 4 6 7 1 6 1 7 7 4 7 7 7]\n",
      " [3 7 0 0 7 1 0 0 4 2 2 7 1 0 5]\n",
      " [3 0 2 0 7 0 5 7 3 0 2 0 3 5 5]\n",
      " [1 3 4 0 1 3 6 1 3 1 4 0 0 6 1]]\n",
      "Rounded list of XP operators A\n"
     ]
    },
    {
     "ename": "TypeError",
     "evalue": "list indices must be integers or slices, not tuple",
     "output_type": "error",
     "traceback": [
      "\u001b[1;31m---------------------------------------------------------------------------\u001b[0m",
      "\u001b[1;31mTypeError\u001b[0m                                 Traceback (most recent call last)",
      "\u001b[1;32m~\\AppData\\Local\\Temp/ipykernel_17340/2219296163.py\u001b[0m in \u001b[0;36m<module>\u001b[1;34m\u001b[0m\n\u001b[0;32m     11\u001b[0m \u001b[0mA\u001b[0m \u001b[1;33m=\u001b[0m \u001b[0mXPRound\u001b[0m\u001b[1;33m(\u001b[0m\u001b[0mA\u001b[0m\u001b[1;33m,\u001b[0m\u001b[0mN\u001b[0m\u001b[1;33m)\u001b[0m\u001b[1;33m\u001b[0m\u001b[1;33m\u001b[0m\u001b[0m\n\u001b[0;32m     12\u001b[0m \u001b[0mprint\u001b[0m\u001b[1;33m(\u001b[0m\u001b[1;34m'Rounded list of XP operators A'\u001b[0m\u001b[1;33m)\u001b[0m\u001b[1;33m\u001b[0m\u001b[1;33m\u001b[0m\u001b[0m\n\u001b[1;32m---> 13\u001b[1;33m \u001b[0mprint\u001b[0m\u001b[1;33m(\u001b[0m\u001b[0mXP2Str\u001b[0m\u001b[1;33m(\u001b[0m\u001b[0mA\u001b[0m\u001b[1;33m,\u001b[0m\u001b[0mN\u001b[0m\u001b[1;33m)\u001b[0m\u001b[1;33m)\u001b[0m\u001b[1;33m\u001b[0m\u001b[1;33m\u001b[0m\u001b[0m\n\u001b[0m\u001b[0;32m     14\u001b[0m \u001b[0mprint\u001b[0m\u001b[1;33m(\u001b[0m\u001b[1;34m'Rescaled to precision P='\u001b[0m\u001b[1;33m,\u001b[0m\u001b[0mP\u001b[0m\u001b[1;33m)\u001b[0m\u001b[1;33m\u001b[0m\u001b[1;33m\u001b[0m\u001b[0m\n\u001b[0;32m     15\u001b[0m \u001b[0mC\u001b[0m \u001b[1;33m=\u001b[0m \u001b[0mXPSetN\u001b[0m\u001b[1;33m(\u001b[0m\u001b[0mA\u001b[0m\u001b[1;33m,\u001b[0m\u001b[0mN\u001b[0m\u001b[1;33m,\u001b[0m\u001b[0mP\u001b[0m\u001b[1;33m)\u001b[0m\u001b[1;33m\u001b[0m\u001b[1;33m\u001b[0m\u001b[0m\n",
      "\u001b[1;32mc:\\Users\\mark\\Dropbox\\PhD\\XPF package\\XPAlgebra.py\u001b[0m in \u001b[0;36mXP2Str\u001b[1;34m(A, N)\u001b[0m\n\u001b[0;32m     86\u001b[0m     \u001b[0mC\u001b[0m \u001b[1;33m=\u001b[0m \u001b[0mXPcomponents\u001b[0m\u001b[1;33m(\u001b[0m\u001b[0mA\u001b[0m\u001b[1;33m)\u001b[0m\u001b[1;33m\u001b[0m\u001b[1;33m\u001b[0m\u001b[0m\n\u001b[0;32m     87\u001b[0m     \u001b[1;32mfor\u001b[0m \u001b[0mi\u001b[0m \u001b[1;32min\u001b[0m \u001b[0mrange\u001b[0m\u001b[1;33m(\u001b[0m\u001b[0mlen\u001b[0m\u001b[1;33m(\u001b[0m\u001b[0mA\u001b[0m\u001b[1;33m)\u001b[0m\u001b[1;33m)\u001b[0m\u001b[1;33m:\u001b[0m\u001b[1;33m\u001b[0m\u001b[1;33m\u001b[0m\u001b[0m\n\u001b[1;32m---> 88\u001b[1;33m         \u001b[0mp\u001b[0m \u001b[1;33m=\u001b[0m \u001b[0mstr\u001b[0m\u001b[1;33m(\u001b[0m\u001b[0mC\u001b[0m\u001b[1;33m[\u001b[0m\u001b[1;36m0\u001b[0m\u001b[1;33m,\u001b[0m\u001b[0mi\u001b[0m\u001b[1;33m]\u001b[0m\u001b[1;33m)\u001b[0m\u001b[1;33m.\u001b[0m\u001b[0mrjust\u001b[0m\u001b[1;33m(\u001b[0m\u001b[0mN2w\u001b[0m\u001b[1;33m)\u001b[0m\u001b[1;33m\u001b[0m\u001b[1;33m\u001b[0m\u001b[0m\n\u001b[0m\u001b[0;32m     89\u001b[0m         \u001b[0mx\u001b[0m \u001b[1;33m=\u001b[0m \u001b[0msep\u001b[0m\u001b[1;33m.\u001b[0m\u001b[0mjoin\u001b[0m\u001b[1;33m(\u001b[0m\u001b[1;33m[\u001b[0m\u001b[0mstr\u001b[0m\u001b[1;33m(\u001b[0m\u001b[0ma\u001b[0m\u001b[1;33m)\u001b[0m \u001b[1;32mfor\u001b[0m \u001b[0ma\u001b[0m \u001b[1;32min\u001b[0m \u001b[0mC\u001b[0m\u001b[1;33m[\u001b[0m\u001b[1;36m1\u001b[0m\u001b[1;33m,\u001b[0m\u001b[0mi\u001b[0m\u001b[1;33m]\u001b[0m\u001b[1;33m]\u001b[0m\u001b[1;33m)\u001b[0m\u001b[1;33m\u001b[0m\u001b[1;33m\u001b[0m\u001b[0m\n\u001b[0;32m     90\u001b[0m         \u001b[0mz\u001b[0m \u001b[1;33m=\u001b[0m \u001b[0msep\u001b[0m\u001b[1;33m.\u001b[0m\u001b[0mjoin\u001b[0m\u001b[1;33m(\u001b[0m\u001b[1;33m[\u001b[0m\u001b[0mstr\u001b[0m\u001b[1;33m(\u001b[0m\u001b[0ma\u001b[0m\u001b[1;33m)\u001b[0m\u001b[1;33m.\u001b[0m\u001b[0mrjust\u001b[0m\u001b[1;33m(\u001b[0m\u001b[0mNw\u001b[0m\u001b[1;33m)\u001b[0m \u001b[1;32mfor\u001b[0m \u001b[0ma\u001b[0m \u001b[1;32min\u001b[0m \u001b[0mC\u001b[0m\u001b[1;33m[\u001b[0m\u001b[1;36m2\u001b[0m\u001b[1;33m,\u001b[0m\u001b[0mi\u001b[0m\u001b[1;33m]\u001b[0m\u001b[1;33m]\u001b[0m\u001b[1;33m)\u001b[0m\u001b[1;33m\u001b[0m\u001b[1;33m\u001b[0m\u001b[0m\n",
      "\u001b[1;31mTypeError\u001b[0m: list indices must be integers or slices, not tuple"
     ]
    }
   ],
   "source": [
    "from XPAlgebra import *\r\n",
    "from common import *\r\n",
    "n = 7\r\n",
    "m = 10\r\n",
    "N = 8\r\n",
    "P = 16\r\n",
    "Q = 2\r\n",
    "A = XPRandom(n,N,m)\r\n",
    "print('List of Random XP operators A')\r\n",
    "print(XP2Str(A,N))\r\n",
    "print('Ensure all have +1 as Eigenvalue')\r\n",
    "A = XPSetEval(A,N)\r\n",
    "print(XP2Str(A,N))\r\n",
    "print('Rescale to precision P=',P)\r\n",
    "C = XPSetN(A,N,P)\r\n",
    "print(XP2Str(C,P))\r\n",
    "print('Rescaled to precision Q=',Q)\r\n",
    "C = XPSetN(A,N,Q)\r\n",
    "print(C)\r\n",
    "print('Operator Distance')\r\n",
    "print(XPdistance(A))\r\n",
    "print('Operator Diagonal')\r\n",
    "print(XPisDiag(A))\r\n",
    "print('Operator Degree')\r\n",
    "print(XPDegree(A,N))\r\n",
    "print('Operator Fundamental Phase')\r\n",
    "print(XPFundamentalPhase(A,N))\r\n",
    "\r\n",
    "print('Test Inverse: A^-1')\r\n",
    "res = XPInverse(A,N)\r\n",
    "print(XP2Str(res,N))\r\n",
    "print(XPInverse(A,N,res))\r\n",
    "print('Test Square: A^2')\r\n",
    "res= XPSquare(A,N)\r\n",
    "print(XP2Str(res,N))\r\n",
    "print(XPSquare(A,N,res))\r\n",
    "d = 5\r\n",
    "print('Test Power: A^',d)\r\n",
    "res = XPPower(A,N,d)\r\n",
    "print(XP2Str(res,N))\r\n",
    "print(XPPower(A,N,d,res))\r\n",
    "\r\n",
    "B = XPRandom(n,N)\r\n",
    "print('B:',XP2Str(B,N))\r\n",
    "print('Test Multiplication: A*B')\r\n",
    "res = XPMul(A,B,N)\r\n",
    "print(XP2Str(res,N))\r\n",
    "print(XPMul(A,B,N,res))\r\n",
    "\r\n",
    "print('Test Conjugation: ABA^-1')\r\n",
    "res = XPConjugate(A,B,N)\r\n",
    "print(XP2Str(res,N))\r\n",
    "print(XPConjugate(A,B,N,res))\r\n",
    "\r\n",
    "print('Test Commutator: [A,B]')\r\n",
    "res = XPCommutator(A,B,N)\r\n",
    "print(XP2Str(res,N))\r\n",
    "print(XPCommutator(A,B,N,res))"
   ]
  }
 ],
 "metadata": {
  "interpreter": {
   "hash": "93ecb60a143ed30565570a02d52015ed90447d3eff4b9455aa518b6d74e1c511"
  },
  "kernelspec": {
   "display_name": "Python 3.8.7 64-bit",
   "name": "python3"
  },
  "language_info": {
   "name": "python",
   "version": ""
  },
  "orig_nbformat": 4
 },
 "nbformat": 4,
 "nbformat_minor": 2
}