__author__ = 'Tom'

import numpy as np
import sys
from numpy import linalg as la
import matplotlib.pyplot as plt
from math import *

def f(x):
    y = np.zeros((3, 1))
    y[0, 0] = 16.*x[0, 0]**4. + 16.*x[1, 0]**4. + x[2, 0]**4. - 16.
    y[1, 0] = x[0, 0]**2. + x[1, 0]**2. + x[2, 0]**2. - 3.
    y[2, 0] = x[0, 0]**3. - x[1, 0]
    return y

def Jf(x):
    J = np.zeros((3, 3))
    J[0, 0] = 64.*x[0]**3.
    J[0, 1] = 64.*x[1]**3.
    J[0, 2] = 4.*x[2]**3.

    J[1, 0] = 2.*x[0]
    J[1, 1] = 2.*x[1]
    J[1, 2] = 2.*x[2]

    J[2, 0] = 3.*x[0]**2.
    J[2, 1] = -1.
    return J

def ForwardElimination(A,b):
    n = len(b[:,0])
    for k in range(len(b[:,0])-1):
        for i in range(k+1,n):
            xmult  = A[i,k]/A[k,k]
            A[i,k] = 0.0 #xmult
            for j in range(k+1,n):
                A[i,j] = A[i,j] - xmult*A[k,j]
            b[i,0] = b[i,0] - xmult*b[k,0]
    return (A,b)

def BackwardSubstitution(A,b):
    n = len(b)
    x = np.zeros((n,1))
    x[n-1,0] = b[n-1,0]/A[n-1,n-1]
    for i in range(n-2,-1,-1):
        s = b[i,0]
    for j in range(i+1,n):
        s = s - A[i,j]*x[j,0]
        x[i,0] = s/A[i,i]
    return x

def GaussianElimination(A,b):
    A,b = ForwardElimination(A, b)
    return BackwardSubstitution(A, b)

def NewtonsMethod(x0, Jf=Jf, f=f, tol=1e-8):
    diff = tol + 1.
    while diff > tol:
        A = Jf(x0)
        b = -f(x0)
        # s = GaussianElimination(A, b)
        s = la.solve(A, b)
        xnew = x0 + s
        diff = la.norm(s)
        x0 = xnew
    return x0

x0 = np.ones((3, 1))

print "f(x0)= ", f(x0)
print "Jf(x0)= ", Jf(x0)
ans = NewtonsMethod(x0)
print "Solution= ", ans
print "Check: ", f(ans)