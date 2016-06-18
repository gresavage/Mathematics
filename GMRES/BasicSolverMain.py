import sys
import os
# import Gnuplot
import math
from numpy import *
from numpy.oldnumeric.linear_algebra import *
from KSPSolvers import *

# ========================================================================
# Initialize the needed vectors.                         
# ========================================================================
A     = ones(16, float);      A.shape=(4,4)
x     = ones(4, float);      x.shape=(4,1)
b     = ones(4, float);      b.shape=(4,1)
y     = ones(4, float);      y.shape=(4,1)

# ========================================================================
# Populate the System Matrix and RHS.                         
# ========================================================================
#for i in range(len(x)):
#    for j in range(len(x)):
#        A[i,j] = (j+1)+3*(i);
#       b[i,0] = 1; 

A[0,0] = 1.0;
A[0,1] = 2.0;
A[0,2] = 3.0;
A[0,3] = 4.0;

A[1,0] = 8.0;
A[1,1] = 7.0;
A[1,2] = 6.0;
A[1,3] = 5.0;

A[2,0] = 2.0;
A[2,1] = 3.0;
A[2,2] = 5.0;
A[2,3] = 4.0;

A[3,0] = 3.0;
A[3,1] = 4.0;
A[3,2] = 1.0;
A[3,3] = 6.0;

b[0,0] = 3.0;
b[1,0] = 1.0;
b[2,0] = 7.0;
b[3,0] = 9.0;

print "A = ", A
print "b = ", b
# ========================================================================
# Traditional Solver
# ========================================================================
Ainv = linalg.inv(A)
y = Ainv.dot(b)

# ========================================================================
# Our GMRES Solver
# ========================================================================
System = gmres(Matrix=A,RHS=b,x=x,Tol=1e-12,maxIts=len(b))
x, error, totalIters = System.solve()

# ========================================================================
# Pring out some information. 
# ========================================================================

print "GMRES: "
Answer = "Solution: x = [" + str(x[0,0]) + ", " + str(x[1,0]) + ", " + str(x[2,0])+ ", " + str(x[3,0]) + "]^T \n"
print Answer
print x
print A.dot(x)

print "Error"
print error
