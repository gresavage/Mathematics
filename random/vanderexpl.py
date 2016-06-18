__author__ = 'Tom'

import sys
import os
import numpy as np
from math import *
import matplotlib.pyplot as plt
# from IPython.html import widgets
# from IPython.html.widgets import interact
# from IPython.display import display
# %matplotlib inline

# Let n be the size of our test system

n = 3


# Matrix Multiplication
# print (A*b) # is Cell by Cell
# print A.dot(b)

# Define a function that creates system matrices
def create_A_and_b(n):
    A = np.zeros((n, n), dtype=float)
    b = np.zeros((n, 1), dtype=float)
    for i in range(n):
        # Put entries in b
        b[i, 0] = (1.0/(i+1))*((1.0+(i+1))**n - 1.0)

        # Put entries in A
        for j in range(n):
            A[i, j] = (i+2)**j

    return A, b

A, b = create_A_and_b(5)

print A, b

# Gaussian Elimination
def forward_elim(A, b):
    n = len(b)
    for k in range(n-2):
        for i in range(k+1, n):
            print "i = ", i
            print "k = ", k
            print "A[%r, %r] = %r" % (i, k, A[i, k])
            print "A[%r, %r] = %r" % (k, k, A[k, k])
            xmult = A[i, k]/A[k, k]
            print "xmult = ", xmult
            # A[i, k] = xmult
            for j in range(k, n):
                print "A[%r, %r] = %r" % (i, j, A[i, j])
                print "xmult*A[%r, %r]= %r\n" %(k, j, xmult*A[k, j])

                A[i, j] = A[i, j] - xmult*A[k, j]
            b[i, 0] = b[i, 0] - xmult*b[k, 0]
            print "\n"
    print "A = %r\nb = %r" % (A, b)
    plt.spy(A)
    plt.show()

forward_elim(A, b)