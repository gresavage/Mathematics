__author__ = 'Tom'

import sys
import os
from math import *
import matplotlib.pyplot as plt

# Use limit definition to calculate d[sin(x)]/dx at x=5

error = [1]
x=0.5
dx = [16]
iter = [0]

while error[-1] > 1e-8:
    deriv = (sin(x+dx[-1])-sin(x))/dx[-1]
    error.append(abs(cos(x)-deriv)/abs(cos(x)))
    dx.append(dx[-1]/2.0)
    iter.append(iter[-1]+1)

print "Limit Approximation: ", deriv
print "Exact value: ", cos(x)
print "dx: %r       Relative Error: %r" %(dx[-1], error[-1])
print len(error)

plt.plot(iter, error)
plt.show()