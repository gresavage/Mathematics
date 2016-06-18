__author__ = 'Tom'

import sys
import os
import platform
import math
import re
import time
import timeit
import profile
from termcolor import colored
import itertools
import numpy as np
from pylab import *
from collections import deque
import matplotlib.pyplot as plt
from custom_library import *

def weak2norm(x):
    x = np.array(x)
    return sqrt(sum(x**2.))

def strong2norm(x):
    r = sys.float_info.min  # Smallest encodable float (smallest model number)
    R = sys.float_info.max  # Largest encodable float (largest model number)
    b = sys.float_info.radix**np.ceil((sys.float_info.min_exp-1)/2.0)   # Ceil of sqrt of smallest number algorithm can handle
    B = sys.float_info.radix**np.floor((sys.float_info.max_exp+1)/2.0)  # Floor of sqrt of largest number algortithm can handle
    s = sys.float_info.radix**np.floor((sys.float_info.min_exp-1)/2.0)  # Floor of sqrt of smallest number algorithm can handle
    S = sys.float_info.radix**np.round((sys.float_info.max_exp+1)/2.0)  # Rounded sqert of largest number algorithm can handle
    N = sys.float_info.radix**(sys.float_info.mant_dig-1)-1 # Maximum size of x for which accuracy can be guaranteed
    epsilon = np.finfo(float).eps # Machine epsilon

    x = np.array(x)
    a = [0, 0, 0]

    if len(x) > N:
        print "length of input (%r) too large to guarantee accuracy. Must be less than %r" %(len(x), N)

    # Loop to scale and store values of x if they are near machine precision
    for i in range(len(x)):
        if np.abs(x[i]) > B:
            a[2] += (x[i]/S)**2.0
        elif np.abs(x[i]) < b:
            a[0] += a[0] + (x[i]/s)**2.0
        else:
            a[1] += x[i]**2.0

    # Loop to calculate 2-norm avoiding overflow and underflow
    if a[2] != 0:
        if sqrt(a[2]) > R/S:
            raise OverflowError
            return R
        if a[1] != 0:
            ymin = min(sqrt(a[1]), S*sqrt(a[2]))
            ymax = max(sqrt(a[1]), S*sqrt(a[2]))
        else:
            return S*sqrt(a[2])
    elif a[0] != 0:
        if a[1] != 0:
            ymin = min(sqrt(a[1]), s*sqrt(a[0]))
            ymax = max(sqrt(a[1]), s*sqrt(a[0]))
        else:
            return s*sqrt(a[0])
    else:
        return sqrt(a[1])
    if ymin < epsilon*ymax/2.:
        return ymax
    else:
        return ymax*sqrt(1+(ymin/ymax)**2.0)

print "*" * 20, "Machine Precision", "*" * 20
print "Min float: ", sys.float_info.min
print "Max float: ", sys.float_info.max, "\n"
weak2 = []
strong2 = []
error = []

for i in range(1000):
    x = [np.random.random()*np.random.random_integers(0, 1000000000) for j in range(1000)]
    weak2.append(weak2norm(x))
    strong2.append(strong2norm(x))
    error.append(np.abs(weak2[i] - strong2[i])/strong2[i])
error.sort(reverse=True)

# Plotting the Error
fig = plt.figure()
fig.suptitle("Sorted relative error in $||x||_2$ for random $x$")

ax = fig.add_subplot(111)
ax.hold(True)
ax.plot(error, color='r')

ax.annotate(s="Maximum Error: %.4e" % max(error), xy=(error.index(max(error)), max(error)), xytext=(0.125, 0.8), textcoords="axes fraction", arrowprops=dict(arrowstyle="->", connectionstyle="arc"))
ax.annotate(s="Mean Error: %.4e" % mean(error), xy=(500, mean(error)), xytext=(0.5, 0.5), textcoords="axes fraction", arrowprops=dict(arrowstyle="->", connectionstyle="arc"))


ax.axhline(mean(error), linestyle="--", color="c")
ax.axhline(max(error), linestyle="--", color="b")

plt.show()

print "*" * 20, "Code Timing", "*" * 20
weak2time = timeit.timeit('weak2norm(%r)' %x, setup="from __main__ import weak2norm", number=1000)
strong2time = timeit.timeit("strong2norm(%r)" %x, setup="from __main__ import strong2norm", number=1000)
print "Weak time: %r\nStrong time: %r" % (weak2time, strong2time)
print "Speed ratio: %r" % (strong2time/weak2time)