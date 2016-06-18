__author__ = 'Tom'

import math

def nCr(n, r):
    return round(math.exp(math.lgamma(n+1)-math.lgamma(n-r+1)-math.lgamma(r+1))).__int__()

def nPr(n, r):
    return round(math.exp(math.lgamma(n+1)-math.lgamma(n-r+1))).__int__()

print nPr(4, 2)

print nCr(61440, 1)

n1  = 2             # Number of internal connections at each input switch
n2  = 2             # Number of internal connections at each output switch
r1  = 2             # Number of input switches
r2  = 2             # Number of output switches
x   = 3             # Maximum number of middle switches used to establish multicast connection
# m   = (x-1)*n1+n2   # Number of middle stage links
m = 2*n1-1

print n1, n2, r1, r2, m
# r = 16
# m = 31
# n = 31
in_combin = 0
mid_combin = 0

total_paths = n1*n2*r1*r2*m
print total_paths
print m
for i in range(r2+1):
    print i
    print float(nCr(total_paths, i))
    print nCr(r2, i)
    print float((n1*r1*m).__pow__(i))
    print float(nPr(n1*r1*m, i))
    print float((n1*r1*m).__pow__(i))/float(nPr(n1*r1*m, i))
    # mid_combin += float(nCr(n1*r1*m*r2*n2, r2-i))
    # mid_combin += float((n1*r1*m).__pow__(i))
    # in_combin += nCr(total_paths, i)/float(mid_combin)
    in_combin += nCr(total_paths, i)/float((n1*r1*m).__pow__(i))
    print "mid_combin: ", mid_combin
    print "in_combin: ", in_combin, '\n'
print m*r2*n2
# for i in range(32):
#     in_combin += nCr(31, i)
#     print i
#
# for i in range(17):
#     mid_combin += nCr(16, i)
#     print i
#
# combin = r*m*mid_combin*(in_combin**2.0)

combin = in_combin

print "Number of configurations: ", float(combin)
fivedays = float(5*24*60*60)
print "Pheasible configurations/Sec: ", combin/fivedays, "tests/sec"

a, b = 4, 2
c = a
print nCr(a, b)
for i in range(a-b, a):
    c *= i
print c
print "r!:", math.factorial(b)
print "n!/r!: ", math.factorial(a)/math.factorial(b)
