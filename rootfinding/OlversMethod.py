__author__ = 'Tom'

import numpy as np
import matplotlib.pyplot as plt

def polynomial(coefficients, name="f"):
    # Creates a polynomial whose coefficients of increasing orders of x are given by "coefficients"
    def function(x):
        y = 0
        for i in range(len(coefficients)):
            y += coefficients[i]*np.power(x, i)
        return y
    function.__name__ = name
    function.coefficients = coefficients

    return function

def polyderiv(polynomial, d=1):
    # Calculates the "d"th derivative of "polynomial"
    __ = polynomial.coefficients
    for j in range(d):
        coefficients = []
        for i in range(len(__)-1):
            coefficients.append(__[i+1]*(i+1))
        __ = coefficients

    def function(x):
        y = 0
        for i in range(len(coefficients)):
            y += coefficients[i]*np.power(x, i)
        return y
    function.__name__ = polynomial.__name__ + "p"*d
    function.coefficients = coefficients
    return function

def olvers(x0, f, tol=1e-8):
    # Finds a root of "f" starting from "x0" using Olver's method and calculates convergence, error, and number of iterations
    fp = polyderiv(f)
    fpp = polyderiv(f, 2)
    error = []
    dist = 1
    xi = x0
    iterations = 0
    while dist>tol:
        xipo = xi-(f(xi)/fp(xi)) - 0.5*(fpp(xi)/fp(xi))*(f(xi)/fp(xi))**2.
        dist = abs(xipo - xi)
        xi = xipo
        iterations+=1
        if abs(xi-root1)!=0:
            error.append(abs(xi-root1))

    # Calculate the convergence rate
    n = len(error)
    A = np.ones((n-1, 2))
    b = np.ones((n-1, 1))

    for j in range(n-1):
        A[j, 0] = np.log(error[j])
        b[j, 0] = np.log(error[j+1])

    convergence = np.linalg.solve(np.transpose(A).dot(A), np.transpose(A).dot(b))

    print "Root: ", xipo
    print "%r(x): " %f.__name__, f(xipo)
    print "Iterations: ", i
    print "Convergence: ", convergence[0, 0]
    return xipo, convergence[0, 0], iterations, error

coefficients = [-9, 0, 1]
f = polynomial(coefficients)

root1 = (-coefficients[1]+np.sqrt(coefficients[1]**2.-4.*coefficients[0])*coefficients[2])/(2.*coefficients[2])
root2 = (-coefficients[1]-np.sqrt(coefficients[1]**2.-4.*coefficients[0])*coefficients[2])/(2.*coefficients[2])

x = [2.**i for i in range(100)]
y = []

for i in range(100):
    root, convergence, iterations, error =  olvers(x[i], f)
    y.append(convergence)

print y.index(min(y))

plt.plot(x, y, "r")
plt.title("Convergence Rates")
plt.xlim([0, 100])
plt.ylim([np.min(y), np.max(y)])
plt.show()