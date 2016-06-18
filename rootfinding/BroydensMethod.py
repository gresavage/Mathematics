__author__ = 'Tom Gresavage'

import numpy as np
import numpy.linalg as la


def fastJacobiInv(xn, xnmo, fn, fnmo, Jinv):
    """
    Fast prediction of the inverse Jacobian.
    :param xn:      x vector at iteration n
    :param xnmo:    x vector at iteration n-1
    :param fn:      f vector at iteration n
    :param fnmo:    f vector at iteration n-1
    :param Jinv:    inverse Jacobian at iteration n-1
    :return:    fastJacobiInv, the computationally inexpensive estimation of the Jacobian inverse.
    See http://link.springer.com/article/10.1007%2FBF01931297
    """
    dx = xn - xnmo
    df = fn - fnmo
    return Jinv + np.dot(dx - Jinv.dot(df), np.transpose(df)) / la.norm(df)


def shermorJacobiInv(xn, xnmo, fn, fnmo, Jinv):
    """
    Calculates the Sherman-Morrison Jacobi Inverse
    :param xn:      x vector at iteration n
    :param xnmo:    x vector at iteration n-1
    :param fn:      f vector at iteration n
    :param fnmo:    f vector at iteration n-1
    :param Jinv:    inverse Jacobian at iteration n-1
    :return:        shermorJacobiInv, the estimated inverse Jacobian calculated using the Sherman-Morrison method
    """
    dx = xn - xnmo
    df = fn - fnmo
    return Jinv + np.dot(dx - Jinv.dot(df), np.dot(np.transpose(dx), Jinv)) / (np.dot(np.transpose(dx), Jinv.dot(df)))


def updatex(xn, fn, Jinv):
    """
    Updates x using the inverse Jacobian matrix and f(x_n)
    :param xn:      x vector at iteration n
    :param fn:      f vector at iteration n
    :param Jinv:    inverse Jacobian at iteration n
    :return:        xnpo, prediction of x at iteration n+1
    """
    return xn - Jinv.dot(fn)


def f(x):
    """
    An arbitrary 4D vector function of x for testing
    :param x:       a 4D vector
    :return:        f(x), the function evaluated at x
    """
    return np.array([x[0] + x[1] - 2., x[0] * x[2] + x[2] * x[3], x[0] * (x[2] ** 2.) + x[1] * (x[3] ** 2.) - 2. / 3.,
                     x[0] * (x[2] ** 3.) + x[1] * (x[3] ** 3.)])


def Jacobi(x):
    """
    The true Jacobian of f(x)
    :param x:
    :return: J(x), the Jacobian of x
    """
    return np.array(
        [[1, 1, 0, 0], [x[2], x[3], x[0], x[1]], [x[2] ** 2., x[3] ** 2., 2. * x[0] * x[2], 2. * x[1] * x[3]],
         [x[2] ** 3., x[3] ** 3., 3. * x[0] * (x[2] ** 2.), 3. * x[1] * (x[3] ** 2.)]])


# 1 Calculate Jacobi_0 using x_0
# 2 Calculate Jacobi inverse
# 3 Calculate f(x_0)
# 4 Update x_n --> x_n+1
# 5 Calculate Jinv_n+1
# 6 Repeat 2-5 until error is below tolerance

def broydens(x0, f, jacobi, tol=1e-8, limit=1e5):
    # Calculates the solution to a nonlinear function 'f' given intial guess 'x0' using Broyden's root finding method in
    # conjunction with the Sherman-Morrison method of updating the Jacobian inverse.
    # Jacobian of 'f' must be non-singular for 'x0'

    error = [1.1 * tol]
    Jinvgood = la.inv(jacobi)
    xnmogood = x0
    iter = 0

    while error[-1] > tol and iter < limit:
        # Calculate f_n-1
        fnmogood = f(xnmogood)

        # Calculate x_n
        xngood = updatex(xnmogood, fnmogood, Jinvgood)

        # Calculate Jinv_n+1
        Jinvgood = shermorJacobiInv(xngood, xnmogood, f(xngood), fnmogood, Jinvgood)

        # Calculate Error
        error.append(la.norm(xngood - xnmogood))

        # Update variables for next iteration
        xnmogood = xngood
        iter += 1

    if iter >= limit:
        print "No convergence or convergence too slow"
        print "Error at %r iterations: %r" % (iter, error[-1])
        if np.sign(error[-1] - error[-2]) < 0.:
            print "Decreasing error"
            print
        else:
            print "Increasing error"
            print

    # Calculate the convergence rate
    n = len(error)
    A = np.ones((n - 1, 2))
    b = np.ones((n - 1, 1))

    for j in range(n - 1):
        A[j, 0] = np.log(error[j])
        b[j, 0] = np.log(error[j + 1])

    convergence = np.linalg.solve(np.transpose(A).dot(A), np.transpose(A).dot(b))
    print "Convergence Rate: ", convergence[0]
    print
    return xngood, error, iter, convergence


x0 = np.array([[2.], [2.], [2. / 3.], [4.]])
j0 = Jacobi(x0)
x1 = np.array([[2.], [2.], [2. / 3.], [5.]])
j1 = Jacobi(x1)
x2 = np.array([[2.], [2.], [1.], [4.]])
j2 = Jacobi(x2)
x3 = np.array([[20.], [-102.], [-1.], [4.]])
j3 = Jacobi(x3)
x4 = np.array([[-2.], [2.], [2. / 3.], [-4.]])
j4 = Jacobi(x4)

y0, error0, iter0, _ = broydens(x0, f, j0, 1e-10, 1e5)
y1, error1, iter1, _ = broydens(x1, f, j1, 1e-10, 1e5)
y2, error2, iter2, _ = broydens(x2, f, j2, 1e-10, 1e5)
y3, error3, iter3, _ = broydens(x3, f, j3, 1e-10, 1e5)
y4, error4, iter4, _ = broydens(x4, f, j4, 1e-10, 1e5)
print y0
print iter0
print
print y1
print iter1
print error1[-1]
print
print y2
print iter2
print
print y3
print iter3
print
print y4
print iter4
print
print la.norm(f(y0))
print la.norm(f(y1))
print la.norm(f(y2))
print la.norm(f(y3))
print la.norm(f(y4))
