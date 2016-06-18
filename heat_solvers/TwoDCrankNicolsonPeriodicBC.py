__author__ = 'Tom'
import sys
import os
import matplotlib.pyplot as plt
import pylab
import math
from mpl_toolkits.mplot3d.axes3d import Axes3D
from numpy import *
from KSPSolvers import *

# Solves the 2D differential equation
# u_t - c \Delta u = f on \Omega = [x0, xEnd]X[y0, yEnd]
# Using periodic boundary condiation

def utrue(x, y, t):
    return exp(-t)*sin(x)*cos(y)

def forcing(x, y, t, c):
    return exp(-t)*sin(x)*cos(y)*(2.0*c-1)

def worldToCode(i, j, grdcellsx, grdcellsy):
    return grdcellsx*pmod(j, grdcellsy) + pmod(i, grdcellsx)

def pmod(i, n):
    value = i % n
    if value < 0:
        value += n
    return value

def codeToWorld(world, grdcellsx, grdcellsy):
    return world % grdcellsx, world/grdcellsy

def worldToCodeVector(V, grdcellsx, grdcellsy):
    v = zeros((grdcellsx*grdcellsy, 1), float)
    for j in range(grdcellsy):
        for i in range(grdcellsx):
            cdi = worldToCode(i, j, grdcellsx, grdcellsy)
            v[cdi] = V[i, j]
    return v

def codeToWorldVector(u, grdcellsx, grdcellsy):
    U = zeros((grdcellsx, grdcellsy), float)
    for k in range(grdcellsx*grdcellsy):
        i, j = codeToWorld(k, grdcellsx, grdcellsy)
        U[i, j] = u[k]
    return U

def customWirePlot(X, Y, U, t, label=None, overlay=None):
    if isinstance(X, ndarray) and isinstance(Y, ndarray) and isinstance(U, ndarray):
        wireframefig = plt.figure(figsize=(10.0, 10.0))
        ax = wireframefig.add_subplot(111, projection="3d")
        ax.plot_wireframe(X, Y, U, label="u(x,y,t)")
        strTitle = "Wire Frame Plot of u(x,y,%r)" % t
        wireframefig.suptitle(strTitle)
        plt.xlabel("x-axis")
        plt.ylabel("y-axis")
        if label == None:
            plt.legend()
        else:
            plt.legend(str(label))
        plt.show()
    elif isinstance(U, (list, tuple)):
        try:
            numPlts = min(len(U), len(X), len(Y))
            wireframefig = plt.figure(figsize=(10.0, 10.0))
            strTitle = "Wire Frame Plots of u%r(x,y,%r) - u%r(x,y,%r)" % (1, t, numPlts, t)
            wireframefig.suptitle(strTitle)
            if overlay:
                ax = wireframefig.add_subplot(111, projection="3d")
                plt.xlabel("x-axis")
                plt.ylabel("y-axis")
                for i in range(numPlts):
                    ax.plot_wireframe(X[i], Y[i], U[i], label="u%r(x,y,%r)" % (i, t), color=(random.random(), random.random(), random.random()))
                if label == None:
                    plt.legend()
                else:
                    plt.legend(label)
                plt.show()
            elif not overlay:
                for i in range(numPlts):
                    ax.append(wireframefig.add_subplot(ceil(sqrt(numPlts), floor(sqrt(numPlts))), i, label="u%r(x,y,%r)" % (i, t), projection="3d"))
                    ax[i].plot_wireframe(X[i], Y[i], U[i])
                    if label == None:
                        plt.legend()
                    else:
                        plt.legend(label)
                plt.show()
        except UnboundLocalError:
            print "Overlay cannot be type 'None' if input is of type 'list' or 'tuple'"
            raise

def computeTimeStepError(U,t,x0,xEnd,y0,yEnd,grdcellsx,grdcellsy):
    # ==========================================================================*/
    # Loop over the cell centers and calculate the Integration of a function    */
    # ==========================================================================*/
    IntLocalGrid = 0.0;
    hx = ((xEnd - x0)/grdcellsx);
    hy = ((yEnd - y0)/grdcellsy);
    for j in range(grdcellsy):
        ypt = y0 + (0.5 + j)*hy;
        for i in range(grdcellsx):
            xpt = x0 + (0.5 + i)*hx;
            # ==================================================================*/
            # Get the approximated solution over the given cell                 */
            # ==================================================================*/
            IntCell = ((U[i,j] - utrue(xpt,ypt,t))**2.0)*(hx*hy);
            IntLocalGrid = IntLocalGrid + IntCell;
    return IntLocalGrid;


########################################################################################################################
########################################################################################################################
#--------------------------------------------- Declaration of Variables -----------------------------------------------#
########################################################################################################################
########################################################################################################################

#----------------------------------------------------------------------------------------------------------------------#
# Set up a grid to compute on
#----------------------------------------------------------------------------------------------------------------------#
# n = 2**input("Enter discretization factor: ")
n = 2**6
grdcellsx = n  # Number of cells in domain x-dir
grdcellsy = n  # Number of cells in domain y-dir
TimeSteps = n  # Number of time steps
PlotEvery = TimeSteps/4   # Decide when to plot
c = 1.0         # Diffusion Coefficient
print "n = ", n

#----------------------------------------------------------------------------------------------------------------------#
# Set bounds of computational domain
#----------------------------------------------------------------------------------------------------------------------#
x0 = 0.0
y0 = 0.0
xEnd = 2.0*pi
yEnd = 2.0*pi
t0 = 0.0
T = 2.0

#----------------------------------------------------------------------------------------------------------------------#
# Mesh Spacings
#----------------------------------------------------------------------------------------------------------------------#
dx = (xEnd - x0)/grdcellsx
dy = (yEnd - x0)/grdcellsy
dt = (T-t0)/TimeSteps
gammax = (c*dt)/(2.0*dx*dx)
gammay = (c*dt)/(2.0*dy*dy)

#----------------------------------------------------------------------------------------------------------------------#
# Position of the cell centers
#----------------------------------------------------------------------------------------------------------------------#
x = zeros((grdcellsx, 1), float)
y = zeros((grdcellsy, 1), float)
x = linspace(x0+0.5*dx, xEnd-0.5*dx, grdcellsx)
y = linspace(y0+0.5*dy, yEnd-0.5*dy, grdcellsy)
Y, X = meshgrid(x, y)
U0 = zeros((grdcellsx, grdcellsy), float)

#----------------------------------------------------------------------------------------------------------------------#
# Set up the Matrix System for Forward Euler.
#----------------------------------------------------------------------------------------------------------------------#
FullError = 0.0
unknowns = grdcellsx*grdcellsy
Af = zeros((unknowns, unknowns), float)     # Coefficient Matrix of u_n
Bf = zeros((unknowns, unknowns), float)     # Coefficient Matrix of u_n+1
un = zeros((unknowns, 1), float)            # Solution vector
unpo = zeros((unknowns, 1), float)          # Solution vector at next time step
fnpoh = zeros((unknowns, 1), float)         # Forcing function f_n+1/2: Used instead of f_n, and f_n+1 to save on memory
FNPOH = zeros((grdcellsx, grdcellsy), float)# Forcing function f_n+1/2 matrix
b = zeros((unknowns, 1), float)

#----------------------------------------------------------------------------------------------------------------------#
# Initial Solution and Coefficient Matrices
#----------------------------------------------------------------------------------------------------------------------#
U0 = utrue(X, Y, t0)
for j in range(grdcellsy):
    for i in range(grdcellsx):
        # set entry posiion
        ipo = (i+1) % grdcellsx
        jpo = (j+1) % grdcellsy
        imo = (i-1) % grdcellsx
        jmo = (j-1) % grdcellsy

        cdi = worldToCode(i, j, grdcellsx, grdcellsy)
        cdj = worldToCode(i, j, grdcellsx, grdcellsy)
        cdipo = worldToCode(ipo, j, grdcellsx, grdcellsy)
        cdimo = worldToCode(imo, j, grdcellsx, grdcellsy)
        cdjpo = worldToCode(i, jpo, grdcellsx, grdcellsy)
        cdjmo = worldToCode(i, jmo, grdcellsx, grdcellsy)

        Af[cdi, cdj] = 1.0-2.0*gammax-2.0*gammay
        Af[cdipo, cdj] = gammax
        Af[cdimo, cdj] = gammax
        Af[cdi, cdjpo] = gammay
        Af[cdi, cdjmo] = gammay

        Bf[cdi, cdj]   = 1.0+2.0*gammax+2.0*gammay
        Bf[cdipo, cdj] = -gammax
        Bf[cdimo, cdj] = -gammax
        Bf[cdi, cdjpo] = -gammay
        Bf[cdi, cdjmo] = -gammay

#----------------------------------------------------------------------------------------------------------------------#
# Plot Initial Solution
#----------------------------------------------------------------------------------------------------------------------#
U0 = utrue(X, Y, t0)
# Wire Frame of Initial Solution
customWirePlot(X, Y, U0, t0)

# Map U(0) --> un
un = worldToCodeVector(U0, grdcellsx, grdcellsy)

########################################################################################################################
########################################################################################################################
#------------------------------------------- Crank Nicolson Time-Stepping ---------------------------------------------#
########################################################################################################################
########################################################################################################################

for k in range(TimeSteps):
    t = (k+1)*dt
#----------------------------------------------------------------------------------------------------------------------#
# Set a forcing function for our code
#----------------------------------------------------------------------------------------------------------------------#
    FNPOH = forcing(X, Y, t-dt/2.0, c)

#----------------------------------------------------------------------------------------------------------------------#
# Set a forcing function for our code
#----------------------------------------------------------------------------------------------------------------------#
    fnpoh = worldToCodeVector(FNPOH, grdcellsx, grdcellsy)

#----------------------------------------------------------------------------------------------------------------------#
# Advance the solution using Crank-Nicolson
#----------------------------------------------------------------------------------------------------------------------#
    b = Af.dot(un)+dt*fnpoh
    System = gmres(Matrix=Bf, RHS=b, x=un, Tol=1e-12, maxIts=len(b))
    unpo, error, totalIters = System.solve()
    UNPO = codeToWorldVector(unpo, grdcellsx, grdcellsy)

    StepError = computeTimeStepError(UNPO, t, x0, xEnd, y0, yEnd, grdcellsx, grdcellsy)*dt
    FullError = FullError + StepError
    # Look at the solution if desired
    if((k+1)%PlotEvery == 0):
    # Plot the solution
        customWirePlot([X, X], [Y, Y], [utrue(X, Y, t), UNPO], t, label=["U-True", "U-Numeric"], overlay=True)
        print "error = ", FullError**0.5
        print "time = ", t
    un = unpo