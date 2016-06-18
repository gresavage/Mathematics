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
    return exp(-2.0*t)*sin(x)*cos(y)

def worldToCode(i, j, grdcellsx, grdcellsy):
    return grdcellsx*pmod(j, grdcellsy) + pmod(i, grdcellsx)

def pmod(i,n):
    value = i%n
    if value < 0:
        value = n + value
    return value

def codeToWorld(world, grdcellsx, grdcellsy):
    return world % grdcellsx, world/grdcellsy

def forcing(x, y, t, c):
    return cos(y)*exp(-2.0*t)*sin(x)*(2.0*c-1.0)

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

def computeTimeStepError(U, t, x0, xEnd, y0, yEnd, grdcellsx, grdcellsy):
    # ==========================================================================*/
    # Set the quad points and weights for the 5 point quadrature rule on [-1,1] */
    # ==========================================================================*/
    quadpt       = zeros((5, 1), float)
    quadwt       = zeros((5, 1), float)

    quadpt[0] = -(1.0/3.0)*((5.0 + 2.0*((10.0/7.0)**0.5))**0.5)
    quadpt[1] = -(1.0/3.0)*((5.0 - 2.0*((10.0/7.0)**0.5))**0.5)
    quadpt[2] =  0.0
    quadpt[3] =  (1.0/3.0)*((5.0 - 2.0*((10.0/7.0)**0.5))**0.5)
    quadpt[4] =  (1.0/3.0)*((5.0 + 2.0*((10.0/7.0)**0.5))**0.5)

    quadwt[0] = (322.0 - 13.0*((70.0)**0.5))/900.0
    quadwt[1] = (322.0 + 13.0*((70.0)**0.5))/900.0
    quadwt[2] = (128.0/225.0);
    quadwt[3] = (322.0 + 13.0*((70.0)**0.5))/900.0
    quadwt[4] = (322.0 - 13.0*((70.0)**0.5))/900.0
    # ==========================================================================*/
    # Loop over the cell centers and calculate the Integration of a function    */
    # ==========================================================================*/
    IntLocalGrid = 0.0

    for j in range(grdcellsy):
        for i in range(grdcellsx):
            # ==================================================================*/
            # Get the approximated solution over the given cell                 */
            # ==================================================================*/
            Uh = U[i,j]

            # ==================================================================*/
            # Get cell corners: cxl, cxr, cyb, cyt assuming fully periodic Grid */
            # ==================================================================*/
            cxl  = x0 + (i  )*((xEnd - x0)/grdcellsx)
            cxr  = x0 + (i+1)*((xEnd - x0)/grdcellsx)
            cyt  = y0 + (j+1)*((yEnd - y0)/grdcellsy)
            cyb  = y0 + (j  )*((yEnd - y0)/grdcellsy)
            detJ = 0.25*(cxr-cxl)*(cyt-cyb)

            IntCell = 0.0

            # This is a 25 point gaussian quadrature rule.
            for indx in range(5):
                xvalue = 0.5*((cxr-cxl)*quadpt[indx] + (cxl + cxr))
                for indy in range(5):
                    yvalue = 0.5*((cyt-cyb)*quadpt[indy] + (cyt + cyb))
                    IntCell = IntCell + ((utrue(xvalue,yvalue,t) - Uh )**2.0)*(quadwt[indx]*quadwt[indy]*detJ)
            IntLocalGrid = IntLocalGrid + IntCell
    return IntLocalGrid


########################################################################################################################
########################################################################################################################
#--------------------------------------------- Declaration of Variables -----------------------------------------------#
########################################################################################################################
########################################################################################################################

#----------------------------------------------------------------------------------------------------------------------#
# Set up a grid to compute on
#----------------------------------------------------------------------------------------------------------------------#
n = 2**5
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
T = 0.002

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
Af = zeros((unknowns, unknowns), float)
Bf = zeros((unknowns, unknowns), float)
un = zeros((unknowns, 1), float)
unpo = zeros((unknowns, 1), float)
fn = zeros((unknowns, 1), float)
fnpo = zeros((unknowns, 1), float)
FN = zeros((grdcellsx, grdcellsy), float)
FNPO = zeros((grdcellsx, grdcellsy), float)
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
        Af[cdi, cdjpo] = gammay

        Bf[cdi, cdj] = 1.0+2.0*gammax+2.0*gammay
        Bf[cdipo, cdj] = -gammax
        Bf[cdimo, cdj] = -gammax
        Bf[cdi, cdjpo] = -gammay
        Bf[cdi, cdjpo] = -gammay


# print "Af = ", Af
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
    FN = forcing(X, Y, (t-dt), c)

#----------------------------------------------------------------------------------------------------------------------#
# Set a forcing function for our code
#----------------------------------------------------------------------------------------------------------------------#
    fn = worldToCodeVector(FN, grdcellsx, grdcellsy)

#----------------------------------------------------------------------------------------------------------------------#
# Advance the solution using Crank-Nicolson
#----------------------------------------------------------------------------------------------------------------------#
    unpo = Af.dot(un)+dt*fn
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