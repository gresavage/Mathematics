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

def worldToCode(i, j, grdcellsx):
    return grdcellsx*j+i

def codeToWorld(world, grdcellsx, grdcellsy):
    return world % grdcellsx, world/grdcellsy

def forcing(x, y, t, c):
    return -2.0*c*cos(y)*exp(-t)*sin(x)-cos(y)*exp(-t)*sin(x)

def worldToCodeVector(FN, grdcellsx, grdcellsy):
    fn = zeros((grdcellsx*grdcellsy, 1), float)
    for j in range(grdcellsy):
        for i in range(grdcellsx):
            cdi = worldToCode(i, j, grdcellsx)
            fn[cdi] = FN[i, j]
    return fn

def codeToWorldVector(u, grdcellsx, grdcellsy):
    U = zeros((grdcellsx, grdcellsy), float)
    for k in range(grdcellsx*grdcellsy):
        i, j = codeToWorld(k, grdcellsx, grdcellsy)
        U[i, j] = u[k]
    return U

def customWirePlot(U, X, Y, t):
    wireframefig = plt.figure(figsize=(10.0, 10.0))
    ax = wireframefig.add_subplot(1, 1, 1, projection="3d")
    p = ax.plot_wireframe(X, Y, U, rstride=1, cstride=1)
    strTitle = "Wire Frame Plot of u(x,y,"+str(t)+")"
    wireframefig.suptitle(strTitle)
    plt.xlabel("x-axis")
    plt.ylabel("y-axis")
    plt.show()

def computeTimeStepError(U,t,x0,xEnd,y0,yEnd,grdcellsx,grdcellsy):
    # ==========================================================================*/
    # Set the quad points and weights for the 5 point quadrature rule on [-1,1] */
    # ==========================================================================*/
    quadpt       = zeros(5,float);
    quadpt.shape = (5,1);
    quadwt       = zeros(5,float);
    quadwt.shape = (5,1);

    quadpt[0] = -(1.0/3.0)*((5.0 + 2.0*((10.0/7.0)**0.5))**0.5);
    quadpt[1] = -(1.0/3.0)*((5.0 - 2.0*((10.0/7.0)**0.5))**0.5);
    quadpt[2] =  0.0;
    quadpt[3] =  (1.0/3.0)*((5.0 - 2.0*((10.0/7.0)**0.5))**0.5);
    quadpt[4] =  (1.0/3.0)*((5.0 + 2.0*((10.0/7.0)**0.5))**0.5);

    quadwt[0] = (322.0 - 13.0*((70.0)**0.5))/900.0;
    quadwt[1] = (322.0 + 13.0*((70.0)**0.5))/900.0;
    quadwt[2] = (128.0/225.0);
    quadwt[3] = (322.0 + 13.0*((70.0)**0.5))/900.0;
    quadwt[4] = (322.0 - 13.0*((70.0)**0.5))/900.0;

# Set up a grid to compute on
grdcellsx = 10  # Number of cells in domain x-dir
grdcellsy = 10  # Number of cells in domain y-dir
TimeSteps = 10  # Number of time steps
PlotEvery = 2   # Decide when to plot
c = 1.0         # Diffusion Coefficient


# Set bounds of computational domain
x0 = 0.0
y0 = 0.0
xEnd = 2.0*pi
yEnd = 2.0*pi
t0 = 0.0
T = 2.0

# Mesh Spacings
dx = (xEnd - x0)/grdcellsx
dy = (yEnd - x0)/grdcellsy
dt = (T-t0)/TimeSteps
gammax = (c*dt)/(dx*dx)
gammay = (c*dt)/(dy*dy)

# Position of the cell centers
x = zeros((grdcellsx, 1), float)
y = zeros((grdcellsy, 1), float)
x = linspace(x0+0.5*dx, xEnd-0.5*dx, grdcellsx)
y = linspace(y0+0.5*dy, yEnd-0.5*dy, grdcellsy)
X, Y = meshgrid(x, y)
U0 = zeros((grdcellsx, grdcellsy), float)
# animate = zeros((grdcellsx, grdcellsy, TimeSteps), float)


U0 = utrue(X, Y, t0)

# Set up the Matrix System for Forward Euler.
unknowns = grdcellsx*grdcellsy
Af = zeros((unknowns, unknowns), float)
un = zeros((unknowns, 1), float)
unpo = zeros((unknowns, 1), float)
fn = zeros((unknowns, 1), float)
fnpo = zeros((unknowns, 1), float)
FN = zeros((grdcellsx, grdcellsy), float)
FNPO = zeros((grdcellsx, grdcellsy), float)

for j in range(grdcellsy):
    for i in range(grdcellsx):
        # set entry posiion
        ipo = (i+1) % grdcellsx
        jpo = (j+1) % grdcellsy
        imo = (i-1) % grdcellsx
        jmo = (j-1) % grdcellsy

        cdi = worldToCode(i, j, grdcellsx)
        cdj = worldToCode(i, j, grdcellsy)
        cdipo = worldToCode(ipo, j, grdcellsx)
        cdimo = worldToCode(imo, j, grdcellsx)
        cdjpo = worldToCode(i, jpo, grdcellsy)
        cdjmo = worldToCode(i, jmo, grdcellsy)

        Af[cdi, cdj] = 1.0-2.0*gammax-2.0*gammay
        Af[cdipo, cdj] = gammax
        Af[cdimo, cdj] = gammax
        Af[cdi, cdjpo] = gammay
        Af[cdi, cdjpo] = gammay

# print "Af = ", Af
U0 = utrue(X, Y, t0)
# animate[:, :, 0] = U0
# Wire Frame of Initial Solution
customWirePlot(U0, X, Y, t0)

# Map U(0) --> un
un = worldToCodeVector(U0, grdcellsx, grdcellsy)

# Move forward in loop
for k in range(TimeSteps):
    t = (k+1)*dt
    
    # Sset a forcing function for our code
    FN = forcing(X, Y, (t-dt), c)
    
    # change the chape of our forcing function
    fn = worldToCodeVector(FN, grdcellsx, grdcellsy)
    
    # Advance the solution in time using a Forward Euler Technique
    unpo = Af.dot(un)+dt*fn
    UNPO = codeToWorldVector(unpo, grdcellsx, grdcellsy)

    # animate[:, :, TimeSteps+1] = UNPO
    # Look at the solution if desired
    if((k+1)%PlotEvery == 0):
    # Reshape the solution for plotting
    
    # Plot the solution
        customWirePlot(UNPO, X, Y, t)
    un = unpo
# # Test worldToCode
# for j in range(grdcellsy):
#     for i in range(grdcellsx):
#         print "(i, j) = (", i, ", ", j, ")\n Code = ", worldToCode(i, j, grdcellsx), "\n"
# 
# # Test codeToWorld
# for k in range(grdcellsx*grdcellsy):
#     print "(i, j) = ", codeToWorld(k, grdcellsx, grdcellsy), "\n Code = ", k, "\n"
