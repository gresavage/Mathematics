import sys
import os
import math
from numpy import *
from pylab import *
from KSPSolvers import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D

##############################################################################
# This is a Projection Method Solver for the Navier Stokes Equations.
# Follows the solver from the 85 Kim and Moin Paper using a cell centered
# finite difference method. Adams Bashforth for the advective transport, and
# Crank-Nicholson Scheme for the diffusion terms.
##############################################################################
def pmod(i,n):
    value = i%n
    if value < 0:
        value = n + value
    return value;

def u1true(x,y,t):
    return -cos(pi*x)*sin(pi*y)*exp(-2.0*t);

def u2true(x,y,t):
    return sin(pi*x)*cos(pi*y)*exp(-2.0*t);

def forcing1(x,y,t,Re):
    # Needs to be checked.
    return (0.5*exp(-4.0*t)*pi*sin(2.0*pi*x) - 2.0*exp(-2.0*t)*pi*pi*cos(pi*x)*sin(pi*y)
             + Re*(-exp(-4.0*t)*pi*cos(pi*x)*cos(pi*y)*cos(pi*y)*sin(pi*x) + 2.0*exp(-2.0*t)*cos(pi*x)*sin(pi*y)
                   - exp(-4.0*t)*pi*cos(pi*x)*sin(pi*x)*sin(pi*y)*sin(pi*y)   )) ;

def forcing2(x,y,t,Re):
    # Needs to be checked.
    return (2.0*exp(-2.0*t)*pi*pi*cos(pi*y)*sin(pi*x) + 0.5*exp(-4.0*t)*pi*sin(2.0*pi*y)
             + Re*(-2.0*exp(-2.0*t)*cos(pi*y)*sin(pi*x) - exp(-4.0*t)*pi*cos(pi*x)*cos(pi*x)*cos(pi*y)*sin(pi*y)
                   - exp(-4.0*t)*pi*cos(pi*y)*sin(pi*x)*sin(pi*x)*sin(pi*y)   )) ;

def WorldToCode(i,j,grdcellsx,grdcellsy):
    return grdcellsx*pmod(j,grdcellsy) + pmod(i,grdcellsx);

def CodeToWorld(k,grdcellsx,grdcellsy):
    i = k%grdcellsx;
    j = k/grdcellsy;
    return i,j;

def WorldToCodeVector(FN,grdcellsx,grdcellsy):
    fn       = zeros(grdcellsx*grdcellsy,float);
    fn.shape = (grdcellsx*grdcellsy,1);

    for j in range(grdcellsy):
        for i in range(grdcellsx):
            cdi = WorldToCode(i,j,grdcellsx,grdcellsy);
            fn[cdi] = FN[i,j];
    return fn;

def CodeToWorldVector(u, grdcellsx,grdcellsy):
    U       = zeros(grdcellsx*grdcellsy,float);
    U.shape = (grdcellsx,grdcellsy);
    for k in range(grdcellsx*grdcellsy):
        i,j = CodeToWorld(k,grdcellsx,grdcellsy);
        U[i,j] = u[k];
    return U;

def Plot(U,X,Y,t,componentName):
    # Plot the a solution as a wireframe.
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(1,1,1,projection='3d')
    p = ax.plot_wireframe(X,Y,U, rstride=1,cstride=1)
    strTitle = "Wire Frame Plot of " +str(componentName) + "(x,y," + str(t) + ")";
    fig.suptitle(strTitle)
    plt.xlabel("$x$-axis")
    plt.ylabel("$y$-axis")
    plt.show()

def ComputeTimeStepError(U1,U2,t,x0,xEnd,y0,yEnd,grdcellsx,grdcellsy):
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


    # ==========================================================================*/
    # Loop over the cell centers and calculate the Integration of a function    */
    # ==========================================================================*/
    IntLocalGrid = 0.0;

    for j in range(grdcellsy):
        for i in range(grdcellsx):

            # ==================================================================*/
            # Get the approximated solution over the given cell                 */
            # ==================================================================*/
            U1h = U1[i,j]; U2h = U2[i,j];

            # ==================================================================*/
            # Get cell corners: cxl, cxr, cyb, cyt assuming fully periodic Grid */
            # ==================================================================*/
            cxl  = x0 + (i  )*((xEnd - x0)/grdcellsx);
            cxr  = x0 + (i+1)*((xEnd - x0)/grdcellsx);
            cyt  = y0 + (j+1)*((yEnd - y0)/grdcellsy);
            cyb  = y0 + (j  )*((yEnd - y0)/grdcellsy);
            detJ = 0.25*(cxr-cxl)*(cyt-cyb);

            IntCell = 0.0;

            # This is a 25 point gaussian quadrature rule.
            for indx in range(5):
                xvalue = 0.5*((cxr-cxl)*quadpt[indx] + (cxl + cxr));
                for indy in range(5):
                    yvalue = 0.5*((cyt-cyb)*quadpt[indy] + (cyt + cyb));
                    IntCell = IntCell + ((u1true(xvalue,yvalue,t) - U1h)**2.0 + (u2true(xvalue,yvalue,t) - U2h)**2.0)*(quadwt[indx]*quadwt[indy]*detJ);

            IntLocalGrid = IntLocalGrid + IntCell;

    return IntLocalGrid;


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
# =-=-=-=-=-                 EDIT CODE PARAMETERS BELOW                 =-=-=-=-=-=- #
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #


# Set up a grid to comput on:

grdcellsx = 8;    # Number of cells in domain x-direction
grdcellsy = 8;    # Number of cells in domain y-direction
TimeSteps = 16;   # Time steps the code will take.
PlotEvery = 4;    # Decide how often to print your work.

# Fluid Parameters:
Re = 1.0;         # The Reynolds number.

# Set bounds on the computational domain (OMEGA).
x0   = 0.0;
y0   = 0.0;
xEnd = 2.0*pi;
yEnd = 2.0*pi;
t0   = 0.0;
T    = 0.002;

# Position of the Cell Centers.
x = zeros(grdcellsx, float);  x.shape=(grdcellsx,1)
y = zeros(grdcellsy, float);  y.shape=(grdcellsy,1)

# Mesh spacings
hx = (xEnd - x0)/grdcellsx
hy = (yEnd - y0)/grdcellsy
dt = (T-t0)/TimeSteps

# Some Derived Constants
gammax = (dt)/(2.0*Re*hx*hx)
gammax = (dt)/(2.0*Re*hy*hy)

if(hx != hy):
    print "Warning !!!! Mesh spacings in x and y are different! "

# Cell center positions.
x = linspace(x0+0.5*hx,xEnd-0.5*hx,grdcellsx)
y = linspace(y0+0.5*hy,yEnd-0.5*hy,grdcellsy)

# Grid of cell centers.
X,Y = meshgrid(x,y)

# Initial contition at time t0.
U10   = u1true(X,Y,t0)
U20   = u2true(X,Y,t0)
U1nmo = u1true(X,Y,t0-dt)
U2nmo = u2true(X,Y,t0-dt)

# Plot the initial solution as a Wire Frame.
Plot(U10,X,Y,t0,"u1")
Plot(U20,X,Y,t0,"u2")

# Set up the Matrix Systems that will be needed. Note we will use Au twice each time step.
unknowns = grdcellsx*grdcellsy
Au    = zeros((unknowns, unknowns), float)
Aphi  = zeros((unknowns, unknowns), float)

# time "n-1"
u1nmo = zeros((unknowns,1),float)
u2nmo = zeros((unknowns,1),float)
# time "n"
u1n   = zeros((unknowns,1),float)
u2n   = zeros((unknowns,1),float)
# time "n+1"
u1npo = zeros((unknowns,1),float)
u2npo = zeros((unknowns,1),float)
phi   = zeros((unknowns,1),float)
# Forcing functions for velocity!
f1nph   = zeros((unknowns,1),float)
f2nph   = zeros((unknowns,1),float)
# RHS vectors for velocity and pressure!
RHSu1   = zeros((unknowns,1),float)
RHSu2   = zeros((unknowns,1),float)
RHSphi  = zeros((unknowns,1),float)
# Intermediate step for velocity!
u1star  = zeros((unknowns,1),float)
u2star  = zeros((unknowns,1),float)


# Populate the system matrix Au1, Au2, and phi.
for j in range (grdcellsy):
    for i in range (grdcellsx):
        # set entry position.
        ipo = (i+1)%grdcellsx
        jpo = (j+1)%grdcellsy
        imo = (i-1)%grdcellsx
        jmo = (j-1)%grdcellsy

        cdi   = WorldToCode(i,j,grdcellsx,grdcellsy)
        cdj   = WorldToCode(i,j,grdcellsx,grdcellsy)
        cdipo = WorldToCode(ipo,j,grdcellsx,grdcellsy)
        cdimo = WorldToCode(imo,j,grdcellsx,grdcellsy)
        cdjpo = WorldToCode(i,jpo,grdcellsx,grdcellsy)
        cdjmo = WorldToCode(i,jmo,grdcellsx,grdcellsy)

        # Set up each of the needed system Matricies.
        # Velocity Matrix For LHS (used twice).
        Au[cdi  ,cdj  ] = 1.0 + 2.0*gammax + 2.0*gammay
        Au[cdipo,cdj  ] = -gammax
        Au[cdimo,cdj  ] = -gammax
        Au[cdi  ,cdjpo] = -gammay
        Au[cdi  ,cdjmo] = -gammay

        # Pressure Solve (used once).
        Aphi[cdi  ,cdj  ] = (-2.0)(1.0/(hx*hx)+1.0/(hy*hy))
        Aphi[cdipo,cdj  ] = 1.0/(hx*hx)
        Aphi[cdimo,cdj  ] = 1.0/(hx*hx)
        Aphi[cdi  ,cdjpo] = 1.0/(hy*hy)
        Aphi[cdi  ,cdjmo] = 1.0/(hy*hy)



# Map U(0) --> un;
u1nmo = WorldToCodeVector(U1nmo,grdcellsx,grdcellsy)
u2nmo = WorldToCodeVector(U2nmo,grdcellsx,grdcellsy)
u1n   = WorldToCodeVector(U10,grdcellsx,grdcellsy)
u2n   = WorldToCodeVector(U20,grdcellsx,grdcellsy)

# Initialize the computational error:
FullError = 0.0

# Move forward in a loop for the time Steps.
for k in range(TimeSteps):
    # Update time step.
    t = (k+1)*dt

    # Set a forcing function for our code.
    F1NPH = forcing1(X,Y,(t-0.5*dt),Re)
    F2NPH = forcing2(X,Y,(t-0.5*dt),Re)

    # Change the shape of our forcing function.
    f1nph = WorldToCodeVector(F1NPH,grdcellsx,grdcellsy)
    f2nph = WorldToCodeVector(F2NPH,grdcellsx,grdcellsy)

    #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
    # STEP ONE: find u*                      #
    #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#

    # Put together RHS of systems.
    #RHSu1 = # Fill in
    #RHSu2 = # Fill in

    # Solve Systems.
    # u1
    System = gmres(Matrix=Au, RHS=RHSu1, x=u1n, Tol=1e-12, maxIts=len(u1n))
    u1star, error, totalIters = System.solve()
    # u2
    System = gmres(Matrix=Au, RHS=RHSu2, x=u2n, Tol=1e-12, maxIts=len(u2n))
    u2star, error, totalIters = System.solve()

    #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
    # STEP TWO: div free pressure            #
    #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#

    # Use phi to find velocity.
    #RHSphi = # Fill in

    # Solve Systems.
    # phi
    System = gmres(Matrix=Aphi, RHS=RHSphi, x=phi, Tol=1e-12, maxIts=len(phi))
    phi, error, totalIters = System.solve()

    # Update systems to back out the velocity at time n+1.
    #u1npo = # Fill in
    #u2npo = # Fill in

    # Reshape the solution for error computation and plotting.
    U1NPO = CodeToWorldVector(u1npo, grdcellsx,grdcellsy)
    U2NPO = CodeToWorldVector(u2npo, grdcellsx,grdcellsy)

    # Compute the error difference between a true solution and our approximations.
    StepError = ComputeTimeStepError(U1NPO,U2NPO,t,x0,xEnd,y0,yEnd,grdcellsx,grdcellsy)*dt
    FullError = FullError + StepError

    # Look at the solution if desired.
    if((k+1)%PlotEvery == 0):
         # Plot the solution.
         Plot(U1NPO,X,Y,t,"u1")
         # Report the Error.
         print "Error = ", FullError**0.5
         print "time  = ", t

    # Lag the found time for the next step.
    u1nmo = u1n
    u2nmo = u2n
    u1n   = u1npo
    u2n   = u2npo


"""
# Test WorldToCode
for j in range (-1,grdcellsy+1):
    for i in range (-1,grdcellsx+1):
        print "i, j =", i,j , "Code = ", WorldToCode(i,j,grdcellsx);

# Test CodeToWorld
for k in range(grdcellsx*grdcellsy):
    i,j = CodeToWorld(k,grdcellsx,grdcellsy);
    print "Code = ", k, "i,j = ", i, j;

"""