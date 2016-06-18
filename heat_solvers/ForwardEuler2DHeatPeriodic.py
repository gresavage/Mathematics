import sys
import os
import math
from numpy import *
from pylab import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D

##############################################################################
# This is a forward Euler Solver for a 2D heat equation. 
#  u_t - c \Delta u = f  on \Omega = [x0, xEnd]X[y0, yEnd] 
#  using periodic boundary conditions. 
##############################################################################
def pmod(i,n):
    value = i%n
    if value < 0:
        value = n + value
    return value;

def utrue(x,y,t):
    return sin(x)*cos(y)*exp(-2.0*t);
    
def forcing(x,y,t,c):
    return 2.0*(c - 1.0)*sin(x)*cos(y)*exp(-2.0*t);

def forceLong(x,y,t,c):
    xpts = size(x);
    print "xpts = ", xpts
    F       = zeros(xpts,float);
    F.shape = (xpts**0.5,xpts**0.5);

    for j in range(int(xpts**0.5)):
        for i in range(int(xpts**0.5)):
            print "i,j = ", i,j;
            print "  x, y = ", X[i,j], " ", Y[i,j];
            F[i,j] = forcing(X[i,j],Y[i,j],t,c);
    return F; 

def WorldToCode(i,j,grdcellsx,grdcellsy):
    return grdcellsx*pmod(j,grdcellsy) + pmod(i,grdcellsx);
    
def CodeToWorld(k,grdcellsx,grdcellsy):
    i = k%grdcellsx; 
    j = k/grdcellsy;
    return i,j;  
    
def WorldToCodeVector(V,grdcellsx,grdcellsy):
    v       = zeros(grdcellsx*grdcellsy,float);
    v.shape = (grdcellsx*grdcellsy,1);

    for j in range(grdcellsy):
        for i in range(grdcellsx):
            cdi = WorldToCode(i,j,grdcellsx,grdcellsy);
            v[cdi] = V[i,j];
    return v;

def CodeToWorldVector(u, grdcellsx,grdcellsy):     
    U       = zeros(grdcellsx*grdcellsy,float);           
    U.shape = (grdcellsx,grdcellsy);
    for k in range(grdcellsx*grdcellsy):
        i,j = CodeToWorld(k,grdcellsx,grdcellsy);
        U[i,j] = u[k]; 
    return U; 
         
def Plot(U,X,Y,t):
    # Plot the a solution as a wireframe. 
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(1,1,1,projection='3d')
    p = ax.plot_wireframe(X,Y,U, rstride=1,cstride=1)
    strTitle = "Wire Frame Plot of u(x,y," + str(t) + ")"; 
    fig.suptitle(strTitle)
    plt.xlabel("$x$-axis")
    plt.ylabel("$y$-axis")
    plt.show()         

def ComputeTimeStepError(U,t,x0,xEnd,y0,yEnd,grdcellsx,grdcellsy):
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
            Uh = U[i,j];
            
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
                    IntCell = IntCell + ((utrue(xvalue,yvalue,t) - Uh )**2.0)*(quadwt[indx]*quadwt[indy]*detJ); 
                    #print "CELL ", i, ",", j;
                    #print "cell x = ", xvalue, "cell y = ", yvalue
                    #print "Uh = ", Uh
                    #print "quad wtx: ", quadwt[indx_x], "quadptx: ", quadpt[indx_x];
                    #print "quad wty: ", quadwt[indx_y], "quadpty: ", quadpt[indx_y]
        
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
TimeSteps = 8;  # Time steps the code will take.
PlotEvery = 2;   # Decide how often to print your work.
# Diffusion Coefficent "c": 
c = 1.0;

# Set bounds on the computational domain (OMEGA). 
x0   = 0.0; 
y0   = 0.0; 
xEnd = 2.0*pi;
yEnd = 2.0*pi;
t0   = 0.0;
T    = 0.002;

# Position of the Cell Centers.     
x = zeros(grdcellsx, float);  x.shape=(grdcellsx,1); 
y = zeros(grdcellsy, float);  y.shape=(grdcellsy,1); 

# Mesh spacings 
hx = (xEnd - x0)/grdcellsx; 
hy = (yEnd - y0)/grdcellsy; 
dt = (T-t0)/TimeSteps;
gamma = (c*dt)/(hx*hx); 

if(hx != hy): 
    print "Warning !!!! Mesh spacings in x and y are different! "

# Cell center positions.
x = linspace(x0+0.5*hx,xEnd-0.5*hx,grdcellsx); 
y = linspace(y0+0.5*hy,yEnd-0.5*hy,grdcellsy); 

# Grid of cell centers.
Y, X = meshgrid(x,y);

# Initial contition at time t0.
U0 = utrue(X,Y,t0);

# Plot the initial solution as a Wire Frame.
Plot(U0,X,Y,t0);

# Set up the Matrix System for Forward Euler. 
unknowns = grdcellsx*grdcellsy; 
Af   = zeros(unknowns*unknowns, float); Af.shape   = (unknowns,unknowns); 
un   = zeros(unknowns,float);           un.shape   = (unknowns,1); 
unpo = zeros(unknowns,float);           unpo.shape = (unknowns,1); 

# Populate the system matrix Af 
for j in range (grdcellsy): 
    for i in range (grdcellsx): 
        # set entry position. 
        ipo = (i+1)%grdcellsx; 
        jpo = (j+1)%grdcellsy; 
        imo = (i-1)%grdcellsx; 
        jmo = (j-1)%grdcellsy;
        
        cdi   = WorldToCode(i,j,grdcellsx,grdcellsy);
        cdj   = WorldToCode(i,j,grdcellsx,grdcellsy);
        cdipo = WorldToCode(ipo,j,grdcellsx,grdcellsy);
        cdimo = WorldToCode(imo,j,grdcellsx,grdcellsy);
        cdjpo = WorldToCode(i,jpo,grdcellsx,grdcellsy);
        cdjmo = WorldToCode(i,jmo,grdcellsx,grdcellsy);
        
        Af[cdi  ,cdj  ] = 1.0 - 4.0*gamma;
        Af[cdipo,cdj  ] = gamma;
        Af[cdimo,cdj  ] = gamma;
        Af[cdi  ,cdjpo] = gamma;
        Af[cdi  ,cdjmo] = gamma;


# Map U(0) --> un; 
un = WorldToCodeVector(U0,grdcellsx,grdcellsy);

# Initialize the computational error:
FullError = 0.0;

# Move forward in a loop for the time Steps. 
for k in range(TimeSteps):
    t = (k+1)*dt

    # Set a forcing function for our code.  
    FN = forcing(X,Y,(t-dt),c);
    #F2 = forceLong(X,Y,(t-dt),c);
    
    #print "Check ==== "
    #print "  FN = ", FN;
    #print " --- ",
    #print "  F2 = ", F2;
    
    # Change the shape of our forcing function. 
    fn = WorldToCodeVector(FN,grdcellsx,grdcellsy);  
    
    # Advance the solution in time using a Forward Euler Technique.   
    unpo = Af.dot(un) + dt*fn;
    
    # Reshape the solution for error computation and plotting.
    UNPO = CodeToWorldVector(unpo, grdcellsx,grdcellsy);

    # Compute the error difference between a true solution and our approximations.
    StepError = ComputeTimeStepError(UNPO,t,x0,xEnd,y0,yEnd,grdcellsx,grdcellsy)*dt;
    FullError = FullError + StepError;

    # Look at the solution if desired. 
    if((k+1)%PlotEvery == 0):
         # Plot the solution. 
         Plot(UNPO,X,Y,t);
         # Report the Error.
         print "Error = ", FullError**0.5;
         print "time  = ", t;

    # Lag the found time for the next step. 
    un = unpo;

"""
print "X = ", X
print "Y = ", Y

# Test WorldToCode
for j in range (-1,grdcellsy+1):
    for i in range (-1,grdcellsx+1):
        print "i, j =", i,j , "Code = ", WorldToCode(i,j,grdcellsx,grdcellsy);
        
# Test CodeToWorld
for k in range(grdcellsx*grdcellsy):
    i,j = CodeToWorld(k,grdcellsx,grdcellsy); 
    print "Code = ", k, "i,j = ", i, j; 
"""