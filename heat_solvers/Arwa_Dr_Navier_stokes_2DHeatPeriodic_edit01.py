__author__ = 'Arwa'
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
# the code was completed by Arwa Ashi 
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

# def forcing3(x, y, t, c):
#     r = sqrt(x**2.0 + y**2.0)
#     exptr = exp(-t*r)
#     factor1 = (t*x)**2.0*exptr/(r**2.0)+(t*x**2.0*exptr)/r**3.0-(t*exptr)/r
#     factor2 = ((t*y)**2.0*exptr)/r**2.0+(t*y**2.0*exptr)/r**3.0-(t*exptr)/r
#     return -c*((r*factor1-(2.0*t*x**2.0*exp(-t*r))/(x**2.0+y**2.0)+(1.0/r-x**2.0/r**3.0)*exp(-t*r))*sin(x*y)*cos(x*y)+(r*factor2-(2.0*t*y**2.0*exptr)/r**2.0+(1.0/r-y**2.0/r**3.0)*exptr)*sin(x*y)*cos(x*y)+2.0*((y*exptr)/r-t*y*exptr)*(x*cos(x*y)**2.0-x*sin(x*y)**2.0)+2.0*((x*exptr)/r-t*x*exptr)*(y*cos(x*y)**2.0-y*sin(x*y)**2.0)-4.0*x**2.0*r*exptr*sin(x*y)*cos(x*y)-4.0*y**2.0*r*exptr*sin(x*y)*cos(x*y))-r**2.0*exptr*sin(2.0*x*y)/2.0

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
    # Plot the stream line 
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
n = 16
grdcellsx = n;    # Number of cells in domain x-direction
grdcellsy = n;    # Number of cells in domain y-direction
TimeSteps = n;   # Time steps the code will take.
PlotEvery = n/2.0;    # Decide how often to print your work.

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
x = zeros(grdcellsx, float);  x.shape=(grdcellsx,1); 
y = zeros(grdcellsy, float);  y.shape=(grdcellsy,1); 

# Mesh spacings 
hx = (xEnd - x0)/grdcellsx; 
hy = (yEnd - y0)/grdcellsy; 
dt = (T-t0)/TimeSteps;

# Some Derived Constants
gamma = (dt)/(2.0*Re*hx*hx);

if(hx != hy): 
    print "Warning !!!! Mesh spacings in x and y are different! "

# Cell center positions.
x = linspace(x0+0.5*hx,xEnd-0.5*hx,grdcellsx); 
y = linspace(y0+0.5*hy,yEnd-0.5*hy,grdcellsy); 

# Grid of cell centers.
Y, X = meshgrid(x,y);

# Initial contition at time t0.
U10   = u1true(X,Y,t0);
U20   = u2true(X,Y,t0);
# U30   = u3true(X,Y,t0);
U1nmo = u1true(X,Y,t0-dt);
U2nmo = u2true(X,Y,t0-dt);
# U3nmo = u3true(X,Y,t0-dt);

# Plot the initial solution as a Wire Frame.
Plot(U10,X,Y,t0,"u1");
# Plot(U20,X,Y,t0,"u2");
# Plot(U30,X,Y,t0,"u3");

# Set up the Matrix Systems that will be needed. Note we will use Au twice each time step. 
unknowns = grdcellsx*grdcellsy; 
Au    = zeros(unknowns*unknowns, float); Au.shape    = (unknowns,unknowns);
Au1n  = zeros(unknowns*unknowns, float); Au1n.shape  = (unknowns,unknowns);
Au2n  = zeros(unknowns*unknowns, float); Au2n.shape  = (unknowns,unknowns);
# Au3n  = zeros(unknowns*unknowns, float); Au3n.shape  = (unknowns,unknowns);
Aphi  = zeros(unknowns*unknowns, float); Aphi.shape  = (unknowns,unknowns);
Hu1n  = zeros(unknowns*unknowns, float); Hu1n.shape  = (unknowns,unknowns);
Hu2n  = zeros(unknowns*unknowns, float); Hu2n.shape  = (unknowns,unknowns);
# Hu3n  = zeros(unknowns*unknowns, float); Hu3n.shape  = (unknowns,unknowns);
Hu1nmo= zeros(unknowns*unknowns, float);Hu1nmo.shape = (unknowns,unknowns);
Hu2nmo= zeros(unknowns*unknowns, float);Hu2nmo.shape = (unknowns,unknowns);
# Hu3nmo= zeros(unknowns*unknowns, float);Hu3nmo.shape = (unknowns,unknowns);

# code viwe 
# time "n-1"
u1nmo = zeros(unknowns,float);           u1nmo.shape = (unknowns,1);
u2nmo = zeros(unknowns,float);           u2nmo.shape = (unknowns,1);
# u3nmo = zeros(unknowns,float);           u3nmo.shape = (unknowns,1);

# time "n"
Au1nv = zeros(unknowns,float);           Au1nv.shape = (unknowns,1);
u1n   = zeros(unknowns,float);           u1n.shape   = (unknowns,1);
u2n   = zeros(unknowns,float);           u2n.shape   = (unknowns,1);
# u3n   = zeros(unknowns,float);           u3n.shape   = (unknowns,1);

# time "n+1"
u1npo = zeros(unknowns,float);           u1npo.shape = (unknowns,1);
u2npo = zeros(unknowns,float);           u2npo.shape = (unknowns,1);
# u3npo = zeros(unknowns,float);           u2npo.shape = (unknowns,1);
phi   = zeros(unknowns,float);           phi.shape   = (unknowns,1);

# Forcing functions for velocity!
f1nph   = zeros(unknowns,float);         f1nph.shape = (unknowns,1);
f2nph   = zeros(unknowns,float);         f2nph.shape = (unknowns,1);
# f3nph   = zeros(unknowns,float);         f2nph.shape = (unknowns,1);

# RHS vectors for velocity and pressure!
RHSu1   = zeros(unknowns,float);         RHSu1.shape = (unknowns,1);
RHSu2   = zeros(unknowns,float);         RHSu2.shape = (unknowns,1);
# RHSu3   = zeros(unknowns,float);         RHSu2.shape = (unknowns,1);
RHSphi  = zeros(unknowns,float);         RHSphi.shape= (unknowns,1);

# Intermediate step for velocity!
u1star  = zeros(unknowns,float);         u1star.shape= (unknowns,1);
u2star  = zeros(unknowns,float);         u2star.shape= (unknowns,1);
# u2star  = zeros(unknowns,float);         u2star.shape= (unknowns,1);


# Populate the system matrix Au1, Au2, and phi. 
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
        
        # Set up each of the needed system Matricies.
        # Velocity Matrix For LHS (used twice). 
        Au[cdi  ,cdj  ] = 1.0 + 4.0*gamma;
        Au[cdipo,cdj  ] = -gamma;
        Au[cdimo,cdj  ] = -gamma;
        Au[cdi  ,cdjpo] = -gamma;
        Au[cdi  ,cdjmo] = -gamma;
                
        Au1n[cdi  ,cdj  ] = 1.0 - (4.0/(hx*hx));
        Au1n[cdipo,cdj  ] = 1.0/(hx*hx);
        Au1n[cdimo,cdj  ] = 1.0/(hx*hx);
        Au1n[cdi  ,cdjpo] = 1.0/(hx*hx);
        Au1n[cdi  ,cdjmo] = 1.0/(hx*hx);

        # Au3n[cdi  ,cdj  ] = 1.0 - (4.0/(hx*hx));
        # Au3n[cdipo,cdj  ] = 1.0/(hx*hx);
        # Au3n[cdimo,cdj  ] = 1.0/(hx*hx);
        # Au3n[cdi  ,cdjpo] = 1.0/(hx*hx);
        # Au3n[cdi  ,cdjmo] = 1.0/(hx*hx);
        
        Au2n[cdi  ,cdj  ] = 1.0 - 4.0/(hx*hx);
        Au2n[cdipo,cdj  ] = 1.0/(hx*hx);
        Au2n[cdimo,cdj  ] = 1.0/(hx*hx);
        Au2n[cdi  ,cdjpo] = 1.0/(hx*hx);
        Au2n[cdi  ,cdjmo] = 1.0/(hx*hx);
        
        
        # Pressure Solve (used once).
        Aphi[cdi  ,cdj  ] = (-4.0)/(hx*hx);
        Aphi[cdipo,cdj  ] = 1.0/(hx*hx);
        Aphi[cdimo,cdj  ] = 1.0/(hx*hx);
        Aphi[cdi  ,cdjpo] = 1.0/(hx*hx);
        Aphi[cdi  ,cdjmo] = 1.0/(hx*hx);
        
# Map U(0) --> un;
u1nmo = WorldToCodeVector(U1nmo,grdcellsx,grdcellsy);
u2nmo = WorldToCodeVector(U2nmo,grdcellsx,grdcellsy);
# u3nmo = WorldToCodeVector(U3nmo,grdcellsx,grdcellsy);
u1n   = WorldToCodeVector(U10,grdcellsx,grdcellsy);
u2n   = WorldToCodeVector(U20,grdcellsx,grdcellsy);
# u3n   = WorldToCodeVector(U30,grdcellsx,grdcellsy);

# Initialize the computational error:
FullError = 0.0;

# Move forward in a loop for the time Steps. 
for k in range(TimeSteps):
    # Update time step.
    t = (k+1)*dt
    for j in range (0,grdcellsy,grdcellsy-1): 
     for i in range (grdcellsx-1): 
         # ====================================================================
         # Note that the alpha values are now dependent on the solution
         # ====================================================================
         #u1n    = u1true(X[i,j],Y[i,j],t);
         #u2n    = u2true(X[i,j],Y[i,j],t);
         alpha  = u1true(X[i, j], Y[i, j], t)/(2.0*hx)
         alphanmo = u1true(X[i, j], Y[i, j], t-dt)/(2.0*hx)
         alpha2  = u2true(X[i, j], Y[i, j], t)/(2.0*hx)
         alpha2nmo = u2true(X[i, j], Y[i, j], t-dt)/(2.0*hx)

         # alpha  = u1true(X[i, j], Y[i, j], t)*(u1true(X[i+1,j],Y[i+1,j],t)-u1true(X[i-1,j],Y[i-1,j],t))/(2.0*hx)+u2true(X[i, j], Y[i, j], t)*(u1true(X[i+1,j],Y[i+1,j],t)-u1true(X[i-1,j],Y[i-1,j],t))/(2.0*hx)
         # alphanmo = u1true(X[i, j], Y[i, j], t-dt)*(u1true(X[i+1,j],Y[i+1,j],t-dt)-u1true(X[i-1,j],Y[i-1,j],t-dt))/(2.0*hx)+u2true(X[i, j], Y[i, j], t-dt)*(u1true(X[i+1,j],Y[i+1,j],t-dt)-u1true(X[i-1,j],Y[i-1,j],t-dt))/(2.0*hx)
         # alpha2  = u1true(X[i, j], Y[i, j], t)*(u2true(X[i+1,j],Y[i+1,j],t)-u2true(X[i-1,j],Y[i-1,j],t))/(2.0*hx)+u2true(X[i, j], Y[i, j], t)*(u2true(X[i+1,j],Y[i+1,j],t)-u2true(X[i-1,j],Y[i-1,j],t))/(2.0*hx)
         # alpha2nmo = u1true(X[i, j], Y[i, j], t-dt)*(u2true(X[i+1,j],Y[i+1,j],t-dt)-u2true(X[i-1,j],Y[i-1,j],t-dt))/(2.0*hx)+u2true(X[i, j], Y[i, j], t-dt)*(u2true(X[i+1,j],Y[i+1,j],t-dt)-u2true(X[i-1,j],Y[i-1,j],t-dt))/(2.0*hx)
         # alpha3 = (u3true(X[i,j],Y[i,j],t)*dt)/(2.0*hx)
         # Hu1n[cdi,cdj]   = alpha;  #u1n[cdi]*(1.0/(2*hx));
         Hu1n[cdipo,cdj]   = alpha;  #u1n[cdi]*(1.0/(2*hx));
         Hu1n[cdimo,cdj]   = -alpha; #u1n[cdi]*(-1.0/(2*hx));
         Hu1n[cdi,cdjpo]   = alpha2; #u2n[cdi]*(1.0/(2*hx));
         Hu1n[cdi,cdjmo]   = -alpha2;#u2n[cdi]*(-1.0/(2*hx));

         # Hu2n[cdi,cdj]   = alpha2;  #u1n[cdi]*(1.0/(2*hx));
         Hu2n[cdipo,cdj]   = alpha;  #u1n[cdi]*(1.0/(2*hx));
         Hu2n[cdimo,cdj]   = -alpha; #u1n[cdi]*(-1.0/(2*hx));
         Hu2n[cdi,cdjpo]   = alpha2; #u2n[cdi]*(1.0/(2*hx));
         Hu2n[cdi,cdjmo]   = -alpha2;#u2n[cdi]*(-1.0/(2*hx));

         # Hu1nmo[cdi,cdj] = alphanmo;  #u1n[cdi]*(1.0/(2*hx));
         Hu1nmo[cdipo,cdj] = alphanmo;  #u1n[cdi]*(1.0/(2*hx));
         Hu1nmo[cdimo,cdj] = -alphanmo; #u1n[cdi]*(-1.0/(2*hx));
         Hu1nmo[cdi,cdjpo] = alpha2nmo; #u2n[cdi]*(1.0/(2*hx));
         Hu1nmo[cdi,cdjmo] = -alpha2nmo;#u2n[cdi]*(-1.0/(2*hx));
        
         # Hu2nmo[cdi,cdj] = alpha2nmo;  #u1n[cdi]*(1.0/(2*hx));
         Hu2nmo[cdipo,cdj] = alphanmo;  #u1n[cdi]*(1.0/(2*hx));
         Hu2nmo[cdimo,cdj] = -alphanmo; #u1n[cdi]*(-1.0/(2*hx));
         Hu2nmo[cdi,cdjpo] = alpha2nmo; #u2n[cdi]*(1.0/(2*hx));
         Hu2nmo[cdi,cdjmo] = -alpha2nmo;#u2n[cdi]*(-1.0/(2*hx));
        
    # Set a forcing function for our code.  
    F1NPH = forcing1(X,Y,(t+0.5*dt),Re);
    F2NPH = forcing2(X,Y,(t+0.5*dt),Re);
    
    # Change the shape of our forcing function. 
    f1nph = WorldToCodeVector(F1NPH,grdcellsx,grdcellsy);
    f2nph = WorldToCodeVector(F2NPH,grdcellsx,grdcellsy);
    
    #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
    # STEP ONE: find u*                      #
    #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
    
    #Au1nv  = Au1n.dot(u1n);
    
    # Put together RHS of systems.
    RHSu1 = (dt)/(2.0*Re)*Au1n.dot(u1n)+(1/2*Re)*(3.0*Hu1n.dot(u1n)-Hu1nmo.dot(u1n))+(dt/Re)*f1nph;
    # RHSu1 = (dt)/(2.0*Re)*Au1n.dot(u1n)+(1/2*Re)*(3.0*(Hu1n.dot(u1n)+Hu1n.dot(u2n))-Hu1nmo.dot(u1n)-Hu1nmo.dot(u2n))+(dt/Re)*f1nph;
    
    RHSu2 = (dt)/(2.0*Re)*Au2n.dot(u2n)+(1/2*Re)*(3.0*Hu2n.dot(u2n)-Hu2nmo.dot(u2n))+(dt/Re)*f2nph;
    # RHSu1 = (dt)/(2.0*Re)*Au1n.dot(u2n)+(1/2*Re)*(3.0*(Hu1n.dot(u1n)+Hu1n.dot(u2n))-Hu1nmo.dot(u1n)-Hu1nmo.dot(u2n))+(dt/Re)*f1nph;

    # Solve Systems.
    # u1
    System = gmres(Matrix=Au, RHS=RHSu1, x=u1n, Tol=1e-12, maxIts=len(u1n))
    u1star, error, totalIters = System.solve();
    # u2
    System = gmres(Matrix=Au, RHS=RHSu2, x=u2n, Tol=1e-12, maxIts=len(u2n))
    u2star, error, totalIters = System.solve();

    #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
    # STEP TWO: div free pressure            #
    #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#

    # Use phi to find velocity.
    RHSphi = (u1star+u2star)/dt;
    # RHSphi = u1star/dt
    
    # Solve Systems.
    # phi
    System = gmres(Matrix=Aphi, RHS=RHSphi, x=phi, Tol=1e-12, maxIts=len(phi))
    phi, error, totalIters = System.solve();
    
    # Update systems to back out the velocity at time n+1.
    u1npo = u1star - dt*phi;
    u2npo = u2star - dt*phi;
    
    # Reshape the solution for error computation and plotting.
    U1NPO = CodeToWorldVector(u1npo, grdcellsx,grdcellsy);
    U2NPO = CodeToWorldVector(u2npo, grdcellsx,grdcellsy);
    
    # Compute the error difference between a true solution and our approximations.
    StepError = ComputeTimeStepError(U1NPO, U2NPO, t, x0, xEnd, y0, yEnd, grdcellsx, grdcellsy)*dt;
    FullError = FullError + StepError;

    # Look at the solution if desired. 
    if((k+1)%PlotEvery == 0):
         # Plot the solution. 
         Plot(U1NPO,X,Y,t,"u1");
         # Report the Error.
         print "Error = ", FullError**0.5;
         print "time  = ", t;

    # Lag the found time for the next step.
    u1nmo = u1n;
    u2nmo = u2n;
    u1n   = u1npo;
    u2n   = u2npo;


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




