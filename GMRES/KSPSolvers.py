__author__ = "Dr. John Chrispell, IUP"
from numpy import *
# from numpy.oldnumeric.linear_algebra import *
###############################################################################
##  Python implementations of the GMRES. 
##                                  John Chrispell IUP MATH 450/550 (1/28/2015)
##                                                           Modified:   Today
###############################################################################
class Solver_base:
    """ 
    This is the base class for some KSP solvers.  
    
    members:
       bcgstab
    """
    def __init__(self):
        pass

class gmres(Solver_base):
    """
    This is a function that solves Ax = b using the gmres routine. 
    """	
    def __init__(self,Matrix,RHS, x, Tol=1e-11, maxIts=50):
	Solver_base.__init__(self)
	self.Matrix = Matrix
	self.RHS = RHS
	self.x = x
        self.Tol = Tol
	self.maxIts = maxIts

    def backsub(self,U,b):
        n = len(b)
        x = b
        
        j = n-1
        while(j >= 0):
            if(U[j,j] == 0):
                print "Singular Matrix Error!"
            x[j] = b[j]/U[j,j]
            for i in range(j):
                b[i] = b[i] - U[i,j]*x[j]
            j = j-1
        return x
            

    def givapp(self,c,s,vin,k):
        # Apply a sequence of givens rotations. 

        self.k = k 
        if(self.k <= 1):
            self.c = []
            self.c.append(c)
            self.s = [] 
            self.s.append(s)
        else:
            self.c = c
            self.s = s

        self.vrot = vin

        for i in range(self.k):
            w1 = self.c[i]*self.vrot[i] - self.s[i]*self.vrot[i+1]
            w2 = self.s[i]*self.vrot[i] + self.c[i]*self.vrot[i+1]

            self.vrot[i] = w1 
            self.vrot[i+1] = w2

        return self.vrot
        
    def solve(self):
        # Initilizeation 
        n = len(self.RHS)
        h = zeros((self.maxIts+1)*(self.maxIts+1)); h.shape=(self.maxIts+1,self.maxIts+1)
        v = zeros((self.maxIts+1)*(n));             v.shape=(n,self.maxIts+1)
        c = zeros((self.maxIts + 1));               c.shape=((self.maxIts+1),1)
        s = zeros((self.maxIts + 1));               s.shape=((self.maxIts+1),1)
        
         
        if( sqrt(float( (transpose(self.x)).dot(self.x) )) != 0.0): 
            r = self.RHS - ( (self.Matrix).dot(self.x) ) 
        else:
            r = self.RHS 
        
        rho = sqrt(float((transpose(r).dot(r))) )

        g = zeros((self.maxIts + 1),float); g.shape=(self.maxIts+1,1)
        g[0,0] = rho
        errtol = self.Tol*(sqrt(float((transpose(self.RHS).dot(self.RHS)))))
        error = []
        error.append(rho)
        totalIters = 0

        if(rho < errtol):
            return (self.x)

        v[:,0] = (1.0/rho)*r[:,0] 
        beta = rho
        k = 0

        # The GMRES Iteration 
        while((rho > errtol) and (k < self.maxIts)):
            k = k + 1
            v[:,k] = self.Matrix.dot(v[:,(k-1)])
            normav = sqrt( (transpose(v[:,k])).dot(v[:,k]) )

            # Modified Gram-Schmidt
            for j in range(k):
                h[j,k-1] = (transpose(v[:,j])).dot(v[:,k])
                v[:,k] = v[:,k] - h[j,k-1]*(v[:,j])

            h[k,k-1] = sqrt( (transpose(v[:,k])).dot(v[:,k]))
            normav2 = h[k,k-1]

            # Reorthoganlaize (Brown-Hindmarsh condition)
            if( (normav + 0.001*normav2) == normav ):
                for j in range(k):
                    hr = (transpose(v[:,j])).dot(v[:,k])
                    h[j,k-1] = h[j,k-1] + hr
                    v[:,k] = v[:,k] - hr*v[:,j]
                h[k,k-1] = sqrt((transpose(v[:,k])).dot(v[:,k]))
                
            if(h[k,k-1] != 0):
                v[:,k] = v[:,k]/h[k,k-1]

            # Form the new Givens Rotation 
            if(k > 1): 
                h[0:k,k-1] = self.givapp(c=c[0:k-1,0],s=s[0:k-1,0],vin = h[0:k,k-1],k=k-1)

            nu = sqrt( (transpose(h[k-1:k+1,k-1])).dot(h[k-1:k+1,k-1]) )

            if(nu != 0):
                c[k-1,0] = float(h[k-1,k-1]/nu)
                s[k-1,0] = float(-h[k,k-1]/nu)
                h[k-1,k-1] = float(c[k-1]*h[k-1,k-1] - s[k-1]*h[k,k-1])
                if(k < self.maxIts):
                    h[k,k-1] = 0
                g[k-1:k+1,0] = self.givapp(c=c[k-1,0],s=s[k-1,0],vin=g[k-1:k+1,0],k=1)

            # Update the residual 
            rho = abs(g[k,0])
            error.append(float(rho))

        # Compute x 
        y = self.backsub(U=h[0:k,0:k],b=g[0:k,0])
        totalIters = k
        self.x[:,0] = self.x[:,0] + (v[0:n,0:k]).dot(y)

	return self.x, error, totalIters