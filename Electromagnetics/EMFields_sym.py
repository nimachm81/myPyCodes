## EMFields_sym.py  
## electromagnetic field, green functions, potentials ...
## symbolic

__all__ = ["CoordSys", "GetScalarGreen2D", "GetScalarGreen3D", "EMFields"]

from sympy import Symbol, I, pi, exp, hankel2

from Electromagnetics.VectorCalculus import *


from enum import Enum
class CoordSys(Enum):
    rectangular = 1
    cylindrical = 2
    spherical = 3
    

def GetScalarGreen2D(rho=Symbol(r'\rho'), k=Symbol('k')):
    ## (del_T^2 + k^2)G= -delta(rho)
    G = 1/(4*I)*hankel2(0, k*rho)
    return G
    
def GetScalarGreen3D(r=Symbol('r'), k=Symbol('k')):
    ## (del^2 + k^2)G= -delta(r)
    G = 1/(4*pi)*exp(-I*k*r)/r
    return G


class EMFields:

    def __init__(self, omega=Symbol(r'\omega'), k=Symbol(r'k'), 
                       eps=Symbol(r'\epsilon'), mu=Symbol(r'\mu'),
                       coord=CoordSys.rectangular):
        self.eps = eps
        self.mu = mu
        self.k = k
        self.omega = omega
        self.coord = coord
        if coord==CoordSys.rectangular:
            self.grad = gradient_r
            self.div = divergence_r
            self.curl = curl_r
        elif coord==CoordSys.cylindrical:
            self.grad = gradient_cy
            self.div = divergence_cy
            self.curl = curl_cy

        return

    def GetE_VecPotElec(self, A):
        E = -I*self.omega*A - I/(self.omega*self.mu*self.eps)*self.grad(self.div(A))
        return E

    def GetH_VecPotElec(self, A):
        H = 1/self.mu*self.curl(A)
        return H

    def GetE_VecPotMag(self, F):
        E = -1/self.eps*self.curl(F)
        return E

    def GetH_VecPotMag(self, F):
        H = -I*self.omega*F - I/(self.omega*self.mu*self.eps)*self.grad(self.div(F))
        return H
        
        
    
        
        






