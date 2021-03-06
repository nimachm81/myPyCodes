## solidstate.py

import numpy
import scipy
from scipy import constants
import math

__all__ = ["FermiDistribution",
            "FerriteMedia"]



def FermiDistribution(E, mu_c, T):
    if T>0.0:
        alpha = (E-mu_c)/(constants.k*T)
        if abs(alpha)<300.0:
            return 1.0/(1.0 + math.exp(alpha))
        elif alpha<=-300.0:
            return 1.0
        else:
            return 0.0
    elif T==0.0:
        if E>mu_c:
            return 0.0
        else:
            return 1.0
    else:
        raise ValueError('negative temperature')
        
        


##------------------ Ferrites ------------


class FerriteMedia:
    
    def __init__(self, B0=None, omega_0=None, M0=None, omega_m=None, alpha=0.0):
        """ B0: magnetic bias
            omega_0: gamma*B0 (if B0 is not provided)
            M0: saturation magnetization
            omega_m: gamma*mu_0*M0
            alpha: loss factor
            only one of M0 and omega_m or one of B0 and omega_0 should be provided
        """
        gyromag_ratio = constants.physical_constants['electron gyromag. ratio']
        print('Gyromagnetic ratio = ', gyromag_ratio)
        self.gamma = gyromag_ratio[0]
        self.omega_0 = None
        if B0!=None:
            self.omega_0 = self.gamma*B0
        if omega_0!=None:
            assert B0==None
            self.omega_0 = omega_0
            print('B_0 = ', omega_0/self.gamma)
        self.omega_m = None
        if M0!=None:
            self.omega_m = self.gamma*constants.mu_0*M0
        if omega_m!=None:
            assert M0==None
            self.omega_m = omega_m
            print('mu_0*M_s = ', omega_m/self.gamma)
        self.alpha = alpha
        return
            
            
    def getRelativePermeabilities(self, omega):
        chi_xx = -self.omega_m*(1j*self.alpha*omega + self.omega_0)/(omega**2 - (1j*self.alpha*omega + self.omega_0)**2)
        chi_xy = -1j*omega*self.omega_m/(omega**2 - (1j*self.alpha*omega + self.omega_0)**2)
        return [1.0+chi_xx, chi_xy]
        
        
    def getRLHandedRelativePermeabilities(self, omega):
        mu_xx, mu_xy = self.getRelativePermeabilities(omega)      
        mu_RH = mu_xx + 1j*mu_xy
        mu_LH = mu_xx - 1j*mu_xy
        return [mu_RH, mu_LH]
        




