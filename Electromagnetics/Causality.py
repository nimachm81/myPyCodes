## Causality.py
## Causality relations, Hilbert transforms, Kramers Kronig relations...


__all__ = ["KKIntegralR", "KKIntegralI",
            "DrudeModel", "LorentzModel"]

import numpy as np
import scipy as sp
import math
import numpy
from scipy.integrate import quad, quadrature
from scipy.fftpack import fft, ifft, fftshift, ifftshift

from scipy import constants


def KKIntegralR__(X2, omega):
    """ Takes the imaginary part and returns the real part
        X2 is a function representing the imaginary part
    """
    if isinstance(omega, float):
        w = omega
        def f(w_p):
            return w_p*(X2(w_p))/(w_p**2 - w**2)
        eps = w*1.0e-5
        X1 = -2.0/math.pi*(quad(f, 0.0, w-eps)[0] + quad(f, w+eps, 10.0*w)[0]
            + quad(f, 10.0*w, 100.0*w)[0])
        return X1
    else:
        assert isinstance(omega, numpy.ndarray)
        N = len(omega)
        X1 = np.zeros(N, dtype=float)
        for i in range(N):
            X1[i] = KKIntegralR(X2, omega[i])
        return X1
    
        
def KKIntegralI__(X1, omega):
    """ Takes the real part and returns the imaginary part
        X1 is a function representing the real part
    """
    if isinstance(omega, float):
        w = omega
        def f(w_p):
            return X1(w_p)/(w_p**2 - w**2)
        
        eps = w*1.0e-5
        X2 = 2.0/math.pi*w*(quad(f, 0.0, w-eps)[0] + quad(f, w+eps, 10.0*w)[0]
            + quad(f, 10.0*w, np.inf)[0])
        return X2
    else:
        assert isinstance(omega, numpy.ndarray)
        N = len(omega)
        X2 = np.zeros(N, dtype=float)
        for i in range(N):
            X2[i] = KKIntegralI(X1, omega[i])
        return X2
       
def KKIntegralR(X2, omega):
    """ Takes the imaginary part and returns the real part
        X2 is a function representing the imaginary part
    """
    if isinstance(omega, float):
        w = omega
        def f(w_p):
            return (X2(w_p)-X2(w))/(w_p - w)
        eps = w*1.0e-10
        X1 = -1.0/math.pi*(quad(f, -np.inf, 0.0)[0] + quad(f, 0.0, w-eps)[0] 
            + quad(f, w+eps, 10.0*w)[0] + quad(f, 10.0*w, np.inf)[0])
        return X1
    else:
        assert isinstance(omega, numpy.ndarray)
        N = len(omega)
        X1 = np.zeros(N, dtype=float)
        for i in range(N):
            X1[i] = KKIntegralR(X2, omega[i])
        return X1
       
def KKIntegralI(X1, omega):
    """ Takes the real part and returns the imaginary part
        X1 is a function representing the real part
    """
    if isinstance(omega, float):
        w = omega
        def f(w_p):
            return (X1(w_p)-X1(w))/(w_p - w)
        
        eps = w*1.0e-10
        X2 = 1.0/math.pi*(quad(f, -np.inf, 0.0)[0] + quad(f, 0.0, w-eps)[0] 
            + quad(f, w+eps, 10.0*w)[0] + quad(f, 10.0*w, np.inf)[0])
        return X2
    else:
        assert isinstance(omega, numpy.ndarray)
        N = len(omega)
        X2 = np.zeros(N, dtype=float)
        for i in range(N):
            X2[i] = KKIntegralI(X1, omega[i])
        return X2


def HilbertTransform(X, N, f_max):
    ## TODO: test
    """ X: function  
        N: number of freq points   
        omega_max: frequency where X(f_max) is essentially zero
        Hilbert transform as
        H = IF(-j*sign(f)*x(f))       IF-->inverse Fourier
    """
    f = np.linspace(-f_max, f_max, 2*N+1)
    f = f[0:2*N]
    x = X(2.0*math.pi*f)
    H_inv = -1j*(f>0.0 - f<0.0)*x
    H = ifft(ifftshift(H_inv))/(2*N)
    return H

from scipy import constants


class DrudeModel:

    def __init__(self, tau, sigma_0=None, n_e=None):
        e = constants.e
        m = constants.m_e
        self.tau = tau
        if sigma_0!=None:
            self.sigma_0 = sigma_0
        elif n_e!=None:
            self.sigma_0 = tau*e**2*n_e/m
        else:
            raise ValueError('at least one of sigma_0 or ne is required')
        return
        
        
    def GetConductivity(self, omega):
        return self.sigma_0/(1j*omega*self.tau + 1)
        
        
        
class LorentzModel:
    def __init__(self, gamma, omega_0=None, chi_0=None, k=None, n_e=None):
        """ n_e: density of atoms
        """
        e = constants.e
        m = constants.m_e
        self.gamma = gamma
        if omega_0!=None:
            self.omega_0 = omega_0
        elif k!=None:
            self.omega_0 = math.sqrt(k/m)
        else:
            raise ValueError('at least one of omega_0 or k is required')
        if chi_0!=None:
            self.chi_0 = chi_0
        elif n_e!=None:
            self.chi_0 = (e**2*n_e/(m*self.omega_0**2))
        else:
            raise ValueError('at least one of chi_0 or ne is required')
        return

        
    def GetSusciptibility(self, omega):
        return self.chi_0/(1j*self.gamma*omega/self.omega_0**2 - omega**2/self.omega_0**2 + 1.0)

    def GetPermittivity(self, omega):
        return 1.0 + self.GetSusciptibility(omega)





        
        
        
