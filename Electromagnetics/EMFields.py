## EMFields.py  
## electromagnetic field, green functions, potentials ...


__all__ = ["EMFields2DHg", "EMFields3DHg"]


import numpy as np
from scipy import constants
from scipy.special import hankel2


class EMFields2DHg:
    """ homogeneous space
    """
    def __init__(self, eps_r=1.0, mu_r=1.0):
        self.JsourceE = []
        self.JsourceM = []
        self.eps_r = eps_r
        self.mu_r = mu_r
        return
        
    def setFreqency(self, f):
        self.omega = 2.0*np.pi*f

    def AddSourceE(self, rho_e, J_e):
        """ re: position 2-vector
            Je: 3-vector amplitude
        """
        self.JsourceE.append([rho_e, J_e])

    def AddSourceM(self, rho_m, J_m):
        """ rm: position
            Jm: 3-vector amplitude
        """
        self.JsourceM.append([rho_e, J_m])
        
    def getEH(self, rho_vec):
        
        ehshape = [3] + list(rho_vec[0].shape)
        E = np.zeros(ehshape, dtype=complex)
        H = np.zeros(ehshape, dtype=complex)
        
        epsilon = self.eps_r*constants.epsilon_0
        mu = self.mu_r*constants.mu_0
        omega = self.omega
        k = omega*np.sqrt(mu*epsilon)
        
        for i in range(len(self.JsourceE)):
            rho_e, J_e = self.JsourceE[i]
            
            J_x, J_y, J_z = J_e
            
            x = rho_vec[0] - rho_e[0]
            y = rho_vec[1] - rho_e[1]
            rho = np.sqrt(x**2 + y**2)

            E[0] +=  (-(4*J_x*epsilon*mu*omega**2*rho**6*hankel2(0, rho*k) + k*(2*J_x*rho**5*(hankel2(-1, rho*k) - hankel2(1, rho*k)) + rho**4*k*x*(J_x*x + J_y*y)*(hankel2(-2, rho*k) - 2*hankel2(0, rho*k) + hankel2(2, rho*k)) - 2*rho**3*x*(J_x*x + J_y*y)*(hankel2(-1, rho*k) - hankel2(1, rho*k))))/(16*epsilon*omega*rho**6))
            E[1] +=  (-(4*J_y*epsilon*mu*omega**2*rho**6*hankel2(0, rho*k) + k*(2*J_y*rho**5*(hankel2(-1, rho*k) - hankel2(1, rho*k)) + rho**4*k*y*(J_x*x + J_y*y)*(hankel2(-2, rho*k) - 2*hankel2(0, rho*k) + hankel2(2, rho*k)) - 2*rho**3*y*(J_x*x + J_y*y)*(hankel2(-1, rho*k) - hankel2(1, rho*k))))/(16*epsilon*omega*rho**6))
            E[2] +=  (-J_z*mu*omega*hankel2(0, rho*k)/4)        
            H[0] +=  (-1j*J_z*k*y*(hankel2(-1, rho*k) - hankel2(1, rho*k))/(8*rho))
            H[1] +=  (1j*J_z*k*x*(hankel2(-1, rho*k) - hankel2(1, rho*k))/(8*rho))
            H[2] +=  (1j*k*(J_x*y - J_y*x)*(hankel2(-1, rho*k) - hankel2(1, rho*k))/(8*rho))        
        
        for i in range(len(self.JsourceM)):
            rho_m, J_m = self.JsourceM[i]
            
            M_x, M_y, M_z = J_m
            
            x = rho_vec[0] - rho_m[0]
            y = rho_vec[1] - rho_m[1]
            rho = np.sqrt(x**2 + y**2)

            E[0] +=  (1j*M_z*mu*k*y*(hankel2(-1, rho*k) - hankel2(1, rho*k))/(8*epsilon*rho))
            E[1] +=  (-1j*M_z*mu*k*x*(hankel2(-1, rho*k) - hankel2(1, rho*k))/(8*epsilon*rho))
            E[2] +=  (1j*mu*k*(-M_x*y + M_y*x)*(hankel2(-1, rho*k) - hankel2(1, rho*k))/(8*epsilon*rho))
            H[0] +=  (-(4*M_x*epsilon*mu*omega**2*rho**6*hankel2(0, rho*k) + k*(2*M_x*rho**5*(hankel2(-1, rho*k) - hankel2(1, rho*k)) + rho**4*k*x*(M_x*x + M_y*y)*(hankel2(-2, rho*k) - 2*hankel2(0, rho*k) + hankel2(2, rho*k)) - 2*rho**3*x*(M_x*x + M_y*y)*(hankel2(-1, rho*k) - hankel2(1, rho*k))))/(16*epsilon*omega*rho**6))
            H[1] +=  (-(4*M_y*epsilon*mu*omega**2*rho**6*hankel2(0, rho*k) + k*(2*M_y*rho**5*(hankel2(-1, rho*k) - hankel2(1, rho*k)) + rho**4*k*y*(M_x*x + M_y*y)*(hankel2(-2, rho*k) - 2*hankel2(0, rho*k) + hankel2(2, rho*k)) - 2*rho**3*y*(M_x*x + M_y*y)*(hankel2(-1, rho*k) - hankel2(1, rho*k))))/(16*epsilon*omega*rho**6))
            H[2] +=  (-M_z*mu*omega*hankel2(0, rho*k)/4)

        return E, H
        


class EMFields3DHg:
    """ homogeneous space
    """
    def __init__(self, eps_r, mu_r):
        self.JsourceE = []
        self.JsourceM = []
        self.eps_r = eps_r
        self.mu_r = mu_r
        return
        
    def setFreqency(self, f):
        self.omega = 2.0*np.pi*f

    def AddSourceE(self, r_e, J_e):
        """ re: position
            Je: 3-vector amplitude
        """
        self.JsourceE.append([r_e, J_e])

    def AddSourceM(self, r_m, J_m):
        """ rm: position
            Jm: 3-vector amplitude
        """
        self.JsourceM.append([r_e, J_m])
        
    def getEH(self, r_vec):
        
        ehshape = [3] + list(rho_vec[0].shape)
        E = np.zeros(ehshape, dtype=complex)
        H = np.zeros(ehshape, dtype=complex)
        
        epsilon = self.eps_r*constants.epsilon_0
        mu = self.mu_r*constants.mu_0
        omega = self.omega
        k = omega*np.sqrt(mu*epsilon)

        for i in range(len(self.JsourceE)):
            r_e, J_e = self.JsourceE[i]
            
            J_x, J_y, J_z = J_e
            
            x = r_vec[0] - r_e[0]
            y = r_vec[1] - r_e[1]
            z = r_vec[2] - r_e[2]
            r = np.sqrt(x**2 + y**2 + z**2)

            E[0] +=  (1j*(-J_x*epsilon*mu*omega**2*r**4 + 1j*J_x*k*r**3 - 3*1j*k*r*x*(J_x*x + J_y*y + J_z*z) + r**2*(J_x*k**2*x**2 + J_x + J_y*k**2*x*y + J_z*k**2*x*z) - 3*x*(J_x*x + J_y*y + J_z*z))*exp(-1j*k*r)/(4*pi*epsilon*omega*r**5))
            E[1] +=  (1j*(-J_y*epsilon*mu*omega**2*r**4 + 1j*J_y*k*r**3 - 3*1j*k*r*y*(J_x*x + J_y*y + J_z*z) + r**2*(J_x*k**2*x*y + J_y*k**2*y**2 + J_y + J_z*k**2*y*z) - 3*y*(J_x*x + J_y*y + J_z*z))*exp(-1j*k*r)/(4*pi*epsilon*omega*r**5))
            E[2] +=  (1j*(-J_z*epsilon*mu*omega**2*r**4 + 1j*J_z*k*r**3 - 3*1j*k*r*z*(J_x*x + J_y*y + J_z*z) + r**2*(J_x*k**2*x*z + J_y*k**2*y*z + J_z*k**2*z**2 + J_z) - 3*z*(J_x*x + J_y*y + J_z*z))*exp(-1j*k*r)/(4*pi*epsilon*omega*r**5))
            H[0] +=  ((J_y*z - J_z*y)*(1j*k*r + 1)*exp(-1j*k*r)/(4*pi*r**3))
            H[1] +=  (-(J_x*z - J_z*x)*(1j*k*r + 1)*exp(-1j*k*r)/(4*pi*r**3))
            H[2] +=  ((J_x*y - J_y*x)*(1j*k*r + 1)*exp(-1j*k*r)/(4*pi*r**3))
        
        for i in range(len(self.JsourceM)):
            r_m, J_m = self.JsourceM[i]
            
            M_x, M_y, M_z = J_m
            
            x = r_vec[0] - r_m[0]
            y = r_vec[1] - r_m[1]
            z = r_vec[2] - r_m[2]
            r = np.sqrt(x**2 + y**2 + z**2)

            E[0] +=  (-mu*(M_y*z - M_z*y)*(1j*k*r + 1)*exp(-1j*k*r)/(4*pi*epsilon*r**3))
            E[1] +=  (mu*(M_x*z - M_z*x)*(1j*k*r + 1)*exp(-1j*k*r)/(4*pi*epsilon*r**3))
            E[2] +=  (-mu*(M_x*y - M_y*x)*(1j*k*r + 1)*exp(-1j*k*r)/(4*pi*epsilon*r**3))
            H[0] +=  (1j*(-M_x*epsilon*mu*omega**2*r**4 + 1j*M_x*k*r**3 - 3*1j*k*r*x*(M_x*x + M_y*y + M_z*z) + r**2*(M_x*k**2*x**2 + M_x + M_y*k**2*x*y + M_z*k**2*x*z) - 3*x*(M_x*x + M_y*y + M_z*z))*exp(-1j*k*r)/(4*pi*epsilon*omega*r**5))
            H[1] +=  (1j*(-M_y*epsilon*mu*omega**2*r**4 + 1j*M_y*k*r**3 - 3*1j*k*r*y*(M_x*x + M_y*y + M_z*z) + r**2*(M_x*k**2*x*y + M_y*k**2*y**2 + M_y + M_z*k**2*y*z) - 3*y*(M_x*x + M_y*y + M_z*z))*exp(-1j*k*r)/(4*pi*epsilon*omega*r**5))
            H[2] +=  (1j*(-M_z*epsilon*mu*omega**2*r**4 + 1j*M_z*k*r**3 - 3*1j*k*r*z*(M_x*x + M_y*y + M_z*z) + r**2*(M_x*k**2*x*z + M_y*k**2*y*z + M_z*k**2*z**2 + M_z) - 3*z*(M_x*x + M_y*y + M_z*z))*exp(-1j*k*r)/(4*pi*epsilon*omega*r**5))
            
        return E, H
            
            
            
        
        
        
