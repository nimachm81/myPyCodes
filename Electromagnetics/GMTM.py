## GMTM.py
# Gradient Metamaterials



__all__ = ["GMTM1D", "GMTM1D_ML"]


from scipy import constants
import cmath
import math

import numpy as np


class GMTM1D:
    ##WARNING: exp(-I*omega*t) convention

    def __init__(self):
        return
        
    def _SetProfileParameters_(self, s1, s2, L1, L2, vbose=False):
        self.s1 = s1
        self.s2 = s2
        self.L1 = L1
        self.L2 = L2
        if vbose:
            y = (L2/(2.0*L1))
            print('y = L2/2L1 = ', y)
            print('U_m = ', 1.0/(1.0+self.s1*y**2))
            print('n0*U_m = ', self.n0_barrier/(1.0+self.s1*y**2))
        
    def SetProfileParameters(self, n0_barrier, y, d, s1, s2, vbose=False):
        # d: barrier length
        self.n0_barrier = n0_barrier
        self.y = y
        self.d = d
        L1 = d/(4.0*y**2)
        L2 = d/(2.0*y)
        self._SetProfileParameters_(s1, s2, L1, L2, vbose)
        return [L1, L2]
        
    def SetMediumParams(self, N_layer, n_0, n_1):
        self.N_layer = N_layer
        self.n_0 = n_0
        self.n_1 = n_1
        
    def GetP(self):
        p_sq = self.s1**2/(4.0*self.L1**2) - self.s2/(self.L2**2)
        return np.sqrt(p_sq)
        
    def GetOmega(self):
        p = self.GetP()
        Omega = constants.c*p/self.n0_barrier
        return Omega
        
    def GetOmega__(self):
        if self.s2<0 and self.s1>0:     ## concave profile
            Omega = 2.0*constants.c*self.y*cmath.sqrt(1.0+self.y**2)/(self.n0_barrier*self.d)
            return Omega
        elif self.s2>0 and self.s1<0:   ## convex profile
            Omega = 2.0*constants.c*self.y*cmath.sqrt(1.0-self.y**2)/(self.n0_barrier*self.d)
            return Omega
            
    def omegaTou(self, omega):
        return abs(self.GetOmega())/omega
        
    def uTomega(self, u):
        return abs(self.GetOmega())/u
        
            
    def GetN(self, omega):
        Omega = self.GetOmega()
        N = np.sqrt(-Omega**2/omega**2 + 1.0+0j)
        return N
        
    def Getq(self, omega):
        N = self.GetN(omega)
        q = (N*omega*self.n0_barrier/constants.c)
        return q
        
        
    def GetU(self, z):
        U = 1.0/(1.0 + self.s1*z/self.L1 + self.s2*z**2/self.L2**2)
        return U     
        
    def GetdU(self, z):
        dU = (-2*self.s2*z/self.L2**2 - self.s1/self.L1)*self.GetU(z)**2
        return dU
        
    def GetUIntegral(self, z_0, z_1, vbose=False):
        from scipy.integrate import quad
        U_int = quad(lambda z: self.GetU(z), z_0, z_1)
        if vbose:
            print('U_int : ', U_int)
        return U_int[0]
        
    def Geteta(self):
        eta = self.GetUIntegral(0.0, self.d)
        return eta
                    
    def GetProfile(self, z):
        return self.n0_barrier*self.GetU(z)
        
        
    def GetReflection__(self, f_0, f_1, N_pts):
        ## one layer
        n_0, n_1 = self.n_0, self.n_1
        
        omega = np.linspace(f_0, f_1, N_pts)*2.0*math.pi
        k_0 = omega/constants.c
        k_i = k_0*n_0
        k_r = k_0*n_0
        k_t = k_0*n_1
        
        s_1, s_2, L_1, L_2 = self.s1, self.s2, self.L1, self.L2
        
        z = self.d
        U = 1/(1 + s_2*z**2/L_2**2 + s_1*z/L_1)
        Up_U2 = (-2*s_2*z/L_2**2 - s_1/L_1)
        eta_d = self.Geteta()
        q = self.Getq(omega)
        
        Q = ((-2*k_t/U + 2*q + 1j*Up_U2)*np.exp(2*1j*eta_d*q)/(2*k_t/U + 2*q - 1j*Up_U2))
        
        R = ((2*1j*L_1*Q*k_i + 2*1j*L_1*Q*q + 2*1j*L_1*k_i - 2*1j*L_1*q - Q*s_1 - s_1)/(2*1j*L_1*Q*k_r - 2*1j*L_1*Q*q + 2*1j*L_1*k_r + 2*1j*L_1*q + Q*s_1 + s_1))
    
        return R
    
    
    def GetReflection(self, f_0, f_1, N_pts):
        
        n_0, n_1 = self.n_0, self.n_1
        N_layer = self.N_layer        
        
        omega = None
        if N_pts==None:
            omega = f_0*2.0*math.pi
        else:
            omega = np.linspace(f_0, f_1, N_pts)*2.0*math.pi
            
        k_0 = omega/constants.c
        k_i = k_0*n_0
        k_r = k_0*n_0
        k_t = k_0*n_1
        
        s_1, s_2, L_1, L_2 = self.s1, self.s2, self.L1, self.L2
        
        z = self.d
        U = 1/(1 + s_2*z**2/L_2**2 + s_1*z/L_1)
        Up_U2 = (-2*s_2*z/L_2**2 - s_1/L_1)
        eta_d = self.Geteta()
        q = self.Getq(omega)
        
        Q = ((-2*k_t/U + 2*q + 1j*Up_U2)*np.exp(2*1j*eta_d*q)/(2*k_t/U + 2*q - 1j*Up_U2))
        
        Q_np1 = Q
        Q_n = Q
        
        U_0 = self.GetU(0.0)
        U_d = self.GetU(self.d)
        dU_0 = self.GetdU(0.0)
        dU_d = self.GetdU(self.d)
        for i in range(N_layer-1):
            Q_n = ((2*1j*Q_np1*q*U_0**2*U_d + 2*1j*Q_np1*q*U_0*U_d**2 - Q_np1*U_0*dU_d + Q_np1*U_d*dU_0 - 2*1j*q*U_0**2*U_d + 2*1j*q*U_0*U_d**2 - U_0*dU_d + U_d*dU_0)*np.exp(2*1j*eta_d*q)/(-2*1j*Q_np1*q*U_0**2*U_d + 2*1j*Q_np1*q*U_0*U_d**2 + Q_np1*U_0*dU_d - Q_np1*U_d*dU_0 + 2*1j*q*U_0**2*U_d + 2*1j*q*U_0*U_d**2 + U_0*dU_d - U_d*dU_0))
            Q_np1 = Q_n
        
        Q = Q_n
        R = ((2*1j*L_1*Q*k_i + 2*1j*L_1*Q*q + 2*1j*L_1*k_i - 2*1j*L_1*q - Q*s_1 - s_1)/(2*1j*L_1*Q*k_r - 2*1j*L_1*Q*q + 2*1j*L_1*k_r + 2*1j*L_1*q + Q*s_1 + s_1))
    
        return R
    

    def GetTransmissionReflection(self, f_0, f_1=None, N_pts=None):
        n_0, n_1 = self.n_0, self.n_1
        N_layer = self.N_layer        

        R = self.GetReflection(f_0, f_1, N_pts)
        
        omega = None
        if N_pts==None:
            omega = f_0*2.0*math.pi
        else:
            omega = np.linspace(f_0, f_1, N_pts)*2.0*math.pi

        k_0 = omega/constants.c
        k_i = k_0*n_0
        k_r = k_0*n_0
        k_t = k_0*n_1
        
        s_1, s_2, L_1, L_2 = self.s1, self.s2, self.L1, self.L2
        
        z = self.d
        U = 1.0/(1 + s_2*z**2/L_2**2 + s_1*z/L_1)
        Up_U2 = (-2*s_2*z/L_2**2 - s_1/L_1)
        eta_d = self.Geteta()
        q = self.Getq(omega)
        
        U_0 = self.GetU(0.0)
        U_d = self.GetU(self.d)
        dU_0 = self.GetdU(0.0)
        dU_d = self.GetdU(self.d)

        Q = ((-2*1j*L_1*R*k_r - 2*1j*L_1*R*q + 2*1j*L_1*k_i - 2*1j*L_1*q - R*s_1 - s_1)/(2*1j*L_1*R*k_r - 2*1j*L_1*R*q - 2*1j*L_1*k_i - 2*1j*L_1*q + R*s_1 + s_1))
        
        A = ((R + 1)*cmath.sqrt(U_0)/(Q + 1))
        
        Q_n = Q
        A_n = A

        Q_np1 = None
        A_np1 = None
        
        AQ_arr = [(A, Q)]
        
        for i in range(1, N_layer):
            Q_np1 = ((2*1j*Q_n*q*U_0**2*U_d + 2*1j*Q_n*q*U_0*U_d**2 + Q_n*U_0*dU_d - Q_n*U_d*dU_0 + 2*1j*q*U_0**2*U_d*np.exp(2*1j*eta_d*q) - 2*1j*q*U_0*U_d**2*np.exp(2*1j*eta_d*q) + U_0*np.exp(2*1j*eta_d*q)*dU_d - U_d*np.exp(2*1j*eta_d*q)*dU_0)/(2*1j*Q_n*q*U_0**2*U_d - 2*1j*Q_n*q*U_0*U_d**2 - Q_n*U_0*dU_d + Q_n*U_d*dU_0 + 2*1j*q*U_0**2*U_d*np.exp(2*1j*eta_d*q) + 2*1j*q*U_0*U_d**2*np.exp(2*1j*eta_d*q) - U_0*np.exp(2*1j*eta_d*q)*dU_d + U_d*np.exp(2*1j*eta_d*q)*dU_0))
            
            A_np1 = (A_n*(Q_n + np.exp(2*1j*eta_d*q))*cmath.sqrt(U_0)*np.exp(-1j*eta_d*q)/((Q_np1 + 1)*cmath.sqrt(U_d)))
            
            Q_n = Q_np1
            A_n = A_np1
            
            AQ_arr.append((A_n, Q_n))
        
        Q = Q_n
        A = A_n
        T = (A*(Q + np.exp(2*1j*eta_d*q))*np.exp(-1j*eta_d*q)/cmath.sqrt(U_d))
    
    
        if N_pts!=None:
            return [omega/(2.0*math.pi), T, R]
        else:
            return [T, R, AQ_arr]
    
        
    def GetFieldPlot(self, freq, n_pts_0, n_pts_i, n_pts_1, d0=None, d1=None, verbose=False):
        """ freq: frequency
            n_pts_i: number of points in each slab
            n_0: number of points in the incident (reflection) medium
            n_1: number of points in the transmission medium
        """
        n_0, n_1 = self.n_0, self.n_1

        lambda_0_z = constants.c/freq
        if d0==None:
            d0 = 2.0*lambda_0_z
        if d1==None:
            d1 = 2.0*lambda_0_z

        T, R, AQ_arr = self.GetTransmissionReflection(freq)
        
        z_pts_0 = np.linspace(-d0, 0.0, n_pts_0, endpoint=False)
        z_pts_i = np.linspace(0.0, self.N_layer*self.d, self.N_layer*n_pts_i, endpoint=False)
        z_pts_1 = np.linspace(self.N_layer*self.d, self.N_layer*self.d+d1,  n_pts_1, endpoint=True)

        omega = 2.0*math.pi*freq
        k_0 = omega/constants.c
        
        ##WARNING: e^-i*w*t time convention
        CONV = -1
        k_z = n_0*k_0
        a_b_i = np.array([1.0, R])
        E_FWD_0 = a_b_i[0]*np.exp(-CONV*1j*k_z*(z_pts_0-z_pts_i[0]))
        E_BWD_0 = a_b_i[1]*np.exp(+CONV*1j*k_z*(z_pts_0-z_pts_i[0]))
                
        E_FWD_i = np.zeros(z_pts_i.shape, dtype=complex)
        E_BWD_i = np.zeros(z_pts_i.shape, dtype=complex)

        for i in range(self.N_layer):
            q = self.Getq(omega)
            z_i = z_pts_i[i*n_pts_i: (i+1)*n_pts_i] - z_pts_i[i*n_pts_i]
            eta_z_i = np.zeros(n_pts_i)
            
            for j in range(n_pts_i):
                eta_z_i[j] = self.GetUIntegral(0.0, z_i[j])
            
            A_i, Q_i = AQ_arr[i]
            a_b_i = [A_i, A_i*Q_i]
            E_FWD_i[i*n_pts_i: (i+1)*n_pts_i] = a_b_i[0]*np.exp(-CONV*1j*q*(eta_z_i))/np.sqrt(self.GetU(z_i))
            E_BWD_i[i*n_pts_i: (i+1)*n_pts_i] = a_b_i[1]*np.exp(+CONV*1j*q*(eta_z_i))/np.sqrt(self.GetU(z_i))


        k_z = n_1*k_0
        a_b_i = np.array([T, 0.0])
        E_FWD_1 = a_b_i[0]*np.exp(-CONV*1j*k_z*(z_pts_1-z_pts_1[0]))
        E_BWD_1 = a_b_i[1]*np.exp(+CONV*1j*k_z*(z_pts_1-z_pts_1[0]))

        z_pts = np.concatenate((z_pts_0, z_pts_i, z_pts_1))
        E_FWD = np.concatenate((E_FWD_0, E_FWD_i, E_FWD_1))
        E_BWD = np.concatenate((E_BWD_0, E_BWD_i, E_BWD_1))
        
        return [z_pts, E_FWD, E_BWD]


    def GetMediumPlot(self, n_pts_0, n_pts_i, n_pts_1, d0=None, d1=None, freq=None):
        n_0, n_1 = self.n_0, self.n_1
        
        if d0==None:
            assert freq!=None
            lambda_0_z = constants.c/freq
            d0 = 2.0*lambda_0_z
        if d1==None:
            assert freq!=None
            lambda_0_z = constants.c/freq
            d1 = 2.0*lambda_0_z
        
        
        z_pts_0 = np.linspace(-d0, 0.0, n_pts_0, endpoint=False)
        z_pts_i = np.linspace(0.0, self.N_layer*self.d, self.N_layer*n_pts_i, endpoint=False)
        z_pts_1 = np.linspace(self.N_layer*self.d, self.N_layer*self.d+d1,  n_pts_1, endpoint=True)
        
        n_0_vec = n_0*np.ones(n_pts_0)
        n_1_vec = n_1*np.ones(n_pts_1)
        
        n_i_vec = np.zeros(z_pts_i.shape, dtype=complex)

        for i in range(self.N_layer):
            z_i = z_pts_i[i*n_pts_i: (i+1)*n_pts_i] - z_pts_i[i*n_pts_i]
            n_i_vec[i*n_pts_i: (i+1)*n_pts_i] = self.GetProfile(z_i)

        z_pts = np.concatenate((z_pts_0, z_pts_i, z_pts_1))
        n_vec = np.concatenate((n_0_vec, n_i_vec, n_1_vec))
        
        return [z_pts, n_vec]



class GMTM1D_ML:
    ##WARNING: exp(-I*omega*t) convention
    """ multilayer gradient metamaterial where each layer can take different
    material properties
    """
    def __init__(self):
        return
        
    def _SetProfileParameters_(self, s1, s2, L1, L2, vbose=False):
        self.s1 = s1
        self.s2 = s2
        self.L1 = L1
        self.L2 = L2
        if vbose:
            y = (L2/(2.0*L1))
            print('y = L2/2L1 = ', y)
            print('U_m = ', 1.0/(1.0+self.s1*y**2))
            print('n0*U_m = ', self.n0_barrier/(1.0+self.s1*y**2))
        
    def SetProfileParameters(self, n0_barrier, y, d, s1, s2, vbose=False):
        # d: barrier length
        self.N_layer = len(n0_barrier)
        self.n0_barrier = n0_barrier
        self.y = y
        self.d = d
        L1 = d/(4.0*y**2)
        L2 = d/(2.0*y)
        self._SetProfileParameters_(s1, s2, L1, L2, vbose)
        return [L1, L2]
        
    def SetMediumParams(self, n_0, n_1):
        self.n_0 = n_0
        self.n_1 = n_1
        
    def GetP(self):
        p_sq = self.s1**2/(4.0*self.L1**2) - self.s2/(self.L2**2)
        return np.sqrt(p_sq)
        
    def GetOmega(self):
        p = self.GetP()
        Omega = constants.c*p/self.n0_barrier
        return Omega
                    
    def omegaTou(self, i, omega):
        return abs(self.GetOmega()[i])/omega
        
    def uTomega(self, i, u):
        return abs(self.GetOmega()[i])/u
        
            
    def GetN(self, i, omega):
        Omega = self.GetOmega()
        N = np.sqrt(-Omega[i]**2/omega**2 + 1.0+0j)
        return N
        
    def Getq(self, i, omega):
        N = self.GetN(i, omega)
        q = (N*omega*self.n0_barrier[i]/constants.c)
        return q
        
        
    def GetU(self, i, z):
        U = 1.0/(1.0 + self.s1[i]*z/self.L1[i] + self.s2[i]*z**2/self.L2[i]**2)
        return U     
        
    def GetdU(self, i, z):
        dU = (-2*self.s2[i]*z/self.L2[i]**2 - self.s1[i]/self.L1[i])*self.GetU(i, z)**2
        return dU
        
    def GetUIntegral(self, i, z_0, z_1, vbose=False):
        from scipy.integrate import quad
        U_int = quad(lambda z: self.GetU(i, z), z_0, z_1)
        if vbose:
            print('U_int : ', U_int)
        return U_int[0]
        
    def Geteta(self, i):
        eta = self.GetUIntegral(i, 0.0, self.d[i])
        return eta
                    
    def GetProfile(self, i, z):
        return self.n0_barrier[i]*self.GetU(i, z)
        
    def GetnAverage(self):
        avg = 0.0
        D = 0.0
        for i in range(self.N_layer):
            avg += self.n0_barrier[i]*self.GetUIntegral(i, 0.0, self.d[i])
            D += self.d[i]
        return avg/D
        
    def GetReflection(self, f_0, f_1, N_pts):
        n_0, n_1 = self.n_0, self.n_1
        N_layer = self.N_layer        
        
        omega = None
        if N_pts==None:
            omega = f_0*2.0*math.pi
        else:
            omega = np.linspace(f_0, f_1, N_pts)*2.0*math.pi
            
        k_0 = omega/constants.c
        k_i = k_0*n_0
        k_r = k_0*n_0
        k_t = k_0*n_1
        
        s_1, s_2, L_1, L_2 = self.s1[N_layer-1], self.s2[N_layer-1], self.L1[N_layer-1], self.L2[N_layer-1]
        
        z = self.d[N_layer-1]
        U = self.GetU(N_layer-1, z) #1.0/(1 + s_2*z**2/L_2**2 + s_1*z/L_1)
        Up_U2 = (-2*s_2*z/L_2**2 - s_1/L_1)
        eta_d = self.Geteta(N_layer-1)
        q = self.Getq(N_layer-1, omega)
        
        Q = ((-2*k_t/U + 2*q + 1j*Up_U2)*np.exp(2*1j*eta_d*q)/(2*k_t/U + 2*q - 1j*Up_U2))  ## final layer
        
        Q_np1 = Q
        Q_n = Q
        
        for i in range(N_layer-2, -1, -1):
            q_n = self.Getq(i, omega)
            q_np1 = self.Getq(i+1, omega)
            eta_nd = self.Geteta(i)

            Unp1_0 = self.GetU(i+1, 0.0)
            Un_d = self.GetU(i, self.d[i])
            dUnp1_0 = self.GetdU(i+1, 0.0)
            dUn_d = self.GetdU(i, self.d[i])

            Q_n = ((2*1j*Q_np1*q_n*Un_d**2*Unp1_0 + 2*1j*Q_np1*q_np1*Un_d*Unp1_0**2 + Q_np1*Un_d*dUnp1_0 - Q_np1*Unp1_0*dUn_d + 2*1j*q_n*Un_d**2*Unp1_0 - 2*1j*q_np1*Un_d*Unp1_0**2 + Un_d*dUnp1_0 - Unp1_0*dUn_d)*np.exp(2*1j*eta_nd*q_n)/(2*1j*Q_np1*q_n*Un_d**2*Unp1_0 - 2*1j*Q_np1*q_np1*Un_d*Unp1_0**2 - Q_np1*Un_d*dUnp1_0 + Q_np1*Unp1_0*dUn_d + 2*1j*q_n*Un_d**2*Unp1_0 + 2*1j*q_np1*Un_d*Unp1_0**2 - Un_d*dUnp1_0 + Unp1_0*dUn_d))
            Q_np1 = Q_n
        
        s_1, s_2, L_1, L_2 = self.s1[0], self.s2[0], self.L1[0], self.L2[0]

        z = self.d[0]
        U = self.GetU(0, z) #1.0/(1 + s_2*z**2/L_2**2 + s_1*z/L_1)
        Up_U2 = (-2*s_2*z/L_2**2 - s_1/L_1)
        eta_d = self.Geteta(0)
        q = self.Getq(0, omega)

        Q = Q_n
        R = ((2*1j*L_1*Q*k_i + 2*1j*L_1*Q*q + 2*1j*L_1*k_i - 2*1j*L_1*q - Q*s_1 - s_1)/(2*1j*L_1*Q*k_r - 2*1j*L_1*Q*q + 2*1j*L_1*k_r + 2*1j*L_1*q + Q*s_1 + s_1))
    
        return R

    def GetTransmissionReflection(self, f_0, f_1=None, N_pts=None):
        n_0, n_1 = self.n_0, self.n_1
        N_layer = self.N_layer        

        R = self.GetReflection(f_0, f_1, N_pts)
        
        omega = None
        if N_pts==None:
            omega = f_0*2.0*math.pi
        else:
            omega = np.linspace(f_0, f_1, N_pts)*2.0*math.pi

        k_0 = omega/constants.c
        k_i = k_0*n_0
        k_r = k_0*n_0
        k_t = k_0*n_1
        
        s_1, s_2, L_1, L_2 = self.s1[0], self.s2[0], self.L1[0], self.L2[0]
        
        s_1, s_2, L_1, L_2 = self.s1[0], self.s2[0], self.L1[0], self.L2[0]

        z = self.d[0]
        U = self.GetU(0, z) #1.0/(1 + s_2*z**2/L_2**2 + s_1*z/L_1)
        Up_U2 = (-2*s_2*z/L_2**2 - s_1/L_1)
        eta_d = self.Geteta(0)
        q = self.Getq(0, omega)
        
        U_0 = self.GetU(0, 0.0)
        U_d = self.GetU(0, self.d[0])
        dU_0 = self.GetdU(0, 0.0)
        dU_d = self.GetdU(0, self.d[0])

        Q = ((-2*1j*L_1*R*k_r - 2*1j*L_1*R*q + 2*1j*L_1*k_i - 2*1j*L_1*q - R*s_1 - s_1)/(2*1j*L_1*R*k_r - 2*1j*L_1*R*q - 2*1j*L_1*k_i - 2*1j*L_1*q + R*s_1 + s_1))
        
        A = ((R + 1)*cmath.sqrt(U_0)/(Q + 1))
        
        Q_n = Q
        A_n = A

        Q_np1 = None
        A_np1 = None
        
        AQ_arr = [(A, Q)]
        
        for i in range(N_layer-1):
            q_n = self.Getq(i, omega)
            q_np1 = self.Getq(i+1, omega)
            eta_nd = self.Geteta(i)

            Unp1_0 = self.GetU(i+1, 0.0)
            Un_d = self.GetU(i, self.d[i])
            dUnp1_0 = self.GetdU(i+1, 0.0)
            dUn_d = self.GetdU(i, self.d[i])
        
            Q_np1 = ((2*1j*Q_n*q_n*Un_d**2*Unp1_0 + 2*1j*Q_n*q_np1*Un_d*Unp1_0**2 - Q_n*Un_d*dUnp1_0 + Q_n*Unp1_0*dUn_d - 2*1j*q_n*Un_d**2*Unp1_0*np.exp(2*1j*eta_nd*q_n) + 2*1j*q_np1*Un_d*Unp1_0**2*np.exp(2*1j*eta_nd*q_n) - Un_d*np.exp(2*1j*eta_nd*q_n)*dUnp1_0 + Unp1_0*np.exp(2*1j*eta_nd*q_n)*dUn_d)/(-2*1j*Q_n*q_n*Un_d**2*Unp1_0 + 2*1j*Q_n*q_np1*Un_d*Unp1_0**2 + Q_n*Un_d*dUnp1_0 - Q_n*Unp1_0*dUn_d + 2*1j*q_n*Un_d**2*Unp1_0*np.exp(2*1j*eta_nd*q_n) + 2*1j*q_np1*Un_d*Unp1_0**2*np.exp(2*1j*eta_nd*q_n) + Un_d*np.exp(2*1j*eta_nd*q_n)*dUnp1_0 - Unp1_0*np.exp(2*1j*eta_nd*q_n)*dUn_d))
            
            A_np1 =(A_n*(Q_n + np.exp(2*1j*eta_nd*q_n))*np.sqrt(Unp1_0)*np.exp(-1j*eta_nd*q_n)/((Q_np1 + 1)*np.sqrt(Un_d)))
            
            Q_n = Q_np1
            A_n = A_np1
            
            AQ_arr.append((A_n, Q_n))
        
        Q = Q_n
        A = A_n

        eta_d = self.Geteta(N_layer-1)
        q = self.Getq(N_layer-1, omega)
        U_d = self.GetU(N_layer-1, self.d[N_layer-1])
        
        #print('A:', A, '\nQ: ', Q, '\neta_d: ', eta_d, '\nq:', q, '\nU_d: ', U_d)
        T = (A*(Q + np.exp(2*1j*eta_d*q))*np.exp(-1j*eta_d*q)/cmath.sqrt(U_d))
    
        if N_pts!=None:
            return [omega/(2.0*math.pi), T, R]
        else:
            return [T, R, AQ_arr]
    
        
    def GetFieldPlot(self, freq, n_pts_0, n_pts_i, n_pts_1, d0=None, d1=None, verbose=False):
        """ freq: frequency
            n_pts_i: number of points in each slab
            n_0: number of points in the incident (reflection) medium
            n_1: number of points in the transmission medium
        """
        n_0, n_1 = self.n_0, self.n_1

        lambda_0_z = constants.c/freq
        if d0==None:
            d0 = 2.0*lambda_0_z
        if d1==None:
            d1 = 2.0*lambda_0_z

        T, R, AQ_arr = self.GetTransmissionReflection(freq)
        
        z_pts_0 = np.linspace(-d0, 0.0, n_pts_0, endpoint=False)
        z_start = 0.0
        z_pts_i = np.linspace(z_start, z_start+self.d[0], n_pts_i, endpoint=False)
        for i in range(1, self.N_layer):
            z_start += self.d[i-1]
            z_pts_i = np.concatenate((z_pts_i, np.linspace(z_start, z_start+self.d[i], n_pts_i, endpoint=False)))
        z_pts_1 = np.linspace(sum(self.d), sum(self.d)+d1,  n_pts_1, endpoint=True)

        omega = 2.0*math.pi*freq
        k_0 = omega/constants.c
        
        ##WARNING: e^-i*w*t time convention
        CONV = -1
        k_z = n_0*k_0
        a_b_i = np.array([1.0, R])
        E_FWD_0 = a_b_i[0]*np.exp(-CONV*1j*k_z*(z_pts_0-z_pts_i[0]))
        E_BWD_0 = a_b_i[1]*np.exp(+CONV*1j*k_z*(z_pts_0-z_pts_i[0]))
                
        E_FWD_i = np.zeros(z_pts_i.shape, dtype=complex)
        E_BWD_i = np.zeros(z_pts_i.shape, dtype=complex)

        for i in range(self.N_layer):
            q = self.Getq(i, omega)
            z_i = z_pts_i[i*n_pts_i: (i+1)*n_pts_i] - z_pts_i[i*n_pts_i]
            eta_z_i = np.zeros(n_pts_i)
            
            for j in range(n_pts_i):
                eta_z_i[j] = self.GetUIntegral(i, 0.0, z_i[j])
            
            A_i, Q_i = AQ_arr[i]
            a_b_i = [A_i, A_i*Q_i]
            E_FWD_i[i*n_pts_i: (i+1)*n_pts_i] = a_b_i[0]*np.exp(-CONV*1j*q*(eta_z_i))/np.sqrt(self.GetU(i, z_i))
            E_BWD_i[i*n_pts_i: (i+1)*n_pts_i] = a_b_i[1]*np.exp(+CONV*1j*q*(eta_z_i))/np.sqrt(self.GetU(i, z_i))


        k_z = n_1*k_0
        a_b_i = np.array([T, 0.0])
        E_FWD_1 = a_b_i[0]*np.exp(-CONV*1j*k_z*(z_pts_1-z_pts_1[0]))
        E_BWD_1 = a_b_i[1]*np.exp(+CONV*1j*k_z*(z_pts_1-z_pts_1[0]))

        z_pts = np.concatenate((z_pts_0, z_pts_i, z_pts_1))
        E_FWD = np.concatenate((E_FWD_0, E_FWD_i, E_FWD_1))
        E_BWD = np.concatenate((E_BWD_0, E_BWD_i, E_BWD_1))
        
        return [z_pts, E_FWD, E_BWD]



    def GetMediumPlot(self, n_pts_0, n_pts_i, n_pts_1, d0=None, d1=None, freq=None):
        n_0, n_1 = self.n_0, self.n_1
        
        if d0==None:
            assert freq!=None
            lambda_0_z = constants.c/freq
            d0 = 2.0*lambda_0_z
        if d1==None:
            assert freq!=None
            lambda_0_z = constants.c/freq
            d1 = 2.0*lambda_0_z
        
        
        z_pts_0 = np.linspace(-d0, 0.0, n_pts_0, endpoint=False)
        z_start = 0.0
        z_pts_i = np.linspace(z_start, z_start+self.d[0], n_pts_i, endpoint=False)
        for i in range(1, self.N_layer):
            z_start += self.d[i-1]
            z_pts_i = np.concatenate((z_pts_i, np.linspace(z_start, z_start+self.d[i], n_pts_i, endpoint=False)))
        z_pts_1 = np.linspace(sum(self.d), sum(self.d)+d1,  n_pts_1, endpoint=True)
        
        n_0_vec = n_0*np.ones(n_pts_0)
        n_1_vec = n_1*np.ones(n_pts_1)
        
        n_i_vec = np.zeros(z_pts_i.shape, dtype=complex)

        for i in range(self.N_layer):
            z_i = z_pts_i[i*n_pts_i: (i+1)*n_pts_i] - z_pts_i[i*n_pts_i]
            n_i_vec[i*n_pts_i: (i+1)*n_pts_i] = self.GetProfile(i, z_i)

        z_pts = np.concatenate((z_pts_0, z_pts_i, z_pts_1))
        n_vec = np.concatenate((n_0_vec, n_i_vec, n_1_vec))
        
        return [z_pts, n_vec]





class GMTM1D_T:
    ##exp(-I*k*z) convention
    """ Gradient time varying material
    """
    def __init__(self):
        return
        
    def _SetProfileParameters_(self, s_1, s_2, T_1, T_2, vbose=False):
        self.s_1 = s_1
        self.s_2 = s_2
        self.T_1 = T_1
        self.T_2 = T_2
        if vbose:
            y = (T_2/(2*T_1))
            print('y = T2/2T1 = ', y)
            U_m = (-s_1**2*y**2/s_2 + 1)
            print('U_m = ', U_m)
            print('n0*U_m = ', self.n0_barrier*U_m)
            t_d = (-T_2**2*s_1/(T_1*s_2))
            print('t_d = ', t_d)
        
    def SetProfileParameters(self, n0_barrier, y, t_d, s_1, s_2, vbose=False):
        # d: barrier thickness
        self.n0_barrier = n0_barrier
        self.y = y
        self.t_d = t_d
        T_1 = (-s_2*t_d/(4*s_1*y**2))
        T_2 = (-s_2*t_d/(2*s_1*y))
        self._SetProfileParameters_(s_1, s_2, T_1, T_2, vbose)
        return [T_1, T_2]
        
    def SetMediumParams(self, N_layer, n_0, n_1):
        self.N_layer = N_layer
        self.n_0 = n_0
        self.n_1 = n_1
        assert N_layer==1
        
    def GetPsq(self):
        s_1, s_2, T_1, T_2 = self.s_1, self.s_2, self.T_1, self.T_2
        p_sq = (-s_2/T_2**2 + s_1**2/(4*T_1**2))
        return p_sq
        
    def GetTransitionK(self):
        # transition from evanescence to propagating
        p_sq = self.GetPsq()
        p_abs = np.sqrt(np.abs(p_sq))
        return self.n0_barrier*p_abs/constants.c
        
    def GetOmega(self, k):
        p_sq = self.GetPsq()
        Omega_sq = (constants.c**2*k**2/self.n0_barrier**2 - p_sq)
        return np.sqrt(Omega_sq+0j)
                        
    def GetU(self, t):
        s_1, s_2, T_1, T_2 = self.s_1, self.s_2, self.T_1, self.T_2
        U = (1.0 + s_2*t**2/T_2**2 + s_1*t/T_1)
        return U     
        
    def GetdU(self, t):
        s_1, s_2, T_1, T_2 = self.s_1, self.s_2, self.T_1, self.T_2
        dU = (2*s_2*t/T_2**2 + s_1/T_1)
        return dU
        
    def GetUinvIntegral(self, t_0, t_1, vbose=False):
        from scipy.integrate import quad
        U_int = quad(lambda t: 1.0/self.GetU(t), t_0, t_1)
        if vbose:
            print('U_int : ', U_int)
        return U_int[0]
        
    def GetUIntegral(self, t_0, t_1, vbose=False):
        from scipy.integrate import quad
        U_int = quad(lambda t: self.GetU(t), t_0, t_1)
        if vbose:
            print('U_int : ', U_int)
        return U_int[0]

    def Geteta(self, t):
        eta = self.GetUinvIntegral(0.0, t)
        return eta
                    
    def GetProfile(self, t):
        return self.n0_barrier*self.GetU(t)
        
        
    def GetTransmissionPlusMinus(self, k_0, k_1=None, N_pts=None):
        ## one layer
        n_0, n_1 = self.n_0, self.n_1
        
        k = None
        if k_1!=None:
            k = np.linspace(k_0, k_1, N_pts)
        else:
            k = k_0

        w_0 = k*constants.c
        omega_1 = w_0/n_0
        omega_3 = w_0/n_1
        
        Omega = self.GetOmega(k)
        
        s_1, s_2, T_1, T_2 = self.s_1, self.s_2, self.T_1, self.T_2

        U_0 = self.GetU(0.0) + 0j
        
        A_p = ((2*T_1*(Omega + omega_1*U_0) + 1j*s_1)/(4*T_1*Omega*np.sqrt(U_0)))
        A_m = ((2*T_1*(Omega - omega_1*U_0) - 1j*s_1)/(4*T_1*Omega*np.sqrt(U_0)))
        
        A_pm_arr = []        
        
        T = self.t_d
        U_T = self.GetU(T) + 0j
        dU_T = self.GetdU(T)
        eta_T = self.Geteta(T)
        
        
        A_p3_p = ((2*Omega + 2*omega_3*U_T - 1j*dU_T)*np.exp(1j*Omega*eta_T)/(4*omega_3*np.sqrt(U_T)))
        A_p3_m = ((-2*Omega + 2*omega_3*U_T - 1j*dU_T)*np.exp(-1j*Omega*eta_T)/(4*omega_3*np.sqrt(U_T)))
        A_m3_p = ((-2*Omega + 2*omega_3*U_T + 1j*dU_T)*np.exp(1j*Omega*eta_T)/(4*omega_3*np.sqrt(U_T)))
        A_m3_m = ((2*Omega + 2*omega_3*U_T + 1j*dU_T)*np.exp(-1j*Omega*eta_T)/(4*omega_3*np.sqrt(U_T)))
        
        TM = np.array([[A_p3_p, A_p3_m],
                       [A_m3_p, A_m3_m]])
                       
                       
        A_pm = np.array([A_p, A_m])
    
        T_p = TM[0,0]*A_pm[0] + TM[0,1]*A_pm[1]
        T_m = TM[1,0]*A_pm[0] + TM[1,1]*A_pm[1]
        
        A_pm_arr.append(A_pm)
        
        if k_1!=None:
            return [k, T_p, T_m]
        else:
            return [[T_p, T_m], A_pm_arr]
    
    
        
    def GetFieldPlot(self, k, n_pts_0, n_pts_i, n_pts_1, t_d0=None, t_d1=None, verbose=False):
        """ freq: frequency
            n_pts_i: number of points in each slab
            n_0: number of points in the incident (reflection) medium
            n_1: number of points in the transmission medium
        """
        n_0, n_1 = self.n_0, self.n_1

        T__0 = 2.0*math.pi/(k*constants.c)
        if t_d0==None:
            t_d0 = 2.0*T__0
        if t_d1==None:
            t_d1 = 2.0*T__0

        T_pm, A_pm = self.GetTransmissionPlusMinus(k)
        
        t_pts_0 = np.linspace(-t_d0, 0.0, n_pts_0, endpoint=False)
        t_pts_i = np.linspace(0.0, self.N_layer*self.t_d, self.N_layer*n_pts_i, endpoint=False)
        t_pts_1 = np.linspace(self.N_layer*self.t_d, self.N_layer*self.t_d+t_d1,  n_pts_1, endpoint=True)

        w_0 = k*constants.c
        omega_1 = w_0/n_0
        omega_3 = w_0/n_1
        
        Omega = self.GetOmega(k)

        ##e^-i*k*z convention
        a_b_i = np.array([1.0, 0.0])
        E_FWD_0 = a_b_i[0]*np.exp(+1j*omega_1*(t_pts_0-t_pts_i[0]))
        E_BWD_0 = a_b_i[1]*np.exp(-1j*omega_1*(t_pts_0-t_pts_i[0]))
                
        E_FWD_i = np.zeros(t_pts_i.shape, dtype=complex)
        E_BWD_i = np.zeros(t_pts_i.shape, dtype=complex)

        for i in range(self.N_layer):
            t_i = t_pts_i[i*n_pts_i: (i+1)*n_pts_i] - t_pts_i[i*n_pts_i]
            eta_t_i = np.zeros(n_pts_i)
            
            for j in range(n_pts_i):
                eta_t_i[j] = self.Geteta(t_i[j])
            
            A_pi, A_mi = A_pm[i]
            a_b_i = [A_pi, A_mi]
            E_FWD_i[i*n_pts_i: (i+1)*n_pts_i] = a_b_i[0]*np.exp(+1j*Omega*(eta_t_i))*np.sqrt(self.GetU(t_i) + 0j)
            E_BWD_i[i*n_pts_i: (i+1)*n_pts_i] = a_b_i[1]*np.exp(-1j*Omega*(eta_t_i))*np.sqrt(self.GetU(t_i) + 0j)


        a_b_i = T_pm
        E_FWD_1 = a_b_i[0]*np.exp(+1j*omega_3*(t_pts_1-t_pts_1[0]))
        E_BWD_1 = a_b_i[1]*np.exp(-1j*omega_3*(t_pts_1-t_pts_1[0]))

        t_pts = np.concatenate((t_pts_0, t_pts_i, t_pts_1))
        E_FWD = np.concatenate((E_FWD_0, E_FWD_i, E_FWD_1))
        E_BWD = np.concatenate((E_BWD_0, E_BWD_i, E_BWD_1))
        
        return [t_pts, E_FWD, E_BWD]


    def GetMediumPlot(self, n_pts_0, n_pts_i, n_pts_1, t_d0=None, t_d1=None, k=None):
        n_0, n_1 = self.n_0, self.n_1
        
        if t_d0==None:
            assert k!=None
            T__0 = 2.0*math.pi/(k*constants.c)
            t_d0 = 2.0*T__0
        if t_d1==None:
            assert k!=None
            T__0 = 2.0*math.pi/(k*constants.c)
            t_d1 = 2.0*T__0
        
        
        t_pts_0 = np.linspace(-t_d0, 0.0, n_pts_0, endpoint=False)
        t_pts_i = np.linspace(0.0, self.N_layer*self.t_d, self.N_layer*n_pts_i, endpoint=False)
        t_pts_1 = np.linspace(self.N_layer*self.t_d, self.N_layer*self.t_d+t_d1,  n_pts_1, endpoint=True)
        
        n_0_vec = n_0*np.ones(n_pts_0)
        n_1_vec = n_1*np.ones(n_pts_1)
        
        n_i_vec = np.zeros(t_pts_i.shape, dtype=complex)

        for i in range(self.N_layer):
            t_i = t_pts_i[i*n_pts_i: (i+1)*n_pts_i] - t_pts_i[i*n_pts_i]
            n_i_vec[i*n_pts_i: (i+1)*n_pts_i] = self.GetProfile(t_i)

        t_pts = np.concatenate((t_pts_0, t_pts_i, t_pts_1))
        n_vec = np.concatenate((n_0_vec, n_i_vec, n_1_vec))
        
        return [t_pts, n_vec]




