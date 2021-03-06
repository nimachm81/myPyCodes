## GSTCFD Freqquency domain GSTC 


__all__ = ["MSTypes", "PlaneWaveFD", "GSTCFD_2D", "GSTCFDSpectralAnalyzer"]

import numpy as np
from enum import Enum
from scipy import constants
import scipy
from scipy import linalg
import matplotlib.pyplot as plt


class MSTypes(Enum):
    scalar = 0
    
    
def IndexVariables(arrs):
    # arrs: list of array shapes
    arrs_inds = []
    ind_st = 0
    for i in range(len(arrs)):
        shape_i = arrs[i]
        ind_last = ind_st+np.prod(shape_i)
        x_ind_i = np.arange(ind_st, ind_last).reshape(shape_i)
        
        arrs_inds.append(x_ind_i)
        ind_st = ind_last
    return arrs_inds


class PlaneWaveFD:
    """ plane wave transformations
        the metasurface is on z=0 plane
    """
    
    def __init__(self, frequency, vbose=False):
        self.vbose = vbose
        self.PWTs = []  ## PW transformations
        self.f = frequency
        
        self.eta_0 = np.sqrt(constants.mu_0/constants.epsilon_0)
        self.eps_r1 = 1.0
        self.mu_r1 = 1.0
        self.eps_r2 = 1.0
        self.mu_r2 = 1.0
        return
        
    @staticmethod
    def AngleToVec(theta, phi):
        a_z = np.array([0.0, 0.0, 1.0])
        a_r = np.array([np.cos(phi)*np.sin(theta), np.sin(phi)*np.sin(theta), np.cos(theta)])
        a_phi = np.cross(a_z, a_r)
        if theta==0.0:
            a_phi = np.array([-np.sin(phi), np.cos(phi), 0.0])
        a_phi /= np.linalg.norm(a_phi)
        a_theta = np.cross(a_phi, a_r)
        return [a_r, a_theta, a_phi]
            
        
    def SetPW(self, E_i, E_r, E_t, k_i, k_r, k_t):
        k_i = k_i/np.linalg.norm(k_i)
        k_r = k_r/np.linalg.norm(k_r)
        k_t = k_t/np.linalg.norm(k_t)
        
        assert k_i[2]*k_r[2] <= 0.0 
        assert k_i[2]*k_t[2] >= 0.0 
        
        E_i = E_i - k_i.dot(E_i)*k_i
        E_r = E_r - k_r.dot(E_r)*k_r
        E_t = E_t - k_t.dot(E_t)*k_t
        
        H_i = np.cross(k_i, E_i)/self.eta_0
        H_r = np.cross(k_r, E_r)/self.eta_0
        H_t = np.cross(k_t, E_t)/self.eta_0
        
        pwt = [[E_i, E_r, E_t], [H_i, H_r, H_t], [k_i, k_r, k_t]]
        self.PWTs.append(pwt)
        
        
    def GetIRTFieldsOnMS(self, N_y=200, nl_y=1.0, l_y_max=None):
        """ electric and magnetic fields on the metasurface (z=0)
            the metasurface extends from -nl_y*l_y_max to nl_y*l_y_max
            if l_y_max==None: it is calculated automatically
        """
        omega = 2.0*np.pi*self.f
        k_0 = omega/constants.c

        a_z = np.array([0.0, 0.0, 1.0])
        E_arr, H_arr, k_arr = self.PWTs[0]
        E_i, E_r, E_t = E_arr
        H_i, H_r, H_t = H_arr 
        ak_i, ak_r, ak_t = k_arr
        assert np.abs(ak_i[0])<1.0e-15 and np.abs(ak_r[0])<1.0e-15 and np.abs(ak_t[0])<1.0e-15

        k_iy = ak_i[1]*k_0
        k_ry = ak_r[1]*k_0
        k_ty = ak_t[1]*k_0
        
        ky_irt = []
        if np.abs(ak_i[1])>1.0e-15:
            ky_irt.append(k_iy)
        if np.abs(ak_r[1])>1.0e-15:
            ky_irt.append(k_ry)
        if np.abs(ak_t[1])>1.0e-15:
            ky_irt.append(k_ty)
        
        if l_y_max==None:
            l_y_max = 1.0
            if len(ky_irt)>0: 
                k_y_min = np.min(np.abs(np.array(ky_irt)))
                if k_y_min>0.0:
                    l_y_max = 2.0*np.pi/k_y_min
        
        y_max = nl_y*l_y_max
        y = np.linspace(-y_max/2, y_max/2, N_y, endpoint=False)
        
        E_i_vec = np.zeros((3, N_y), dtype=complex)
        E_r_vec = np.zeros((3, N_y), dtype=complex)
        E_t_vec = np.zeros((3, N_y), dtype=complex)
        H_i_vec = np.zeros((3, N_y), dtype=complex)
        H_r_vec = np.zeros((3, N_y), dtype=complex)
        H_t_vec = np.zeros((3, N_y), dtype=complex)
        
        assert ak_i[2]>=0 ## excitation from other side not implemented

        for i in range(3):
            E_i_vec[i,:] = E_i[i]*np.exp(-1j*k_iy*y)
            E_r_vec[i,:] = E_r[i]*np.exp(-1j*k_ry*y)
            E_t_vec[i,:] = E_t[i]*np.exp(-1j*k_ty*y)
            
            H_i_vec[i,:] = H_i[i]*np.exp(-1j*k_iy*y)
            H_r_vec[i,:] = H_r[i]*np.exp(-1j*k_ry*y)
            H_t_vec[i,:] = H_t[i]*np.exp(-1j*k_ty*y)
        return y, [E_i_vec, E_r_vec, E_t_vec], [H_i_vec, H_r_vec, H_t_vec]


    def GetFieldsOnMS(self, y, IRT='ALL'):
        """ electric and magnetic fields on the metasurface (z=0)
        IRT: ALL/I/R/T
        """
        N_y = len(y)
        omega = 2.0*np.pi*self.f
        k_0 = omega/constants.c

        a_z = np.array([0.0, 0.0, 1.0])
        E_arr, H_arr, k_arr = self.PWTs[0]
        E_i, E_r, E_t = E_arr
        H_i, H_r, H_t = H_arr 
        ak_i, ak_r, ak_t = k_arr
        assert np.abs(ak_i[0])<1.0e-15 and np.abs(ak_r[0])<1.0e-15 and np.abs(ak_t[0])<1.0e-15

        k_iy = ak_i[1]*k_0
        k_ry = ak_r[1]*k_0
        k_ty = ak_t[1]*k_0
                
        E_i_vec = np.zeros((3, N_y), dtype=complex)
        E_r_vec = np.zeros((3, N_y), dtype=complex)
        E_t_vec = np.zeros((3, N_y), dtype=complex)
        H_i_vec = np.zeros((3, N_y), dtype=complex)
        H_r_vec = np.zeros((3, N_y), dtype=complex)
        H_t_vec = np.zeros((3, N_y), dtype=complex)
        
        assert ak_i[2]>=0 ## excitation from other side not implemented

        for i in range(3):
            E_i_vec[i,:] = E_i[i]*np.exp(-1j*k_iy*y)
            E_r_vec[i,:] = E_r[i]*np.exp(-1j*k_ry*y)
            E_t_vec[i,:] = E_t[i]*np.exp(-1j*k_ty*y)
            
            H_i_vec[i,:] = H_i[i]*np.exp(-1j*k_iy*y)
            H_r_vec[i,:] = H_r[i]*np.exp(-1j*k_ry*y)
            H_t_vec[i,:] = H_t[i]*np.exp(-1j*k_ty*y)
        if IRT=='ALL':
            return [E_i_vec, E_r_vec, E_t_vec], [H_i_vec, H_r_vec, H_t_vec]
        elif IRT=='I':
            return E_i_vec, H_i_vec
        elif IRT=='R':
            return E_r_vec, H_r_vec
        elif IRT=='T':
            return E_t_vec, H_t_vec
        else:
            raise ValueError()
            
    

class GSTCFD_2D:
    """ 2D problem, or 1D metasurface on z=0 plane
    """
    
    def __init__(self, frequency, msType, vbose=False):
        self.vbose = vbose
        self.PWTs = []  ## PW transformations
        self.f = frequency
        self.msType = msType
        
        self.eta_0 = np.sqrt(constants.mu_0/constants.epsilon_0)
        self.eps_r1 = 1.0
        self.mu_r1 = 1.0
        self.eps_r2 = 1.0
        self.mu_r2 = 1.0
        return
        
    def SetFields(self, y, E_irt, H_irt):
        self.y = y
        self.E_i, self.E_r, self.E_t = E_irt
        self.H_i, self.H_r, self.H_t = H_irt
        return
        
    def GetSusceptibilities(self):
        N_y = len(self.y)
        omega = 2.0*np.pi*self.f
        k_0 = omega/constants.c
        if self.msType==MSTypes.scalar:
            a_z = np.array([0.0, 0.0, 1.0])
            
            E_p = self.E_t
            E_m = self.E_i + self.E_r 
            zxdE = np.zeros((3, N_y), dtype=complex)
            E_avg = np.zeros((3, N_y), dtype=complex)
            H_p = self.H_t
            H_m = self.H_i + self.H_r
            zxdH = np.zeros((3, N_y), dtype=complex)
            H_avg = np.zeros((3, N_y), dtype=complex)

            for j in range(N_y):
                zxdH[:, j] = np.cross(a_z, H_p[:,j]-H_m[:,j])
                zxdE[:, j] = np.cross(a_z, E_p[:,j]-E_m[:,j])
                
            P = zxdH/(1j*omega)
            M = -zxdE/(1j*omega*constants.mu_0)
            
            H_avg = (H_p + H_m)/2.0
            E_avg = (E_p + E_m)/2.0
            
            X_ee, X_mmm = None, None
            
            if np.sum(np.abs(self.E_i[0]))>np.sum(np.abs(self.E_i[1])):   ##TODO: fix
                if self.vbose:
                    print('X polarized incidence..')
                X_ee = P[0,:]/(E_avg[0,:]*constants.epsilon_0)
                X_mm = M[1,:]/H_avg[1,:]
            else:
                if self.vbose:
                    print('Y polarized incidence..')
                X_ee = P[1,:]/(E_avg[1,:]*constants.epsilon_0)
                X_mm = M[0,:]/H_avg[0,:]
                
            return [X_ee, X_mm], [P, M], [E_avg, H_avg]
        else:
            raise NotImplementedError()
    
    
class GSTCFDSpectralAnalyzer:
    """ spectral domain analysis
        the metasurface is on z=0 plane
    """
    
    def __init__(self, frequency, msType, vbose=False):
        self.vbose = vbose
        self.f = frequency
        self.msType = msType
        
        self.eta_0 = np.sqrt(constants.mu_0/constants.epsilon_0)
        self.eps_r1 = 1.0
        self.mu_r1 = 1.0
        self.eps_r2 = 1.0
        self.mu_r2 = 1.0
        return

    def getKz(self, k, k_x, k_y):
        k_z = np.sqrt(k**2 - k_x**2 - k_y**2 + 0j)
        proper = (np.imag(k_z)<=0.0)
        k_z = proper*k_z - np.logical_not(proper)*k_z 
        return k_z
    
    def getKy(self, y):
        N_y = len(y)
        dy = y[1]-y[0]
        k_y = 2.0*np.pi*np.fft.fftfreq(N_y, d=dy)
        k_y = np.fft.fftshift(k_y)
        #print('k_y:', k_y)
        return k_y
        
    def getEzf(self, Ex_f, Ey_f, k_x, k_y, k_z):
        Ez_f = -(Ex_f*k_x + Ey_f*k_y)/k_z
        return Ez_f

    def AnalyzeSpectral1DScalar(self, y, X_ee, X_mm, E_i, compare_iterative=False):
        omega = 2.0*np.pi*self.f
        N_y = len(y)
        dy = y[1]-y[0]
        X_ee_f = np.fft.ifftshift(np.fft.ifft(X_ee))#/N_y
        X_mm_f = np.fft.ifftshift(np.fft.ifft(X_mm))#/N_y
        ##TODO: multiply by exp(j*k*y_min) to set the origin at y_min 
        ## and make consistent changes throught the whole code
        assert len(E_i)>=2
        E_i_f = np.zeros((len(E_i),N_y), dtype=complex)
        for i in range(len(E_i)):
            E_i_f[i,:] = np.fft.ifftshift(np.fft.ifft(E_i[i,:]))#/N_y
        
        k_y = self.getKy(y)
        
        E_r_Tf_size = (2,N_y)
        E_t_Tf_size = (2,N_y)
        N_tot = 4*N_y
        
        E_r_Tf_inds, E_t_Tf_inds = IndexVariables([E_r_Tf_size, E_t_Tf_size])
       
        if self.vbose:
            plt.plot(E_r_Tf_inds[0,:])
            plt.plot(E_r_Tf_inds[1,:])
            plt.plot(E_t_Tf_inds[0,:])
            plt.plot(E_t_Tf_inds[1,:])
            plt.ylabel('indices', fontsize=18)
            plt.show()
        
        A = np.zeros((N_tot, N_tot), dtype=complex)
        rhs = np.zeros(N_tot, dtype=complex)
        
        eps_r1, mu_r1 = self.eps_r1, self.mu_r1
        eps_r2, mu_r2 = self.eps_r2, self.mu_r2
        if eps_r1!=eps_r2 or mu_r1!=mu_r2:
            raise NotImplementedError()
        
        epsilon_0 = constants.epsilon_0
        mu_0 = constants.mu_0
        mu_1 = mu_r1*mu_0
        mu_2 = mu_r2*mu_0
        
        k_0 = omega/constants.c
        k_1 = k_0*np.sqrt(mu_r1*eps_r1)
        k_2 = k_0*np.sqrt(mu_r2*eps_r2)
        
        j0 = int(N_y/2)
        
        k_x_i = 0.0
        for j1 in range(N_y):
            k_y_j1 = k_y[j1]
            k_z1_j1 = self.getKz(k_1, k_x_i, k_y_j1)
            k_z2_j1 = self.getKz(k_2, k_x_i, k_y_j1)
            
            ind_erx_j1 = E_r_Tf_inds[0, j1]
            ind_ery_j1 = E_r_Tf_inds[1, j1]
            ind_etx_j1 = E_t_Tf_inds[0, j1]
            ind_ety_j1 = E_t_Tf_inds[1, j1]
            
            E_xi_j1 = E_i_f[0, j1]
            E_yi_j1 = E_i_f[1, j1]
            
            if np.abs(k_z1_j1)<1.0e-5:
                print('k_z1_j1==0.0 ', k_z1_j1)
                k_z1_j1 = dk_y/10.0

            A[ind_erx_j1, ind_erx_j1] += (-k_z1_j1/(mu_1*omega))
            A[ind_erx_j1, ind_ery_j1] += (0)
            A[ind_erx_j1, ind_etx_j1] += (-k_z1_j1/(mu_1*omega))
            A[ind_erx_j1, ind_ety_j1] += (0)
            A[ind_ery_j1, ind_erx_j1] += (0)
            A[ind_ery_j1, ind_ery_j1] += (-(k_y_j1**2 + k_z1_j1**2)/(mu_1*omega*k_z1_j1))
            A[ind_ery_j1, ind_etx_j1] += (0)
            A[ind_ery_j1, ind_ety_j1] += (-(k_y_j1**2 + k_z1_j1**2)/(mu_1*omega*k_z1_j1))
            A[ind_etx_j1, ind_erx_j1] += (0)
            A[ind_etx_j1, ind_ery_j1] += (1)
            A[ind_etx_j1, ind_etx_j1] += (0)
            A[ind_etx_j1, ind_ety_j1] += (-1)
            A[ind_ety_j1, ind_erx_j1] += (-1)
            A[ind_ety_j1, ind_ery_j1] += (0)
            A[ind_ety_j1, ind_etx_j1] += (1)
            A[ind_ety_j1, ind_ety_j1] += (0)

            rhs[ind_erx_j1] += (-E_xi_j1*k_z1_j1/(mu_1*omega))
            rhs[ind_ery_j1] += (-E_yi_j1*(k_y_j1**2 + k_z1_j1**2)/(mu_1*omega*k_z1_j1))
            rhs[ind_etx_j1] += (-E_yi_j1)
            rhs[ind_ety_j1] += (E_xi_j1)
            
            
            for j2 in range(N_y):
                J1 = j1 - j0
                J2 = j2 - j0
                if 0<=J1-J2+j0<N_y and 0<=J2+j0<N_y:
                    k_y_j2 = k_y[J2+j0]
                    k_z1_j2 = self.getKz(k_1, k_x_i, k_y_j2)
                    k_z2_j2 = self.getKz(k_2, k_x_i, k_y_j2)

                    if np.abs(k_z1_j2)<1.0e-6:
                        ##TODO: fix
                        dk_y = k_y[1]-k_y[0]
                        if self.vbose:
                            print('k_z1_j2==0.0 ', k_z1_j2)
                        k_z1_j2 = dk_y/10.0


                    ind_erx_j2 = E_r_Tf_inds[0, J2+j0]
                    ind_ery_j2 = E_r_Tf_inds[1, J2+j0]
                    ind_etx_j2 = E_t_Tf_inds[0, J2+j0]
                    ind_ety_j2 = E_t_Tf_inds[1, J2+j0]
                    
                    E_xi_j2 = E_i_f[0, J2+j0]
                    E_yi_j2 = E_i_f[1, J2+j0]
                    ind_X_ = J1-J2+j0

                    A[ind_erx_j1, ind_erx_j2] += (-1j*X_ee_f[ind_X_]*epsilon_0*omega/2)
                    A[ind_erx_j1, ind_ery_j2] += (0)
                    A[ind_erx_j1, ind_etx_j2] += (-1j*X_ee_f[ind_X_]*epsilon_0*omega/2)
                    A[ind_erx_j1, ind_ety_j2] += (0)
                    A[ind_ery_j1, ind_erx_j2] += (0)
                    A[ind_ery_j1, ind_ery_j2] += (-1j*X_ee_f[ind_X_]*epsilon_0*omega/2)
                    A[ind_ery_j1, ind_etx_j2] += (0)
                    A[ind_ery_j1, ind_ety_j2] += (-1j*X_ee_f[ind_X_]*epsilon_0*omega/2)
                    A[ind_etx_j1, ind_erx_j2] += (0)
                    A[ind_etx_j1, ind_ery_j2] += (1j*X_mm_f[ind_X_]*mu_0*(k_y_j2**2 + k_z1_j2**2)/(2*mu_1*k_z1_j2))
                    A[ind_etx_j1, ind_etx_j2] += (0)
                    A[ind_etx_j1, ind_ety_j2] += (-1j*X_mm_f[ind_X_]*mu_0*(k_y_j2**2 + k_z1_j2**2)/(2*mu_1*k_z1_j2))
                    A[ind_ety_j1, ind_erx_j2] += (-1j*X_mm_f[ind_X_]*mu_0*k_z1_j2/(2*mu_1))
                    A[ind_ety_j1, ind_ery_j2] += (0)
                    A[ind_ety_j1, ind_etx_j2] += (1j*X_mm_f[ind_X_]*mu_0*k_z1_j2/(2*mu_1))
                    A[ind_ety_j1, ind_ety_j2] += (0)

                    rhs[ind_erx_j1] += (1j*E_xi_j2*X_ee_f[ind_X_]*epsilon_0*omega/2)
                    rhs[ind_ery_j1] += (1j*E_yi_j2*X_ee_f[ind_X_]*epsilon_0*omega/2)
                    rhs[ind_etx_j1] += (1j*E_yi_j2*X_mm_f[ind_X_]*mu_0*(k_y_j2**2 + k_z1_j2**2)/(2*mu_1*k_z1_j2))
                    rhs[ind_ety_j1] += (-1j*E_xi_j2*X_mm_f[ind_X_]*mu_0*k_z1_j2/(2*mu_1))
                    
                            
        A_d = np.diagonal(A).copy()
        #A_d_inv = 1.0/A_d
        #from scipy.sparse import dia_matrix
        #A_d_inv = dia_matrix((A_d_inv, 0), shape=A.shape)
        #A = A_d_inv*A
        #rhs = A_d_inv*rhs
        
        if self.vbose:
            plt.plot(np.abs(A_d))
            plt.ylabel('diagonals', fontsize=18)
            plt.show()
        
        
        self.A = A
        self.rhs = rhs
        self.N_tot = N_tot
        
        if self.vbose:
            plt.matshow(np.abs(A), cmap=plt.cm.Blues)
            plt.colorbar()
            plt.show()        
            print('A_max: ', np.max(np.abs(A)))


        if compare_iterative:
            from scipy.sparse.linalg import gmres
            res = gmres(A, rhs, x0=np.random.rand(N_tot), restart=100)
            print('info: ', res[1])
            x_vec = res[0]
            self.x_vec = x_vec

            if self.vbose:
                plt.plot(np.abs(x_vec), 'g')
                plt.ylabel(r'x_{vec} iterative', fontsize=18)
                plt.show()
            
                plt.plot(np.abs(rhs), 'b')
                plt.plot(np.abs(A.dot(x_vec)).T, 'r-.')
                #plt.gca().set_ylim([0.0, 2000])
                plt.ylabel('rhs iterative', fontsize=18)
                plt.show()

        x_vec = linalg.solve(A, rhs)
        self.x_vec = x_vec

        if self.vbose:
            plt.plot(np.abs(x_vec), 'g')
            plt.ylabel(r'x_{vec}', fontsize=18)
            plt.show()
        
            plt.plot(np.abs(rhs), 'b')
            plt.plot(np.abs(A.dot(x_vec)).T, 'r-.')
            #plt.gca().set_ylim([0.0, 2000])
            plt.ylabel('rhs', fontsize=18)
            plt.show()


                
        E_r_Tf = np.zeros((2, N_y), dtype=complex)
        E_t_Tf = np.zeros((2, N_y), dtype=complex)
        E_r_Tf[0, :] = x_vec[E_r_Tf_inds[0,:]]
        E_r_Tf[1, :] = x_vec[E_r_Tf_inds[1,:]]
        E_t_Tf[0, :] = x_vec[E_t_Tf_inds[0,:]]
        E_t_Tf[1, :] = x_vec[E_t_Tf_inds[1,:]]

        if self.vbose:
            plt.plot(np.abs(E_r_Tf[0, :]), 'b', label='rx')
            plt.plot(np.abs(E_r_Tf[1, :]), 'b-.', label='ry')
            plt.plot(np.abs(E_t_Tf[0, :]), 'r', label='tx')
            plt.plot(np.abs(E_t_Tf[1, :]), 'r-.', label='ty')
            plt.ylabel('$E_{rf}, E_{tf}$', fontsize=18)
            plt.legend(fontsize=18)
            plt.show()
            
            
            plt.plot(np.abs(E_i_f[0, :]), 'r', label='x')
            plt.plot(np.abs(E_i_f[1, :]), 'b', label='y')
            if len(E_i)>2:
                plt.plot(np.abs(E_i_f[2, :]), 'g', label='z')
            plt.ylabel('$E_{if}$', fontsize=18)
            plt.legend(fontsize=18)
            plt.show()

        E_r_T = np.zeros((2, N_y), dtype=complex)
        E_t_T = np.zeros((2, N_y), dtype=complex)
        for i in range(2):
            E_r_T[i,:] = np.fft.fft(np.fft.ifftshift(E_r_Tf[i,:]))#*N_y
            E_t_T[i,:] = np.fft.fft(np.fft.ifftshift(E_t_Tf[i,:]))#*N_y
        
        return E_r_T, E_t_T
        
        
    def GetCoeffsSparse1DScalar(self, y, X_ee, X_mm, E_i):
        omega = 2.0*np.pi*self.f
        N_y = len(y)
        dy = y[1]-y[0]
        X_ee_f = np.fft.ifftshift(np.fft.ifft(X_ee))#/N_y
        X_mm_f = np.fft.ifftshift(np.fft.ifft(X_mm))#/N_y
        assert len(E_i)>=2
        E_i_f = np.zeros((len(E_i),N_y), dtype=complex)
        for i in range(len(E_i)):
            E_i_f[i,:] = np.fft.ifftshift(np.fft.ifft(E_i[i,:]))#/N_y
        
        k_y = self.getKy(y)
        
        E_r_Tf_size = (2,N_y)
        E_t_Tf_size = (2,N_y)
        N_tot = 4*N_y
        
        E_r_Tf_inds, E_t_Tf_inds = IndexVariables([E_r_Tf_size, E_t_Tf_size])
       
        if self.vbose:
            plt.plot(E_r_Tf_inds[0,:])
            plt.plot(E_r_Tf_inds[1,:])
            plt.plot(E_t_Tf_inds[0,:])
            plt.plot(E_t_Tf_inds[1,:])
            plt.ylabel('indices', fontsize=18)
            plt.show()
        
        #A = np.zeros((N_tot, N_tot), dtype=complex)
        row, col, data = [], [], []
        rhs = np.zeros(N_tot, dtype=complex)
        
        eps_r1, mu_r1 = self.eps_r1, self.mu_r1
        eps_r2, mu_r2 = self.eps_r2, self.mu_r2
        if eps_r1!=eps_r2 or mu_r1!=mu_r2:
            raise NotImplementedError()
        
        epsilon_0 = constants.epsilon_0
        mu_0 = constants.mu_0
        mu_1 = mu_r1*mu_0
        mu_2 = mu_r2*mu_0
        
        k_0 = omega/constants.c
        k_1 = k_0*np.sqrt(mu_r1*eps_r1)
        k_2 = k_0*np.sqrt(mu_r2*eps_r2)
        
        j0 = int(N_y/2)
        
        k_x_i = 0.0
        for j1 in range(N_y):
            k_y_j1 = k_y[j1]
            k_z1_j1 = self.getKz(k_1, k_x_i, k_y_j1)
            k_z2_j1 = self.getKz(k_2, k_x_i, k_y_j1)
            
            ind_erx_j1 = E_r_Tf_inds[0, j1]
            ind_ery_j1 = E_r_Tf_inds[1, j1]
            ind_etx_j1 = E_t_Tf_inds[0, j1]
            ind_ety_j1 = E_t_Tf_inds[1, j1]
            
            E_xi_j1 = E_i_f[0, j1]
            E_yi_j1 = E_i_f[1, j1]
            
            if np.abs(k_z1_j1)<1.0e-5:
                print('k_z1_j1==0.0 ', k_z1_j1)
                k_z1_j1 = (k_y[1]-k_y[0])/10.0

            row.append(ind_erx_j1); col.append(ind_erx_j1); data.append((-k_z1_j1/(mu_1*omega)))
            row.append(ind_erx_j1); col.append(ind_ery_j1); data.append((0))
            row.append(ind_erx_j1); col.append(ind_etx_j1); data.append((-k_z1_j1/(mu_1*omega)))
            row.append(ind_erx_j1); col.append(ind_ety_j1); data.append((0))
            row.append(ind_ery_j1); col.append(ind_erx_j1); data.append((0))
            row.append(ind_ery_j1); col.append(ind_ery_j1); data.append((-(k_y_j1**2 + k_z1_j1**2)/(mu_1*omega*k_z1_j1)))
            row.append(ind_ery_j1); col.append(ind_etx_j1); data.append((0))
            row.append(ind_ery_j1); col.append(ind_ety_j1); data.append((-(k_y_j1**2 + k_z1_j1**2)/(mu_1*omega*k_z1_j1)))
            row.append(ind_etx_j1); col.append(ind_erx_j1); data.append((0))
            row.append(ind_etx_j1); col.append(ind_ery_j1); data.append((1))
            row.append(ind_etx_j1); col.append(ind_etx_j1); data.append((0))
            row.append(ind_etx_j1); col.append(ind_ety_j1); data.append((-1))
            row.append(ind_ety_j1); col.append(ind_erx_j1); data.append((-1))
            row.append(ind_ety_j1); col.append(ind_ery_j1); data.append((0))
            row.append(ind_ety_j1); col.append(ind_etx_j1); data.append((1))
            row.append(ind_ety_j1); col.append(ind_ety_j1); data.append((0))

            rhs[ind_erx_j1] += (-E_xi_j1*k_z1_j1/(mu_1*omega))
            rhs[ind_ery_j1] += (-E_yi_j1*(k_y_j1**2 + k_z1_j1**2)/(mu_1*omega*k_z1_j1))
            rhs[ind_etx_j1] += (-E_yi_j1)
            rhs[ind_ety_j1] += (E_xi_j1)
            
        ##get H_i
        k_y = self.getKy(y)
        j0 = int(N_y/2)

        k_x = 0.0
        k_z = +self.getKz(k_1, k_x, k_y)
        mu = mu_r1*constants.mu_0

        H_i_f = np.zeros((3,N_y), dtype=complex)
        for j in range(N_y):
            k_vec_j = np.array([k_x, k_y[j], k_z[j]])
            H_i_f[:,j] = np.cross(k_vec_j, E_i_f[:,j])/(omega*mu)
    
        H_ix = np.zeros(N_y, dtype=complex)
        H_iy = np.zeros(N_y, dtype=complex)
        H_iz = np.zeros(N_y, dtype=complex)
        for j in range(N_y):
            k_y_j = k_y[j]
            k_z_j = k_z[j]
            k_dot_r = k_y_j*y
            H_ix += H_i_f[0, j]*np.exp(-1j*k_dot_r)
            H_iy += H_i_f[1, j]*np.exp(-1j*k_dot_r)
            H_iz += H_i_f[2, j]*np.exp(-1j*k_dot_r)
        
        H_ix = np.fft.fftshift(H_ix, axes=(0,))
        H_iy = np.fft.fftshift(H_iy, axes=(0,))
        H_iz = np.fft.fftshift(H_iz, axes=(0,))
        
        H_i = np.array([H_ix, H_iy, H_iz])
        ##
        
        Ptx_i = 0.5*epsilon_0*X_ee*E_i[0]
        Pty_i = 0.5*epsilon_0*X_ee*E_i[1]
        Mtx_i = 0.5*X_mm*H_i[0]
        Mty_i = 0.5*X_mm*H_i[1]
        
        inds_erx = E_r_Tf_inds[0, np.arange(N_y)]
        inds_ery = E_r_Tf_inds[1, np.arange(N_y)]
        inds_etx = E_t_Tf_inds[0, np.arange(N_y)]
        inds_ety = E_t_Tf_inds[1, np.arange(N_y)]
        
        rhs[inds_erx] += 1j*omega*np.fft.ifftshift(np.fft.ifft(Ptx_i))
        rhs[inds_ery] += 1j*omega*np.fft.ifftshift(np.fft.ifft(Pty_i))
        rhs[inds_etx] += -1j*omega*mu_0*np.fft.ifftshift(np.fft.ifft(Mtx_i))
        rhs[inds_ety] += -1j*omega*mu_0*np.fft.ifftshift(np.fft.ifft(Mty_i))
            
        if self.vbose:
            plt.plot(np.abs(rhs), 'b')
            if self.rhs!=None:
                plt.plot(np.abs(self.rhs), 'r')
            plt.ylabel('sparse rhs', fontsize=18)
            plt.show()
        
        self.A_sp = scipy.sparse.coo_matrix((data, (row, col)), shape=(N_tot, N_tot)).tocsr()
        self.rhs = rhs

        self.N_tot = N_tot
        self.y = y
        self.X_ee = X_ee
        self.X_mm = X_mm
        self.E_r_Tf_inds, self.E_t_Tf_inds = E_r_Tf_inds, E_t_Tf_inds
        return
        
    def GetMatrixVectorProduct(self, Ert_f_vec):
        y = self.y
        X_ee = self.X_ee
        X_mm = self.X_mm
        E_r_Tf_inds, E_t_Tf_inds = self.E_r_Tf_inds, self.E_t_Tf_inds
        
        omega = 2.0*np.pi*self.f
        N_y = len(y)

        inds_erx = E_r_Tf_inds[0, np.arange(N_y)]
        inds_ery = E_r_Tf_inds[1, np.arange(N_y)]
        inds_etx = E_t_Tf_inds[0, np.arange(N_y)]
        inds_ety = E_t_Tf_inds[1, np.arange(N_y)]
        
        Er_f = np.zeros((3,N_y), dtype=complex)
        Et_f = np.zeros((3,N_y), dtype=complex)
        Er_f[0] = Ert_f_vec[inds_erx]
        Er_f[1] = Ert_f_vec[inds_ery]
        Et_f[0] = Ert_f_vec[inds_etx]
        Et_f[1] = Ert_f_vec[inds_ety]

        ## E_f[2]
        epsilon_0 = constants.epsilon_0
        mu_0 = constants.mu_0
        eps_r1, mu_r1 = self.eps_r1, self.mu_r1
        eps_r2, mu_r2 = self.eps_r2, self.mu_r2
        k_0 = omega/constants.c
        k_1 = k_0*np.sqrt(mu_r1*eps_r1)
        k_2 = k_0*np.sqrt(mu_r2*eps_r2)

        k_y = self.getKy(y)
        j0 = int(N_y/2)

        k_x = 0.0
        k_zr = -self.getKz(k_1, k_x, k_y)
        k_zt = +self.getKz(k_2, k_x, k_y)
        mu_1 = mu_r1*constants.mu_0
        mu_2 = mu_r2*constants.mu_0

        Er_f[2] = self.getEzf(Er_f[0], Er_f[1], k_x, k_y, k_zr)
        Et_f[2] = self.getEzf(Et_f[0], Et_f[1], k_x, k_y, k_zt)
        
        ##get E_r, E_t
        E_r = np.zeros((3,N_y), dtype=complex)
        E_t = np.zeros((3,N_y), dtype=complex)
        for j in range(N_y):
            k_y_j = k_y[j]
            k_dot_r = k_y_j*y
            E_r[0] += Er_f[0, j]*np.exp(-1j*k_dot_r)
            E_r[1] += Er_f[1, j]*np.exp(-1j*k_dot_r)
            E_r[2] += Er_f[2, j]*np.exp(-1j*k_dot_r)
            E_t[0] += Et_f[0, j]*np.exp(-1j*k_dot_r)
            E_t[1] += Et_f[1, j]*np.exp(-1j*k_dot_r)
            E_t[2] += Et_f[2, j]*np.exp(-1j*k_dot_r)
        
        E_r[0] = np.fft.fftshift(E_r[0], axes=(0,))
        E_r[1] = np.fft.fftshift(E_r[1], axes=(0,))
        E_r[2] = np.fft.fftshift(E_r[2], axes=(0,))
        E_t[0] = np.fft.fftshift(E_t[0], axes=(0,))
        E_t[1] = np.fft.fftshift(E_t[1], axes=(0,))
        E_t[2] = np.fft.fftshift(E_t[2], axes=(0,))
                
        ##get H_r, H_t
        H_r_f = np.zeros((3,N_y), dtype=complex)
        H_t_f = np.zeros((3,N_y), dtype=complex)
        for j in range(N_y):
            kr_vec_j = np.array([k_x, k_y[j], k_zr[j]])
            kt_vec_j = np.array([k_x, k_y[j], k_zt[j]])
            H_r_f[:,j] = np.cross(kr_vec_j, Er_f[:,j])/(omega*mu_1)
            H_t_f[:,j] = np.cross(kt_vec_j, Et_f[:,j])/(omega*mu_2)
    
        H_r = np.zeros((3,N_y), dtype=complex)
        H_t = np.zeros((3,N_y), dtype=complex)
        for j in range(N_y):
            k_y_j = k_y[j]
            k_dot_r = k_y_j*y
            H_r[0] += H_r_f[0, j]*np.exp(-1j*k_dot_r)
            H_r[1] += H_r_f[1, j]*np.exp(-1j*k_dot_r)
            H_r[2] += H_r_f[2, j]*np.exp(-1j*k_dot_r)
            H_t[0] += H_t_f[0, j]*np.exp(-1j*k_dot_r)
            H_t[1] += H_t_f[1, j]*np.exp(-1j*k_dot_r)
            H_t[2] += H_t_f[2, j]*np.exp(-1j*k_dot_r)
        
        H_r[0] = np.fft.fftshift(H_r[0], axes=(0,))
        H_r[1] = np.fft.fftshift(H_r[1], axes=(0,))
        H_r[2] = np.fft.fftshift(H_r[2], axes=(0,))
        H_t[0] = np.fft.fftshift(H_t[0], axes=(0,))
        H_t[1] = np.fft.fftshift(H_t[1], axes=(0,))
        H_t[2] = np.fft.fftshift(H_t[2], axes=(0,))
        
        ##
        res = self.A_sp*Ert_f_vec

        Ptx_r = 0.5*epsilon_0*X_ee*E_r[0]
        Pty_r = 0.5*epsilon_0*X_ee*E_r[1]
        Mtx_r = 0.5*X_mm*H_r[0]
        Mty_r = 0.5*X_mm*H_r[1]
        Ptx_t = 0.5*epsilon_0*X_ee*E_t[0]
        Pty_t = 0.5*epsilon_0*X_ee*E_t[1]
        Mtx_t = 0.5*X_mm*H_t[0]
        Mty_t = 0.5*X_mm*H_t[1]

        res[inds_erx] += -1j*omega*np.fft.ifftshift(np.fft.ifft(Ptx_r + Ptx_t))
        res[inds_ery] += -1j*omega*np.fft.ifftshift(np.fft.ifft(Pty_r + Pty_t))
        res[inds_etx] += +1j*omega*mu_0*np.fft.ifftshift(np.fft.ifft(Mtx_r + Mtx_t))
        res[inds_ety] += +1j*omega*mu_0*np.fft.ifftshift(np.fft.ifft(Mty_r + Mty_t))

        """
        plt.plot(np.abs(res), 'b')
        if np.all(self.A!=None):
            res_2 = self.A.dot(Ert_f_vec)
            plt.plot(np.abs(res_2), 'r-.')
            print(res.shape, res_2.shape)
        plt.show()
        """
            
        return res
        

    def SolveIterative(self, tol=1.0e-3, maxiter=1000):
        N_tot = self.N_tot
        
        def matvec(v):
            return self.GetMatrixVectorProduct(v)
        
        from scipy.sparse.linalg import LinearOperator
        
        A = LinearOperator(dtype=complex, shape=(N_tot,N_tot), matvec=matvec)
        A = self.A
        
        from scipy.sparse.linalg import gmres
        
        res = gmres(A, self.rhs, tol=tol, maxiter=maxiter)
        return  res
        
        
    def GetYZField(self, y, E_T, Z=1.0, N_z=100, n_y_repeat=0, getH=False):
        """ y: y samples on the metasurface
            E_T: electric fild transverse vector on the y samples
            Z: samples (y, z=0..Z) plane (Z>0 for transmission, Z<0 for reflection)
            N_z = number of z samples
        """
        assert Z!=0.0
        omega = 2.0*np.pi*self.f
        N_y = len(y)
        dy = y[1]-y[0]
        E_f = np.zeros((3,N_y), dtype=complex)
        for i in range(2):
            E_f[i,:] = np.fft.ifftshift(np.fft.ifft(E_T[i,:]))#/N_y
            
        eps_r1, mu_r1 = self.eps_r1, self.mu_r1
        eps_r2, mu_r2 = self.eps_r2, self.mu_r2
        k_0 = omega/constants.c
        k_1 = k_0*np.sqrt(mu_r1*eps_r1)
        k_2 = k_0*np.sqrt(mu_r2*eps_r2)

        k_y = self.getKy(y)
        j0 = int(N_y/2)

        k_x = 0.0
        k_z = None
        z = None
        mu = None
        if Z>0.0:
            z = np.linspace(0.0, Z, N_z)
            k_z = self.getKz(k_2, k_x, k_y)
            mu = mu_r2*constants.mu_0
        else:    
            z = np.linspace(Z, 0.0, N_z)
            k_z = -self.getKz(k_1, k_x, k_y)
            mu = mu_r1*constants.mu_0
        
        if n_y_repeat>0:
            D_y = y[-1]-y[0]+dy
            y_plus  = np.linspace(y[-1]+dy, y[-1]+dy+n_y_repeat*D_y, N_y*n_y_repeat, endpoint=False) 
            y_minus = np.linspace(y[0]-n_y_repeat*D_y, y[0], N_y*n_y_repeat, endpoint=False) 
            y = np.concatenate((y_minus, y, y_plus))
                
        
        E_f[2] = self.getEzf(E_f[0], E_f[1], k_x, k_y, k_z)
        
        
        Z, Y = np.meshgrid(z, y)
        Ex = np.zeros(Y.shape, dtype=complex)
        Ey = np.zeros(Y.shape, dtype=complex)
        Ez = np.zeros(Y.shape, dtype=complex)
        for j in range(N_y):
            k_y_j = k_y[j]
            k_z_j = k_z[j]
            k_dot_r = k_y_j*Y + k_z_j*Z
            Ex += E_f[0, j]*np.exp(-1j*k_dot_r)
            Ey += E_f[1, j]*np.exp(-1j*k_dot_r)
            Ez += E_f[2, j]*np.exp(-1j*k_dot_r)
        
        Ex = np.fft.fftshift(Ex, axes=(0,))
        Ey = np.fft.fftshift(Ey, axes=(0,))
        Ez = np.fft.fftshift(Ez, axes=(0,))

        if getH:
            H_f = np.zeros((3,N_y), dtype=complex)
            for j in range(N_y):
                k_vec_j = np.array([k_x, k_y[j], k_z[j]])
                H_f[:,j] = np.cross(k_vec_j, E_f[:,j])/(omega*mu)
        
            Hx = np.zeros(Y.shape, dtype=complex)
            Hy = np.zeros(Y.shape, dtype=complex)
            Hz = np.zeros(Y.shape, dtype=complex)
            for j in range(N_y):
                k_y_j = k_y[j]
                k_z_j = k_z[j]
                k_dot_r = k_y_j*Y + k_z_j*Z
                Hx += H_f[0, j]*np.exp(-1j*k_dot_r)
                Hy += H_f[1, j]*np.exp(-1j*k_dot_r)
                Hz += H_f[2, j]*np.exp(-1j*k_dot_r)
            
            Hx = np.fft.fftshift(Hx, axes=(0,))
            Hy = np.fft.fftshift(Hy, axes=(0,))
            Hz = np.fft.fftshift(Hz, axes=(0,))

            return [[Y, Z], [Ex, Ey, Ez], [Hx, Hy, Hz]]

        else:
            return [[Y, Z], [Ex, Ey, Ez]]
        
    
    
    
    
    
