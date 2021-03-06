##  TMM.py
##  Transmission Matrix Method


__all__ = ["TMM", "TMMGen", "TMM_EM_N", "TMM_EM_EPSMU", "TMM_EM_N_Time"]


from scipy import constants

import numpy as np
import math, cmath


""" e^(+j*omega*t) time convention
    WARNING: n = sqrt(eps_r*mu_r) ---> results presented here are valid for mu_r=1 only
    otherwise the formulation for reflection and transmission should be developed 
    using eps_r and mu_r instead of n
"""


'''
class TMM:
    def __init__(self, freq, theta=0.0):
        """ theta: angle with respect to the interface normal in rad
        """
        if theta!=0.0:
            print('Normal incidence only. Use TMM_EM_N instead.')
            raise ValueError('theta=0 only!')
            
        self.theta = theta
        self.SetFrequency(freq)
        return
        
    def SetFrequency(self, freq):
        self.freq = freq
        self.omega = 2.0*math.pi*freq
        self.k_0 = self.omega/constants.c
        self.k_0z = self.k_0*math.cos(self.theta)
        
    def SetupMultilayer(self, n_0, n_1, n_arr, d_arr):
        """ 
        """
        self.multilayer = [n_0, n_1, (n_arr, d_arr)]
        return
        
    def MaxwellInterfaceMatrix(self, n_1, n_2):
        """ Electric field transmission matrix
            a:forward wave  b: backward wave
            a1 + b1 = a2 + b2
            (a1 - b1)/n1 = (a2 - b2)/n2
        """
        TM = np.array([[n_1/(2.0*n_2) + 1.0/2.0, -n_1/(2.0*n_2) + 1.0/2.0], [-n_1/(2.0*n_2) + 1.0/2.0, n_1/(2.0*n_2) + 1.0/2.0]])
        return TM

    def MaxwellPropagMatrix(self, d, k_z):
        """ propagation matrix n:refractive index   d: thickness   k_z:propagation number (normal to interface)
            a2 = a1*exp(-j*k*d)
            b2 = b1*exp(+j*k*d)
        """
        PM = np.array([[cmath.exp(-1j*d*k_z), 0], [0, cmath.exp(1j*d*k_z)]])
        return PM
        
    def StaircaseContinuousProfile(self, f_n, n_slabs, z_0, z_1):
        """ It descretizes the refractive index profile defined by the function
            f_n to n_slabs slabs between z_0 and z_1 
            returns n_slabs slabs with a given thickness
        """
        z = np.linspace(z_0, z_1, n_slabs, endpoint=False)
        d = (z[1] - z[0])
        z = z + d/2.0
        
        n = f_n(z)
        d = np.ones(n_slabs)*d
        return [n, d]
        
    def RepeatProfile(self, n, d, n_repeat):
        ## n_repeat: final number of layers
        assert n_repeat > 0
        n_p = n.copy()
        d_p = d.copy()
        for i in range(n_repeat-1):
            n_p = np.concatenate((n_p, n))
            d_p = np.concatenate((d_p, d))
        return [n_p, d_p]
        
    
    def GetKz(self, n):
        return self.k_0z*n
    
    def GetTMTotal(self):
        """ get the total cascaded transmission matrix
            n_0: incident medium (start)
            n_arr__d: slabs in the middle and their thickness
            n_1: transmission medium (final)
        """
        n_0, n_1, n_arr__d = self.multilayer
        n_arr, d = n_arr__d
        len_arr = len(n_arr)
        assert len_arr>0
        IM = self.MaxwellInterfaceMatrix(n_0,n_arr[0])
        TM = IM
        
        for i in range(len(n_arr)-1):
            n_i = n_arr[i]
            k_z = self.GetKz(n_i)
            PM = self.MaxwellPropagMatrix(d[i], k_z)
            IM = self.MaxwellInterfaceMatrix(n_i, n_arr[i+1])
            TM = IM.dot(PM).dot(TM)
        

        PM = self.MaxwellPropagMatrix(d[len_arr-1], self.GetKz(n_arr[len_arr-1]))
        IM = self.MaxwellInterfaceMatrix(n_arr[len_arr-1], n_1)
        TM = IM.dot(PM).dot(TM)
        
        return TM
        
        
    def GetTransmissionReflection(self):
        TM = self.GetTMTotal()
        if TM.shape==(2, 2):
            # A_10*I + A_11*R = 0
            R = -TM[1, 0]/TM[1, 1]
            # A_00*I + A_01*R = T
            T = TM[0, 0] + R*TM[0, 1]
            return [T, R]
        else:
            return None
    
        
    def GetFieldPlot(self, n_pts_0, n_pts_i, n_pts_1, d0=None, d1=None, verbose=False):
        """ n_pts_i: number of points in each slab
            n_0: number of points in the incident (reflection) medium
            n_1: number of points in the transmission medium
        """
        n_0, n_1, n_arr__d = self.multilayer
        
        n_arr, d = n_arr__d
        
        lambda_0_z = 2.0*math.pi/self.k_0z
        if d0==None:
            d0 = 2.0*lambda_0_z
        if d1==None:
            d1 = 2.0*lambda_0_z
        
        T, R = self.GetTransmissionReflection()
        
        if verbose:
            print('T:', T, ' --- ', abs(T))
            print('R:', R, ' --- ', abs(R))
        
        n_slabs = len(n_arr) #number of slabs
        assert n_slabs>0
        if verbose:
            print('n_slabs: ', n_slabs)
        
        z_pts_0 = np.linspace(-d0, 0.0, n_pts_0, endpoint=False)
        z_start = 0.0
        z_pts_i = np.linspace(z_start, z_start+d[0], n_pts_i, endpoint=False)
        for i in range(1, len(d)):
            z_start += d[i-1]
            z_pts_i = np.concatenate((z_pts_i, np.linspace(z_start, z_start+d[i], n_pts_i, endpoint=False)))
        z_pts_1 = np.linspace(sum(d), sum(d)+d1,  n_pts_1, endpoint=True)
        
        k_z = self.GetKz(n_0)
        a_b_0 = np.array([1.0, R])
        a_b_i = a_b_0

        E_FWD_0 = a_b_i[0]*np.exp(-1j*k_z*(z_pts_0-z_pts_i[0]))
        E_BWD_0 = a_b_i[1]*np.exp(+1j*k_z*(z_pts_0-z_pts_i[0]))
        
        IM = self.MaxwellInterfaceMatrix(n_0,n_arr[0])
        TM = IM
        
        E_FWD_i = np.zeros(z_pts_i.shape, dtype=complex)
        E_BWD_i = np.zeros(z_pts_i.shape, dtype=complex)

        for i in range(n_slabs-1):
            n_i = n_arr[i]
            k_z = self.GetKz(n_i)
            a_b_i = TM.dot(a_b_0)
            
            
            E_FWD_i[i*n_pts_i: (i+1)*n_pts_i] = a_b_i[0]*np.exp(-1j*k_z*(z_pts_i[i*n_pts_i: (i+1)*n_pts_i] - z_pts_i[i*n_pts_i]))
            E_BWD_i[i*n_pts_i: (i+1)*n_pts_i] = a_b_i[1]*np.exp(+1j*k_z*(z_pts_i[i*n_pts_i: (i+1)*n_pts_i] - z_pts_i[i*n_pts_i]))
            
            
            PM = self.MaxwellPropagMatrix(d[i], k_z)
            IM = self.MaxwellInterfaceMatrix(n_i, n_arr[i+1])
            
            #print(E_FWD_i[(i+1)*n_pts_i-1], PM.dot(a_b_i)[0])
            #print(E_BWD_i[(i+1)*n_pts_i-1], PM.dot(a_b_i)[1])
            
            TM = IM.dot(PM).dot(TM)
        
        k_z = self.GetKz(n_arr[len(n_arr)-1])
        a_b_i = TM.dot(a_b_0)

        i = n_slabs-1
        E_FWD_i[i*n_pts_i: (i+1)*n_pts_i] = a_b_i[0]*np.exp(-1j*k_z*(z_pts_i[i*n_pts_i: (i+1)*n_pts_i] - z_pts_i[i*n_pts_i]))
        E_BWD_i[i*n_pts_i: (i+1)*n_pts_i] = a_b_i[1]*np.exp(+1j*k_z*(z_pts_i[i*n_pts_i: (i+1)*n_pts_i] - z_pts_i[i*n_pts_i]))
        
        PM = self.MaxwellPropagMatrix(d[len(d)-1], k_z)
        IM = self.MaxwellInterfaceMatrix(n_arr[len(n_arr)-1], n_1)
        TM = IM.dot(PM).dot(TM)
        
        k_z = self.GetKz(n_1)
        a_b_i = TM.dot(a_b_0)
        
        if verbose:
            print(a_b_i)

        E_FWD_1 = a_b_i[0]*np.exp(-1j*k_z*(z_pts_1-z_pts_1[0]))
        E_BWD_1 = a_b_i[1]*np.exp(+1j*k_z*(z_pts_1-z_pts_1[0]))
        
        z_pts = np.concatenate((z_pts_0, z_pts_i, z_pts_1))
        E_FWD = np.concatenate((E_FWD_0, E_FWD_i, E_FWD_1))
        E_BWD = np.concatenate((E_BWD_0, E_BWD_i, E_BWD_1))
        
        return [z_pts, E_FWD, E_BWD]
        
        
    def GetMediumPlot(self, n_pts_0, n_pts_i, n_pts_1, d0=None, d1=None):
        n_0, n_1, n_arr__d = self.multilayer
        
        n_arr, d = n_arr__d
        
        lambda_0_z = 2.0*math.pi/self.k_0z
        if d0==None:
            d0 = 2.0*lambda_0_z
        if d1==None:
            d1 = 2.0*lambda_0_z
        
        n_slabs = len(n_arr) #number of slabs
        assert n_slabs>0
        
        z_pts_0 = np.linspace(-d0, 0.0, n_pts_0, endpoint=False)
        z_start = 0.0
        z_pts_i = np.linspace(z_start, z_start+d[0], n_pts_i, endpoint=False)
        for i in range(1, len(d)):
            z_start += d[i-1]
            z_pts_i = np.concatenate((z_pts_i, np.linspace(z_start, z_start+d[i], n_pts_i, endpoint=False)))
        z_pts_1 = np.linspace(sum(d), sum(d)+d1,  n_pts_1, endpoint=True)
        
        
        n_0_vec = n_0*np.ones(n_pts_0)
        n_1_vec = n_1*np.ones(n_pts_1)
        
        n_i_vec = np.zeros(z_pts_i.shape, dtype=complex)

        for i in range(n_slabs):
            n_i = n_arr[i]
            n_i_vec[i*n_pts_i: (i+1)*n_pts_i] = n_i*np.ones(n_pts_i)

        z_pts = np.concatenate((z_pts_0, z_pts_i, z_pts_1))
        n_vec = np.concatenate((n_0_vec, n_i_vec, n_1_vec))
        
        return [z_pts, n_vec]
    
        
    def GetTransmissionReflectionFreqBand(self, f_0, f_1, N):
        f = np.linspace(f_0, f_1, N)
        T = np.zeros(N, dtype=complex)
        R = np.zeros(N, dtype=complex)
        for i in range(N):
            self.SetFrequency(f[i])
            T[i], R[i] = self.GetTransmissionReflection()
            
        return [f, T, R]
    

    def StaircaseGradientBarrier1D(self, n0, y, d, s1, s2, n_slabs):
        # d: width of the barrier
        from Electromagnetics.GMTM import GMTM1D
        
        gmtm = GMTM1D()
        gmtm.SetProfileParameters(n0, y, d, s1, s2)
        
        def f_n(z):
            return gmtm.GetProfile(z)
        
        z_0 = 0.0
        z_1 = d
        return self.StaircaseContinuousProfile(f_n, n_slabs, z_0, z_1)
        
'''


class TMMGen:
    """ General TMM method, the problem is modelled as a general cascading of
        layers defined by associated matrices
        the parameters and the associated matrices will be provided by user
    """
    def __init__(self, mats, params, n_var):
        """ params: parameters that define the matrices, it is a list of lists
                each list representing the parameters needed for the associated
                matrix as strings
            mats: a set of matrices, representing different types of layers
                for exapmle for a series of cascaded slabs for an electromagnetic 
                transmission/reflection problem there are two types of matrices, one
                representing the interface between 2 layers, the othr one representing
                the propagation matrix, the non-zero components of the matrices are 
                functions in terms of the parameters in the params list
            n_var: number of input/ouptt variables. The matrices have dimensions
                n_var*n_var
        """
        self.params = params
        self.mats = mats
        self.n_var = n_var


    def setLayers(self, layers):
        """ layers: list of lists describing each layer with elements
                [i, params] where params is the list of parameters required for
                matrix i
                i refers to a matrix in self.mats
                each layer element describes a matrix. For example in a photonic 
                crystal, interfaces and slabs are elements of the later
        """  
        self.layers = layers
        
        
    def getTotalTransmissionMatrix(self):
        """ multiplies all the matrices in self.layers and returns the resulting 
            matrix
        """
        layers = self.layers
        mats = self.mats
        N = self.n_var
        M_tot = np.eye(N, dtype=complex)
        for i in range(len(layers)):
            ind_mat, params = layers[i]
            
            M_i_func = mats[ind_mat]
            M_i = M_i_func(*tuple(params))
            
            M_tot = M_i.dot(M_tot)
        
        self.M_tot = M_tot
        return M_tot
        
        
    def setInputsOutputs(self, io_vals, io_vars):
        """ it presets some of the i/o variables and calculates the transmission
        matrix for the rest 
        io_vals: [['in'/'out', index, value], ...]
        io_vars: [['in'/'out', index], ...]  (index: associated index in the 
        input or output vectors 
        io_vars only sets the indices for variables not included in io_vals
        they both have length N, where N is the total number of inputs (outputs)  
        """
        ##TODO: some more conclusive tests to be done
        N = self.n_var
        assert len(io_vars)==N and len(io_vals)==N
        
        self.getTotalTransmissionMatrix()
        
        ## construct a matrix of all input and output variables 
        ## A_tot*[x_in, x_out]^T = 0
        A_tot = np.zeros((N, 2*N), dtype=complex)
        A_tot[0:N, 0:N] = self.M_tot
        A_tot[0:N, N:2*N] = -np.eye(N)
        
        ## assign new indices to all variables i.e. [x_in, x_out]
        vars_order = [-1]*(2*N)
        for i in range(len(io_vars)):
            io, ind = io_vars[i]
            if io=='in':
                vars_order[ind] = i
            else:
                assert io=='out'
                vars_order[ind+N] = i
        
        ind_last = N
        for i in range(len(vars_order)):
            if vars_order[i]<0:
                vars_order[i] = ind_last
                ind_last += 1
        assert ind_last==2*N
        
        ## get the permutation matrix  x_old = P*x_new
        P = np.zeros((2*N, 2*N))
        for i in range(2*N):
            P[i, vars_order[i]] = 1
        
        ## set A_tot_new = A_tot*P  --->   A_tot_new*x_new = 0
        A_tot_new = A_tot.dot(P)
        
        ## set A_new and B_new such that --> A_new*x_new[0:N] = B_new*x_new[N:2N]
        A_new = A_tot_new[0:N, 0:N]
        B_new = -A_tot_new[0:N, N:2*N]

        ## set y=x_new[N:2N]
        y = np.zeros(N, dtype=complex)
        for i in range(N):
            io, ind, val = io_vals[i]
            if io=='in':
                assert vars_order[ind]>=N
                y[vars_order[ind]-N] = val
            else:
                assert io=='out'
                assert vars_order[ind+N]>=N
                y[vars_order[ind+N]-N] = val
                
        ## solve for x = x_new[0:N]
        b = B_new.dot(y)
        #x = np.linalg.inv(A_new).dot(b)
        x = np.linalg.solve(A_new, b)
        return x


    def GetLayerByLayerOutputs(self, in_vals):
        """ it takes the input values and returns the output after each layer 
        (after each multiplication)
        in_vals: input values (numpy vector)
        """
        layers = self.layers
        mats = self.mats
        N = self.n_var
        M_tot = np.eye(N, dtype=complex)
        outs = [None]*len(layers)
        for i in range(len(layers)):
            ind_mat, params = layers[i]
            
            M_i_func = mats[ind_mat]
            M_i = M_i_func(*tuple(params))
            
            M_tot = M_i.dot(M_tot)
            outs[i] = M_tot.dot(in_vals)
            
        return outs




class TMM_EM_N:
    """ TMM for electromagnetic wave, described by refractive index
    
    for theta!=0.0:
        for TE, T and R refer to electric field transmission and reflection
        for TM, T and R refer to magnetic field transmission and reflection
    for theta==0.0:
        T and R refer to electric field transmission and reflection
    """
    
    def __init__(self, freq, theta=0.0, TETM='TE'):
        """ theta: angle with respect to the interface normal in rad
        """
        self.multilayer = None
        self.theta = theta
        self.SetFrequency(freq)
        assert TETM in ['TE', 'TM']
        self.TETM = TETM
        self.eta_0 = np.sqrt(constants.mu_0/constants.epsilon_0)
        self.k_imag_max = 100.0
        return
        
    def SetFrequency(self, freq):
        self.freq = freq
        self.omega = 2.0*math.pi*freq
        self.k_0 = self.omega/constants.c
        if self.multilayer!=None:
            self.SetKt()
        #self.k_0z = self.k_0*math.cos(self.theta)
        
    def SetupMultilayer(self, n_0, n_1, n_arr, d_arr):
        """ 
        """
        self.multilayer = [n_0, n_1, (n_arr, d_arr)]
        self.SetKt()
        self.SetupTMMGen()
        return
        
    def SetKt(self):
        n_0 = self.multilayer[0]
        k_t = self.k_0*n_0*math.sin(self.theta) 
        self.k_t = k_t
        #print('n_0:', n_0, ' k_t:', k_t, 'k_0:', self.k_0, 'k_z:', self.GetKz(n_0))
        
    def GetKz(self, n):
        k_0z_sq = self.k_0**2*n**2 - self.k_t**2
        k_0z = np.sqrt(k_0z_sq + 0j)
        imag_pos = (np.imag(k_0z)>=0)
        k_0z = k_0z*np.logical_not(imag_pos) + np.conjugate(k_0z)*imag_pos
        imag_big = np.abs(np.imag(k_0z))<self.k_imag_max
        k_0z = k_0z*imag_big + (np.real(k_0z)-1j*self.k_imag_max)*np.logical_not(imag_big)
        return k_0z
        
    def MaxwellInterfaceMatrix(self, n_1, n_2):
        """ Electric field transmission matrix
            a:forward wave  b: backward wave
            a1 + b1 = a2 + b2
            (a1 - b1)/n1 = (a2 - b2)/n2
        """
        #if self.theta==0.0:
        #    TM = np.array([[n_1/(2.0*n_2) + 1.0/2.0, -n_1/(2.0*n_2) + 1.0/2.0], [-n_1/(2.0*n_2) + 1.0/2.0, n_1/(2.0*n_2) + 1.0/2.0]])
        #    return TM
        if self.TETM=='TE':
            k_z1 = self.GetKz(n_1)
            k_z2 = self.GetKz(n_2)
            TM = np.array([[(k_z1 + k_z2)/(2*k_z2), (-k_z1 + k_z2)/(2*k_z2)], [(-k_z1 + k_z2)/(2*k_z2), (k_z1 + k_z2)/(2*k_z2)]])
            return TM
        else:
            assert self.TETM=='TM'
            k_z1 = self.GetKz(n_1)
            k_z2 = self.GetKz(n_2)
            TM = np.array([[(k_z1*n_2**2 + k_z2*n_1**2)/(2*k_z2*n_1**2), (-k_z1*n_2**2 + k_z2*n_1**2)/(2*k_z2*n_1**2)], \
                    [(-k_z1*n_2**2 + k_z2*n_1**2)/(2*k_z2*n_1**2), (k_z1*n_2**2 + k_z2*n_1**2)/(2*k_z2*n_1**2)]])
            return TM


    def MaxwellPropagMatrix(self, d, k_z):
        """ propagation matrix n:refractive index   d: thickness   k_z:propagation number (normal to interface)
            a2 = a1*exp(-j*k*d)
            b2 = b1*exp(+j*k*d)
        """
        PM = np.array([[np.exp(-1j*d*k_z), 0], [0, np.exp(1j*d*k_z)]])
        return PM
        
    def StaircaseContinuousProfile(self, f_n, n_slabs, z_0, z_1):
        """ It descretizes the refractive index profile defined by the function
            f_n to n_slabs slabs between z_0 and z_1 
            returns n_slabs slabs with a given thickness
        """
        z = np.linspace(z_0, z_1, n_slabs, endpoint=False)
        d = (z[1] - z[0])
        z = z + d/2.0
        
        n = f_n(z)
        d = np.ones(n_slabs)*d
        return [n, d]
        
    def RepeatProfile(self, n, d, n_repeat):
        ## n_repeat: final number of layers
        assert n_repeat > 0
        n_p = n.copy()
        d_p = d.copy()
        for i in range(n_repeat-1):
            n_p = np.concatenate((n_p, n))
            d_p = np.concatenate((d_p, d))
        return [n_p, d_p]
        
    
    def SetupTMMGen(self):
        """ get the total cascaded transmission matrix
            n_0: incident medium (start)
            n_arr__d: slabs in the middle and their thickness
            n_1: transmission medium (final)
        """
        n_0, n_1, n_arr__d = self.multilayer
        n_arr, d = n_arr__d
        n_slabs = len(n_arr)
        assert n_slabs>0
        
        
        param_interface = ['n_1', 'n_2']
        mat_interface = self.MaxwellInterfaceMatrix
        param_prop = ['d', 'k_z']
        mat_prop = self.MaxwellPropagMatrix
        
        params = [param_interface, param_prop]
        mats = [mat_interface, mat_prop]
        tmmgen = TMMGen(mats, params, n_var=2)
        
        mattype_interface = 0
        mattype_propagation = 1

        layers = [None]*(2*n_slabs+1)
        layers[0] = [mattype_interface, [n_0, n_arr[0]]]
        ind_layer = 1
        for i in range(n_slabs-1):
            k_zi = self.GetKz(n_arr[i])
            layers[ind_layer] = [mattype_propagation, [d[i], k_zi]]
            ind_layer += 1
            layers[ind_layer] = [mattype_interface, [n_arr[i], n_arr[i+1]]]
            ind_layer += 1

        i = n_slabs-1
        k_zi = self.GetKz(n_arr[i])
        layers[ind_layer] = [mattype_propagation, [d[i], k_zi]]
        ind_layer += 1
        layers[ind_layer] = [mattype_interface, [n_arr[i], n_1]]
        ind_layer += 1
        
        tmmgen.setLayers(layers)
        self.tmmgen = tmmgen

    def GetTMTotal(self):
        """ get the total cascaded transmission matrix
            n_0: incident medium (start)
            n_arr__d: slabs in the middle and their thickness
            n_1: transmission medium (final)
        """
        self.SetupTMMGen()
        TM = self.tmmgen.getTotalTransmissionMatrix()     
        return TM
        
        
    def GetTransmissionReflection(self):
        self.SetupTMMGen()
        io_vals = [['in', 0, 1.0], ['out', 1, 0.0]]
        io_vars = [['in', 1], ['out', 0]]
        [R, T] = self.tmmgen.setInputsOutputs(io_vals, io_vars)

        return [T, R]    
        
        
    def SetGtoEH(self, EH, k_z, n, G_FWD, G_BWD):
        eta = self.eta_0/n
        a_kz = k_z/(self.k_0*n)
        a_kt = self.k_t/(self.k_0*n)
        if self.TETM=='TE': 
            if EH=='Et':
                pass
            elif EH=='Ez':
                G_FWD *= 0.0
                G_BWD *= 0.0
            elif EH=='Ht':
                G_FWD *= a_kz/eta 
                G_BWD *= -a_kz/eta
            elif EH=='Hz':
                G_FWD *= -a_kt/eta 
                G_BWD *= -a_kt/eta
        else:
            assert self.TETM=='TM' 
            if EH=='Et':
                G_FWD *= -a_kz*eta 
                G_BWD *= a_kz*eta
            elif EH=='Ez':
                G_FWD *= a_kt*eta 
                G_BWD *= a_kt*eta
            elif EH=='Ht':
                pass
            elif EH=='Hz':
                G_FWD *= 0.0
                G_BWD *= 0.0
        
    def GetFieldPlot(self, n_pts_0, n_pts_i, n_pts_1, d0=None, d1=None, EH='Et', verbose=False):
        """ n_pts_i: number of points in each slab
            n_0: number of points in the incident (reflection) medium
            n_1: number of points in the transmission medium
        """
        assert EH in ['Et', 'Ht', 'Ez', 'Hz']
        n_0, n_1, n_arr__d = self.multilayer
        
        n_arr, d = n_arr__d
        
        assert np.real(self.GetKz(n_0))!=0.0
        lambda_0_z = 2.0*math.pi/np.real(self.GetKz(n_0))
        if d0==None:
            d0 = 2.0*lambda_0_z
        if d1==None:
            d1 = 2.0*lambda_0_z
        #print('d0: ', d0, 'n_0:', n_0, '  k_z:', self.GetKz(n_0), 'k_t:', self.k_t, \
        #    'k_0:', self.k_0, ' lambda_0_z: ', lambda_0_z)
        
        T, R = self.GetTransmissionReflection()
        
        in_vals = np.array([1.0, R])
        layers_outpt = self.tmmgen.GetLayerByLayerOutputs(in_vals)
        
        if verbose:
            print('T:', T, ' --- ', abs(T))
            print('R:', R, ' --- ', abs(R))
        
        n_slabs = len(n_arr) #number of slabs
        assert n_slabs>0
        if verbose:
            print('n_slabs: ', n_slabs)
        
        z_pts_0 = np.linspace(-d0, 0.0, n_pts_0, endpoint=False)
        z_start = 0.0
        z_pts_i = np.linspace(z_start, z_start+d[0], n_pts_i, endpoint=False)
        for i in range(1, len(d)):
            z_start += d[i-1]
            z_pts_i = np.concatenate((z_pts_i, np.linspace(z_start, z_start+d[i], n_pts_i, endpoint=False)))
        z_pts_1 = np.linspace(sum(d), sum(d)+d1,  n_pts_1, endpoint=True)
        
        
        # G: E or H TODO: set G  t: x
        k_z = self.GetKz(n_0)
        G_FWD_0 = in_vals[0]*np.exp(-1j*k_z*(z_pts_0-z_pts_i[0]))
        G_BWD_0 = in_vals[1]*np.exp(+1j*k_z*(z_pts_0-z_pts_i[0]))
        self.SetGtoEH(EH, k_z, n_0, G_FWD_0, G_BWD_0)
        
        G_FWD_i = np.zeros(z_pts_i.shape, dtype=complex)
        G_BWD_i = np.zeros(z_pts_i.shape, dtype=complex)

        for i in range(n_slabs):
            n_i = n_arr[i]
            k_z = self.GetKz(n_i)
            a_b_i = layers_outpt[2*i]
            
            
            G_FWD_i[i*n_pts_i: (i+1)*n_pts_i] = a_b_i[0]*np.exp(-1j*k_z*(z_pts_i[i*n_pts_i: (i+1)*n_pts_i] - z_pts_i[i*n_pts_i]))
            G_BWD_i[i*n_pts_i: (i+1)*n_pts_i] = a_b_i[1]*np.exp(+1j*k_z*(z_pts_i[i*n_pts_i: (i+1)*n_pts_i] - z_pts_i[i*n_pts_i]))
            self.SetGtoEH(EH, k_z, n_i, G_FWD_i[i*n_pts_i: (i+1)*n_pts_i], G_BWD_i[i*n_pts_i: (i+1)*n_pts_i])
            
                
        k_z = self.GetKz(n_1)
        a_b_i = layers_outpt[-1]
        
        if verbose:
            print(a_b_i)

        G_FWD_1 = a_b_i[0]*np.exp(-1j*k_z*(z_pts_1-z_pts_1[0]))
        G_BWD_1 = a_b_i[1]*np.exp(+1j*k_z*(z_pts_1-z_pts_1[0]))
        self.SetGtoEH(EH, k_z, n_1, G_FWD_1, G_BWD_1)
        
        z_pts = np.concatenate((z_pts_0, z_pts_i, z_pts_1))
        G_FWD = np.concatenate((G_FWD_0, G_FWD_i, G_FWD_1))
        G_BWD = np.concatenate((G_BWD_0, G_BWD_i, G_BWD_1))
        
        return [z_pts, G_FWD, G_BWD]
        
        
    def GetMediumPlot(self, n_pts_0, n_pts_i, n_pts_1, d0=None, d1=None):
        n_0, n_1, n_arr__d = self.multilayer
        
        n_arr, d = n_arr__d
        
        lambda_0_z = 2.0*math.pi/self.GetKz(n_0)
        if d0==None:
            d0 = 2.0*lambda_0_z
        if d1==None:
            d1 = 2.0*lambda_0_z
        
        n_slabs = len(n_arr) #number of slabs
        assert n_slabs>0
        
        z_pts_0 = np.linspace(-d0, 0.0, n_pts_0, endpoint=False)
        z_start = 0.0
        z_pts_i = np.linspace(z_start, z_start+d[0], n_pts_i, endpoint=False)
        for i in range(1, len(d)):
            z_start += d[i-1]
            z_pts_i = np.concatenate((z_pts_i, np.linspace(z_start, z_start+d[i], n_pts_i, endpoint=False)))
        z_pts_1 = np.linspace(sum(d), sum(d)+d1,  n_pts_1, endpoint=True)
        
        
        n_0_vec = n_0*np.ones(n_pts_0)
        n_1_vec = n_1*np.ones(n_pts_1)
        
        n_i_vec = np.zeros(z_pts_i.shape, dtype=complex)

        for i in range(n_slabs):
            n_i = n_arr[i]
            n_i_vec[i*n_pts_i: (i+1)*n_pts_i] = n_i*np.ones(n_pts_i)

        z_pts = np.concatenate((z_pts_0, z_pts_i, z_pts_1))
        n_vec = np.concatenate((n_0_vec, n_i_vec, n_1_vec))
        
        return [z_pts, n_vec]
    
        
    def GetTransmissionReflectionFreqBand(self, f_0, f_1, N):
        f = np.linspace(f_0, f_1, N)
        T = np.zeros(N, dtype=complex)
        R = np.zeros(N, dtype=complex)
        for i in range(N):
            self.SetFrequency(f[i])
            self.SetupTMMGen()
            T[i], R[i] = self.GetTransmissionReflection()
            
        return [f, T, R]
    

    def StaircaseGradientBarrier1D(self, n0, y, d, s1, s2, n_slabs):
        # d: width of the barrier
        from Electromagnetics.GMTM import GMTM1D
        
        gmtm = GMTM1D()
        gmtm.SetProfileParameters(n0, y, d, s1, s2)
        
        def f_n(z):
            return gmtm.GetProfile(z)
        
        z_0 = 0.0
        z_1 = d
        return self.StaircaseContinuousProfile(f_n, n_slabs, z_0, z_1)


TMM = TMM_EM_N



class TMM_EM_EPSMU:
    """ TMM for electromagnetic wave, described by eps_r and mu_r
    
    for theta!=0.0:
        for TE, T and R refer to electric field transmission and reflection
        for TM, T and R refer to magnetic field transmission and reflection
    for theta==0.0:
        T and R refer to electric field transmission and reflection
    """
    
    def __init__(self, freq, theta=0.0, TETM='TE'):
        """ theta: angle with respect to the interface normal in rad
        """
        self.multilayer = None
        self.theta = theta
        self.SetFrequency(freq)
        assert TETM in ['TE', 'TM']
        self.TETM = TETM
        self.eta_0 = np.sqrt(constants.mu_0/constants.epsilon_0)
        self.k_imag_max = 100.0
        return
        
    def SetFrequency(self, freq):
        self.freq = freq
        self.omega = 2.0*math.pi*freq
        self.k_0 = self.omega/constants.c
        if self.multilayer!=None:
            self.SetKt()
        #self.k_0z = self.k_0*math.cos(self.theta)
        
    def SetupMultilayer(self, epsmu_0, epsmu_1, eps_arr, mu_arr, d_arr):
        """ 
        """
        self.multilayer = [epsmu_0, epsmu_1, (eps_arr, mu_arr, d_arr)]
        self.SetKt()
        self.SetupTMMGen()
        return
      
    def GetSqrtEpsMu(self, eps, mu):
        n_0 = np.sqrt(eps*mu)
        if np.real(eps)<0.0 and np.real(mu)<0.0:
            if np.real(n_0)>0.0:
                n_0 *= -1.0
        if np.imag(n_0)>0.0:
            n_0 = np.conjugate(n_0)
        return n_0
        
    def SetKt(self):
        eps_0, mu_0 = self.multilayer[0]
        n_0 = self.GetSqrtEpsMu(eps_0, mu_0)
        k_t = self.k_0*n_0*math.sin(self.theta) 
        self.k_t = k_t
        #print('n_0:', n_0, ' k_t:', k_t, 'k_0:', self.k_0, 'k_z:', self.GetKz(n_0))
        
    def GetKz(self, eps, mu):
        n__2 = eps*mu
        k_0z_sq = self.k_0**2*n__2 - self.k_t**2
        k_0z = np.sqrt(k_0z_sq + 0j)
        imag_pos = (np.imag(k_0z)>=0)
        k_0z = k_0z*np.logical_not(imag_pos) + np.conjugate(k_0z)*imag_pos
        imag_big = np.abs(np.imag(k_0z))<self.k_imag_max
        k_0z = k_0z*imag_big + (np.real(k_0z)-1j*self.k_imag_max)*np.logical_not(imag_big)
        return k_0z
        
    def MaxwellInterfaceMatrix(self, eps_1, mu_1, eps_2, mu_2):
        """ Electric field transmission matrix
            a:forward wave  b: backward wave
            a1 + b1 = a2 + b2
            (a1 - b1)/n1 = (a2 - b2)/n2
        """
        if self.TETM=='TE':
            k_z1 = self.GetKz(eps_1, mu_1)
            k_z2 = self.GetKz(eps_2, mu_2)
            TM = np.array([[(mu_1*k_z2 + mu_2*k_z1)/(2*mu_1*k_z2), (mu_1*k_z2 - mu_2*k_z1)/(2*mu_1*k_z2)], \
                    [(mu_1*k_z2 - mu_2*k_z1)/(2*mu_1*k_z2), (mu_1*k_z2 + mu_2*k_z1)/(2*mu_1*k_z2)]])
            return TM
        else:
            assert self.TETM=='TM'
            k_z1 = self.GetKz(eps_1, mu_1)
            k_z2 = self.GetKz(eps_2, mu_2)
            TM = np.array([[(eps_1*k_z2 + eps_2*k_z1)/(2*eps_1*k_z2), \
                (eps_1*k_z2 - eps_2*k_z1)/(2*eps_1*k_z2)], \
                [(eps_1*k_z2 - eps_2*k_z1)/(2*eps_1*k_z2), \
                (eps_1*k_z2 + eps_2*k_z1)/(2*eps_1*k_z2)]])
            return TM


    def MaxwellPropagMatrix(self, d, k_z):
        """ propagation matrix n:refractive index   d: thickness   k_z:propagation number (normal to interface)
            a2 = a1*exp(-j*k*d)
            b2 = b1*exp(+j*k*d)
        """
        PM = np.array([[np.exp(-1j*d*k_z), 0], [0, np.exp(1j*d*k_z)]])
        return PM
        
    def StaircaseContinuousProfile(self, f_n, n_slabs, z_0, z_1):
        """ It descretizes the refractive index profile defined by the function
            f_n to n_slabs slabs between z_0 and z_1 
            returns n_slabs slabs with a given thickness
        """
        z = np.linspace(z_0, z_1, n_slabs, endpoint=False)
        d = (z[1] - z[0])
        z = z + d/2.0
        
        n = f_n(z)
        d = np.ones(n_slabs)*d
        return [n, d]
        
    def RepeatProfile(self, eps, mu, d, n_repeat):
        ## n_repeat: final number of layers
        assert n_repeat > 0
        eps_p = eps.copy()
        mu_p = mu.copy()
        d_p = d.copy()
        for i in range(n_repeat-1):
            eps_p = np.concatenate((eps_p, eps))
            mu_p = np.concatenate((mu_p, mu))
            d_p = np.concatenate((d_p, d))
        return [eps_p, mu_p, d_p]
        
    
    def SetupTMMGen(self):
        """ get the total cascaded transmission matrix
            n_0: incident medium (start)
            n_arr__d: slabs in the middle and their thickness
            n_1: transmission medium (final)
        """
        epsmu_0, epsmu_1, epsmu_arr__d = self.multilayer
        eps_0, mu_0 = epsmu_0
        eps_1, mu_1 = epsmu_1
        eps_arr, mu_arr, d = epsmu_arr__d
        n_slabs = len(eps_arr)
        assert n_slabs>0
        
        
        param_interface = ['eps_1', 'mu_1', 'eps_2', 'mu_2']
        mat_interface = self.MaxwellInterfaceMatrix
        param_prop = ['d', 'k_z']
        mat_prop = self.MaxwellPropagMatrix
        
        params = [param_interface, param_prop]
        mats = [mat_interface, mat_prop]
        tmmgen = TMMGen(mats, params, n_var=2)
        
        mattype_interface = 0
        mattype_propagation = 1

        layers = [None]*(2*n_slabs+1)
        layers[0] = [mattype_interface, [eps_0, mu_0, eps_arr[0], mu_arr[0]]]
        ind_layer = 1
        for i in range(n_slabs-1):
            k_zi = self.GetKz(eps_arr[i], mu_arr[i])
            layers[ind_layer] = [mattype_propagation, [d[i], k_zi]]
            ind_layer += 1
            layers[ind_layer] = [mattype_interface, [eps_arr[i], mu_arr[i], eps_arr[i+1], mu_arr[i+1]]]
            ind_layer += 1

        i = n_slabs-1
        k_zi = self.GetKz(eps_arr[i], mu_arr[i])
        layers[ind_layer] = [mattype_propagation, [d[i], k_zi]]
        ind_layer += 1
        layers[ind_layer] = [mattype_interface, [eps_arr[i], mu_arr[i], eps_1, mu_1]]
        ind_layer += 1
        
        tmmgen.setLayers(layers)
        self.tmmgen = tmmgen

    def GetTMTotal(self):
        """ get the total cascaded transmission matrix
            n_0: incident medium (start)
            n_arr__d: slabs in the middle and their thickness
            n_1: transmission medium (final)
        """
        self.SetupTMMGen()
        TM = self.tmmgen.getTotalTransmissionMatrix()     
        return TM
        
        
    def GetTransmissionReflection(self):
        self.SetupTMMGen()
        io_vals = [['in', 0, 1.0], ['out', 1, 0.0]]
        io_vars = [['in', 1], ['out', 0]]
        [R, T] = self.tmmgen.setInputsOutputs(io_vals, io_vars)

        return [T, R]    
        
        
    def SetGtoEH(self, EH, k_z, eps, mu, G_FWD, G_BWD):
        ##TODO: confirm the sqrt imaginary sign
        eta = self.eta_0*np.sqrt(mu/eps)
        a_kz = k_z/(self.k_0*n)
        a_kt = self.k_t/(self.k_0*n)
        if self.TETM=='TE': 
            if EH=='Et':
                pass
            elif EH=='Ez':
                G_FWD *= 0.0
                G_BWD *= 0.0
            elif EH=='Ht':
                G_FWD *= a_kz/eta 
                G_BWD *= -a_kz/eta
            elif EH=='Hz':
                G_FWD *= -a_kt/eta 
                G_BWD *= -a_kt/eta
        else:
            assert self.TETM=='TM' 
            if EH=='Et':
                G_FWD *= -a_kz*eta 
                G_BWD *= a_kz*eta
            elif EH=='Ez':
                G_FWD *= a_kt*eta 
                G_BWD *= a_kt*eta
            elif EH=='Ht':
                pass
            elif EH=='Hz':
                G_FWD *= 0.0
                G_BWD *= 0.0
        
    def GetFieldPlot(self, n_pts_0, n_pts_i, n_pts_1, d0=None, d1=None, EH='Et', verbose=False):
        """ n_pts_i: number of points in each slab
            n_0: number of points in the incident (reflection) medium
            n_1: number of points in the transmission medium
        """
        assert EH in ['Et', 'Ht', 'Ez', 'Hz']
        epsmu_0, epsmu_1, epsmu_arr__d = self.multilayer
        eps_0, mu_0 = epsmu_0
        eps_1, mu_1 = epsmu_1
        eps_arr, mu_arr, d = epsmu_arr__d
        
        assert np.real(self.GetKz(eps_0, mu_0))!=0.0
        lambda_0_z = 2.0*math.pi/np.real(self.GetKz(eps_0, mu_0))
        if d0==None:
            d0 = 2.0*lambda_0_z
        if d1==None:
            d1 = 2.0*lambda_0_z
        #print('d0: ', d0, 'n_0:', n_0, '  k_z:', self.GetKz(n_0), 'k_t:', self.k_t, \
        #    'k_0:', self.k_0, ' lambda_0_z: ', lambda_0_z)
        
        T, R = self.GetTransmissionReflection()
        
        in_vals = np.array([1.0, R])
        layers_outpt = self.tmmgen.GetLayerByLayerOutputs(in_vals)
        
        if verbose:
            print('T:', T, ' --- ', abs(T))
            print('R:', R, ' --- ', abs(R))
        
        n_slabs = len(eps_arr) #number of slabs
        assert n_slabs>0
        if verbose:
            print('n_slabs: ', n_slabs)
        
        z_pts_0 = np.linspace(-d0, 0.0, n_pts_0, endpoint=False)
        z_start = 0.0
        z_pts_i = np.linspace(z_start, z_start+d[0], n_pts_i, endpoint=False)
        for i in range(1, len(d)):
            z_start += d[i-1]
            z_pts_i = np.concatenate((z_pts_i, np.linspace(z_start, z_start+d[i], n_pts_i, endpoint=False)))
        z_pts_1 = np.linspace(sum(d), sum(d)+d1,  n_pts_1, endpoint=True)
        
        
        # G: E or H TODO: set G  t: x
        k_z = self.GetKz(eps_0, mu_0)
        G_FWD_0 = in_vals[0]*np.exp(-1j*k_z*(z_pts_0-z_pts_i[0]))
        G_BWD_0 = in_vals[1]*np.exp(+1j*k_z*(z_pts_0-z_pts_i[0]))
        self.SetGtoEH(EH, k_z, eps_0, mu_0, G_FWD_0, G_BWD_0)
        
        G_FWD_i = np.zeros(z_pts_i.shape, dtype=complex)
        G_BWD_i = np.zeros(z_pts_i.shape, dtype=complex)

        for i in range(n_slabs):
            eps_i = eps_arr[i]
            mu_i = mu_arr[i]
            k_z = self.GetKz(eps_i, mu_i)
            a_b_i = layers_outpt[2*i]
            
            
            G_FWD_i[i*n_pts_i: (i+1)*n_pts_i] = a_b_i[0]*np.exp(-1j*k_z*(z_pts_i[i*n_pts_i: (i+1)*n_pts_i] - z_pts_i[i*n_pts_i]))
            G_BWD_i[i*n_pts_i: (i+1)*n_pts_i] = a_b_i[1]*np.exp(+1j*k_z*(z_pts_i[i*n_pts_i: (i+1)*n_pts_i] - z_pts_i[i*n_pts_i]))
            self.SetGtoEH(EH, k_z, eps_i, mu_i, G_FWD_i[i*n_pts_i: (i+1)*n_pts_i], G_BWD_i[i*n_pts_i: (i+1)*n_pts_i])
            
                
        k_z = self.GetKz(eps_1, mu_1)
        a_b_i = layers_outpt[-1]
        
        if verbose:
            print(a_b_i)

        G_FWD_1 = a_b_i[0]*np.exp(-1j*k_z*(z_pts_1-z_pts_1[0]))
        G_BWD_1 = a_b_i[1]*np.exp(+1j*k_z*(z_pts_1-z_pts_1[0]))
        self.SetGtoEH(EH, k_z, eps_1, mu_1, G_FWD_1, G_BWD_1)
        
        z_pts = np.concatenate((z_pts_0, z_pts_i, z_pts_1))
        G_FWD = np.concatenate((G_FWD_0, G_FWD_i, G_FWD_1))
        G_BWD = np.concatenate((G_BWD_0, G_BWD_i, G_BWD_1))
        
        return [z_pts, G_FWD, G_BWD]
        
        
    def GetMediumPlot(self, n_pts_0, n_pts_i, n_pts_1, d0=None, d1=None):
        epsmu_0, epsmu_1, epsmu_arr__d = self.multilayer
        eps_0, mu_0 = epsmu_0
        eps_1, mu_1 = epsmu_1
        eps_arr, mu_arr, d = epsmu_arr__d
        
        lambda_0_z = 2.0*math.pi/self.GetKz(eps_0, mu_0)
        if d0==None:
            d0 = 2.0*lambda_0_z
        if d1==None:
            d1 = 2.0*lambda_0_z
        
        n_slabs = len(eps_arr) #number of slabs
        assert n_slabs>0
        
        z_pts_0 = np.linspace(-d0, 0.0, n_pts_0, endpoint=False)
        z_start = 0.0
        z_pts_i = np.linspace(z_start, z_start+d[0], n_pts_i, endpoint=False)
        for i in range(1, len(d)):
            z_start += d[i-1]
            z_pts_i = np.concatenate((z_pts_i, np.linspace(z_start, z_start+d[i], n_pts_i, endpoint=False)))
        z_pts_1 = np.linspace(sum(d), sum(d)+d1,  n_pts_1, endpoint=True)
        
        
        eps_0_vec = eps_0*np.ones(n_pts_0)
        eps_1_vec = eps_1*np.ones(n_pts_1)
        mu_0_vec = mu_0*np.ones(n_pts_0)
        mu_1_vec = mu_1*np.ones(n_pts_1)
        
        eps_i_vec = np.zeros(z_pts_i.shape, dtype=complex)
        mu_i_vec = np.zeros(z_pts_i.shape, dtype=complex)

        for i in range(n_slabs):
            eps_i = eps_arr[i]
            eps_i_vec[i*n_pts_i: (i+1)*n_pts_i] = eps_i*np.ones(n_pts_i)
            mu_i = mu_arr[i]
            mu_i_vec[i*n_pts_i: (i+1)*n_pts_i] = mu_i*np.ones(n_pts_i)

        z_pts = np.concatenate((z_pts_0, z_pts_i, z_pts_1))
        eps_vec = np.concatenate((eps_0_vec, eps_i_vec, eps_1_vec))
        mu_vec = np.concatenate((mu_0_vec, mu_i_vec, mu_1_vec))
        
        return [z_pts, eps_vec, mu_vec]
    
        
    def GetTransmissionReflectionFreqBand(self, f_0, f_1, N):
        f = np.linspace(f_0, f_1, N)
        T = np.zeros(N, dtype=complex)
        R = np.zeros(N, dtype=complex)
        for i in range(N):
            self.SetFrequency(f[i])
            self.SetupTMMGen()
            T[i], R[i] = self.GetTransmissionReflection()
            
        return [f, T, R]
    



class TMM_EM_N_Time:
    """ TMM for electromagnetic wave in time varying media, described by 
    refractive index
    """
    
    def __init__(self, k):
        """ theta: angle with respect to the interface normal in rad
        """
        self.SetK(k)
        return
        
    def SetK(self, k):
        self.k = k
        self.omega_0 = self.k*constants.c
        
    def SetupMultilayer(self, n_0, n_1, n_arr, t_arr):
        """ 
        """
        self.multilayer = [n_0, n_1, (n_arr, t_arr)]
        self.SetupTMMGen()
        return
        
    def MaxwellInterfaceMatrix(self, n_1, n_2):
        """ Electric field transmission matrix
            a:forward wave  b: backward wave
            a1 + b1 = a2 + b2
            (a1 - b1)/n1 = (a2 - b2)/n2
        """
        ##TODO: to be modified for oblique incidence
        TM = np.array([[(n_1 + n_2)/(2.0*n_1), (n_1 - n_2)/(2.0*n_1)], [(n_1 - n_2)/(2.0*n_1), (n_1 + n_2)/(2.0*n_1)]])
        return TM

    def MaxwellPropagMatrix(self, T, omega):
        """ propagation matrix n:refractive index   d: thickness   k_z:propagation number (normal to interface)
            a2 = a1*exp(+j*omega*T)
            b2 = b1*exp(+j*omega*T)
        """
        PM = np.array([[np.exp(1j*T*omega), 0], [0, np.exp(-1j*T*omega)]])
        return PM

    def StaircaseContinuousProfile(self, f_n, n_slabs, z_0, z_1):
        """ It descretizes the refractive index profile defined by the function
            f_n to n_slabs slabs between z_0 and z_1 
            returns n_slabs slabs with a given thickness
        """
        z = np.linspace(z_0, z_1, n_slabs, endpoint=False)
        d = (z[1] - z[0])
        z = z + d/2.0
        
        n = f_n(z)
        d = np.ones(n_slabs)*d
        return [n, d]
        
    def RepeatProfile(self, n, d, n_repeat):
        ## n_repeat: final number of layers
        assert n_repeat > 0
        n_p = n.copy()
        d_p = d.copy()
        for i in range(n_repeat-1):
            n_p = np.concatenate((n_p, n))
            d_p = np.concatenate((d_p, d))
        return [n_p, d_p]
        
    
    def GetOmega(self, n):
        return self.omega_0/n
    
    def SetupTMMGen(self):
        """ get the total cascaded transmission matrix
            n_0: incident medium (start)
            n_arr__d: slabs in the middle and their thickness
            n_1: transmission medium (final)
        """
        n_0, n_1, n_arr__d = self.multilayer
        n_arr, d = n_arr__d
        n_slabs = len(n_arr)
        assert n_slabs>0
        
        
        param_interface = ['n_1', 'n_2']
        mat_interface = self.MaxwellInterfaceMatrix
        param_prop = ['T', 'omega']
        mat_prop = self.MaxwellPropagMatrix
        
        params = [param_interface, param_prop]
        mats = [mat_interface, mat_prop]
        tmmgen = TMMGen(mats, params, n_var=2)
        
        mattype_interface = 0
        mattype_propagation = 1

        layers = [None]*(2*n_slabs+1)
        layers[0] = [mattype_interface, [n_0, n_arr[0]]]
        ind_layer = 1
        for i in range(n_slabs-1):
            omega_i = self.GetOmega(n_arr[i])
            layers[ind_layer] = [mattype_propagation, [d[i], omega_i]]
            ind_layer += 1
            layers[ind_layer] = [mattype_interface, [n_arr[i], n_arr[i+1]]]
            ind_layer += 1

        i = n_slabs-1
        omega_i = self.GetOmega(n_arr[i])
        layers[ind_layer] = [mattype_propagation, [d[i], omega_i]]
        ind_layer += 1
        layers[ind_layer] = [mattype_interface, [n_arr[i], n_1]]
        ind_layer += 1
        
        tmmgen.setLayers(layers)
        self.tmmgen = tmmgen

    def GetTMTotal(self):
        """ get the total cascaded transmission matrix
            n_0: incident medium (start)
            n_arr__d: slabs in the middle and their thickness
            n_1: transmission medium (final)
        """
        self.SetupTMMGen()
        TM = self.tmmgen.getTotalTransmissionMatrix()     
        return TM
        
        
    def GetTransmissionPlusMinus(self):
        self.SetupTMMGen()
        io_vals = [['in', 0, 1.0], ['in', 1, 0.0]]
        io_vars = [['out', 0], ['out', 1]]
        [Tp, Tm] = self.tmmgen.setInputsOutputs(io_vals, io_vars)

        return [Tp, Tm]    
        
    def GetFieldPlot(self, n_pts_0, n_pts_i, n_pts_1, d0=None, d1=None, verbose=False):
        """ n_pts_i: number of points in each slab
            n_0: number of points in the incident (reflection) medium
            n_1: number of points in the transmission medium
        """
        n_0, n_1, n_arr__d = self.multilayer
        
        n_arr, d = n_arr__d
        
        T_0 = 2.0*math.pi/self.omega_0
        if d0==None:
            d0 = 2.0*T_0
        if d1==None:
            d1 = 2.0*T_0
        
        Tp, Tm = self.GetTransmissionPlusMinus()
        
        in_vals = np.array([1.0, 0.0])
        layers_outpt = self.tmmgen.GetLayerByLayerOutputs(in_vals)
        
        if verbose:
            print('T+:', Tp, ' --- ', abs(Tp))
            print('T-:', Tm, ' --- ', abs(Tm))
        
        n_slabs = len(n_arr) #number of slabs
        assert n_slabs>0
        if verbose:
            print('n_slabs: ', n_slabs)
        
        z_pts_0 = np.linspace(-d0, 0.0, n_pts_0, endpoint=False)
        z_start = 0.0
        z_pts_i = np.linspace(z_start, z_start+d[0], n_pts_i, endpoint=False)
        for i in range(1, len(d)):
            z_start += d[i-1]
            z_pts_i = np.concatenate((z_pts_i, np.linspace(z_start, z_start+d[i], n_pts_i, endpoint=False)))
        z_pts_1 = np.linspace(sum(d), sum(d)+d1,  n_pts_1, endpoint=True)
        
        
        omega = self.GetOmega(n_0)
        E_FWD_0 = in_vals[0]*np.exp(+1j*omega*(z_pts_0-z_pts_i[0]))
        E_BWD_0 = in_vals[1]*np.exp(-1j*omega*(z_pts_0-z_pts_i[0]))
        
        E_FWD_i = np.zeros(z_pts_i.shape, dtype=complex)
        E_BWD_i = np.zeros(z_pts_i.shape, dtype=complex)

        for i in range(n_slabs):
            n_i = n_arr[i]
            omega_i = self.GetOmega(n_i)
            a_b_i = layers_outpt[2*i]
            
            E_FWD_i[i*n_pts_i: (i+1)*n_pts_i] = a_b_i[0]*np.exp(+1j*omega_i*(z_pts_i[i*n_pts_i: (i+1)*n_pts_i] - z_pts_i[i*n_pts_i]))
            E_BWD_i[i*n_pts_i: (i+1)*n_pts_i] = a_b_i[1]*np.exp(-1j*omega_i*(z_pts_i[i*n_pts_i: (i+1)*n_pts_i] - z_pts_i[i*n_pts_i]))
            
                
        omega = self.GetOmega(n_1)
        a_b_i = layers_outpt[-1]
        
        if verbose:
            print(a_b_i)

        E_FWD_1 = a_b_i[0]*np.exp(+1j*omega*(z_pts_1-z_pts_1[0]))
        E_BWD_1 = a_b_i[1]*np.exp(-1j*omega*(z_pts_1-z_pts_1[0]))
        
        z_pts = np.concatenate((z_pts_0, z_pts_i, z_pts_1))
        E_FWD = np.concatenate((E_FWD_0, E_FWD_i, E_FWD_1))
        E_BWD = np.concatenate((E_BWD_0, E_BWD_i, E_BWD_1))
        
        return [z_pts, E_FWD, E_BWD]
        
        
    def GetMediumPlot(self, n_pts_0, n_pts_i, n_pts_1, d0=None, d1=None):
        n_0, n_1, n_arr__d = self.multilayer
        
        n_arr, d = n_arr__d
        
        T_0 = 2.0*math.pi/self.omega_0
        if d0==None:
            d0 = 2.0*T_0
        if d1==None:
            d1 = 2.0*T_0
        
        n_slabs = len(n_arr) #number of slabs
        assert n_slabs>0
        
        z_pts_0 = np.linspace(-d0, 0.0, n_pts_0, endpoint=False)
        z_start = 0.0
        z_pts_i = np.linspace(z_start, z_start+d[0], n_pts_i, endpoint=False)
        for i in range(1, len(d)):
            z_start += d[i-1]
            z_pts_i = np.concatenate((z_pts_i, np.linspace(z_start, z_start+d[i], n_pts_i, endpoint=False)))
        z_pts_1 = np.linspace(sum(d), sum(d)+d1,  n_pts_1, endpoint=True)
        
        
        n_0_vec = n_0*np.ones(n_pts_0)
        n_1_vec = n_1*np.ones(n_pts_1)
        
        n_i_vec = np.zeros(z_pts_i.shape, dtype=complex)

        for i in range(n_slabs):
            n_i = n_arr[i]
            n_i_vec[i*n_pts_i: (i+1)*n_pts_i] = n_i*np.ones(n_pts_i)

        z_pts = np.concatenate((z_pts_0, z_pts_i, z_pts_1))
        n_vec = np.concatenate((n_0_vec, n_i_vec, n_1_vec))
        
        return [z_pts, n_vec]
    
        
    def GetTransmissionPlusMinusKBand(self, k_0, k_1, N):
        k = np.linspace(k_0, k_1, N)
        Tp = np.zeros(N, dtype=complex)
        Tm = np.zeros(N, dtype=complex)
        for i in range(N):
            self.SetK(k[i])
            self.SetupTMMGen()
            Tp[i], Tm[i] = self.GetTransmissionPlusMinus()
            
        return [k, Tp, Tm]


    def StaircaseGradientBarrier1D(self, n0, y, d, s1, s2, n_slabs):
        # d: width of the barrier
        from Electromagnetics.GMTM import GMTM1D_T
        
        gmtm = GMTM1D_T()
        gmtm.SetProfileParameters(n0, y, d, s1, s2)
        
        def f_n(t):
            return gmtm.GetProfile(t)
        
        t_0 = 0.0
        t_1 = d
        return self.StaircaseContinuousProfile(f_n, n_slabs, t_0, t_1)





        
        
