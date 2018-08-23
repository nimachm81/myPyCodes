##  TMM_sym.py
##  Transmission Matrix Method - Symbolic


__all__ = ["TMMGen_sym", "TMM_sym_EM_N"]


from sympy import eye, zeros, solve, sin, Symbol, sqrt, Matrix, exp, I

import math, cmath



class TMMGen_sym:
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
        M_tot = eye(N)
        for i in range(len(layers)):
            ind_mat, params = layers[i]
            
            M_i_func = mats[ind_mat]
            M_i = M_i_func(*tuple(params))
            
            M_tot = M_i*M_tot
        
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
        A_tot = zeros(N, 2*N)
        A_tot[0:N, 0:N] = self.M_tot
        A_tot[0:N, N:2*N] = -eye(N)
        
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
        P = zeros(2*N, 2*N)
        for i in range(2*N):
            P[i, vars_order[i]] = 1
        
        ## set A_tot_new = A_tot*P  --->   A_tot_new*x_new = 0
        A_tot_new = A_tot*P
        
        ## set A_new and B_new such that --> A_new*x_new[0:N] = B_new*x_new[N:2N]
        A_new = A_tot_new[0:N, 0:N]
        B_new = -A_tot_new[0:N, N:2*N]

        ## set y=x_new[N:2N]
        y = zeros(N)
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
        b = B_new*y
        #x = np.linalg.inv(A_new).dot(b)
        x = A_new.inv()*b
        return x


    def GetLayerByLayerOutputs(self, in_vals):
        """ it takes the input values and returns the output after each layer 
        (after each multiplication)
        in_vals: input values (numpy vector)
        """
        layers = self.layers
        mats = self.mats
        N = self.n_var
        M_tot = eye(N)
        outs = [None]*len(layers)
        for i in range(len(layers)):
            ind_mat, params = layers[i]
            
            M_i_func = mats[ind_mat]
            M_i = M_i_func(*tuple(params))
            
            M_tot = M_i*M_tot
            outs[i] = M_tot*in_vals
            
        return outs






class TMM_sym_EM_N:
    """ TMM for electromagnetic wave, described by refractive index
    
    for theta!=0.0:
        for TE, T and R refer to electric field transmission and reflection
        for TM, T and R refer to magnetic field transmission and reflection
    for theta==0.0:
        T and R refer to electric field transmission and reflection
    """
    
    def __init__(self, TETM='TE'):
        """ theta: angle with respect to the interface normal in rad
        """
        self.multilayer = None
        self.theta = Symbol(r'\theta')
        self.TETM = TETM

        self.freq = Symbol('f')
        self.omega = Symbol(r'\omega')
        self.k_0 = Symbol('k_0')
        
        self.k_t = Symbol('k_t')
        self.k_z = Symbol('k_z')

        return
                
    def SetupMultilayer(self, n_0_kz0, n_1_kz1, n_arr_kzarr, d_arr):
        """ 
        """
        n_0, kz_n0 = n_0_kz0
        n_1, kz_n1 = n_0_kz0
        n_arr, kz_narr = n_arr_kzarr
        self.multilayer = [n_0, n_1, (n_arr, d_arr)]
        self.SetupMultilayer_Kz(kz_n0, kz_n1, kz_narr)
        self.SetupTMMGen()
        return

    def SetupMultilayer_Kz(self, k_n0, k_n1, k_narr):
        """ k_n0   -->  n_0
            k_n1   -->  n_1
            k_narr -->  n_arr
        """
        self.multilayer_kz = [k_n0, k_n1, k_narr]
        return

        
    def GetKt_sub(self):
        n_0 = self.multilayer[0]
        k_t_sub = self.k_0*n_0*sin(self.theta) 
        return k_t_sub
        
    def GetKz(self, n):
        n_0, n_1, n_arr__d = self.multilayer
        n_arr, d = n_arr__d
        n_slabs = len(n_arr)
        
        kz_n0, kz_n1, kz_narr = self.multilayer_kz
        if n==n_0:
            return kz_n0
        elif n==n_1:
            return kz_n1
        else:
            ind = n_arr.index(n)
            assert ind>=0
            return kz_narr[ind]          
        
    def GetKz_sub(self, n):
        k_0z_sq = self.k_0**2*n**2 - self.k_t**2
        k_0z = sqrt(k_0z_sq)
        return k_0z

        
    def MaxwellInterfaceMatrix(self, n_1, n_2):
        """ Electric field transmission matrix
            a:forward wave  b: backward wave
            a1 + b1 = a2 + b2
            (a1 - b1)/n1 = (a2 - b2)/n2
        """
        if self.TETM=='TE':
            k_z1 = self.GetKz(n_1)
            k_z2 = self.GetKz(n_2)
            TM = Matrix([[(k_z1 + k_z2)/(2*k_z2), (-k_z1 + k_z2)/(2*k_z2)], [(-k_z1 + k_z2)/(2*k_z2), (k_z1 + k_z2)/(2*k_z2)]])
            return TM
        else:
            assert self.TETM=='TM'
            k_z1 = self.GetKz(n_1)
            k_z2 = self.GetKz(n_2)
            TM = Matrix([[(k_z1*n_2**2 + k_z2*n_1**2)/(2*k_z2*n_1**2), (-k_z1*n_2**2 + k_z2*n_1**2)/(2*k_z2*n_1**2)], \
                    [(-k_z1*n_2**2 + k_z2*n_1**2)/(2*k_z2*n_1**2), (k_z1*n_2**2 + k_z2*n_1**2)/(2*k_z2*n_1**2)]])
            return TM


    def MaxwellPropagMatrix(self, d, k_z):
        """ propagation matrix n:refractive index   d: thickness   k_z:propagation number (normal to interface)
            a2 = a1*exp(-j*k*d)
            b2 = b1*exp(+j*k*d)
        """
        PM = Matrix([[exp(-I*d*k_z), 0], [0, exp(I*d*k_z)]])
        return PM
        
               
    
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
        tmmgen = TMMGen_sym(mats, params, n_var=2)
        
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
        io_vals = [['in', 0, 1], ['out', 1, 0]]
        io_vars = [['in', 1], ['out', 0]]
        RT = self.tmmgen.setInputsOutputs(io_vals, io_vars)

        return RT    
        

