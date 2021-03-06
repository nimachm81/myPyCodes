## NetworkAnalysis_sym.py
## Electromagnetic network analysis 
## scattering matrix analysis


from sympy import *

__all__ = ["EmNetwork", "TL",
           "Tl_reflection", "TL_junc_Smat", "TL_Tjunc_Smat"]


class EmNetwork:
    """
    Defining an EM network as a collection of components
    described by scattering matrices, whose ports can become
    interconnected using the class methods
    """

    def __init__(self, smat_list=[], name_list=[]):
        self.components = []
        self.connections = []
        if len(smat_list)==len(name_list):
            for i in range(len(smat_list)):
                self.addComponent(smat_list[i], name_list[i])
        else:
            raise ValueError("S-matrix list and name_list should have the same length.")
        return



    def addComponent(self, s_mat, name, n=1):
        """
        adds the component to the network
        if n=1: only one component with the given name is added
        for n>1: n components with the names name_0, name_1, ... are added
        """
        if n==1:
            if isinstance(s_mat, Matrix):
                if s_mat.rows==s_mat.cols:
                    if self.__check_name_(name)[0]==True:
                        self.components.append([name, s_mat, [0]*s_mat.rows])
                    else:
                        raise ValueError(self.__check_name_(name)[1])
                else:
                    raise ShapeError("scattering matrics should be sqaure matrics")
            else:
                raise TypeError("argument should be a list of scattering matrices")
        else:
            for i in range(n):
                self.addComponent(s_mat, name+'_'+str(i), 1)
        return



    def addConnection(self, comp1_name, port_1, comp2_name, port_2):
        """
        Connects port 1 of component 1 to port 2 of component 2
        """
        chk_1 = self.__check_ports_(comp1_name, port_1)
        chk_2 = self.__check_ports_(comp2_name, port_2)
        if chk_1[0]==True:
            if chk_2[0]==True:
                self.connections.append([[comp1_name, port_1], [comp2_name, port_2]])
                self.closePort(comp1_name, port_1)
                self.closePort(comp2_name, port_2)
            else:
                raise ValueError(chk_2[1])
        else:
            raise ValueError(chk_1[1])
        return


    def getTotalScatteringMatrix_ignoreConnections(self):
        """
        Returns a big scattering matrix representing all the ports in
        a single matrix - all the ports are assumed open. port connections
        are then applied in the next step
        """
        n_comps = len(self.components)
        n_ports_i = [len(self.components[i][2]) for i in range(n_comps)]
        n_ports_tot = sum(n_ports_i)
        S = zeros(n_ports_tot)
        start = 0
        for i in range(n_comps):
            for j in range(n_ports_i[i]):
                for k in range(n_ports_i[i]):
                    S[start+j, start+k] = self.components[i][1][j, k]
            start += n_ports_i[i]
        return S

    
    def getTotalScatteringMatrix(self):
        """
        It creates the total scattering matrix, taking into account all the 
        port connections
        """
        S = self.getTotalScatteringMatrix_ignoreConnections()
        n_connects = len(self.connections)
        to_delete = []
        for i in range(n_connects):
            p_0_tot = self.__port_to_port_total(self.connections[i][0][0], self.connections[i][0][1])
            p_1_tot = self.__port_to_port_total(self.connections[i][1][0], self.connections[i][1][1])
            S = self.connectPorts(S, p_0_tot, p_1_tot)
            to_delete.append(p_0_tot)
            to_delete.append(p_1_tot)
        to_delete.sort(reverse=True)
        for i in to_delete:
            S.row_del(i)
            S.col_del(i)
        return S


    def connectPorts(self, S, p_i, p_j):
        """
        It connects ports p_i and p_j of the given scattering matrix S
        and returns the new scattering matrix
        """
        n = S.rows
        if S.cols != n:
            raise ValueError('S Should be an square matrix.')
        row_pi = S.row(p_i)
        for i in range(n):
            if row_pi[p_j]==1:
                print('warning: matrix ', row_pi, ' was 1 at index i=', p_j)
                print('to ensure power conservation all other components must be zero!')
                break
            temp = S[i,p_j]/(1-row_pi[p_j])
            #S[i*n] = S[i*n] + temp*S[p_i*n]
            for j in range(n):
                S[i, j] = S[i, j] + temp*row_pi[j]
        for i in range(n):
            S[i, p_j] = 0
            S[p_i, i] = 0
        row_pj = S.row(p_j)
        for i in range(n):
            if row_pj[p_i]==1:
                print('warning: matrix ', row_pj, ' was 1 at index i=', p_i)
                print('to ensure power conservation all other components must be zero!')
                break
            temp = S[i,p_i]/(1-row_pj[p_i])
            #S[i*n] = S[i*n] + temp*S[p_j*n]
            for j in range(n):
                S[i, j] = S[i, j] + temp*row_pj[j]
        for i in range(n):
            S[i, p_i] = 0
            S[p_j, i] = 0
        return S


    def closePort(self, comp_name, port):
        ## marks port as connected   0:open   1:connected
        chk = self.__check_ports_(comp_name, port)
        if chk[0]==True:
            self.components[chk[1]][2][port] = 1
        return

    def __port_to_port_total(self, comp_name, port):
        """
        Gets component name and port number and returns port index among 
        all the ports
        """
        n_comps = len(self.components)
        port_ind_start = 0
        for i in range(0, n_comps):
            if comp_name==self.components[i][0]:
                port_tot = port + port_ind_start
                return port_tot
            port_ind_start += len(self.components[i][2])
        raise ValueError('port not found')
        return


    def __check_name_(self, name):
        for comp in self.components:
            if comp[0]==name:
                return [False, "component name already exists!"]
        return [True]


    def __check_ports_(self, comp_name, port):
        for i, comp in enumerate(self.components):
            if comp[0]==comp_name:
                if port<len(comp[2]):
                    if comp[2][port]==0:
                        return [True, i]
                    else:
                        return [False, 'port closed']
                else:
                    return [False, 'port does not exist']
        return [False, 'component does not exist']





class TL:
    """
    Transmission Line Parameters
    """
    ##TODO: for index=n add n as subscript to parameters R_n, C_n...
    def __init__(self, lossless=False, index=None):
        self.omega=Symbol('\\omega', real=True)
        self.R=Symbol('R', real=True)
        self.L=Symbol('L', real=True)
        self.G=Symbol('G', real=True)
        self.C=Symbol('C', real=True)
        self.Z_c = Symbol('Z_c')
        self.gamma = Symbol('\\gamma')
        self.alpha = Symbol('\\alpha', real=True)
        self.beta = Symbol('\\beta', real=True)
        self.Z_L = Symbol('Z_L')
        self.l = Symbol('l')
        if lossless:
            self.R = 0
            self.G = 0
            self.Z_c = Symbol('Z_c', real=True)
        return
        
    def isLossless(self, lossless=True):
        if lossless:
            self.R = 0
            self.G = 0
            self.Z_c = Symbol('Z_c', real=True)
        else:
            self.R=Symbol('R', real=True)
            self.G=Symbol('G', real=True)
            self.Z_c = Symbol('Z_c')
        return    
        
    def propConst(self):
        return sqrt((self.R+I*self.omega*self.L)*(self.G+I*self.omega*self.C))
        
    def charImp(self):
        return sqrt((self.R+I*self.omega*self.L)/(self.G+I*self.omega*self.C))
        
    def reflection(self):
        if self.Z_L==oo or self.Z_L==-oo:
            return 1
        return (self.Z_L - self.Z_c)/(self.Z_L + self.Z_c)
        
    def transmission(self):
        return 2*self.Z_L/(self.Z_L + self.Z_c)        
    
    def inputImp(self):
        if self.R==0 and self.G==0:
            if self.Z_L!=oo or self.Z_L==-oo:
                return self.Z_c*(self.Z_L + I*self.Z_c*tan(self.beta*self.l))/\
                    (self.Z_c + I*self.Z_L*tan(self.beta*self.l))
            else:
                return -I*self.Z_c/tan(self.beta*self.l)
        else:
            if self.Z_L!=oo or self.Z_L==-oo:
                return self.Z_c*(self.Z_L + self.Z_c*tanh(self.gamma*self.l))/\
                    (self.Z_c + self.Z_L*tanh(self.gamma*self.l))
            else:
                return self.Z_c/tanh(self.gamma*self.l)
        
        
def Tl_reflection(Z_c, Z_L):
    """ Voltage reflection coefficient
    """
    if Z_L=='open':
        return 1.0
    return (Z_L-Z_c)/(Z_L+Z_c)

def TL_junc_Smat(Z_c, Z_1):
    """ returns the 2*2 scattering matrix of a transmission line junstion where the 
    char impedance of the input line is Z_c and the char impedance of output line is Z_1 
    """
    r = Tl_reflection(Z_c, Z_1)
    #t = TL_transmission(Z_c, Z_1)
    t = cmath.sqrt(1.0 - abs(r)**2)
    
    S_junc = Matrix([[r, t],
                     [t, -conjugate(r)]])
    return S_junc


def TL_Tjunc_Smat(Z_c, Z_1):
    """ returns the 3*3 scattering matrix of a symmetric T junstion where the 
    char impedance of the main line is Z_c and the char impedance of branches 
    are Z_1
    """
    x = Tl_reflection(Z_c, Z_1/2.0)
    y = (1.0/sqrt(2.0)*sqrt(1.0-abs(x)**2))
    Zc_Z1_par = Z_c*Z_1/(Z_c + Z_1)
    z = (Zc_Z1_par - Z_1)/(Zc_Z1_par + Z_1)
    w = -x - z
    
    S_junc = Matrix([[x, y, y],
                     [y, z, w],
                     [y, w, z]])
 
    return S_junc


