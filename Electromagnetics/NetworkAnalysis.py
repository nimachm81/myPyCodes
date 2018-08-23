## NetworkAnalysis_sym.py
## Electromagnetic network analysis 
## scattering matrix analysis

from numpy import *
from cmath import *
from math import *

import scipy as sp
import numpy as np
import cmath
import math

__all__ = ["EmNetwork",
           "Tl_reflection", "TL_transmission", "TL_inpImp", "TL_inpImp_Lossy",
           "TL_Smat", "TL_junc_Smat", "TL_Tjunc_Smat"]


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
            if isinstance(s_mat, np.ndarray):
                s_size = s_mat.shape
                if len(s_size)==2 and s_size[0]==s_size[1]:
                    n_rows = s_size[0]
                    if self.__check_name_(name)[0]==True:
                        self.components.append([name, s_mat, [0]*n_rows])
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
        S = np.zeros((n_ports_tot, n_ports_tot), dtype=complex)
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
            S = np.delete(S, i, 0)
            S = np.delete(S, i, 1)
        return S


    def connectPorts(self, S, p_i, p_j):
        """
        It connects ports p_i and p_j of the given scattering matrix S
        and returns the new scattering matrix
        """
        n = S.shape[0]
        if S.shape[1] != n:
            raise ValueError('S Should be an square matrix.')
        row_pi = S[p_i,:].copy()
        for i in range(n):
            if row_pi[p_j]==1:
                print('warning: matrix ', row_pi, ' was 1 at index i=', p_j)
                print('to ensure power conservation all other components must be zero!')
                break
            temp = S[i,p_j]/(1.0-row_pi[p_j])
            #S[i*n] = S[i*n] + temp*S[p_i*n]
            for j in range(n):
                S[i, j] = S[i, j] + temp*row_pi[j]
        for i in range(n):
            S[i, p_j] = 0.0
            S[p_i, i] = 0.0
        row_pj = S[p_j,:].copy()
        for i in range(n):
            if row_pj[p_i]==1:
                print('warning: matrix ', row_pj, ' was 1 at index i=', p_i)
                print('to ensure power conservation all other components must be zero!')
                break
            temp = S[i,p_i]/(1.0-row_pj[p_i])
            #S[i*n] = S[i*n] + temp*S[p_j*n]
            for j in range(n):
                S[i, j] = S[i, j] + temp*row_pj[j]
        for i in range(n):
            S[i, p_i] = 0.0
            S[p_j, i] = 0.0
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



##------------  Transmisssion lines
""" Notice: the components inside the scattering matrices are related to power 
    waves (V/sqrt(Z)), and not the voltage waves (V)
    V: the voltage wave   Z:characteristic impedance of the port
"""


def Tl_reflection(Z_c, Z_L):
    """ Voltage reflection coefficient
    """
    if Z_L=='open':
        return 1.0
    return (Z_L-Z_c)/(Z_L+Z_c)
    
def TL_transmission(Z_c, Z_L):
    """ Voltage transmission coefficient T = 1 + R
    """
    R = Tl_reflection(Z_c, Z_L)
    return 1.0 + R
    
def TL_inpImp(theta, Z_c, Z_L='open'):
    """ input impedance of a lossless tranmission line of electrical length theta, 
    terminated by a load Z_L
    """
    if Z_L == 'open':
        tan_t = tan(theta)
        if tan_t==0.0:
            return 'open'
        elif tan_t==float('inf') or tan_t==-float('inf'):
            return 0.0
        else:
            return -1j*Z_c / tan(theta)
    return Z_c*(Z_L + 1j*Z_c*tan(theta))/(Z_c + 1j*Z_L*tan(theta))
    
    
def TL_inpImp_Lossy(gamma_l, Z_c, Z_L='open'):
    """ input impedance of a lossy tranmission line, 
    terminated by a load Z_L ---  gamma_l = gamma*l
    """
    if Z_L == 'open':
        return Z_c / tanh(gamma_l)
    return Z_c*(Z_L + Z_c*tanh(gamma_l))/(Z_c + Z_L*tanh(gamma_l))



def TL_Smat(theta, Z_c=None, Z_L=None):
    """ if Z_L==None it gives the scattering matrix of a 2 port transmission line
    of length theta. Otherwise it gives the scattering matrix of a 1 port loaded
    transmission line with characteristics impedance Z_c and load Z_L
    """
    if Z_c==None or Z_L==None:
        S = np.zeros((2,2), dtype=complex)
        S[1,0] = cmath.exp(-1j*theta)
        S[0,1] = cmath.exp(-1j*theta)
        return S
    else:
        Z_in = TL_inpImp(theta, Z_c, Z_L)
        R = Tl_reflection(Z_c, Z_in)
        S = np.zeros((1,1), dtype=complex)
        S[0,0] = R
        return S

def TL_junc_Smat(Z_c, Z_1):
    """ returns the 2*2 scattering matrix of a transmission line junstion where the 
    char impedance of the input line is Z_c and the char impedance of output line is Z_1 
    """
    r = Tl_reflection(Z_c, Z_1)
    #t = TL_transmission(Z_c, Z_1)
    t = cmath.sqrt(1.0 - abs(r)**2)
    
    S_junc = np.array([[r, t],
                       [t, -np.conjugate(r)]])
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
    
    S_junc = np.array([[x, y, y],
                       [y, z, w],
                       [y, w, z]])
 
    return S_junc
    
    





