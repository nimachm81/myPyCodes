## PhotonicCristals.py


__all__ = ["PC2D_EPS"]

from Electromagnetics.VectorCalculus import *
from Electromagnetics.FourierBlochND import *
from IPython.display import display, Math, Latex
from scipy import constants
import math

from Electromagnetics.Misc import Point2D, Point3D, GetReiprocalVecs

from sympy import Symbol, symbols, Matrix, I, latex, simplify, Mul,\
                 Add, Sum, Pow, Integer, Function, oo
import sympy

from scipy.sparse.linalg import LinearOperator, eigs, eigsh
import numpy as np

import pickle
import os


A_CONVS_IND = 1


class PC2D_EPS:
    """ 2D photonic crystal lossless TM using H components
        the wave propagates in the xy plane
    """
    def __init__(self, mode='TE', Ns=[16, 16], vb=False, symm=False, mu_r=1.0, 
            DX=[1.0, 1.0], D_phi=math.pi/2.0):
        """ D_phi the angle between the 2 axis
            direct cell: -DX[0]:DX[0] along the first axis and
                         -DX[1]:DX[1] along the second axis
        """
        assert mode=='TM' or mode=='TE'
        
        self.mode = mode
        
        self.x, self.y, self.z = symbols('x y z')
        self.Hz = Symbol('H_{z}')
        self.Ez = Symbol('E_{z}')
        self.k, self.k0 = symbols('k k_0')
        self.omega = Symbol('\\omega')
        self.eps_r, self.mu_r = symbols('\\epsilon_r \\mu_r')
        self.eps_0, self.mu_0 = symbols('\\epsilon_0 \\mu_0')
        self.eps_r_inv = Symbol('\\epsilon_r^{-1}')
        
        self.Ns = Ns
        
        self.mu_r_val = mu_r        
        self.vbose = vb
        
        self.f_eps_inv = None
        
        # -DX[i]/2.0:DX[i]/2.0 domain size along axix i, generally assumed oriented (prallelogram)
        self.DX = DX
        self.D_phi = D_phi
        
        v_dir_0 = Point3D(self.DX[0], 0.0, 0.0)
        v_dir_1 = Point3D(self.DX[1]*math.cos(D_phi), self.DX[1]*math.sin(D_phi), 0.0)
        self.vecDirec_vals = [v_dir_0, v_dir_1]
        print('v_dir_0: ', v_dir_0)
        print('v_dir_1: ', v_dir_1)

        v_rec_0, v_rec_1 = GetReiprocalVecs([v_dir_0, v_dir_1])
        self.vecRecip_vals = [v_rec_0, v_rec_1]
        print('v_rec_0: ', v_rec_0)
        print('v_rec_1: ', v_rec_1)
        
        
        self.b0x, self.b0y, self.b0z = symbols('b_{0x} b_{0y} b_{0z}')
        self.b1x, self.b1y, self.b1z = symbols('b_{1x} b_{1y} b_{1z}')
        self.vecRecip = [Matrix([[self.b0x, self.b0y, 0]]), Matrix([[self.b1x, self.b1y, 0]])]
        
        self.theta = Symbol('\\theta')
        
        self.useFFT = True
        self.n_mat_vec_fft = 0
        self.n_mat_vec_dir = 0
        
        
        
    def set_f_eps_rect(self, eps_0, eps_1, dx1, dy1):
        ## dx1: width    dy1:height
        assert eps_0.imag==0 and eps_1.imag==0
        #assert dx1<self.DX[0] and dy1<self.DX[1]
        def f(r):
            return 1.0/self.f_eps_rect(r, eps_0, eps_1, dx1, dy1)
        self.f_eps_inv = f
        
    def set_f_eps_rect_Fermi(self, eps_0, eps_1, dx1, dy1, beta=10):
        ## dx1: width    dy1:height
        assert eps_0.imag==0 and eps_1.imag==0
        #assert dx1<self.DX[0] and dy1<self.DX[1]
        def f(r):
            return 1.0/self.f_eps_rect_Fermi(r, eps_0, eps_1, dx1, dy1, beta=beta)
        self.f_eps_inv = f

    def set_f_eps_circ(self, eps_0, eps_1, dr):
        assert eps_0.imag==0 and eps_1.imag==0
        #assert dr<self.DX[0]/2.0 and dr<self.DX[1]/2.0
        def f(r):
            return 1.0/self.f_eps_circ(r, eps_0, eps_1, dr)
        self.f_eps_inv = f

    def set_f_eps_circ_Fermi(self, eps_0, eps_1, dr, beta=10.0):
        assert eps_0.imag==0 and eps_1.imag==0
        #assert dr<self.DX[0]/2.0 and dr<self.DX[1]/2.0
        def f(r):
            return 1.0/self.f_eps_circ_Fermi(r, eps_0, eps_1, dr, beta=beta)
        self.f_eps_inv = f

    def set_f_eps_circ_honeycomb_supercell_Fermi(self, eps_0, eps_1, dr, beta=10.0, SC_scales=[3, 3]):
        assert eps_0.imag==0 and eps_1.imag==0
        assert SC_scales[0]==SC_scales[1] and SC_scales[0]>1
        #assert dr<self.DX[0]/2.0 and dr<self.DX[1]/2.0
        def f(r):
            return 1.0/self.f_eps_circ_honeycomb_supercell_Fermi(r, eps_0, eps_1, dr, beta=beta, SC_scales=SC_scales)
        self.f_eps_inv = f

    def set_f_eps_circ_honeycomb_defect_supercell_Fermi(self, eps_0, eps_1, eps_2, dr_1, dr_2, beta=10.0, SC_scales=[3, 3]):
        assert eps_0.imag==0 and eps_1.imag==0
        assert SC_scales[0]==SC_scales[1] and SC_scales[0]>1
        #assert dr<self.DX[0]/2.0 and dr<self.DX[1]/2.0
        def f(r):
            return 1.0/self.f_eps_circ_honeycomb_defect_supercell_Fermi(r, eps_0, eps_1, eps_2, dr_1, dr_2, beta=beta, SC_scales=SC_scales)
        self.f_eps_inv = f

    def set_f_eps_circ_defect_supercell_Fermi(self, eps_0, eps_1, eps_2, dr_1, dr_2, beta, SC_scales, mat_defect):
        assert eps_0.imag==0 and eps_1.imag==0
        assert SC_scales[0]==SC_scales[1] and SC_scales[0]>1
        #assert dr<self.DX[0]/2.0 and dr<self.DX[1]/2.0
        def f(r):
            return 1.0/self.f_eps_circ_defect_supercell_Fermi(r, eps_0, eps_1,
             eps_2, dr_1, dr_2, beta, SC_scales, mat_defect)
        self.f_eps_inv = f

    def set_f_eps_circ_gaussian(self, eps_0, eps_1, dr):
        assert eps_0.imag==0 and eps_1.imag==0
        #assert dr<self.DX[0]/2.0 and dr<self.DX[1]/2.0
        def f(r):
            return 1.0/self.f_eps_circ_gaussian(r, eps_0, eps_1, dr)
        self.f_eps_inv = f

    def f_eps_rect(self, r, eps_0, eps_1, dx1, dy1):
        ## 0:background index     1:rod index   dx1, dy1: rod width and height
        r_x = r[0] + r[1]*math.cos(self.D_phi)
        r_y = r[1]*math.sin(self.D_phi)
        eps = ((abs(r_x)<=dx1/2.0)*(abs(r_y)<=dy1/2.0))*eps_1 + \
              ((abs(r_x)>dx1/2.0) + (abs(r_y)>dy1/2.0))*eps_0
        return eps              
        
    def f_eps_rect_Fermi(self, r, eps_0, eps_1, dx1, dy1, beta):
        ## 0:background index     1:rod index   dx1, dy1: rod width and height
        r_x = r[0] + r[1]*math.cos(self.D_phi)
        r_y = r[1]*math.sin(self.D_phi)
        eps = 1.0/(1.0+np.exp(beta*(np.abs(r_x)-dx1/2.0)))/(1.0+np.exp(beta*(np.abs(r_y)-dy1/2.0)))*(eps_1-eps_0) + eps_0
        return eps              


    def f_eps_circ(self, r, eps_0, eps_1, dr):
        ## 0:background index     1:rod index   dr: rod radius
        r_x = r[0] + r[1]*math.cos(self.D_phi)
        r_y = r[1]*math.sin(self.D_phi)
        eps = ((r_x**2 + r_y**2)<=dr**2)*eps_1 + \
              ((r_x**2 + r_y**2)>dr**2)*eps_0
        return eps              
        
    def f_eps_circ_Fermi(self, r, eps_0, eps_1, dr, beta):
        ## 0:background index     1:rod index   dr: rod radius
        #assert eps_0<eps_1
        r_x = r[0] + r[1]*math.cos(self.D_phi)
        r_y = r[1]*math.sin(self.D_phi)
        r_x, r_y = self.moveRecCoordInsideBZ([r_x, r_y])
        eps = 1.0/(1.0+np.exp(beta*(np.sqrt((r_x**2 + r_y**2))-dr)))*(eps_1-eps_0) + eps_0
        return eps              

    def f_eps_circ_honeycomb_supercell_Fermi(self, r, eps_0, eps_1, dr, beta, SC_scales):
        ## 0:background index     1:rod index   dr: rod radius
        #assert eps_0<eps_1
        r_x = r[0] + r[1]*math.cos(self.D_phi)
        r_y = r[1]*math.sin(self.D_phi)
        r_xy, an_bn = self.moveRecCoordInsideBZGetDistanceToFBZ([r_x, r_y], SC_scales)
        r_x, r_y = r_xy
        a_n, b_n = an_bn
        assert SC_scales[0]==SC_scales[1] and SC_scales[0]%3==0
        scale = SC_scales[0]
        # honeycomb lattice points
        HC_atoms = ((a_n%scale==0)*(b_n%scale==0) + (a_n%scale==2*int(scale/3))*(b_n%scale==0) +
                    (a_n%scale==0)*(b_n%scale==1*int(scale/3)) + (a_n%scale==1*int(scale/3))*(b_n%scale==1*int(scale/3)) +
                    (a_n%scale==1*int(scale/3))*(b_n%scale==2*int(scale/3)) + (a_n%scale==2*int(scale/3))*(b_n%scale==2*int(scale/3))).astype(bool)
        eps_hc_atoms = (1.0/(1.0+np.exp(beta*(np.sqrt((r_x**2 + r_y**2))-dr)))*(eps_1-eps_0) + eps_0)*HC_atoms
        #background esilon
        eps_back = eps_0*(np.logical_not(HC_atoms))
        eps = eps_hc_atoms + eps_back
        return eps              

    def f_eps_circ_honeycomb_defect_supercell_Fermi(self, r, eps_0, eps_1, eps_2, dr_1, dr_2, beta, SC_scales):
        ## 0:background index     1:rod index   2: defect index   dr_1: rod radius dr_2: defect radius
        #assert eps_0<eps_1
        r_x = r[0] + r[1]*math.cos(self.D_phi)
        r_y = r[1]*math.sin(self.D_phi)
        r_xy, an_bn = self.moveRecCoordInsideBZGetDistanceToFBZ([r_x, r_y], SC_scales)
        r_x, r_y = r_xy
        a_n, b_n = an_bn
        assert SC_scales[0]==SC_scales[1] and SC_scales[0]%3==0
        scale = SC_scales[0]
        # honeycomb lattice points
        HC_atoms = ((a_n%scale==0)*(b_n%scale==0) + (a_n%scale==2*int(scale/3))*(b_n%scale==0) +
                    (a_n%scale==0)*(b_n%scale==1*int(scale/3)) + (a_n%scale==1*int(scale/3))*(b_n%scale==1*int(scale/3)) +
                    (a_n%scale==1*int(scale/3))*(b_n%scale==2*int(scale/3)) + (a_n%scale==2*int(scale/3))*(b_n%scale==2*int(scale/3))).astype(bool)
        eps_hc_atoms = (1.0/(1.0+np.exp(beta*(np.sqrt((r_x**2 + r_y**2))-dr_2)))*(eps_2-eps_0) + eps_0)*HC_atoms
        #background esilon
        eps_back = (1.0/(1.0+np.exp(beta*(np.sqrt((r_x**2 + r_y**2))-dr_1)))*(eps_1-eps_0) + eps_0)*(np.logical_not(HC_atoms))
        eps = eps_hc_atoms + eps_back
        return eps              
        
    def f_eps_circ_defect_supercell_Fermi(self, r, eps_0, eps_1, eps_2, dr_1, dr_2, beta, SC_scales, mat_defect):
        ## 0:background index     1:rod index   2: defect index   dr_1: rod radius dr_2: defect radius
        # mat_defect: SC_scales[0]*SC_scales[1] matrix of 0s and 1s, where 1s denote defects
        assert mat_defect.shape == tuple(SC_scales)
        r_x = r[0] + r[1]*math.cos(self.D_phi)
        r_y = r[1]*math.sin(self.D_phi)
        r_xy, an_bn = self.moveRecCoordInsideBZGetDistanceToFBZ([r_x, r_y], SC_scales)
        r_x, r_y = r_xy
        a_n, b_n = an_bn
        scale_0, scale_1 = SC_scales
        # defect points
        HC_atoms = np.zeros(a_n.shape, dtype=int)
        for i in range(scale_0):
            for j in range(scale_1):
                if mat_defect[i, j] == 1:
                    HC_atoms += (a_n%scale_0==i)*(b_n%scale_1==j)
        HC_atoms = HC_atoms.astype(bool)
        eps_hc_atoms = (1.0/(1.0+np.exp(beta*(np.sqrt((r_x**2 + r_y**2))-dr_2)))*(eps_2-eps_0) + eps_0)*HC_atoms
        #background esilon
        eps_back = (1.0/(1.0+np.exp(beta*(np.sqrt((r_x**2 + r_y**2))-dr_1)))*(eps_1-eps_0) + eps_0)*(np.logical_not(HC_atoms))
        eps = eps_hc_atoms + eps_back
        return eps              
        

    def f_eps_circ_gaussian(self, r, eps_0, eps_1, dr):
        ## 0:background index     1:rod index   dr: rod radius
        #assert eps_0<eps_1
        r_x = r[0] + r[1]*math.cos(self.D_phi)
        r_y = r[1]*math.sin(self.D_phi)
        r_x, r_y = self.moveRecCoordInsideBZ([r_x, r_y])
        eps = np.exp(-(r_x**2 + r_y**2)/dr**2)*(eps_1-eps_0) + eps_0
        return eps              

    def get_eps_r_inv_harmonics(self, X):
        # Domain:[-X[0]/2:X[0]/2, -X[1]/2:X[1]/2] the axis are in genreral oriented
        eps_r_inv_vec = self.pc2d.getFourierCoeffs(self.f_eps_inv, [-X[0]/2.0, -X[1]/2.0], 
                [X[0]/2.0, X[1]/2.0], self.Ns)
        return eps_r_inv_vec
        
    def moveRecCoordInsideBZ(self, r):
        # BZ for the direct cell (the smallest volume closest to the center 
        # with all the information)
        # r is in rectangular coordinates
        r = self.moveToFirstUnitCell(r)
        v0 = Point2D(self.vecDirec_vals[0].x, self.vecDirec_vals[0].y)
        v1 = Point2D(self.vecDirec_vals[1].x, self.vecDirec_vals[1].y)
        v2, v3 = v0+v1, v0-v1
        vecs_unit = [v0, v1, v2, v3, -v0, -v1, -v2, -v3]
        r_ = r.copy()
        r_x_list = [r[0]-vecs_unit[i].x for i in range(len(vecs_unit))]
        r_y_list = [r[1]-vecs_unit[i].y for i in range(len(vecs_unit))]
        r_abs_list = [r_x_list[i]**2+r_y_list[i]**2 for i in range(len(r_x_list))]
        r_abs = r_[0]**2 + r_[1]**2
        
        for i in range(len(r_abs_list)):
            r_[0] = (r_abs<=r_abs_list[i])*r_[0] + (r_abs>r_abs_list[i])*r_x_list[i]
            r_[1] = (r_abs<=r_abs_list[i])*r_[1] + (r_abs>r_abs_list[i])*r_y_list[i]
            r_abs = r_[0]**2 + r_[1]**2
        
        return r_
        
    def moveRecCoordInsideBZGetDistanceToFBZ(self, r, SC_scales):
        """ BZ for the direct cell (the smallest volume closest to the center 
            with all the information)
            r is in rectangular coordinates
            gets the distance to first BZ of the amall unit cells inside the supercell
            self.vecDirec_vals represents the super cell unit vectors
            SC_scales: super cell scales i.e: self.vecDirec_vals[i]/SC_scales[i] 
            are the smallest units inside the supercell
        """
        r, an_bn = self.moveToFirstUnitCellGetDistanceToFUC(r, SC_scales)
        a_n, b_n = an_bn
        v0 = Point2D(self.vecDirec_vals[0].x, self.vecDirec_vals[0].y)/float(SC_scales[0])
        v1 = Point2D(self.vecDirec_vals[1].x, self.vecDirec_vals[1].y)/float(SC_scales[1])
        v2, v3 = v0+v1, v0-v1
        vecs_unit = [v0, v1, v2, v3, -v0, -v1, -v2, -v3]
        BZ_distance_0 = np.array([1, 0, 1, 1, -1, 0, -1, -1])
        BZ_distance_1 = np.array([0, 1, 1, -1, 0, -1, -1, 1])
        BZ_index = -np.ones(r[0].shape, dtype=int)
        vec_ones = np.ones(r[0].shape, dtype=int)
        r_ = r.copy()
        r_x_list = [r[0]-vecs_unit[i].x for i in range(len(vecs_unit))]
        r_y_list = [r[1]-vecs_unit[i].y for i in range(len(vecs_unit))]
        r_abs_list = [r_x_list[i]**2+r_y_list[i]**2 for i in range(len(r_x_list))]
        r_abs = r_[0]**2 + r_[1]**2
        
        for i in range(len(r_abs_list)):
            r_[0] = (r_abs<=r_abs_list[i])*r_[0] + (r_abs>r_abs_list[i])*r_x_list[i]
            r_[1] = (r_abs<=r_abs_list[i])*r_[1] + (r_abs>r_abs_list[i])*r_y_list[i]
            BZ_index = (r_abs<=r_abs_list[i])*BZ_index + (r_abs>r_abs_list[i])*i
            r_abs = r_[0]**2 + r_[1]**2
            
        #print(a_n.shape, BZ_index.shape, BZ_distance.shape)
        
        a_n = a_n + BZ_distance_0[(BZ_index>=0)*BZ_index]*(BZ_index>=0)
        b_n = b_n + BZ_distance_1[(BZ_index>=0)*BZ_index]*(BZ_index>=0)
        
        return [r_, [a_n, b_n]]
        

    def rectangularToDirectVectorCoords(self, r):
        # returns (a, b) where r = a*v0 + b*v1 where v0 and v1 are the unit cell vectors
        # x = a + b*cos(theta)
        # y = b*sin(theta)
        # --->  b = y/sin(theta)  ,  a = x - y*cot(theta)        
        a = (r[0] - r[1]/np.tan(self.D_phi))/self.DX[0]
        b = r[1]/np.sin(self.D_phi)/self.DX[1]
        return [a, b]
        
    def moveToFirstUnitCell(self, r):
        a, b = self.rectangularToDirectVectorCoords(r)
        a = a - np.floor(a)
        b = b - np.floor(b)
        a = a - (a>=0.5)*1.0
        b = b - (b>=0.5)*1.0
        x = a*self.DX[0] + b*self.DX[1]*np.cos(self.D_phi)
        y = b*self.DX[1]*np.sin(self.D_phi)
        return [x, y]

    def moveToFirstUnitCellGetDistanceToFUC(self, r, SC_scales):
        # how much a cell is displaced from the FUC first unit cell
        # SC_scales: supercell scales
        a, b = self.rectangularToDirectVectorCoords(r)
        a = a*SC_scales[0]
        b = b*SC_scales[1]
        a_n = np.floor(a)
        b_n = np.floor(b)
        a = a - a_n
        b = b - b_n
        a_n = a_n + (a>=0.5)*1.0
        b_n = b_n + (b>=0.5)*1.0
        a = a - (a>=0.5)*1.0
        b = b - (b>=0.5)*1.0
        x = a*self.DX[0]/SC_scales[0] + b*self.DX[1]/SC_scales[1]*np.cos(self.D_phi)
        y = b*self.DX[1]/SC_scales[1]*np.sin(self.D_phi)
        return [[x, y], [a_n, b_n]]

    def setupEquations(self, regenerateEQs=True):
        # theta: propagation angle (with respect to axis 0)
        Bloch_fact = sympy.exp(-I*(self.k*sympy.cos(self.theta)*self.x + self.k*sympy.sin(self.theta)*self.y))

        EQ = None
        if self.mode=='TE':
            H = Matrix([[0, 0, self.Hz]])*Bloch_fact
            E = 1/(I*self.omega*self.eps_0)*self.eps_r_inv*curl_r(H)
            
            #EQ = curl_r(E) + I*self.omega*self.mu_0*self.mu_r*H
            #EQ = self.omega*self.eps_0*curl_r(E) + I*self.omega*self.mu_0*self.mu_r*H
            EQ = I*self.omega/(self.mu_0*self.mu_r)*curl_r(E) - self.omega**2*H

            self.Gz = self.Hz           ## Gz = Ez or Hz depending on the mode
        else:
            E = Matrix([[0, 0, self.Ez]])*Bloch_fact
            H = -1/(I*self.omega*self.mu_0*self.mu_r)*curl_r(E)
            
            #EQ = curl_r(H) - I*self.omega*self.eps_0*self.eps_r*E
            EQ = -I*self.omega*self.eps_r_inv/self.eps_0*curl_r(H) - self.omega**2*E
        
            self.Gz = self.Ez
        
        if self.vbose:
            display(Math(latex(EQ.T)))
        
        
        pc2d = PDEFourierSeriesND(EQ, vars_fourier=[self.Gz, self.eps_r_inv], n_dim=2, vecRecip=self.vecRecip)
        if regenerateEQs==True:
            EQ_F = pc2d.putSums()
            EQ_F = Matrix([[EQ_F[i,j].doit() for j in range(EQ_F.cols)] for i in range(EQ_F.rows)])

            if self.vbose:
                display(Math(latex(EQ_F.T)))
                display(Math(latex(pc2d.varsHarm)))

            EQ_F = [EQ_F[2]]

            EQ_F = [pc2d.applyConvolutions(EQ_F[i]) for i in range(len(EQ_F))]
            if self.vbose:
                display(Math(latex(Matrix(EQ_F))))

            EQ_F = [simplify(pc2d.applyOrthogonalities(EQ_F[i])/Bloch_fact) for i in range(len(EQ_F))]
            if self.vbose:
                display(Math(latex(Matrix(EQ_F))))
            
            self.pc2d = pc2d
            self.EQ_F = EQ_F

            from sympy import srepr
            for i in range(len(self.EQ_F)):
                print('equation ', i, ': ', srepr(self.EQ_F[i]))
            
            
        elif self.mode=='TM':
            pc2d.createVarsHarm()
            self.pc2d = pc2d
            from sympy import sin, cos, exp, Tuple

            EQ_F_0 = Mul(Pow(Symbol('\\epsilon_0'), Integer(-1)), Pow(Symbol('\\mu_0'), Integer(-1)), Pow(Symbol('\\mu_r'), Integer(-1)), Add(Mul(Integer(-1), Symbol('\\epsilon_0'), Symbol('\\mu_0'), Symbol('\\mu_r'), Pow(Symbol('\\omega'), Integer(2)), Function('\\tilde{E_{z}}')(Symbol('n_0'), Symbol('n_1'))), Mul(Pow(Symbol('k'), Integer(2)), Pow(sin(Symbol('\\theta')), Integer(2)), Sum(Mul(Function('\\tilde{E_{z}}')(Symbol('m_0'), Symbol('m_1')), Function('\\tilde{\\epsilon_r^{-1}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1')))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo))), Mul(Pow(Symbol('k'), Integer(2)), Pow(cos(Symbol('\\theta')), Integer(2)), Sum(Mul(Function('\\tilde{E_{z}}')(Symbol('m_0'), Symbol('m_1')), Function('\\tilde{\\epsilon_r^{-1}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1')))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo))), Mul(Symbol('k'), sin(Symbol('\\theta')), Sum(Mul(Integer(2), Symbol('b_{0y}'), Symbol('m_0'), Function('\\tilde{E_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo))), Mul(Symbol('k'), sin(Symbol('\\theta')), Sum(Mul(Integer(-1), Integer(2), Symbol('b_{0y}'), Symbol('n_0'), Function('\\tilde{E_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo))), Mul(Symbol('k'), sin(Symbol('\\theta')), Sum(Mul(Integer(2), Symbol('b_{1y}'), Symbol('m_1'), Function('\\tilde{E_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo))), Mul(Symbol('k'), sin(Symbol('\\theta')), Sum(Mul(Integer(-1), Integer(2), Symbol('b_{1y}'), Symbol('n_1'), Function('\\tilde{E_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo))), Mul(Symbol('k'), cos(Symbol('\\theta')), Sum(Mul(Integer(2), Symbol('b_{0x}'), Symbol('m_0'), Function('\\tilde{E_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo))), Mul(Symbol('k'), cos(Symbol('\\theta')), Sum(Mul(Integer(-1), Integer(2), Symbol('b_{0x}'), Symbol('n_0'), Function('\\tilde{E_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo))), Mul(Symbol('k'), cos(Symbol('\\theta')), Sum(Mul(Integer(2), Symbol('b_{1x}'), Symbol('m_1'), Function('\\tilde{E_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo))), Mul(Symbol('k'), cos(Symbol('\\theta')), Sum(Mul(Integer(-1), Integer(2), Symbol('b_{1x}'), Symbol('n_1'), Function('\\tilde{E_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo))), Sum(Mul(Pow(Symbol('b_{0x}'), Integer(2)), Pow(Symbol('m_0'), Integer(2)), Function('\\tilde{E_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Pow(Symbol('b_{0x}'), Integer(2)), Pow(Symbol('n_0'), Integer(2)), Function('\\tilde{E_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Pow(Symbol('b_{0y}'), Integer(2)), Pow(Symbol('m_0'), Integer(2)), Function('\\tilde{E_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Pow(Symbol('b_{0y}'), Integer(2)), Pow(Symbol('n_0'), Integer(2)), Function('\\tilde{E_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Pow(Symbol('b_{1x}'), Integer(2)), Pow(Symbol('m_1'), Integer(2)), Function('\\tilde{E_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Pow(Symbol('b_{1x}'), Integer(2)), Pow(Symbol('n_1'), Integer(2)), Function('\\tilde{E_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Pow(Symbol('b_{1y}'), Integer(2)), Pow(Symbol('m_1'), Integer(2)), Function('\\tilde{E_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Pow(Symbol('b_{1y}'), Integer(2)), Pow(Symbol('n_1'), Integer(2)), Function('\\tilde{E_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Integer(-1), Integer(2), Pow(Symbol('b_{0x}'), Integer(2)), Symbol('m_0'), Symbol('n_0'), Function('\\tilde{E_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Integer(-1), Integer(2), Pow(Symbol('b_{0y}'), Integer(2)), Symbol('m_0'), Symbol('n_0'), Function('\\tilde{E_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Integer(-1), Integer(2), Pow(Symbol('b_{1x}'), Integer(2)), Symbol('m_1'), Symbol('n_1'), Function('\\tilde{E_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Integer(-1), Integer(2), Pow(Symbol('b_{1y}'), Integer(2)), Symbol('m_1'), Symbol('n_1'), Function('\\tilde{E_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Integer(2), Symbol('b_{0x}'), Symbol('b_{1x}'), Symbol('m_0'), Symbol('m_1'), Function('\\tilde{E_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Integer(-1), Integer(2), Symbol('b_{0x}'), Symbol('b_{1x}'), Symbol('m_0'), Symbol('n_1'), Function('\\tilde{E_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Integer(-1), Integer(2), Symbol('b_{0x}'), Symbol('b_{1x}'), Symbol('m_1'), Symbol('n_0'), Function('\\tilde{E_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Integer(2), Symbol('b_{0x}'), Symbol('b_{1x}'), Symbol('n_0'), Symbol('n_1'), Function('\\tilde{E_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Integer(2), Symbol('b_{0y}'), Symbol('b_{1y}'), Symbol('m_0'), Symbol('m_1'), Function('\\tilde{E_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Integer(-1), Integer(2), Symbol('b_{0y}'), Symbol('b_{1y}'), Symbol('m_0'), Symbol('n_1'), Function('\\tilde{E_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Integer(-1), Integer(2), Symbol('b_{0y}'), Symbol('b_{1y}'), Symbol('m_1'), Symbol('n_0'), Function('\\tilde{E_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Integer(2), Symbol('b_{0y}'), Symbol('b_{1y}'), Symbol('n_0'), Symbol('n_1'), Function('\\tilde{E_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo))))
        
            scale_fact = self.mu_0*self.eps_0
            self.EQ_F = [EQ_F_0*scale_fact]
            
            if self.vbose:
                display(Math(latex(Matrix(self.EQ_F))))
        else:
            assert self.mode=='TE'
            pc2d.createVarsHarm()
            self.pc2d = pc2d
            from sympy import sin, cos, exp, Tuple

            EQ_F_0 = Mul(Pow(Symbol('\\epsilon_0'), Integer(-1)), Pow(Symbol('\\mu_0'), Integer(-1)), Pow(Symbol('\\mu_r'), Integer(-1)), Add(Mul(Integer(-1), Symbol('\\epsilon_0'), Symbol('\\mu_0'), Symbol('\\mu_r'), Pow(Symbol('\\omega'), Integer(2)), Function('\\tilde{H_{z}}')(Symbol('n_0'), Symbol('n_1'))), Mul(Pow(Symbol('k'), Integer(2)), Pow(sin(Symbol('\\theta')), Integer(2)), Sum(Mul(Function('\\tilde{H_{z}}')(Symbol('m_0'), Symbol('m_1')), Function('\\tilde{\\epsilon_r^{-1}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1')))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo))), Mul(Pow(Symbol('k'), Integer(2)), Pow(cos(Symbol('\\theta')), Integer(2)), Sum(Mul(Function('\\tilde{H_{z}}')(Symbol('m_0'), Symbol('m_1')), Function('\\tilde{\\epsilon_r^{-1}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1')))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo))), Mul(Symbol('k'), sin(Symbol('\\theta')), Sum(Mul(Symbol('b_{0y}'), Symbol('m_0'), Function('\\tilde{H_{z}}')(Symbol('m_0'), Symbol('m_1')), Function('\\tilde{\\epsilon_r^{-1}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1')))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo))), Mul(Symbol('k'), sin(Symbol('\\theta')), Sum(Mul(Integer(2), Symbol('b_{0y}'), Symbol('m_0'), Function('\\tilde{H_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo))), Mul(Symbol('k'), sin(Symbol('\\theta')), Sum(Mul(Integer(-1), Symbol('b_{0y}'), Symbol('n_0'), Function('\\tilde{H_{z}}')(Symbol('m_0'), Symbol('m_1')), Function('\\tilde{\\epsilon_r^{-1}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1')))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo))), Mul(Symbol('k'), sin(Symbol('\\theta')), Sum(Mul(Integer(-1), Integer(2), Symbol('b_{0y}'), Symbol('n_0'), Function('\\tilde{H_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo))), Mul(Symbol('k'), sin(Symbol('\\theta')), Sum(Mul(Symbol('b_{1y}'), Symbol('m_1'), Function('\\tilde{H_{z}}')(Symbol('m_0'), Symbol('m_1')), Function('\\tilde{\\epsilon_r^{-1}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1')))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo))), Mul(Symbol('k'), sin(Symbol('\\theta')), Sum(Mul(Integer(2), Symbol('b_{1y}'), Symbol('m_1'), Function('\\tilde{H_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo))), Mul(Symbol('k'), sin(Symbol('\\theta')), Sum(Mul(Integer(-1), Symbol('b_{1y}'), Symbol('n_1'), Function('\\tilde{H_{z}}')(Symbol('m_0'), Symbol('m_1')), Function('\\tilde{\\epsilon_r^{-1}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1')))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo))), Mul(Symbol('k'), sin(Symbol('\\theta')), Sum(Mul(Integer(-1), Integer(2), Symbol('b_{1y}'), Symbol('n_1'), Function('\\tilde{H_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo))), Mul(Symbol('k'), cos(Symbol('\\theta')), Sum(Mul(Symbol('b_{0x}'), Symbol('m_0'), Function('\\tilde{H_{z}}')(Symbol('m_0'), Symbol('m_1')), Function('\\tilde{\\epsilon_r^{-1}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1')))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo))), Mul(Symbol('k'), cos(Symbol('\\theta')), Sum(Mul(Integer(2), Symbol('b_{0x}'), Symbol('m_0'), Function('\\tilde{H_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo))), Mul(Symbol('k'), cos(Symbol('\\theta')), Sum(Mul(Integer(-1), Symbol('b_{0x}'), Symbol('n_0'), Function('\\tilde{H_{z}}')(Symbol('m_0'), Symbol('m_1')), Function('\\tilde{\\epsilon_r^{-1}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1')))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo))), Mul(Symbol('k'), cos(Symbol('\\theta')), Sum(Mul(Integer(-1), Integer(2), Symbol('b_{0x}'), Symbol('n_0'), Function('\\tilde{H_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo))), Mul(Symbol('k'), cos(Symbol('\\theta')), Sum(Mul(Symbol('b_{1x}'), Symbol('m_1'), Function('\\tilde{H_{z}}')(Symbol('m_0'), Symbol('m_1')), Function('\\tilde{\\epsilon_r^{-1}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1')))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo))), Mul(Symbol('k'), cos(Symbol('\\theta')), Sum(Mul(Integer(2), Symbol('b_{1x}'), Symbol('m_1'), Function('\\tilde{H_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo))), Mul(Symbol('k'), cos(Symbol('\\theta')), Sum(Mul(Integer(-1), Symbol('b_{1x}'), Symbol('n_1'), Function('\\tilde{H_{z}}')(Symbol('m_0'), Symbol('m_1')), Function('\\tilde{\\epsilon_r^{-1}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1')))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo))), Mul(Symbol('k'), cos(Symbol('\\theta')), Sum(Mul(Integer(-1), Integer(2), Symbol('b_{1x}'), Symbol('n_1'), Function('\\tilde{H_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo))), Sum(Mul(Integer(-1), Pow(Symbol('b_{0x}'), Integer(2)), Pow(Symbol('m_0'), Integer(2)), Function('\\tilde{H_{z}}')(Symbol('m_0'), Symbol('m_1')), Function('\\tilde{\\epsilon_r^{-1}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1')))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Pow(Symbol('b_{0x}'), Integer(2)), Pow(Symbol('m_0'), Integer(2)), Function('\\tilde{H_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Pow(Symbol('b_{0x}'), Integer(2)), Pow(Symbol('n_0'), Integer(2)), Function('\\tilde{H_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Integer(-1), Pow(Symbol('b_{0y}'), Integer(2)), Pow(Symbol('m_0'), Integer(2)), Function('\\tilde{H_{z}}')(Symbol('m_0'), Symbol('m_1')), Function('\\tilde{\\epsilon_r^{-1}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1')))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Pow(Symbol('b_{0y}'), Integer(2)), Pow(Symbol('m_0'), Integer(2)), Function('\\tilde{H_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Pow(Symbol('b_{0y}'), Integer(2)), Pow(Symbol('n_0'), Integer(2)), Function('\\tilde{H_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Integer(-1), Pow(Symbol('b_{1x}'), Integer(2)), Pow(Symbol('m_1'), Integer(2)), Function('\\tilde{H_{z}}')(Symbol('m_0'), Symbol('m_1')), Function('\\tilde{\\epsilon_r^{-1}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1')))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Pow(Symbol('b_{1x}'), Integer(2)), Pow(Symbol('m_1'), Integer(2)), Function('\\tilde{H_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Pow(Symbol('b_{1x}'), Integer(2)), Pow(Symbol('n_1'), Integer(2)), Function('\\tilde{H_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Integer(-1), Pow(Symbol('b_{1y}'), Integer(2)), Pow(Symbol('m_1'), Integer(2)), Function('\\tilde{H_{z}}')(Symbol('m_0'), Symbol('m_1')), Function('\\tilde{\\epsilon_r^{-1}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1')))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Pow(Symbol('b_{1y}'), Integer(2)), Pow(Symbol('m_1'), Integer(2)), Function('\\tilde{H_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Pow(Symbol('b_{1y}'), Integer(2)), Pow(Symbol('n_1'), Integer(2)), Function('\\tilde{H_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Pow(Symbol('b_{0x}'), Integer(2)), Symbol('m_0'), Symbol('n_0'), Function('\\tilde{H_{z}}')(Symbol('m_0'), Symbol('m_1')), Function('\\tilde{\\epsilon_r^{-1}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1')))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Integer(-1), Integer(2), Pow(Symbol('b_{0x}'), Integer(2)), Symbol('m_0'), Symbol('n_0'), Function('\\tilde{H_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Pow(Symbol('b_{0y}'), Integer(2)), Symbol('m_0'), Symbol('n_0'), Function('\\tilde{H_{z}}')(Symbol('m_0'), Symbol('m_1')), Function('\\tilde{\\epsilon_r^{-1}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1')))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Integer(-1), Integer(2), Pow(Symbol('b_{0y}'), Integer(2)), Symbol('m_0'), Symbol('n_0'), Function('\\tilde{H_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Pow(Symbol('b_{1x}'), Integer(2)), Symbol('m_1'), Symbol('n_1'), Function('\\tilde{H_{z}}')(Symbol('m_0'), Symbol('m_1')), Function('\\tilde{\\epsilon_r^{-1}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1')))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Integer(-1), Integer(2), Pow(Symbol('b_{1x}'), Integer(2)), Symbol('m_1'), Symbol('n_1'), Function('\\tilde{H_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Pow(Symbol('b_{1y}'), Integer(2)), Symbol('m_1'), Symbol('n_1'), Function('\\tilde{H_{z}}')(Symbol('m_0'), Symbol('m_1')), Function('\\tilde{\\epsilon_r^{-1}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1')))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Integer(-1), Integer(2), Pow(Symbol('b_{1y}'), Integer(2)), Symbol('m_1'), Symbol('n_1'), Function('\\tilde{H_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Integer(-1), Symbol('b_{0x}'), Symbol('b_{1x}'), Symbol('m_0'), Symbol('m_1'), Function('\\tilde{H_{z}}')(Symbol('m_0'), Symbol('m_1')), Function('\\tilde{\\epsilon_r^{-1}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1')))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Integer(-1), Symbol('b_{0x}'), Symbol('b_{1x}'), Symbol('m_0'), Symbol('m_1'), Function('\\tilde{H_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Integer(2), Symbol('b_{0x}'), Symbol('b_{1x}'), Symbol('m_0'), Symbol('m_1'), Function('\\tilde{H_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Symbol('b_{0x}'), Symbol('b_{1x}'), Symbol('m_0'), Symbol('n_1'), Function('\\tilde{H_{z}}')(Symbol('m_0'), Symbol('m_1')), Function('\\tilde{\\epsilon_r^{-1}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1')))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Integer(-1), Integer(2), Symbol('b_{0x}'), Symbol('b_{1x}'), Symbol('m_0'), Symbol('n_1'), Function('\\tilde{H_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Symbol('b_{0x}'), Symbol('b_{1x}'), Symbol('m_0'), Symbol('n_1'), Function('\\tilde{H_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Integer(-1), Integer(2), Symbol('b_{0x}'), Symbol('b_{1x}'), Symbol('m_1'), Symbol('n_0'), Function('\\tilde{H_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Integer(2), Symbol('b_{0x}'), Symbol('b_{1x}'), Symbol('n_0'), Symbol('n_1'), Function('\\tilde{H_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Integer(-1), Symbol('b_{0y}'), Symbol('b_{1y}'), Symbol('m_0'), Symbol('m_1'), Function('\\tilde{H_{z}}')(Symbol('m_0'), Symbol('m_1')), Function('\\tilde{\\epsilon_r^{-1}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1')))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Integer(-1), Symbol('b_{0y}'), Symbol('b_{1y}'), Symbol('m_0'), Symbol('m_1'), Function('\\tilde{H_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Integer(2), Symbol('b_{0y}'), Symbol('b_{1y}'), Symbol('m_0'), Symbol('m_1'), Function('\\tilde{H_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Symbol('b_{0y}'), Symbol('b_{1y}'), Symbol('m_0'), Symbol('n_1'), Function('\\tilde{H_{z}}')(Symbol('m_0'), Symbol('m_1')), Function('\\tilde{\\epsilon_r^{-1}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1')))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Integer(-1), Integer(2), Symbol('b_{0y}'), Symbol('b_{1y}'), Symbol('m_0'), Symbol('n_1'), Function('\\tilde{H_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Symbol('b_{0y}'), Symbol('b_{1y}'), Symbol('m_0'), Symbol('n_1'), Function('\\tilde{H_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Integer(-1), Integer(2), Symbol('b_{0y}'), Symbol('b_{1y}'), Symbol('m_1'), Symbol('n_0'), Function('\\tilde{H_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo)), Sum(Mul(Integer(2), Symbol('b_{0y}'), Symbol('b_{1y}'), Symbol('n_0'), Symbol('n_1'), Function('\\tilde{H_{z}}')(Add(Mul(Integer(-1), Symbol('m_0')), Symbol('n_0')), Add(Mul(Integer(-1), Symbol('m_1')), Symbol('n_1'))), Function('\\tilde{\\epsilon_r^{-1}}')(Symbol('m_0'), Symbol('m_1'))), Tuple(Symbol('m_0'), -oo, oo), Tuple(Symbol('m_1'), -oo, oo))))
            
            scale_fact = self.mu_0*self.eps_0
            self.EQ_F = [EQ_F_0*scale_fact]
            
            if self.vbose:
                display(Math(latex(Matrix(self.EQ_F))))
           
        
    def initializeEqsk(self):
        self.Gz_tilde, self.eps_r_inv_tilde = self.pc2d.varsHarm

        self.eps_r_inv_vec = self.get_eps_r_inv_harmonics(self.DX)#TODO: double check --> *math.sin(self.D_phi)
        print(self.eps_r_inv_vec)

        EQ_F_fin_k = [self.EQ_F[i].subs([(self.eps_0, constants.epsilon_0),
             (self.mu_0, constants.mu_0), (self.mu_r, self.mu_r_val), 
             (self.b0x, self.vecRecip_vals[0][0]), (self.b0y, self.vecRecip_vals[0][1]),
             (self.b1x, self.vecRecip_vals[1][0]), (self.b1y, self.vecRecip_vals[1][1])]) for i in range(len(self.EQ_F))]
             
        self.EQ_F_fin_k = EQ_F_fin_k
        
        if self.vbose:
            display(Math(latex(Matrix(EQ_F_fin_k))))
            
    def set_eps_r_vec_for_test(self):
        self.eps_r_inv_vec = np.random.random((2*self.Ns[0], 2*self.Ns[1]))
        for i in range(2*self.Ns[0]):
            for j in range(2*self.Ns[1]):
                if i<=self.Ns[0]/2 or j<=self.Ns[1]/2 or i>=3*self.Ns[0]/2-1 or j>=3*self.Ns[1]/2-1:
                    self.eps_r_inv_vec[i,j] = 0

    def get_randon_vec_for_test(self):
        x_vec = np.random.random((2*self.Ns[0], 2*self.Ns[1])) + 1j*np.random.random((2*self.Ns[0], 2*self.Ns[1]))
        for i in range(2*self.Ns[0]):
            for j in range(2*self.Ns[1]):
                if i<=self.Ns[0]/2 or j<=self.Ns[1]/2 or i>=3*self.Ns[0]/2-1 or j>=3*self.Ns[1]/2-1:
                    x_vec[i,j] = 0
        x_vec = x_vec.reshape((self.pc2d.n_total,))        
        return x_vec


    def setEqskinit(self, k_init, theta_init=0.0):
        EQ_F_fin = [self.EQ_F_fin_k[i].subs([(self.k, k_init), (self.theta, theta_init)]) for i in range(len(self.EQ_F_fin_k))]
        self.EQ_F_fin = EQ_F_fin
        if self.vbose:
            display(Math(latex(Matrix(EQ_F_fin))))

        self.pc2d.setupNumericalParameters(EQ_F_fin, Ns=self.Ns, vars=[self.Gz_tilde], 
            pars=[self.eps_r_inv_tilde], pars_vecs=[self.eps_r_inv_vec], eig_vars=[self.omega])

        self.pc2d.save_n_as_powers = False
        self.pc2d.save_m_as_powers = False
        self.pc2d.useFFTforConvs = self.useFFT
        self.A_mat = self.pc2d.orthogonalToNumpyMatrix()
    

    def matVecProductFFT(self, x_vec):
        [y_tot, y_tot_eigs] = self.pc2d.getConvProductsFFT(x_vec)
        assert y_tot_eigs==[]
        self.n_mat_vec_fft += 1
        return y_tot
        
    def solve_A__sig_I_FFT(self, b_rhs, eig_0):
        """ calculates (A - sigma*I)*x = b_rhs using iterative solvers
        """
        def matVecProductFFT_shifted(x_vec):
            return self.matVecProductFFT(x_vec) - eig_0*x_vec

        n_tot = self.pc2d.n_total
        A_shape = (n_tot, n_tot)
        A = LinearOperator(A_shape, matvec=matVecProductFFT_shifted, dtype=complex)        
        
        from scipy.sparse.linalg import gmres, cg
        
        x = cg(A, b_rhs, tol=1.0e-2, maxiter=None)
        #print(x)
        print('converged..', b_rhs.shape)
        return x[0]
    

    def matVecProductDirect(self, x_vec):
        A = self.A_mat[A_CONVS_IND]
        y_tot = A.dot(x_vec)
        self.n_mat_vec_dir += 1
        return y_tot

    def solveEigenvalues(self, n_eigs=5, eig_0=None, eig_vec_0=None, 
                which='LM', tol=1.0e-5, ncv=None, maxiter=None, shift_manually=True):
        """ shift_manually: manually shifts matrix to A - eig_0*I (to avoid negative 
            eigenvalues if necessary (optional) ) and calculates the 
            smallest magnitude eigenvalues and shifts back the calculated eigenvalues
        """
        def matVecProductFFT_shifted(x_vec):
            return self.matVecProductFFT(x_vec) - eig_0*x_vec
                
        def matVecProductDirect_shifted(x_vec):
            return self.matVecProductDirect(x_vec) - eig_0*x_vec

        def OPinv_A_sig_I_FFT(b):
            return self.solve_A__sig_I_FFT(b, eig_0)

        A = None
        A__sig_I_inv = None
        self.n_mat_vec_fft = 0
        self.n_mat_vec_dir = 0
        if self.useFFT:
            n_tot = self.pc2d.n_total
            A_shape = (n_tot, n_tot)
            if shift_manually==True and eig_0!=None:
                A = LinearOperator(A_shape, matvec=matVecProductFFT_shifted, dtype=complex)
            elif eig_0!=None:
                print('eig_0!=None ... A__sig_I_inv')
                A__sig_I_inv = LinearOperator(A_shape, matvec=OPinv_A_sig_I_FFT, dtype=complex)
                A = LinearOperator(A_shape, matvec=self.matVecProductFFT, dtype=complex)
            else:
                A = LinearOperator(A_shape, matvec=self.matVecProductFFT, dtype=complex)
        else:
            n_tot = self.pc2d.n_total
            A_shape = (n_tot, n_tot)
            if shift_manually==True and eig_0!=None:
                A = LinearOperator(A_shape, matvec=matVecProductDirect_shifted, dtype=complex)        
            else:
                A = LinearOperator(A_shape, matvec=self.matVecProductDirect, dtype=complex)        
            #A = self.A_mat[A_CONVS_IND]
    

        if shift_manually==True and eig_0!=None:
            omega_eigs = eigs(A, k=n_eigs, sigma=None, which='SM', v0=eig_vec_0, tol=tol, 
                    ncv=ncv, maxiter=maxiter)
            for i in range(len(omega_eigs[0])):
                omega_eigs[0][i] += eig_0
            return omega_eigs
        else:
            omega_eigs = eigs(A, k=n_eigs, sigma=eig_0, which=which, v0=eig_vec_0, tol=tol, 
                    ncv=ncv, maxiter=maxiter, OPinv=A__sig_I_inv)
            return omega_eigs
            

    def solveEigenvalues_RQ_NLCG(self, eig_0=None, eig_vec_0=None, eigs_prev=None, eigs_vec_prev=None):
        """ find the eigenvalues using Rayleigh Quotient minimization with 
            non-linear conjugate gradient method
        """
        n_total = self.pc2d.n_total
        
        n_prev = 0
        if eigs_prev!=None:
            n_prev = len(eigs_prev)
            assert len(eigs_vec_prev)==n_prev
            assert eig_0!=None
            assert eig_0.imag==0.0 and eig_0.real>=0.0
            for i in range(n_prev):
                assert eig_0 >= eigs_prev[i]
                assert eigs_prev[i].imag==0.0 and eigs_prev[i].real>=0.0
                if i<n_prev-1:
                    assert eigs_prev[i]<=eigs_prev[i+1]
                    
        c__ = 4.0
        
        def RQ(x_vec_2):
            x_vec = x_vec_2[0:n_total] + 1j*x_vec_2[n_total:2*n_total]
            Ax = self.matVecProductFFT(x_vec)
            
            xAx = np.conjugate(x_vec).dot(Ax)
            xx = np.conjugate(x_vec).dot(x_vec)
            rq = (xAx/xx).real
            print(rq)
            
            for i in range(n_prev):
                x_i = eigs_vec_prev[i]
                xx_i = np.conjugate(x_vec).dot(eigs_vec_prev[i])
                xixi = np.conjugate(x_i).dot(x_i)
                rq += c__*(eig_0 - eigs_prev[i])*abs(xx_i)**2/(xx*xixi)
            return rq.real
            
        def RQ_prime(x_vec_2):
            x_vec = x_vec_2[0:n_total] + 1j*x_vec_2[n_total:2*n_total]
            Ax = self.matVecProductFFT(x_vec)
            
            xAx = np.conjugate(x_vec).dot(Ax)
            xx = np.conjugate(x_vec).dot(x_vec)
            rq = (xAx/xx).real
            
            rq_prime = 2.0*(Ax - rq*x_vec)/xx
            
            for i in range(n_prev):
                x_i = eigs_vec_prev[i]
                xx_i = np.conjugate(x_vec).dot(eigs_vec_prev[i])
                xixi = np.conjugate(x_i).dot(x_i)
                N_prime = 2.0*np.conjugate(xx_i)*x_i 
                D_prime = 2.0*xixi*x_vec
                N_ = abs(xx_i)**2
                D_ = (xx*xixi)
                rq_prime += c__*(eig_0 - eigs_prev[i])*(N_prime-N_/D_*D_prime)/D_
            
            rq_prime_2 = np.zeros(2*n_total)
            rq_prime_2[0:n_total] = np.real(rq_prime)
            rq_prime_2[n_total:2*n_total] = np.imag(rq_prime)
            return rq_prime_2

        from scipy.optimize import fmin_cg
        
        x0 = None
        if eig_vec_0==None:
            x0 = np.random.rand(2*n_total)
        else:
            x0 = np.zeros(2*n_total)
            x0[0:n_total] = np.real(eig_vec_0)
            x0[n_total:2*n_total] = np.imag(eig_vec_0)
        res = fmin_cg(RQ, x0, fprime=RQ_prime, epsilon=1.0e-5, gtol=1.0e-3, maxiter=1000)
        #print(res)
        return [RQ(res), res]
        
        


    def getVarsFields(self, x_vec):
        """ take the IFFT and get the fields in the direct space
            here Ez or Hz depending on the mode (TE/TM)
        """
        vars_fields = self.pc2d.getVarsFields(x_vec)
        assert len(vars_fields)==1
        return vars_fields[0]#/self.pc2d.getInverseFourierCoeffs(self.eps_r_inv_vec)



    def sweepK(self, k_0_2D, k_1_2D, n_pts=30, n_bands=3, eigs_vecs_init=None, 
                    saveToFile=None, readFromFile=None, endpoint=False, params={"traceBands":False}):
        ##TODO: there might be a problem with identifying igenvalues using orthogonalities
        # if they are degenerate
        eig_shift = -10.0
        
        k_2D_pts = [None]*n_pts
        dk_2D = k_1_2D - k_0_2D
        for i in range(n_pts):
            if endpoint==True:
                k_2D_pts[i] = k_0_2D + dk_2D*float(i)/(n_pts-1)
            else:
                k_2D_pts[i] = k_0_2D + dk_2D*float(i)/n_pts
        
        omega_eigs = [None]*n_bands
        omega_eigs_vecs = [None]*n_bands
        for i in range(n_bands):
            omega_eigs[i] = [None]*n_pts
            omega_eigs_vecs[i] = [None]*n_pts
        
        def saveDataToFile(i_k):
            if saveToFile!=None:
                f_name = saveToFile
                f_dir = os.path.dirname(f_name)
                if not os.path.exists(f_dir):
                    os.makedirs(f_dir)
                f = open(f_name, 'wb')
                params_dic = {'k_0_2D': k_0_2D, 'k_1_2D':k_1_2D, 'n_pts':n_pts, 
                    'n_bands':n_bands, 'eigs_vecs_init': eigs_vecs_init,
                    'eig_shift':eig_shift, 'i_k_start':i_k, 'endpoint':endpoint, 
                    'params':params}
                pickle.dump(params_dic, f)
                pickle.dump(k_2D_pts, f)
                pickle.dump(omega_eigs, f)
                pickle.dump(omega_eigs_vecs, f)
                f.close()

        i_k_start = 0
        if readFromFile!=None:
            f_name = readFromFile
            f = open(f_name, 'rb')
            params_dic = pickle.load(f)
            k_2D_pts = pickle.load(f)
            omega_eigs = pickle.load(f)
            omega_eigs_vecs = pickle.load(f)
            k_0_2D = params_dic['k_0_2D']
            k_1_2D = params_dic['k_1_2D']
            n_pts = params_dic['n_pts']
            n_bands = params_dic['n_bands']
            eigs_vecs_init = params_dic['eigs_vecs_init']
            eig_shift = params_dic['eig_shift']
            i_k_start = params_dic['i_k_start']
            endpoint = params_dic['endpoint']
            params = params_dic['params']
            f.close()
                

        n_solve = n_bands+2
        if params["traceBands"]==False:
            n_solve = n_bands

        for i_k in range(i_k_start, n_pts):
            k_i_2D = k_2D_pts[i_k]
            k_i, theta_i = self.getMagAngleK2D(k_i_2D)
            
            print('i:', i_k,  'k_i_2D: ', k_i_2D, 'k_i:', k_i, 'theta_i:', theta_i)
            self.setEqskinit(k_i, theta_i)
            
            omega_eigs_and_vecs = self.solveEigenvalues(n_eigs=n_solve, eig_0=eig_shift, eig_vec_0=None, which='SM', tol=1.0e-3, 
                                   ncv=3*n_solve+1, maxiter=None, shift_manually=True)
            
            omega_eigs_i = omega_eigs_and_vecs[0]
            omega_eigs_vecs_i = omega_eigs_and_vecs[1].T
            
            print('omega_eigs_i: ', np.sqrt(omega_eigs_i)/(2.0*math.pi))
            print('number of mat vec products: ', self.n_mat_vec_fft)
            assert n_solve>=n_bands
            
            inds_omega_sorted = None
            omega_eigs_i_abs = [abs(omega_eigs_i[j]) for j in range(n_solve)]
            inds_omega_sorted = sorted(range(n_solve), key=lambda j: omega_eigs_i_abs[j])
            print('inds_omega_sorted: ', inds_omega_sorted)

            if params["traceBands"]==True:
                for i_b in range(n_bands):
                    vec_ib_prev = None
                    if i_k>0:
                        vec_ib_prev = omega_eigs_vecs[i_b][i_k-1]
                    elif eigs_vecs_init!=None:
                        vec_ib_prev = eigs_vecs_init[i_b]
                    else:
                        omega_eigs[i_b][i_k] = omega_eigs_i[inds_omega_sorted[i_b]]
                        omega_eigs_vecs[i_b][i_k] = omega_eigs_vecs_i[inds_omega_sorted[i_b]]
                        continue
                        
                    corr = [0.0]*n_solve
                    vec_ib_prev /= np.linalg.norm(vec_ib_prev)
                    for i_s in range(n_solve):
                        vec_is = omega_eigs_vecs_i[i_s]
                        vec_is /= np.linalg.norm(vec_is)
                        corr[i_s] = np.abs(np.conjugate(vec_ib_prev).dot(vec_is))
                    
                    print(i_b, ' correlations: \n', corr)
                    ind_corr = 0
                    corr_max = corr[0]
                    for i_s in range(1, n_solve):
                        if corr[i_s] > corr_max:
                            ind_corr = i_s
                            corr_max = corr[i_s]
                    
                    if corr_max>0.9:
                        omega_eigs[i_b][i_k] = omega_eigs_i[ind_corr]
                        omega_eigs_vecs[i_b][i_k] = omega_eigs_vecs_i[ind_corr]
                    elif i_k>0:
                        omega_eigs[i_b][i_k] = None
                        omega_eigs_vecs[i_b][i_k] = omega_eigs_vecs[i_b][i_k-1]
                    else:
                        assert False
                        
                    
                    sq_omeg_2pi = omega_eigs[i_b].copy()
                    for j__ in range(len(sq_omeg_2pi)):
                        if sq_omeg_2pi[j__]!=None:
                            #print(sq_omeg_2pi[j__])
                            sq_omeg_2pi[j__] = np.sqrt(sq_omeg_2pi[j__])/(2.0*math.pi)
                    print(sq_omeg_2pi, '\n')
            else:
                for i_b in range(n_bands):
                    omega_eigs[i_b][i_k] = omega_eigs_i[inds_omega_sorted[i_b]]
                    omega_eigs_vecs[i_b][i_k] = omega_eigs_vecs_i[inds_omega_sorted[i_b]]

                    sq_omeg_2pi = omega_eigs[i_b].copy()
                    for j__ in range(len(sq_omeg_2pi)):
                        if sq_omeg_2pi[j__]!=None:
                            #print(sq_omeg_2pi[j__])
                            sq_omeg_2pi[j__] = np.sqrt(sq_omeg_2pi[j__])/(2.0*math.pi)
                    print(sq_omeg_2pi, '\n')

            if i_k>0 and i_k%5==0:
                saveDataToFile(i_k+1)
        saveDataToFile(n_pts)
        return [k_2D_pts, omega_eigs, omega_eigs_vecs]

    def sweepK_RQ(self, k_0_2D, k_1_2D, n_pts=30, n_bands=3, eigs_init=None, eigs_vecs_init=None, 
                    saveToFile=None, readFromFile=None, endpoint=False, params={"traceBands":False}):
        ##TODO: degenerate eigenvalues might pollute the interpolation process (to be checked)
        eig_shift = -10.0
        
        k_2D_pts = [None]*n_pts
        dk_2D = k_1_2D - k_0_2D
        for i in range(n_pts):
            if endpoint==True:
                k_2D_pts[i] = k_0_2D + dk_2D*float(i)/(n_pts-1)
            else:
                k_2D_pts[i] = k_0_2D + dk_2D*float(i)/n_pts
        
        omega_eigs = [None]*n_bands
        omega_eigs_vecs = [None]*n_bands
        for i in range(n_bands):
            omega_eigs[i] = [None]*n_pts
            omega_eigs_vecs[i] = [None]*n_pts
                    
        def saveDataToFile(i_k):
            if saveToFile!=None:
                f_name = saveToFile
                f_dir = os.path.dirname(f_name)
                if not os.path.exists(f_dir):
                    os.makedirs(f_dir)
                f = open(f_name, 'wb')
                params_dic = {'k_0_2D': k_0_2D, 'k_1_2D':k_1_2D, 'n_pts':n_pts, 
                    'n_bands':n_bands, 'eigs_init': eigs_init, 'eigs_vecs_init': eigs_vecs_init,
                    'eig_shift':eig_shift, 'i_k_start':i_k, 'endpoint':endpoint, 
                    'params':params}
                pickle.dump(params_dic, f)
                pickle.dump(k_2D_pts, f)
                pickle.dump(omega_eigs, f)
                pickle.dump(omega_eigs_vecs, f)
                f.close()

        i_k_start = 0
        if readFromFile!=None:
            f_name = readFromFile
            f = open(f_name, 'rb')
            params_dic = pickle.load(f)
            k_2D_pts = pickle.load(f)
            omega_eigs = pickle.load(f)
            omega_eigs_vecs = pickle.load(f)
            k_0_2D = params_dic['k_0_2D']
            k_1_2D = params_dic['k_1_2D']
            n_pts = params_dic['n_pts']
            n_bands = params_dic['n_bands']
            eigs_init = params_dic['eigs_init']
            eigs_vecs_init = params_dic['eigs_vecs_init']
            eig_shift = params_dic['eig_shift']
            i_k_start = params_dic['i_k_start']
            endpoint = params_dic['endpoint']
            params = params_dic['params']
            f.close()
                

        def extrap_eigs_eigsvecs(i_b, i_k, N_interp=3):
            #i_b: band index   i_k: k index 
            coeff_lagrange = None
            N_lagr = None
            k_samp = None
            eig_samp = None
            eig_vec_samp = None
            k = None
            if i_k==0:
                inds_init_sorted = None
                eigs_init_abs = [abs(eigs_init[j]) for j in range(len(eigs_init))]
                ind_eigs_init_sorted = sorted(range(len(eigs_init)), key=lambda j: eigs_init_abs[j])
                return [eigs_init[ind_eigs_init_sorted[i_b]], eigs_vecs_init[ind_eigs_init_sorted[i_b]]]
            elif i_k<=N_interp:
                k_samp = [(k_2D_pts[i]-k_0_2D).norm() for i in range(i_k)]
                eig_samp = [omega_eigs[i_b][i] for i in range(i_k)]
                eig_vec_samp = [omega_eigs_vecs[i_b][i] for i in range(i_k)]
                k = (k_2D_pts[i_k]-k_0_2D).norm()
            else:
                k_samp = [(k_2D_pts[i]-k_0_2D).norm() for i in range(i_k-N_interp, i_k)]
                eig_samp = [omega_eigs[i_b][i] for i in range(i_k-N_interp, i_k)]
                eig_vec_samp = [omega_eigs_vecs[i_b][i] for i in range(i_k-N_interp, i_k)]
                k = (k_2D_pts[i_k]-k_0_2D).norm()
            if k_samp!=None:
                coeff_lagrange = lagrange1d_coeffs(k_samp, k)
                
            if coeff_lagrange!=None:
                eig_xtrap = coeff_lagrange[0]*eig_samp[0]
                eig_vec_xtrap = coeff_lagrange[0]*eig_vec_samp[0]
                for i in range(1, len(k_samp)):
                    eig_xtrap += coeff_lagrange[i]*eig_samp[i]
                    eig_vec_xtrap += coeff_lagrange[i]*eig_vec_samp[i]
                return [eig_xtrap, eig_vec_xtrap]
            else:
                return [None, None]
            

        n_solve = n_bands            

        for i_k in range(i_k_start, n_pts):
            k_i_2D = k_2D_pts[i_k]
            k_i, theta_i = self.getMagAngleK2D(k_i_2D)
            
            print('i:', i_k,  'k_i_2D: ', k_i_2D, 'k_i:', k_i, 'theta_i:', theta_i)
            self.setEqskinit(k_i, theta_i)
            
            omega_eigs_vecs_ik = None
            omega_eigs_ik = None
            for i_b in range(n_solve):
                [eig_xtrap, eig_vec_xtrap] = extrap_eigs_eigsvecs(i_b, i_k, N_interp=3)

                omega_eigs_and_vecs = self.solveEigenvalues_RQ_NLCG(eig_0=eig_xtrap, eig_vec_0=eig_vec_xtrap, 
                    eigs_prev=omega_eigs_ik, eigs_vec_prev=omega_eigs_vecs_ik)
            
                omega_eigs_i_b = omega_eigs_and_vecs[0]
                omega_eigs_vecs_i_b = omega_eigs_and_vecs[1].T
                
                if eigs_i_k_arr==None:
                    omega_eigs_ik = [omega_eigs_i_b]
                    omega_eigs_vecs_ik = [omega_eigs_vecs_i_b]
                else:
                    omega_eigs_ik.append(omega_eigs_i_b)
                    omega_eigs_vecs_ik.append(omega_eigs_vecs_i_b)
            
            print('omega_eigs_ik: ', np.sqrt(omega_eigs_ik)/(2.0*math.pi))
            print('number of mat vec products: ', self.n_mat_vec_fft)
            #assert n_solve>=n_bands
        
            inds_omega_sorted = None
            omega_eigs_ik_abs = [abs(omega_eigs_ik[j]) for j in range(n_solve)]
            inds_omega_sorted = sorted(range(n_solve), key=lambda j: omega_eigs_ik_abs[j])
            print('inds_omega_sorted: ', inds_omega_sorted)

            if params["traceBands"]==True:
                ## problem with degenerate modes (orthogonality condition might not work)
                for i_b in range(n_bands):
                    vec_ib_prev = None
                    if i_k>0:
                        vec_ib_prev = omega_eigs_vecs[i_b][i_k-1]
                    elif eigs_vecs_init!=None:
                        vec_ib_prev = eigs_vecs_init[i_b]
                    else:
                        omega_eigs[i_b][i_k] = omega_eigs_ik[inds_omega_sorted[i_b]]
                        omega_eigs_vecs[i_b][i_k] = omega_eigs_vecs_ik[inds_omega_sorted[i_b]]
                        continue
                        
                    corr = [0.0]*n_solve
                    vec_ib_prev /= np.linalg.norm(vec_ib_prev)
                    for i_s in range(n_solve):
                        vec_is = omega_eigs_vecs_ik[i_s]
                        vec_is /= np.linalg.norm(vec_is)
                        corr[i_s] = np.abs(np.conjugate(vec_ib_prev).dot(vec_is))
                    
                    print(i_b, ' correlations: \n', corr)
                    ind_corr = 0
                    corr_max = corr[0]
                    for i_s in range(1, n_solve):
                        if corr[i_s] > corr_max:
                            ind_corr = i_s
                            corr_max = corr[i_s]
                    
                    if corr_max>0.9:
                        omega_eigs[i_b][i_k] = omega_eigs_ik[ind_corr]
                        omega_eigs_vecs[i_b][i_k] = omega_eigs_vecs_ik[ind_corr]
                    elif i_k>0:
                        omega_eigs[i_b][i_k] = None
                        omega_eigs_vecs[i_b][i_k] = omega_eigs_vecs[i_b][i_k-1]
                    else:
                        assert False
                        
                    
                    sq_omeg_2pi = omega_eigs[i_b].copy()
                    for j__ in range(len(sq_omeg_2pi)):
                        if sq_omeg_2pi[j__]!=None:
                            #print(sq_omeg_2pi[j__])
                            sq_omeg_2pi[j__] = np.sqrt(sq_omeg_2pi[j__])/(2.0*math.pi)
                    print(sq_omeg_2pi, '\n')
            else:
                for i_b in range(n_bands):
                    omega_eigs[i_b][i_k] = omega_eigs_ik[inds_omega_sorted[i_b]]
                    omega_eigs_vecs[i_b][i_k] = omega_eigs_vecs_ik[inds_omega_sorted[i_b]]

                    sq_omeg_2pi = omega_eigs[i_b].copy()
                    for j__ in range(len(sq_omeg_2pi)):
                        if sq_omeg_2pi[j__]!=None:
                            #print(sq_omeg_2pi[j__])
                            sq_omeg_2pi[j__] = np.sqrt(sq_omeg_2pi[j__])/(2.0*math.pi)
                    print(sq_omeg_2pi, '\n')

            if i_k>0 and i_k%5==0:
                saveDataToFile(i_k+1)
        saveDataToFile(n_pts)
        return [k_2D_pts, omega_eigs, omega_eigs_vecs]

    def getSymmetryPointsSquareLattice(self):
        Gamma = Point2D(0.0, 0.0)
        X = Point2D(self.vecRecip_vals[0].x/2.0, self.vecRecip_vals[0].y/2.0)
        M = X + Point2D(self.vecRecip_vals[1].x/2.0, self.vecRecip_vals[1].y/2.0)
        return [Gamma, X, M]

    def getSymmetryPointsHexagonalLattice(self):
        Gamma = Point2D(0.0, 0.0)
        v0, v1 = self.vecRecip_vals
        M = Point2D(v1.x/2.0, v1.y/2.0)
        cos_theta = v0&v1/(v0.norm()*v1.norm())
        ang = None
        if cos_theta>0.0:
            ang = math.acos(cos_theta)/2.0
        else:
            ang = math.acos(cos_theta)/4.0
        M_pp = Point2D(M.y, -M.x)   #perpandicular to M
        #M_pp = M_pp/M_pp.norm()
        K = M + M_pp*math.tan(ang)
        return [Gamma, M, K]

    def getMagAngleK2D(self, k_i_2D):
        k_i = k_i_2D.norm()
        theta_i = math.atan2(k_i_2D.y, k_i_2D.x)
        if abs(theta_i)<1.0e-15:
            theta_i = 0.0
        return [k_i, theta_i]


    def getBZLines(self):
        v0 = Point2D(self.vecRecip_vals[0].x, self.vecRecip_vals[0].y)
        v1 = Point2D(self.vecRecip_vals[1].x, self.vecRecip_vals[1].y)
        v2, v3 = v0+v1, v0-v1
        vecs_unit = [v0, v1, v2, v3, -v0, -v1, -v2, -v3]
        len_max = v0.norm()        
        if len_max<v1.norm():
            len_max = v1.norm()
        
        n_vec = len(vecs_unit)
        endpoints = []
        for i in range(n_vec):
            v_pp = Point2D(vecs_unit[i].y, -vecs_unit[i].x)
            v_pp = v_pp/v_pp.norm()
            v_start = vecs_unit[i]/2.0 - v_pp*len_max/2.0
            v_end   = vecs_unit[i]/2.0 + v_pp*len_max/2.0
            endpoints.append([v_start, v_end])
        return endpoints
        
        
        
        
            
            
            


