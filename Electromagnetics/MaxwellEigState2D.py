## MaxwellEigState2D.py
## solves for the eigenstates of 2D maxwell equations


__all__ = ["MaxEigState2DSolver", "materialDic"]


from sympy import *
from Electromagnetics.EMEqs_sym import *
from IPython.display import display, Math

import math
import numpy as np
from scipy import constants

from Electromagnetics.RecGrid2D import *
from Electromagnetics.Misc import *

import pickle
import os



from enum import Enum
class MaterialTypes(Enum):
    ems_sss = 1     #e:epsilon m:mu s:sigma, _s:scalar _t:tensorial
    ems_sst = 2     ## scalar epsilon and mu, tensorial sigma


materialDic = {'all_scalar':MaterialTypes.ems_sss, 'sigma_tensor':MaterialTypes.ems_sst}

MAT_NAME = 0
MAT_REC = 1
MAT_TYPE = 2
MAT_PARAMS = 3

FUN_EPS = 0
FUN_MU = 1
FUN_SIG = 2

MESH_REC = 0
MESH_NREF = 1

class MaxEigState2DSolver:
    """ it takes the domain rectangle, the structures (rectangles for the moment)
        and their material parameters as a function of frequency, and the rectangles
        for grid refinement and sets up an RecGrid2D structure for solving the 
        eigenvalues,
        It can sweep over the given frequency band
    """
    def __init__(self, materials, meshparams, vb=False):
        """ materials[i] = [name, rectangle, material_type, material_parameters]
                name: string
                rectangle: [point2D, point2D]
                material_type: picked from materialDic ('all_scalar', ...)
                material_params: [fun_epsilon, fun_mu, fun_sigma]
            depending on material type, if the type is scalar
            it includes three functions returning 3 scalars. 
            epsilon = fun_epsilon(freq)
            mu = fun_mu(freq)
            sigma = fun_sigma(freq)

                meshparams[i] = [rectangle, n_refines]
                n_refines: number of refinements
        """
        self.materials = materials
        self.meshparams = meshparams
        self.vbose = vb
        self.GetMaterialList()
        self.GenerateSymEquations()
        self.CreateGrid()
        
        
    def GetMaterialList(self):
        self.materialList = []
        print('possible material types : ', materialDic.keys())
        for mat in self.materials:
            print('material name: ', mat[MAT_NAME], '   type :', mat[MAT_TYPE])
            assert mat[MAT_TYPE] in list(materialDic.keys())
            if materialDic[mat[MAT_TYPE]] not in self.materialList:
                self.materialList.append(materialDic[mat[MAT_TYPE]])
            
    def setFreqRange(self, f0, f1, fstep):
        self.freq0 = f0
        self.freq1 = f1
        self.fstep = fstep
        
    def GenerateSymEquations(self):
        if MaterialTypes.ems_sss in self.materialList:
            if self.vbose:
                print("-"*20, "all scalar material", "-"*20)
            maxwell = MaxwellEqs_r()
            self.maxwell_sss = maxwell
            MF = maxwell.getMaxwellFaraday()
            MA = maxwell.getMaxwellAmper(hasSigma=True)
            MF = MF[0] - MF[1]
            MA = MA[0] - MA[1]
            if self.vbose:
                display(Math(latex(MF)))
                display(Math(latex(MA)))

            #assuming z propagation        
            MFz = maxwell.getZPropEqs(MF)
            MAz = maxwell.getZPropEqs(MA)
            if self.vbose:
                display(Math(latex(MFz)))
                display(Math(latex(MAz)))

            #supressing the harmonic term
            for i in range(MFz.shape[0]):
                for j in range(MFz.shape[1]):
                    MFz[i,j] = (MFz[i,j].doit()*exp(I*maxwell.k*maxwell.z)).simplify()
            for i in range(MAz.shape[0]):
                for j in range(MAz.shape[1]):
                    MAz[i,j] = (MAz[i,j].doit()*exp(I*maxwell.k*maxwell.z)).simplify()

            ## MF:maxwell faraday -- MA:maxwell ampere        
            if self.vbose:
                display(Math('MF = ' + latex(MFz.T)))
                display(Math('MA = ' + latex(MAz.T)))

            MFz = maxwell.ZPropEqs_ReplFunBySym(MFz)
            MAz = maxwell.ZPropEqs_ReplFunBySym(MAz)
            if self.vbose:
                display(Math('MF = ' + latex(MFz.T)))
                display(Math('MA = ' + latex(MAz.T)))

            k_0 = Symbol('k_0')
            self.k_0 = k_0          #@
            if self.vbose:
                print('normalizing to k_0')
            MFz[0] = (MFz[0]/(-I*k_0)).simplify()
            MFz[1] = (MFz[1]/(I*k_0)).simplify()
            MAz[0] = (MAz[0]/(-I*k_0)).simplify()
            MAz[1] = (MAz[1]/(I*k_0)).simplify()
            if self.vbose:
                display(Math('MF = ' + latex(MFz.T)))
                display(Math('MA = ' + latex(MAz.T)))

            MFz = MFz.subs(maxwell.k, 0)
            MAz = MAz.subs(maxwell.k, 0)
            if self.vbose:
                display(Math('MF = ' + latex(MFz.T)))
                display(Math('MA = ' + latex(MAz.T)))
            
            Eqs = [MFz[1], MFz[0], MAz[1], MAz[0], MAz[2], MFz[2]]
            self.Eqs_sss = Eqs     #@ sss: scalar_epsilon, scalar_mu, scalar_sigma
            if self.vbose:
                print('scalar epsilon mu sigma equations: ')
                display(Math('Eqs = ' + latex(Matrix(Eqs))))
                
            self.var_list = [maxwell.Ex, maxwell.Ey, maxwell.Hx, maxwell.Hy, maxwell.Ez, maxwell.Hz]
            if self.vbose:
                display(Math('Vars = ' + latex(Matrix(self.var_list))))


            ## regions boundary equations
            if self.vbose:
                print('disintegrating equations')
            Eqs_b = [self.DisintegrateLineaEquation(Eqs[0], [maxwell.Hy, maxwell.Ez], [maxwell.x, maxwell.y, maxwell.z]),
                     self.DisintegrateLineaEquation(Eqs[1], [maxwell.Hx, maxwell.Ez], [maxwell.x, maxwell.y, maxwell.z]),
     self.DisintegrateLineaEquation(Eqs[2], [maxwell.Ey, maxwell.Hz], [maxwell.x, maxwell.y, maxwell.z]),
     self.DisintegrateLineaEquation(Eqs[3], [maxwell.Ex, maxwell.Hz], [maxwell.x, maxwell.y, maxwell.z]),
     self.DisintegrateLineaEquation(Eqs[4], [maxwell.Ez, maxwell.Hx, maxwell.Hy], [maxwell.x, maxwell.y, maxwell.z]),
     self.DisintegrateLineaEquation(Eqs[5], [maxwell.Hz, maxwell.Ex, maxwell.Ey], [maxwell.x, maxwell.y, maxwell.z])]

            """
            Eqs_b = [[Eqs[0].subs(maxwell.Ez, 0).doit().simplify(), Eqs[0].subs(maxwell.Hy, 0).simplify()],
                         [Eqs[1].subs(maxwell.Ez, 0).doit().simplify(), Eqs[1].subs(maxwell.Hx, 0).simplify()],
                         [Eqs[2].subs(maxwell.Hz, 0).doit().simplify(), Eqs[2].subs(maxwell.Ey, 0).simplify()],
                         [Eqs[3].subs(maxwell.Hz, 0).doit().simplify(), Eqs[3].subs(maxwell.Ex, 0).simplify()],
                         [Eqs[4].subs([(maxwell.Hx, 0), (maxwell.Hy, 0)]).doit().simplify(), Eqs[4].subs([(maxwell.Ez, 0),\
              (maxwell.Hy, 0)]).subs(Derivative(0, maxwell.x), 0).simplify(), Eqs[4].subs([(maxwell.Ez, 0), (maxwell.Hx, 0)])\
              .subs(Derivative(0, maxwell.y), 0).simplify()],
                         [Eqs[5].subs([(maxwell.Ex, 0), (maxwell.Ey, 0)]).doit().simplify(), Eqs[5].subs([(maxwell.Hz, 0),\
              (maxwell.Ey, 0)]).subs(Derivative(0, maxwell.x), 0).simplify(), Eqs[5].subs([(maxwell.Hz, 0), (maxwell.Ex, 0)])\
              .subs(Derivative(0, maxwell.y), 0).simplify()]]
            """

            self.Eqs_b_sss = Eqs_b      #@
            if self.vbose:
                for i in range(len(Eqs_b)):
                    display(Math('Eqs_b[{}] = '.format(i) + latex(Matrix(Eqs_b[i]).T)))
        
        if MaterialTypes.ems_sst in self.materialList:
            if self.vbose:
                print("-"*20, "Sigma tensor material", "-"*20)
            maxwell = MaxwellEqs_r(sigma_type='tensor')
            self.maxwell_sst = maxwell
            MF = maxwell.getMaxwellFaraday()
            MA = maxwell.getMaxwellAmper(hasSigma=True)
            MF = MF[0] - MF[1]
            MA = MA[0] - MA[1]
            if self.vbose:
                display(Math(latex(MF)))
                display(Math(latex(MA)))

            #assuming z propagation        
            MFz = maxwell.getZPropEqs(MF)
            MAz = maxwell.getZPropEqs(MA)
            if self.vbose:
                display(Math(latex(MFz)))
                display(Math(latex(MAz)))

            #supressing the harmonic term
            for i in range(MFz.shape[0]):
                for j in range(MFz.shape[1]):
                    MFz[i,j] = (MFz[i,j].doit()*exp(I*maxwell.k*maxwell.z)).simplify()
            for i in range(MAz.shape[0]):
                for j in range(MAz.shape[1]):
                    MAz[i,j] = (MAz[i,j].doit()*exp(I*maxwell.k*maxwell.z)).simplify()

            ## MF:maxwell faraday -- MA:maxwell ampere        
            if self.vbose:
                display(Math('MF = ' + latex(MFz.T)))
                display(Math('MA = ' + latex(MAz.T)))

            MFz = maxwell.ZPropEqs_ReplFunBySym(MFz)
            MAz = maxwell.ZPropEqs_ReplFunBySym(MAz)
            if self.vbose:
                display(Math('MF = ' + latex(MFz.T)))
                display(Math('MA = ' + latex(MAz.T)))

            k_0 = Symbol('k_0')
            self.k_0 = k_0          #@
            if self.vbose:
                print('normalizing to k_0')
            MFz[0] = (MFz[0]/(-I*k_0)).simplify()
            MFz[1] = (MFz[1]/(I*k_0)).simplify()
            MAz[0] = (MAz[0]/(-I*k_0)).simplify()
            MAz[1] = (MAz[1]/(I*k_0)).simplify()
            if self.vbose:
                display(Math('MF = ' + latex(MFz.T)))
                display(Math('MA = ' + latex(MAz.T)))

            MFz = MFz.subs(maxwell.k, 0)
            MAz = MAz.subs(maxwell.k, 0)
            if self.vbose:
                display(Math('MF = ' + latex(MFz.T)))
                display(Math('MA = ' + latex(MAz.T)))
            
            Eqs = [MFz[1], MFz[0], MAz[1], MAz[0], MAz[2], MFz[2]]
            self.Eqs_sst = Eqs     #@ sss: scalar_epsilon, scalar_mu, scalar_sigma
            if self.vbose:
                print('scalar epsilon mu sigma equations: ')
                display(Math('Eqs = ' + latex(Matrix(Eqs))))
                
            self.var_list = [maxwell.Ex, maxwell.Ey, maxwell.Hx, maxwell.Hy, maxwell.Ez, maxwell.Hz]
            if self.vbose:
                display(Math('Vars = ' + latex(Matrix(self.var_list))))


            ## regions boundary equations
            if self.vbose:
                print('disintegrating equations')
            Eqs_b = [self.DisintegrateLineaEquation(Eqs[0], [maxwell.Hy, maxwell.Ez], [maxwell.x, maxwell.y, maxwell.z]),
                     self.DisintegrateLineaEquation(Eqs[1], [maxwell.Hx, maxwell.Ez], [maxwell.x, maxwell.y, maxwell.z]),
     self.DisintegrateLineaEquation(Eqs[2], [maxwell.Ex, maxwell.Ey, maxwell.Ez, maxwell.Hz], [maxwell.x, maxwell.y, maxwell.z]),
     self.DisintegrateLineaEquation(Eqs[3], [maxwell.Ex, maxwell.Ey, maxwell.Ez, maxwell.Hz], [maxwell.x, maxwell.y, maxwell.z]),
     self.DisintegrateLineaEquation(Eqs[4], [maxwell.Ex, maxwell.Ey, maxwell.Ez, maxwell.Hx, maxwell.Hy], [maxwell.x, maxwell.y, maxwell.z]),
     self.DisintegrateLineaEquation(Eqs[5], [maxwell.Hz, maxwell.Ex, maxwell.Ey], [maxwell.x, maxwell.y, maxwell.z])]

            self.Eqs_b_sst = Eqs_b      #@
            if self.vbose:
                for i in range(len(Eqs_b)):
                    display(Math('Eqs_b[{}] = '.format(i) + latex(Matrix(Eqs_b[i]).T)))

    def DisintegrateLineaEquation(self, Eq, eq_vars, ind_vars):
        """ eq_vars: list of vasiables inside the equation
            ind_vars: list of indipendant variables (x, y, z...)
        """
        eqs_dis = []
        for i in range(len(eq_vars)):
            Eq_i = Eq
            for j in range(len(eq_vars)):
                if j != i:
                    Eq_i = Eq_i.subs(eq_vars[j], 0)
            for j in range(len(ind_vars)):
                Eq_i = Eq_i.subs(Derivative(0, ind_vars[j]), 0)
            Eq_i = Eq_i.simplify()
            eqs_dis.append(Eq_i)
        return eqs_dis

        
    def CreateGrid(self):
        rec_D = self.materials[0][MAT_REC]        ## main domain region
        p0_D = rec_D[0]
        p1_D = rec_D[1]
        W_x = p1_D.x - p0_D.x
        W_y = p1_D.y - p0_D.y
        assert W_x>0.0 and W_y>0.0
        
        self.rg = RG2D(0.0, W_x, 0.0, W_y, W_x/10.0, W_y/10.0)
        
        meshparams = self.meshparams
        rg = self.rg
        if self.vbose:
            print("refining mesh")
        for i in range(len(meshparams)):
            rec_G = meshparams[i][MESH_REC]    ## grid rectangle
            n_refine = meshparams[i][MESH_NREF]
            for n in range(n_refine):
                rg.RefineRectangle(rec_G)
        
        if self.vbose:
            print("refining materials")
        for i in range(1, len(self.materials)):
            rect_m = self.materials[i][MAT_REC]    ## material rectangle
            #rg.RefineRectangle(rect_m)
            rg.RefineThinRectangle(rect_m)
            if self.vbose:
                print("refining border of {}".format(self.materials[i][MAT_NAME]))
            rg.RefineRectangleBorders(rect_m)
            rg.RefineRectangleBorders(rect_m)
            rg.RefineRectangleBorders(rect_m)
        if self.vbose:        
            print('n_cells : ', len(rg.cells), ' n_sides: ', len(rg.sides), ' n_nodes: ', len(rg.nodes))


    def SetupEquations(self, freq):
        self.Eqs_list = [None]*len(self.materials)      #@
        self.Eqs_b_list = [None]*len(self.materials)      #@

        omega_0 = 2.0*math.pi*freq
        eps_0 = constants.epsilon_0
        mu_0 = constants.mu_0
        k_0_val = omega_0*sqrt(mu_0*eps_0)
        k_0 = self.k_0

        for i_mat in range(len(self.materials)):
            mat = self.materials[i_mat]
            if materialDic[mat[MAT_TYPE]]==MaterialTypes.ems_sss:
                maxwell = self.maxwell_sss
                mat_name = mat[MAT_NAME]
                fun_epsilon = mat[MAT_PARAMS][FUN_EPS]
                fun_mu = mat[MAT_PARAMS][FUN_MU]
                fun_sig = mat[MAT_PARAMS][FUN_SIG]
                eps = fun_epsilon(freq)
                mu = fun_mu(freq)
                sig = fun_sig(freq)
                
                Eqs = self.Eqs_sss
                Eqs_i = [Eqs[i].subs([(maxwell.mu, mu), (maxwell.eps, eps), (maxwell.sigma, sig), \
                    (maxwell.omega, omega_0), (k_0, k_0_val)]).simplify() for i in range(len(Eqs))]
                if self.vbose:
                    display(Math('Eqs_\\text{{ {} }} = '.format(mat_name) + latex(Matrix(Eqs_i))))
                
                self.Eqs_list[i_mat] = Eqs_i
 
                ## boundary equations
                Eqs_b = self.Eqs_b_sss
                if self.vbose:
                    print('Eq1 boundary equations ')
                Eqs_b_i = [[Eqs_b[i][j].subs([(maxwell.mu, mu), (maxwell.eps, eps), (maxwell.sigma, sig), \
                    (maxwell.omega, omega_0), (k_0, k_0_val)]) for j in range(len(Eqs_b[i]))] for i in range(len(Eqs_b))]
                if self.vbose:
                    for i in range(len(Eqs_b_i)):
                        display(Math('EQs_\\text{{ {} }}-b[{}] = '.format(mat_name, i) + latex(Matrix(Eqs_b_i[i]).T)))
                        
                self.Eqs_b_list[i_mat] = Eqs_b_i
            elif materialDic[mat[MAT_TYPE]]==MaterialTypes.ems_sst:
                maxwell = self.maxwell_sst
                mat_name = mat[MAT_NAME]
                fun_epsilon = mat[MAT_PARAMS][FUN_EPS]
                fun_mu = mat[MAT_PARAMS][FUN_MU]
                fun_sig = mat[MAT_PARAMS][FUN_SIG]
                eps = fun_epsilon(freq)
                mu = fun_mu(freq)
                sig = fun_sig(freq)
                
                print(sig)
                
                sig_xx = sig[0][0]
                sig_xy = sig[0][1]
                sig_xz = sig[0][2]
                sig_yx = sig[1][0]
                sig_yy = sig[1][1]
                sig_yz = sig[1][2]
                sig_zx = sig[2][0]
                sig_zy = sig[2][1]
                sig_zz = sig[2][2]
                
                Eqs = self.Eqs_sst
                Eqs_i = [Eqs[i].subs([(maxwell.mu, mu), (maxwell.eps, eps), \
                    (maxwell.sigma_xx, sig_xx), (maxwell.sigma_xy, sig_xy), (maxwell.sigma_xz, sig_xz),\
                    (maxwell.sigma_yx, sig_yx), (maxwell.sigma_yy, sig_yy), (maxwell.sigma_yz, sig_yz),\
                    (maxwell.sigma_zx, sig_zx), (maxwell.sigma_zy, sig_zy), (maxwell.sigma_zz, sig_zz),\
                    (maxwell.omega, omega_0), (k_0, k_0_val)]).simplify() for i in range(len(Eqs))]
                if self.vbose:
                    display(Math('Eqs_\\text{{ {} }} = '.format(mat_name) + latex(Matrix(Eqs_i))))
                
                self.Eqs_list[i_mat] = Eqs_i
 
                ## boundary equations
                Eqs_b = self.Eqs_b_sst
                if self.vbose:
                    print('Eq1 boundary equations ')
                Eqs_b_i = [[Eqs_b[i][j].subs([(maxwell.mu, mu), (maxwell.eps, eps), \
                    (maxwell.sigma_xx, sig_xx), (maxwell.sigma_xy, sig_xy), (maxwell.sigma_xz, sig_xz),\
                    (maxwell.sigma_yx, sig_yx), (maxwell.sigma_yy, sig_yy), (maxwell.sigma_yz, sig_yz),\
                    (maxwell.sigma_zx, sig_zx), (maxwell.sigma_zy, sig_zy), (maxwell.sigma_zz, sig_zz),\
                    (maxwell.omega, omega_0), (k_0, k_0_val)]) for j in range(len(Eqs_b[i]))] for i in range(len(Eqs_b))]
                if self.vbose:
                    for i in range(len(Eqs_b_i)):
                        display(Math('EQs_\\text{{ {} }}-b[{}] = '.format(mat_name, i) + latex(Matrix(Eqs_b_i[i]).T)))
                        
                self.Eqs_b_list[i_mat] = Eqs_b_i
            else:
                raise NotImplementedError()
        
        
    def AddEquationsToGrid(self, freq, firstPass=False):
        cells_nm_list = [None]*len(self.materials) 
        for i_mat in range(len(self.materials)):
            mat = self.materials[i_mat]
            if materialDic[mat[MAT_TYPE]]==MaterialTypes.ems_sss:
                mat_name = mat[MAT_NAME]
                region_nm = mat_name    ## region name
                self.rg.AddExpression(region_nm, self.Eqs_list[i_mat], self.var_list)
            elif materialDic[mat[MAT_TYPE]]==MaterialTypes.ems_sst:
                mat_name = mat[MAT_NAME]
                region_nm = mat_name    ## region name
                self.rg.AddExpression(region_nm, self.Eqs_list[i_mat], self.var_list)
            else:
                raise NotImplementedError()

        for i_mat in range(len(self.materials)):
            mat = self.materials[i_mat]
            if materialDic[mat[MAT_TYPE]]==MaterialTypes.ems_sss:
                mat_name = mat[MAT_NAME]
                region_nm = mat_name    ## region name
                Eqs_b = self.Eqs_b_list[i_mat]
                eq_list =  [{'zero':Float(0), 'one':(Eqs_b[0][0]+Eqs_b[0][1])/2.0, 'two':(Eqs_b[0][0]+Eqs_b[0][1])},
                            {'zero':Float(0), 'one':(Eqs_b[1][0]+Eqs_b[1][1])/2.0, 'two':(Eqs_b[1][0]+Eqs_b[1][1])},
                            {'zero':Float(0), 'one':(Eqs_b[2][0]/2.0+Eqs_b[2][1]), 'two':(Eqs_b[2][0]+Eqs_b[2][1])},
                            {'zero':Float(0), 'one':(Eqs_b[3][0]/2.0+Eqs_b[3][1]), 'two':(Eqs_b[3][0]+Eqs_b[3][1])},
                            {'zero':Float(0), 'one':(Eqs_b[4][0]/4.0), 'two':(Eqs_b[4][0]/2.0), \
                                'two_x':(Eqs_b[4][1]+Eqs_b[4][2]/2.0), 'two_y':(Eqs_b[4][1]/2.0+Eqs_b[4][2]), \
                                'three':(Eqs_b[4][0]*3.0/4.0+Eqs_b[4][1]+Eqs_b[4][2]), \
                                'four':(Eqs_b[4][0]+Eqs_b[4][1]+Eqs_b[4][2])},
                           {}]
                var_type = ['SX', 'SY', 'SY', 'SX', 'N', 'C']
                self.rg.SetRegionBoundaryEquations(region_nm, eq_list, self.var_list, var_type)
            elif materialDic[mat[MAT_TYPE]]==MaterialTypes.ems_sst:
                mat_name = mat[MAT_NAME]
                region_nm = mat_name    ## region name
                Eqs_b = self.Eqs_b_list[i_mat]
                eq_list =  [{'zero':Float(0), 'one':(Eqs_b[0][0]+Eqs_b[0][1])/2.0, 'two':(Eqs_b[0][0]+Eqs_b[0][1])},
                            {'zero':Float(0), 'one':(Eqs_b[1][0]+Eqs_b[1][1])/2.0, 'two':(Eqs_b[1][0]+Eqs_b[1][1])},
                            {'zero':Float(0), 'one':(Eqs_b[2][0]/2.0+Eqs_b[2][1]/2.0+Eqs_b[2][2]/2.0+Eqs_b[2][3]),\
                                 'two':(Eqs_b[2][0]+Eqs_b[2][1]+Eqs_b[2][2]+Eqs_b[2][3])},
                            {'zero':Float(0), 'one':(Eqs_b[3][0]/2.0+Eqs_b[3][1]/2.0+Eqs_b[3][2]/2.0+Eqs_b[3][3]),\
                                 'two':(Eqs_b[3][0]+Eqs_b[3][1]+Eqs_b[3][2]+Eqs_b[3][3])},
                            {'zero':Float(0), 'one':(Eqs_b[4][0]/4.0+Eqs_b[4][1]/4.0+Eqs_b[4][2]/4.0),\
                                'two':(Eqs_b[4][0]/2.0+Eqs_b[4][1]/2.0+Eqs_b[4][2]/2.0), \
                                'two_x':(Eqs_b[4][3]+Eqs_b[4][4]/2.0), 'two_y':(Eqs_b[4][3]/2.0+Eqs_b[4][4]), \
                                'three':(Eqs_b[4][0]*3.0/4.0+Eqs_b[4][1]*3.0/4.0+Eqs_b[4][2]*3.0/4.0+Eqs_b[4][3]+Eqs_b[4][4]), \
                                'four':(Eqs_b[4][0]+Eqs_b[4][1]+Eqs_b[4][2]+Eqs_b[4][3]+Eqs_b[4][4])},
                           {}]
                var_type = ['SX', 'SY', 'SY', 'SX', 'N', 'C']
                self.rg.SetRegionBoundaryEquations(region_nm, eq_list, self.var_list, var_type)
            else:
                raise NotImplementedError()
        
        
        if firstPass:   ## assuming for other passes the grid is unchanged        
            for i_mat in range(len(self.materials)):
                mat = self.materials[i_mat]
                mat_name = mat[MAT_NAME]
                region_nm = mat_name    ## region name
                mat_rec = mat[MAT_REC]
                if i_mat==0:
                    cells_nm = list(range(len(self.rg.cells)))
                    cells_nm_list[i_mat] = cells_nm
                else:
                    cells_nm, nodes_nm, sidesX_nm, sidesY_nm = self.rg.FindElementsMarkedByRectangle(mat_rec[0], mat_rec[1])
                    cells_nm_list[i_mat] = cells_nm

            cells_nm_list = self.rg.RegionsRemoveSharedCells(cells_nm_list)

            rg = self.rg
            maxwell = self.maxwell_sss
            for i_mat in range(len(self.materials)):
                mat = self.materials[i_mat]
                mat_name = mat[MAT_NAME]
                region_nm = mat_name    ## region name
                cells_nm = cells_nm_list[i_mat]
                nodes_nm, sidesX_nm, sidesY_nm = self.rg.Cells_GetNodesSide(cells_nm)

                rg.AssignVariableToSides(maxwell.Ex, sidesX_nm, region_nm, orient='X')
                rg.AssignVariableToSides(maxwell.Ey, sidesY_nm, region_nm, orient='Y')
                rg.AssignVariableToSides(maxwell.Hx, sidesY_nm, region_nm, orient='Y')
                rg.AssignVariableToSides(maxwell.Hy, sidesX_nm, region_nm, orient='X')
                rg.AssignVariableToNodes(maxwell.Ez, nodes_nm, region_nm)
                rg.AssignVariableToCells(maxwell.Hz, cells_nm, region_nm)

            rg.SetTotalIndicesForVariables()
            if self.vbose:
                print('rg.CellsToIndexTotalShared: ', len(rg.CellsToIndexTotalShared))
                print('rg.SidesToIndexTotalShared: ', len(rg.SidesToIndexTotalShared))
                print('rg.NodesToIndexTotalShared: ', len(rg.NodesToIndexTotalShared))
            
            self.cells_nm_list = cells_nm_list

    
    def GetCellsNodesAndSidesForMaterial(self, material_name):
        for i_mat in range(len(self.materials)):
            mat = self.materials[i_mat]
            mat_name = mat[MAT_NAME]
            print('material {}: '.format(i_mat), mat_name)
            if mat_name == material_name:
                cells_nm = self.cells_nm_list[i_mat]
                nodes_nm, sidesX_nm, sidesY_nm = self.rg.Cells_GetNodesSide(cells_nm)
                return [cells_nm, nodes_nm, sidesX_nm, sidesY_nm]
        return None
        

        
    def ConstructMatrixAndSolve(self, n_eig, lambda_0, vec_0=None, maxiter=None, tol=0):
        rg = self.rg
        maxwell = self.maxwell_sss
        tic()
        A_coo = rg.ConstructInitialMatrix(save=True)
        if self.vbose:
            print("initial matrix nnz", A_coo.nnz, toc())
        tic()
        A_coo = rg.ConstructInitialMatrix_Shared(A_coo, save=True)
        if self.vbose:
            print("initial matrix shared nnz", A_coo.nnz, toc())
        tic()
        A_coo = rg.ApplyDirichletBoundaryCondition(A_coo, save=True)
        if self.vbose:
            print("final matrix (application of boundary conditions) nnz: ", A_coo.nnz, toc())
        self.A_coo_final = A_coo
        tic()
        A_fin = rg.EliminateVarible_Direct(A_coo, [maxwell.Ez, maxwell.Hz])
        if self.vbose:
            print("eliminating z components", toc())
        tic()
        k_eigs, k_vecs = rg.SolveEigenvalues(A_fin, n_eig, lambda_0=lambda_0, vec_0=vec_0, tol=tol, maxiter=maxiter)
        if self.vbose:
            print("A_reduced nnz: ", A_fin.nnz)
            print("solving eigenvalues", toc())
        if self.vbose:
            print(k_eigs)
            print(k_vecs.shape)
        return [k_eigs, k_vecs]
    
    
    def GetEtInCellsFromFinalSolutionVector(self, k_vec, A_coo_final):
        rg = self.rg
        maxwell = self.maxwell_sss
        tic()
        rg.GetFromSolutionIndexToInitialIndex(k_vec, A_coo_final)
        if self.vbose:
            print('GetFromSolutionIndexToInitialIndex : ', toc())

        tic()
        exy_cell = rg.InterpolateVariablesValuesToCellCenters([maxwell.Ex, maxwell.Ey])
        ex = np.array(exy_cell[0])
        ey = np.array(exy_cell[1])
        if self.vbose:
            print('InterpolateVariablesValuesToCellCenters : ', toc())

        tic()
        ex_sq = rg.CalculateSurfaceIntegral(ex, np.conjugate(ex))
        ey_sq = rg.CalculateSurfaceIntegral(ey, np.conjugate(ey))
        exy_norm = np.sqrt(ex_sq + ey_sq)
        ex /= exy_norm
        ey /= exy_norm
        if self.vbose:
                print('CalculateSurfaceIntegral : ', toc())
        return [ex, ey]


    def GetHtInCellsFromFinalSolutionVector(self, k_vec, A_coo_final):
        rg = self.rg
        maxwell = self.maxwell_sss
        tic()
        rg.GetFromSolutionIndexToInitialIndex(k_vec, A_coo_final)
        if self.vbose:
            print('GetFromSolutionIndexToInitialIndex : ', toc())

        tic()
        hxy_cell = rg.InterpolateVariablesValuesToCellCenters([maxwell.Hx, maxwell.Hy])
        hx = np.array(hxy_cell[0])
        hy = np.array(hxy_cell[1])
        if self.vbose:
            print('InterpolateVariablesValuesToCellCenters : ', toc())

        tic()
        hx_sq = rg.CalculateSurfaceIntegral(hx, np.conjugate(hx))
        hy_sq = rg.CalculateSurfaceIntegral(hy, np.conjugate(hy))
        hxy_norm = np.sqrt(hx_sq + hy_sq)
        hx /= hxy_norm
        hy /= hxy_norm
        if self.vbose:
                print('CalculateSurfaceIntegral : ', toc())
        return [hx, hy]

        
    
    def GetTransverseCrossProductEt1Ht2(self, k_vec_1, k_vec_2, A_coo_final, conjHt2=False):
        rg = self.rg
        maxwell = self.maxwell_sss
        conj = ''
        if conjHt2:
            conj = 'np.conjugate'

        rg.GetFromSolutionIndexToInitialIndex(k_vec_1, A_coo_final)
        exy_cell = rg.InterpolateVariablesValuesToCellCenters([maxwell.Ex, maxwell.Ey])
        ex_1 = np.array(exy_cell[0])
        ey_1 = np.array(exy_cell[1])

        ex1_sq = rg.CalculateSurfaceIntegral(ex_1, np.conjugate(ex_1))
        ey1_sq = rg.CalculateSurfaceIntegral(ey_1, np.conjugate(ey_1))
        exy1_norm = np.sqrt(ex1_sq + ey1_sq)
        ex_1 /= exy1_norm
        ey_1 /= exy1_norm

        rg.GetFromSolutionIndexToInitialIndex(k_vec_2, A_coo_final)
        hxy_cell = rg.InterpolateVariablesValuesToCellCenters([maxwell.Hx, maxwell.Hy])
        hx_2 = np.array(hxy_cell[0])
        hy_2 = np.array(hxy_cell[1])

        hx2_sq = rg.CalculateSurfaceIntegral(hx_2, np.conjugate(hx_2))
        hy2_sq = rg.CalculateSurfaceIntegral(hy_2, np.conjugate(hy_2))
        hxy2_norm = np.sqrt(hx2_sq + hy2_sq)
        hx_2 /= hxy2_norm
        hy_2 /= hxy2_norm
        
        ex1_hy2 = rg.CalculateSurfaceIntegral(ex_1, eval(conj + '(hy_2)'))
        ey1_hx2 = rg.CalculateSurfaceIntegral(ey_1, eval(conj + '(hx_2)'))
        
        e1_h2 = ex1_hy2 - ey1_hx2
        print(ex1_hy2, ey1_hx2)
        
        return e1_h2
        
        
    def Solve(self, freq, n_eig, lambda_0, vec_0=None, maxiter=None, tol=0, firstPass=False):
        self.SetupEquations(freq)
        self.AddEquationsToGrid(freq, firstPass=firstPass)
        return self.ConstructMatrixAndSolve(n_eig=n_eig, lambda_0=lambda_0, vec_0=vec_0, tol=tol, maxiter=maxiter)
        
        
    def sweepFrequency(self, eig_val0, eig_vec0, f_start, f_end, f_step, n_eig=5, maxiter=None, tol=0,
            stop_cond=None, record_data=False):
        self.setFreqRange(f_start, f_end, f_step)
        freq_sweep = ParamSweeper()
        freq_sweep.setLimitPoints(f_start, f_end)
        freq_sweep.setMaxStepSize(f_step)
        freq_sweep.addNext()
        eigvar_vals_start = [eig_val0]
        freq_sweep.setLastY(eigvar_vals_start)
        if self.vbose:
            freq_sweep.printXY()
        
        x_null_prev = eig_vec0
        x_null_prev /= np.linalg.norm(x_null_prev)
        
        x_null_sweep = []
        if record_data:
            x_null_sweep.append(x_null_prev)
            
        n_diverge = 0
        n_diverg_stat = []
        while(freq_sweep.addNext()):
            freq = freq_sweep.getLastX()
            k_init = freq_sweep.extrapolateLast()[0]
            
            k_eigs, k_vecs = self.Solve(freq, n_eig=n_eig, lambda_0=k_init, tol=tol, 
                    maxiter=maxiter, firstPass=False)
            
            solvefailure = False
            if solvefailure:
                freq_sweep.refineLast()
                n_diverge += 1
                if n_diverge>3:
                    freq_sweep.decreaseMaxStepSize()
                    n_diverge = 0
                if self.vbose:
                    print('k_eigs: ', k_eigs)
                    print('making refinement...')
                    print('=-'*50)
                print('making refinement...')
                continue
            
            ind_corr = -1
            ind_x_corr = 0
            for i in range(n_eig):
                x_null = k_vecs[:,i]
                x_null /= np.linalg.norm(x_null)

                x_corr = x_null.dot(np.conj(x_null_prev))
            
                if self.vbose:
                    print('x_correlation: ', x_corr, '    abs(x_correlation): ', abs(x_corr))
                k_corr_r = abs(k_eigs[i].real - k_init.real)/max(abs(k_eigs[i].real), abs(k_init.real))
                k_corr_i = abs(k_eigs[i].imag - k_init.imag)/max(abs(k_eigs[i].imag), abs(k_init.imag))
                if self.vbose:
                    print('k_correlation_r: ', k_corr_r, '    k_correlation_i: ', k_corr_i)
                if abs(k_init.real/k_init.imag)>1.0e4:
                    k_corr_i = 0.0
                    if self.vbose:
                        print('**** k_correlation_i set to zero ****')
                if ind_corr<0:
                    if abs(x_corr)>0.9 and abs(k_corr_r)<0.1 and abs(k_corr_i)<0.1:
                        ind_corr = i
                        ind_x_corr = x_corr
                else:
                    if abs(x_corr)>0.9 and abs(k_corr_r)<0.1 and abs(k_corr_i)<0.1:
                        if x_corr > ind_x_corr:
                            ind_corr = i
                            ind_x_corr = x_corr
            
            n_diverg_stat.append(n_diverge)
            if ind_corr<0:
                freq_sweep.refineLast()
                n_diverge += 1
                if n_diverge>3:
                    freq_sweep.decreaseMaxStepSize()
                    n_diverge = 0
                if self.vbose:
                    print('making refinement...')
                    print('=-'*50)
                print('making refinement...')
                continue
            
            x_null = k_vecs[:,ind_corr]
            x_null /= np.linalg.norm(x_null)

            x_null_prev = x_null
            freq_sweep.setLastY([k_eigs[ind_corr]])
            if self.vbose:
                freq_sweep.printXY()
                print('--'*50)

            if record_data:
                x_null_sweep.append(x_null)
                
            if stop_cond!=None:
                cond__0 = "Re(k)<"
                if stop_cond.startswith(cond__0):
                    cond__1 = complex(stop_cond[len(cond__0):])
                    assert cond__1.imag == 0.0
                    if k_eigs[ind_corr].real < cond__1:
                        break
                    

        print('divergence stats: ', n_diverg_stat)
        return [freq_sweep, x_null_sweep]
        
        
    def sweepFrequencyMultiple(self, eig_val0, eig_vec0, f_start, f_end, f_step, n_eig=10, maxiter=None, tol=0,
            interp_ord=3, stop_cond=None, record_data=False, saveToFile=None, readFromFile=None):
        n_modes = len(eig_val0)
        assert len(eig_vec0) == n_modes
        assert n_eig>n_modes
        self.setFreqRange(f_start, f_end, f_step)
        freq_sweep = [0]*n_modes
        for i in range(n_modes):
            freq_sweep[i] = ParamSweeper()
            freq_sweep[i].setLimitPoints(f_start, f_end)
            freq_sweep[i].setMaxStepSize(f_step)
            freq_sweep[i].addNext()     #1st element
            eigvar_vals_start = [eig_val0[i]]
            freq_sweep[i].setLastY(eigvar_vals_start)
            if self.vbose:
                freq_sweep[i].printXY()
        
        x_null_prev = eig_vec0
        for i in range(n_modes):
            x_null_prev[i] /= np.linalg.norm(x_null_prev[i])
        
        x_null_sweep = []
        if record_data:
            x_null_sweep.append(x_null_prev)
            
        n_diverge = 0
        n_diverg_stat = []
        
        ind_lead = 0
        mode_active = [True]*n_modes  #0:active 1:terminated
        
        if saveToFile!=None:
            f_name = saveToFile
            f_dir = os.path.dirname(f_name)
            if not os.path.exists(f_dir):
                os.makedirs(f_dir)
            f = open(f_name, 'wb')
            params_dic = {'n_modes': n_modes, 'mode_active':mode_active, 
                'ind_lead':ind_lead, 'n_eig':n_eig, 'x_null_prev': x_null_prev}
            pickle.dump(params_dic, f)
            pickle.dump(freq_sweep, f)
            f.close()
        
        if readFromFile!=None:
            f_name = readFromFile
            f = open(f_name, 'rb')
            params_dic = pickle.load(f)
            freq_sweep = pickle.load(f)
            n_modes = params_dic['n_modes']
            mode_active = params_dic['mode_active']
            ind_lead = params_dic['ind_lead']
            n_eig = params_dic['n_eig']
            x_null_prev = params_dic['x_null_prev']

        while(True):

            final_f_reached = False
            for i in range(n_modes):
                if mode_active[i]:
                    fin_f_ = not freq_sweep[i].addNext()
                    if fin_f_:
                        final_f_reached = True
                        mode_active[i] = False
                        freq_sweep[i].RemoveLast()

            if final_f_reached:
                #-- reset leader
                if mode_active[ind_lead]==False:
                    for i in range(n_modes):
                        if mode_active[i]:
                            ind_lead = i
                            break
                if mode_active[ind_lead]==False:
                    for i in range(n_modes):
                        if freq_sweep[i].Y[-1]==None:
                            ind_lead = i
                            freq_sweep[i].RemoveLast()
                            mode_active[i] = True
                            n_eig = n_modes     ##to avoid long convergence times
                                                # only one eigenvalue is looed for in this case
                            break
                #---
                for i in range(n_modes):
                    if mode_active[i]:
                        freq_sweep[i].addNext()
            
            if mode_active[ind_lead]==False:
                break

            freq = freq_sweep[ind_lead].getLastX()
            k_init = freq_sweep[ind_lead].extrapolateLast(n=interp_ord)[0]
            
            n_mode_inactive = 0
            for i in range(n_modes):
                if mode_active[i]==False:
                    n_mode_inactive += 1
            
            n_eig_act = n_eig - n_mode_inactive
            
            if self.vbose:
                print('frequency:', freq, '  k_init:', k_init)
            
            k_eigs, k_vecs = self.Solve(freq, n_eig=n_eig_act, lambda_0=k_init, tol=tol, 
                    maxiter=maxiter, firstPass=False)
            
            solvefailure = False
            if solvefailure:
                for i in range(n_modes):
                    if mode_active[i]:
                        freq_sweep[i].refineLast()
                n_diverge += 1
                if n_diverge>3:
                    for i in range(n_modes):
                        if mode_active[i]:
                            freq_sweep[i].decreaseMaxStepSize()
                        n_diverge = 0
                if self.vbose:
                    print('k_eigs: ', k_eigs)
                    print('making refinement...')
                    print('=-'*50)
                print('making refinement...')
                continue
            
            ind_corr = [-1]*n_modes
            ind_x_corr = [-1]*n_modes
            for i in range(n_eig_act):
                x_null = k_vecs[:,i]
                x_null /= np.linalg.norm(x_null)

                for j in range(n_modes):
                    if mode_active[j]==False:
                        continue
                    x_corr = x_null.dot(np.conj(x_null_prev[j]))
                    if self.vbose:
                        print(i, j, 'x_correlation: ', x_corr, '    abs(x_correlation): ', abs(x_corr))

                    k_init_j = freq_sweep[j].extrapolateLast(n=interp_ord)[0]

                    k_corr_r = abs(k_eigs[i].real - k_init_j.real)/max(abs(k_eigs[i].real), abs(k_init_j.real))
                    k_corr_i = abs(k_eigs[i].imag - k_init_j.imag)/max(abs(k_eigs[i].imag), abs(k_init_j.imag))
                    if self.vbose:
                        print('k_correlation_r: ', k_corr_r, '    k_correlation_i: ', k_corr_i)
                    if abs(k_init_j.real/k_init_j.imag)>1.0e4 or k_init_j.imag==0.0:
                        k_corr_i = 0.0
                        if self.vbose:
                            print('**** k_correlation_i set to zero ****')

                    if ind_corr[j]<0:
                        if abs(x_corr)>0.9 and k_corr_r<0.1 and k_corr_i<0.1:
                            ind_corr[j] = i
                            ind_x_corr[j] = x_corr
                    else:
                        if abs(x_corr)>0.9 and k_corr_r<0.1 and k_corr_i<0.1:
                            if x_corr > ind_x_corr[j]:
                                ind_corr[j] = i
                                ind_x_corr[j] = x_corr
            
            n_diverg_stat.append(n_diverge)
            if ind_corr[ind_lead]<0:
                for i in range(n_modes):
                    if mode_active[i]:
                        freq_sweep[i].refineLast()
                n_diverge += 1
                if n_diverge>3:
                    for i in range(n_modes):
                        if mode_active[i]:
                            freq_sweep[i].decreaseMaxStepSize()
                    n_diverge = 0
                if self.vbose:
                    print('making refinement...')
                    print('=-'*50)
                print('making refinement...')
                continue
                
            for i in range(n_modes):
                if ind_corr[i]<0:
                    mode_active[i] = False
            
            x_null = [0]*n_modes
            for i in range(n_modes):
                if mode_active[i]:
                    x_null[i] = k_vecs[:,ind_corr[i]]
                    x_null[i] /= np.linalg.norm(x_null[i])
                    x_null_prev[i] = x_null[i]

            #x_null_prev = x_null
            for i in range(n_modes):
                if mode_active[i]:
                    freq_sweep[i].setLastY([k_eigs[ind_corr[i]]])
            if self.vbose:
                for i in range(n_modes):
                    freq_sweep[i].printXY()
                print('--'*50)

            if record_data:
                x_null_sweep.append(x_null)
                
            if stop_cond!=None:
                cond__0 = "Re(k)<"
                if stop_cond.startswith(cond__0):
                    cond__1 = complex(stop_cond[len(cond__0):])
                    assert cond__1.imag == 0.0
                    for i in range(n_modes):
                        if mode_active[i]:
                            if k_eigs[ind_corr[i]].real < cond__1:
                                mode_active[i] = False
            #reset leader                            
            if mode_active[ind_lead]==False:
                for i in range(n_modes):
                    if mode_active[i]:
                        ind_lead = i
                        break
            if mode_active[ind_lead]==False:
                for i in range(n_modes):
                    if freq_sweep[i].Y[-1]==None:
                        ind_lead = i
                        freq_sweep[i].RemoveLast()
                        mode_active[i] = True
                        n_eig = n_modes     ## to avoid long convergence times
                        break
            #---
            if mode_active[ind_lead]==False:
                break

            if saveToFile!=None:
                f_name = saveToFile
                f = open(f_name, 'wb')
                params_dic = {'n_modes': n_modes, 'mode_active':mode_active, 
                    'ind_lead':ind_lead, 'n_eig':n_eig, 'x_null_prev': x_null_prev}
                pickle.dump(params_dic, f)
                pickle.dump(freq_sweep, f)
                f.close()

        print('divergence stats: ', n_diverg_stat)
        return [freq_sweep, x_null_sweep]
        
        
