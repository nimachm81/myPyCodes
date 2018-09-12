## FDTD 1D 2D 3D

__all__ = ["SimSpecToFDTDModel", "FDTDSimulator", "FDTDMPYEE", "FDTDYEE", 
            "FdtdCmds", "FVTypes"]

import numpy as np
import inspect
from multiprocessing import Process, Pipe, Queue
import multiprocessing, logging, sys
multiprocessing.log_to_stderr(logging.DEBUG)
#logger = multiprocessing.log_to_stderr()
#logger.setLevel(logging.DEBUG)
import time
import copy
import os

from sympy import Symbol, symbols, lambdify
from sympy.parsing.sympy_parser import parse_expr


from enum import Enum
class FdtdCmds(Enum):
    quit = 0
    stringcmd = 1
    getArrSlice = 2
    getArrSliceResp = 3

class FVTypes(Enum):
    ##fields
    E = 0
    H = 1
    D = 2
    B = 3
    ##sources
    JePoint = 10
    JmPoint = 11
    JeSheet = 12
    JmSheet = 13
    JeBox = 14
    JmBox = 15
    ##materials
    EpsIsoBox = 20
    EpsIsoSTvarBox = 21
    EpsDiagBox = 22
    EpsDiagSTvarBox = 23
    MuIsoBox = 24
    MuIsoSTvarBox = 25
    SigmaEDrudeIsoBox = 26
    SigmaEDrudeIsoSTvarBox = 27
    SigmaMDrudeIsoBox = 28
    SigmaMDrudeIsoSTvarBox = 29
    UPML = 30
    ##Boundary conditions
    BCond = 40
    PEC = 41
    PMC = 42
    PBC = 43
    ##view planes
    VPSide = 50
    
class SideFaceEnum(Enum):
    side = 0
    face = 1
    

TypeSideOrFace = {FVTypes.E:SideFaceEnum.side, 
                  FVTypes.H:SideFaceEnum.face,  
                  FVTypes.D:SideFaceEnum.side,
                  FVTypes.B:SideFaceEnum.face,
                  FVTypes.JePoint:SideFaceEnum.side,
                  FVTypes.JmPoint:SideFaceEnum.face,
                  FVTypes.JeSheet:SideFaceEnum.side,
                  FVTypes.JmSheet:SideFaceEnum.face,
                  FVTypes.JeBox:SideFaceEnum.side,
                  FVTypes.JmBox:SideFaceEnum.face,
                  FVTypes.EpsIsoBox:SideFaceEnum.side,
                  FVTypes.EpsIsoSTvarBox:SideFaceEnum.side,
                  FVTypes.EpsDiagBox:SideFaceEnum.side,
                  FVTypes.EpsDiagSTvarBox:SideFaceEnum.side,
                  FVTypes.MuIsoBox:SideFaceEnum.face,
                  FVTypes.MuIsoSTvarBox:SideFaceEnum.face,
                  FVTypes.SigmaEDrudeIsoBox:SideFaceEnum.side,
                  FVTypes.SigmaEDrudeIsoSTvarBox:SideFaceEnum.side,
                  FVTypes.SigmaMDrudeIsoBox:SideFaceEnum.face,
                  FVTypes.SigmaMDrudeIsoSTvarBox:SideFaceEnum.face }




class SimSpecToFDTDModel:
    """ takes the simulation specifications and parameters and runs an FDTD 
        simulation.
    """
    def __init__(self, simulationParameters, elements, globalvariables):
        self.simulationParameters = simulationParameters
        self.elements = elements
        self.globalvariables = globalvariables
        
        self.fdtd_sources = None
        self.fdtd_box = None
        self.fdtd_materials = None
        self.fdtd_ndims = 3
        self.fdtd_gridparams = None
        self.fdtd_outputs = None
        self.fdtd_boundaries = None
        return
        
    
    def validateSpecs(self):
        if "Global" in self.simulationParameters:
            globalParams = self.simulationParameters["Global"]
            if 'NumDimensions' in globalParams:
                n_dims = globalParams['NumDimensions']
                self.fdtd_ndims = n_dims
            else:
                print("Error: number of dimensions not set.")
                #return False
        else:
            print("Error: Global parameters not set.")
            #return False
        
        if "Simulation" in self.simulationParameters:
            grid_params = self.simulationParameters["Simulation"]
            dx_str = grid_params['dx']
                        
            dx,dy,dz = [None]*3
            dx = self.StrExpressionToFloat(dx_str)
            if self.fdtd_ndims>1:
                dy_str = grid_params['dy']
                dy = self.StrExpressionToFloat(dy_str)
            if self.fdtd_ndims>2:
                dz_str = grid_params['dz']
                dz = self.StrExpressionToFloat(dz_str)
            S_str = grid_params['S']
            S = self.StrExpressionToFloat(S_str)
            
            dr = None
            if self.fdtd_ndims==1:
                dr = np.array([dx])
            elif self.fdtd_ndims==2:
                dr = np.array([dx, dy])
            elif self.fdtd_ndims==3:
                dr = np.array([dx, dy, dz])
                
            
            N_save_str = grid_params['N_save']
            N_save = int(self.StrExpressionToFloat(N_save_str))
            
            
            self.fdtd_gridparams = {'type':'Yee', 'dr':dr, 'S':S, 'N_save':N_save}
            
            print('fdtd_gridparams : ', self.fdtd_gridparams)
            
        else:
            print("Error: FDTD simulation (grid) parameters not defind.")
        
        if "SimBox" in self.simulationParameters:
            box_name = self.simulationParameters["SimBox"]
            simbox = None
            print('box_name: ', box_name)
            
            assert "Objects" in self.simulationParameters
            objs = self.simulationParameters["Objects"]
            for obj in objs:
                if obj["name"] == box_name:
                    simbox = obj
                    break
            
            assert simbox!=None
            print(simbox)
            assert simbox['type']=='cube'
            
            r0 = [self.StrExpressionToFloat(r0_i) for r0_i in simbox['r0'] ]
            r1 = [self.StrExpressionToFloat(r1_i) for r1_i in simbox['r1'] ]
            
            self.fdtd_box = {'r0':r0, 'r1':r1}
            print('self.fdtd_box : ', self.fdtd_box)
        else:
            print("Error: Simulation box not set.")
            #return False


        if "Sources" in self.simulationParameters:
            fdtd_sources = []
            for i_s, source in enumerate(self.simulationParameters["Sources"]):
                if source['type'] == "electric delta" or source['type'] == "magnetic delta":
                    x = source['x']
                    y = source['y']
                    z = source['z']
                    r = np.array([x, y, z])
                    
                    Jx = source['Jx']
                    Jy = source['Jy']
                    Jz = source['Jz']
                    
                    Jx_t = self.validateTimeFunc(Jx)
                    Jy_t = self.validateTimeFunc(Jy)
                    Jz_t = self.validateTimeFunc(Jz)
                    
                    j_type = FVTypes.JePoint
                    if source['type'] == "magnetic delta":
                        j_type = FVTypes.JmPoint
                    
                    if Jx_t is None or Jy_t is None or Jz_t is None:
                        print("Error in source definitions.")
                        #return False
                        
                    if not (isinstance(Jx_t, float) and Jx_t==0.0):
                        print("Jx non-zero.")
                        je_args = {'r0':r, 'mag':1.0, 'f_t':Jx_t, 'src_dir':'x'}
                        
                        je_source = {'type':j_type, 'name':'Je-{}-x'.format(i_s), 'args':je_args}
                        
                        fdtd_sources.append(je_source)

                    if not (isinstance(Jy_t, float) and Jy_t==0.0):
                        print("Jy non-zero.")
                        je_args = {'r0':r, 'mag':1.0, 'f_t':Jx_t, 'src_dir':'y'}
                        
                        je_source = {'type':j_type, 'name':'Je-{}-y'.format(i_s), 'args':je_args}
                        
                        fdtd_sources.append(je_source)

                    if not (isinstance(Jz_t, float) and Jz_t==0.0):
                        print("Jz non-zero.")
                        je_args = {'r0':r, 'mag':1.0, 'f_t':Jx_t, 'src_dir':'z'}
                        
                        je_source = {'type':j_type, 'name':'Je-{}-z'.format(i_s), 'args':je_args}
                        
                        fdtd_sources.append(je_source)
            self.fdtd_sources = fdtd_sources
                    
        else:
            print("No sources defined.")
            #return False
            
        if "Boundaries" in self.simulationParameters:
            bound_conds = self.simulationParameters['Boundaries']
            
            bound_type = bound_conds["type"]
            
            fdtd_boundaries = {'type':bound_type}
            if bound_type=="PML":
                pml_thickness_str = bound_conds['pml_thickness']
                pml_stretch_fact_str = bound_conds['pml_stretch_fact']
                
                pml_thickness = self.StrExpressionToFloat(pml_thickness_str) 
                pml_stretch_fact = self.StrExpressionToComplex(pml_stretch_fact_str)
                
                fdtd_boundaries['pml_thickness'] = pml_thickness
                fdtd_boundaries['pml_stretch_fact'] = pml_stretch_fact
            
            self.fdtd_boundaries = fdtd_boundaries
            print('fdtd_boundaries : ', fdtd_boundaries)
            
        else:
            print("Error: Boundary conditions not set.")
            #return False
            
        if "Media" in self.simulationParameters:
            media = self.simulationParameters["Media"]
            materials_names = {}
            for i in range(len(media)):
                med_i = media[i]
                materials_names[med_i['name']] = i
                
            print("materials_names: ", materials_names)

            obj_materials = []
            fdtd_materials = []
            if "Objects" in self.simulationParameters:
                objs = self.simulationParameters["Objects"]
                for obj in objs:
                    if 'material' in obj and obj['material']!=None:
                        obj_mat_name = obj['material']
                        
                        assert obj_mat_name in materials_names
                        print('name: {},  material:{} '.format(obj['name'], obj_mat_name))
                        
                        mat_params = media[materials_names[obj_mat_name]]
                        print('medium parameters: ', mat_params)
                        
                        if mat_params['type']=='diagonal':
                            eps_r_xx_str = mat_params['eps_r_xx']
                            st_vars_eps_xx, eps_xx_func = self.getSpaceTimeDependanceAndFunc( eps_r_xx_str, ret='rt' )
                            eps_r_yy_str = mat_params['eps_r_yy']
                            st_vars_eps_yy, eps_yy_func = self.getSpaceTimeDependanceAndFunc( eps_r_yy_str, ret='rt' )
                            eps_r_zz_str = mat_params['eps_r_zz']
                            st_vars_eps_zz, eps_zz_func = self.getSpaceTimeDependanceAndFunc( eps_r_zz_str, ret='rt' )
                            
                            mu_r_xx_str = mat_params['mu_r_xx']
                            st_vars_mu_xx, mu_xx_func = self.getSpaceTimeDependanceAndFunc( mu_r_xx_str, ret='rt' )
                            mu_r_yy_str = mat_params['mu_r_yy']
                            st_vars_mu_yy, mu_yy_func = self.getSpaceTimeDependanceAndFunc( mu_r_yy_str, ret='rt' )
                            mu_r_zz_str = mat_params['mu_r_zz']
                            st_vars_mu_zz, mu_zz_func = self.getSpaceTimeDependanceAndFunc( mu_r_zz_str, ret='rt' )
                                                        
                            if obj['type']=='cube':
                                r0_str = obj['r0']
                                r1_str = obj['r1']
                                
                                r0 = [self.StrExpressionToFloat( r0_str[i] ) for i in range(3)]
                                r1 = [self.StrExpressionToFloat( r1_str[i] ) for i in range(3)]
                                
                                eps_r_diag = [eps_xx_func, eps_yy_func, eps_zz_func]
                                eps_args = {'r0':r0, 'r1':r1, 'mag_in':eps_r_diag, 'mag_out':0.0}
                                fdtd_mat = {'type':FVTypes.EpsDiagSTvarBox, 'name':'eps_{}'.format(len(fdtd_materials)), 'args':eps_args}
                                fdtd_materials.append(fdtd_mat)
                            else:
                                raise NotImplementedError()
                        else:
                            raise NotImplementedError()

            self.fdtd_materials = fdtd_materials
            print('fdtd_materials : ', fdtd_materials)

        else:
            print("Error: No medium defined.")
            
        if "Outputs" in self.simulationParameters:
            output_planes = self.simulationParameters["Outputs"]
            n_dims = self.fdtd_ndims
            
            fdtd_outputs = []
            for out in output_planes:
                out_type = out["type"]
                out_field = out["field"]
                out_r_str = out["r"]
                out_name = out["name"]
                out_plane = out["plane"]
                
                out_r = [self.StrExpressionToFloat(out_r_str[i]) for i in range(3)]
                
                vp_type = None
                if out_field[0]=='E':
                    vp_type = FVTypes.VPSide
                elif out_field[0]=='H':
                    raise NotImplementedError()
                    
                
                out_args = None
                if out_type=="entire":
                    out_args = {'type':vp_type, 'r':out_r, 'args':{'A':out_field[0], 'A_dir':out_field[1], 'O_dir':None, 'name':out_name}}
                elif out_type=="cut":
                    out_args = {'type':vp_type, 'r':out_r, 'args':{'A':out_field[0], 'A_dir':out_field[1], 'O_dir':out_plane, 'name':out_name}}
                
                fdtd_outputs.append(out_args)
            self.fdtd_outputs = fdtd_outputs
            print('fdtd_outputs: ', self.fdtd_outputs)
            
        else:
            print("Error: No output plane defined.")
            #return False
             
        
        return True


    def getSpaceTimeDependanceAndFunc(self, f_str, ret='txyz'):
        """ ret = 'txyz' or 'rt'
            ret = 'txyz' --->   returns f(t, x, y, z)
            ret = 'rt'   --->   returns f(r, t)
        """
        f_sym = parse_expr(f_str)
        for var in self.globalvariables:
            f_sym = f_sym.subs(Symbol(var), self.globalvariables[var])
        
        #print(f_sym)

        t, x, y, z = symbols('t x y z')
        f_vars = f_sym.free_symbols
        st_vars = []
        if t in f_vars:
            st_vars.append(t)
        if x in f_vars:
            st_vars.append(x)
        if y in f_vars:
            st_vars.append(y)
        if z in f_vars:
            st_vars.append(z)
        
        f_func = None
        if len(st_vars)==0:
            f_func = float(f_sym)
        elif len(st_vars)==1:
            f_func = lambdify(*tuple(st_vars), f_sym, modules="numpy")
        else:
            f_func = lambdify(st_vars, f_sym, modules="numpy")
            
        if ret=='txyz':
            return st_vars, f_func
        elif ret=='rt':
            f_func_std = None
            if st_vars==[t]:
                f_func_std = lambda r,t : f_func(t)
            elif st_vars==[x]:
                f_func_std = lambda r,t : f_func(r[0])
            elif st_vars==[y]:
                f_func_std = lambda r,t : f_func(r[1])
            elif st_vars==[z]:
                f_func_std = lambda r,t : f_func(r[2])
            elif st_vars==[t, x]:
                f_func_std = lambda r,t : f_func(t, r[0])
            elif st_vars==[t, y]:
                f_func_std = lambda r,t : f_func(t, r[1])
            elif st_vars==[t, z]:
                f_func_std = lambda r,t : f_func(t, r[2])
            elif st_vars==[x, y]:
                f_func_std = lambda r,t : f_func(r[0], r[1])
            elif st_vars==[x, z]:
                f_func_std = lambda r,t : f_func(r[0], r[2])
            elif st_vars==[y, z]:
                f_func_std = lambda r,t : f_func(r[1], r[2])
            elif st_vars==[t, x, y]:
                f_func_std = lambda r,t : f_func(t, r[0], r[1])
            elif st_vars==[t, x, z]:
                f_func_std = lambda r,t : f_func(t, r[0], r[2])
            elif st_vars==[t, y, z]:
                f_func_std = lambda r,t : f_func(t, r[1], r[2])
            elif st_vars==[t, x, y, z]:
                f_func_std = lambda r,t : f_func(t, r[0], r[1], r[2])
            else:
                assert len(st_vars)==0
                f_func_std = f_func


            return st_vars, f_func_std

            
                                        


    def validateTimeFunc(self, f_str):
        f_sym = parse_expr(f_str)
        for var in self.globalvariables:
            f_sym = f_sym.subs(Symbol(var), self.globalvariables[var])
        
        print(f_sym)

        t = Symbol('t')
        if t not in f_sym.free_symbols:
            return float(f_sym)
        else:
            f_t = lambdify(t, f_sym, modules="numpy")
            if isinstance(f_t(1.0), float):
                return f_t
                
            else:
                return None
        

    def StrExpressionToNumpyFunc(self, expr, vars=Symbol('t')):
        eq = parse_expr(expr)
        for var in self.globalvariables:
            eq = eq.subs(Symbol(var), self.globalvariables[var])
            
        np_func = lambdify(vars, eq, modules='numpy')  
        return np_func
        

    def StrExpressionToFloat(self, expr):
        eq = parse_expr(expr)
        for var in self.globalvariables:
            eq = eq.subs(Symbol(var), self.globalvariables[var])
        
        return float(eq.evalf())
        

    def StrExpressionToComplex(self, expr):
        eq = parse_expr(expr)
        for var in self.globalvariables:
            eq = eq.subs(Symbol(var), self.globalvariables[var])
            
        return complex(eq.evalf())



class FDTDSimulator:
    """ sets up FDTD parameters and runs the simulation
    """
    def __init__(self, pr_grid=None, dtype=float, vbose=False):
        """ pr_grid: None ---> runs single process
                    tuple ---> runs MP version
        """
        self.pr_grid = pr_grid
        self.materials = []
        self.sources = []
        self.pmls = []
        self.viewPlans = []
        self.save_every = 50
        self.bcs = None
        self.dtype = dtype
        self.vbose = vbose
        self.k_vec = None

    def SetSimulationBox(self, r0, r1, dr, dt):
        self.r0 = r0
        self.r1 = r1
        self.dr = dr
        self.dt = dt

    def AddMaterial(self, mat_params):
        self.materials.append(mat_params)

    def AddSources(self, source_params):
        self.sources.append(source_params)

    def AddPML(self, pml_params):
        self.pmls.append(pml_params)
        
    def AddViewPlane(self, vp_params):
        self.viewPlans.append(vp_params)
        
    def SetSaveEvery(self, save_every):
        self.save_every = save_every
        
    def SetBCs(self, bcs):
        self.bcs = bcs
        
    def StepFields(self, n_t):
        if self.pr_grid is None:
            fdtd = FDTDYEE(self.r0, self.r1, self.dr, self.dt, dtype=self.dtype, vbose=self.vbose)
            fdtd.AllocateSideArr_list(["E", "D"])
            fdtd.AllocateFaceArr_list(["B", "H"])            

            fdtd.SetSpatialPointsSide()
            fdtd.SetSpatialPointsFace()

            fdtd.SetupFdtdVars(FVTypes.E, None, var=fdtd.E)
            fdtd.SetupFdtdVars(FVTypes.H, None, var=fdtd.H)
            fdtd.SetupFdtdVars(FVTypes.D, None, var=fdtd.D)
            fdtd.SetupFdtdVars(FVTypes.B, None, var=fdtd.B)
            
            fdtd.save_every = self.save_every 

            for src_params in self.sources:
                src_type = src_params['type']
                src_name = src_params['name']
                src_args = src_params['args']
                if TypeSideOrFace[src_type] == SideFaceEnum.side:
                    fdtd.AllocateSideArr_list([src_name])
                elif TypeSideOrFace[src_type] == SideFaceEnum.face:
                    fdtd.AllocateFaceArr_list([src_name])
                else:
                    raise NotImplementedError()
                src_var = [None]
                exec('src_var[0] = fdtd.{}'.format(src_name))
                src_var = src_var[0]
                fdtd.SetupFdtdVars(src_type, src_args, var=src_var)
                
                
            for mat_params in self.materials:
                mat_type = mat_params['type']
                mat_name = mat_params['name']
                mat_args = mat_params['args']
                if TypeSideOrFace[mat_type] == SideFaceEnum.side:
                    fdtd.AllocateSideArr_list([mat_name])
                elif TypeSideOrFace[mat_type] == SideFaceEnum.face:
                    fdtd.AllocateFaceArr_list([mat_name])
                else:
                    raise NotImplementedError()
                mat_var = [None]
                exec('mat_var[0] = fdtd.{}'.format(mat_name))
                mat_var = mat_var[0]
                fdtd.SetupFdtdVars(mat_type, mat_args, var=mat_var)
            
            assert len(self.pmls)<=1
            for pml_params in self.pmls:
                pml_type = pml_params['type']
                #pml_name = pml_params['name']
                pml_args = pml_params['args']
                fdtd.SetupFdtdVars(pml_type, pml_args, var=None)
                

            if self.bcs!=None:
                bc_args = {'bcs':self.bcs}
                fdtd.SetupFdtdVars(FVTypes.BCond, bc_args, var=None)

            for v in self.viewPlans:
                #Example:   v = {'type':'side', 'r':r_0, 
                #               'args':{'A':'E', 'A_dir':'z', 'O_dir':'x', 'name':'E'}}
                # 'A' --> variable name
                # 'name': vp name
                vp_type = v['type']
                vp_r = v['r']
                vp_args = v['args'] 
                vp_args_var = [None]
                exec("vp_args_var[0] = fdtd.{}".format(vp_args['A']))
                
                vp_args_cp = copy.deepcopy(vp_args)
                vp_args_cp['A'] = vp_args_var[0]
                
                if vp_type==FVTypes.VPSide:
                    fdtd.SetViewPlane_Side(vp_r, [vp_args_cp])
                else:
                    raise NotImplementedError()
            
            if self.k_vec!=None:
                fdtd.SetKVec(self.k_vec)
            tic = time.perf_counter()
            self.fdtd = fdtd
            self.fdtd.StepFields(n_t)

            toc = time.perf_counter()
            print('simulation time: {}:{}'.format(int((toc-tic)/60), int((toc-tic)%60)))
            self.n_saved = self.fdtd.n_saved
            return

        else:
            raise NotImplementedError()
        
    def GetOutputs(self, name, file_ind=0):
        return self.fdtd.GetOutputs(name, file_ind)        
        
    def SetKVec(self, k_vec):
        self.k_vec = k_vec
    


##TODO: sources may be counted multiple times if they are located on the boundary
class FDTDMPYEE:
    def __init__(self, r0, r1, dr, dt, pr_grid):
        """ pr_grid: process grid  Ex: (2, 3, 4) in 3D --> 24 processes
        """
        assert len(r0)==len(r1)==len(dr)==len(pr_grid)
        self.n_dim=len(r0)
        self.pr_grid = pr_grid
        self.processGrid = np.empty(pr_grid, dtype=object)

        self.r0 = np.copy(r0)
        self.r1 = np.copy(r1)
        self.W = self.r1-self.r0
        self.dt = dt
        
        self.vbose = False
                
        self.tag = 0

        dr_chunk = (r1-r0)/pr_grid
        self.N_chunk = None
        self.dr = None
        if self.n_dim==3:
            for i in range(pr_grid[0]):
                for j in range(pr_grid[1]):
                    for k in range(pr_grid[2]):
                        r0_ch = r0+np.array([i, j, k])*dr_chunk
                        r1_ch = r0+np.array([i+1, j+1, k+1])*dr_chunk #to make corners exactly match
                        #r1_ch = r0_ch + dr_chunk
                        if self.N_chunk is None:
                            self.processGrid[i,j,k] =  FDTDYEE(r0_ch, r1_ch, dr, dt)
                            self.N_chunk = np.copy(self.processGrid[i,j,k].N)
                            self.dr = np.copy(self.processGrid[i,j,k].dr)
                        else:
                            self.processGrid[i,j,k] =  FDTDYEE(r0_ch, r1_ch, dr, dt, force_N=self.N_chunk)
                            assert np.all(self.processGrid[i,j,k].N==self.N_chunk)
                            assert np.all(self.processGrid[i,j,k].dr==self.dr)
                        self.processGrid[i,j,k].SetCartesianIndex(ChunkInd=np.array([i, j, k]), ChunkIndTot=np.copy(pr_grid))
                        self.processGrid[i,j,k].SetCartesianIndexNeighbors()
                            
        elif self.n_dim==2:
            for i in range(pr_grid[0]):
                for j in range(pr_grid[1]):
                    r0_ch = r0+np.array([i, j])*dr_chunk
                    r1_ch = r0+np.array([i+1, j+1])*dr_chunk #to make corners exactly match
                    #r1_ch = r0_ch + dr_chunk
                    if self.N_chunk is None:
                        self.processGrid[i,j] =  FDTDYEE(r0_ch, r1_ch, dr, dt)
                        self.N_chunk = np.copy(self.processGrid[i,j].N)
                        self.dr = np.copy(self.processGrid[i,j].dr)
                    else:
                        self.processGrid[i,j] =  FDTDYEE(r0_ch, r1_ch, dr, dt, force_N=self.N_chunk)
                        assert np.all(self.processGrid[i,j].N==self.N_chunk)
                        assert np.all(self.processGrid[i,j].dr==self.dr)
                    self.processGrid[i,j].SetCartesianIndex(ChunkInd=np.array([i, j]), ChunkIndTot=np.copy(pr_grid))
                    self.processGrid[i,j].SetCartesianIndexNeighbors()
                            
        elif self.n_dim==1:
            for i in range(pr_grid[0]):
                r0_ch = r0+np.array([i])*dr_chunk
                r1_ch = r0+np.array([i+1])*dr_chunk #to make corners exactly match
                #r1_ch = r0_ch + dr_chunk
                if self.N_chunk is None:
                    self.processGrid[i] =  FDTDYEE(r0_ch, r1_ch, dr, dt)
                    self.N_chunk = np.copy(self.processGrid[i].N)
                    self.dr = np.copy(self.processGrid[i].dr)
                else:
                    self.processGrid[i] =  FDTDYEE(r0_ch, r1_ch, dr, dt, force_N=self.N_chunk)
                    assert np.all(self.processGrid[i].N==self.N_chunk)
                    assert np.all(self.processGrid[i].dr==self.dr)
                self.processGrid[i].SetCartesianIndex(ChunkInd=np.array([i]), ChunkIndTot=np.copy(pr_grid))
                self.processGrid[i].SetCartesianIndexNeighbors()

        print('self.pr_grid:', self.pr_grid, 'self.N_chunk:', self.N_chunk)
        self.SetInterCommunicators()
        self.SetCommandPipes()
        self.StartProcesses()

        self.logProcContents()
        return


    def getProcessPMLs(self, d_pml_n, s_pml_n, d_pml_p=None, s_pml_p=None):
        if d_pml_p is None:
            d_pml_p = d_pml_n
        if s_pml_p is None:
            s_pml_p = s_pml_n
        pr_grid = self.pr_grid
        self.pml_grid = np.empty(pr_grid, dtype=object)
        r0_nopml = self.r0 + d_pml_n
        r1_nopml = self.r1 - d_pml_p
        if self.n_dim==3:
            for i in range(pr_grid[0]):
                for j in range(pr_grid[1]):
                    for k in range(pr_grid[2]):
                        r0_ch = self.processGrid[i,j,k].r0
                        r1_ch = self.processGrid[i,j,k].r1
                        d_pml_n_ch = (r0_nopml-r0_ch)*(r0_nopml>r0_ch)
                        d_pml_p_ch = (r1_ch-r1_nopml)*(r1_nopml<r1_ch)
                        self.pml_grid[i,j,k] = [d_pml_n_ch, s_pml_n, d_pml_p_ch, s_pml_p]
                        
        elif self.n_dim==2:
            for i in range(pr_grid[0]):
                for j in range(pr_grid[1]):
                    r0_ch = self.processGrid[i,j].r0
                    r1_ch = self.processGrid[i,j].r1
                    d_pml_n_ch = (r0_nopml-r0_ch)*(r0_nopml>r0_ch)
                    d_pml_p_ch = (r1_ch-r1_nopml)*(r1_nopml<r1_ch)
                    self.pml_grid[i,j] = [d_pml_n_ch, s_pml_n, d_pml_p_ch, s_pml_p]
                        
        elif self.n_dim==1:
            for i in range(pr_grid[0]):
                r0_ch = self.processGrid[i].r0
                r1_ch = self.processGrid[i].r1
                d_pml_n_ch = (r0_nopml-r0_ch)*(r0_nopml>r0_ch)
                d_pml_p_ch = (r1_ch-r1_nopml)*(r1_nopml<r1_ch)
                self.pml_grid[i] = [d_pml_n_ch, s_pml_n, d_pml_p_ch, s_pml_p]
                        
    
    ##TODO: set bcs for each process
    def setProcessBCs(self, bcs):
        raise NotImplementedError()
        return
        

    def logProcContents(self):
        pr_grid = self.pr_grid
        if self.n_dim==3:
            for i in range(pr_grid[0]):
                for j in range(pr_grid[1]):
                    for k in range(pr_grid[2]):
                        self.processGrid[i,j,k].LogInfo()
        if self.n_dim==2:
            for i in range(pr_grid[0]):
                for j in range(pr_grid[1]):
                        self.processGrid[i,j].LogInfo()
        if self.n_dim==1:
            for i in range(pr_grid[0]):
                self.processGrid[i].LogInfo()


                            
    def SetCommandPipes(self):
        pr_grid = self.pr_grid
        self.cmdGrid = np.empty(pr_grid, dtype=object)
        if self.n_dim==3:
            for i in range(pr_grid[0]):
                for j in range(pr_grid[1]):
                    for k in range(pr_grid[2]):
                        cmd, cmd_res = Pipe()
                        self.cmdGrid[i,j,k] = cmd
                        self.processGrid[i,j,k].SetCommandPipe(cmd_res)
        
        elif self.n_dim==2:
            for i in range(pr_grid[0]):
                for j in range(pr_grid[1]):
                    cmd, cmd_res = Pipe()
                    self.cmdGrid[i,j] = cmd
                    self.processGrid[i,j].SetCommandPipe(cmd_res)

        elif self.n_dim==1:
            for i in range(pr_grid[0]):
                cmd, cmd_res = Pipe()
                self.cmdGrid[i] = cmd
                self.processGrid[i].SetCommandPipe(cmd_res)
                                

    def SetInterCommunicators(self):
        pr_grid = self.pr_grid
        if self.n_dim==3:
            for i in range(pr_grid[0]):
                for j in range(pr_grid[1]):
                    for k in range(pr_grid[2]):
                        ChunkIndNB = self.processGrid[i,j,k].ChunkIndNB
                        for n_dir in range(3):
                            if ChunkIndNB['p'][n_dir]!=None:
                                comm, comm_rec = Pipe()
                                grid_ind_nb = ChunkIndNB['p'][n_dir]
                                self.processGrid[i,j,k].AddCommPipesNBSend(grid_ind_nb, comm)
                                self.processGrid[tuple(grid_ind_nb)].AddCommPipesNBRecv((i,j,k), comm_rec)
                            if ChunkIndNB['n'][n_dir]!=None:
                                comm, comm_rec = Pipe()
                                grid_ind_nb = ChunkIndNB['n'][n_dir]
                                self.processGrid[i,j,k].AddCommPipesNBSend(grid_ind_nb, comm)
                                self.processGrid[tuple(grid_ind_nb)].AddCommPipesNBRecv((i,j,k), comm_rec)
        
        elif self.n_dim==2:
            for i in range(pr_grid[0]):
                for j in range(pr_grid[1]):
                    ChunkIndNB = self.processGrid[i,j].ChunkIndNB
                    for n_dir in range(2):
                        if ChunkIndNB['p'][n_dir]!=None:
                            comm, comm_rec = Pipe()
                            grid_ind_nb = ChunkIndNB['p'][n_dir]
                            self.processGrid[i,j].AddCommPipesNBSend(grid_ind_nb, comm)
                            self.processGrid[tuple(grid_ind_nb)].AddCommPipesNBRecv((i,j), comm_rec)
                        if ChunkIndNB['n'][n_dir]!=None:
                            comm, comm_rec = Pipe()
                            grid_ind_nb = ChunkIndNB['n'][n_dir]
                            self.processGrid[i,j].AddCommPipesNBSend(grid_ind_nb, comm)
                            self.processGrid[tuple(grid_ind_nb)].AddCommPipesNBRecv((i,j), comm_rec)

        elif self.n_dim==1:
            for i in range(pr_grid[0]):
                ChunkIndNB = self.processGrid[i].ChunkIndNB
                for n_dir in range(1):
                    if ChunkIndNB['p'][n_dir]!=None:
                        comm, comm_rec = Pipe()
                        grid_ind_nb = ChunkIndNB['p'][n_dir]
                        self.processGrid[i].AddCommPipesNBSend(grid_ind_nb, comm)
                        self.processGrid[tuple(grid_ind_nb)].AddCommPipesNBRecv((i,), comm_rec)
                    if ChunkIndNB['n'][n_dir]!=None:
                        comm, comm_rec = Pipe()
                        grid_ind_nb = ChunkIndNB['n'][n_dir]
                        self.processGrid[i].AddCommPipesNBSend(grid_ind_nb, comm)
                        self.processGrid[tuple(grid_ind_nb)].AddCommPipesNBRecv((i,), comm_rec)


    def StartProcesses(self):
        pr_grid = self.pr_grid
        if self.n_dim==3:
            for i in range(pr_grid[0]):
                for j in range(pr_grid[1]):
                    for k in range(pr_grid[2]):
                        self.processGrid[i,j,k].start()
        elif self.n_dim==2:
            for i in range(pr_grid[0]):
                for j in range(pr_grid[1]):
                    self.processGrid[i,j].start()
        elif self.n_dim==1:
            for i in range(pr_grid[0]):
                self.processGrid[i].start()
        
    def TerminateProcesses(self):
        pr_grid = self.pr_grid
        if self.n_dim==3:
            for i in range(pr_grid[0]):
                for j in range(pr_grid[1]):
                    for k in range(pr_grid[2]):
                        #self.processGrid[i,j,k].terminate()
                        #assert self.processGrid[i, j, k].quitLoop==True
                        self.processGrid[i,j,k].join()
        elif self.n_dim==2:
            for i in range(pr_grid[0]):
                for j in range(pr_grid[1]):
                    #self.processGrid[i,j].terminate()
                    #assert self.processGrid[i, j].quitLoop==True
                    self.processGrid[i,j].join()
        elif self.n_dim==1:
            for i in range(pr_grid[0]):
                #self.processGrid[i].terminate()
                #assert self.processGrid[i].quitLoop==True
                self.processGrid[i].join()
        
              
    def getUniqueTag(self):
        self.tag += 1
        return self.tag

    def SendQuitMessage(self):
        pr_grid = self.pr_grid
        tag = self.getUniqueTag()
        if self.n_dim==3:
            for i in range(pr_grid[0]):
                for j in range(pr_grid[1]):
                    for k in range(pr_grid[2]):
                        self.cmdGrid[i,j,k].send([FdtdCmds.quit, tag, None, None])
            for i in range(pr_grid[0]):
                for j in range(pr_grid[1]):
                    for k in range(pr_grid[2]):
                        resp = self.cmdGrid[i,j,k].recv()
                        assert resp[2][0]=='ok'
        elif self.n_dim==2:
            for i in range(pr_grid[0]):
                for j in range(pr_grid[1]):
                    self.cmdGrid[i,j].send([FdtdCmds.quit, None, None, None])
            for i in range(pr_grid[0]):
                for j in range(pr_grid[1]):
                    resp = self.cmdGrid[i,j].recv()
                    assert resp[2][0]=='ok'
        elif self.n_dim==1:
            for i in range(pr_grid[0]):
                self.cmdGrid[i].send([FdtdCmds.quit, None, None, None])
            for i in range(pr_grid[0]):
                resp = self.cmdGrid[i].recv()
                assert resp[2][0]=='ok'
    
    def SendStrMessage(self, msg_str, msg_ret=None):
        pr_grid = self.pr_grid
        tag = self.getUniqueTag()
        return_dic = {}
        if self.n_dim==3:
            for i in range(pr_grid[0]):
                for j in range(pr_grid[1]):
                    for k in range(pr_grid[2]):
                        self.cmdGrid[i,j,k].send([FdtdCmds.stringcmd, tag, msg_str, msg_ret])
            for i in range(pr_grid[0]):
                for j in range(pr_grid[1]):
                    for k in range(pr_grid[2]):
                        resp = self.cmdGrid[i,j,k].recv()
                        assert resp[1]==tag and resp[2][0]=='ok'
                        return_dic[(i,j,k)] = resp[2][1]
        elif self.n_dim==2:
            for i in range(pr_grid[0]):
                for j in range(pr_grid[1]):
                    self.cmdGrid[i,j].send([FdtdCmds.stringcmd, tag, msg_str, msg_ret])
            for i in range(pr_grid[0]):
                for j in range(pr_grid[1]):
                    resp = self.cmdGrid[i,j].recv()
                    assert resp[1]==tag and resp[2][0]=='ok'
                    return_dic[(i,j)] = resp[2][1]
        elif self.n_dim==1:
            for i in range(pr_grid[0]):
                self.cmdGrid[i].send([FdtdCmds.stringcmd, tag, msg_str, msg_ret])
            for i in range(pr_grid[0]):
                resp = self.cmdGrid[i].recv()
                assert resp[1]==tag and resp[2][0]=='ok'
                return_dic[(i,)] = resp[2][1]
        return return_dic


    def GetOutputs(self, name):
        ##TODO: not working.. perhaps scope issue
        pr_grid = self.pr_grid
        outs = np.empty(pr_grid, dtype=object)
        if self.n_dim==3:
            for i in range(pr_grid[0]):
                for j in range(pr_grid[1]):
                    for k in range(pr_grid[2]):
                        outs[i,j,k] = self.processGrid[i,j,k].GetOutputs(name)
        if self.n_dim==2:
            for i in range(pr_grid[0]):
                for j in range(pr_grid[1]):
                    outs[i,j] = self.processGrid[i,j].GetOutputs(name)
    
        if self.n_dim==1:
            for i in range(pr_grid[0]):
                outs[i] = self.processGrid[i].GetOutputs(name)
        return outs



class FDTDYEE(Process):
    
    def __init__(self, r0, r1, dr, dt, N_is_even=None, force_N=None, dtype=float, vbose=False):
        super(FDTDYEE, self).__init__()
        self.daemon = True

        assert len(r0)==len(r1)==len(dr)
        self.n_dim=len(r0)
        
        self.vbose = vbose
    
        self.r0 = np.copy(r0)   #3 vectors
        self.r1 = np.copy(r1)
        self.W = self.r1-self.r0
        if force_N is None:
            self.N = (self.W/dr).astype(int)
            if N_is_even!=None:
                assert len(N_is_even)==self.n_dim
                for i in range(self.n_dim):
                    if N_is_even[i]==True:
                        if self.N[i]%2==1:
                            self.N[i] += 1
                    else:
                        if self.N[i]%2==0:
                            self.N[i] += 1
        else:
            self.N = np.copy(force_N)
                
        self.dr = self.W/self.N
        self.dt = dt

        self.save_every = 10
        self.n_saved = 0

        self.dtype = dtype
        
        if self.n_dim==3:
            self.Nsx = self.N+np.array([0, 1, 1])
            self.Nsy = self.N+np.array([1, 0, 1])
            self.Nsz = self.N+np.array([1, 1, 0])

            self.Nfx = self.N+np.array([1, 0, 0])
            self.Nfy = self.N+np.array([0, 1, 0])
            self.Nfz = self.N+np.array([0, 0, 1])
        elif self.n_dim==2:
            self.Nsx = self.N+np.array([0, 1])
            self.Nsy = self.N+np.array([1, 0])
            self.Nsz = self.N+np.array([1, 1])

            self.Nfx = self.N+np.array([1, 0])
            self.Nfy = self.N+np.array([0, 1])
            self.Nfz = self.N+np.array([0, 0])
        elif self.n_dim==1:
            self.Nsx = self.N+np.array([0])
            self.Nsy = self.N+np.array([1])
            self.Nsz = self.N+np.array([1])

            self.Nfx = self.N+np.array([1])
            self.Nfy = self.N+np.array([0])
            self.Nfz = self.N+np.array([0])
        
        
        self.ViewPlanes = {}

        self.ChunkInd = None
        self.ChunkIndTot = None
        
        self.tag = 0
        
        self.cmdPipe = None
        self.commPipesNBDicSend = {}
        self.commPipesNBDicRecv = {}
        self.quitLoop = False
        
        self.Sources = []    ## {'type':xx, 'args':{}, 'var':xx}
        self.Materials = []
        self.Fields = []
        
        self.bcs = None
        self.k_vec = None
        return

    def SetKVec(self, k_vec):
        self.k_vec = k_vec

        
    def SetupFdtdVars(self, vtype, args, var):
        if vtype==FVTypes.E or vtype==FVTypes.H or vtype==FVTypes.D or vtype==FVTypes.B:
            self.Fields.append({'type':vtype, 'args':args, 'var':var})
        elif vtype==FVTypes.JePoint:
            r0 = args['r0']
            mag = args['mag']
            f_t = args['f_t']
            src_dir = args['src_dir']
            
            f_r = self.GetPointSourceFunc(r0, mag, src_dir, em_type='e')
            J_0 = self.UpdateFuncSide_space(var, f_r, getCopy=True)
            
            jargs = {'r0':r0, 'mag':mag, 'f_t':f_t, 'src_dir':src_dir, 'f_r':f_r, 'V_0':J_0}
            
            jdic = {'type':vtype, 'args':jargs, 'var':var}
            self.Sources.append(jdic)
        elif vtype==FVTypes.JmPoint:
            r0 = args['r0']
            mag = args['mag']
            f_t = args['f_t']
            src_dir = args['src_dir']
            
            f_r = self.GetPointSourceFunc(r0, mag, src_dir, em_type='m')
            J_0 = self.UpdateFuncFace_space(var, f_r, getCopy=True)
            
            jargs = {'r0':r0, 'mag':mag, 'f_t':f_t, 'src_dir':src_dir, 'f_r':f_r, 'V_0':J_0}
            
            jdic = {'type':vtype, 'args':jargs, 'var':var}
            self.Sources.append(jdic)
        elif vtype==FVTypes.JeSheet:
            r0 = args['r0']
            mag = args['mag']
            f_t = args['f_t']
            src_dir = args['src_dir']
            norm_dir = args['norm_dir']
            
            f_r = self.GetSheetSourceFunc(r0, mag, src_dir, em_type='e')
            J_0 = self.UpdateFuncSide_space(var, f_r, getCopy=True)
            
            jargs = {'r0':r0, 'mag':mag, 'f_t':f_t, 'src_dir':src_dir, \
                     'norm_dir':norm_dir, 'f_r':f_r, 'V_0':J_0}
            
            jdic = {'type':vtype, 'args':jargs, 'var':var}
            self.Sources.append(jdic)
        elif vtype==FVTypes.JmSheet:
            r0 = args['r0']
            mag = args['mag']
            f_t = args['f_t']
            src_dir = args['src_dir']
            norm_dir = args['norm_dir']
            
            f_r = self.GetSheetSourceFunc(r0, mag, src_dir, em_type='m')
            J_0 = self.UpdateFuncFace_space(var, f_r, getCopy=True)
            
            jargs = {'r0':r0, 'mag':mag, 'f_t':f_t, 'src_dir':src_dir, \
                     'norm_dir':norm_dir, 'f_r':f_r, 'V_0':J_0}
            
            jdic = {'type':vtype, 'args':jargs, 'var':var}
            self.Sources.append(jdic)
        elif vtype==FVTypes.JeBox:
            r0 = args['r0']
            r1 = args['r1']
            mag = args['mag']
            f_t = args['f_t']
            src_dir = args['src_dir']
            
            f_r = self.GetBoxSourceFunc(r0, r1, mag, src_dir, em_type='e')
            J_0 = self.UpdateFuncSide_space(var, f_r, getCopy=True)
            
            jargs = {'r0':r0, 'r1':r1, 'mag':mag, 'f_t':f_t, 'src_dir':src_dir, \
                     'f_r':f_r, 'V_0':J_0}
            
            jdic = {'type':vtype, 'args':jargs, 'var':var}
            self.Sources.append(jdic)
        elif vtype==FVTypes.JmBox:
            r0 = args['r0']
            r1 = args['r1']
            mag = args['mag']
            f_t = args['f_t']
            src_dir = args['src_dir']
            
            f_r = self.GetBoxSourceFunc(r0, r1, mag, src_dir, em_type='e')
            J_0 = self.UpdateFuncFace_space(var, f_r, getCopy=True)
            
            jargs = {'r0':r0, 'r1':r1, 'mag':mag, 'f_t':f_t, 'src_dir':src_dir, \
                     'f_r':f_r, 'V_0':J_0}
            
            jdic = {'type':vtype, 'args':jargs, 'var':var}
            self.Sources.append(jdic)
        elif vtype==FVTypes.EpsIsoBox:
            ##Warning: if the boundaries of 2 contiguent boxes are located on a side element
            ## one of the boxes should exclude the side
            r0 = args['r0']
            r1 = args['r1']
            a_eps = args['mag_in']
            b_eps = args['mag_out']
            exclude = False
            if 'exclude' in args:
                exclude = args['exclude']
            
            f_eps = self.FunctionBox(r0, r1, a_eps, b_eps, exclude)
            f_eps_iso = [f_eps]*3            
            self.UpdateFuncSide_space(var, f_eps_iso)
            
            epsargs = {'r0':r0, 'r1':r1, 'mag_in':a_eps, 'mag_out':b_eps, 'f_eps':f_eps_iso}
            
            epsdic = {'type':vtype, 'args':epsargs, 'var':var}
            self.Materials.append(epsdic)
        elif vtype==FVTypes.EpsIsoSTvarBox:
            r0 = args['r0']
            r1 = args['r1']
            a_eps = args['mag_in']
            b_eps = args['mag_out']
            exclude = False
            if 'exclude' in args:
                exclude = args['exclude']
            
            f_eps = self.FunctionBox(r0, r1, a_eps, b_eps, exclude)
            f_eps_iso = [f_eps]*3            
            self.UpdateFuncSide_spacetime(var, f_eps_iso, t=0.0)
            
            epsargs = {'r0':r0, 'r1':r1, 'mag_in':a_eps, 'mag_out':b_eps, 'f_eps':f_eps_iso}
            epsdic = {'type':vtype, 'args':epsargs, 'var':var}
            self.Materials.append(epsdic)
        elif vtype==FVTypes.EpsDiagBox: ##diagonal tensor
            r0 = args['r0']
            r1 = args['r1']
            a_eps = args['mag_in']
            b_eps = args['mag_out']
            exclude = False
            if 'exclude' in args:
                exclude = args['exclude']
            
            f_eps_0 = self.FunctionBox(r0, r1, a_eps[0], b_eps[0], exclude)
            f_eps_1 = self.FunctionBox(r0, r1, a_eps[1], b_eps[1], exclude)
            f_eps_2 = self.FunctionBox(r0, r1, a_eps[2], b_eps[2], exclude)
            f_eps_diag = [f_eps_0, f_eps_1, f_eps_2]            
            self.UpdateFuncSide_space(var, f_eps_diag)
            
            epsargs = {'r0':r0, 'r1':r1, 'mag_in':a_eps, 'mag_out':b_eps, 'f_eps':f_eps_diag}
            
            epsdic = {'type':vtype, 'args':epsargs, 'var':var}
            self.Materials.append(epsdic)
        elif vtype==FVTypes.EpsDiagSTvarBox: ##diagonal tensor
            r0 = args['r0']
            r1 = args['r1']
            a_eps = args['mag_in']
            b_eps = args['mag_out']
            exclude = False
            if 'exclude' in args:
                exclude = args['exclude']
            
            f_eps_0 = self.FunctionBox(r0, r1, a_eps[0], b_eps[0], exclude)
            f_eps_1 = self.FunctionBox(r0, r1, a_eps[1], b_eps[1], exclude)
            f_eps_2 = self.FunctionBox(r0, r1, a_eps[2], b_eps[2], exclude)
            f_eps_diag = [f_eps_0, f_eps_1, f_eps_2]            
            self.UpdateFuncSide_space(var, f_eps_diag)
            
            epsargs = {'r0':r0, 'r1':r1, 'mag_in':a_eps, 'mag_out':b_eps, 'f_eps':f_eps_diag}
            
            epsdic = {'type':vtype, 'args':epsargs, 'var':var}
            self.Materials.append(epsdic)
        elif vtype==FVTypes.MuIsoBox:
            r0 = args['r0']
            r1 = args['r1']
            a_mu = args['mag_in']
            b_mu = args['mag_out']
            exclude = False
            if 'exclude' in args:
                exclude = args['exclude']
            
            f_mu = self.FunctionBox(r0, r1, a_mu, b_mu, exclude)
            f_mu_iso = [f_mu]*3            
            self.UpdateFuncFace_space(var, f_mu_iso)
            
            muargs = {'r0':r0, 'r1':r1, 'mag_in':a_mu, 'mag_out':b_mu, 'f_mu':f_mu_iso}
            mudic = {'type':vtype, 'args':muargs, 'var':var}
            self.Materials.append(mudic)
        elif vtype==FVTypes.MuIsoSTvarBox:
            r0 = args['r0']
            r1 = args['r1']
            a_mu = args['mag_in']
            b_mu = args['mag_out']
            exclude = False
            if 'exclude' in args:
                exclude = args['exclude']
            
            f_mu = self.FunctionBox(r0, r1, a_mu, b_mu, exclude)
            f_mu_iso = [f_mu]*3            
            self.UpdateFuncFace_spacetime(var, f_mu_iso, t=0.0)
            
            muargs = {'r0':r0, 'r1':r1, 'mag_in':a_mu, 'mag_out':b_mu, 'f_mu':f_mu_iso}
            mudic = {'type':vtype, 'args':muargs, 'var':var}
            self.Materials.append(mudic)
        elif vtype==FVTypes.SigmaEDrudeIsoBox:
            r0 = args['r0']
            r1 = args['r1']
            sig_0 = args['sig_0']
            tau = args['tau']
            exclude = False
            if 'exclude' in args:
                exclude = args['exclude']

            f_sig = self.FunctionBox(r0, r1, sig_0, 0.0, exclude)
            f_sig_iso = [f_sig]*3            
            self.UpdateFuncSide_space(var, f_sig_iso)
            
            sigargs = {'r0':r0, 'r1':r1, 'sig_0':sig_0, 'tau':tau, 'f_sig':f_sig_iso}
            
            sigdic = {'type':vtype, 'args':sigargs, 'var':var}
            self.Materials.append(sigdic)
        elif vtype==FVTypes.SigmaEDrudeIsoSTvarBox:
            r0 = args['r0']
            r1 = args['r1']
            sig_0 = args['sig_0']
            tau = args['tau']
            exclude = False
            if 'exclude' in args:
                exclude = args['exclude']
            
            f_sig = self.FunctionBox(r0, r1, sig_0, 0.0, exclude)
            f_sig_iso = [f_sig]*3            
            self.UpdateFuncSide_spacetime(var, f_sig_iso, t=0.0)
            
            sigargs = {'r0':r0, 'r1':r1, 'sig_0':sig_0, 'tau':tau, 'f_sig':f_sig_iso}
            
            sigdic = {'type':vtype, 'args':sigargs, 'var':var}
            self.Materials.append(sigdic)
        elif vtype==FVTypes.SigmaMDrudeIsoBox:
            r0 = args['r0']
            r1 = args['r1']
            sig_0 = args['sig_0']
            tau = args['tau']
            exclude = False
            if 'exclude' in args:
                exclude = args['exclude']

            f_sig = self.FunctionBox(r0, r1, sig_0, 0.0, exclude)
            f_sig_iso = [f_sig]*3            
            self.UpdateFuncFace_space(var, f_sig_iso)
            
            sigargs = {'r0':r0, 'r1':r1, 'sig_0':sig_0, 'tau':tau, 'f_sig':f_sig_iso}
            
            sigdic = {'type':vtype, 'args':sigargs, 'var':var}
            self.Materials.append(sigdic)
        elif vtype==FVTypes.SigmaMDrudeIsoSTvarBox:
            r0 = args['r0']
            r1 = args['r1']
            sig_0 = args['sig_0']
            tau = args['tau']
            exclude = False
            if 'exclude' in args:
                exclude = args['exclude']
            
            f_sig = self.FunctionBox(r0, r1, sig_0, 0.0, exclude)
            f_sig_iso = [f_sig]*3            
            self.UpdateFuncFace_spacetime(var, f_sig_iso, t=0.0)
            
            sigargs = {'r0':r0, 'r1':r1, 'sig_0':sig_0, 'tau':tau, 'f_sig':f_sig_iso}
            
            sigdic = {'type':vtype, 'args':sigargs, 'var':var}
            self.Materials.append(sigdic)
        elif vtype==FVTypes.UPML:
            d_pml = args['d_pml']
            s_pml = args['s_pml']
            walls = self.GetWallsAllDic(d_pml, s_pml)

            sx_arr_s, sy_arr_s, sz_arr_s = self.GetUPMLFactor(walls, eps_or_s='s', side_or_face='side')
            sx_arr_f, sy_arr_f, sz_arr_f = self.GetUPMLFactor(walls, eps_or_s='s', side_or_face='face')
                    
            self.upml_c_fd0 = np.array([-np.imag(sy_arr_s[0]), -np.imag(sz_arr_s[1]), -np.imag(sx_arr_s[2])])
            self.upml_c_fd1 = np.array([np.imag(sx_arr_s[0]), np.imag(sy_arr_s[1]), np.imag(sz_arr_s[2])])
            self.upml_c_d = np.array([-np.imag(sz_arr_s[0]), -np.imag(sx_arr_s[1]), -np.imag(sy_arr_s[2])])

            self.upml_c_fh0 = np.array([-np.imag(sy_arr_f[0]), -np.imag(sz_arr_f[1]), -np.imag(sx_arr_f[2])])
            self.upml_c_fh1 = np.array([np.imag(sx_arr_f[0]), np.imag(sy_arr_f[1]), np.imag(sz_arr_f[2])])
            self.upml_c_h = np.array([-np.imag(sz_arr_f[0]), -np.imag(sx_arr_f[1]), -np.imag(sy_arr_f[2])])

            upmlargs = {'d_pml':d_pml, 's_pml':s_pml, 'walls':walls}
            upmlvars = {'c_fd0':self.upml_c_fd0, 'c_fd1':self.upml_c_fd1, 'c_d':self.upml_c_d, \
                        'c_fh0':self.upml_c_fh0, 'c_fh1':self.upml_c_fh1, 'c_h':self.upml_c_h}
            upmldic = {'type':vtype, 'args':upmlargs, 'var':upmlvars}
            self.Materials.append(upmldic)
        
        elif vtype==FVTypes.BCond:
            """Example: 
            bcs = FVTypes.PEC  --> everywhere PEC
            bcs = FVTypes.PBC  --> everywhere PBC
            bcs[0]=[FVTypes.PEC, FVTypes.PEC]  ---> PEC at x- ad x+
            bcs[1]=[FVTypes.PBC, FVTypes.PBC]  ---> PBC at y
            """
            bcs = args['bcs']
            if isinstance(bcs, list):
                assert len(bcs)==self.n_dim
            else:
                assert bcs in [FVTypes.PEC, FVTypes.PMC, FVTypes.PBC]
            self.bcs = bcs
        return        



    def GetMaterialTypes(self):
        AllMaterialTypes = []
        for M in self.Materials:
            if M['type'] not in AllMaterialTypes:
                AllMaterialTypes.append(M['type'])
        self.AllMaterialTypes = AllMaterialTypes
        

    def AssociateFieldTypesToVars(self):
        fieldTypeVar = {}
        for F in self.Fields:
            assert F['type'] not in fieldTypeVar
            fieldTypeVar[F['type']] = F['var']
        self.fieldTypeVar = fieldTypeVar


    def UpdateSources(self, t , em_type='e'):
        """ em_type: e/m/h
        """
        if em_type=='e':
            for jdic in self.Sources:
                if jdic['type'] in [FVTypes.JePoint, FVTypes.JeSheet, FVTypes.JeBox]:
                    jargs = jdic['args']
                    var = jdic['var']
                    f_t = jargs['f_t']
                    J_0 = jargs['V_0']
                    self.UpdateSeperableFunc_Time(var, J_0, f_t, t)
        elif em_type in ['m', 'h']:
            for jdic in self.Sources:
                if jdic['type'] in [FVTypes.JmPoint, FVTypes.JmSheet, FVTypes.JmBox]:
                    jargs = jdic['args']
                    var = jdic['var']
                    f_t = jargs['f_t']
                    J_0 = jargs['V_0']
                    self.UpdateSeperableFunc_Time(var, J_0, f_t, t)

    def sumSources(self, em_type='e'):
        if em_type=='e':
            j_tot = [0.0]*3
            for jdic in self.Sources:
                if jdic['type'] in [FVTypes.JePoint, FVTypes.JeSheet, FVTypes.JeBox]:
                    var = jdic['var']
                    for i in range(3):
                        j_tot[i] = var[i] + j_tot[i]
            return j_tot
        elif em_type in ['m', 'h']:
            j_tot = [0.0]*3
            for jdic in self.Sources:
                if jdic['type'] in [FVTypes.JmPoint, FVTypes.JmSheet, FVTypes.JmBox]:
                    var = jdic['var']
                    for i in range(3):
                        j_tot[i] = var[i] + j_tot[i]
            return j_tot
        

    def updateTvarMaterials(self, t, em_type='e'):
        if em_type=='e':
            for matdic in self.Materials:
                if matdic['type']==FVTypes.EpsIsoSTvarBox:
                    var = matdic['var']
                    f_eps_iso = matdic['args']['f_eps']
                    self.UpdateFuncSide_spacetime(var, f_eps_iso, t)
                if matdic['type']==FVTypes.SigmaEDrudeIsoSTvarBox:
                    var = matdic['var']
                    f_sig_iso = matdic['args']['f_sig']
                    self.UpdateFuncSide_spacetime(var, f_sig_iso, t)
                if matdic['type']==FVTypes.EpsDiagSTvarBox:
                    var = matdic['var']
                    f_eps_diag = matdic['args']['f_eps']
                    self.UpdateFuncSide_spacetime(var, f_eps_diag, t)
        else:
            for matdic in self.Materials:
                if matdic['type']==FVTypes.MuIsoSTvarBox:
                    var = matdic['var']
                    f_mu_iso = matdic['args']['f_mu']
                    self.UpdateFuncFace_spacetime(var, f_mu_iso, t)
                if matdic['type']==FVTypes.SigmaMDrudeIsoSTvarBox:
                    var = matdic['var']
                    f_sig_iso = matdic['args']['f_sig']
                    self.UpdateFuncFace_spacetime(var, f_sig_iso, t)
            

    def sumMaterials(self, mat_type):
        ## vacuum should be defined in only one of the box definitions 
        ## TODO: this condition to be relaxed
        assert mat_type in ['eps', 'mu']
        if mat_type=='eps':
            has_eps = False
            mat_tot = [0.0]*3
            for matdic in self.Materials:
                if matdic['type'] in [FVTypes.EpsIsoBox, FVTypes.EpsIsoSTvarBox, \
                        FVTypes.EpsDiagBox, FVTypes.EpsDiagSTvarBox]:
                    has_eps = True
                    var = matdic['var']
                    for i in range(3):
                        mat_tot[i] = var[i] + mat_tot[i]
            if has_eps:
                return mat_tot
            else:
                return [1.0]*3
        elif mat_type=='mu':
            has_mu = False
            mat_tot = [0.0]*3
            for matdic in self.Materials:
                if matdic['type'] in [FVTypes.MuIsoBox, FVTypes.MuIsoSTvarBox]:
                    has_mu = True
                    var = matdic['var']
                    for i in range(3):
                        mat_tot[i] = var[i] + mat_tot[i]
            if has_mu:
                return mat_tot
            else:
                return [1.0]*3
    

    ##TODO: update conduction currents
    def UpdateE(self, t):
        H = self.fieldTypeVar[FVTypes.H]
        curl_H = self.GetCurlFace(H)
        self.curl_H = curl_H
        self.sumCurlFaceWalls('self.curl_H')
        self.curl_H = None
        
        ## dD/dt = curl H - Je
        D = self.fieldTypeVar[FVTypes.D]
        Je_tot = self.sumSources(em_type='e')
        if self.bcs==FVTypes.PBC:
            self.ResetPBC(D)
            self.ResetPBC(Je_tot)
        self.Update_adAdt_bB(D, curl_H, C_list=[Je_tot], c_list=[-np.ones(3)])

        if self.bcs==FVTypes.PBC:
            self.ApplyPBC(D)
        
        ## eps*E = D
        E = self.fieldTypeVar[FVTypes.E]
        eps_tot = self.sumMaterials(mat_type='eps')
        self.Update_aA_bB(E, D, a=eps_tot)
                
        ## PEC
        if self.bcs is None or self.bcs==FVTypes.PEC:
            self.ResetSideWalls(E)
        return


    def UpdateH(self, t):
        E = self.fieldTypeVar[FVTypes.E]
        curl_E = self.GetCurlSide(E)
        
        ## dB/dt = -curl E - Jm
        B = self.fieldTypeVar[FVTypes.B]
        Jm_tot = self.sumSources(em_type='m')
        self.Update_adAdt_bB(B, curl_E, C_list=[Jm_tot], b=-np.ones(3), c_list=[-np.ones(3)])
        
        ## mu*H = B
        H = self.fieldTypeVar[FVTypes.H]
        mu_tot = self.sumMaterials(mat_type='mu')
        self.Update_aA_bB(H, B, a=mu_tot)
        return


    def UpdateE_UPML(self, t):
        H = self.fieldTypeVar[FVTypes.H]
        curl_H = self.GetCurlFace(H)
        self.curl_H = curl_H
        self.sumCurlFaceWalls('self.curl_H')
        self.curl_H = None
        
        ## dF_D/dt + s_fd*F_D = curl H - Je
        F_D = self.upmlvars['F_D']
        dF_D = self.upmlvars['dF_D']
        c_fd0 = self.upmlvars['c_fd0']
        Je_tot = self.sumSources(em_type='e')
        self.Update_aA_bB(dF_D, curl_H, C_list=[F_D, Je_tot], c_list=[c_fd0, -np.ones(3)])
        ## dD/dt + s_d*D = dF_D/dt + s_fd*F_D
        D = self.fieldTypeVar[FVTypes.D]
        c_fd1 = self.upmlvars['c_fd1']
        c_d = self.upmlvars['c_d']
        self.Update_adAdt_bB(D, dF_D, C_list=[D, F_D], c_list=[c_d, c_fd1])
        self.Update_adAdt_bB(F_D, dF_D)
        ## eps*E = D
        E = self.fieldTypeVar[FVTypes.E]
        eps_tot = self.sumMaterials(mat_type='eps')
        self.Update_aA_bB(E, D, a=eps_tot)
        ## PEC
        if self.bcs is None or self.bcs==FVTypes.PEC:
            self.ResetSideWalls(E)
            self.ResetSideWalls(D)
        else:
            raise NotImplementedError()
        return


    def UpdateH_UPML(self, t):
        E = self.fieldTypeVar[FVTypes.E]
        curl_E = self.GetCurlSide(E)
        
        ## dF_M/dt + s_fm*F_M = -curl E - Jm
        F_H = self.upmlvars['F_H']
        dF_H = self.upmlvars['dF_H']
        c_fh0 = self.upmlvars['c_fh0']
        Jm_tot = self.sumSources(em_type='m')
        self.Update_aA_bB(dF_H, curl_E, C_list=[F_H, Jm_tot], b=-np.ones(3), c_list=[c_fh0, -np.ones(3)])
        ## dH/dt + s_h*H = dF_H/dt + s_fd*F_H
        B = self.fieldTypeVar[FVTypes.B]
        c_fh1 = self.upmlvars['c_fh1']
        c_h = self.upmlvars['c_h']
        self.Update_adAdt_bB(B, dF_H, C_list=[B, F_H], c_list=[c_h, c_fh1])
        self.Update_adAdt_bB(F_H, dF_H)
        ## mu*H = B
        H = self.fieldTypeVar[FVTypes.H]
        mu_tot = self.sumMaterials(mat_type='mu')
        self.Update_aA_bB(H, B, a=mu_tot)

    
    #TODO
    def UpdateE_UPML_Lossy(self, t):
        return

    #TODO
    def UpdateH_UPML_Lossy(self, t):
        return


    def SetUPmlVars(self):
        LossyPml = False
        
        for M in self.Materials:
            if M['type']==FVTypes.UPML:
                upmlvars = M['var']
                if not LossyPml:
                    self.AllocateSideArr_list(["upml_dF_D", "upml_F_D"])
                    self.AllocateFaceArr_list(["upml_dF_H", "upml_F_H"])     
                    upmlvars['dF_D'] = self.upml_dF_D
                    upmlvars['dF_H'] = self.upml_dF_H
                    upmlvars['F_D'] = self.upml_F_D
                    upmlvars['F_H'] = self.upml_F_H
                    
                    self.upmlvars = upmlvars
                else:
                    raise NotImplementedError()



    def StepFields(self, n_t):
        dt = self.dt
        if self.vbose:
            print('self.ChunkInd:', self.ChunkInd, 'dt:', dt, 'n_t:', n_t, 't_final:', n_t*dt)
            print('starting time: ', time.perf_counter())
        self.GetMaterialTypes()
        self.AssociateFieldTypesToVars()
        
        has_pml = False
        if FVTypes.UPML in self.AllMaterialTypes:
            assert not has_pml
            has_pml = True
            self.SetUPmlVars()
        if self.vbose:
            print('has_pml:', has_pml)
                
        if not has_pml:
            for i in range(n_t):
                t = float(i*dt)
                self.UpdateSources(t , em_type='e')
                self.updateTvarMaterials(t, em_type='e')
                self.UpdateE(t)

                t += dt/2
                self.UpdateSources(t , em_type='m')
                self.UpdateH(t)

                if i%self.save_every==0:
                    self.SaveOutputs()
                    self.n_saved += 1
        else:
            for i in range(n_t):
                t = float(i*dt)
                self.UpdateSources(t , em_type='e')
                self.updateTvarMaterials(t, em_type='e')
                self.UpdateE_UPML(t)

                t += dt/2
                self.UpdateSources(t , em_type='m')
                self.UpdateH_UPML(t)

                if i%self.save_every==0:
                    self.SaveOutputs()
                    self.n_saved += 1
                    
        self.DumpRemainingOutputsToFile()

        if self.vbose:
            print('finishing time: ', time.perf_counter())

        
    def LogInfo(self):
        print('self.ChunkInd:', self.ChunkInd, 'self.ChunkIndTot:', self.ChunkIndTot)
        print('self.ChunkIndNB:', self.ChunkIndNB)
        print('self.N:', self.N, 'self.dr:', self.dr, 'self.r0:', self.r0, 'self.r1:', self.r1)
        print('self.cmdPipe:', self.cmdPipe)
        print('self.commPipesNBDicSend: ', self.commPipesNBDicSend)
        print('self.commPipesNBDicRecv: ', self.commPipesNBDicRecv)
        print('self.quitLoop: ', self.quitLoop)
    

    def run(self):
        if self.vbose:
            print('starting process: ', self.ChunkInd)
        while not self.quitLoop:
            ## wait for commands from the main process
            cmd = self.cmdPipe.recv()
            ##process command
            self.runCmd(cmd)
            sys.stdout.flush()
        if self.vbose:
            print('Exiting process: ', self.ChunkInd)
        sys.stdout.flush()
            
    
    def runCmd(self, cmd):
        if cmd[0]==FdtdCmds.stringcmd:
            tag = cmd[1]
            str_cmd = cmd[2]
            exec(str_cmd)
            str_return = cmd[3]
            if str_return!=None:
                var_return = [None]
                exec('var_return[0] = '+str_return)
                self.cmdPipe.send([self.ChunkInd, tag, ['ok', var_return[0]]])
            else:
                self.cmdPipe.send([self.ChunkInd, tag, ['ok', None]])
        elif cmd[0]==FdtdCmds.quit:
            self.quitLoop = True
            tag = cmd[1]
            self.cmdPipe.send([self.ChunkInd, tag, ['ok', None]])
            if self.vbose:
                print('FdtdCmds.quit message received.', self.ChunkInd, 'self.quitLoop:', self.quitLoop)


    def SetCartesianIndex(self, ChunkInd=None, ChunkIndTot=None):
        """ ChunkInd: np.array(n_dim)
        """
        self.ChunkInd = ChunkInd
        self.ChunkIndTot = ChunkIndTot
        
    def SetCartesianIndexNeighbors(self):
        ChunkIndNB = {'p':[None]*self.n_dim, 'n':[None]*self.n_dim}
        ChunkInd0 = np.zeros(self.n_dim)
        if self.n_dim==3:
            chind_xp = self.ChunkInd + np.array([1, 0, 0], dtype=int)
            if np.all(chind_xp < self.ChunkIndTot):
                ChunkIndNB['p'][0] = chind_xp
            chind_xn = self.ChunkInd - np.array([1, 0, 0], dtype=int)
            if np.all(chind_xn >= ChunkInd0):
                ChunkIndNB['n'][0] = chind_xn

            chind_yp = self.ChunkInd + np.array([0, 1, 0], dtype=int)
            if np.all(chind_yp < self.ChunkIndTot):
                ChunkIndNB['p'][1] = chind_yp
            chind_yn = self.ChunkInd - np.array([0, 1, 0], dtype=int)
            if np.all(chind_yn >= ChunkInd0):
                ChunkIndNB['n'][1] = chind_yn

            chind_zp = self.ChunkInd + np.array([0, 0, 1], dtype=int)
            if np.all(chind_zp < self.ChunkIndTot):
                ChunkIndNB['p'][2] = chind_zp
            chind_zn = self.ChunkInd - np.array([0, 0, 1], dtype=int)
            if np.all(chind_zn >= ChunkInd0):
                ChunkIndNB['n'][2] = chind_zn
                
        elif self.n_dim==2:
            chind_xp = self.ChunkInd + np.array([1, 0], dtype=int)
            if np.all(chind_xp < self.ChunkIndTot):
                ChunkIndNB['p'][0] = chind_xp
            chind_xn = self.ChunkInd - np.array([1, 0], dtype=int)
            if np.all(chind_xn >= ChunkInd0):
                ChunkIndNB['n'][0] = chind_xn

            chind_yp = self.ChunkInd + np.array([0, 1], dtype=int)
            if np.all(chind_yp < self.ChunkIndTot):
                ChunkIndNB['p'][1] = chind_yp
            chind_yn = self.ChunkInd - np.array([0, 1], dtype=int)
            if np.all(chind_yn >= ChunkInd0):
                ChunkIndNB['n'][1] = chind_yn
        
        elif self.n_dim==1:
            chind_xp = self.ChunkInd + np.array([1], dtype=int)
            if np.all(chind_xp < self.ChunkIndTot):
                ChunkIndNB['p'][0] = chind_xp
            chind_xn = self.ChunkInd - np.array([1], dtype=int)
            if np.all(chind_xn >= ChunkInd0):
                ChunkIndNB['n'][0] = chind_xn

        self.ChunkIndNB = ChunkIndNB
        

    def AddCommPipesNBSend(self, grid_ind_nb, commPipeSend):
        """ neighbor to neighbor communication pipes
        """
        self.commPipesNBDicSend[tuple(grid_ind_nb)] = commPipeSend

    def AddCommPipesNBRecv(self, grid_ind_nb, commPipeRecv):
        """ neighbor to neighbor communication pipes
        """
        self.commPipesNBDicRecv[tuple(grid_ind_nb)] = commPipeRecv


    def SetCommandPipe(self, commandPipe):
        """ communication between the main process and this cartesian chunk
        """
        self.cmdPipe = commandPipe
        

    def AllocateSideArr_list(self, str_list):
        for A_str in str_list:
            self.allocateSideArrays(A_str)

    def AllocateFaceArr_list(self, str_list):
        for A_str in str_list:
            self.allocateFaceArrays(A_str)

    def allocateSideArrays(self, A_str):
        cmd_x = "%sx = np.zeros(tuple(self.Nsx), dtype=self.dtype)"%A_str
        cmd_y = "%sy = np.zeros(tuple(self.Nsy), dtype=self.dtype)"%A_str
        cmd_z = "%sz = np.zeros(tuple(self.Nsz), dtype=self.dtype)"%A_str
        exec(cmd_x)
        exec(cmd_y)
        exec(cmd_z)
        
        cmd = "self.%s = [%sx, %sy, %sz]"%(A_str, A_str, A_str, A_str)
        exec(cmd)
        if self.vbose:
            exec("print('self.%s[0].shape:', self.%s[0].shape)"%(A_str, A_str))

    def allocateFaceArrays(self, A_str):
        cmd_x = "%sx = np.zeros(tuple(self.Nfx), dtype=self.dtype)"%A_str
        cmd_y = "%sy = np.zeros(tuple(self.Nfy), dtype=self.dtype)"%A_str
        cmd_z = "%sz = np.zeros(tuple(self.Nfz), dtype=self.dtype)"%A_str
        exec(cmd_x)
        exec(cmd_y)
        exec(cmd_z)

        cmd = "self.%s = [%sx, %sy, %sz]"%(A_str, A_str, A_str, A_str)
        exec(cmd)
        if self.vbose:
            exec("print('self.%s[0].shape:', self.%s[0].shape)"%(A_str, A_str))


    def UpdateSide_dAdt_CurlB_C(self, A, B, C_list=None, a=None, b=None, c_list=None):
        """ Wall indices for A are not updated, it will be later forced to 0
            a*d/dt A = b*curl(B) + c*C
            a, A, b, c, C are side elems
            B is face elem
            it starts with D then C then B
            A is on the left hand side.. the rest on the sight hand side
            
            This function should be replaced with GetCurlFace in the parallel version
        """
        if a is None:
            a = np.ones(3)
        if b is None:
            b = np.ones(3)
        Nfx, Nfy, Nfz = self.Nfx, self.Nfy, self.Nfz
        
        ax, ay, az = a
        bx, by, bz = b

        Ax, Ay, Az = A
        Bx, By, Bz = B

        #print('side', '-'*50)
        #print(np.max(np.abs(Ax)), np.max(np.abs(Ay)), np.max(np.abs(Az)))
        #print(np.max(np.abs(A[0])), np.max(np.abs(A[1])), np.max(np.abs(A[2])))


        dx, dy, dz = [None]*3
        if self.n_dim==3:
            dx, dy, dz = self.dr[0], self.dr[1], self.dr[2]
        elif self.n_dim==2:
            dx, dy, dz = self.dr[0], self.dr[1], None
        elif self.n_dim==1:
            dx, dy, dz = self.dr[0], None, None
        dt = self.dt
        
        if C_list!=None:
            if c_list is None:
                c_list = [None]*len(C_list)
            for i in range(len(C_list)):
                C = C_list[i]
                c = c_list[i]
                if c is None:
                    c = np.ones(3)
                    
                Cx, Cy, Cz = C
                cx, cy, cz = c
                Ax += Cx*(cx*dt/ax)
                Ay += Cy*(cy*dt/ay)
                Az += Cz*(cz*dt/az)

        ax_, ay_, az_ = ax, ay, az
        bx_, by_, bz_ = bx, by, bz

        if self.n_dim==3:
            if type(ax)==np.ndarray:
                ax_ = ax[:, 1:-1, 1:-1]
            if type(bx)==np.ndarray:
                bx_ = bx[:, 1:-1, 1:-1]

            if type(ay)==np.ndarray:
                ay_ = ay[1:-1, :, 1:-1]
            if type(by)==np.ndarray:
                by_ = by[1:-1, :, 1:-1]

            if type(az)==np.ndarray:
                az_ = az[1:-1, 1:-1, :]
            if type(bz)==np.ndarray:
                bz_ = bz[1:-1, 1:-1, :]
        elif self.n_dim==2:
            if type(ax)==np.ndarray:
                ax_ = ax[:, 1:-1]
            if type(bx)==np.ndarray:
                bx_ = bx[:, 1:-1]

            if type(ay)==np.ndarray:
                ay_ = ay[1:-1, :]
            if type(by)==np.ndarray:
                by_ = by[1:-1, :]

            if type(az)==np.ndarray:
                az_ = az[1:-1, 1:-1]
            if type(bz)==np.ndarray:
                bz_ = bz[1:-1, 1:-1]
        elif self.n_dim==1:
            if type(ax)==np.ndarray:
                ax_ = ax[:]
            if type(bx)==np.ndarray:
                bx_ = bx[:]

            if type(ay)==np.ndarray:
                ay_ = ay[1:-1]
            if type(by)==np.ndarray:
                by_ = by[1:-1]

            if type(az)==np.ndarray:
                az_ = az[1:-1]
            if type(bz)==np.ndarray:
                bz_ = bz[1:-1]

        if self.n_dim==3:
            Ax[:, 1:-1, 1:-1] -= (By[:, 1:-1, 1:Nfy[2]] - By[:, 1:-1, 0:Nfy[2]-1])*(bx_*dt/(dz*ax_)) 
            Ax[:, 1:-1, 1:-1] += (Bz[:, 1:Nfz[1], 1:-1] - Bz[:, 0:Nfz[1]-1, 1:-1])*(bx_*dt/(dy*ax_)) 

            Ay[1:-1, :, 1:-1] -= (Bz[1:Nfz[0], :, 1:-1] - Bz[0:Nfz[0]-1, :, 1:-1])*(by_*dt/(dx*ay_)) 
            Ay[1:-1, :, 1:-1] += (Bx[1:-1, :, 1:Nfx[2]] - Bx[1:-1, :, 0:Nfx[2]-1])*(by_*dt/(dz*ay_)) 

            Az[1:-1, 1:-1, :] -= (Bx[1:-1, 1:Nfx[1], :] - Bx[1:-1, 0:Nfx[1]-1, :])*(bz_*dt/(dy*az_)) 
            Az[1:-1, 1:-1, :] += (By[1:Nfy[0], 1:-1, :] - By[0:Nfy[0]-1, 1:-1, :])*(bz_*dt/(dx*az_)) 
        elif self.n_dim==2:
            #Ax[:, 1:-1] -= (By[:, 1:-1    ] - By[:, 1:-1      ])*(bx_*dt/(dz*ax_)) 
            Ax[:, 1:-1] += (Bz[:, 1:Nfz[1]] - Bz[:, 0:Nfz[1]-1])*(bx_*dt/(dy*ax_)) 

            Ay[1:-1, :] -= (Bz[1:Nfz[0], :] - Bz[0:Nfz[0]-1, :])*(by_*dt/(dx*ay_)) 
            #Ay[1:-1, :] += (Bx[1:-1, :    ] - Bx[1:-1, :      ])*(by_*dt/(dz*ay_)) 

            Az[1:-1, 1:-1] -= (Bx[1:-1, 1:Nfx[1]] - Bx[1:-1, 0:Nfx[1]-1])*(bz_*dt/(dy*az_)) 
            Az[1:-1, 1:-1] += (By[1:Nfy[0], 1:-1] - By[0:Nfy[0]-1, 1:-1])*(bz_*dt/(dx*az_)) 
        elif self.n_dim==1:
            #Ax[:] -= (By[:] - By[:])*(bx_*dt/(dz*ax_)) 
            #Ax[:] += (Bz[:] - Bz[:])*(bx_*dt/(dy*ax_)) 

            Ay[1:-1] -= (Bz[1:Nfz[0]] - Bz[0:Nfz[0]-1])*(by_*dt/(dx*ay_)) 
            #Ay[1:-1] += (Bx[1:-1    ] - Bx[1:-1      ])*(by_*dt/(dz*ay_)) 

            #Az[1:-1] -= (Bx[1:-1    ] - Bx[1:-1      ])*(bz_*dt/(dy*az_)) 
            Az[1:-1] += (By[1:Nfy[0]] - By[0:Nfy[0]-1])*(bz_*dt/(dx*az_)) 
                
        ##walls
        ax_ii0, ax_ii1, ax_i0i, ax_i1i = [ax]*4
        bx_ii0, bx_ii1, bx_i0i, bx_i1i = [bx]*4
        ay_0ii, ay_1ii, ay_ii0, ay_ii1 = [ay]*4
        by_0ii, by_1ii, by_ii0, by_ii1 = [by]*4
        az_i0i, az_i1i, az_0ii, az_1ii = [az]*4
        bz_i0i, bz_i1i, bz_0ii, bz_1ii = [bz]*4
        ax_i0, ax_i1 = [ax]*2
        bx_i0, bx_i1 = [bx]*2
        ay_0i, ay_1i = [ay]*2
        by_0i, by_1i = [by]*2
        az_i0, az_i1, az_0i, az_1i = [az]*4
        bz_i0, bz_i1, bz_0i, bz_1i = [bz]*4
        ay_0, ay_1 = [ay]*2
        by_0, by_1 = [by]*2
        az_0, az_1 = [az]*2
        bz_0, bz_1 = [bz]*2
        if self.n_dim==3:
            if type(ax)==np.ndarray:
                ax_ii0 = ax[:, :, 0]
                ax_ii1 = ax[:, :, 1]
                ax_i0i = ax[:, 0, :]
                ax_i1i = ax[:, 1, :]
            if type(bx)==np.ndarray:
                bx_ii0 = bx[:, :, 0]
                bx_ii1 = bx[:, :, 1]
                bx_i0i = bx[:, 0, :]
                bx_i1i = bx[:, 1, :]

            if type(ay)==np.ndarray:
                ay_0ii = ay[0, :, :]
                ay_1ii = ay[1, :, :]
                ay_ii0 = ay[:, :, 0]
                ay_ii1 = ay[:, :, 1]
            if type(by)==np.ndarray:
                by_0ii = by[0, :, :]
                by_1ii = by[1, :, :]
                by_ii0 = by[:, :, 0]
                by_ii1 = by[:, :, 1]

            if type(az)==np.ndarray:
                az_i0i = az[:, 0, :]
                az_i1i = az[:, 1, :]
                az_0ii = az[0, :, :]
                az_1ii = az[1, :, :]
            if type(bz)==np.ndarray:
                bz_i0i = bz[:, 0, :]
                bz_i1i = bz[:, 1, :]
                bz_0ii = bz[0, :, :]
                bz_1ii = bz[1, :, :]

        elif self.n_dim==2:
            if type(ax)==np.ndarray:
                ax_i0 = ax[:, 0]
                ax_i1 = ax[:, 1]
            if type(bx)==np.ndarray:
                bx_i0 = bx[:, 0]
                bx_i1 = bx[:, 1]

            if type(ay)==np.ndarray:
                ay_0i = ay[0, :]
                ay_1i = ay[1, :]
            if type(by)==np.ndarray:
                by_0i = by[0, :]
                by_1i = by[1, :]

            if type(az)==np.ndarray:
                az_i0 = az[:, 0]
                az_i1 = az[:, 1]
                az_0i = az[0, :]
                az_1i = az[1, :]
            if type(bz)==np.ndarray:
                bz_i0 = bz[:, 0]
                bz_i1 = bz[:, 1]
                bz_0i = bz[0, :]
                bz_1i = bz[1, :]
        elif self.n_dim==1:
            if type(ay)==np.ndarray:
                ay_0 = ay[0]
                ay_1 = ay[1]
            if type(by)==np.ndarray:
                by_0 = by[0]
                by_1 = by[1]

            if type(az)==np.ndarray:
                az_0 = az[0]
                az_1 = az[1]
            if type(bz)==np.ndarray:
                bz_0 = bz[0]
                bz_1 = bz[1]
                
        if self.n_dim==3:
            Ax[:, :,  0] -= ( By[:, :, 0       ])*(bx_ii0*dt/(dz*ax_ii0)) 
            Ax[:, :, -1] -= (-By[:, :, Nfy[2]-1])*(bx_ii1*dt/(dz*ax_ii1)) 
            Ax[:, 0 , :] += ( Bz[:, 0,        :])*(bx_i0i*dt/(dy*ax_i0i)) 
            Ax[:,-1 , :] += (-Bz[:, Nfz[1]-1, :])*(bx_i1i*dt/(dy*ax_i1i)) 


            Ay[ 0, :, :] -= ( Bz[0,        :, :])*(by_0ii*dt/(dx*ay_0ii)) 
            Ay[-1, :, :] -= (-Bz[Nfz[0]-1, :, :])*(by_1ii*dt/(dx*ay_1ii)) 
            Ay[:, :,  0] += ( Bx[:, :,        0])*(by_ii0*dt/(dz*ay_ii0)) 
            Ay[:, :, -1] += (-Bx[:, :, Nfx[2]-1])*(by_ii1*dt/(dz*ay_ii1)) 

            Az[:, 0 , :] -= ( Bx[:, 0,        :])*(bz_i0i*dt/(dy*az_i0i)) 
            Az[:, -1, :] -= (-Bx[:, Nfx[1]-1, :])*(bz_i1i*dt/(dy*az_i1i)) 
            Az[0 , :, :] += ( By[0,        :, :])*(bz_0ii*dt/(dx*az_0ii)) 
            Az[-1, :, :] += (-By[Nfy[0]-1, :, :])*(bz_1ii*dt/(dx*az_1ii)) 
        elif self.n_dim==2:
            #Ax[:, 1:-1] -= (By[:, 1:-1    ] - By[:, 1:-1      ])*(bx_*dt/(dz*ax_)) 
            Ax[:, 0] += ( Bz[:, 0       ])*(bx_i0*dt/(dy*ax_i0)) 
            Ax[:,-1] += (-Bz[:, Nfz[1]-1])*(bx_i1*dt/(dy*ax_i1)) 

            Ay[0 , :] -= ( Bz[0,        :])*(by_0i*dt/(dx*ay_0i)) 
            Ay[-1, :] -= (-Bz[Nfz[0]-1, :])*(by_1i*dt/(dx*ay_1i)) 
            #Ay[1:-1, :] += (Bx[1:-1, :    ] - Bx[1:-1, :      ])*(by_*dt/(dz*ay_)) 

            Az[:,  0] -= ( Bx[:, 0       ])*(bz_i0*dt/(dy*az_i0)) 
            Az[:, -1] -= (-Bx[:, Nfx[1]-1])*(bz_i1*dt/(dy*az_i1)) 
            Az[0 , :] += ( By[0,        :])*(bz_0i*dt/(dx*az_0i)) 
            Az[-1, :] += (-By[Nfy[0]-1, :])*(bz_1i*dt/(dx*az_1i)) 
        elif self.n_dim==1:
            #Ax[:] -= (By[:] - By[:])*(bx_*dt/(dz*ax_)) 
            #Ax[:] += (Bz[:] - Bz[:])*(bx_*dt/(dy*ax_)) 

            Ay[0 ] -= ( Bz[0       ])*(by_0*dt/(dx*ay_0)) 
            Ay[-1] -= (-Bz[Nfz[0]-1])*(by_1*dt/(dx*ay_1)) 
            #Ay[1:-1] += (Bx[1:-1    ] - Bx[1:-1      ])*(by_*dt/(dz*ay_)) 

            #Az[1:-1] -= (Bx[1:-1    ] - Bx[1:-1      ])*(bz_*dt/(dy*az_)) 
            Az[0 ] += ( By[0       ])*(bz_0*dt/(dx*az_0)) 
            Az[-1] += (-By[Nfy[0]-1])*(bz_1*dt/(dx*az_1)) 

    
    def GetCurlFace(self, B, b=None):
        """ b*curl(B) B is Face element
            returns side element
        """
        if b is None:
            b = np.ones(3)

        Nfx, Nfy, Nfz = self.Nfx, self.Nfy, self.Nfz
        Bx, By, Bz = B
        bx, by, bz = b

        bx_, by_, bz_ = bx, by, bz

        if self.n_dim==3:
            if type(bx)==np.ndarray:
                bx_ = bx[:, 1:-1, 1:-1]
            if type(by)==np.ndarray:
                by_ = by[1:-1, :, 1:-1]
            if type(bz)==np.ndarray:
                bz_ = bz[1:-1, 1:-1, :]
        elif self.n_dim==2:
            if type(bx)==np.ndarray:
                bx_ = bx[:, 1:-1]
            if type(by)==np.ndarray:
                by_ = by[1:-1, :]
            if type(bz)==np.ndarray:
                bz_ = bz[1:-1, 1:-1]
        elif self.n_dim==1:
            if type(bx)==np.ndarray:
                bx_ = bx[:, 1:-1]
            if type(by)==np.ndarray:
                by_ = by[1:-1, :]
            if type(bz)==np.ndarray:
                bz_ = bz[1:-1, 1:-1]

        dx, dy, dz = [None]*3
        if self.n_dim==3:
            dx, dy, dz = self.dr[0], self.dr[1], self.dr[2]
        elif self.n_dim==2:
            dx, dy, dz = self.dr[0], self.dr[1], None
        elif self.n_dim==1:
            dx, dy, dz = self.dr[0], None, None

        Ax = np.zeros(self.Nsx, dtype=self.dtype)
        Ay = np.zeros(self.Nsy, dtype=self.dtype)
        Az = np.zeros(self.Nsz, dtype=self.dtype)
        
        if self.n_dim==3:
            Ax[:, 1:-1, 1:-1] -= (By[:, 1:-1, 1:Nfy[2]] - By[:, 1:-1, 0:Nfy[2]-1])*(bx_/dz) 
            Ax[:, 1:-1, 1:-1] += (Bz[:, 1:Nfz[1], 1:-1] - Bz[:, 0:Nfz[1]-1, 1:-1])*(bx_/dy) 

            Ay[1:-1, :, 1:-1] -= (Bz[1:Nfz[0], :, 1:-1] - Bz[0:Nfz[0]-1, :, 1:-1])*(by_/dx) 
            Ay[1:-1, :, 1:-1] += (Bx[1:-1, :, 1:Nfx[2]] - Bx[1:-1, :, 0:Nfx[2]-1])*(by_/dz) 

            Az[1:-1, 1:-1, :] -= (Bx[1:-1, 1:Nfx[1], :] - Bx[1:-1, 0:Nfx[1]-1, :])*(bz_/dy) 
            Az[1:-1, 1:-1, :] += (By[1:Nfy[0], 1:-1, :] - By[0:Nfy[0]-1, 1:-1, :])*(bz_/dx)
        elif self.n_dim==2:
            #Ax[:, 1:-1] -= (By[:, 1:-1    ] - By[:, 1:-1      ])*(bx_/dz) 
            Ax[:, 1:-1] += (Bz[:, 1:Nfz[1]] - Bz[:, 0:Nfz[1]-1])*(bx_/dy) 

            Ay[1:-1, :] -= (Bz[1:Nfz[0], :] - Bz[0:Nfz[0]-1, :])*(by_/dx) 
            #Ay[1:-1, :] += (Bx[1:-1, :    ] - Bx[1:-1, :      ])*(by_/dz) 

            Az[1:-1, 1:-1] -= (Bx[1:-1, 1:Nfx[1]] - Bx[1:-1, 0:Nfx[1]-1])*(bz_/dy) 
            Az[1:-1, 1:-1] += (By[1:Nfy[0], 1:-1] - By[0:Nfy[0]-1, 1:-1])*(bz_/dx)
        elif self.n_dim==1:
            #Ax[:] -= (By[:] - By[:])*(bx_/dz) 
            #Ax[:] += (Bz[:] - Bz[:])*(bx_/dy) 

            Ay[1:-1] -= (Bz[1:Nfz[0]] - Bz[0:Nfz[0]-1])*(by_/dx) 
            #Ay[1:-1] += (Bx[1:-1    ] - Bx[1:-1      ])*(by_/dz) 

            #Az[1:-1] -= (Bx[1:-1    ] - Bx[1:-1      ])*(bz_/dy) 
            Az[1:-1] += (By[1:Nfy[0]] - By[0:Nfy[0]-1])*(bz_/dx)
        
        ##walls
        bx_ii0, bx_ii1, bx_i0i, bx_i1i = [bx]*4
        by_0ii, by_1ii, by_ii0, by_ii1 = [by]*4
        bz_i0i, bz_i1i, bz_0ii, bz_1ii = [bz]*4
        bx_i0, bx_i1 = [bx]*2
        by_0i, by_1i = [by]*2
        bz_i0, bz_i1, bz_0i, bz_1i = [bz]*4
        by_0, by_1 = [by]*2
        bz_0, bz_1 = [bz]*2
        if self.n_dim==3:
            if type(bx)==np.ndarray:
                bx_ii0 = bx[:, :, 0]
                bx_ii1 = bx[:, :, 1]
                bx_i0i = bx[:, 0, :]
                bx_i1i = bx[:, 1, :]

            if type(by)==np.ndarray:
                by_0ii = by[0, :, :]
                by_1ii = by[1, :, :]
                by_ii0 = by[:, :, 0]
                by_ii1 = by[:, :, 1]

            if type(bz)==np.ndarray:
                bz_i0i = bz[:, 0, :]
                bz_i1i = bz[:, 1, :]
                bz_0ii = bz[0, :, :]
                bz_1ii = bz[1, :, :]

        elif self.n_dim==2:
            if type(bx)==np.ndarray:
                bx_i0 = bx[:, 0]
                bx_i1 = bx[:, 1]

            if type(by)==np.ndarray:
                by_0i = by[0, :]
                by_1i = by[1, :]

            if type(bz)==np.ndarray:
                bz_i0 = bz[:, 0]
                bz_i1 = bz[:, 1]
                bz_0i = bz[0, :]
                bz_1i = bz[1, :]
        elif self.n_dim==1:
            if type(by)==np.ndarray:
                by_0 = by[0]
                by_1 = by[1]

            if type(bz)==np.ndarray:
                bz_0 = bz[0]
                bz_1 = bz[1]

        if self.n_dim==3:
            Ax[:, :,  0] -= ( By[:, :, 0       ])*(bx_ii0/dz) 
            Ax[:, :, -1] -= (-By[:, :, Nfy[2]-1])*(bx_ii1/dz) 
            Ax[:, 0 , :] += ( Bz[:, 0,        :])*(bx_i0i/dy) 
            Ax[:,-1 , :] += (-Bz[:, Nfz[1]-1, :])*(bx_i1i/dy) 


            Ay[ 0, :, :] -= ( Bz[0,        :, :])*(by_0ii/dx) 
            Ay[-1, :, :] -= (-Bz[Nfz[0]-1, :, :])*(by_1ii/dx) 
            Ay[:, :,  0] += ( Bx[:, :,        0])*(by_ii0/dz) 
            Ay[:, :, -1] += (-Bx[:, :, Nfx[2]-1])*(by_ii1/dz) 

            Az[:, 0 , :] -= ( Bx[:, 0,        :])*(bz_i0i/dy) 
            Az[:, -1, :] -= (-Bx[:, Nfx[1]-1, :])*(bz_i1i/dy) 
            Az[0 , :, :] += ( By[0,        :, :])*(bz_0ii/dx) 
            Az[-1, :, :] += (-By[Nfy[0]-1, :, :])*(bz_1ii/dx) 
        elif self.n_dim==2:
            #Ax[:, 1:-1] -= (By[:, 1:-1    ] - By[:, 1:-1      ])*(bx_*dt/(dz*ax_)) 
            Ax[:, 0] += ( Bz[:, 0       ])*(bx_i0/dy) 
            Ax[:,-1] += (-Bz[:, Nfz[1]-1])*(bx_i1/dy) 

            Ay[0 , :] -= ( Bz[0,        :])*(by_0i/dx) 
            Ay[-1, :] -= (-Bz[Nfz[0]-1, :])*(by_1i/dx) 
            #Ay[1:-1, :] += (Bx[1:-1, :    ] - Bx[1:-1, :      ])*(by_*dt/(dz*ay_)) 

            Az[:,  0] -= ( Bx[:, 0       ])*(bz_i0/dy) 
            Az[:, -1] -= (-Bx[:, Nfx[1]-1])*(bz_i1/dy) 
            Az[0 , :] += ( By[0,        :])*(bz_0i/dx) 
            Az[-1, :] += (-By[Nfy[0]-1, :])*(bz_1i/dx) 
        elif self.n_dim==1:
            #Ax[:] -= (By[:] - By[:])*(bx_*dt/(dz*ax_)) 
            #Ax[:] += (Bz[:] - Bz[:])*(bx_*dt/(dy*ax_)) 

            Ay[0 ] -= ( Bz[0       ])*(by_0/dx) 
            Ay[-1] -= (-Bz[Nfz[0]-1])*(by_1/dx) 
            #Ay[1:-1] += (Bx[1:-1    ] - Bx[1:-1      ])*(by_*dt/(dz*ay_)) 

            #Az[1:-1] -= (Bx[1:-1    ] - Bx[1:-1      ])*(bz_*dt/(dy*az_)) 
            Az[0 ] += ( By[0       ])*(bz_0/dx) 
            Az[-1] += (-By[Nfy[0]-1])*(bz_1/dx) 

        return [Ax, Ay, Az]
        

    ##TODO: corners will be correct?
    def sumCurlFaceWalls(self, A_str):
        """ A_str is the curl and is assumed to have the same name among all processes
        for example self.A
        """
        if self.ChunkInd is None:
            return
        elif np.all(self.ChunkInd==self.ChunkIndTot):
            return
        Nsx, Nsy, Nsz = self.Nsx, self.Nsy, self.Nsz
        ChunkIndNB = self.ChunkIndNB
        A = [None]*3
        exec('A[0], A[1], A[2] = %s'%(A_str))
        Ax, Ay, Az = A
        if self.n_dim==3:
            run_later = []
            rl_vars = {}    ## rl_ : run later
            ##------  X 
            ind_nb_xp = ChunkIndNB['p'][0]
            ind_nb_xn = ChunkIndNB['n'][0]
            if ind_nb_xp is None:
                ##now all chunks respond to same messages sent by neighbors on the opposite side
                if ind_nb_xn != None:
                    self.processNeighbArrSliceRequest(ind_nb_xn)
            else:
                tag = self.getUniqueTag()
                self.sendGetArrSliceCommand(ind_nb_xp, tag, A_str_list=[A_str+'[1][0,:,:]', A_str+'[2][0,:,:]'])
                ##now all chunks respond to same messages sent by neighbors on the opposite side
                if ind_nb_xn != None:
                    self.processNeighbArrSliceRequest(ind_nb_xn)
                Anby_0ii, Anbz_0ii = self.GetArrSliceCommandResp(ind_nb_xp, tag)
                #Ay[Nsy[0]-1, :       , :       ] += Anby_0ii
                #Az[Nsz[0]-1, :       , :       ] += Anbz_0ii
                rl_vars["Anby_0ii"] = Anby_0ii
                rl_vars["Anbz_0ii"] = Anbz_0ii
                run_later.append('Ay[Nsy[0]-1, :       , :       ] += rl_vars["Anby_0ii"]')
                run_later.append('Az[Nsz[0]-1, :       , :       ] += rl_vars["Anbz_0ii"]')
            if ind_nb_xn is None:
                if ind_nb_xp != None:
                    self.processNeighbArrSliceRequest(ind_nb_xp)
            else:
                tag = self.getUniqueTag()
                self.sendGetArrSliceCommand(ind_nb_xn, tag, A_str_list=[A_str+'[1][self.Nsy[0]-1,:,:]', A_str+'[2][self.Nsz[0]-1,:,:]'])
                if ind_nb_xp != None:
                    self.processNeighbArrSliceRequest(ind_nb_xp)
                Anby_1ii, Anbz_1ii = self.GetArrSliceCommandResp(ind_nb_xn, tag)
                #Ay[0       , :       , :       ] += Anby_1ii
                #Az[0       , :       , :       ] += Anbz_1ii
                rl_vars["Anby_1ii"] = Anby_1ii
                rl_vars["Anbz_1ii"] = Anbz_1ii
                run_later.append('Ay[0       , :       , :       ] += rl_vars["Anby_1ii"]')
                run_later.append('Az[0       , :       , :       ] += rl_vars["Anbz_1ii"]')
            ##----- Y 
            ind_nb_yp = ChunkIndNB['p'][1]
            ind_nb_yn = ChunkIndNB['n'][1]
            if ind_nb_yp is None:
                if ind_nb_yn != None:
                    self.processNeighbArrSliceRequest(ind_nb_yn)
            else:
                tag = self.getUniqueTag()
                self.sendGetArrSliceCommand(ind_nb_yp, tag, A_str_list=[A_str+'[0][:,0,:]', A_str+'[2][:,0,:]'])
                if ind_nb_yn != None:
                    self.processNeighbArrSliceRequest(ind_nb_yn)
                Anbx_i0i, Anbz_i0i = self.GetArrSliceCommandResp(ind_nb_yp, tag)
                #Ax[:       , Nsx[1]-1, :       ] += Anbx_i0i
                #Az[:       , Nsz[1]-1, :       ] += Anbz_i0i
                rl_vars["Anbx_i0i"] = Anbx_i0i
                rl_vars["Anbz_i0i"] = Anbz_i0i
                run_later.append('Ax[:       , Nsx[1]-1, :       ] += rl_vars["Anbx_i0i"]')
                run_later.append('Az[:       , Nsz[1]-1, :       ] += rl_vars["Anbz_i0i"]')
            if ind_nb_yn is None:
                if ind_nb_yp != None:
                    self.processNeighbArrSliceRequest(ind_nb_yp)
            else:
                tag = self.getUniqueTag()
                self.sendGetArrSliceCommand(ind_nb_yn, tag, A_str_list=[A_str+'[0][:,self.Nsx[1]-1,:]', A_str+'[2][:,self.Nsz[1]-1,:]'])
                if ind_nb_yp != None:
                    self.processNeighbArrSliceRequest(ind_nb_yp)
                Anbx_i1i, Anbz_i1i = self.GetArrSliceCommandResp(ind_nb_yn, tag)
                #Ax[:       , 0        , :       ] += Anbx_i1i
                #Az[:       , 0        , :       ] += Anbz_i1i
                rl_vars["Anbx_i1i"] = Anbx_i1i
                rl_vars["Anbz_i1i"] = Anbz_i1i
                run_later.append('Ax[:       , 0        , :       ] += rl_vars["Anbx_i1i"]')
                run_later.append('Az[:       , 0        , :       ] += rl_vars["Anbz_i1i"]')
            ##----- Z 
            ind_nb_zp = ChunkIndNB['p'][2]
            ind_nb_zn = ChunkIndNB['n'][2]
            if ind_nb_zp is None:
                if ind_nb_zn != None:
                    self.processNeighbArrSliceRequest(ind_nb_zn)
            else:
                tag = self.getUniqueTag()
                self.sendGetArrSliceCommand(ind_nb_zp, tag, A_str_list=[A_str+'[0][:,:,0]', A_str+'[1][:,:,0]'])
                if ind_nb_zn != None:
                    self.processNeighbArrSliceRequest(ind_nb_zn)
                Anbx_ii0, Anby_ii0 = self.GetArrSliceCommandResp(ind_nb_zp, tag)
                #Ax[:       , :       , Nsx[2]-1] += Anbx_ii0
                #Ay[:       , :       , Nsy[2]-1] += Anby_ii0
                rl_vars["Anbx_ii0"] = Anbx_ii0
                rl_vars["Anby_ii0"] = Anby_ii0
                run_later.append('Ax[:       , :       , Nsx[2]-1] += rl_vars["Anbx_ii0"]')
                run_later.append('Ay[:       , :       , Nsy[2]-1] += rl_vars["Anby_ii0"]')
            if ind_nb_zn is None:
                if ind_nb_zp != None:
                    self.processNeighbArrSliceRequest(ind_nb_zp)
            else:
                tag = self.getUniqueTag()
                self.sendGetArrSliceCommand(ind_nb_zn, tag, A_str_list=[A_str+'[0][:,:,self.Nsx[2]-1]', A_str+'[1][:,:,self.Nsy[2]-1]'])
                if ind_nb_zp != None:
                    self.processNeighbArrSliceRequest(ind_nb_zp)
                Anbx_ii1, Anby_ii1 = self.GetArrSliceCommandResp(ind_nb_zn, tag)
                #Ax[:       , :       , 0       ] += Anbx_ii1
                #Ay[:       , :       , 0       ] += Anby_ii1
                rl_vars["Anbx_ii1"] = Anbx_ii1
                rl_vars["Anby_ii1"] = Anby_ii1
                run_later.append('Ax[:       , :       , 0       ] += rl_vars["Anbx_ii1"]')
                run_later.append('Ay[:       , :       , 0       ] += rl_vars["Anby_ii1"]')
                
            for expr in run_later:
                exec(expr)

        elif self.n_dim==2:
            run_later = []
            rl_vars = {}    ## rl_ : run later
            ##------  X 
            ind_nb_xp = ChunkIndNB['p'][0]
            ind_nb_xn = ChunkIndNB['n'][0]
            if ind_nb_xp is None:
                ##now all chunks respond to same messages sent by neighbors on the opposite side
                if ind_nb_xn != None:
                    self.processNeighbArrSliceRequest(ind_nb_xn)
            else:
                tag = self.getUniqueTag()
                self.sendGetArrSliceCommand(ind_nb_xp, tag, A_str_list=[A_str+'[1][0,:]', A_str+'[2][0,:]'])
                ##now all chunks respond to same messages sent by neighbors on the opposite side
                if ind_nb_xn != None:
                    self.processNeighbArrSliceRequest(ind_nb_xn)
                Anby_0i, Anbz_0i = self.GetArrSliceCommandResp(ind_nb_xp, tag)
                #Ay[Nsy[0]-1, :       ] += Anby_0i
                #Az[Nsz[0]-1, :       ] += Anbz_0i
                rl_vars["Anby_0i"] = Anby_0i
                rl_vars["Anbz_0i"] = Anbz_0i
                run_later.append('Ay[Nsy[0]-1, :       ] += rl_vars["Anby_0i"]')
                run_later.append('Az[Nsz[0]-1, :       ] += rl_vars["Anbz_0i"]')
            if ind_nb_xn is None:
                if ind_nb_xp != None:
                    self.processNeighbArrSliceRequest(ind_nb_xp)
            else:
                tag = self.getUniqueTag()
                self.sendGetArrSliceCommand(ind_nb_xn, tag, A_str_list=[A_str+'[1][self.Nsy[0]-1,:]', A_str+'[2][self.Nsz[0]-1,:]'])
                if ind_nb_xp != None:
                    self.processNeighbArrSliceRequest(ind_nb_xp)
                Anby_1i, Anbz_1i = self.GetArrSliceCommandResp(ind_nb_xn, tag)
                #Ay[0       , :       ] += Anby_1i
                #Az[0       , :       ] += Anbz_1i
                rl_vars["Anby_1i"] = Anby_1i
                rl_vars["Anbz_1i"] = Anbz_1i
                run_later.append('Ay[0       , :       ] += rl_vars["Anby_1i"]')
                run_later.append('Az[0       , :       ] += rl_vars["Anbz_1i"]')
            ##----- Y 
            ind_nb_yp = ChunkIndNB['p'][1]
            ind_nb_yn = ChunkIndNB['n'][1]
            if ind_nb_yp is None:
                if ind_nb_yn != None:
                    self.processNeighbArrSliceRequest(ind_nb_yn)
            else:
                tag = self.getUniqueTag()
                self.sendGetArrSliceCommand(ind_nb_yp, tag, A_str_list=[A_str+'[0][:,0]', A_str+'[2][:,0]'])
                if ind_nb_yn != None:
                    self.processNeighbArrSliceRequest(ind_nb_yn)
                Anbx_i0, Anbz_i0 = self.GetArrSliceCommandResp(ind_nb_yp, tag)
                #Ax[:       , Nsx[1]-1] += Anbx_i0
                #Az[:       , Nsz[1]-1] += Anbz_i0
                rl_vars["Anbx_i0"] = Anbx_i0
                rl_vars["Anbz_i0"] = Anbz_i0
                run_later.append('Ax[:       , Nsx[1]-1] += rl_vars["Anbx_i0"]')
                run_later.append('Az[:       , Nsz[1]-1] += rl_vars["Anbz_i0"]')
            if ind_nb_yn is None:
                if ind_nb_yp != None:
                    self.processNeighbArrSliceRequest(ind_nb_yp)
            else:
                tag = self.getUniqueTag()
                self.sendGetArrSliceCommand(ind_nb_yn, tag, A_str_list=[A_str+'[0][:,self.Nsx[1]-1]', A_str+'[2][:,self.Nsz[1]-1]'])
                if ind_nb_yp != None:
                    self.processNeighbArrSliceRequest(ind_nb_yp)
                Anbx_i1, Anbz_i1 = self.GetArrSliceCommandResp(ind_nb_yn, tag)
                #Ax[:       , 0        ] += Anbx_i1
                #Az[:       , 0        ] += Anbz_i1
                rl_vars["Anbx_i1"] = Anbx_i1
                rl_vars["Anbz_i1"] = Anbz_i1
                run_later.append('Ax[:       , 0        ] += rl_vars["Anbx_i1"]')
                run_later.append('Az[:       , 0        ] += rl_vars["Anbz_i1"]')
                
            for expr in run_later:
                exec(expr)
                
        elif self.n_dim==1:
            run_later = []
            rl_vars = {}    ## rl_ : run later
            ##------  X 
            ind_nb_xp = ChunkIndNB['p'][0]
            ind_nb_xn = ChunkIndNB['n'][0]
            if ind_nb_xp is None:
                ##now all chunks respond to same messages sent by neighbors on the opposite side
                if ind_nb_xn != None:
                    self.processNeighbArrSliceRequest(ind_nb_xn)
            else:
                tag = self.getUniqueTag()
                self.sendGetArrSliceCommand(ind_nb_xp, tag, A_str_list=[A_str+'[1][0]', A_str+'[2][0]'])
                ##now all chunks respond to same messages sent by neighbors on the opposite side
                if ind_nb_xn != None:
                    self.processNeighbArrSliceRequest(ind_nb_xn)
                Anby_0, Anbz_0 = self.GetArrSliceCommandResp(ind_nb_xp, tag)
                #Ay[Nsy[0]-1] += Anby_0
                #Az[Nsz[0]-1] += Anbz_0
                rl_vars["Anby_0"] = Anby_0
                rl_vars["Anbz_0"] = Anbz_0
                run_later.append('Ay[Nsy[0]-1] += rl_vars["Anby_0"]')
                run_later.append('Az[Nsz[0]-1] += rl_vars["Anbz_0"]')
            if ind_nb_xn is None:
                if ind_nb_xp != None:
                    self.processNeighbArrSliceRequest(ind_nb_xp)
            else:
                tag = self.getUniqueTag()
                self.sendGetArrSliceCommand(ind_nb_xn, tag, A_str_list=[A_str+'[1][self.Nsy[0]-1]', A_str+'[2][self.Nsz[0]-1]'])
                if ind_nb_xp != None:
                    self.processNeighbArrSliceRequest(ind_nb_xp)
                Anby_1, Anbz_1 = self.GetArrSliceCommandResp(ind_nb_xn, tag)
                #Ay[0       ] += Anby_1
                #Az[0       ] += Anbz_1
                rl_vars["Anby_1"] = Anby_1
                rl_vars["Anbz_1"] = Anbz_1
                run_later.append('Ay[0       ] += rl_vars["Anby_1"]')
                run_later.append('Az[0       ] += rl_vars["Anbz_1"]')
                
            for expr in run_later:
                exec(expr)
    
    def getUniqueTag(self):
        self.tag += 1
        return self.tag
        
    def processNeighbArrSliceRequest(self, ind_nb):
        req_nb = self.commPipesNBDicRecv[tuple(ind_nb)].recv()
        assert req_nb[0]==FdtdCmds.getArrSlice
        tag_nb = req_nb[1]
        arr_slice_inds = req_nb[2]
        arr_nb = [None]*len(arr_slice_inds)
        for i in range(len(arr_slice_inds)):
            exec('arr_nb[i] = '+arr_slice_inds[i])
        req_nb_resp = [FdtdCmds.getArrSliceResp, tag_nb, arr_nb]
        self.commPipesNBDicRecv[tuple(ind_nb)].send(req_nb_resp)
        
    
    def sendGetArrSliceCommand(self, ind_nb, tag, A_str_list):
        cmd = [FdtdCmds.getArrSlice , tag, A_str_list]
        self.commPipesNBDicSend[tuple(ind_nb)].send(cmd)
            

    def GetArrSliceCommandResp(self, ind_nb, tag):
        resp = self.commPipesNBDicSend[tuple(ind_nb)].recv()
        assert resp[0]==FdtdCmds.getArrSliceResp
        assert resp[1]==tag
        return resp[2]
        
        
    def UpdateFace_dAdt_CurlB_C(self, A, B, C_list=None, a=None, b=None, c_list=None):
        """ Wall indices for A are not updated, it will be later forced to 0
            a*d/dt A = b*curl(B) + c*C
            a, A, b, c, C are face elems
            B is side elem
            
            Should be replaced with the GetCurlSide in the parallel version
        """
        if a is None:
            a = np.ones(3)
        if b is None:
            b = np.ones(3)

        Nsx, Nsy, Nsz = self.Nsx, self.Nsy, self.Nsz
        ax, ay, az = a
        bx, by, bz = b

        Ax, Ay, Az = A
        Bx, By, Bz = B
        dx, dy, dz = [None]*3
        if self.n_dim==3:
            dx, dy, dz = self.dr[0], self.dr[1], self.dr[2]
        elif self.n_dim==2:
            dx, dy, dz = self.dr[0], self.dr[1], None
        elif self.n_dim==1:
            dx, dy, dz = self.dr[0], None, None
        dt = self.dt

        if C_list!=None:
            if c_list is None:
                c_list = [None]*len(C_list)
            for i in range(len(C_list)):
                C = C_list[i]
                c = c_list[i]
                if c is None:
                    c = np.ones(3)
                    
                Cx, Cy, Cz = C
                cx, cy, cz = c
                Ax += Cx*(cx*dt/ax)
                Ay += Cy*(cy*dt/ay)
                Az += Cz*(cz*dt/az)

        if self.n_dim==3:
            Ax -= (By[:, :, 1:Nsy[2]] - By[:, :, 0:Nsy[2]-1])*(bx*dt/(dz*ax)) 
            Ax += (Bz[:, 1:Nsz[1], :] - Bz[:, 0:Nsz[1]-1, :])*(bx*dt/(dy*ax)) 

            Ay -= (Bz[1:Nsz[0], :, :] - Bz[0:Nsz[0]-1, :, :])*(by*dt/(dx*ay)) 
            Ay += (Bx[:, :, 1:Nsx[2]] - Bx[:, :, 0:Nsx[2]-1])*(by*dt/(dz*ay)) 

            Az -= (Bx[:, 1:Nsx[1], :] - Bx[:, 0:Nsx[1]-1, :])*(bz*dt/(dy*az)) 
            Az += (By[1:Nsy[0], :, :] - By[0:Nsy[0]-1, :, :])*(bz*dt/(dx*az)) 
        elif self.n_dim==2:
            #Ax -= (By[:, :       ] - By[:, :         ])*(bx*dt/(dz*ax)) 
            Ax += (Bz[:, 1:Nsz[1]] - Bz[:, 0:Nsz[1]-1])*(bx*dt/(dy*ax)) 

            Ay -= (Bz[1:Nsz[0], :] - Bz[0:Nsz[0]-1, :])*(by*dt/(dx*ay)) 
            #Ay += (Bx[:, :       ] - Bx[:, :         ])*(by*dt/(dz*ay)) 

            Az -= (Bx[:, 1:Nsx[1]] - Bx[:, 0:Nsx[1]-1])*(bz*dt/(dy*az)) 
            Az += (By[1:Nsy[0], :] - By[0:Nsy[0]-1, :])*(bz*dt/(dx*az)) 
        elif self.n_dim==1:
            #Ax -= (By[:] - By[:])*(bx*dt/(dz*ax)) 
            #Ax += (Bz[:] - Bz[:])*(bx*dt/(dy*ax)) 

            Ay -= (Bz[1:Nsz[0]] - Bz[0:Nsz[0]-1])*(by*dt/(dx*ay)) 
            #Ay += (Bx[:       ] - Bx[:         ])*(by*dt/(dz*ay)) 

            #Az -= (Bx[:       ] - Bx[:         ])*(bz*dt/(dy*az)) 
            Az += (By[1:Nsy[0]] - By[0:Nsy[0]-1])*(bz*dt/(dx*az)) 


    def GetCurlSide(self, B, b=None):
        """ b*curl(B) B is side element
            returns face element
        """
        if b is None:
            b = np.ones(3)
            
        Nsx, Nsy, Nsz = self.Nsx, self.Nsy, self.Nsz
        Bx, By, Bz = B
        bx, by, bz = b

        dx, dy, dz = [None]*3
        if self.n_dim==3:
            dx, dy, dz = self.dr[0], self.dr[1], self.dr[2]
        elif self.n_dim==2:
            dx, dy, dz = self.dr[0], self.dr[1], None
        elif self.n_dim==1:
            dx, dy, dz = self.dr[0], None, None
        dt = self.dt

        Ax, Ay, Az = np.zeros(self.Nfx, dtype=self.dtype), np.zeros(self.Nfy, dtype=self.dtype), \
                        np.zeros(self.Nfz, dtype=self.dtype)
        if self.n_dim==3:
            Ax -= (By[:, :, 1:Nsy[2]] - By[:, :, 0:Nsy[2]-1])*(bx/dz) 
            Ax += (Bz[:, 1:Nsz[1], :] - Bz[:, 0:Nsz[1]-1, :])*(bx/dy) 

            Ay -= (Bz[1:Nsz[0], :, :] - Bz[0:Nsz[0]-1, :, :])*(by/dx) 
            Ay += (Bx[:, :, 1:Nsx[2]] - Bx[:, :, 0:Nsx[2]-1])*(by/dz) 

            Az -= (Bx[:, 1:Nsx[1], :] - Bx[:, 0:Nsx[1]-1, :])*(bz/dy) 
            Az += (By[1:Nsy[0], :, :] - By[0:Nsy[0]-1, :, :])*(bz/dx) 
        if self.n_dim==2:
            #Ax -= (By[:, :       ] - By[:, :         ])*(bx/dz) 
            Ax += (Bz[:, 1:Nsz[1]] - Bz[:, 0:Nsz[1]-1])*(bx/dy) 

            Ay -= (Bz[1:Nsz[0], :] - Bz[0:Nsz[0]-1, :])*(by/dx) 
            #Ay += (Bx[:, :       ] - Bx[:, :         ])*(by/dz) 

            Az -= (Bx[:, 1:Nsx[1]] - Bx[:, 0:Nsx[1]-1])*(bz/dy) 
            Az += (By[1:Nsy[0], :] - By[0:Nsy[0]-1, :])*(bz/dx) 
        if self.n_dim==1:
            #Ax -= (By[:] - By[:])*(bx/dz) 
            #Ax += (Bz[:] - Bz[:])*(bx/dy) 

            Ay -= (Bz[1:Nsz[0]] - Bz[0:Nsz[0]-1])*(by/dx) 
            #Ay += (Bx[:       ] - Bx[:         ])*(by/dz) 

            #Az -= (Bx[:       ] - Bx[:         ])*(bz/dy) 
            Az += (By[1:Nsy[0]] - By[0:Nsy[0]-1])*(bz/dx) 
        return [Ax, Ay, Az]


    def Update_adAdt_bB(self, A, B, C_list=None, a=None, b=None, c_list=None):
        """ a*d/dt A = b*B
            all are side/face elems
        """
        if a is None:
            a = np.ones(3)
        if b is None:
            b = np.ones(3)

        ax, ay, az = a
        bx, by, bz = b

        Ax, Ay, Az = A
        Bx, By, Bz = B
        
        dt = self.dt
        
        if C_list!=None:
            if c_list is None:
                c_list = [None]*len(C_list)
            for i in range(len(C_list)):
                C = C_list[i]
                c = c_list[i]
                if c is None:
                    c = np.ones(3)
                    
                Cx, Cy, Cz = C
                cx, cy, cz = c
                Ax += Cx*(cx*dt/ax)
                Ay += Cy*(cy*dt/ay)
                Az += Cz*(cz*dt/az)

        Ax += Bx*(bx*dt/ax)
        Ay += By*(by*dt/ay)
        Az += Bz*(bz*dt/az)
        

    def Update_aA_bB(self, A, B, C_list=None, a=None, b=None, c_list=None):
        """ a*A = b*B
            all are side/face elems
        """
        if a is None:
            a = np.ones(3)
        if b is None:
            b = np.ones(3)

        ax, ay, az = a
        bx, by, bz = b

        Bx, By, Bz = B
                
        A[0] = Bx*(bx/ax)
        A[1] = By*(by/ay)
        A[2] = Bz*(bz/az)

        Ax, Ay, Az = A

        if C_list!=None:
            if c_list is None:
                c_list = [None]*len(C_list)
            for i in range(len(C_list)):
                C = C_list[i]
                c = c_list[i]
                if c is None:
                    c = np.ones(3)
                    
                Cx, Cy, Cz = C
                cx, cy, cz = c
                Ax += Cx*(cx/ax)
                Ay += Cy*(cy/ay)
                Az += Cz*(cz/az)


    def ResetSideWalls(self, A):
        """ A is side element
        """
        if self.ChunkInd is None or np.all(self.ChunkInd==self.ChunkIndTot):
            Nsx, Nsy, Nsz = self.Nsx, self.Nsy, self.Nsz
            Ax, Ay, Az = A
            
            if self.n_dim==3:
                Ax[:       , 0       , :       ] = 0.0
                Ax[:       , Nsx[1]-1, :       ] = 0.0
                Ax[:       , :       , 0       ] = 0.0
                Ax[:       , :       , Nsx[2]-1] = 0.0
                
                Ay[0       , :       , :       ] = 0.0
                Ay[Nsy[0]-1, :       , :       ] = 0.0
                Ay[:       , :       , 0       ] = 0.0
                Ay[:       , :       , Nsy[2]-1] = 0.0

                Az[0       , :       , :       ] = 0.0
                Az[Nsz[0]-1, :       , :       ] = 0.0
                Az[:       , 0       , :       ] = 0.0
                Az[:       , Nsz[1]-1, :       ] = 0.0
            elif self.n_dim==2:
                Ax[:       , 0       ] = 0.0
                Ax[:       , Nsx[1]-1] = 0.0
                
                Ay[0       , :       ] = 0.0
                Ay[Nsy[0]-1, :       ] = 0.0

                Az[0       , :       ] = 0.0
                Az[Nsz[0]-1, :       ] = 0.0
                Az[:       , 0       ] = 0.0
                Az[:       , Nsz[1]-1] = 0.0
            elif self.n_dim==1:            
                Ay[0       ] = 0.0
                Ay[Nsy[0]-1] = 0.0

                Az[0       ] = 0.0
                Az[Nsz[0]-1] = 0.0
        else:
            ChunkIndNB = self.ChunkIndNB
            Nsx, Nsy, Nsz = self.Nsx, self.Nsy, self.Nsz
            Ax, Ay, Az = A
            if self.n_dim==3:
                ind_nb_xp, ind_nb_yp, ind_nb_zp = ChunkIndNB['p']
                ind_nb_xn, ind_nb_yn, ind_nb_zn = ChunkIndNB['n']
                if ind_nb_yn is None:
                    Ax[:       , 0       , :       ] = 0.0
                if ind_nb_yp is None:
                    Ax[:       , Nsx[1]-1, :       ] = 0.0
                if ind_nb_zn is None:
                    Ax[:       , :       , 0       ] = 0.0
                if ind_nb_zp is None:
                    Ax[:       , :       , Nsx[2]-1] = 0.0
                
                if ind_nb_xn is None:
                    Ay[0       , :       , :       ] = 0.0
                if ind_nb_xp is None:
                    Ay[Nsy[0]-1, :       , :       ] = 0.0
                if ind_nb_zn is None:
                    Ay[:       , :       , 0       ] = 0.0
                if ind_nb_zp is None:
                    Ay[:       , :       , Nsy[2]-1] = 0.0

                if ind_nb_xn is None:
                    Az[0       , :       , :       ] = 0.0
                if ind_nb_xp is None:
                    Az[Nsz[0]-1, :       , :       ] = 0.0
                if ind_nb_yn is None:
                    Az[:       , 0       , :       ] = 0.0
                if ind_nb_yp is None:
                    Az[:       , Nsz[1]-1, :       ] = 0.0
        
            elif self.n_dim==2:
                ind_nb_xp, ind_nb_yp = ChunkIndNB['p']
                ind_nb_xn, ind_nb_yn = ChunkIndNB['n']
                if ind_nb_yn is None:
                    Ax[:       , 0       ] = 0.0
                if ind_nb_yp is None:
                    Ax[:       , Nsx[1]-1] = 0.0
                
                if ind_nb_xn is None:
                    Ay[0       , :       ] = 0.0
                if ind_nb_xp is None:
                    Ay[Nsy[0]-1, :       ] = 0.0

                if ind_nb_xn is None:
                    Az[0       , :       ] = 0.0
                if ind_nb_xp is None:
                    Az[Nsz[0]-1, :       ] = 0.0
                if ind_nb_yn is None:
                    Az[:       , 0       ] = 0.0
                if ind_nb_yp is None:
                    Az[:       , Nsz[1]-1] = 0.0
            elif self.n_dim==1:            
                ind_nb_xp = ChunkIndNB['p'][0]
                ind_nb_xn = ChunkIndNB['n'][0]
                if ind_nb_xn is None:
                    Ay[0       ] = 0.0
                if ind_nb_xp is None:
                    Ay[Nsy[0]-1] = 0.0

                if ind_nb_xn is None:
                    Az[0       ] = 0.0
                if ind_nb_xp is None:
                    Az[Nsz[0]-1] = 0.0


    def ResetSideWalls_inDir(self, A, dirs):
        """ reset at the given direction
            A is side element
            dirs: list of strings: 'x+', 'x-', 'y+', 'y-', 'z+', 'z-'
        """
        if self.ChunkInd is None or np.all(self.ChunkInd==self.ChunkIndTot):
            Nsx, Nsy, Nsz = self.Nsx, self.Nsy, self.Nsz
            Ax, Ay, Az = A
            
            if self.n_dim==3:
                if 'y-' in dirs:
                    Ax[:       , 0       , :       ] = 0.0
                if 'y+' in dirs:
                    Ax[:       , Nsx[1]-1, :       ] = 0.0
                if 'z-' in dirs:
                    Ax[:       , :       , 0       ] = 0.0
                if 'z+' in dirs:
                    Ax[:       , :       , Nsx[2]-1] = 0.0
                
                if 'x-' in dirs:
                    Ay[0       , :       , :       ] = 0.0
                if 'x+' in dirs:
                    Ay[Nsy[0]-1, :       , :       ] = 0.0
                if 'z-' in dirs:
                    Ay[:       , :       , 0       ] = 0.0
                if 'z+' in dirs:
                    Ay[:       , :       , Nsy[2]-1] = 0.0

                if 'x-' in dirs:
                    Az[0       , :       , :       ] = 0.0
                if 'x+' in dirs:
                    Az[Nsz[0]-1, :       , :       ] = 0.0
                if 'y-' in dirs:
                    Az[:       , 0       , :       ] = 0.0
                if 'y+' in dirs:
                    Az[:       , Nsz[1]-1, :       ] = 0.0
            elif self.n_dim==2:
                if 'y-' in dirs:
                    Ax[:       , 0       ] = 0.0
                if 'y+' in dirs:
                    Ax[:       , Nsx[1]-1] = 0.0
                
                if 'x-' in dirs:
                    Ay[0       , :       ] = 0.0
                if 'x+' in dirs:
                    Ay[Nsy[0]-1, :       ] = 0.0

                if 'x-' in dirs:
                    Az[0       , :       ] = 0.0
                if 'x+' in dirs:
                    Az[Nsz[0]-1, :       ] = 0.0
                if 'y-' in dirs:
                    Az[:       , 0       ] = 0.0
                if 'y+' in dirs:
                    Az[:       , Nsz[1]-1] = 0.0
            elif self.n_dim==1:            
                if 'x-' in dirs:
                    Ay[0       ] = 0.0
                if 'x+' in dirs:
                    Ay[Nsy[0]-1] = 0.0

                if 'x-' in dirs:
                    Az[0       ] = 0.0
                if 'x+' in dirs:
                    Az[Nsz[0]-1] = 0.0
        else:
            ChunkIndNB = self.ChunkIndNB
            Nsx, Nsy, Nsz = self.Nsx, self.Nsy, self.Nsz
            Ax, Ay, Az = A
            if self.n_dim==3:
                ind_nb_xp, ind_nb_yp, ind_nb_zp = ChunkIndNB['p']
                ind_nb_xn, ind_nb_yn, ind_nb_zn = ChunkIndNB['n']
                if ind_nb_yn is None:
                    if 'y-' in dirs:
                        Ax[:       , 0       , :       ] = 0.0
                if ind_nb_yp is None:
                    if 'y+' in dirs:
                        Ax[:       , Nsx[1]-1, :       ] = 0.0
                if ind_nb_zn is None:
                    if 'z-' in dirs:
                        Ax[:       , :       , 0       ] = 0.0
                if ind_nb_zp is None:
                    if 'z+' in dirs:
                        Ax[:       , :       , Nsx[2]-1] = 0.0
                
                if ind_nb_xn is None:
                    if 'x-' in dirs:
                        Ay[0       , :       , :       ] = 0.0
                if ind_nb_xp is None:
                    if 'x+' in dirs:
                        Ay[Nsy[0]-1, :       , :       ] = 0.0
                if ind_nb_zn is None:
                    if 'z-' in dirs:
                        Ay[:       , :       , 0       ] = 0.0
                if ind_nb_zp is None:
                    if 'z+' in dirs:
                        Ay[:       , :       , Nsy[2]-1] = 0.0

                if ind_nb_xn is None:
                    if 'x-' in dirs:
                        Az[0       , :       , :       ] = 0.0
                if ind_nb_xp is None:
                    if 'x+' in dirs:
                        Az[Nsz[0]-1, :       , :       ] = 0.0
                if ind_nb_yn is None:
                    if 'y-' in dirs:
                        Az[:       , 0       , :       ] = 0.0
                if ind_nb_yp is None:
                    if 'y+' in dirs:
                        Az[:       , Nsz[1]-1, :       ] = 0.0
        
            elif self.n_dim==2:
                ind_nb_xp, ind_nb_yp = ChunkIndNB['p']
                ind_nb_xn, ind_nb_yn = ChunkIndNB['n']
                if ind_nb_yn is None:
                    if 'y-' in dirs:
                        Ax[:       , 0       ] = 0.0
                if ind_nb_yp is None:
                    if 'y+' in dirs:
                        Ax[:       , Nsx[1]-1] = 0.0
                
                if ind_nb_xn is None:
                    if 'x-' in dirs:
                        Ay[0       , :       ] = 0.0
                if ind_nb_xp is None:
                    if 'x+' in dirs:
                        Ay[Nsy[0]-1, :       ] = 0.0

                if ind_nb_xn is None:
                    if 'x-' in dirs:
                        Az[0       , :       ] = 0.0
                if ind_nb_xp is None:
                    if 'x+' in dirs:
                        Az[Nsz[0]-1, :       ] = 0.0
                if ind_nb_yn is None:
                    if 'y-' in dirs:
                        Az[:       , 0       ] = 0.0
                if ind_nb_yp is None:
                    if 'y+' in dirs:
                        Az[:       , Nsz[1]-1] = 0.0
            elif self.n_dim==1:            
                ind_nb_xp = ChunkIndNB['p'][0]
                ind_nb_xn = ChunkIndNB['n'][0]
                if ind_nb_xn is None:
                    if 'x-' in dirs:
                        Ay[0       ] = 0.0
                if ind_nb_xp is None:
                    if 'x+' in dirs:
                        Ay[Nsy[0]-1] = 0.0

                if ind_nb_xn is None:
                    if 'x-' in dirs:
                        Az[0       ] = 0.0
                if ind_nb_xp is None:
                    if 'x+' in dirs:
                        Az[Nsz[0]-1] = 0.0

    def ResetPBC(self, A):
        """ zeros the last side element, as if it does not exist
        """
        if self.ChunkInd is None or np.all(self.ChunkInd==self.ChunkIndTot):
            Nsx, Nsy, Nsz = self.Nsx, self.Nsy, self.Nsz
            Ax, Ay, Az = A
            
            if not isinstance(Ax, np.ndarray):
                assert Ax == 0.0
            if not isinstance(Ay, np.ndarray):
                assert Ay == 0.0
            if not isinstance(Az, np.ndarray):
                assert Az == 0.0
            
            if self.n_dim==3:
                if isinstance(Ax, np.ndarray):
                    Ax[:       , Nsx[1]-1, :       ] *= 0.0
                    Ax[:       , :       , Nsx[2]-1] *= 0.0
                
                if isinstance(Ay, np.ndarray):
                    Ay[Nsy[0]-1, :       , :       ] *= 0.0
                    Ay[:       , :       , Nsy[2]-1] *= 0.0

                if isinstance(Az, np.ndarray):
                    Az[Nsz[0]-1, :       , :       ] *= 0.0
                    Az[:       , Nsz[1]-1, :       ] *= 0.0
            elif self.n_dim==2:
                if isinstance(Ax, np.ndarray):
                    Ax[:       , Nsx[1]-1] *= 0.0
                
                if isinstance(Ay, np.ndarray):
                    Ay[Nsy[0]-1, :       ] *= 0.0

                if isinstance(Az, np.ndarray):
                    Az[Nsz[0]-1, :       ] *= 0.0
                    Az[:       , Nsz[1]-1] *= 0.0
            elif self.n_dim==1:            
                if isinstance(Ay, np.ndarray):
                    Ay[Nsy[0]-1] *= 0.0

                if isinstance(Az, np.ndarray):
                    Az[Nsz[0]-1] *= 0.0
        else:
            ##TODO
            raise NotImplementedError()


    def ApplyPBC(self, A):
        """ apply periodic boundary conditions in all directions
            A is a side element
            ASSUMPTION: border elements are updated one sidedly
        """
        if self.ChunkInd is None or np.all(self.ChunkInd==self.ChunkIndTot):
            Nsx, Nsy, Nsz = self.Nsx, self.Nsy, self.Nsz
            Ax, Ay, Az = A
            
            exp_mjkr = np.exp(-1j*self.W*self.k_vec)
            exp_pjkr = np.exp(+1j*self.W*self.k_vec)
            
            if self.n_dim==3:
                Ax[:       , 0       , :       ] += Ax[:       , Nsx[1]-1, :       ]*exp_pjkr[1]
                Ax[:       , Nsx[1]-1, :       ] =  Ax[:       , 0       , :       ]*exp_mjkr[1]
                Ax[:       , :       , 0       ] += Ax[:       , :       , Nsx[2]-1]*exp_pjkr[2]
                Ax[:       , :       , Nsx[2]-1] =  Ax[:       , :       , 0       ]*exp_mjkr[2]
                
                Ay[0       , :       , :       ] += Ay[Nsy[0]-1, :       , :       ]*exp_pjkr[0]
                Ay[Nsy[0]-1, :       , :       ] =  Ay[0       , :       , :       ]*exp_mjkr[0]
                Ay[:       , :       , 0       ] += Ay[:       , :       , Nsy[2]-1]*exp_pjkr[2]
                Ay[:       , :       , Nsy[2]-1] =  Ay[:       , :       , 0       ]*exp_mjkr[2]

                Az[0       , :       , :       ] += Az[Nsz[0]-1, :       , :       ]*exp_pjkr[0]
                Az[Nsz[0]-1, :       , :       ] =  Az[0       , :       , :       ]*exp_mjkr[0]
                Az[:       , 0       , :       ] += Az[:       , Nsz[1]-1, :       ]*exp_pjkr[1]
                Az[:       , Nsz[1]-1, :       ] =  Az[:       , 0       , :       ]*exp_mjkr[1]
            elif self.n_dim==2:
                Ax[:       , 0       ] += Ax[:       , Nsx[1]-1]*exp_pjkr[1]
                Ax[:       , Nsx[1]-1] =  Ax[:       , 0       ]*exp_mjkr[1]
                
                Ay[0       , :       ] += Ay[Nsy[0]-1, :       ]*exp_pjkr[0]
                Ay[Nsy[0]-1, :       ] =  Ay[0       , :       ]*exp_mjkr[0]

                Az[0       , :       ] += Az[Nsz[0]-1, :       ]*exp_pjkr[0]
                Az[Nsz[0]-1, :       ] =  Az[0       , :       ]*exp_mjkr[0]
                Az[:       , 0       ] += Az[:       , Nsz[1]-1]*exp_pjkr[1]
                Az[:       , Nsz[1]-1] =  Az[:       , 0       ]*exp_mjkr[1]
            elif self.n_dim==1:            
                Ay[0       ] += Ay[Nsy[0]-1]*exp_pjkr[0]
                Ay[Nsy[0]-1] =  Ay[0       ]*exp_mjkr[0]

                Az[0       ] += Az[Nsz[0]-1]*exp_pjkr[0]
                Az[Nsz[0]-1] =  Az[0       ]*exp_mjkr[0]
        else:
            ##TODO
            raise NotImplementedError()
        

    ##TODO:
    def ApplyPBC_inDir(self, A, dirs):
        """ apply periodic boundary conditions in the given directions
            A is a side element
            dirs: list of strings: 'x+', 'x-', 'y+', 'y-', 'z+', 'z-'
            ASSUMPTION: border elements are updated one sidedly
        """
        raise NotImplementedError()
        return


    def SetSpatialPointsSide(self):     
        if self.n_dim==3:   
            dx, dy, dz = self.dr[0], self.dr[1], self.dr[2]
            x0, y0, z0 = self.r0[0], self.r0[1], self.r0[2]
            x1, y1, z1 = self.r1[0], self.r1[1], self.r1[2]
            nx, ny, nz = self.N[0],  self.N[1],  self.N[2]
            
            x_1d = np.linspace(x0+dx/2, x1-dx/2, nx, endpoint=True)
            y_1d = np.linspace(y0+dy/2, y1-dy/2, ny, endpoint=True)
            z_1d = np.linspace(z0+dz/2, z1-dz/2, nz, endpoint=True)

            x_1dp = np.linspace(x0, x1, nx+1, endpoint=True)
            y_1dp = np.linspace(y0, y1, ny+1, endpoint=True)
            z_1dp = np.linspace(z0, z1, nz+1, endpoint=True)
             
            rx_x, rx_y, rx_z = np.meshgrid(x_1d, y_1dp, z_1dp, indexing='ij')
            self.rsx = [rx_x, rx_y, rx_z]

            ry_x, ry_y, ry_z = np.meshgrid(x_1dp, y_1d, z_1dp, indexing='ij')
            self.rsy = [ry_x, ry_y, ry_z]

            rz_x, rz_y, rz_z = np.meshgrid(x_1dp, y_1dp, z_1d, indexing='ij')
            self.rsz = [rz_x, rz_y, rz_z]
        elif self.n_dim==2:   
            dx, dy, dz = self.dr[0], self.dr[1], None
            x0, y0, z0 = self.r0[0], self.r0[1], None
            x1, y1, z1 = self.r1[0], self.r1[1], None
            nx, ny, nz = self.N[0],  self.N[1],  None
            
            x_1d = np.linspace(x0+dx/2, x1-dx/2, nx, endpoint=True)
            y_1d = np.linspace(y0+dy/2, y1-dy/2, ny, endpoint=True)

            x_1dp = np.linspace(x0, x1, nx+1, endpoint=True)
            y_1dp = np.linspace(y0, y1, ny+1, endpoint=True)
             
            rx_x, rx_y = np.meshgrid(x_1d, y_1dp, indexing='ij')
            self.rsx = [rx_x, rx_y]

            ry_x, ry_y = np.meshgrid(x_1dp, y_1d, indexing='ij')
            self.rsy = [ry_x, ry_y]

            rz_x, rz_y = np.meshgrid(x_1dp, y_1dp, indexing='ij')
            self.rsz = [rz_x, rz_y]
        elif self.n_dim==1:   
            dx, dy, dz = self.dr[0], None, None
            x0, y0, z0 = self.r0[0], None, None
            x1, y1, z1 = self.r1[0], None, None
            nx, ny, nz = self.N[0],  None,  None
            
            x_1d = np.linspace(x0+dx/2, x1-dx/2, nx, endpoint=True)

            x_1dp = np.linspace(x0, x1, nx+1, endpoint=True)
             
            rx_x = x_1d
            self.rsx = [rx_x]

            ry_x = x_1dp
            self.rsy = [ry_x]

            rz_x = x_1dp
            self.rsz = [rz_x]


    def SetSpatialPointsFace(self):        
        if self.n_dim==3:
            dx, dy, dz = self.dr[0], self.dr[1], self.dr[2]
            x0, y0, z0 = self.r0[0], self.r0[1], self.r0[2]
            x1, y1, z1 = self.r1[0], self.r1[1], self.r1[2]
            nx, ny, nz = self.N[0],  self.N[1],  self.N[2]
            
            x_1d = np.linspace(x0+dx/2, x1-dx/2, nx, endpoint=True)
            y_1d = np.linspace(y0+dy/2, y1-dy/2, ny, endpoint=True)
            z_1d = np.linspace(z0+dz/2, z1-dz/2, nz, endpoint=True)

            x_1dp = np.linspace(x0, x1, nx+1, endpoint=True)
            y_1dp = np.linspace(y0, y1, ny+1, endpoint=True)
            z_1dp = np.linspace(z0, z1, nz+1, endpoint=True)
            
            rx_x, rx_y, rx_z = np.meshgrid(x_1dp, y_1d, z_1d, indexing='ij')
            self.rfx = [rx_x, rx_y, rx_z]
             
            ry_x, ry_y, ry_z = np.meshgrid(x_1d, y_1dp, z_1d, indexing='ij')
            self.rfy = [ry_x, ry_y, ry_z]

            rz_x, rz_y, rz_z = np.meshgrid(x_1d, y_1d, z_1dp, indexing='ij')
            self.rfz = [rz_x, rz_y, rz_z]
        elif self.n_dim==2:
            dx, dy, dz = self.dr[0], self.dr[1], None
            x0, y0, z0 = self.r0[0], self.r0[1], None
            x1, y1, z1 = self.r1[0], self.r1[1], None
            nx, ny, nz = self.N[0],  self.N[1],  None
            
            x_1d = np.linspace(x0+dx/2, x1-dx/2, nx, endpoint=True)
            y_1d = np.linspace(y0+dy/2, y1-dy/2, ny, endpoint=True)

            x_1dp = np.linspace(x0, x1, nx+1, endpoint=True)
            y_1dp = np.linspace(y0, y1, ny+1, endpoint=True)
            
            rx_x, rx_y = np.meshgrid(x_1dp, y_1d, indexing='ij')
            self.rfx = [rx_x, rx_y]
             
            ry_x, ry_y = np.meshgrid(x_1d, y_1dp, indexing='ij')
            self.rfy = [ry_x, ry_y]

            rz_x, rz_y = np.meshgrid(x_1d, y_1d, indexing='ij')
            self.rfz = [rz_x, rz_y]
        elif self.n_dim==1:
            dx, dy, dz = self.dr[0], None, None
            x0, y0, z0 = self.r0[0], None, None
            x1, y1, z1 = self.r1[0], None, None
            nx, ny, nz = self.N[0],  None,  None
            
            x_1d = np.linspace(x0+dx/2, x1-dx/2, nx, endpoint=True)

            x_1dp = np.linspace(x0, x1, nx+1, endpoint=True)
            
            rx_x = x_1dp
            self.rfx = [rx_x]
             
            ry_x = x_1d
            self.rfy = [ry_x]

            rz_x = x_1d
            self.rfz = [rz_x]
    

    def UpdateFuncSide_space(self, A, f, getCopy=False):
        fx, fy, fz = f
        rsx, rsy, rsz = self.rsx, self.rsy, self.rsz
        if self.n_dim==1:
            rsx, rsy, rsz = self.rsx[0], self.rsy[0], self.rsz[0]
        
        if fx!=None:                
            A[0] = fx(rsx)
        if fy!=None:                
            A[1] = fy(rsy)
        if fz!=None:                
            A[2] = fz(rsz)
        if getCopy:
            A_copy = [np.copy(A[0]), np.copy(A[1]), np.copy(A[2])]
            return A_copy
            

    def UpdateFuncSide_spacetime(self, A, f, t, getCopy=False):
        fx, fy, fz = f
        rsx, rsy, rsz = self.rsx, self.rsy, self.rsz
        if self.n_dim==1:
            rsx, rsy, rsz = self.rsx[0], self.rsy[0], self.rsz[0]
        
        if fx!=None:                
            A[0] = fx(rsx, t)
        if fy!=None:                
            A[1] = fy(rsy, t)
        if fz!=None:                
            A[2] = fz(rsz, t)
        if getCopy:
            A_copy = [np.copy(A[0]), np.copy(A[1]), np.copy(A[2])]
            return A_copy

    
    def UpdateFuncFace_space(self, A, f, getCopy=False):
        fx, fy, fz = f
        rfx, rfy, rfz = self.rfx, self.rfy, self.rfz
        if self.n_dim==1:
            rfx, rfy, rfz = self.rfx[0], self.rfy[0], self.rfz[0]

        if fx!=None:                
            A[0] = fx(rfx)
        if fy!=None:                
            A[1] = fy(rfy)
        if fz!=None:                
            A[2] = fz(rfz)
        if getCopy:
            A_copy = [np.copy(A[0]), np.copy(A[1]), np.copy(A[2])]
            return A_copy
        

    def UpdateFuncFace_spacetime(self, A, f, t, getCopy=False):
        fx, fy, fz = f
        rfx, rfy, rfz = self.rfx, self.rfy, self.rfz
        if self.n_dim==1:
            rfx, rfy, rfz = self.rfx[0], self.rfy[0], self.rfz[0]

        if fx!=None:                
            A[0] = fx(rfx, t)
        if fy!=None:                
            A[1] = fy(rfy, t)
        if fz!=None:                
            A[2] = fz(rfz, t)
        if getCopy:
            A_copy = [np.copy(A[0]), np.copy(A[1]), np.copy(A[2])]
            return A_copy


    def UpdateSeperableFunc_Time(self, A, A_0, f, t):
        A[0] = f(t)*A_0[0]
        A[1] = f(t)*A_0[1]
        A[2] = f(t)*A_0[2]

        
    def ComplexEpsToEpsSigmaE(self, eps_c, omega):
        eps = np.real(eps_c)
        sig = -np.imag(eps_c)*omega
        return [eps, sig]


    def ComplexMuToMuSigmaM(self, mu_c, omega):
        mu = np.real(mu_c)
        sig = -np.imag(mu_c)*omega
        return [eps, sig]


    def GetWallsAllDic__ver_0(self, dw, s):
        if self.n_dim==3:
            dx, dy, dz = dw[0], dw[1], dw[2]
            sx, sy, sz = s
            walls = {('x', 'n'):[sx, dx], ('x', 'p'):[sx, dx],
                     ('y', 'n'):[sy, dy], ('y', 'p'):[sy, dy],
                     ('z', 'n'):[sz, dz], ('z', 'p'):[sz, dz]}
            return walls
        elif self.n_dim==2:
            dx, dy = dw[0], dw[1]
            sx, sy = s
            walls = {('x', 'n'):[sx, dx], ('x', 'p'):[sx, dx],
                     ('y', 'n'):[sy, dy], ('y', 'p'):[sy, dy]}
            return walls
        elif self.n_dim==1:
            dx = dw[0]
            sx = s
            walls = {('x', 'n'):[sx, dx], ('x', 'p'):[sx, dx]}
            return walls


    def GetWallsAllDic(self, dwn, sn, dwp=None, sp=None):
        ##symmetric
        if dwp is None:
            dw = dwn
            s = sn
            if self.n_dim==3:
                dx, dy, dz = dw[0], dw[1], dw[2]
                sx, sy, sz = s
                walls = {}
                if dx>0.0:
                    walls[('x', 'n')] = [sx, dx]
                    walls[('x', 'p')] = [sx, dx]
                if dy>0.0:
                    walls[('y', 'n')] = [sy, dy]
                    walls[('y', 'p')] = [sy, dy]
                if dz>0.0:
                    walls[('z', 'n')] = [sz, dz]
                    walls[('z', 'p')] = [sz, dz]
                return walls
            elif self.n_dim==2:
                dx, dy = dw[0], dw[1]
                sx, sy = s
                walls = {}
                if dx>0.0:
                    walls[('x', 'n')] = [sx, dx]
                    walls[('x', 'p')] = [sx, dx]
                if dy>0.0:
                    walls[('y', 'n')] = [sy, dy]
                    walls[('y', 'p')] = [sy, dy]
                return walls
            elif self.n_dim==1:
                dx = dw[0]
                sx = s          ##TODO: s or s[0]
                walls = {}
                if dx>0.0:
                    walls[('x', 'n')] = [sx, dx]
                    walls[('x', 'p')] = [sx, dx]
                return walls
        ##asymmetric
        else:
            if self.n_dim==3:
                dxp, dyp, dzp = dwp[0], dwp[1], dwp[2]
                sxp, syp, szp = sp
                dxn, dyn, dzn = dwn[0], dwn[1], dwn[2]
                sxn, syn, szn = sn
                walls = {}
                if dxn>0.0:
                    walls[('x', 'n')] = [sxn, dxn]
                if dxp>0.0:
                    walls[('x', 'p')] = [sxp, dxp]
                if dyn>0.0:
                    walls[('y', 'n')] = [syn, dyn]
                if dyp>0.0:
                    walls[('y', 'p')] = [syp, dyp]
                if dzn>0.0:
                    walls[('z', 'n')] = [szn, dzn]
                if dzp>0.0:
                    walls[('z', 'p')] = [szp, dzp]
                return walls
            elif self.n_dim==2:
                dxn, dyn = dwn[0], dwn[1]
                sxn, syn = sn
                dxp, dyp = dwp[0], dwp[1]
                sxp, syp = sp
                walls = {}
                if dxn>0.0:
                    walls[('x', 'n')] = [sxn, dxn]
                if dxp>0.0:
                    walls[('x', 'p')] = [sxp, dxp]
                if dyn>0.0:
                    walls[('y', 'n')] = [syn, dyn]
                if dyp>0.0:
                    walls[('y', 'p')] = [syp, dyp]
                return walls
            elif self.n_dim==1:
                dxn = dwn[0]
                sxn = sn
                dxp = dwp[0]
                sxp = sp
                walls = {}
                if dxn>0.0:
                    walls[('x', 'n')] = [sxn, dxn]
                if dxp>0.0:
                    walls[('x', 'p')] = [sxp, dxp]
                return walls

    def GetUPMLFactor(self, walls, eps_or_s='s', side_or_face='side'):
        """ walls:{('x', 'n'):[sx, dx], ('y', 'p'):[sy, dy]..}
            ('x', 'n'):[sx, dx] --> plane normal to x direction, negative side, 
                stretch factor sx, thickness dx
            eps_or_s='eps'/'s'
            returns boxes with complex upml permittivity/permeability factors
            if eps_or_s=='s' returns where sx/sy/sz are non-zero
        """
        Nx, Ny, Nz = [None]*3
        if side_or_face=='side':
            Nx, Ny, Nz = self.Nsx, self.Nsy, self.Nsz
        else:
            assert side_or_face=='face'
            Nx, Ny, Nz = self.Nfx, self.Nfy, self.Nfz
            
        if eps_or_s=='eps':
            eps_xx = np.ones(Nx)
            eps_yy = np.ones(Ny)
            eps_zz = np.ones(Nz)
            for wall_dir in walls:
                w, n_dir = wall_dir
                f = None
                if w=='x':
                    r1 = np.copy(self.r1)
                    sx, dx = walls[wall_dir]
                    if n_dir=='n':
                        r1[0] = self.r0[0] + dx
                        r0 = np.copy(self.r0)
                        f_xx = self.FunctionBox(r0, r1, a=1.0/sx, b=1.0)    
                        f_yy = self.FunctionBox(r0, r1, a=sx, b=1.0)    
                        f_zz = self.FunctionBox(r0, r1, a=sx, b=1.0)    
                        f = [f_xx, f_yy, f_zz]
                    else:
                        assert n_dir=='p'
                        r0 = np.copy(self.r0)
                        r0[0] = self.r1[0] - dx
                        f_xx = self.FunctionBox(r0, r1, a=1.0/sx, b=1.0)    
                        f_yy = self.FunctionBox(r0, r1, a=sx, b=1.0)    
                        f_zz = self.FunctionBox(r0, r1, a=sx, b=1.0)    
                        f = [f_xx, f_yy, f_zz]
                elif w=='y':
                    assert self.n_dim>=2
                    r1 = np.copy(self.r1)
                    sy, dy = walls[wall_dir]
                    if n_dir=='n':
                        r1[1] = self.r0[1] + dy
                        r0 = np.copy(self.r0)
                        f_xx = self.FunctionBox(r0, r1, a=sy, b=1.0)    
                        f_yy = self.FunctionBox(r0, r1, a=1.0/sy, b=1.0)    
                        f_zz = self.FunctionBox(r0, r1, a=sy, b=1.0)    
                        f = [f_xx, f_yy, f_zz]
                    else:
                        assert n_dir=='p'
                        r0 = np.copy(self.r0)
                        r0[1] = self.r1[1] - dy
                        f_xx = self.FunctionBox(r0, r1, a=sy, b=1.0)    
                        f_yy = self.FunctionBox(r0, r1, a=1.0/sy, b=1.0)    
                        f_zz = self.FunctionBox(r0, r1, a=sy, b=1.0)    
                        f = [f_xx, f_yy, f_zz]
                else:
                    assert w=='z'
                    assert self.n_dim==3
                    r1 = np.copy(self.r1)
                    sz, dz = walls[wall_dir]
                    if n_dir=='n':
                        r1[2] = self.r0[2] + dz
                        r0 = np.copy(self.r0)
                        f_xx = self.FunctionBox(r0, r1, a=sz, b=1.0)    
                        f_yy = self.FunctionBox(r0, r1, a=sz, b=1.0)    
                        f_zz = self.FunctionBox(r0, r1, a=1.0/sz, b=1.0)    
                        f = [f_xx, f_yy, f_zz]
                    else:
                        assert n_dir=='p'
                        r0 = np.copy(self.r0)
                        r0[2] = self.r1[2] - dz
                        f_xx = self.FunctionBox(r0, r1, a=sz, b=1.0)    
                        f_yy = self.FunctionBox(r0, r1, a=sz, b=1.0)    
                        f_zz = self.FunctionBox(r0, r1, a=1.0/sz, b=1.0)    
                        f = [f_xx, f_yy, f_zz]
                
                A = [None]*3
                if side_or_face=='side':
                    self.UpdateFuncSide_space(A, f)
                else:
                    self.UpdateFuncFace_space(A, f)
                eps_xx *= A[0]
                eps_yy *= A[1]
                eps_zz *= A[2]
            return [eps_xx, eps_yy, eps_zz]
        else:
            assert eps_or_s == 's'
            sx_arr, sy_arr, sz_arr = [None]*3, [None]*3, [None]*3
            sx_arr[0] = np.ones(Nx, dtype=complex)   ## s_x for sides x
            sx_arr[1] = np.ones(Ny, dtype=complex)   ## s_x for sides y
            sx_arr[2] = np.ones(Nz, dtype=complex)

            sy_arr[0] = np.ones(Nx, dtype=complex)   ## s_y for sides x
            sy_arr[1] = np.ones(Ny, dtype=complex)   ## s_y for sides y
            sy_arr[2] = np.ones(Nz, dtype=complex)

            sz_arr[0] = np.ones(Nx, dtype=complex)   ## s_z for sides x
            sz_arr[1] = np.ones(Ny, dtype=complex)   ## s_z for sides y
            sz_arr[2] = np.ones(Nz, dtype=complex)

            for wall_dir in walls:
                w, n_dir = wall_dir
                f = None
                if w=='x':
                    r1 = np.copy(self.r1)
                    sx, dx = walls[wall_dir]
                    if n_dir=='n':
                        r1[0] = self.r0[0] + dx
                        r0 = np.copy(self.r0)
                        f = self.FunctionBox(r0, r1, a=sx-1.0, b=0.0)    
                    else:
                        assert n_dir=='p'
                        r0 = np.copy(self.r0)
                        r0[0] = self.r1[0] - dx
                        f = self.FunctionBox(r0, r1, a=sx-1.0, b=0.0)    
                    A = [None]*3
                    if side_or_face=='side':
                        self.UpdateFuncSide_space(A, [f, f, f])
                    else:
                        self.UpdateFuncFace_space(A, [f, f, f])
                    sx_arr[0] += A[0]
                    sx_arr[1] += A[1]
                    sx_arr[2] += A[2]
                elif w=='y':
                    assert self.n_dim>=2
                    r1 = np.copy(self.r1)
                    sy, dy = walls[wall_dir]
                    if n_dir=='n':
                        r1[1] = self.r0[1] + dy
                        r0 = np.copy(self.r0)
                        f = self.FunctionBox(r0, r1, a=sy-1.0, b=0.0)    
                    else:
                        assert n_dir=='p'
                        r0 = np.copy(self.r0)
                        r0[1] = self.r1[1] - dy
                        f = self.FunctionBox(r0, r1, a=sy-1.0, b=0.0)    
                    A = [None]*3
                    if side_or_face=='side':
                        self.UpdateFuncSide_space(A, [f, f, f])
                    else:
                        self.UpdateFuncFace_space(A, [f, f, f])
                    sy_arr[0] += A[0]
                    sy_arr[1] += A[1]
                    sy_arr[2] += A[2]
                else:
                    assert w=='z'
                    assert self.n_dim==3
                    r1 = np.copy(self.r1)
                    sz, dz = walls[wall_dir]
                    if n_dir=='n':
                        r1[2] = self.r0[2] + dz
                        r0 = np.copy(self.r0)
                        f = self.FunctionBox(r0, r1, a=sz-1.0, b=0.0)    
                    else:
                        assert n_dir=='p'
                        r0 = np.copy(self.r0)
                        r0[2] = self.r1[2] - dz
                        f = self.FunctionBox(r0, r1, a=sz-1.0, b=0.0)    
                    A = [None]*3
                    if side_or_face=='side':
                        self.UpdateFuncSide_space(A, [f, f, f])
                    else:
                        self.UpdateFuncFace_space(A, [f, f, f])
                    sz_arr[0] += A[0]
                    sz_arr[1] += A[1]
                    sz_arr[2] += A[2]
            return [sx_arr, sy_arr, sz_arr]
                    


    def FunctionSphere(self, r0, rad, a, b=0.0):
        """ r0: center
            rad: radius
            a: amplitude on shpere
            b: amplitude outside sphere
        """
        if self.n_dim==3:
            x0, y0, z0 = r0[0], r0[1], r0[2]        
            def f(r):
                x, y, z = r
                R2 = (x-x0)**2 + (y-y0)**2 + (z-z0)**2
                isin = (R2 <= rad**2)
                isout = np.logical_not(isin)
                res = None
                if inspect.isfunction(a):
                    res = a(r)*isin
                else:
                    res = a*isin
                if inspect.isfunction(b):
                    res += b(r)*isout
                else:
                    res += b*isout
                return res
            return f
        if self.n_dim==2:
            x0, y0 = r0[0], r0[1]        
            def f(r):
                x, y = r
                R2 = (x-x0)**2 + (y-y0)**2
                isin = (R2 <= rad**2)
                isout = np.logical_not(isin)
                res = None
                if inspect.isfunction(a):
                    res = a(r)*isin
                else:
                    res = a*isin
                if inspect.isfunction(b):
                    res += b(r)*isout
                else:
                    res += b*isout
                return res
            return f
        if self.n_dim==1:
            x0 = r0[0]        
            def f(r):
                x = r   ##TODO: make consistent --> x=r[0]
                R2 = (x-x0)**2
                isin = (R2 <= rad**2)
                isout = np.logical_not(isin)
                res = None
                if inspect.isfunction(a):
                    res = a(r)*isin
                else:
                    res = a*isin
                if inspect.isfunction(b):
                    res += b(r)*isout
                else:
                    res += b*isout
                return res
            return f

    def FunctionBox(self, r0, r1, a, b=0.0, exclude=False):
        """ r0: lowe left corner
            r1: upper right corner
            a: amplitude
            exclude: exclude boundary from interior
        """
        if self.n_dim==3:
            x0, y0, z0 = r0[0], r0[1], r0[2]        
            x1, y1, z1 = r1[0], r1[1], r1[2]        
            assert x0<x1 and y0<y1 and z0<z1
            def f(r, t=0):
                x, y, z = r
                isin = (x>=x0) & (x<=x1) & (y>=y0) & (y<=y1) & (z>=z0) & (z<=z1)
                if exclude:
                    isin = (x>x0) & (x<x1) & (y>y0) & (y<y1) & (z>z0) & (z<z1)
                isout = np.logical_not(isin)
                res = None
                if inspect.isfunction(a):
                    a_args = inspect.getargspec(a)[0]
                    if len(a_args)==2:
                        res = a(r, t)*isin
                    else:
                        assert len(a_args)==1
                        res = a(r)*isin
                else:
                    res = a*isin
                if inspect.isfunction(b):
                    b_args = inspect.getargspec(b)[0]
                    if len(b_args)==2:
                        res += b(r, t)*isout
                    else:
                        assert len(b_args)==1
                        res += b(r)*isout
                else:
                    res += b*isout
                return res
            return f
        elif self.n_dim==2:
            x0, y0 = r0[0], r0[1]        
            x1, y1 = r1[0], r1[1]        
            assert x0<x1 and y0<y1
            def f(r, t=0):
                x, y = r
                isin = (x>=x0) & (x<=x1) & (y>=y0) & (y<=y1)
                if exclude:
                    isin = (x>x0) & (x<x1) & (y>y0) & (y<y1)
                isout = np.logical_not(isin)
                res = None
                if inspect.isfunction(a):
                    a_args = inspect.getargspec(a)[0]
                    if len(a_args)==2:
                        res = a(r, t)*isin
                    else:
                        assert len(a_args)==1
                        res = a(r)*isin
                else:
                    res = a*isin
                if inspect.isfunction(b):
                    b_args = inspect.getargspec(b)[0]
                    if len(b_args)==2:
                        res += b(r, t)*isout
                    else:
                        assert len(b_args)==1
                        res += b(r)*isout
                else:
                    res += b*isout
                return res
            return f
        elif self.n_dim==1:
            x0 = r0[0]        
            x1 = r1[0]        
            assert x0<x1
            def f(r, t=0):
                x = r
                isin = (x>=x0) & (x<=x1)
                if exclude:
                    isin = (x>x0) & (x<x1)
                isout = np.logical_not(isin)
                res = None
                if inspect.isfunction(a):
                    a_args = inspect.getargspec(a)[0]
                    if len(a_args)==2:
                        res = a(r, t)*isin
                    else:
                        assert len(a_args)==1
                        res = a(r)*isin
                else:
                    res = a*isin
                if inspect.isfunction(b):
                    b_args = inspect.getargspec(b)[0]
                    if len(b_args)==2:
                        res += b(r, t)*isout
                    else:
                        assert len(b_args)==1
                        res += b(r)*isout
                else:
                    res += b*isout
                return res
            return f
        
            

    def GetPointSourceFunc(self, r_0, mag=1.0, src_dir='x', em_type='e'):
        ##TODO: do the same for other functions
        if not (np.all(r_0>=self.r0) and np.all(r_0<=self.r1)):
            return [lambda r : 0.0]*3
        if em_type=='e':
            i_xyz_min, r_xyz_min, d_xyz_min = self.FindClosestSides(r_0)
            f_je = None
            if src_dir=='x':
                i_j, r_j, d_j = i_xyz_min[0], r_xyz_min[0], d_xyz_min[0]

                f_jex = self.FunctionSphere(r_j, d_j/2.0, mag)
                f_jey = lambda r : 0.0
                f_jez = lambda r : 0.0
                f_je = [f_jex, f_jey, f_jez]
            elif src_dir=='y':
                i_j, r_j, d_j = i_xyz_min[1], r_xyz_min[1], d_xyz_min[1]

                f_jex = lambda r : 0.0
                f_jey = self.FunctionSphere(r_j, d_j/2.0, mag)
                f_jez = lambda r : 0.0
                f_je = [f_jex, f_jey, f_jez]
            elif src_dir=='z':
                i_j, r_j, d_j = i_xyz_min[2], r_xyz_min[2], d_xyz_min[2]

                f_jex = lambda r : 0.0
                f_jey = lambda r : 0.0
                f_jez = self.FunctionSphere(r_j, d_j/2.0, mag)
                f_je = [f_jex, f_jey, f_jez]
            else:
                raise ValueError()
            
            return f_je
        elif em_type=='m' or em_type=='h':
            i_xyz_min, r_xyz_min, d_xyz_min = self.FindClosestFaces(r_0)
            f_jm = None
            if src_dir=='x':
                i_j, r_j, d_j = i_xyz_min[0], r_xyz_min[0], d_xyz_min[0]

                f_jmx = self.FunctionSphere(r_j, d_j/2.0, mag)
                f_jmy = lambda r : 0.0
                f_jmz = lambda r : 0.0
                f_jm = [f_jmx, f_jmy, f_jmz]
            elif src_dir=='y':
                i_j, r_j, d_j = i_xyz_min[1], r_xyz_min[1], d_xyz_min[1]

                f_jmx = lambda r : 0.0
                f_jmy = self.FunctionSphere(r_j, d_j/2.0, mag)
                f_jmz = lambda r : 0.0
                f_jm = [f_jmx, f_jmy, f_jmz]
            elif src_dir=='z':
                i_j, r_j, d_j = i_xyz_min[2], r_xyz_min[2], d_xyz_min[2]

                f_jmx = lambda r : 0.0
                f_jmy = lambda r : 0.0
                f_jmz = self.FunctionSphere(r_j, d_j/2.0, mag)
                f_jm = [f_jmx, f_jmy, f_jmz]
            else:
                raise ValueError()

            return f_jm
        else:
            raise ValueError()


    def GetSheetSourceFunc(self, r_0, mag=1.0, src_dir='x', norm_dir='x', em_type='e'):
        assert src_dir in ['x', 'y', 'z']
        assert norm_dir in ['x', 'y', 'z']
        assert em_type in ['e', 'm', 'h']
        if em_type=='e':
            i_xyz_min, r_xyz_min, d_xyz_min = self.FindClosestSides(r_0)
            f_je = None
            if src_dir=='x':
                i_j, r_j, d_j = i_xyz_min[0], r_xyz_min[0], d_xyz_min[0]
                d_min = min(self.dr)

                r0 = np.copy(self.r0)
                r1 = np.copy(self.r1)
                if norm_dir=='x':
                    r0[0] = r_j[0]-d_min/3
                    r1[0] = r_j[0]+d_min/3
                elif norm_dir=='y':
                    assert self.n_dim>=2
                    r0[1] = r_j[1]-d_min/3
                    r1[1] = r_j[1]+d_min/3
                elif norm_dir=='z':
                    assert self.n_dim==3
                    r0[2] = r_j[2]-d_min/3
                    r1[2] = r_j[2]+d_min/3
                
                f_jex = self.FunctionBox(r0, r1, mag)
                f_jey = lambda r : 0.0
                f_jez = lambda r : 0.0
                f_je = [f_jex, f_jey, f_jez]
            elif src_dir=='y':
                i_j, r_j, d_j = i_xyz_min[1], r_xyz_min[1], d_xyz_min[1]
                d_min = min(self.dr)

                r0 = np.copy(self.r0)
                r1 = np.copy(self.r1)
                if norm_dir=='x':
                    r0[0] = r_j[0]-d_min/3
                    r1[0] = r_j[0]+d_min/3
                elif norm_dir=='y':
                    assert self.n_dim>=2
                    r0[1] = r_j[1]-d_min/3
                    r1[1] = r_j[1]+d_min/3
                elif norm_dir=='z':
                    assert self.n_dim==3
                    r0[2] = r_j[2]-d_min/3
                    r1[2] = r_j[2]+d_min/3
                

                f_jex = lambda r : 0.0
                f_jey = self.FunctionBox(r0, r1, mag)
                f_jez = lambda r : 0.0
                f_je = [f_jex, f_jey, f_jez]
            elif src_dir=='z':
                i_j, r_j, d_j = i_xyz_min[2], r_xyz_min[2], d_xyz_min[2]
                d_min = min(self.dr)

                r0 = np.copy(self.r0)
                r1 = np.copy(self.r1)
                if norm_dir=='x':
                    r0[0] = r_j[0]-d_min/3
                    r1[0] = r_j[0]+d_min/3
                elif norm_dir=='y':
                    assert self.n_dim>=2
                    r0[1] = r_j[1]-d_min/3
                    r1[1] = r_j[1]+d_min/3
                elif norm_dir=='z':
                    assert self.n_dim==3
                    r0[2] = r_j[2]-d_min/3
                    r1[2] = r_j[2]+d_min/3

                f_jex = lambda r : 0.0
                f_jey = lambda r : 0.0
                f_jez = self.FunctionBox(r0, r1, mag)
                f_je = [f_jex, f_jey, f_jez]
            else:
                raise ValueError()
            
            return f_je
        elif em_type=='m' or em_type=='h':
            i_xyz_min, r_xyz_min, d_xyz_min = self.FindClosestFaces(r_0)
            f_jm = None
            if src_dir=='x':
                i_j, r_j, d_j = i_xyz_min[0], r_xyz_min[0], d_xyz_min[0]
                d_min = min(self.dr)

                r0 = np.copy(self.r0)
                r1 = np.copy(self.r1)
                if norm_dir=='x':
                    r0[0] = r_j[0]-d_min/3
                    r1[0] = r_j[0]+d_min/3
                elif norm_dir=='y':
                    assert self.n_dim>=2
                    r0[1] = r_j[1]-d_min/3
                    r1[1] = r_j[1]+d_min/3
                elif norm_dir=='z':
                    assert self.n_dim==3
                    r0[2] = r_j[2]-d_min/3
                    r1[2] = r_j[2]+d_min/3

                f_jmx = self.FunctionBox(r0, r1, mag)
                f_jmy = lambda r : 0.0
                f_jmz = lambda r : 0.0
                f_jm = [f_jmx, f_jmy, f_jmz]
            elif src_dir=='y':
                i_j, r_j, d_j = i_xyz_min[1], r_xyz_min[1], d_xyz_min[1]
                d_min = min(self.dr)

                r0 = np.copy(self.r0)
                r1 = np.copy(self.r1)
                if norm_dir=='x':
                    r0[0] = r_j[0]-d_min/3
                    r1[0] = r_j[0]+d_min/3
                elif norm_dir=='y':
                    assert self.n_dim>=2
                    r0[1] = r_j[1]-d_min/3
                    r1[1] = r_j[1]+d_min/3
                elif norm_dir=='z':
                    assert self.n_dim==3
                    r0[2] = r_j[2]-d_min/3
                    r1[2] = r_j[2]+d_min/3
                
                f_jmx = lambda r : 0.0
                f_jmy = self.FunctionBox(r0, r1, mag)
                f_jmz = lambda r : 0.0
                f_jm = [f_jmx, f_jmy, f_jmz]
            elif src_dir=='z':
                i_j, r_j, d_j = i_xyz_min[2], r_xyz_min[2], d_xyz_min[2]
                d_min = min(self.dr)

                r0 = np.copy(self.r0)
                r1 = np.copy(self.r1)
                if norm_dir=='x':
                    r0[0] = r_j[0]-d_min/3
                    r1[0] = r_j[0]+d_min/3
                elif norm_dir=='y':
                    assert self.n_dim>=2
                    r0[1] = r_j[1]-d_min/3
                    r1[1] = r_j[1]+d_min/3
                elif norm_dir=='z':
                    assert self.n_dim==3
                    r0[2] = r_j[2]-d_min/3
                    r1[2] = r_j[2]+d_min/3

                f_jmx = lambda r : 0.0
                f_jmy = lambda r : 0.0
                f_jmz = self.FunctionBox(r0, r1, mag)
                f_jm = [f_jmx, f_jmy, f_jmz]
            else:
                raise ValueError()

            return f_jm
        else:
            raise ValueError()


    def GetBoxSourceFunc(self, r0, r1, mag=1.0, src_dir='x', em_type='e'):
        assert src_dir in ['x', 'y', 'z']
        assert em_type in ['e', 'm', 'h']
        if em_type=='e':
            f_je = None
            if src_dir=='x':
                f_jex = self.FunctionBox(r0, r1, mag)
                f_jey = lambda r : 0.0
                f_jez = lambda r : 0.0
                f_je = [f_jex, f_jey, f_jez]
            elif src_dir=='y':
                f_jex = lambda r : 0.0
                f_jey = self.FunctionBox(r0, r1, mag)
                f_jez = lambda r : 0.0
                f_je = [f_jex, f_jey, f_jez]
            elif src_dir=='z':
                f_jex = lambda r : 0.0
                f_jey = lambda r : 0.0
                f_jez = self.FunctionBox(r0, r1, mag)
                f_je = [f_jex, f_jey, f_jez]
            else:
                raise ValueError()
            
            return f_je
        elif em_type=='m' or em_type=='h':
            f_jm = None
            if src_dir=='x':
                f_jmx = self.FunctionBox(r0, r1, mag)
                f_jmy = lambda r : 0.0
                f_jmz = lambda r : 0.0
                f_jm = [f_jmx, f_jmy, f_jmz]
            elif src_dir=='y':
                f_jmx = lambda r : 0.0
                f_jmy = self.FunctionBox(r0, r1, mag)
                f_jmz = lambda r : 0.0
                f_jm = [f_jmx, f_jmy, f_jmz]
            elif src_dir=='z':
                f_jmx = lambda r : 0.0
                f_jmy = lambda r : 0.0
                f_jmz = self.FunctionBox(r0, r1, mag)
                f_jm = [f_jmx, f_jmy, f_jmz]
            else:
                raise ValueError()

            return f_jm
        else:
            raise ValueError()

            

    def FindClosestSides(self, r0):
        """ r0:point
        """
        #self.vbose=True
        drx2=None
        if self.n_dim==3:
            drx2 = (self.rsx[0]-r0[0])**2 + (self.rsx[1]-r0[1])**2 + (self.rsx[2]-r0[2])**2
        elif self.n_dim==2:
            drx2 = (self.rsx[0]-r0[0])**2 + (self.rsx[1]-r0[1])**2
        elif self.n_dim==1:
            drx2 = (self.rsx[0]-r0[0])**2
        inds_minx_F = np.argmin(drx2)
        inds_min_x = np.unravel_index(inds_minx_F, self.Nsx)
        x_min = np.array([self.rsx[i][inds_min_x] for i in range(self.n_dim)])
        dx_min = np.sqrt(np.min(drx2))
        if self.vbose:
            print('sides inds_min_x:', inds_min_x, 'd:', dx_min)

        dry2 = None
        if self.n_dim==3:
            dry2 = (self.rsy[0]-r0[0])**2 + (self.rsy[1]-r0[1])**2 + (self.rsy[2]-r0[2])**2
        elif self.n_dim==2:
            dry2 = (self.rsy[0]-r0[0])**2 + (self.rsy[1]-r0[1])**2
        elif self.n_dim==1:
            dry2 = (self.rsy[0]-r0[0])**2
        inds_miny_F = np.argmin(dry2)
        inds_min_y = np.unravel_index(inds_miny_F, self.Nsy)
        y_min = np.array([self.rsy[i][inds_min_y] for  i in range(self.n_dim)])
        dy_min = np.sqrt(np.min(dry2))
        if self.vbose:
            print('sides inds_min_y:', inds_min_y, 'd:', np.sqrt(np.min(dry2)))

        drz2 = None
        if self.n_dim==3:
            drz2 = (self.rsz[0]-r0[0])**2 + (self.rsz[1]-r0[1])**2 + (self.rsz[2]-r0[2])**2
        elif self.n_dim==2:
            drz2 = (self.rsz[0]-r0[0])**2 + (self.rsz[1]-r0[1])**2 
        elif self.n_dim==1:
            drz2 = (self.rsz[0]-r0[0])**2 
        inds_minz_F = np.argmin(drz2)
        inds_min_z = np.unravel_index(inds_minz_F, self.Nsz)
        z_min = np.array([self.rsz[i][inds_min_z] for i in range(self.n_dim)])
        dz_min = np.sqrt(np.min(drz2))
        if self.vbose:
            print('sides inds_min_z:', inds_min_z, 'd:', np.sqrt(np.min(drz2)))
        
        return [inds_min_x, inds_min_y, inds_min_z], [x_min, y_min, z_min],\
            np.array([dx_min, dy_min, dz_min])


    def FindClosestFaces(self, r0):
        """ r0:point
        """
        #self.vbose=True
        drx2=None
        if self.n_dim==3:
            drx2 = (self.rfx[0]-r0[0])**2 + (self.rfx[1]-r0[1])**2 + (self.rfx[2]-r0[2])**2
        elif self.n_dim==2:
            drx2 = (self.rfx[0]-r0[0])**2 + (self.rfx[1]-r0[1])**2
        elif self.n_dim==1:
            drx2 = (self.rfx[0]-r0[0])**2
        inds_minx_F = np.argmin(drx2)
        inds_min_x = np.unravel_index(inds_minx_F, self.Nfx)
        x_min = np.array([self.rfx[i][inds_min_x] for i in range(self.n_dim)])
        dx_min = np.sqrt(np.min(drx2))
        if self.vbose:
            print('face inds_min_x:', inds_min_x, 'd:', dx_min)

        dry2 = None
        if self.n_dim==3:
            dry2 = (self.rfy[0]-r0[0])**2 + (self.rfy[1]-r0[1])**2 + (self.rfy[2]-r0[2])**2
        elif self.n_dim==2:
            dry2 = (self.rfy[0]-r0[0])**2 + (self.rfy[1]-r0[1])**2
        elif self.n_dim==1:
            dry2 = (self.rfy[0]-r0[0])**2
        inds_miny_F = np.argmin(dry2)
        inds_min_y = np.unravel_index(inds_miny_F, self.Nfy)
        y_min = np.array([self.rfy[i][inds_min_y] for  i in range(self.n_dim)])
        dy_min = np.sqrt(np.min(dry2))
        if self.vbose:
            print('face inds_min_y:', inds_min_y, 'd:', np.sqrt(np.min(dry2)))

        drz2 = None
        if self.n_dim==3:
            drz2 = (self.rfz[0]-r0[0])**2 + (self.rfz[1]-r0[1])**2 + (self.rfz[2]-r0[2])**2
        elif self.n_dim==2:
            drz2 = (self.rfz[0]-r0[0])**2 + (self.rfz[1]-r0[1])**2 
        elif self.n_dim==1:
            drz2 = (self.rfz[0]-r0[0])**2 
        inds_minz_F = np.argmin(drz2)
        inds_min_z = np.unravel_index(inds_minz_F, self.Nfz)
        z_min = np.array([self.rfz[i][inds_min_z] for i in range(self.n_dim)])
        dz_min = np.sqrt(np.min(drz2))
        if self.vbose:
            print('face inds_min_z:', inds_min_z, 'd:', np.sqrt(np.min(drz2)))
        
        return [inds_min_x, inds_min_y, inds_min_z], [x_min, y_min, z_min],\
            np.array([dx_min, dy_min, dz_min])


    def SetViewPlane_Side(self, r0, A_dir_list):
        """ 3 normal planes passing through r0 will be returned 
            A: the fields to output for example [{'A':self.Ex, 'A_dir:''x', 
            'O_dir':'y', 'v_out':v), ..]
            var_out: save output to this variable (a list) when demanded
        """
        inds = self.FindClosestSides(r0)[0]
        ind_x, ind_y, ind_z = inds
        
        Ellipsis = slice(None, None, None)
        
        Ax_r = [None]*self.n_dim
        ix_p = [None]*self.n_dim
        for n in range(self.n_dim):
            ix_p[n] = [Ellipsis]*self.n_dim    ##plane normal to n
            ix_p[n][n] = ind_x[n]
            Ax_r[n] = [self.rsx[i][ix_p[n]] for i in range(self.n_dim)]
        
        ix_all = [Ellipsis]*self.n_dim   
        Ax_r_all = self.rsx
        if self.n_dim==0:
            Ax_r_all = self.rsx[0]


        Ay_r = [None]*self.n_dim
        iy_p = [None]*self.n_dim
        for n in range(self.n_dim):
            iy_p[n] = [Ellipsis]*self.n_dim    ##plane normal to n
            iy_p[n][n] = ind_y[n]
            Ay_r[n] = [self.rsy[i][iy_p[n]] for i in range(self.n_dim)]
        
        iy_all = [Ellipsis]*self.n_dim   
        Ay_r_all = self.rsy
        if self.n_dim==0:
            Ay_r_all = self.rsy[0]


        Az_r = [None]*self.n_dim
        iz_p = [None]*self.n_dim
        for n in range(self.n_dim):
            iz_p[n] = [Ellipsis]*self.n_dim    ##plane normal to n
            iz_p[n][n] = ind_z[n]
            Az_r[n] = [self.rsz[i][iz_p[n]] for i in range(self.n_dim)]
        
        iz_all = [Ellipsis]*self.n_dim   
        Az_r_all = self.rsz
        if self.n_dim==0:
            Az_r_all = self.rsz[0]


        for view_dic in A_dir_list:
            #print('view_dic:', view_dic)
            A = view_dic['A']           ##field
            A_dir = view_dic['A_dir']   ##field direction
            O_dir = view_dic['O_dir']   ##output plane is normal to this direction
            name = view_dic['name']   ##output list
            v_out = []
            save_folder = None
            save_size = None
            
            if "save_folder" in view_dic:
                save_folder = view_dic["save_folder"]
                if save_folder is not None:
                    assert os.path.exists(save_folder)
                    save_size = view_dic["save_size"]
                    assert isinstance(save_size, int) and save_size>0
                
            
            assert name not in self.ViewPlanes
            
            if A_dir=='x':
                n_dir=0
                if O_dir=='x':
                    self.ViewPlanes[name]=  {'v_out':v_out, 'A':A, 'A_dir':n_dir, \
                        'ind':ix_p[0], 'r':Ax_r[0], "save_folder":save_folder, "save_size":save_size}
                elif O_dir=='y':
                    assert self.n_dim>=2
                    self.ViewPlanes[name] = {'v_out':v_out, 'A':A, 'A_dir':n_dir, \
                    'ind':ix_p[1], 'r':Ax_r[1], "save_folder":save_folder, "save_size":save_size}
                elif O_dir=='z':
                    assert self.n_dim==3
                    self.ViewPlanes[name] = {'v_out':v_out, 'A':A, 'A_dir':n_dir, \
                    'ind':ix_p[2], 'r':Ax_r[2], "save_folder":save_folder, "save_size":save_size}
                elif O_dir is None:
                    self.ViewPlanes[name] = {'v_out':v_out, 'A':A, 'A_dir':n_dir, \
                    'ind':ix_all, 'r':Ax_r_all, "save_folder":save_folder, "save_size":save_size}
                else:
                    raise ValueError()
            elif A_dir=='y':
                n_dir=1
                if O_dir=='x':
                    self.ViewPlanes[name] = {'v_out':v_out, 'A':A, 'A_dir':n_dir, \
                    'ind':iy_p[0], 'r':Ay_r[0], "save_folder":save_folder, "save_size":save_size}
                elif O_dir=='y':
                    assert self.n_dim>=2
                    self.ViewPlanes[name] = {'v_out':v_out, 'A':A, 'A_dir':n_dir, \
                    'ind':iy_p[1], 'r':Ay_r[1], "save_folder":save_folder, "save_size":save_size}
                elif O_dir=='z':
                    assert self.n_dim==3
                    self.ViewPlanes[name] = {'v_out':v_out, 'A':A, 'A_dir':n_dir, \
                    'ind':iy_p[2], 'r':Ay_r[2], "save_folder":save_folder, "save_size":save_size}
                elif O_dir is None:
                    self.ViewPlanes[name] = {'v_out':v_out, 'A':A, 'A_dir':n_dir, \
                    'ind':iy_all, 'r':Ay_r_all, "save_folder":save_folder, "save_size":save_size}
                else:
                    raise ValueError()
            elif A_dir=='z':
                n_dir=2
                if O_dir=='x':
                    self.ViewPlanes[name] = {'v_out':v_out, 'A':A, 'A_dir':n_dir, \
                    'ind':iz_p[0], 'r':Az_r[0], "save_folder":save_folder, "save_size":save_size}
                elif O_dir=='y':
                    assert self.n_dim>=2
                    self.ViewPlanes[name] = {'v_out':v_out, 'A':A, 'A_dir':n_dir, \
                    'ind':iz_p[1], 'r':Az_r[1], "save_folder":save_folder, "save_size":save_size}
                elif O_dir=='z':
                    assert self.n_dim==3
                    self.ViewPlanes[name] = {'v_out':v_out, 'A':A, 'A_dir':n_dir, \
                    'ind':iz_p[2], 'r':Az_r[2], "save_folder":save_folder, "save_size":save_size}
                elif O_dir is None:
                    self.ViewPlanes[name] = {'v_out':v_out, 'A':A, 'A_dir':n_dir, \
                    'ind':iz_all, 'r':Az_r_all, "save_folder":save_folder, "save_size":save_size}
                else:
                    raise ValueError()
            else:
                raise NotImplementedError()
        if self.vbose:
            print('self.ViewPlanes: ', self.ViewPlanes.keys())
                    
    
    def SaveOutputs(self):
        for name in self.ViewPlanes:
            vp = self.ViewPlanes[name]
            v_out = vp['v_out']
            A = vp['A']
            n_dir = vp['A_dir']
            ind = vp['ind']
            r = vp['r']
            
            ## TODO: ind should be modified...
            ## FutureWarning: Using a non-tuple sequence for multidimensional indexing is deprecated; 
            ## use `arr[tuple(seq)]` instead of `arr[seq]`
            v_out.append(np.copy(A[n_dir][ind]))    
            
            save_folder = vp["save_folder"]
            #print("save_folder : ", vp)
            if save_folder is not None:
                save_size = vp["save_size"]
                #print(len(v_out))
                if len(v_out)>=save_size:
                    ind_last = 0
                    if "save_ind_last" in vp:
                        ind_last = vp["save_ind_last"]
                    else:
                        vp["save_ind_last"] = ind_last
                    file_name = os.path.join(save_folder, name+"_"+str(ind_last))
                    np.savez(file_name, arr=np.array(v_out, dtype=object))
                    
                    vp['v_out'] = []
                    vp["save_ind_last"] += 1
                    
                    if self.vbose:
                        print(name+"_"+str(ind_last)+ " saved! ")

    def DumpRemainingOutputsToFile(self):
        for name in self.ViewPlanes:
            vp = self.ViewPlanes[name]
            v_out = vp['v_out']
            A = vp['A']
            n_dir = vp['A_dir']
            ind = vp['ind']
            r = vp['r']
            
            save_folder = vp["save_folder"]
            #print("save_folder : ", vp)
            if save_folder is not None:
                save_size = vp["save_size"]
                #print(len(v_out))
                if len(v_out)>0:
                    ind_last = 0
                    if "save_ind_last" in vp:
                        ind_last = vp["save_ind_last"]
                    else:
                        vp["save_ind_last"] = ind_last
                    file_name = os.path.join(save_folder, name+"_"+str(ind_last))
                    np.savez(file_name, arr=np.array(v_out))
                    
                    vp['v_out'] = []
                    vp["save_ind_last"] += 1
                    
                    if self.vbose:
                        print(name+"_"+str(ind_last)+ " saved! ")
        
        
    def GetOutputs(self, name, file_ind=0):
        if self.vbose:
            print('Saved outputs:', self.ViewPlanes.keys())
        if name in self.ViewPlanes:
            vp = self.ViewPlanes[name]
            v_out = vp['v_out']
            A = vp['A']
            n_dir = vp['A_dir']
            ind = vp['ind']
            r = vp['r']
            
            save_folder = vp["save_folder"]
            if save_folder is None:
                return [r, v_out]
            else:
                ind_max = vp["save_ind_last"]
                assert file_ind<ind_max
                file_name = os.path.join(save_folder, name+"_"+str(file_ind)+'.npz')
                
                file_data = np.load(file_name)
                v_out = file_data['arr'].tolist()
                #print(v_out[0])
                
                return [r, v_out, ind_max]

        else:
            if self.vbose:
                print('Variable %s does not exist!'%(name))
            return None
        
        
    def GetAllOutputs(self):
        out = {}
        if name in self.ViewPlanes:
            vp = self.ViewPlanes[name]
            v_out = vp['v_out']
            A = vp['A']
            n_dir = vp['A_dir']
            ind = vp['ind']
            r = vp['r']
            
            out[name] =  [r, v_out]
        return out


