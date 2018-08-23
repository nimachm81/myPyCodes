## FEM.py finite elements


__all__ = ["BasisTypes", "AssemblyProcess", "CoeffTypes", "BCTypes",
            "FEMElemNode1D", "MeshGen1D", "FEM1D",
            "FEMElemNode2D", "MeshGen2D", "FEM2D"]

import numpy as np
import scipy
from scipy.misc import factorial
from numpy.polynomial import polynomial

from enum import Enum
import copy



class BasisTypes(Enum):
    lagrange_1D = 1
    nodeSimplex_2D = 2

class AssemblyProcess(Enum):
    galerkin = 1
    
class CoeffTypes(Enum):
    constant = 1
    polynomial_1D = 2

class BCTypes(Enum):
    Dirichlet = 1
    Neumann = 2
    Robin = 3


class FEMElemNode1D:
    
    def __init__(self, elemOrder=1, basisType=BasisTypes.lagrange_1D, x0=0.0, x1=1.0):
        self.elemOrder = elemOrder
        self.basisType = basisType
        self.BB_int = {}
        self.BP_int = {}
        
        self.x0 = x0
        self.x1 = x1
        
        self.setupBasisFuncs()
        return


    def setupBasisFuncs(self):
        if self.basisType==BasisTypes.lagrange_1D:
            from scipy.interpolate import lagrange
            n_base = self.elemOrder+1
            bases_poly = [None]*n_base
            for i in range(n_base):
                x_vec = np.linspace(self.x0, self.x1, n_base)
                y_vec = np.zeros(n_base)
                y_vec[i] = 1.0
                b_i_poly = lagrange(x_vec, y_vec)
                bases_poly[i] = b_i_poly
            self.B_lagr = bases_poly
    

    def getOvIntdmBidnBj(self, m=0, n=0):
        if self.basisType==BasisTypes.lagrange_1D:
            n_base = self.elemOrder+1
            B_int = np.zeros((n_base, n_base))
            for i in range(n_base):
                for j in range(n_base):
                    p_ij = self.B_lagr[i].deriv(m)*self.B_lagr[j].deriv(n)
                    p_ij_int = p_ij.integ()
                    B_int[i, j] = p_ij_int(self.x1) - p_ij_int(self.x0)
            self.BB_int[(m,n)] = B_int
            return B_int
        

    def getOvIntdmBiP(self, m=0, poly_order=2):
        if self.basisType==BasisTypes.lagrange_1D:
            n_base = self.elemOrder+1
            B_int = np.zeros((n_base, poly_order+1))
            for i in range(n_base):
                for j in range(poly_order+1):
                    p_j = np.zeros(poly_order+1)
                    p_j[poly_order-j] = 1.0
                    p_ij = self.B_lagr[i].deriv(m)*np.poly1d(p_j)
                    p_ij_int = p_ij.integ()
                    B_int[i, j] = p_ij_int(self.x1) - p_ij_int(self.x0)
            self.BP_int[m] = B_int
            return B_int

        
    def getBasisPlot(self, N=100):
        x = np.linspace(self.x0, self.x1, N)
        n_base = (self.elemOrder+1)
        bases = [None]*n_base
        if self.basisType==BasisTypes.lagrange_1D:
            for i in range(n_base):
                bases[i] = self.B_lagr[i](x)
        else:
            raise NotImplementedError()
        return x, bases
            
            
    def getFieldPlot(self, coeffs, x0, x1, N=100):
        n_base = (self.elemOrder+1)
        assert len(coeffs)==n_base
        x_ = np.linspace(self.x0, self.x1, N)
        x = np.linspace(x0, x1, N)
        fx = np.zeros(N)
        if self.basisType==BasisTypes.lagrange_1D:
            for i in range(n_base):
                fx += coeffs[i]*self.B_lagr[i](x_)
        else:
            raise NotImplementedError()
        return x, fx
        
        

        
class MeshGen1D:

    def __init__(self, x0, x1, dx_max):
        self.x0 = x0
        self.x1 = x1
        self.dx_max = dx_max
        return

    def GenerateUniformMesh(self):
        N = int((self.x1 - self.x0)/self.dx_max + 1)
        nodes = np.linspace(self.x0, self.x1, N)
        sides = np.zeros((N-1, 2), dtype=int)
        sides[:, 0] = np.arange(0, N-1)
        sides[:, 1] = np.arange(1, N)
        nodes_labels = -np.ones(N, dtype=int)
        nodes_labels[0] = 1
        nodes_labels[-1] = 2
        sides_labels = -np.ones(N-1, dtype=int)
        
        self.nodes = nodes
        self.sides = sides
        self.nodes_labels = nodes_labels
        self.sides_labels = sides_labels
        return nodes, sides, nodes_labels, sides_labels




class FEM1D:
    
    def __init__(self):
        return
        
    
    def SetStruct(self, nodes, sides, nodelabels=None, sidelabels=None):
        self.nodes = nodes
        self.sides = sides
        self.nodelabels = nodelabels
        self.sidelabels = sidelabels
        self.dx = nodes[sides[:,1]] - nodes[sides[:,0]]
        
        self.SetLabelsAssociatedSidesOrNodes()
        self.SetNodesConnectedSides()
        
        self.assemblyParams = []
        return
        
                
    def SetLabelsAssociatedSidesOrNodes(self):
        labels_sides = {}
        n_sides = len(self.sides)
        if self.sidelabels is None:
            labels_sides[None] = np.arange(n_sides)
        else:
            labelsAll = set(self.sidelabels)
            for label in labelsAll:
                labels_sides[label] = np.where(self.sidelabels==label)[0]
        self.labels_sides = labels_sides
            
        labels_nodes = {}
        n_nodes = len(self.nodes)
        if self.nodelabels is None:
            labels_nodes[None] = np.arange(n_nodes)
        else:
            labelsAll = set(self.nodelabels)
            for label in labelsAll:
                labels_nodes[label] = np.where(self.nodelabels==label)[0]
        self.labels_nodes = labels_nodes
        return


    def SetNodesConnectedSides(self):
        n_nodes = len(self.nodes)
        n_sides = len(self.sides)
        nodes_connectedSides = -np.ones((n_nodes, 2), dtype=int)
        nodes_connectedSides[self.sides[:, 0], 1] = np.arange(n_sides)
        nodes_connectedSides[self.sides[:, 1], 0] = np.arange(n_sides)
        self.nodes_connectedSides = nodes_connectedSides

        
    def SetGlobalIndices(self, elem_order, elem_type):
        self.elem_order = elem_order
        self.elem_type = elem_type
        if elem_type==BasisTypes.lagrange_1D:
            n_nodes = len(self.nodes)
            n_sides = len(self.sides)
            n_base_T = n_sides*(elem_order-1) + n_nodes
            
            sidesBases = np.zeros((n_sides, elem_order+1), dtype=int)
            
            sidesBases[:, 0] = self.sides[:, 0]
            sidesBases[:, elem_order] = self.sides[:, 1]
            
            n_b_inside = elem_order-1
    
            if n_b_inside>0:
                ind_last = n_nodes
                for s in range(n_sides):
                    sidesBases[s, 1:elem_order] = np.arange(ind_last, ind_last+n_b_inside)
                    ind_last += n_b_inside
            
            ##TODO: reorder based on element locality
            
            self.n_base_T = n_base_T
            self.sidesBases = sidesBases
            return sidesBases
        else:
            raise NotImplementedError()


    def DefineMatrixAssemblyProcess(self, description, label, proc, procParams):
        """ example:description: 'mat'/'rhs'/'bc'
                    proc=AssemblyProcess.galerkin
                    procParams = {'coeffType':CoeffTypes.constant, 'C':2.0, 'dUi':1, 'dUj':1}
            for integral over cell 2*dw/dx*du/dx,   
            i.e. C=2.0, derivative order w=1, derivative order u=1
            w=wighting function ---> U_i
            u=basis function  -----> U_j
            
            'mat' ---> stifness matrix      label --> sidelabel
            'rhs' ---> rhs                  label --> sidelabel
            'bc'  ---> boundary condition   label --> nodelabel
            
            for 'bc'
            procParams = {'bcType':BCTypes.Dirichlet, 'bcParams':bcParams}
                bcParams:
                    Dirichlet : 'rhs'   ---->  u(0) = rhs
                    Neumann   : 'rhs'   ---->  u'(0) = rhs
                    Robin     : 'sigma', 'rhs' ---> u'(0)-sigma*u(0)=rhs
        """
        assert description in ['mat', 'rhs', 'bc']
        
        self.assemblyParams.append({'description':description, 'label':label, \
                'proc':proc, 'procParams':procParams}) 


    def AssembleMatrix(self):
        row, col, data = np.array([], dtype=int), np.array([], dtype=int), np.array([])
        row_rhs, data_rhs = np.array([], dtype=int), np.array([])
        row_preset, data_preset = np.array([], dtype=int), np.array([])
        for assproc in self.assemblyParams:
            description = assproc['description']
            proc = assproc['proc']
            if description=='mat' and proc==AssemblyProcess.galerkin:
                label = assproc['label']
                procParams = assproc['procParams']
                
                assert label in self.labels_sides
                lab_sides = self.labels_sides[label]
                row_, col_, data_ = self.GetGalerkinMatElems(lab_sides, procParams)
                
                row = np.append(row, row_)
                col = np.append(col, col_)
                data = np.append(data, data_)
                
            elif description=='rhs' and proc==AssemblyProcess.galerkin:
                label = assproc['label']
                procParams = assproc['procParams']

                assert label in self.labels_sides
                lab_sides = self.labels_sides[label]

                row_rhs_, data_rhs_ = self.GetGalerkinRhs(lab_sides, procParams)

                row_rhs = np.append(row_rhs, row_rhs_)
                data_rhs = np.append(data_rhs, data_rhs_)
                
            elif description=='bc' and proc==AssemblyProcess.galerkin:
                label = assproc['label']
                procParams = assproc['procParams']

                assert label in self.labels_nodes
                lab_nodes = self.labels_nodes[label]

                bcType = procParams['bcType']
                bcParams = procParams['bcParams']
                
                if bcType==BCTypes.Dirichlet:
                    row_ps, data_ps = self.GetGalerkinDirichlet(lab_nodes, bcParams)

                    row_preset = np.append(row_preset, row_ps)
                    data_preset = np.append(data_preset, data_ps)
                    
                elif bcType==BCTypes.Neumann:
                    raise NotImplementedError()
                elif bcType==BCTypes.Robin:
                    raise NotImplementedError()
                else:
                    raise NotImplementedError()
                
            else:
                raise NotImplementedError()
            
        #print(row, col, data, sep='\n')
        
        ##reset presets
        mask = np.ones(len(row), dtype=bool)
        if len(row_preset)>0:
            mask = np.logical_not(row==row_preset[0])
        for i in range(1, len(row_preset)):
            mask = np.logical_and(mask, np.logical_not(row==row_preset[i]))
        
        row = row[mask]
        col = col[mask]
        data = data[mask]
        
        row = np.append(row, row_preset)
        col = np.append(col, row_preset)
        data = np.append(data, np.ones(len(row_preset)))
        
        mask = np.ones(len(row), dtype=bool)
        if len(row_preset)>0:
            mask_rhs = np.logical_not(row_rhs==row_preset[0])
        for i in range(1, len(row_preset)):
            mask_rhs = np.logical_and(mask_rhs, np.logical_not(row_rhs==row_preset[i]))
        
        row_rhs = row_rhs[mask_rhs]
        data_rhs = data_rhs[mask_rhs]
        
        row_rhs = np.append(row_rhs, row_preset)
        data_rhs = np.append(data_rhs, data_preset)
        
        ## generate sparse matrix, solve
        A_coo = scipy.sparse.coo_matrix((data, (row, col)), shape=(self.n_base_T, self.n_base_T))
        #print('A_coo : \n', A_coo)
        
        b_vec = np.zeros(self.n_base_T)
        for i in range(len(row_rhs)):
            b_vec[row_rhs[i]] += data_rhs[i]
            
        from scipy.sparse import linalg
        A_csc = A_coo.tocsc()
        x_res = linalg.spsolve(A_csc, b_vec)
        
        self.x_res = x_res
        return x_res
        
        
    def GetGalerkinMatElems(self, sides_inds, params):
        sides = self.sides
        nodes = self.nodes
        dx = self.dx
        sidesBases = self.sidesBases
        if params['coeffType']==CoeffTypes.constant:
            coeff = params['C']
            dUi = params['dUi']
            dUj = params['dUj']
            
            elem_order = self.elem_order
            fem_elem = FEMElemNode1D(elem_order, basisType=self.elem_type, x0=0.0, x1=1.0)
            ovInts = fem_elem.getOvIntdmBidnBj(m=dUi, n=dUj)
            
            row, col, data = [], [], []
            for s in sides_inds:
                ## the integral scales as dx and the derivatives scale as 1/dx**m
                mat_s = ovInts/dx[s]**(dUi+dUj-1)
                shape = mat_s.shape
                row_s = np.zeros(shape, dtype=int)
                col_s = np.zeros(shape, dtype=int)
                for i in range(shape[0]):
                    row_s[i,:] = np.ones(shape[1], dtype=int)*sidesBases[s, i]
                    col_s[:,i] = np.ones(shape[1], dtype=int)*sidesBases[s, i]
                
                #print(row_s, col_s, sep='\n')
                
                data_s = mat_s.reshape(shape[0]*shape[1])
                row_s = row_s.reshape(shape[0]*shape[1])
                col_s = col_s.reshape(shape[0]*shape[1])
                
                nnz_s = data_s!=0.0
                
                row_s = row_s[nnz_s]
                col_s = col_s[nnz_s]
                data_s = data_s[nnz_s]
                
                row.extend(row_s)
                col.extend(col_s)
                data.extend(data_s)
                
            return np.array(row), np.array(col), coeff*np.array(data)

        else:
            raise NotImplementedError()
                    
        
    
    def GetGalerkinRhs(self, sides_inds, params):
        sides = self.sides
        nodes = self.nodes
        dx = self.dx
        sidesBases = self.sidesBases
        if params['coeffType']==CoeffTypes.constant:
            coeff = params['C']
            dUi = params['dUi']
            
            elem_order = self.elem_order
            fem_elem = FEMElemNode1D(elem_order, basisType=self.elem_type, x0=0.0, x1=1.0)
            BPInts = fem_elem.getOvIntdmBiP(m=dUi, poly_order=0)
             
            row, data = [], []
            for s in sides_inds:
                rhs_s = BPInts[:,0]/dx[s]**(dUi-1)
                shape = rhs_s.shape
                row_s = np.zeros(shape[0], dtype=int)
                for i in range(shape[0]):
                    row_s[i] = sidesBases[s, i]
                                
                nnz_s = rhs_s!=0.0
                
                row_s = row_s[nnz_s]
                data_s = rhs_s[nnz_s]
                
                row.extend(row_s)
                data.extend(data_s)
                
            return np.array(row), coeff*np.array(data)

        else:
            raise NotImplementedError()
    
    

    def GetGalerkinDirichlet(self, lab_nodes, bcParams):
        sides = self.sides
        nodes = self.nodes
        dx = self.dx
        sidesBases = self.sidesBases
        nodes_connectedSides = self.nodes_connectedSides
        
        bc_rhs = bcParams['rhs']
        row_ps, data_ps = [], []
        if self.elem_type==BasisTypes.lagrange_1D:
            for n in lab_nodes:
                if nodes_connectedSides[n, 0]<0:
                    s = nodes_connectedSides[n, 1]
                    assert s>=0
                    row_ps.append(sidesBases[s, 0])
                    data_ps.append(bc_rhs)
            
                elif nodes_connectedSides[n, 1]<0:
                    s = nodes_connectedSides[n, 0]
                    assert s>=0
                    row_ps.append(sidesBases[s, -1])
                    data_ps.append(bc_rhs)
        else:
            raise NotImplementedError()
        
        return np.array(row_ps), np.array(data_ps)
    
        
    def GetArrangedSides(self):
        sides = self.sides
        nodes = self.nodes
        dx = self.dx
        sidesBases = self.sidesBases
        nodes_connectedSides = self.nodes_connectedSides
        n_sides = len(self.sides)
        
        n0_ind = np.argmax(nodes_connectedSides[:,0]<0)
        assert nodes_connectedSides[n0_ind,0]<0
        
        s0 = nodes_connectedSides[n0_ind,1]
        assert s0>=0
        sides_arranged = -np.ones(n_sides, dtype=int)
        _s_ind_ = 0
        n_last = n0_ind
        while True:
            s_last = nodes_connectedSides[n_last, 1]
            if s_last>=0:
                sides_arranged[_s_ind_] = s_last
                _s_ind_ += 1
            else:
                break
            n_last = sides[s_last, 1]
        assert np.all(sides_arranged>=0)
        self.sides_arranged = sides_arranged
        return sides_arranged
        
    
    def GetFieldPlot(self):
        self.GetArrangedSides()
        shape = self.sidesBases.shape
        x_res_ = np.zeros(shape)
        n_sides = len(self.sides)
        for i in range(n_sides):
            for j in range(self.elem_order+1):
                x_res_[i, j] = self.x_res[self.sidesBases[self.sides_arranged[i], j]]

        femelem = FEMElemNode1D(elemOrder=self.elem_order, basisType=self.elem_type)
        x, fx = [], []
        for i in range(n_sides):
            coeffs_i = x_res_[i, :]
            x0 = self.nodes[self.sides[self.sides_arranged[i]][0]]
            x1 = self.nodes[self.sides[self.sides_arranged[i]][1]]
            x_i, fx_i = femelem.getFieldPlot(coeffs_i, x0, x1, N=10)
            x.extend(x_i)
            fx.extend(fx_i)
        
        return x, fx
    
    


##-------------------------------------  2D --------------------------

class FEMElemNode2D:
    
    def __init__(self, elemOrder=1, basisType=BasisTypes.nodeSimplex_2D, \
            r0=np.array([0.0, 0.0]), r1=np.array([1.0, 0.0]), r2=np.array([0.0, 1.0])):
        self.elemOrder = elemOrder
        self.basisType = basisType
        self.BB_int_divA = {}
        self.BP_int = {}
        
        self.r0 = r0
        self.r1 = r1
        self.r2 = r2
        
        self.setupBasisFuncs()
        return


    def setupBasisFuncs(self):
        if self.basisType==BasisTypes.nodeSimplex_2D:
            n_base = self.elemOrder+1
            bases_poly = [None]*n_base
            bases_poly[0] = polynomial.Polynomial([1.0], domain=[0, 1], window=[0, 1])
            for i in range(1, n_base):
                bases_poly[i] = polynomial.Polynomial([1.0], domain=[0, 1], window=[0, 1])
                for j in range(i):
                    bases_poly[i] *= polynomial.Polynomial([-j, self.elemOrder], domain=[0, 1], window=[0, 1])
                bases_poly[i] /= float(factorial(i))
            self.B_NodeSimplex = bases_poly

    def getBasisPlot(self, b_ind, n_pts=100):
        """ b_ind: basis index example (I, J, K) where I+J+K = n_base
        """
        if self.basisType==BasisTypes.nodeSimplex_2D:
            assert len(b_ind)==3
            I, J, K = b_ind
            N = self.elemOrder
            if I+J+K!=N:
                raise ValueError("I+J+K=N not satisfied.")
            L0_arr = np.linspace(0.0, 1.0, n_pts)
            L1_arr = np.linspace(0.0, 1.0, n_pts)
            L0_mesh, L1_mesh = np.meshgrid(L0_arr, L1_arr)
            L2_mesh = 1.0 - L0_mesh - L1_mesh
            x_mesh = L0_mesh*self.r0[0] + L1_mesh*self.r1[0] + L2_mesh*self.r2[0]
            y_mesh = L0_mesh*self.r0[1] + L1_mesh*self.r1[1] + L2_mesh*self.r2[1]
            
                
            b_mesh = self.B_NodeSimplex[I](L0_mesh)*self.B_NodeSimplex[J](L1_mesh)\
                    *self.B_NodeSimplex[K](L2_mesh)    
                    
            ##mask
            mask = L0_mesh+L1_mesh<=1.0
            b_mesh *= mask
            b_mesh = np.ma.masked_where(np.logical_not(mask), b_mesh)
            
            return x_mesh, y_mesh, b_mesh
        
        
    def getSimplexIntegralL0aL1bL2c_divA(self, a, b, c):
        """ get integral L0**a*L1**b*L2**c/A over triangle with area A
        """
        facts_args = np.array([a, b, c, a+b+c+2])
        facts = factorial(facts_args)
        I = 2.0*facts[0]*facts[1]*facts[2]/facts[3]
        return I


    def getOvIntd0Bid1Bj_divA(self, IJK0, IJK1, d0=(0, 0, 0), d1=(0, 0, 0)):
        """ dB(IJK0)/dd0*dB(IJK1)/dd1
            d0:d3/d0[0]d0[1]d0[2]
        """
        ##TODO: test
        if self.basisType==BasisTypes.nodeSimplex_2D:
            assert len(IJK0)==3 and len(IJK1)==3
            I0, J0, K0 = IJK0
            I1, J1, K1 = IJK1
            N = self.elemOrder
            if I0+J0+K0!=N or I1+J1+K1!=N:
                raise ValueError("I+J+K=N not satisfied.")
            
            #P_0 = self.B_NodeSimplex[I0]*self.B_NodeSimplex[J0]*self.B_NodeSimplex[K0] 
            #P_1 = self.B_NodeSimplex[I1]*self.B_NodeSimplex[J1]*self.B_NodeSimplex[K1] 
            
            dP_0 = [self.B_NodeSimplex[I0].deriv(d0[0]), self.B_NodeSimplex[J0].deriv(d0[1]), \
                self.B_NodeSimplex[K0].deriv(d0[2])]
            dP_1 = [self.B_NodeSimplex[I1].deriv(d1[0]), self.B_NodeSimplex[J1].deriv(d1[1]), \
                self.B_NodeSimplex[K1].deriv(d1[2])]
            
            dP_01 = [(dP_0[i]*dP_1[i]).coef for i in range(3)]
            
            getI_divA = self.getSimplexIntegralL0aL1bL2c_divA
            B_int = 0.0
            
            for i in range(len(dP_01[0])):
                for j in range(len(dP_01[1])):
                    for k in range(len(dP_01[2])):
                        B_int += getI_divA(dP_01[0][i], dP_01[1][j], dP_01[2][k])

            self.BB_int_divA[(IJK0, IJK1)] = B_int
            return B_int


    def getOvIntdBiP_divA(self, IJK, d, P):
        """ d3B(IJK)/dd*P
            d:d3/dd[0]dd[1]dd[2]
            P:P0(l0)*P1(l1)*P2(l2)
        """
        ##TODO: test
        if self.basisType==BasisTypes.nodeSimplex_2D:
            assert len(IJK)==3 and len(d)==3
            I, J, K = IJK
            N = self.elemOrder
            if I+J+K!=N:
                raise ValueError("I+J+K=N not satisfied.")
        
            dB = [self.B_NodeSimplex[I].deriv(d[0]), self.B_NodeSimplex[J].deriv(d[1]), \
                self.B_NodeSimplex[K].deriv(d[2])]
            
            dBP_c = [(dB[i]*P[i]).coef for i in range(3)]
            
            getI_divA = self.getSimplexIntegralL0aL1bL2c_divA
            BP_int = 0.0
            for i in range(len(dBP_c[0])):
                for j in range(len(dBP_c[1])):
                    for k in range(len(dBP_c[2])):
                        BP_int += getI_divA(dBP_c[0][i], dBP_c[1][j], dBP_c[2][k])

            self.BP_int_divA[(IJK, P)] = BP_int
            return BP_int
    
    @staticmethod        
    def mapNodeSimplex2DIJKto1D(N, ret="ijk-to-1d"):
        """ maps (I, J, K) to a unique number
            N: element order
        """
        IJK_to_N = {}
        IJK_arr = [(i, j, k) for i in range(N, -1, -1) for j in range(N, -1, -1) \
             for k in range(N, -1, -1) if i+j+k==N]
        #print(IJK_arr)
        
        if ret=="ijk-to-1d":
            return dict(zip(IJK_arr, range(len(IJK_arr))))
        elif ret=="1d-to-ijk":
            return dict(zip(range(len(IJK_arr)), IJK_arr))
        elif ret=="nodes-1d":
            ijk_1d = dict(zip(IJK_arr, range(len(IJK_arr))))
            nodes_ijk = [ijk for ijk in IJK_arr if np.sum(np.array(ijk)==0)==2]
            nodes_1d = [ijk_1d[ijk] for ijk in nodes_ijk]
            return nodes_1d
        elif ret=="sides-1d":
            ijk_1d = dict(zip(IJK_arr, range(len(IJK_arr))))
            sides_ijk = [ijk for ijk in IJK_arr if np.sum(np.array(ijk)==0)==1]
            sides_ijk_0 = [ijk for ijk in sides_ijk if ijk[0]==0]
            sides_ijk_1 = [ijk for ijk in sides_ijk if ijk[1]==0]
            sides_ijk_2 = [ijk for ijk in sides_ijk if ijk[2]==0]
            sides_1d_0 = [ijk_1d[ijk] for ijk in sides_ijk_0]
            sides_1d_1 = [ijk_1d[ijk] for ijk in sides_ijk_1]
            sides_1d_2 = [ijk_1d[ijk] for ijk in sides_ijk_2]
            sides_1d_1 = list(reversed(sides_1d_1))   ## counter clockwise
            return [sides_1d_0, sides_1d_1, sides_1d_2]
            
    
    @staticmethod
    def NodeSimplex2D_dxyTodL(dx_order, dy_order, r0=np.array([0.0, 0.0]), r1=np.array([1.0, 0.0]), \
            r2=np.array([0.0, 1.0])):
        """ dx_order: derivative order with respect to x
            dy_order: derivative order with respect to y
            dxdy ---> dL0dL1dL2
        """
        M_inv = FEMElemNode2D.NodeSimplex2D_CartesianToBarycentric(r0, r1, r2)
        print('M_inv = \n', M_inv)
        
        ## d/dx transformed to d/dL0 , d/dL1 , d/dL2
        _dx_ = {}
        if dx_order==0:
            _dx_ = {(0, 0, 0): 1.0}
        else:
            for i in range(dx_order):
                _dx_copy = copy.deepcopy(_dx_)
                if i==0:
                    _dx_copy = {(0, 0, 0): 0.0}
                for der_L012 in _dx_copy:
                    ## d/dL0
                    der0_der_L012 = np.array(der_L012, dtype=int)
                    der0_der_L012[0] += 1
                    der0_der_L012 = tuple(der0_der_L012)
                    if der0_der_L012 not in _dx_:
                        _dx_[der0_der_L012] = M_inv[0, 0]
                    else:
                        _dx_[der0_der_L012] += M_inv[0, 0]
                        
                    ## d/dL1
                    der1_der_L012 = np.array(der_L012, dtype=int)
                    der1_der_L012[1] += 1
                    der1_der_L012 = tuple(der1_der_L012)
                    if der1_der_L012 not in _dx_:
                        _dx_[der1_der_L012] = M_inv[1, 0]
                    else:
                        _dx_[der1_der_L012] += M_inv[1, 0]
                    
                    ## d/dL2
                    der2_der_L012 = np.array(der_L012, dtype=int)
                    der2_der_L012[2] += 1
                    der2_der_L012 = tuple(der2_der_L012)
                    if der2_der_L012 not in _dx_:
                        _dx_[der2_der_L012] = M_inv[2, 0]
                    else:
                        _dx_[der2_der_L012] += M_inv[2, 0]
        
        print("_dx_ : ", _dx_)
        
        _dy_ = {}
        if dy_order==0:
            _dy_ = {(0, 0, 0): 1.0}
        else:
            for i in range(dy_order):
                _dy_copy = copy.deepcopy(_dy_)
                if i==0:
                    _dy_copy = {(0, 0, 0): 0.0}
                for der_L012 in _dy_copy:
                    ## d/dL0
                    der0_der_L012 = np.array(der_L012, dtype=int)
                    der0_der_L012[0] += 1
                    der0_der_L012 = tuple(der0_der_L012)
                    if der0_der_L012 not in _dy_:
                        _dy_[der0_der_L012] = M_inv[0, 1]
                    else:
                        _dy_[der0_der_L012] += M_inv[0, 1]
                        
                    ## d/dL1
                    der1_der_L012 = np.array(der_L012, dtype=int)
                    der1_der_L012[1] += 1
                    der1_der_L012 = tuple(der1_der_L012)
                    if der1_der_L012 not in _dy_:
                        _dy_[der1_der_L012] = M_inv[1, 1]
                    else:
                        _dy_[der1_der_L012] += M_inv[1, 1]
                    
                    ## d/dL2
                    der2_der_L012 = np.array(der_L012, dtype=int)
                    der2_der_L012[2] += 1
                    der2_der_L012 = tuple(der2_der_L012)
                    if der2_der_L012 not in _dy_:
                        _dy_[der2_der_L012] = M_inv[2, 1]
                    else:
                        _dy_[der2_der_L012] += M_inv[2, 1]

        print("_dy_ : ", _dy_)
        
        _dx_dy_ = {}
        for der_L012__x in _dx_:
            for der_L012__y in _dy_:
                der_mul_xy =  np.array(der_L012__x, dtype=int) + np.array(der_L012__y, dtype=int)
                der_mul_xy = tuple(der_mul_xy)
                c_mul = _dx_[der_L012__x]*_dy_[der_L012__y]
                
                if der_mul_xy not in _dx_dy_:
                    _dx_dy_[der_mul_xy] = c_mul
                else:
                    _dx_dy_[der_mul_xy] += c_mul
                    
        
        return _dx_dy_
        

    @staticmethod
    def NodeSimplex2D_CartesianToBarycentric(r0=np.array([0.0, 0.0]), r1=np.array([1.0, 0.0]), \
            r2=np.array([0.0, 1.0])):
        
        M = np.array([[r0[0], r1[0], r2[0]],
                      [r0[1], r1[1], r2[1]],
                      [1,     1,     1    ]])
                      
        ## M*[L0, L1, L2]^T = [x, y, 1]^T   --->   [L0, L1, L2]^T = M^-1*[x, y, 1]^T
        return np.linalg.inv(M)
        
        
    def TransXYDerToXYDerOnMasterElem(self, m, n, r0_m=np.array([0.0, 0.0]), \
            r1_m=np.array([1.0, 0.0]), r2_m=np.array([0.0, 1.0])):
        """
           It transforms d/dx^m*d/dy^n derivatives from the triangle defined by self.r0, self.r1, self.r2
        to the equivalent derivatives on the master triangle r0_m, r1_m, r2_m.
        m: order of the x derivative
        n: order of the y derivative
        """
        ## find the linear transformation 
        ## x_p = a_p*x + b_p*y + c_p       (1)
        ## y_p = e_p*x + f_p*y + g_p       (2)
        ## that transforms self.r0, self.r1, self.r2 to the primed coordinate
        ## r0_m, r1_m, r2_m, respectively.
        x0, y0 = self.r0
        x1, y1 = self.r1
        x2, y2 = self.r2
        
        x0_p, y0_p = r0_m
        x1_p, y1_p = r1_m
        x2_p, y2_p = r2_m

        A = np.array([[x0, y0, 1],
                      [x1, y1, 1],
                      [x2, y2, 1]])
                      
        ## rhs for calculating transformation (1) and (2)
        b = np.array([[x0_p, x1_p, x2_p], [y0_p, y1_p, y2_p]]).T
                      
        ## solve for both transformations simultaneously
        abc_efg = np.linalg.solve(A, b)
        
        ## transformation coefficients
        abc_p = abc_efg[:,0]
        efg_p = abc_efg[:,1]
        
        ## differentials : d/dx**m in the prime system
        a_p, e_p = abc_p[0], efg_p[0]
        diffs_dx_m_p = {}
        binom_m = 1.0
        for i in range(m+1):
            diffs_dx_m_p[(m-i, i)] = binom_m * a_p**(m-i) * e_p**i
            binom_m *= (m-i)
            binom_m /= i+1

        ## differentials : d/dy**n in the prime system
        b_p, f_p = abc_p[1], efg_p[1]
        diffs_dy_n_p = {}
        binom_n = 1.0
        for i in range(n+1):
            diffs_dy_n_p[(n-i, i)] = binom_n * b_p**(n-i) * f_p**i
            binom_n *= (n-i)
            binom_n /= i+1

        ## multiply the two multivariate polynomials diffs_dx_m_p*diffs_dy_n_p
        diffs_dxdy_mn_p = {}
        for dx_dxdy_i in diffs_dx_m_p:
            for dy_dxdy_j in diffs_dy_n_p:
                dxy_combined = 
        
        ## treat m=0 or n=0 separately
                
            
        
        
                        

        
from meshpy.triangle import MeshInfo, build, write_gnuplot_mesh
from meshpy.geometry import GeometryBuilder, Marker, make_box

        
class MeshGen2D:

    def __init__(self):
        self.labels = {}
        self.labelsColors = {}
        return
        
    def AddLabel(self, label, color='b'):
        assert label not in self.labels
        self.labels[label] = Marker.FIRST_USER_MARKER+1+len(self.labels)
        self.labelsColors[label] = color
                
    def roundTripConnect(self, points):
        N = len(points)
        sides = []
        for i in range(N):
            sides.append((points[i], points[(i+1)%N]))
        return sides

    def SetDomain(self, nodes, boundary, holes, regions):
        """ boundary:{'points':[], 'sides':[], 'points_labels':[], 'sides_labels':[], 'v_max':v_max}
            holes:[ {'points':[], 'sides':[], 'points_labels':[], 'sides_labels':[], 'point_inside':(x,y)} ]
            regions:[ {'points':[], 'sides':[], 'points_labels':[], 'sides_labels':[], 'region_label':xx, 'point_inside':(x,y), 'v_max':v_max} ]
        """
        self.nodes = nodes
        self.boundary = boundary
        self.holes = holes
        self.regions = regions
        return
        
    
    def BuildMesh(self):
        """
        """
        points = self.nodes
        points_mk = [0]*len(points)
        pts = self.boundary['points']
        pts_labels = self.boundary['points_labels']
        for j in range(len(pts)):
            points_mk[pts[j]] = pts_labels[j]
        for i in range(len(self.holes)):
            pts = self.holes[i]['points']
            pts_labels = self.holes[i]['points_labels']
            for j in range(len(pts)):
                points_mk[pts[j]] = pts_labels[j]
        for i in range(len(self.regions)):
            pts = self.regions[i]['points']
            pts_labels = self.regions[i]['points_labels']
            for j in range(len(pts)):
                points_mk[pts[j]] = pts_labels[j]

        facets = self.boundary['sides']
        facets_mk = self.boundary['sides_labels']
        for i in range(len(self.holes)):
            facets.extend(self.holes[i]['sides'])
            facets_mk.extend(self.holes[i]['sides_labels'])
        for i in range(len(self.regions)):
            facets.extend(self.regions[i]['sides'])
            facets_mk.extend(self.regions[i]['sides_labels'])

        points_holes = []
        for i in range(len(self.holes)):
            points_holes.append(self.holes[i]['point_inside'])
        

        mesh_info = MeshInfo()
        mesh_info.set_points(points, points_mk)
        mesh_info.set_facets(facets, facets_mk)
        mesh_info.set_holes(points_holes)

        if len(self.regions):
            mesh_info.regions.resize(len(self.regions))
            for i in range(len(self.regions)):
                x, y = self.regions[i]['point_inside']
                mesh_info.regions[i] = (x, y, self.regions[i]['region_label'], self.regions[i]['v_max'])

        self.mesh = build(mesh_info, max_volume=self.boundary['v_max'], attributes=True, 
                            volume_constraints=True, generate_faces=True)
        return self.mesh
        

    def PlotMesh(self, showlabels=['nodes', 'sides', 'regions']):
        import matplotlib.pyplot as plt
        mesh = self.mesh
        labelsVals = self.labels.values()
        valToLabel = {}
        for key, val in self.labels.items():
            valToLabel[val] = key
        plt.triplot(np.array(mesh.points)[:, 0], np.array(mesh.points)[:, 1], np.array(mesh.elements))
        if 'nodes' in showlabels:
            for i in range(len(mesh.points)):
                if mesh.point_markers[i] in labelsVals:
                    plt.plot(mesh.points[i][0], mesh.points[i][1], self.labelsColors[valToLabel[mesh.point_markers[i]]]+'o')
        
        if 'sides' in showlabels:
            for i in range(len(mesh.facets)):
                if mesh.facet_markers[i] in labelsVals:
                    p0, p1 = mesh.facets[i]
                    plt.plot([mesh.points[p0][0], mesh.points[p1][0]], [mesh.points[p0][1], mesh.points[p1][1]], \
                        self.labelsColors[valToLabel[mesh.facet_markers[i]]], lw=2)

        if 'regions' in showlabels:
            for ind in range(len(self.regions)):
                reg_elems = [mesh.elements[i] for i in range(len(mesh.elements)) if mesh.element_attributes[i]==self.regions[ind]['region_label']]  
                plt.triplot(np.array(mesh.points)[:, 0], np.array(mesh.points)[:, 1], np.array(reg_elems), \
                    color=self.labelsColors[valToLabel[self.regions[ind]['region_label']]])

        plt.show()




class FEM2D:
    
    def __init__(self):
        return
        
    
    def SetStruct(self, mesh):
        self.mesh = mesh
        ##TODO: confirm elements nodes are counter-clockwise
        
        self.nodes = np.array(mesh.points)
        self.elems = np.array(mesh.elements)
        self.sides = np.array(mesh.faces)

        self.nodelabels = mesh.point_markers
        self.sidelabels = mesh.facet_markers
        self.elemlabels = mesh.element_attributes
        
        self.SetLabelsAssociatedSidesOrNodesOrElems()
        self.SetNodesConnectedSides()
        self.SetElemsConnectedSides()
        
        self.assemblyParams = []
        return
        
                
    def SetLabelsAssociatedSidesOrNodesOrElems(self):
        labels_sides = {}
        n_sides = len(self.sides)
        if self.sidelabels is None:
            labels_sides[None] = np.arange(n_sides)
        else:
            labelsAll = set(self.sidelabels)
            for label in labelsAll:
                labels_sides[label] = np.where(self.sidelabels==label)[0]
        self.labels_sides = labels_sides
            
        labels_nodes = {}
        n_nodes = len(self.nodes)
        if self.nodelabels is None:
            labels_nodes[None] = np.arange(n_nodes)
        else:
            labelsAll = set(self.nodelabels)
            for label in labelsAll:
                labels_nodes[label] = np.where(self.nodelabels==label)[0]
        self.labels_nodes = labels_nodes

        labels_elems = {}
        n_elems = len(self.elems)
        if self.elemlabels is None:
            labels_elems[None] = np.arange(n_elems)
        else:
            labelsAll = set(self.elemlabels)
            for label in labelsAll:
                labels_elems[label] = np.where(self.elemlabels==label)[0]
        self.labels_elems = labels_elems
        return


    def SetNodesConnectedSides(self):
        n_nodes = len(self.nodes)
        n_sides = len(self.sides)
        nodes_connectedSides = [None]*n_nodes
        #print(self.sides)
        for i in range(n_nodes):
            nodes_connectedSides[i] = []
        for i in range(n_sides):
            n0 = self.sides[i, 0]
            if i not in nodes_connectedSides[n0]:
                nodes_connectedSides[n0].append(i)
            n1 = self.sides[i, 1]
            if i not in nodes_connectedSides[n1]:
                nodes_connectedSides[n1].append(i)
        self.nodes_connectedSides = nodes_connectedSides
        #print(self.nodes_connectedSides)
        return


    def SetNodesConnectedElems(self):
        n_nodes = len(self.nodes)
        n_elems = len(self.elems)
        nodes_connectedElems = [None]*n_nodes
        for i in range(n_sides):
            nodes_connectedElems[i] = []
        for i in range(n_elems):
            n0 = self.elems[i, 0]
            if i not in nodes_connectedElems[n0]:
                nodes_connectedElems[n0].append(i)
            n1 = self.elems[i, 1]
            if i not in nodes_connectedElems[n1]:
                nodes_connectedElems[n1].append(i)
            n2 = self.elems[i, 2]
            if i not in nodes_connectedElems[n2]:
                nodes_connectedElems[n2].append(i)
        self.nodes_connectedElems = nodes_connectedElems
        return


    def SetElemsConnectedSides(self):
        n_elems = len(self.elems)
        n_sides = len(self.sides)
        elems_connectedSides = -np.ones((n_elems, 3), dtype=int)
        sides_connectedElems = [None]*n_sides
        for i in range(n_sides):
            sides_connectedElems[i] = []
        nodes_connectedSides = self.nodes_connectedSides
        for i in range(n_elems):
            n0, n1, n2 = self.elems[i, :]
            n0_s = set(nodes_connectedSides[n0])
            n1_s = set(nodes_connectedSides[n1])
            n2_s = set(nodes_connectedSides[n2])
            s0 = list(n1_s & n2_s)
            s1 = list(n0_s & n2_s)
            s2 = list(n0_s & n1_s)
            assert len(s0)==1 and len(s1)==1 and len(s2)==1
            s0, s1, s2 = s0[0], s1[0], s2[0]
            elems_connectedSides[i, 0] = s0
            elems_connectedSides[i, 1] = s1
            elems_connectedSides[i, 2] = s2
            sides_connectedElems[s0].append(i)
            sides_connectedElems[s1].append(i)
            sides_connectedElems[s2].append(i)
        self.elems_connectedSides = elems_connectedSides
        self.sides_connectedElems = sides_connectedElems
        return elems_connectedSides, sides_connectedElems

        
    def SetGlobalIndices(self, elem_order, elem_type):
        self.elem_order = elem_order
        self.elem_type = elem_type
        if elem_type==BasisTypes.nodeSimplex_2D:
            n_nodes = len(self.nodes)
            n_sides = len(self.sides)
            n_elems = len(self.elems)
            
            n_base_elem = (elem_order+1)*(elem_order+2)//2
            n_b_side = (elem_order-1)
            n_b_inside = (elem_order-2)*(elem_order-1)//2

            n_base_T = n_elems*n_b_inside + n_sides*n_b_side + n_nodes
            print("n_base_T : ", n_base_T)

            elemsBases = -np.ones((n_elems, n_base_elem), dtype=int)
            sidesBases = -np.ones((n_sides, n_b_side), dtype=int)
            nodesBases = -np.ones(n_nodes, dtype=int)
                        
            
            ## first set nodes, then set bases on side, finally
            ## set bases inside  
            ind_last = 0
            for i in range(n_elems):
                ns = self.elems[i, :]
                assert len(ns)==3
                for j in range(3):
                    nj = ns[j]
                    if nodesBases[nj]<0:
                        elemsBases[i, j] = ind_last
                        nodesBases[nj] = ind_last
                        ind_last += 1
                    else:
                        elemsBases[i, j] = nodesBases[nj]
                
            if n_b_side>0:
                _i_s_start = 3
                elems_connectedSides = self.elems_connectedSides
                for i in range(n_elems):
                    ns = self.elems[i, :]
                    ss = elems_connectedSides[i,:]
                    for j in range(3):
                        sj = ss[j]
                        sj_n0, sj_n1 = self.sides[sj, :]
                        e_n0, e_n1 = ns[(j+1)%3], ns[(j+2)%3]
                        ccw = True
                        assert abs(sj_n1-sj_n0)==abs(e_n1-e_n0)
                        if (sj_n1-sj_n0)*(e_n1-e_n0)<0:
                            ccw = False
                        if sidesBases[sj,0]<0:
                            for k in range(n_b_side):
                                elemsBases[i, _i_s_start + j*n_b_side + k] = ind_last
                                if ccw:
                                    sidesBases[sj,k] = ind_last
                                else:
                                    sidesBases[sj,n_b_side-1-k] = ind_last
                                ind_last += 1
                        else:
                            for k in range(n_b_side):
                                if ccw:
                                    elemsBases[i, _i_s_start + j*n_b_side + k] = sidesBases[sj,k]
                                else:
                                    elemsBases[i, _i_s_start + j*n_b_side + k] = sidesBases[sj,n_b_side-1-k]
                            
                
            if n_b_inside>0:
                _i_inside_start = _i_s_start + 3*n_b_side
                #print('_i_inside_start: ', _i_inside_start, _i_inside_start+n_b_inside)
                for i in range(n_elems):
                    for k in range(n_b_inside):
                        elemsBases[i, _i_inside_start + k] = ind_last
                        ind_last += 1
            
            ##TODO: reorder based on element closeness
            
            self.n_base_T = n_base_T
            self.elemsBases = elemsBases
            return elemsBases
        else:
            raise NotImplementedError()
        return

    def DefineMatrixAssemblyProcess(self, description, label, proc, procParams):
        """ example:description: 'mat'/'rhs'/'bc'
                    proc=AssemblyProcess.galerkin
                    procParams = {'coeffType':CoeffTypes.constant, 'C':2.0, 'dx_Ui':1, 'dy_Ui':0, 'dx_Uj':1, 'dy_Uj':0}
            for integral over cell 2*dw_dx*du_dx, (w is the testing and u is the basis function)  
            i.e. C=2.0, derivative x order w=1, derivative x order u=1
            w=wighting function ---> U_i
            u=basis function  -----> U_j
            
            'mat' ---> stifness matrix      label --> elementlabel
            'rhs' ---> rhs                  label --> elementlabel
            'bc'  ---> boundary condition   label --> sidelabel
            
            for 'bc'
            procParams = {'bcType':BCTypes.Dirichlet, 'bcParams':bcParams}
                bcParams:
                    Dirichlet : 'rhs'   ---->  u(0) = rhs
                    Neumann   : 'rhs'   ---->  u'(0) = rhs
                    Robin     : 'sigma', 'rhs' ---> u'(0)-sigma*u(0)=rhs
        """
        assert description in ['mat', 'rhs', 'bc']
        
        self.assemblyParams.append({'description':description, 'label':label, \
                'proc':proc, 'procParams':procParams}) 


    def AssembleMatrix(self):
        row, col, data = np.array([], dtype=int), np.array([], dtype=int), np.array([])
        row_rhs, data_rhs = np.array([], dtype=int), np.array([])
        row_preset, data_preset = np.array([], dtype=int), np.array([])
        for assproc in self.assemblyParams:
            description = assproc['description']
            proc = assproc['proc']
            if description=='mat' and proc==AssemblyProcess.galerkin:
                label = assproc['label']
                procParams = assproc['procParams']
                
                assert label in self.labels_elems
                lab_elems = self.labels_elems[label]
                row_, col_, data_ = self.GetGalerkinMatElems(lab_elems, procParams)
                
                row = np.append(row, row_)
                col = np.append(col, col_)
                data = np.append(data, data_)
                
            elif description=='rhs' and proc==AssemblyProcess.galerkin:
                label = assproc['label']
                procParams = assproc['procParams']

                assert label in self.labels_elems
                lab_elems = self.labels_elems[label]

                row_rhs_, data_rhs_ = self.GetGalerkinRhs(lab_elems, procParams)

                row_rhs = np.append(row_rhs, row_rhs_)
                data_rhs = np.append(data_rhs, data_rhs_)
                
            elif description=='bc' and proc==AssemblyProcess.galerkin:
                label = assproc['label']
                procParams = assproc['procParams']

                assert label in self.labels_sides
                lab_sides = self.labels_sides[label]

                bcType = procParams['bcType']
                bcParams = procParams['bcParams']
                
                if bcType==BCTypes.Dirichlet:
                    row_ps, data_ps = self.GetGalerkinDirichlet(lab_sides, bcParams)

                    row_preset = np.append(row_preset, row_ps)
                    data_preset = np.append(data_preset, data_ps)
                    
                elif bcType==BCTypes.Neumann:
                    raise NotImplementedError()
                elif bcType==BCTypes.Robin:
                    raise NotImplementedError()
                else:
                    raise NotImplementedError()
                
            else:
                raise NotImplementedError()
            
        #print(row, col, data, sep='\n')
        
        ##reset presets
        mask = np.ones(len(row), dtype=bool)
        if len(row_preset)>0:
            mask = np.logical_not(row==row_preset[0])
        for i in range(1, len(row_preset)):
            mask = np.logical_and(mask, np.logical_not(row==row_preset[i]))
        
        row = row[mask]
        col = col[mask]
        data = data[mask]
        
        row = np.append(row, row_preset)
        col = np.append(col, row_preset)
        data = np.append(data, np.ones(len(row_preset)))
        
        mask = np.ones(len(row), dtype=bool)
        if len(row_preset)>0:
            mask_rhs = np.logical_not(row_rhs==row_preset[0])
        for i in range(1, len(row_preset)):
            mask_rhs = np.logical_and(mask_rhs, np.logical_not(row_rhs==row_preset[i]))
        
        row_rhs = row_rhs[mask_rhs]
        data_rhs = data_rhs[mask_rhs]
        
        row_rhs = np.append(row_rhs, row_preset)
        data_rhs = np.append(data_rhs, data_preset)
        
        ## generate sparse matrix, solve
        A_coo = scipy.sparse.coo_matrix((data, (row, col)), shape=(self.n_base_T, self.n_base_T))
        #print('A_coo : \n', A_coo)
        
        b_vec = np.zeros(self.n_base_T)
        for i in range(len(row_rhs)):
            b_vec[row_rhs[i]] += data_rhs[i]
            
        from scipy.sparse import linalg
        A_csc = A_coo.tocsc()
        x_res = linalg.spsolve(A_csc, b_vec)
        
        self.x_res = x_res
        return x_res



    def GetGalerkinMatElems(self, elems_inds, params):
        elems = self.elems
        sides = self.sides
        nodes = self.nodes
        dx = self.dx
        elemsBases = self.elemsBases
        if params['coeffType']==CoeffTypes.constant:
            coeff = params['C']
            dxUi = params['dx_Ui']
            dxUj = params['dx_Uj']
            dyUi = params['dy_Ui']
            dyUj = params['dy_Uj']
            
            elem_order = self.elem_order
            fem_elem = FEMElemNode2D(elem_order, basisType=self.elem_type, \
                r0=np.array([0.0, 0.0]), r1=np.array([1.0, 0.0]), r2=np.array([0.0, 1.0]))
                
            dLUi = NodeSimplex2D_dxyTodL(dx_order=dxUi, dy_order=dyUi, r0=np.array([0.0, 0.0]), \
                                        r1=np.array([1.0, 0.0]), r2=np.array([0.0, 1.0]))
                                        
            ### to complete ???
            
            ovInts = fem_elem.getOvIntd0Bid1Bj_divA(IJK0, IJK1, d0=(0, 0, 0), d1=(0, 0, 0))
            
            row, col, data = [], [], []
            for s in sides_inds:
                mat_s = ovInts/dx[s]**(dUi+dUj-1)
                shape = mat_s.shape
                row_s = np.zeros(shape, dtype=int)
                col_s = np.zeros(shape, dtype=int)
                for i in range(shape[0]):
                    row_s[i,:] = np.ones(shape[1], dtype=int)*sidesBases[s, i]
                    col_s[:,i] = np.ones(shape[1], dtype=int)*sidesBases[s, i]
                
                #print(row_s, col_s, sep='\n')
                
                data_s = mat_s.reshape(shape[0]*shape[1])
                row_s = row_s.reshape(shape[0]*shape[1])
                col_s = col_s.reshape(shape[0]*shape[1])
                
                nnz_s = data_s!=0.0
                
                row_s = row_s[nnz_s]
                col_s = col_s[nnz_s]
                data_s = data_s[nnz_s]
                
                row.extend(row_s)
                col.extend(col_s)
                data.extend(data_s)
                
            return np.array(row), np.array(col), coeff*np.array(data)

        else:
            raise NotImplementedError()
                    
        







