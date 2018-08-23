## RecGridND.py
##

__all__ = ["RGND", "GUIMakerND", "ParType"]

import numpy as np
import scipy as sp
import random
import math

from enum import Enum, IntEnum

from Electromagnetics.SymExprTree import symExpr_generate_tree, symExpr_update_tree
import copy
from scipy.sparse import coo_matrix
from scipy.special import binom
from scipy.interpolate import lagrange

from IPython.display import display, Math
from sympy import latex, lambdify, sympify

from collections import deque

##TODO: test sympy .coeff function


class CoeffType(Enum):
    const = 0
    par_int = 1
    nonconst_anal = 2
    par_ext = 3

class DebugMode(IntEnum):
    lev_0 = 0
    lev_1 = 1
    lev_2 = 2
    lev_3 = 3
    
class ParType(IntEnum):
    general = 0
    seperable = 1
    partialseperable = 2
    polyarg = 3


NEXIST = -1
IND_NEXIST = -1
N_level_max = 30
H_CN = 0 ## cell nodes index inside cellsHier
H_CP = 1 ## cell parent index inside cellsHier
H_CC = 2 ## cells child index inside cellsHier
H_CI = 3 ## cell index index inside cellsHier

class RGND:
    """ N-Dimensional structured rectangular grid - each cell can be refined on demand
    - the adjacent cells either are the same size, or one is half the other, if
    this condition does not hold cells are reifined automatically until this 
    condition is satisfied
    """
    def __init__(self, x0, x1, dx, rg=None):
        """ [x0, x1]: intervals
            x0: ND numpy array
            x1: numpy array
            dx: numpy array
            rg: copy grid properties from another RGND 
        """
        if rg==None:
            nx = np.ceil((x1 - x0)/dx).astype(int)
            dx = (x1 - x0)/nx
            
            self.nx = nx
            
            N = x0.shape[0]
            assert x1.shape[0]==N and dx.shape[0]==N
            self.N = N
            self.W = x1 - x0
            self.x0, self.x1 = np.copy(x0), np.copy(x1)
            assert np.all(x0<x1)
            
            self.dx_levels = [None]*N_level_max
            self.dx_levels[0] = dx
            for i in range(1, N_level_max):
                self.dx_levels[i] = self.dx_levels[i-1]/2.0
            
            n_cell = 1
            n_node = 1
            for i in range(N):
                n_cell *= nx[i]
                n_node *= nx[i]+1
            print('--', nx, n_cell)
            self.cellsHier = [None]
            self.cellsHier[0] = [None]*n_cell    #cells in a hierarchical order
            self.nodesPoints = [None]*n_node

            counter = np.zeros(nx.shape, dtype=int)
            
            nx_pow = np.ones(N, dtype=int)      ## nx**i
            for i in range(N-2, -1, -1):
                nx_pow[i] = nx[i+1]*nx_pow[i+1]
            self.nx_pow = nx_pow

            nxp1_pow = np.ones(N, dtype=int)    ## (nx+1)**i
            for i in range(N-2, -1, -1):
                nxp1_pow[i] = (nx[i+1]+1)*nxp1_pow[i+1]
            self.nxp1_pow = nxp1_pow
            
            ## cells on the boundary
            self.cellsWithNodesOnBorder = [{'n':[None]*N, 'p':[None]*N}]  # [neg[Dir]:cells, pos[Dir]:cells] : for each level
            for i in range(N):
                self.cellsWithNodesOnBorder[0]['n'][i] = []
                self.cellsWithNodesOnBorder[0]['p'][i] = []
            cellsWithNodesOnBorder_0_n = self.cellsWithNodesOnBorder[0]['n']
            cellsWithNodesOnBorder_0_p = self.cellsWithNodesOnBorder[0]['p']

            while True:
                ## cell :[nodes, cell_parent, cells_child]
                cell_ind = self._uniformGridCellInd(counter, nx, nx_pow)
                cell_nodes = self._uniformGridGetCellNodes(counter, nx, nxp1_pow)
                cell_parent = NEXIST
                cells_child = NEXIST
                cell_ind_tot = IND_NEXIST
                self.cellsHier[0][cell_ind] = [cell_nodes, cell_parent, cells_child, cell_ind_tot]
                for i in range(N):
                    if counter[i]==0:
                        cellsWithNodesOnBorder_0_n[i].append(cell_ind)
                    elif counter[i]==nx[i]-1:
                        cellsWithNodesOnBorder_0_p[i].append(cell_ind)
                #print(counter, cell_ind)
                if not self._increaseCellCounterIndex(counter, nx):
                    break
            
            counter = np.zeros(nx.shape, dtype=int)
            while True:
                node_ind = self._uniformGridNodeInd(counter, nx, nxp1_pow)
                self.nodesPoints[node_ind] = counter*dx
                if not self._increaseNodeCounterIndex(counter, nx):
                    break
        else:
            self.nx = rg.nx
            self.N = rg.N
            self.W = rg.W
            ##x0 and x1 (reference points) could be different but x1-x0 should be the same
            self.x0, self.x1 = x0, x1
            assert np.all(np.abs(self.x1-self.x0-self.W)/self.W<1.0e-15)
            assert np.all(self.x0<self.x1)
            
            self.dx_levels = rg.dx_levels
            self.cellsHier = rg.cellsHier
            self.nodesPoints = rg.nodesPoints
            self.nx_pow = rg.nx_pow
            self.nxp1_pow = rg.nxp1_pow
            self.cellsWithNodesOnBorder = rg.cellsWithNodesOnBorder
            
        
        self.BCs = None
        self.CCs = None
        self.cells_nb_dic = None
        self.cc_dic = None
        self.cc_eqindex_and_multiplicity = None
        self.cc_nbmul_dic = None
        self.bctoic_map = None
        self.bc_dic = None
        self.bc_eqindex_and_multiplicity = None
        self.cellsCC = None
        self.bc_inds = None
        self.N_indexed_acc = None
        self.pf_dic = None
        
        self.pars_extern = None
        
        self.nx_arr = None
        self.P0P1_M = None
        
        self.verbose = 0
        
        self.FieldsAtSurf = {}
        self.que_maxlen = 10
        
        self.matsolver='SLU'
        return
            
    def CloneGridAndIndexing(self, rg):
        self.nx = rg.nx
        self.N = rg.N
        self.W = rg.W
        ##x0 and x1 (reference points) could be different but x1-x0 should be the same
        assert np.all(np.abs(self.x1-self.x0-self.W)/self.W<1.0e-15)
        assert np.all(self.x0<self.x1)
        
        self.dx_levels = rg.dx_levels
        self.cellsHier = rg.cellsHier
        self.nodesPoints = rg.nodesPoints
        self.nx_pow = rg.nx_pow
        self.nxp1_pow = rg.nxp1_pow
        self.cellsWithNodesOnBorder = rg.cellsWithNodesOnBorder
        return

    def _uniformGridCellInd(self, i, nx, nx_pow):
        #print(i, nx)
        if np.all((0<=i)*(i<nx)):
            return np.sum(i*nx_pow)
        else:
            return NEXIST
        

    def _uniformGridNodeInd(self, i, nx, nxp1_pow):
        if np.all((0<=i)*(i<=nx)):
            return np.sum(i*nxp1_pow)
        else:
            return NEXIST
        
    def _uniformGridGetCellNodes(self, i, nx, nxp1_pow):
        ## i: (i1, i2, i3...) i-n:cell index in the n-th dimension
        N = self.N
        di = self._getBinaryCounter()
        cell_nodes = [-1]*2**N
        ind = 0
        while True:
            cell_nodes[ind] = self._uniformGridNodeInd(i+di, nx, nxp1_pow)
            if not self._binaryCounterIncrease(di):
                break
            ind += 1
        for i in range(len(cell_nodes)):
            assert cell_nodes[i] >= 0
        return cell_nodes
        
    def _getBinaryCounter(self):
        return np.zeros(self.N, dtype=int)

    def _binaryCounterIncrease(self, counter):
        for i in range(self.N-1, -1, -1):
            if counter[i]==0:
                counter[i] = 1
                return True
            else:
                counter[i] = 0
        return False

    def _binaryCounterIncrease_Masked(self, counter, inds_mask):
        for i in range(self.N-1, -1, -1):
            if i not in inds_mask:
                if counter[i]==0:
                    counter[i] = 1
                    return True
                else:
                    counter[i] = 0
        return False
        
        
    def getBinaryCounterValue(self, counter):
        """ counter[0] --> 2**(N-1)
            counter[N-1]--> 2**0
        """
        val = 0
        _2_pn = 1
        N = len(counter)
        for i in range(N-1, -1, -1):
            val += counter[i]*_2_pn
            _2_pn *= 2
        return val
        

    def getBinaryCounterValue_Masked(self, counter, inds_mask):
        """ indices coressponding to inds_mask are removed befor calculating the
        binary value
        """
        val = 0
        _2_pn = 1
        N = len(counter)
        for i in range(N-1, -1, -1):
            if i not in inds_mask:
                val += counter[i]*_2_pn
                _2_pn *= 2
        return val


    def _increaseCellCounterIndex(self, counter, nx):
        for i in range(self.N-1, -1, -1):
            if counter[i]<nx[i]-1:
                counter[i] += 1
                return True
            else:
                counter[i] = 0
        return False
        

    def _increaseCellCounterIndex_Masked(self, counter, nx, inds_mask):
        for i in range(self.N-1, -1, -1):
            if i not in inds_mask:
                if counter[i]<nx[i]-1:
                    counter[i] += 1
                    return True
                else:
                    counter[i] = 0
        return False


    def _increaseNodeCounterIndex(self, counter, nx):
        for i in range(self.N-1, -1, -1):
            if counter[i]<nx[i]:
                counter[i] += 1
                return True
            else:
                counter[i] = 0
        return False

        
    def _increasePolyCounterIndex(self, counter, p_ord):
        for i in range(self.N-1, -1, -1):
            if counter[i]<p_ord[i]:
                counter[i] += 1
                return True
            else:
                counter[i] = 0
        return False
        
        
    def _increasePolyCounterIndex_Masked(self, counter, p_ord, inds_mask):
        for i in range(self.N-1, -1, -1):
            if i not in inds_mask:
                if counter[i]<p_ord[i]:
                    counter[i] += 1
                    return True
                else:
                    counter[i] = 0
        return False
        

    def getPolyCounterValue(self, counter, p_ord):
        """ counter[0] --> (P_ord[1]+1)*(P_ord[2]+1)..
            counter[N-1]--> (p_ord[N-1]+1)**0
        """
        val = 0
        _pi_pim1_ = 1
        N = len(counter)
        for i in range(N-1, -1, -1):
            val += counter[i]*_pi_pim1_
            _pi_pim1_ *= (p_ord[i]+1)
        return val
        
        
    def CellGetSubchilrenInDirection(self, cell_lev, cell_ind, dir_ind, dir_sign='n'):
        """ dir_ind: 0, 1...   0-->x   1-->y   2-->z ...
            dir_sign : 'n'/'p'    'n'-->negative   'p'-->positive
        """ 
        N = self.N
        counter = np.zeros(N, dtype=int)        
        if dir_sign=='p':
            counter[dir_ind] = 1
        inds_mask = [dir_ind]
        cell_children = self.cellsHier[cell_lev][cell_ind][H_CC]
        if cell_children==NEXIST:
            return None
        sub_cells = []
        while True:
            child_ind = self.getBinaryCounterValue(counter)
            sub_cells.append(cell_children[child_ind])
            if not self._binaryCounterIncrease_Masked(counter, inds_mask):
                break
        return sub_cells
        

    def NodesGetConnectedCells(self):
        """ it returns a list specifying the cells connected to each node and 
        their levels
        returns node_conn_cells[lev][node] = [connectedcells] 
        lev: level
        node: node_index
        the node_conn_cells should be updated if the grid is modified in any ways
        """
        # n_ind: node index
        N_nodes = len(self.nodesPoints)
        N_levels = len(self.cellsHier)
        node_conn_cells = [None]*N_nodes
        for n in range(N_nodes):
            node_conn_cells[n] = {}
        
        for lev in range(N_levels):
            cells_lev = self.cellsHier[lev]
            for c_ind in range(len(cells_lev)):
                c_nodes = cells_lev[c_ind][H_CN]    ## cell nodes
                for i in range(len(c_nodes)):
                    node = c_nodes[i]
                    if lev in node_conn_cells[node]:
                        node_conn_cells[node][lev].append(c_ind)
                    else:
                        node_conn_cells[node][lev] = [c_ind]
        self.node_conn_cells = node_conn_cells
        return node_conn_cells
        
        
    def NodesGetConnectedLeafCells(self):
        """ assuming self.node_conn_cells is already constructed and up to date
        it returns a copy of it that includes only the leaf cells (cells withought
         children)
        """
        node_conn_cells = self.node_conn_cells
        cellsHier = self.cellsHier
        N_nodes = len(self.nodesPoints)
        N_levels = len(self.cellsHier)

        #count number of leaf connected cells
        n_node_conn_cells_leaf = [None]*N_nodes
        for n in range(N_nodes):
            n_node_conn_cells_leaf[n] = {}

        for n in range(N_nodes):
            for lev in node_conn_cells[n]:
                node_conn_cells_n_lev = node_conn_cells[n][lev]
                n_node_conn_cells_leaf_n_lev = 0
                cellsHier_lev = cellsHier[lev]
                for c in node_conn_cells_n_lev:
                    if cellsHier_lev[c][H_CC]==NEXIST:
                        n_node_conn_cells_leaf_n_lev += 1
                n_node_conn_cells_leaf[n][lev] = n_node_conn_cells_leaf_n_lev
                
        #set up list
        node_conn_cells_leaf = [None]*N_nodes
        for n in range(N_nodes):
            node_conn_cells_leaf[n] = {}

        for n in range(N_nodes):
            for lev in node_conn_cells[n]:
                node_conn_cells_n_lev = node_conn_cells[n][lev]
                n_cell = n_node_conn_cells_leaf[n][lev]
                if n_cell>0:
                    node_conn_cells_leaf[n][lev] = [-1]*n_cell
                    node_conn_cells_leaf_n_lev = node_conn_cells_leaf[n][lev]
                    cellsHier_lev = cellsHier[lev]
                    _ind = 0
                    for c in node_conn_cells_n_lev:
                        if cellsHier_lev[c][H_CC]==NEXIST:
                            node_conn_cells_leaf_n_lev[_ind] = c
                            _ind += 1
        
        return node_conn_cells_leaf
    

    def CellGetCellsAroundSameNexLev(self, c_lev, c_ind):
        """ it returns the neigboring cells in the same and next level
        The cell itself is included        
        
        c_lev: cell level
        c_ind : cell index
        """
        c_nodes = self.cellsHier[c_lev][c_ind][H_CN]
        cells_neighb_set_SL = set()     ## SL: same level
        for i in range(len(c_nodes)):
            cells_nb_i = self.node_conn_cells[c_nodes[i]][c_lev]
            for j in range(len(cells_nb_i)):
                cells_neighb_set_SL.add(cells_nb_i[j])
        N_nb = len(cells_neighb_set_SL)
        N_levels = len(self.cellsHier)

        """
        cells_neighb_set_NL = set()     ## NL: Next Level
        for c_nb in cells_neighb_set_SL:
            c_childs = self.cellsHier[c_lev][c_nb][H_CC]
            if c_childs!=NEXIST:
                assert c_lev<N_levels-1
                for c_ch in c_childs:
                    cells_neighb_set_NL.add(c_ch)
                    
        """

        ##next level connected cells
        cells_neighb_set_NL = set()     ## NL: next level
        if len(self.cellsHier)>c_lev+1:
            for i in range(len(c_nodes)):
                #print(self.node_conn_cells[c_nodes[i]])
                if c_lev+1 in self.node_conn_cells[c_nodes[i]]: 
                    cells_nb_i = self.node_conn_cells[c_nodes[i]][c_lev+1]
                    for j in range(len(cells_nb_i)):
                        cells_neighb_set_NL.add(cells_nb_i[j])
        
        return [cells_neighb_set_SL, cells_neighb_set_NL]
        
    
    def CellGetFaceConnectedCellsSameLevel(self, c_lev, c_ind):
        """ it returns neighbor cells which have a face in common
        """
        c_nodes = self.cellsHier[c_lev][c_ind][H_CN]
        cells_neighb_set_SL = set()     ## SL: same level
        for i in range(len(c_nodes)):
            cells_nb_i = self.node_conn_cells[c_nodes[i]][c_lev]
            for j in range(len(cells_nb_i)):
                cells_neighb_set_SL.add(cells_nb_i[j])
        
        cells_neighb_set_SL.remove(c_ind)
        
        ##check if there is 2**(N-1) nodes in common
        N = self.N
        dx_lev = self.dx_levels[c_lev]
        nodesPoints = self.nodesPoints
        node_corner = nodesPoints[c_nodes[0]]
        cells_nb_neg = [None]*N
        cells_nb_pos = [None]*N
        for cells_nb_i in cells_neighb_set_SL:
            nodes_nb = self.cellsHier[c_lev][cells_nb_i][H_CN]
            nodes_comm = []
            for n_nb in nodes_nb:
                for n in c_nodes:
                    if n_nb==n:
                        nodes_comm.append(n)
                        break
            if len(nodes_comm)==2**(N-1):
                ##common face: find orientation
                node_ref = nodesPoints[nodes_comm[0]]
                Dir = None  ##direction
                for i in range(N):
                    Dir = i
                    for j in range(1, len(nodes_comm)):
                        if abs(nodesPoints[nodes_comm[j]][i] - node_ref[i]) > dx_lev[i]/10:
                            Dir = None
                            break
                    if Dir!=None:
                        break
                if node_ref[Dir]-node_corner[Dir]>dx_lev[Dir]/10:
                    cells_nb_pos[Dir] = cells_nb_i
                else:
                    cells_nb_neg[Dir] = cells_nb_i
        
        return [cells_nb_neg, cells_nb_pos]
        

    def CellGetFaceConnectedCellsSameNextLevel(self, c_lev, c_ind, getConnNodes=False):
        """ it returns neighbor cells which have a face in common
        c_lev: current cell level
        c_ind: cell index
        """
        c_nodes = self.cellsHier[c_lev][c_ind][H_CN]
        cells_neighb_set_SL = set()     ## SL: same level
        for i in range(len(c_nodes)):
            cells_nb_i = self.node_conn_cells[c_nodes[i]][c_lev]
            for j in range(len(cells_nb_i)):
                cells_neighb_set_SL.add(cells_nb_i[j])
        
        cells_neighb_set_SL.remove(c_ind)
        
        ##check if there is 2**(N-1) nodes in common
        N = self.N
        dx_lev = self.dx_levels[c_lev]
        nodesPoints = self.nodesPoints
        node_corner = nodesPoints[c_nodes[0]]
        cells_nb_neg = [None]*N
        cells_nb_pos = [None]*N
        cells_nb_next_neg = [None]*N
        cells_nb_next_pos = [None]*N
        for i in range(N):
            cells_nb_next_neg[i] = []
            cells_nb_next_pos[i] = []
        for cells_nb_i in cells_neighb_set_SL:
            nodes_nb = self.cellsHier[c_lev][cells_nb_i][H_CN]
            nodesInds_comm = []
            for n_nb in nodes_nb:
                for n in c_nodes:
                    if n_nb==n:
                        nodesInds_comm.append(n)
                        break
            if len(nodesInds_comm)==2**(N-1):
                ##common face: find orientation
                nodeInd_ref = nodesInds_comm[0]
                Dir = None  ##direction
                for i in range(N):
                    Dir = i
                    for j in range(1, len(nodesInds_comm)):
                        if abs(nodesPoints[nodesInds_comm[j]][i] - nodesPoints[nodeInd_ref][i]) > dx_lev[i]/10:
                            Dir = None
                            break
                    if Dir!=None:
                        break
                if nodesPoints[nodeInd_ref][Dir]-node_corner[Dir]>dx_lev[Dir]/10:
                    cells_nb_pos[Dir] = cells_nb_i
                    c_childs = self.cellsHier[c_lev][cells_nb_i][H_CC]
                    if c_childs!=NEXIST:
                        for c_ch_ind in range(len(c_childs)):
                            c_ch = c_childs[c_ch_ind]
                            c_ch_nodes = self.cellsHier[c_lev+1][c_ch][H_CN]
                            for n_comm_ in nodesInds_comm:
                                if n_comm_ in c_ch_nodes:
                                    if getConnNodes:
                                        n_comm__ch_ind_ = c_ch_nodes.index(n_comm_)
                                        n_comm_____ind_ = c_nodes.index(n_comm_)
                                        cells_nb_next_pos[Dir].append([c_ch, (n_comm_____ind_, n_comm__ch_ind_)])
                                    else:
                                        cells_nb_next_pos[Dir].append(c_ch)
                                    break
                else:
                    cells_nb_neg[Dir] = cells_nb_i
                    c_childs = self.cellsHier[c_lev][cells_nb_i][H_CC]
                    if c_childs!=NEXIST:
                        for c_ch_ind in range(len(c_childs)):
                            c_ch = c_childs[c_ch_ind]
                            c_ch_nodes = self.cellsHier[c_lev+1][c_ch][H_CN]
                            for n_comm_ in nodesInds_comm:
                                if n_comm_ in c_ch_nodes:
                                    if getConnNodes:
                                        n_comm__ch_ind_ = c_ch_nodes.index(n_comm_)
                                        n_comm_____ind_ = c_nodes.index(n_comm_)
                                        cells_nb_next_neg[Dir].append([c_ch, (n_comm_____ind_, n_comm__ch_ind_)])
                                    else:
                                        cells_nb_next_neg[Dir].append(c_ch)
                                    break
        
        return [[cells_nb_neg, cells_nb_pos], [cells_nb_next_neg, cells_nb_next_pos]]
        
        

    def CellGetNeighborNodesSameNexLev(self, c_lev, c_ind):
        """ It returns the nodes of the neighboring cells in the same and next 
        level (near neigbor and next near neighbors for the next level)

        c_lev: cell level
        c_ind : cell index
        """
        cells_neighb_set_SL,  cells_neighb_set_NL = self.CellGetCellsAroundSameNexLev(c_lev, c_ind)
        nodes_neighb_set = set()
        for c in cells_neighb_set_SL:
            c_nodes = self.cellsHier[c_lev][c][H_CN]
            for n in c_nodes:
                nodes_neighb_set.add(n)
        for c in cells_neighb_set_NL:
            assert c_lev<len(self.cellsHier)-1
            c_nodes = self.cellsHier[c_lev+1][c][H_CN]
            for n in c_nodes:
                nodes_neighb_set.add(n)
        return nodes_neighb_set
        

    def RefineSingleCell(self, c_lev, c_ind):
        c_no, c_pa, c_ch = H_CN, H_CP, H_CC
        N = self.N
        assert c_lev<len(self.cellsHier)
        assert c_ind<len(self.cellsHier[c_lev])
        assert self.cellsHier[c_lev][c_ind][c_ch]==NEXIST
        ## get all neigbor points (to check later for duplicates)
        nodes_neighb_set = self.CellGetNeighborNodesSameNexLev(c_lev, c_ind)
        ## --        
        self.cellsHier[c_lev][c_ind][c_ch] = []
        cell_children = self.cellsHier[c_lev][c_ind][c_ch]        
        if(len(self.cellsHier)==c_lev+1):
            self.cellsHier.append([])
        if(len(self.cellsWithNodesOnBorder)==c_lev+1):
            self.cellsWithNodesOnBorder.append({'n':[None]*N, 'p':[None]*N})
            for i in range(N):
                self.cellsWithNodesOnBorder[c_lev+1]['n'][i] = []
                self.cellsWithNodesOnBorder[c_lev+1]['p'][i] = []
        is_border_cell_n = [False]*N
        is_border_cell_p = [False]*N
        for i in range(N):
            ##TODO: check search performance (use a preset True/False list to avoid search)
            if c_ind in self.cellsWithNodesOnBorder[c_lev]['n'][i]:
                is_border_cell_n[i] = True 
            if c_ind in self.cellsWithNodesOnBorder[c_lev]['p'][i]:
                is_border_cell_p[i] = True
        is_border_cell = False
        if (True in is_border_cell_n) or (True in is_border_cell_p):
            is_border_cell = True
        cellsWithNodesOnBorder_levp1_n = self.cellsWithNodesOnBorder[c_lev+1]['n']
        cellsWithNodesOnBorder_levp1_p = self.cellsWithNodesOnBorder[c_lev+1]['p']
            
        lev_p1 = c_lev+1
        cells_lev_p1 = self.cellsHier[lev_p1]
        n_cells_lev_p1 = len(cells_lev_p1)
        cell_nodes = self.cellsHier[c_lev][c_ind][c_no]
        nodesPoints = self.nodesPoints
        n_nodes_last = len(nodesPoints)
        dx__min = np.min(self.dx_levels[c_lev+1])
        dx_thrsh = dx__min/10.0
        node_conn_cells = self.node_conn_cells
        ## subcells
        subcells_nodes = [None]*2**N ## contains subcell nodes for each subcell
        for i in range(2**N):
            subcells_nodes[i] = []
        ## scale the cell (scale all cell dimensions by half)
        dx_scaled = [None]*2**N
        node_corner = nodesPoints[cell_nodes[0]]
        for i in range(2**N):
            dx_scaled[i] = (nodesPoints[cell_nodes[i]] - node_corner)/2.0
        ## shift the scaled cell 2**N times
        di = self._getBinaryCounter()
        for i in range(2**N):
            subcell_nodesPoints_i = [-1]*2**N
            subcell_node_corner = node_corner + di*self.dx_levels[c_lev+1]
            #print(di)
            #print(subcell_node_corner)
            for j in range(2**N):
                subcell_nodesPoints_i[j] = subcell_node_corner + dx_scaled[j]
                ## take care of the duplicated nodes
                # check for duplicates: get the neighbor cells and check the nodes
                # of the neighbor cells in the current and next hierarchical level
                node_exists = -1
                for n_nb in nodes_neighb_set:
                    if np.linalg.norm(subcell_nodesPoints_i[j] - nodesPoints[n_nb])<dx_thrsh:
                        ## node already exists
                        node_exists = n_nb
                        break
                if node_exists>=0:
                    subcells_nodes[i].append(node_exists)
                else:
                    n_new = len(nodesPoints)
                    subcells_nodes[i].append(n_new)
                    # add new nodes to self.nodesPoints
                    nodesPoints.append(subcell_nodesPoints_i[j])
                    nodes_neighb_set.add(n_new)
                    node_conn_cells.append({})
            ## update self.cellsWithNodesOnBorder
            if is_border_cell:
                for j in range(N):
                    if is_border_cell_n[j] and di[j]==0:
                        cellsWithNodesOnBorder_levp1_n[j].append(len(cells_lev_p1))
                    if is_border_cell_p[j] and di[j]==1:
                        cellsWithNodesOnBorder_levp1_p[j].append(len(cells_lev_p1))

            #print(subcell_nodesPoints_i)
            #print(subcells_nodes[i])
            #the order of the next lines is important
            assert self.getBinaryCounterValue(di)==len(cell_children)
            cell_children.append(len(cells_lev_p1)) 
            cells_lev_p1.append([subcells_nodes[i], c_ind, NEXIST, IND_NEXIST])
            if not self._binaryCounterIncrease(di):
                break
        
        
        #update self.node_conn_cells
        for c_ind_ in range(n_cells_lev_p1, len(cells_lev_p1)):
            c_nodes = cells_lev_p1[c_ind_][H_CN]    ## cell nodes
            for i in range(len(c_nodes)):
                node = c_nodes[i]
                if lev_p1 in node_conn_cells[node]:
                    node_conn_cells[node][lev_p1].append(c_ind_)
                else:
                    node_conn_cells[node][lev_p1] = [c_ind_]
        
        return
        
        
    def RefineCells(self, cells_ref):
        n_lev = len(self.cellsHier)
        assert len(cells_ref)==n_lev
        for lev in range(n_lev):
            for c_ind in cells_ref[lev]:
                if self.cellsHier[lev][c_ind][H_CC]==NEXIST:
                    self.RefineSingleCell(lev, c_ind)
                        

    def SetupBoundaryCellConnections(self):
        """
        Sets up self.boundaryCellConnections after refinement
        """
        N = self.N
        nx = self.nx
        nx_pow = self.nx_pow
        
        ## first level
        
        ## update cells connected by periodic boundary conditions (border cells)
        # [neg[Dir]:cells, pos[Dir]:cells] --> same neg and pos indices are connected
        self.boundaryCellConnections = [{'n':[None]*N, 'p':[None]*N}]  
        for i in range(N):
            self.boundaryCellConnections[0]['n'][i] = []
            self.boundaryCellConnections[0]['p'][i] = []
        boundaryCellConnections_0_n = self.boundaryCellConnections[0]['n']
        boundaryCellConnections_0_p = self.boundaryCellConnections[0]['p']
        for i in range(N):
            counter_0 = np.zeros(nx.shape, dtype=int)
            counter_1 = np.zeros(nx.shape, dtype=int)
            counter_1[i] = nx[i]-1
            while True:
                cell_ind_0 = self._uniformGridCellInd(counter_0, nx, nx_pow)
                cell_ind_1 = self._uniformGridCellInd(counter_1, nx, nx_pow)
                boundaryCellConnections_0_n[i].append(cell_ind_0)
                boundaryCellConnections_0_p[i].append(cell_ind_1)
                counter_0_increased = self._increaseCellCounterIndex_Masked(counter_0, nx, [i])
                counter_1_increased = self._increaseCellCounterIndex_Masked(counter_1, nx, [i])
                assert counter_0_increased==counter_1_increased
                if not counter_0_increased:
                    break
                    
        ## higher levels
        n_lev = len(self.cellsHier)
        for lev in range(1, n_lev):
            self.boundaryCellConnections.append({'n':[None]*N, 'p':[None]*N})
            for n in range(N):  ##n: direction
                self.boundaryCellConnections[lev]['n'][n] = []
                self.boundaryCellConnections[lev]['p'][n] = []
                boundaryCellConnections_lev_n = self.boundaryCellConnections[lev]['n'][n]
                boundaryCellConnections_lev_p = self.boundaryCellConnections[lev]['p'][n]
                n_prev = len(self.boundaryCellConnections[lev-1]['n'][n])
                boundaryCellConnections__1_neg = self.boundaryCellConnections[lev-1]['n'][n]
                boundaryCellConnections__1_pos = self.boundaryCellConnections[lev-1]['p'][n]
                for i in range(n_prev):
                    cell_1_n = boundaryCellConnections__1_neg[i]    ## _n: neg
                    cell_1_p = boundaryCellConnections__1_pos[i]    ## _p: pos
                    cells_nn = self.CellGetSubchilrenInDirection(lev-1, cell_1_n, n, 'n')
                    cells_pp = self.CellGetSubchilrenInDirection(lev-1, cell_1_p, n, 'p')
                    if cells_nn!=None and cells_pp!=None:
                        boundaryCellConnections_lev_n.extend(cells_nn)
                        boundaryCellConnections_lev_p.extend(cells_pp)
        return



    def GetInternalAndBoundaryCells(self, cells_ref):
        """ cells_ref[lev][cell_ind]
        it decomposes cells_ref to two lists cells_ref_intern[lev][cell_ind] and 
        cells_ref_bound[lev][cell_ind] 
        """
        N = self.N
        cellsHier = self.cellsHier
        n_lev = len(cellsHier)
        cells_bound_marked = self.MarkBoundaryCells(info='simple')
        cells_ref_intern = [None]*n_lev
        cells_ref_bound = [None]*n_lev
        n_cells_ref_intern = [0]*n_lev
        n_cells_ref_bound = [0]*n_lev
        for lev in range(n_lev):
            n_ref = len(cells_ref[lev])
            for i in range(n_ref):
                if cells_bound_marked[lev][cells_ref[i]]:
                    n_cells_ref_bound[lev] += 1
                else:
                    n_cells_ref_intern[lev] += 1
        
        for lev in range(n_lev):
            cells_ref_intern[lev] = [None]*n_cells_ref_intern[lev]
            cells_ref_bound[lev] = [None]*n_cells_ref_bound[lev]
            n_ref = len(cells_ref[lev])
            n_b, n_c = 0, 0
            for i in range(n_ref):
                c_ind = cells_ref[i]
                if cells_bound_marked[lev][c_ind]:
                    cells_ref_bound[lev][n_b] = c_ind
                    n_b += 1
                else:
                    n_cells_ref_intern[lev][n_i] = c_ind
                    n_c += 1
        
        return [cells_ref_intern, cells_ref_bound]
        

    def MarkAdditionalCellsForRefinement(self, cells_ref):
        """ cells_ref[lev][cell_ind]
        Appends additional cells to refine to cells_ref such that in the refined
        structure there is no neighbors with more than one level difference.
        """        
        cellsHier = self.cellsHier
        n_lev = len(cellsHier)
        n_nodes = len(self.nodesPoints)
        assert len(cells_ref)==n_lev
        cells_marked = self.MarkCells(cells_ref)
        ##TODO: construct and keep track of node_conn_cells_leaf after each refinement
        ## for speedup
        node_conn_cells_leaf = self.NodesGetConnectedLeafCells()
        nodes_redo = []
        #print('node_conn_cells_leaf: ', node_conn_cells_leaf)
        for n in range(n_nodes):
            for lev in sorted(node_conn_cells_leaf[n].keys(), reverse=True):
                node_conn_cells_leaf_n_lev = node_conn_cells_leaf[n][lev]
                for c_ind in node_conn_cells_leaf_n_lev: 
                    if cells_marked[lev][c_ind]:
                        if lev-1 in node_conn_cells_leaf[n]:
                            #print(node_conn_cells_leaf[n][lev-1])
                            node_conn_cells_leaf_n_lev_1 = node_conn_cells_leaf[n][lev-1]
                            for c_ind_1 in node_conn_cells_leaf_n_lev_1:
                                if not cells_marked[lev-1][c_ind_1]:
                                    cells_marked[lev-1][c_ind_1] = True
                                    ## add the nodes of the newly marked cells
                                    ## to a list and repeat the process for them
                                    nodes_redo.extend(cellsHier[lev-1][c_ind_1][H_CN])

        nodes_redo = list(set(nodes_redo))
        while len(nodes_redo)>0:
            nodes_redo_next = []
            for n in nodes_redo:
                for lev in sorted(node_conn_cells_leaf[n].keys(), reverse=True):
                    node_conn_cells_leaf_n_lev = node_conn_cells_leaf[n][lev]
                    for c_ind in node_conn_cells_leaf_n_lev: 
                        if cells_marked[lev][c_ind]:
                            if lev-1 in node_conn_cells_leaf[n]:
                                #print(node_conn_cells_leaf[n][lev-1])
                                node_conn_cells_leaf_n_lev_1 = node_conn_cells_leaf[n][lev-1]
                                for c_ind_1 in node_conn_cells_leaf_n_lev_1:
                                    if not cells_marked[lev-1][c_ind_1]:
                                        cells_marked[lev-1][c_ind_1] = True
                                        nodes_redo_next.extend(cellsHier[lev-1][c_ind_1][H_CN])
            nodes_redo = list(set(nodes_redo_next))
                    
        ##construct cells_ref_new
        cells_ref_new = [None]*n_lev
        n_ref = [0]*n_lev
        for lev in range(n_lev):
            n_ref_lev = 0
            for c_ind in range(len(cells_marked[lev])):
                if cells_marked[lev][c_ind]:
                    n_ref_lev += 1
            n_ref[lev] = n_ref_lev
                    
        for lev in range(n_lev):
            cells_ref_new[lev] = [None]*n_ref[lev]
            n_ref_lev = 0
            for c_ind in range(len(cells_marked[lev])):
                if cells_marked[lev][c_ind]:
                    cells_ref_new[lev][n_ref_lev] = c_ind
                    n_ref_lev += 1
        
        return cells_ref_new


    def MarkAdditionalCellsForRefinement_PeriodicBoundary(self, cells_ref):
        """ cells_ref[lev][cell_ind]
        Appends additional cells to refine to cells_ref such that in the refined
        structure there is no neighbors with more than one level difference.
        In addition if two cells are connected through periodic boundary 
        (corressponding cells at the beginning and end boundaries) their refinement
        levels are identical
        i.e. Appends additional cells to refine to cells_ref such that in the refined
        structure, for each boundary cell (to refine) the coressponding boundary 
        cell (through periodic boundary) will be refined as well.
        """        
        cellsHier = self.cellsHier
        n_lev = len(cellsHier)
        N = self.N
        n_nodes = len(self.nodesPoints)
        assert len(cells_ref)==n_lev
        cells_marked = self.MarkCells(cells_ref)

        cells_bound_marked_simple = self.MarkBoundaryCells(info='simple')
        cells_bound_marked_detail = self.MarkBoundaryCells(info='detailed-connected')

        node_conn_cells_leaf = self.NodesGetConnectedLeafCells()
        nodes_border_updated = []
        nodes_redo = []
        for n in range(n_nodes):
            for lev in sorted(node_conn_cells_leaf[n].keys(), reverse=True):
                node_conn_cells_leaf_n_lev = node_conn_cells_leaf[n][lev]
                for c_ind in node_conn_cells_leaf_n_lev: 
                    if cells_marked[lev][c_ind]:
                        if cells_bound_marked_simple[lev][c_ind]:
                            for n_dir in range(N):
                                c_ind_opp = cells_bound_marked_detail[lev]['n'][n_dir][c_ind]
                                if c_ind_opp>=0:
                                    if not cells_marked[lev][c_ind_opp]:
                                        cells_marked[lev][c_ind_opp] = True
                                        nodes_border_updated.extend(cellsHier[lev][c_ind_opp][H_CN])
                                c_ind_opp = cells_bound_marked_detail[lev]['p'][n_dir][c_ind]
                                if c_ind_opp>=0:
                                    if not cells_marked[lev][c_ind_opp]:
                                        cells_marked[lev][c_ind_opp] = True
                                        nodes_border_updated.extend(cellsHier[lev][c_ind_opp][H_CN])
                                    
                        if lev-1 in node_conn_cells_leaf[n]:
                            node_conn_cells_leaf_n_lev_1 = node_conn_cells_leaf[n][lev-1]
                            for c_ind_1 in node_conn_cells_leaf_n_lev_1:
                                if not cells_marked[lev-1][c_ind_1]:
                                    cells_marked[lev-1][c_ind_1] = True
                                    nodes_redo.extend(cellsHier[lev-1][c_ind_1][H_CN])
                                    if cells_bound_marked_simple[lev-1][c_ind_1]:
                                        for n_dir in range(N):
                                            c_ind_opp = cells_bound_marked_detail[lev-1]['n'][n_dir][c_ind_1]
                                            if c_ind_opp>=0:
                                                if not cells_marked[lev-1][c_ind_opp]:
                                                    cells_marked[lev-1][c_ind_opp] = True
                                                    nodes_border_updated.extend(cellsHier[lev-1][c_ind_opp][H_CN])
                                            c_ind_opp = cells_bound_marked_detail[lev-1]['p'][n_dir][c_ind_1]
                                            if c_ind_opp>=0:
                                                if not cells_marked[lev-1][c_ind_opp]:
                                                    cells_marked[lev-1][c_ind_opp] = True
                                                    nodes_border_updated.extend(cellsHier[lev-1][c_ind_opp][H_CN])

        nodes_redo.extend(nodes_border_updated)
        nodes_redo = list(set(nodes_redo))
        if self.verbose:
            print('nodes_redo:', nodes_redo)
        while len(nodes_redo)>0:
            nodes_border_updated_nx = []
            nodes_redo_next = []
            for n in nodes_redo:
                for lev in sorted(node_conn_cells_leaf[n].keys(), reverse=True):
                    node_conn_cells_leaf_n_lev = node_conn_cells_leaf[n][lev]
                    for c_ind in node_conn_cells_leaf_n_lev: 
                        if cells_marked[lev][c_ind]:
                            if cells_bound_marked_simple[lev][c_ind]:
                                for n_dir in range(N):
                                    c_ind_opp = cells_bound_marked_detail[lev]['n'][n_dir][c_ind]
                                    if c_ind_opp>=0:
                                        if not cells_marked[lev][c_ind_opp]:
                                            cells_marked[lev][c_ind_opp] = True
                                            nodes_border_updated_nx.extend(cellsHier[lev][c_ind_opp][H_CN])
                                    c_ind_opp = cells_bound_marked_detail[lev]['p'][n_dir][c_ind]
                                    if c_ind_opp>=0:
                                        if not cells_marked[lev][c_ind_opp]:
                                            cells_marked[lev][c_ind_opp] = True
                                            nodes_border_updated_nx.extend(cellsHier[lev][c_ind_opp][H_CN])
                                        
                            if lev-1 in node_conn_cells_leaf[n]:
                                node_conn_cells_leaf_n_lev_1 = node_conn_cells_leaf[n][lev-1]
                                for c_ind_1 in node_conn_cells_leaf_n_lev_1:
                                    if not cells_marked[lev-1][c_ind_1]:
                                        cells_marked[lev-1][c_ind_1] = True
                                        nodes_redo_next.extend(cellsHier[lev-1][c_ind_1][H_CN])
                                        if cells_bound_marked_simple[lev-1][c_ind_1]:
                                            for n_dir in range(N):
                                                c_ind_opp = cells_bound_marked_detail[lev-1]['n'][n_dir][c_ind_1]
                                                if c_ind_opp>=0:
                                                    if not cells_marked[lev-1][c_ind_opp]:
                                                        cells_marked[lev-1][c_ind_opp] = True
                                                        nodes_border_updated_nx.extend(cellsHier[lev-1][c_ind_opp][H_CN])
                                                c_ind_opp = cells_bound_marked_detail[lev-1]['p'][n_dir][c_ind_1]
                                                if c_ind_opp>=0:
                                                    if not cells_marked[lev-1][c_ind_opp]:
                                                        cells_marked[lev-1][c_ind_opp] = True
                                                        nodes_border_updated_nx.extend(cellsHier[lev-1][c_ind_opp][H_CN])
            
            #nodes_border_updated = list(set(nodes_border_updated_nx))
            nodes_redo_next.extend(nodes_border_updated_nx)
            nodes_redo = list(set(nodes_redo_next))
            if self.verbose:
                print('nodes_border_updated:', len(nodes_border_updated))
                print('nodes_redo:', len(nodes_redo))
            
        ##construct cells_ref_new
        cells_ref_new = [None]*n_lev
        n_ref = [0]*n_lev
        for lev in range(n_lev):
            n_ref_lev = 0
            for c_ind in range(len(cells_marked[lev])):
                if cells_marked[lev][c_ind]:
                    n_ref_lev += 1
            n_ref[lev] = n_ref_lev
                    
        for lev in range(n_lev):
            cells_ref_new[lev] = [None]*n_ref[lev]
            n_ref_lev = 0
            for c_ind in range(len(cells_marked[lev])):
                if cells_marked[lev][c_ind]:
                    cells_ref_new[lev][n_ref_lev] = c_ind
                    n_ref_lev += 1
        
        return cells_ref_new


    def MarkAdditionalCellsForRefinement_BCreflective(self, cells_ref):
        ##TODO
        """ mark cells that are connected through reflection with respect to 'p' 
        side boundary conditions, for example in 2d if x'p' and y'p' boundary conditions
        are specified, a cell on the boundary x'p' is first reflected with respect 
        to the y axis, and then with respect to x axis, to get its mirrored cell.
        """


    def VerifyCellsToRefineContinuity(self, cells_ref):
        """ The difference in level between two leaf neigbor cells is verfied
        to be (will be after refinement) at most one.
        """
        cellsHier = self.cellsHier
        n_lev = len(cellsHier)
        n_nodes = len(self.nodesPoints)
        assert len(cells_ref)==n_lev
        cells_marked = self.MarkCells(cells_ref)
        node_conn_cells_leaf = self.NodesGetConnectedLeafCells()
        for n in range(n_nodes):
            for lev in sorted(node_conn_cells_leaf[n].keys(), reverse=True):
                node_conn_cells_leaf_n_lev = node_conn_cells_leaf[n][lev]
                for c_ind in node_conn_cells_leaf_n_lev: 
                    if cells_marked[lev][c_ind]:
                        if lev-1 in node_conn_cells_leaf[n]:
                            node_conn_cells_leaf_n_lev_1 = node_conn_cells_leaf[n][lev-1]
                            for c_ind_1 in node_conn_cells_leaf_n_lev_1:
                                if not cells_marked[lev-1][c_ind_1]:
                                    return False
        return True


    def VerifyCellsContinuity(self):
        """ The difference in level between two leaf neigbor cells is verfied
        to be at most one. (after refinements have been done)
        """
        cellsHier = self.cellsHier
        n_lev = len(cellsHier)
        n_nodes = len(self.nodesPoints)
        node_conn_cells_leaf = self.NodesGetConnectedLeafCells()
        for n in range(n_nodes):
            for lev in sorted(node_conn_cells_leaf[n].keys(), reverse=True):
                node_conn_cells_leaf_n_lev = node_conn_cells_leaf[n][lev]
                for c_ind in node_conn_cells_leaf_n_lev:
                    for lev_l2 in range(lev-1): 
                        if lev_l2 in node_conn_cells_leaf[n]:
                            if len(node_conn_cells_leaf[n][lev_l2])>0:
                                #print('lev_l2:', lev_l2, 'n', n, \
                                #'node_conn_cells_leaf[n][lev_l2]:', node_conn_cells_leaf[n][lev_l2])
                                return False
        return True

        
    def VerifyInternalAndBorderCellsContinuity(self):
        """ The difference in level between two leaf neigbor cells is verfied
        to be at most one. (after refinements have been done)
        """
        N = self.N
        cellsHier = self.cellsHier
        n_lev = len(cellsHier)
        n_nodes = len(self.nodesPoints)
        node_conn_cells_leaf = self.NodesGetConnectedLeafCells()

        cells_bound_marked_simple = self.MarkBoundaryCells(info='simple')
        cells_bound_marked_detail = self.MarkBoundaryCells(info='detailed-connected')

        for n in range(n_nodes):
            for lev in sorted(node_conn_cells_leaf[n].keys(), reverse=True):
                node_conn_cells_leaf_n_lev = node_conn_cells_leaf[n][lev]
                for c_ind in node_conn_cells_leaf_n_lev:
                    ##internal connections
                    for lev_l2 in range(lev-1): 
                        if lev_l2 in node_conn_cells_leaf[n]:
                            if len(node_conn_cells_leaf[n][lev_l2])>0:
                                #print('lev_l2:', lev_l2, 'n', n, \
                                #'node_conn_cells_leaf[n][lev_l2]:', node_conn_cells_leaf[n][lev_l2])
                                return False
                    ##boundary connections
                    if cells_bound_marked_simple[lev][c_ind]:
                        for n_dir in range(N):
                            ##TODO: verify at the opposite side cell always exists
                            c_ind_opp = cells_bound_marked_detail[lev]['n'][n_dir][c_ind]
                            if c_ind_opp>=0:
                                if cellsHier[lev][c_ind_opp][H_CC]!=NEXIST:
                                    return False
                            c_ind_opp = cells_bound_marked_detail[lev]['p'][n_dir][c_ind]
                            if c_ind_opp>=0:
                                if cellsHier[lev][c_ind_opp][H_CC]!=NEXIST:
                                    return False
        return True
        
        
    def MarkCells(self, cells_ref):
        """ returns an array n_lev*n_cell_lev where each cell in cells_ref is 
        marked as True
        """
        cellsHier = self.cellsHier
        n_lev = len(cellsHier)
        cellsMarked = [None]*n_lev
        for lev in range(n_lev):
            n_cells_lev = len(cellsHier[lev])
            cellsMarked[lev] = [False]*n_cells_lev
            for c_ind in cells_ref[lev]:
                cellsMarked[lev][c_ind] = True
        return cellsMarked
        

    def MarkBoundaryCells(self, info='simple'):
        """ returns an array n_lev*n_cell_lev where each boundary cell is marked
        'simple': [lev][cell] = True/False
        'detailed': [lev]['n'/'p'][n][cell] = True/False
        """
        if info=='simple':
            N = self.N
            cellsHier = self.cellsHier
            n_lev = len(cellsHier)
            cellsMarked = [None]*n_lev
            for lev in range(n_lev):
                n_cells_lev = len(cellsHier[lev])
                cellsMarked[lev] = [False]*n_cells_lev
                for n in range(N):
                    for c_ind in self.cellsWithNodesOnBorder[lev]['n'][n]:
                        cellsMarked[lev][c_ind] = True
                    for c_ind in self.cellsWithNodesOnBorder[lev]['p'][n]:
                        cellsMarked[lev][c_ind] = True
            return cellsMarked
        elif info=='detailed':
            N = self.N
            cellsHier = self.cellsHier
            n_lev = len(cellsHier)
            cellsMarked = [None]*n_lev
            for lev in range(n_lev):
                n_cells_lev = len(cellsHier[lev])
                cellsMarked[lev] = {'n':[None]*N, 'p':[None]*N}
                for n in range(N):
                    cellsMarked[lev]['n'][n] = [False]*n_cells_lev
                    cellsMarked_lev_n_n = cellsMarked[lev]['n'][n]
                    for c_ind in self.cellsWithNodesOnBorder[lev]['n'][n]:
                        cellsMarked_lev_n_n[c_ind] = True
                    cellsMarked[lev]['p'][n] = [False]*n_cells_lev
                    cellsMarked_lev_p_n = cellsMarked[lev]['p'][n]
                    for c_ind in self.cellsWithNodesOnBorder[lev]['p'][n]:
                        cellsMarked_lev_p_n[c_ind] = True
            return cellsMarked
        elif info=='detailed-connected':
            N = self.N
            cellsHier = self.cellsHier
            n_lev = len(cellsHier)
            cellsMarked = [None]*n_lev
            for lev in range(n_lev):
                n_cells_lev = len(cellsHier[lev])
                cellsMarked[lev] = {'n':[None]*N, 'p':[None]*N}
                for n in range(N):
                    cellsMarked[lev]['n'][n] = [-1]*n_cells_lev
                    cellsMarked_lev_n_n = cellsMarked[lev]['n'][n]
                    for _ind, c_ind in enumerate(self.boundaryCellConnections[lev]['n'][n]):
                        cellsMarked_lev_n_n[c_ind] = self.boundaryCellConnections[lev]['p'][n][_ind]
                        #print(cellsMarked_lev_n_n[c_ind])
                    cellsMarked[lev]['p'][n] = [-1]*n_cells_lev
                    cellsMarked_lev_p_n = cellsMarked[lev]['p'][n]
                    for _ind, c_ind in enumerate(self.boundaryCellConnections[lev]['p'][n]):
                        cellsMarked_lev_p_n[c_ind] = self.boundaryCellConnections[lev]['n'][n][_ind]
            return cellsMarked
        

    def GetNeigborCellTypesInTermsOfConnectedNodes(self, direction='neg'):
        """ categorize connected cells (with a face in common) in terms of their
        level, connected nodes. 
        
        cells_nb_dic = {'n':{'sl':[None]*n_lev, 'nl':[None]*n_lev, 'pl':[None]*n_lev}}
        
        For neigbor cells in the same level:
        cells_nb_dic['n']['sl'][lev] = [.., (cell_ind, cell_nb_ind), ..]  #nb:neighbor 
        
        For neigbor cells in the next/previous level:
        cells_nb_dic['n']['nl'/'pl'][lev][(c_ind_node, c_nb_node)] = [.., (cell_ind, cell_nb_ind), ..]
        (c_ind_node, c_nb_node) : common nodes (0..7 in 3 dimensions)
                
        direction='neg': get the neighbor cells in the negative direction only
        (the neighbor cell is in lefr/down/...)
        """
        N = self.N
        cellsHier = self.cellsHier
        n_lev = len(cellsHier)
        if direction=='neg':
            cells_nb_dic = {'n':{'sl':[None]*n_lev, 'nl':[None]*n_lev, 'pl':[None]*n_lev}}
            ## 'n' : negative
            ## 'sl': same level 
            ## 'nl': next level (neighbor cell is in next level)
            ## 'pl': previous level (parent level)

            for lev in range(n_lev):
                cells_nb_dic['n']['sl'][lev] = [None]*N
                cells_nb_dic['n']['nl'][lev] = [None]*N
                cells_nb_dic['n']['pl'][lev] = [None]*N
                for i in range(N):
                    cells_nb_dic['n']['sl'][lev][i] = []
                    cells_nb_dic['n']['nl'][lev][i] = {}
                    cells_nb_dic['n']['pl'][lev][i] = {}
            
            for lev in range(n_lev):
                cells_nb_dic_neg_sl_lev = cells_nb_dic['n']['sl'][lev]
                cells_nb_dic_neg_nl_lev = cells_nb_dic['n']['nl'][lev]
                cells_nb_dic_neg_pl_levp1 = None
                if lev<n_lev-1:
                    cells_nb_dic_neg_pl_levp1 = cells_nb_dic['n']['pl'][lev+1]

                cellsHier_lev = cellsHier[lev]
                cellsHier_levp1 = None
                if lev<n_lev-1:
                    cellsHier_levp1 = cellsHier[lev+1]
                for c_ind in range(len(cellsHier_lev)):
                    if cellsHier_lev[c_ind][H_CC]==NEXIST:
                        assert cellsHier_lev[c_ind][H_CI]!=IND_NEXIST
                        cellsConnSaNx = self.CellGetFaceConnectedCellsSameNextLevel(lev, 
                                c_ind, getConnNodes=True)
                        cellsConnSL, cellsConnNL = cellsConnSaNx
                        cells_nb_SL_neg, cells_nb_SL_pos = cellsConnSL
                        cells_nb_NL_neg, cells_nb_NL_pos = cellsConnNL
                        
                        for i in range(N):
                            if cells_nb_SL_neg[i]!=None:
                                assert cells_nb_SL_neg[i]>=0
                                if cellsHier_lev[cells_nb_SL_neg[i]][H_CC]==NEXIST:
                                    cells_nb_dic_neg_sl_lev[i].append((c_ind, cells_nb_SL_neg[i]))
                                else:
                                    cells_nb_NL_neg_i = cells_nb_NL_neg[i]
                                    for j in range(len(cells_nb_NL_neg_i)):
                                        c_nb_ = cells_nb_NL_neg_i[j][0]
                                        ##TODO: check if the following condition is automatically checked before...
                                        if cellsHier_levp1[c_nb_][H_CC]==NEXIST:
                                            n_conn_ = cells_nb_NL_neg_i[j][1]
                                            if n_conn_ in cells_nb_dic_neg_nl_lev[i]:
                                                cells_nb_dic_neg_nl_lev[i][n_conn_].append((c_ind, c_nb_))
                                            else:
                                                cells_nb_dic_neg_nl_lev[i][n_conn_] = [(c_ind, c_nb_)]
                            if cells_nb_SL_pos[i]!=None:
                                assert cells_nb_SL_pos[i]>=0
                                if cellsHier_lev[cells_nb_SL_pos[i]][H_CC]==NEXIST:
                                    pass
                                else:
                                    cells_nb_NL_pos_i = cells_nb_NL_pos[i]
                                    for j in range(len(cells_nb_NL_pos_i)):
                                        c_nb_ = cells_nb_NL_pos_i[j][0]
                                        if cellsHier_levp1[c_nb_][H_CC]==NEXIST:
                                            n_conn_ = cells_nb_NL_pos_i[j][1]
                                            n_conn_ = (n_conn_[1], n_conn_[0]) ## (c_nb_ node, c_ind_ node)
                                            if n_conn_ in cells_nb_dic_neg_pl_levp1[i]:
                                                cells_nb_dic_neg_pl_levp1[i][n_conn_].append((c_nb_, c_ind))
                                            else:
                                                cells_nb_dic_neg_pl_levp1[i][n_conn_] = [(c_nb_, c_ind)]
            self.cells_nb_dic = cells_nb_dic
            return cells_nb_dic
        else:
            raise NotImplementedError()
        
        

    ##-----------
    ## setting up the differential equations and the basis functions    
        
    def AttachDiffEquations(self, eqs_list, vars_list, pars_list, indepVars_list, pars_extern=None):
        """ eqs_list: list of equations
            vars_list: list of unknowns
            pars_list: list of known parameters such as variable coefficients
            ind_vars_list: list of independant variables 
        """
        self.eqs_list = eqs_list
        self.vars_list = vars_list
        self.pars_list = pars_list
        self.indepVars_list = indepVars_list #independant variables
        assert len(indepVars_list)==self.N
        if pars_extern!=None:
            self.SetParsExtern(pars_extern)
        else:
            self.pars_extern=None
        if self.verbose>0:
            print('self.pars_extern:', self.pars_extern)
            
        self.VarIndex = {}
        for i in range(len(vars_list)):
            self.VarIndex[vars_list[i]] = i
        if self.verbose>0:
            print('self.VarIndex:', self.VarIndex)

        self.IndepVarIndex = {}
        for i in range(len(indepVars_list)):
            self.IndepVarIndex[indepVars_list[i]] = i
        if self.verbose>0:
            print('self.IndepVarIndex:', self.IndepVarIndex)
        return
        

    def SetParsValues(self, pars_values, pars_types=None, pars_periodicity=None):
        """ pars_types[i]==general --->  pars_values[i]:symbolic function
            pars_types[i]==seperable --->  pars_values[i]:[f_x, f_y, f_z]
            pars_types[i]==partialseperable --->  pars_values[i]:[(t,):f_t, (x,y):f_xy]
                    or pars_values[i]:[(t,):f_t, (x,):f_x, (y,z):f_yz]
            pars_periodicity = [t:(2*pi, [t0,t1]), x:None, y:None]
            --> between t0 and t1 the function can pe considered periodic in t with 
            period 2*pi
        """
        assert len(pars_values)==len(self.pars_list)
        self.pars_values = pars_values 
        if self.verbose>0:
            print('self.pars_list:', self.pars_list, 'self.pars_values:', self.pars_values, sep='\n')
        if pars_types==None:
            self.pars_types = [ParType.general]*len(self.pars_list)
        else:
            assert len(pars_types)==len(self.pars_list)
            self.pars_types = pars_types
        if pars_periodicity==None:
            self.pars_periodicity = [None]*len(self.indepVars_list)
        else:
            assert len(pars_periodicity)==len(self.indepVars_list)
            self.pars_periodicity = pars_periodicity
            
                  
    def SetParsExtern(self, pars_extern):
        """ pars_extern:{par_0:rgnd_0, par_1:rgnd_1}
        par_0,.. external parameters
        rgnd_0... another RGND object having par_0 as a variable
        """
        for par in pars_extern:
            assert par not in self.pars_list
        self.pars_extern = pars_extern
        if self.verbose>0:
            print('self.pars_extern:', self.pars_extern)
        return
                                             
                                     
    def DisintegrateEquations(self):
        EQs_parts = [None]*len(self.eqs_list)
        rhs = [None]*len(self.eqs_list)
        for i in range(len(self.eqs_list)):
            EQs_parts[i] = [None]*len(self.vars_list)
            for j in range(len(self.vars_list)):
                EQs_parts[i][j] = []
        for i_eq, eq in enumerate(self.eqs_list):
            eq_tree = symExpr_generate_tree(eq)
            for i_v, var in enumerate(self.vars_list):
                while eq_tree[1].has(var):
                    coeff_sym, der_ord, der_vars = sym_getcoeff_setzero(eq_tree, var)
                    EQs_parts[i_eq][i_v].append([coeff_sym, der_ord, der_vars])
            rhs[i_eq] = -eq_tree[1]
            
        self.EQs_parts = EQs_parts
        self.rhs = rhs
        
        self.DisintegrateSelfRhs()
        return [EQs_parts, rhs]
        

    def DisintegrateBoundaryConds(self):
        BC_parts = [None]*len(self.BCs)
        BC_rhs = [None]*len(self.BCs)
        for i in range(len(self.BCs)):
            BC_parts[i] = {}
        for bc_ind, b_cond in enumerate(self.BCs):
            n_dir, expr, face = b_cond['dir'], b_cond['expr'], b_cond['face']
            #BC_parts[bc_ind][n_dir][face] = expr
            eq_tree = symExpr_generate_tree(expr)
            for i_v, var in enumerate(self.vars_list):
                while eq_tree[1].has(var):
                    coeff_sym, der_ord, der_vars = sym_getcoeff_setzero(eq_tree, var)
                    if var not in BC_parts[bc_ind]:
                        BC_parts[bc_ind][i_v] = []    
                    BC_parts[bc_ind][i_v].append([coeff_sym, der_ord, der_vars])
            BC_rhs[bc_ind] = -eq_tree[1]
        self.BC_parts = BC_parts
        self.BC_rhs = BC_rhs
        return [BC_parts, BC_rhs]
        

    def DisintegrateSelfRhs(self):
        """ for each element in self.rhs disintegrate it into coeff, der, par 
        """
        #determine what parameters are involved in each self.rhs
        #internal pars
        n_rhs = len(self.rhs)   ## n_rhs=n_eq?
        rhs_pars_inside = [None]*n_rhs
        for i in range(n_rhs):
            rhs_pars_inside[i] = []
        for i_rhs, rhs in enumerate(self.rhs):
            for par in self.pars_list:
                if rhs.has(par):
                    rhs_pars_inside[i_rhs].append(par)
        #external pars
        rhs_pars_ex_inside = [None]*n_rhs
        for i in range(n_rhs):
            rhs_pars_ex_inside[i] = []
        if self.pars_extern!=None:
            for i_rhs, rhs in enumerate(self.rhs):
                for par in self.pars_extern:
                    if rhs.has(par):
                        rhs_pars_ex_inside[i_rhs].append(par)
        for i in range(n_rhs):
            if len(rhs_pars_ex_inside[i])>1:
                raise NotImplementedError()
        #disintegrate
        rhs_parts = [None]*n_rhs
        for i in range(n_rhs):
            rhs_parts[i] = []
        for i_rhs, rhs in enumerate(self.rhs):
            eq_tree = symExpr_generate_tree(rhs)
            for par in rhs_pars_ex_inside[i_rhs]:
                while eq_tree[1].has(par):
                    coeff_sym, der_ord, der_vars = sym_getcoeff_setzero(eq_tree, par)
                    rhs_parts[i_rhs].append([par, CoeffType.par_ext, coeff_sym, der_ord, der_vars])
            has_par_in = False
            for par in rhs_pars_inside[i_rhs]:
                if eq_tree[1].has(par):
                    rhs_parts[i_rhs].append([eq_tree[1], CoeffType.par_int, None, None, None])
                    has_par_in = True
                    break
            if not has_par_in:
                has_var_indep = False
                for v_indep in self.indepVars_list:
                    if eq_tree[1].has(v_indep):
                        rhs_parts[i_rhs].append([eq_tree[1], CoeffType.nonconst_anal, None, None, None])
                        has_var_indep = True
                        break
                if not has_var_indep:
                    if eq_tree[1]!=0:
                        rhs_parts[i_rhs].append([eq_tree[1], CoeffType.const, None, None, None])

        self.rhs_parts = rhs_parts
        self.rhs_pars_ex_inside = rhs_pars_ex_inside
        self.rhs_pars_inside = rhs_pars_inside
        
        if self.verbose>0:
            print('self.rhs:', self.rhs, 'self.rhs_parts:', self.rhs_parts, \
                'self.rhs_pars_ex_inside:', self.rhs_pars_ex_inside,\
                'self.rhs_pars_inside:', self.rhs_pars_inside, '-'*30, sep='\n')
        return


    def ResetBoundaryConditions(self):
        self.BCs = None
             
    def DefineBoundaryConds(self, BCs):
        """ BCs: a list [{'dir':dir, 'face':'n'/'p', 'bc':expr}, ...]
        dir:0/1/2..
        face:'n'/'p'   'n':negative face  'p':positive face
        expr=0 boundary condition  
        """
        self.BCs = BCs
        return
        
    def AddBoundaryCondition(self, BC):
        """ BC: {'dir':dir, 'face':'n'/'p', 'expr':expr}
        dir:0/1/2..
        face:'n'/'p'   'n':negative face  'p':positive face
        expr=0 boundary condition  
        """
        if self.BCs==None:
            self.BCs = []
        self.BCs.append(BC)
        
    
    def DefineContinuityCondsAcrossCellBorders(self, CCs):
        """ CCs:{var:{'dir':dir, 'der':der, 'face':'n'/'p'} ...}
            var:direction:derivative:face
            face='n'/'p' negative or positive face at the given direction
            apply the continituity of the n-th derivative of the variable at the
            beginning/end face of each cell in the given direction(x/y/z)
        """
        self.CCs = CCs
        return
    
     
    def SetMaxPolyOrderToKeepForEachEq(self, orders_max):
        """ orders_max = [(n0, n1, n2..), (n0, n1, n2..) ...]
            orders_max[i] describes the i-th equation
            (n0, n1, n2..) keep up to power n0 in the x direction, n1 in the y 
            direction ...
        """
        self.orders_max = orders_max
        return
    
    def SetMaxDerivativeOrderForEachVar(self, der_orders_max):
        self.der_orders_max = der_orders_max


    def SetPolynomialBasisFuncs(self, poly_ords, indexing='children'):
        """ indexing='children' : only index children
        """
        assert self.vars_list!=None
        n_vars = len(self.vars_list)
        assert len(poly_ords)==n_vars
        for v in range(n_vars):
            assert len(poly_ords[v])==self.N            
        self.poly_ords = poly_ords
        arr_shape = [tuple([poly_ords[v][i]+1 for i in range(self.N)]) for v in range(n_vars)]
        
        cellsHier = self.cellsHier
        n_lev = len(self.cellsHier)
        self.polyHier = [None]*n_vars
        self.indexInfo = [None]*n_lev
        for i in range(n_lev):
            self.indexInfo[i] = {}
        n_vars = len(self.vars_list)
        if indexing=='children':
            self.indexing = indexing
            for lev in range(n_lev):
                n_cell_lev = len(cellsHier[lev])
                cellsHier_lev = cellsHier[lev]
                ## set indices
                ind = 0
                for i in range(n_cell_lev):
                    if cellsHier_lev[i][H_CC]==NEXIST:
                        cellsHier_lev[i][H_CI] = ind
                        ind += 1
                self.indexInfo[lev]['N'] = ind      ## 'N': number of indexed cells
                ## set up self.polyHier
            for v in range(n_vars):
                self.polyHier[v] = [None]*n_lev
                for lev in range(n_lev):
                    n_cell_lev = len(cellsHier[lev])
                    self.polyHier[v][lev] = [None]*n_cell_lev
                    polyHier_v_lev = self.polyHier[v][lev]
                    cellsHier_lev = cellsHier[lev]
                    for i in range(n_cell_lev):
                        if cellsHier_lev[i][H_CC]==NEXIST:
                            polyHier_v_lev[i] = None
        else:
            raise NotImplementedError()
        if self.verbose>0:
            print('self.indexInfo : ', self.indexInfo)       
                             

    def GetNumOfIndexedCells(self):
        n_lev = len(self.cellsHier)
        N_indexed = 0
        for lev in range(n_lev):
            N_indexed += self.indexInfo[lev]['N']
        return N_indexed

    def GetAccumukativeNumOfIndexedCells(self):
        n_lev = len(self.cellsHier)
        N_indexed = 0
        N_indexed_acc = [0]
        for lev in range(n_lev):
            N_indexed += self.indexInfo[lev]['N']
            N_indexed_acc.append(N_indexed)
        
        self.N_indexed_acc = N_indexed_acc
        return N_indexed_acc


    def GetTaylorCoeff(self, f, r_0, h, ijk_tup):
        """ r_0: corner
            h : scaling factor
            f : symbolic function
            ijk_tup: taylor index
        """
        N = self.N
        xyz_list = self.indepVars_list
        arg_der = []
        for i in range(N):
            arg_der.append(xyz_list[i])
            arg_der.append(int(ijk_tup[i]))
        f_ijk = Derivative(f, *tuple(arg_der)).doit()
        f_ijk_lambda = lambdify(tuple(xyz_list), f_ijk)
        fact_ijk = 1.0
        for i in range(N):
            fact_ijk *= math.factorial(ijk_tup[i])
        tailor_ijk = f_ijk_lambda(*tuple(r_0.tolist()))/fact_ijk
        h_sc = h**ijk_tup
        scale = 1.0
        for s in h_sc:
            scale *= s
        tailor_ijk_scaled = tailor_ijk*scale
        return tailor_ijk_scaled
        
        
    def GetTaylorCoeffAll(self, f, r_0, h):
        N = self.N
        mask = self.inds_keep
        shape = mask.shape
        
        taylor = np.zeros(shape)
        #inds = self.inds_0_allp        
        counter = np.zeros(N, dtype=int)
        while True:
            cnt_ind = tuple(counter.tolist())
            if mask[cnt_ind]==1:
                taylor[cnt_ind] = self.GetTaylorCoeff(f, r_0, h, cnt_ind)
            if not self._increaseCellCounterIndex(counter, shape):
                break
        return taylor

    def GetTaylorCoeffAllPresetDers(self, f, r_0, h, f_ders_lambda, factorials, scale_factor):
        N = self.N
        mask = self.inds_keep
        shape = mask.shape
        
        taylor = np.zeros(shape)
        #inds = self.inds_0_allp        
        
        counter = np.zeros(N, dtype=int)
        while True:
            cnt_ind = tuple(counter.tolist())
            if mask[cnt_ind]==1:
                tailor_ijk = f_ders_lambda[cnt_ind](*tuple(r_0.tolist()))/factorials[cnt_ind]
                tailor_ijk_scaled = tailor_ijk*scale_factor[cnt_ind]
                taylor[cnt_ind] = tailor_ijk_scaled
            if not self._increaseCellCounterIndex(counter, shape):
                break
        return taylor

        
    ##TODO: get taylor coeffs with respect to the center of cell and then change
    ## the the reference point to cell corner
    def SetTaylorCoeffAlllevels(self, f):
        N = self.N
        cellsHier = self.cellsHier
        n_lev = len(cellsHier)
        xyz_list = self.indepVars_list
        mask = self.inds_keep
        shape = mask.shape
        factorials = np.zeros(shape)
        f_ders_lambda = np.zeros(shape, dtype=object)
        r_corner = self.x0
        ##set derivatives and factoriels
        counter = np.zeros(N, dtype=int)
        #print('f:', f)
        while True:
            cnt_ind = tuple(counter.tolist())
            if mask[cnt_ind]==1:
                fact_ijk = 1.0
                for i in range(N):
                    fact_ijk *= math.factorial(cnt_ind[i])
                factorials[cnt_ind] = fact_ijk
                
                arg_der = []
                for i in range(N):
                    arg_der.append(xyz_list[i])
                    arg_der.append(int(cnt_ind[i]))
                f_ijk = (Derivative(f, *tuple(arg_der))).doit()
                f_ijk_lambda = lambdify(tuple(xyz_list), f_ijk)
                f_ders_lambda[cnt_ind] = f_ijk_lambda
                #print('cnt_ind:', cnt_ind, 'f_ijk:', f_ijk)

            if not self._increaseCellCounterIndex(counter, shape):
                break
        ##set scale factor
        scale_factor = [None]*n_lev
        for lev in range(n_lev):
            scale_factor[lev] = np.zeros(shape)
            h = self.dx_levels[lev]            
            counter = np.zeros(N, dtype=int)
            while True:
                cnt_ind = tuple(counter.tolist())
                if mask[cnt_ind]==1:
                    h_sc = h**counter
                    scale = 1.0
                    for s in h_sc:
                        scale *= s
                    
                    scale_factor[lev][cnt_ind] = scale

                if not self._increaseCellCounterIndex(counter, shape):
                    break
        
        #print('scale_factor:', scale_factor)
        ##set taylor
        TayCo = [None]*n_lev
        for lev in range(n_lev):
            cellsHier_lev = cellsHier[lev]
            n_cell_lev = len(cellsHier_lev)
            TayCo[lev] = [None]*n_cell_lev
            dx_lev = self.dx_levels[lev]
            for c_ind in range(n_cell_lev):
                if cellsHier_lev[c_ind][H_CC]==NEXIST:
                    node_0 = cellsHier_lev[c_ind][H_CN][0]
                    r_0 = self.nodesPoints[node_0]
                    TayCo[lev][c_ind] = self.GetTaylorCoeffAllPresetDers(f, r_0+r_corner, dx_lev, \
                        f_ders_lambda, factorials, scale_factor[lev])
        #print('TayCo:', TayCo)
        return TayCo


    def GetTaylorCoeff_1D_PresetDers(self, n_dir, r_0, f_ders_lambda, factorials, scale_factor):
        n = self.inds_keep.shape[n_dir]
        taylor = np.zeros(n)
        for i in range(n):
            taylor[i] = f_ders_lambda[i](r_0)/factorials[i]*scale_factor[i]
        return taylor


    def GetTaylorCoeffAllPresetDers_SeperableAllDir(self, r_0, f_ders_lambda_list, factorials_list, scale_factor_list):
        """ seperable in all directions
        """
        N = self.N
        taylor_1d = [None]*N
        for n_dir in range(N):
            taylor_1d[n_dir] = self.GetTaylorCoeff_1D_PresetDers(n_dir, r_0[n_dir], f_ders_lambda_list[n_dir], factorials_list[n_dir], scale_factor_list[n_dir])
        
        mask = self.inds_keep
        shape = mask.shape
        
        taylor = np.ones(shape)
        #inds = self.inds_0_allp        
        
        for n_dir in range(N):
            indx = [Ellipsis]*N
            for j in range(shape[n_dir]):
                indx[n_dir] = j
                taylor[indx] *= taylor_1d[n_dir][j]
        taylor *= mask
        return taylor


    def SetTaylorCoeffAlllevels_SeperableAllDir(self, f_list):
        """ seperable in all directions
        f_0: function of 0_th indep variable
        f_1: function of 1_th indep variable
        ...
        """
        N = self.N
        cellsHier = self.cellsHier
        n_lev = len(cellsHier)
        xyz_list = self.indepVars_list
        mask = self.inds_keep
        shape = mask.shape
        r_corner = self.x0
        factorials_list = [None]*N
        f_ders_lambda_list = [None]*N
        
        for n_dir in range(N):
            factorials = np.zeros(shape[n_dir])
            f_ders_lambda = np.zeros(shape[n_dir], dtype=object)
            f = f_list[n_dir]
            for i in range(shape[n_dir]):
                factorials[i] = math.factorial(i)
                f_i = (Derivative(f, xyz_list[n_dir], i)).doit()
                f_i_lambda = lambdify(xyz_list[n_dir], f_i)
                f_ders_lambda[i] = f_i_lambda
            factorials_list[n_dir] = factorials
            f_ders_lambda_list[n_dir] = f_ders_lambda
        
        ##set scale factor
        scale_factor_list = [None]*n_lev
        for lev in range(n_lev):
            scale_factor_list[lev] = [None]*N
            for n_dir in range(N):
                scale_factor = np.zeros(shape[n_dir])
                h = self.dx_levels[lev][n_dir]            
                for i in range(shape[n_dir]):
                    h_sc = h**i
                    scale_factor[i] = h_sc
                scale_factor_list[lev][n_dir] = scale_factor
        
        #print('scale_factor:', scale_factor)
        ##set taylor
        TayCo = [None]*n_lev
        for lev in range(n_lev):
            cellsHier_lev = cellsHier[lev]
            n_cell_lev = len(cellsHier_lev)
            TayCo[lev] = [None]*n_cell_lev
            dx_lev = self.dx_levels[lev]
            for c_ind in range(n_cell_lev):
                if cellsHier_lev[c_ind][H_CC]==NEXIST:
                    node_0 = cellsHier_lev[c_ind][H_CN][0]
                    r_0 = self.nodesPoints[node_0]
                    TayCo[lev][c_ind] = self.GetTaylorCoeffAllPresetDers_SeperableAllDir(r_0+r_corner, \
                        f_ders_lambda_list, factorials_list, scale_factor_list[lev])
        #print('TayCo:', TayCo)
        return TayCo

    def TestSeperableTaylor(self, f, f_list, printall=False):
        """ f is seperable in all directions
        """
        T0 = self.SetTaylorCoeffAlllevels(f)
        T1 = self.SetTaylorCoeffAlllevels_SeperableAllDir(f_list)
        cellsHier = self.cellsHier
        n_lev = len(cellsHier)
        for lev in range(n_lev):
            for c_ind in range(len(cellsHier[lev])):
                if cellsHier[lev][c_ind][H_CC]==NEXIST:
                    diff_rel = np.max(np.abs(T0[lev][c_ind] - T1[lev][c_ind]))/np.max(np.abs(T0[lev][c_ind]))
                    if printall:
                        print('T0[lev][c_ind]:', T0[lev][c_ind], 'T1[lev][c_ind]:', T1[lev][c_ind], sep='\n')
                        print('relative difference: ', diff_rel)
                    if diff_rel>1.0e-8:
                        print('WARNING!!!')
                        print('T0[lev][c_ind]:', T0[lev][c_ind], 'T1[lev][c_ind]:', T1[lev][c_ind], sep='\n')
                        print('relative difference: ', diff_rel)
                        

    def GetTaylorCoeff_Face(self, f, r_0, h, ijk_1_tup, n_dir):
        """ r_0: corner
            h : scaling factor
            f : symbolic function
            ijk_1_tup: taylor index on the face
            n_dir: on which face
        """
        N = self.N
        N_1 = N-1
        xyz_list = self.indepVars_list
        xyz_1_list = [xyz_list[i] for i in range(N) if i!=n_dir]
        arg_der = []
        for i in range(N_1):
            arg_der.append(xyz_1_list[i])
            arg_der.append(int(ijk_1_tup[i]))
        f_ijk = Derivative(f, *tuple(arg_der)).doit()
        f_ijk_lambda = lambdify(tuple(xyz_1_list), f_ijk)
        fact_1_ijk = 1.0
        for i in range(N_1):
            fact_1_ijk *= math.factorial(ijk_1_tup[i])
        r_0f = tuple([r_0[i] for i in range(N) if i!=n_dir])
        #print('f_ijk:', f_ijk, 'xyz_1_list:', xyz_1_list, 'f_ijk_lambda:', f_ijk_lambda, \
        #    'ijk_1_tup:', ijk_1_tup, 'r_0:', r_0, 'r_0f:', r_0f, sep='\n')
        tailor_ijk = f_ijk_lambda(*r_0f)/fact_1_ijk
        h_1 = np.array([h[i] for i in range(N) if i!=n_dir])
        h_1_sc = h_1**ijk_1_tup
        scale = 1.0
        for s in h_1_sc:
            scale *= s
        tailor_ijk_scaled = tailor_ijk*scale
        return tailor_ijk_scaled
        

    def GetTaylorCoeffAll_face(self, f, r_0, h, n_dir, flat=False):
        N = self.N
        N_1 = N-1
        mask = self.inds_keep
        shape = mask.shape
        shape_1 = tuple([shape[i] for i in range(N) if i!=n_dir])
        
        indx = [Ellipsis]*N
        indx[n_dir] = 0
        mask_1 = mask[indx]
        
        inds_0_allp = self.inds_0_allp        
        taylor = np.zeros(shape_1)
        inds = -np.ones(shape_1)
        
        counter = np.zeros(N, dtype=int)
        cnt_mask = [n_dir]
        while True:
            cnt_ind = tuple(counter.tolist())
            cnt_ind_1 = tuple([int(counter[i]) for i in range(N) if i!=n_dir])
            if mask_1[cnt_ind_1]==1:
                taylor[cnt_ind_1] = self.GetTaylorCoeff_Face(f, r_0, h, cnt_ind_1, n_dir)
                inds[cnt_ind_1] = inds_0_allp[cnt_ind]
            if not self._increaseCellCounterIndex_Masked(counter, shape, cnt_mask):
                break
        
        if flat:
            taylor_F, inds_F = [], []
            counter = np.zeros(N, dtype=int)
            cnt_mask = [n_dir]
            while True:
                cnt_ind = tuple(counter.tolist())
                cnt_ind_1 = tuple([int(counter[i]) for i in range(N) if i!=n_dir])
                if mask_1[cnt_ind_1]==1:
                    taylor_F.append(taylor[cnt_ind_1])
                    inds_F.append(inds_0_allp[cnt_ind])
                if not self._increaseCellCounterIndex_Masked(counter, shape, cnt_mask):
                    break
            #print('f:', f, '\n r_0:', r_0, 'h:', h, 'n_dir:', n_dir)
            #print('taylor_F:', taylor_F, 'inds_F:', inds_F, sep='\n')
            return [taylor_F, inds_F]
        else:
            #print('f:', f, '\n r_0:', r_0, 'h:', h, 'n_dir:', n_dir)
            #print('taylor:', taylor, 'inds:', inds, sep='\n')
            return [taylor, inds]
        

    def GetTaylorCoeffAllPresetDers_face(self, f, r_0, n_dir, f_ders_lambda, factorials, scale_factor, flat=False):
        N = self.N
        N_1 = N-1
        mask = self.inds_keep
        shape = mask.shape
        shape_1 = tuple([shape[i] for i in range(N) if i!=n_dir])
        
        indx = [Ellipsis]*N
        indx[n_dir] = 0
        mask_1 = mask[indx]
        
        inds_0_allp = self.inds_0_allp        
        taylor = np.zeros(shape_1)
        inds = -np.ones(shape_1)
        
        r_0f = tuple([r_0[i] for i in range(N) if i!=n_dir])
        counter = np.zeros(N, dtype=int)
        cnt_mask = [n_dir]
        while True:
            cnt_ind = tuple(counter.tolist())
            cnt_ind_1 = tuple([int(counter[i]) for i in range(N) if i!=n_dir])
            if mask_1[cnt_ind_1]==1:
                tailor_ijk = f_ders_lambda[cnt_ind_1](*r_0f)/factorials[cnt_ind_1]
                tailor_ijk_scaled = tailor_ijk*scale_factor[cnt_ind_1]
                taylor[cnt_ind_1] = tailor_ijk_scaled
                inds[cnt_ind_1] = inds_0_allp[cnt_ind]
            if not self._increaseCellCounterIndex_Masked(counter, shape, cnt_mask):
                break
        
        if flat:
            taylor_F, inds_F = [], []
            counter = np.zeros(N, dtype=int)
            cnt_mask = [n_dir]
            while True:
                cnt_ind = tuple(counter.tolist())
                cnt_ind_1 = tuple([int(counter[i]) for i in range(N) if i!=n_dir])
                if mask_1[cnt_ind_1]==1:
                    taylor_F.append(taylor[cnt_ind_1])
                    inds_F.append(inds_0_allp[cnt_ind])
                if not self._increaseCellCounterIndex_Masked(counter, shape, cnt_mask):
                    break
            #print('f:', f, '\n r_0:', r_0, 'h:', h, 'n_dir:', n_dir)
            #print('taylor_F:', taylor_F, 'inds_F:', inds_F, sep='\n')
            return [taylor_F, inds_F]
        else:
            #print('f:', f, '\n r_0:', r_0, 'h:', h, 'n_dir:', n_dir)
            #print('taylor:', taylor, 'inds:', inds, sep='\n')
            return [taylor, inds]


    def SetTaylorCoeffAlllevels_face(self, f, face, n_dir, flat=False):
        N = self.N
        N_1 = N-1

        cellsHier = self.cellsHier
        n_lev = len(cellsHier)
        
        mask = self.inds_keep
        shape = mask.shape
        shape_1 = tuple([shape[i] for i in range(N) if i!=n_dir])

        indx = [Ellipsis]*N
        indx[n_dir] = 0
        mask_1 = mask[indx]
        
        inds_0_allp = self.inds_0_allp        
        taylor = np.zeros(shape_1)
        inds = -np.ones(shape_1)
        
        xyz_list = self.indepVars_list
        xyz_1_list = [xyz_list[i] for i in range(N) if i!=n_dir]
        
        factorials = np.zeros(shape_1)
        f_ders_lambda = np.zeros(shape_1, dtype=object)
        r_corner = self.x0
        ##set derivatives and factoriels
        counter = np.zeros(N, dtype=int)
        cnt_mask = [n_dir]
        while True:
            cnt_ind = tuple(counter.tolist())
            cnt_ind_1 = tuple([int(counter[i]) for i in range(N) if i!=n_dir])
            if mask_1[cnt_ind_1]==1:
                arg_der = []
                for i in range(N_1):
                    arg_der.append(xyz_1_list[i])
                    arg_der.append(int(cnt_ind_1[i]))
                f_ijk = Derivative(f, *tuple(arg_der)).doit()
                f_ijk_lambda = lambdify(tuple(xyz_1_list), f_ijk)
                f_ders_lambda[cnt_ind_1] = f_ijk_lambda

                fact_1_ijk = 1.0
                for i in range(N_1):
                    fact_1_ijk *= math.factorial(cnt_ind_1[i])
                
                factorials[cnt_ind_1] = fact_1_ijk
            if not self._increaseCellCounterIndex_Masked(counter, shape, cnt_mask):
                break
        ##set scale factor
        scale_factor = [None]*n_lev
        for lev in range(n_lev):
            scale_factor[lev] = np.zeros(shape_1)
            h = self.dx_levels[lev]            
            h_1 = np.array([h[i] for i in range(N) if i!=n_dir])
            counter = np.zeros(N, dtype=int)
            while True:
                cnt_ind = tuple(counter.tolist())
                cnt_ind_1 = tuple([int(counter[i]) for i in range(N) if i!=n_dir])
                if mask_1[cnt_ind_1]==1:
                    h_1_sc = h_1**cnt_ind_1
                    scale = 1.0
                    for s in h_1_sc:
                        scale *= s
                    scale_factor[lev][cnt_ind_1] = scale
                if not self._increaseCellCounterIndex_Masked(counter, shape, cnt_mask):
                    break
        ##
        bcTayco = [None]*n_lev
        for lev in range(n_lev):
            cells_lev = cellsHier[lev]
            bc_conns_lev_face_dir = self.boundaryCellConnections[lev][face][n_dir]
            bcTayco[lev] = [None]*len(bc_conns_lev_face_dir)
            for i in range(len(bc_conns_lev_face_dir)):
                c_ind = bc_conns_lev_face_dir[i]
                if cells_lev[c_ind][H_CI]>=0:
                    node_0 = cells_lev[c_ind][H_CN][0]
                    r_0 = self.nodesPoints[node_0]
                    bcTayco[lev][i] = self.GetTaylorCoeffAllPresetDers_face(f, r_0+r_corner, n_dir, \
                        f_ders_lambda, factorials, scale_factor[lev], flat)
        return bcTayco
        

    def MultiplyPolynomialsNDArrays(self, poly_order, P0, P1):
        N = self.N
        mask = self.inds_keep
        shape = mask.shape
        assert P0.shape==P1.shape==shape
        counter = np.zeros(N, dtype=int)
        P0P1 = np.zeros(shape)
        while True:
            if mask[tuple(counter.tolist())]>0:
                a = P1[tuple(counter.tolist())]*P0
                for n in range(N):
                    a = np.roll(a, counter[n], axis=n)
                    for i in range(counter[n]):
                        indx = [Ellipsis]*N
                        indx[n] = i
                        a[indx] *= 0
                P0P1 += a
            if not self._increasePolyCounterIndex(counter, poly_order):
                break
        P0P1*=mask
        return P0P1


    def GetPolyMultTemplate(self, poly_order):
        N = self.N
        mask = self.inds_keep
        shape = mask.shape
        P0 = np.copy(self.inds_0_allp)
        P1 = np.copy(self.inds_0_allp)
        counter = np.zeros(N, dtype=int)
        n_F = len(self.inds_0_all_F)
        P0P1 = [None]*n_F
        for i in range(n_F):
            P0P1[i] = []
        while True:
            a_coeff = P1[tuple(counter.tolist())]
            if a_coeff>=0:
                a = np.copy(P0)
                for n in range(N):
                    a = np.roll(a, counter[n], axis=n)
                    for i in range(counter[n]):
                        indx = [Ellipsis]*N
                        indx[n] = i
                        a[indx] = -1
                for i in range(n_F):
                    i_a = self.inds_0_all_F[i][1]
                    if a[i_a]>=0:
                        P0P1[i].append((a[i_a], a_coeff))
            if not self._increasePolyCounterIndex(counter, poly_order):
                break
        P0P1_M = [None]*n_F     ## multiplicities determined
        for i in range(n_F):
            P0P1_M[i] = []
        for i in range(n_F):
            for j in range(len(P0P1[i])):
                duplicate = False
                for k in range(len(P0P1_M[i])):
                    if P0P1[i][j]==P0P1_M[i][k][0]:
                        P0P1_M[i][k][1] += 1
                        duplicate = True
                        break        
                if not duplicate:
                    P0P1_M[i].append([P0P1[i][j], 1])
        if self.verbose>0:
            print('P0P1:', P0P1, 'P0P1_M:', P0P1_M, sep='\n')
        self.P0P1_M = P0P1_M
        return P0P1_M


    def MultiplyPolynomialsNDArrays_Face(self, poly_order, P0, P1, n_dir):
        N = self.N
        N_1 = N-1
        mask = self.inds_keep
        shape = mask.shape
        shape_1 = tuple([shape[i] for i in range(N) if i!=n_dir])
        assert P0.shape==P1.shape==shape
        
        poly_order_1 = [poly_order[i] for i in range(N) if i!=n_dir]

        indx = [Ellipsis]*N
        indx[n_dir] = 0
        mask_1 = mask[indx]

        counter = np.zeros(N_1, dtype=int)
        P0P1 = np.zeros(shape_1)
        while True:
            if mask_1[tuple(counter.tolist())]:
                a = P1[tuple(counter.tolist())]*P0
                for n in range(N_1):
                    a = np.roll(a, counter[n], axis=n)
                    for i in range(counter[n]):
                        indx = [Ellipsis]*N_1
                        indx[n] = i
                        a[indx] *= 0
                P0P1 += a
            if not self._increasePolyCounterIndex(counter, poly_order_1):
                break
        P0P1*=mask_1
        return P0P1
        

    def GetNumOfUnknownsInEachCell(self):
        n_unkn_v = self.GetNumOfUnknownsInEachCellForEachVar()
        N_unkn = sum(n_unkn_v)
        return N_unkn
        
    def GetNumOfUnknownsInEachCellForEachVar(self):
        n_vars = len(self.vars_list)
        n_unkn_v = [0]*n_vars
        for v in range(n_vars):
            ind_to_rem = self.GetIndsToRemove(self.poly_ords[v], self.orders_max[v])
            ind_to_keep = (ind_to_rem==0)*1
            shape = ind_to_keep.shape
            n_tot = 1 
            for i in range(self.N):
                n_tot *= shape[i]
            ind_to_keep = ind_to_keep.reshape(n_tot, order='C')
            n_unkn_v[v] = sum(ind_to_keep)
        if self.verbose>0:
            print('n_unkn_v: ', n_unkn_v)
        self.n_unkn_v = n_unkn_v
        return n_unkn_v


    def GetNumOfUnknowns(self):
        n_cell_ind = self.GetNumOfIndexedCells()
        n_unkn_cell = self.GetNumOfUnknownsInEachCell()
        return n_cell_ind*n_unkn_cell



    ##TODO: to be extended to non-linear differential equations
    def SetupMatrixEQs(self):
        ## ** Works for ICs but gives singular matrices for BCs **
        """ take the field continuity on the boundaries, as initial conditions for
        the next cell. i.e. the initial condition for each cell is satisfied through
        the continuity condition at x=0, y=0... and the continuity condition at
        x=1 or y=1.. serves as initial condition for the cell on right or above.
        In case of boundary conditions the conditions at x=x0 and x=x1 serve as 
        filling the requirements of the leftmost (or downmost) cells.. 
        
        choice 1: create the final sparce matrix directly
        
        choice 2: put each sparce row inside the corressponding polyHier hierarchical 
        list, and finally gather them in a matrix
            - the indices can be either local (cell index, poly order index) or 
            final (final index number)
            PROS: may be more efficient for multi-grid methods
        """
        ##we assume all variables have the same polynomial order
        ##TODO: possibility of using unequal orders to be analyzed
        N = self.N
        cellsHier = self.cellsHier
        n_lev = len(cellsHier)
        poly_orders = self.poly_ords[0]
        orders_max = self.orders_max[0]     ##max poly orders to keep
        der_orders_max = self.der_orders_max[0]
        
        cell_n_tot__vars = self.GetNumOfUnknownsInEachCellForEachVar()
        if self.verbose>0:
            print('cell_n_tot__vars: ', cell_n_tot__vars)
        cell_n_tot = sum(cell_n_tot__vars)
        vars_list = self.vars_list
        n_vars = len(vars_list)
        cell_n_tot_vi = cell_n_tot__vars[0]
        for i in range(n_vars):
            assert cell_n_tot__vars[i]==cell_n_tot_vi
            
        ##cumulative number of indexed cells        
        if not self.N_indexed_acc:
            self.GetAccumukativeNumOfIndexedCells()
        N_Indexed_Cumul = self.N_indexed_acc
        """
        N_Indexed_Cumul = [0]*(n_lev+1)
        for lev in range(n_lev):
            N_Indexed_Cumul[lev+1] = N_Indexed_Cumul[lev]+self.indexInfo[lev]['N']
        self.N_Indexed_Cumul = N_Indexed_Cumul
        """
        
        ##get the total number of unknowns
        n_total = self.GetNumOfUnknowns()
        
        ##TODO:calculate the total number of elements to preallocate arrays
        n_elem = None 

        ##get the sparse matrix components 
        row, col, data = [], [], []

        ##---
        self.SetEqIndicesToKeep(poly_orders, orders_max)

        
        ##set matrix components from interior interactions
        if self.verbose>0:
            print('self.EQs_parts:', self.EQs_parts)
        for eq_ind, eq_arr in enumerate(self.EQs_parts):
            #eq_arr: a differential eq in the input set of differential equations
            if self.verbose>0:
                print('eq_ind:', eq_ind, '  eq_arr:', eq_arr)
            for v_ind, eq_v in enumerate(eq_arr):
                #eq_v: the differential equation associated with the v-th variable
                if self.verbose>0:
                    print('v_ind:', v_ind, '  eq_v:', eq_v)
                for eq_part in eq_v:
                    #eq_part: an additive part of the differential equation
                    if self.verbose>0:
                        print('eq_part:', eq_part)
                    coeff_sym, der_ord, der_vars = eq_part
                    
                    display(Math(latex(eq_part)+' : ' + latex(coeff_sym) + \
                    ';' + latex(der_ord) + ';' + latex(der_vars)))
                    
                    der_orders = self.ChangeDerOrderAndDerVarsToStandardForm(der_ord, der_vars)
                    
                    poly_diff = [None]*n_lev
                    inds_diff = [None]*n_lev
                    inds_orig = [None]*n_lev
                    for lev in range(n_lev):
                        poly_diff[lev], inds_diff[lev], inds_orig[lev] = \
                        self.GetPolyArrayDiff_scaled_indices_masked_flat(eq_ind, poly_orders, der_orders, orders_max, lev)
                        #print('lev:', lev, 'poly_diff:', poly_diff[lev], 'inds_diff:',\
                        # inds_diff[lev], 'inds_orig:', inds_orig[lev], sep='\n')
                    
                    if self.verbose>0:
                        print('coeff_sym:', coeff_sym)
                    coeff_sym = sympify(coeff_sym)
                    coeff_is_constant = 0
                    for par in self.pars_list:
                        if coeff_sym.has(par):
                            coeff_is_constant = 1
                            break
                    if coeff_is_constant==0:
                        for x_indep in self.indepVars_list:
                            if coeff_sym.has(x_indep):
                                coeff_is_constant = 2
                                break

                    if coeff_is_constant==0:
                        ##constant coefficients
                        if self.verbose>0:
                            print('coeff_sym: ', coeff_sym)
                        coeff = complex(coeff_sym)
                        for lev in range(n_lev):
                            ind_lev_st = N_Indexed_Cumul[lev]
                            cells_lev = self.cellsHier[lev]
                            
                            for i in range(len(cells_lev)):
                                if cells_lev[i][H_CI]>=0:
                                    c_ind_tot_st_row = (ind_lev_st + cells_lev[i][H_CI])*cell_n_tot \
                                        + eq_ind*cell_n_tot_vi
                                    c_ind_tot_st_col = (ind_lev_st + cells_lev[i][H_CI])*cell_n_tot \
                                        + v_ind*cell_n_tot_vi
                                    rows_new = c_ind_tot_st_row + inds_orig[lev]
                                    cols_new = c_ind_tot_st_col + inds_diff[lev]
                                    data_new = coeff*poly_diff[lev]
                                    
                                    row.extend(rows_new)
                                    col.extend(cols_new)
                                    data.extend(data_new)
                                    
                            assert np.all(inds_orig[lev]<cell_n_tot_vi)
                            assert np.all(inds_diff[lev]<cell_n_tot_vi)
                    elif coeff_is_constant==1 or coeff_is_constant==2:
                        if self.P0P1_M==None:
                            self.GetPolyMultTemplate(poly_orders)
                            
                        poly_diff_inds_mul = [None]*n_lev
                        inds_orig_nc = [None]*n_lev         ## _nc : non constant
                        for lev in range(n_lev):
                            poly_diff_inds_mul[lev], inds_orig_nc[lev] = \
                            self.GetPolyArrayDiff_scaled_indices_masked_flat_multTemp(eq_ind, poly_orders, \
                                der_orders, orders_max, lev)
                        
                        ##set taylor coeffs for coeff_sym
                        coeff_sym_sub = coeff_sym
                        if self.verbose>0:
                            print('coeff_sym:', coeff_sym)
                        if coeff_is_constant==1:
                            for i_p in range(len(self.pars_list)):
                                if self.verbose>0:
                                    print('self.pars_list[i_p]:', self.pars_list[i_p], 'self.pars_values[i_p]:', self.pars_values[i_p])
                                coeff_sym_sub = coeff_sym_sub.subs(self.pars_list[i_p], self.pars_values[i_p])
                                assert coeff_sym_sub.has(self.pars_list[i_p])==False
                        if self.verbose>0:
                            print('coeff_sym_sub:', coeff_sym_sub)
                        coeff_sym_taylor = self.SetTaylorCoeffAlllevels(coeff_sym_sub)
                        #print('coeff_sym_sub:', coeff_sym_sub, 'coeff_sym_taylor:', coeff_sym_taylor, sep='\n')

                        for lev in range(n_lev):
                            ind_lev_st = N_Indexed_Cumul[lev]
                            cells_lev = self.cellsHier[lev]
                            
                            for i in range(len(cells_lev)):
                                if cells_lev[i][H_CI]>=0:
                                    c_ind_tot_st_row = (ind_lev_st + cells_lev[i][H_CI])*cell_n_tot \
                                        + eq_ind*cell_n_tot_vi
                                    c_ind_tot_st_col = (ind_lev_st + cells_lev[i][H_CI])*cell_n_tot \
                                        + v_ind*cell_n_tot_vi

                                    P0_tay = coeff_sym_taylor[lev][i]
                                    n_el0 = len(inds_orig_nc[lev])
                                    poly_diff_lev = np.zeros(n_el0) 
                                    poly_mul_lev = poly_diff_inds_mul[lev]
                                    inds_orig_nc_lev = inds_orig_nc[lev]
                                    assert len(poly_mul_lev)==n_el0
                                    col_coeff = [None]*n_el0
                                    #print('P0_F:', P0_F, 'len(P0_F):', len(P0_F))
                                    for j in range(n_el0):
                                        col_j = [None]*len(poly_mul_lev[j])
                                        for k in range(len(poly_mul_lev[j])):
                                            ij_01, multip = poly_mul_lev[j][k]
                                            ind_P0_tup = self.inds_0_all_F[ij_01[0]][1]
                                            col_j[k] = (ij_01[1], P0_tay[ind_P0_tup]*multip)
                                        col_coeff[j] = col_j
                                        
                                    n_el = sum([len(col_coeff[j]) for j in range(len(col_coeff))])
                                    rows_0, cols_0, data_0 = [None]*n_el, [None]*n_el, [None]*n_el
                                    _ind_=0
                                    for j in range(n_el0):
                                        for k in range(len(col_coeff[j])):
                                            rows_0[_ind_] = inds_orig_nc_lev[j]
                                            cols_0[_ind_] = col_coeff[j][k][0]
                                            data_0[_ind_] = col_coeff[j][k][1]
                                            _ind_ += 1
                                    assert _ind_==n_el
                                            
                                    rows_0, cols_0, data_0 = np.array(rows_0), np.array(cols_0), np.array(data_0)
                                    
                                    c_ind_tot_st_row = (ind_lev_st + cells_lev[i][H_CI])*cell_n_tot \
                                        + eq_ind*cell_n_tot_vi
                                    c_ind_tot_st_col = (ind_lev_st + cells_lev[i][H_CI])*cell_n_tot \
                                        + v_ind*cell_n_tot_vi
                                    rows_new = c_ind_tot_st_row + rows_0
                                    cols_new = c_ind_tot_st_col + cols_0
                                    data_new = data_0
                                    
                                    row.extend(rows_new)
                                    col.extend(cols_new)
                                    data.extend(data_new)
                    else:    
                        raise NotImplementedError()


        ##TODO: for simplicity assuming all variables are on equal footing in terms
        ## of the derivative orders, to be generalized (poly_orders may need to be
        ## different for each var otherwise)..

        ## continuity conditions across cells boundaries
        if not self.cells_nb_dic:
            self.GetNeigborCellTypesInTermsOfConnectedNodes()
        assert self.cells_nb_dic!=None
        if not self.cc_dic:
            self.GetCellContCondEqs_scaled_masked(poly_orders, polyorders_max=orders_max, der_ord_max=der_orders_max)
        assert self.cc_dic!=None
        if not self.cc_nbmul_dic:
            self.SetCCCellMultiplicity()
        assert self.cc_nbmul_dic!=None
        if not self.cellsCC:
            self.MapBCsToICsForEachBoundaryCell(poly_orders, self.der_orders_max)
        CCs = self.CCs
        cells_bound_marked = self.MarkBoundaryCells(info='simple')
        for var in CCs:
            v_ind = self.VarIndex[var]    ## variable index (order)
            v_ccs = CCs[var] ## variable's continuity conditions    
            for v_cc in v_ccs:
                n_dir, d_ord, face = v_cc['dir'], v_cc['der'], v_cc['face']
                #print('v_ccs: ', v_ccs)
                if face=='n':
                    ## 'sl'
                    for lev in range(n_lev):
                        ind_lev_st = N_Indexed_Cumul[lev]
                        cells_lev = cellsHier[lev]
                        cells_nb_n_sl_lev_dir = self.cells_nb_dic['n']['sl'][lev][n_dir]
                        if cells_nb_n_sl_lev_dir==None:
                            continue
                        if self.cc_dic['n']['sl'][lev][n_dir]==None:
                            continue
                        arr_inds_dir, arr_inds_dir_nb = self.cc_dic['n']['sl'][lev][n_dir][d_ord]
                        arr_I_FL_f0_dir, inds_0_FL_f0_dir = arr_inds_dir
                        arr_I_FL_f0_nb_dir, inds_0_FL_f0_nb_dir = arr_inds_dir_nb
                        
                        cc_multip, cc_eq_index = self.cc_eqindex_and_multiplicity[n_dir][d_ord]
                        assert len(arr_I_FL_f0_dir)==len(cc_eq_index)
                        
                        ccnb_multip = self.cc_nbmul_dic['n']['sl'][n_dir]   ## neigbor multiplicity
                        cc_multip = [m*ccnb_multip for m in cc_multip]
                        
                        n_el = sum([len(el) for el in arr_I_FL_f0_dir])
                        n_el_nb = sum([len(el) for el in arr_I_FL_f0_nb_dir])
                        rows_0, cols_0, data_0 = [None]*n_el, [None]*n_el, [None]*n_el
                        rows_nb_0, cols_nb_0, data_nb_0 = [None]*n_el_nb, [None]*n_el_nb, [None]*n_el_nb
                        ind_0 = 0
                        for i in range(len(arr_I_FL_f0_dir)):
                            for col_ind in arr_I_FL_f0_dir[i]:
                                rows_0[ind_0] = cc_eq_index[i]
                                cols_0[ind_0] = col_ind
                                data_0[ind_0] = arr_I_FL_f0_dir[i][col_ind]/cc_multip[i]
                                ind_0 += 1
                        assert ind_0==n_el
                        ind_0 = 0
                        for i in range(len(arr_I_FL_f0_nb_dir)):
                            for col_ind in arr_I_FL_f0_nb_dir[i]:
                                ##rows are set according to cell, not cell_nb
                                rows_nb_0[ind_0] = cc_eq_index[i] 
                                cols_nb_0[ind_0] = col_ind
                                data_nb_0[ind_0] = -arr_I_FL_f0_nb_dir[i][col_ind]/cc_multip[i]
                                ind_0 += 1
                        assert ind_0==n_el_nb
                        
                        rows_0, cols_0, data_0 = np.array(rows_0), np.array(cols_0), np.array(data_0)
                        rows_nb_0, cols_nb_0, data_nb_0 = np.array(rows_nb_0), np.array(cols_nb_0), np.array(data_nb_0)
                        
                        #print('rows_0:', rows_0)
                        assert np.all(rows_0<cell_n_tot)
                        
                        for i in range(len(cells_nb_n_sl_lev_dir)):
                            c_ind, c_nb_ind = cells_nb_n_sl_lev_dir[i]
                            assert cells_lev[c_ind][H_CC]==NEXIST and cells_lev[c_nb_ind][H_CC]==NEXIST
                            assert cells_lev[c_ind][H_CI]>=0 and cells_lev[c_nb_ind][H_CI]>=0
                            if not cells_bound_marked[lev][c_ind]:
                                c_ind_tot_st = (ind_lev_st + cells_lev[c_ind][H_CI])*cell_n_tot \
                                    + v_ind*cell_n_tot_vi
                                
                                rows_new = c_ind_tot_st + rows_0
                                cols_new = c_ind_tot_st + cols_0
                                data_new = data_0
                                        
                                row.extend(rows_new)
                                col.extend(cols_new)
                                data.extend(data_new)

                                c_nb_ind_tot_st = (ind_lev_st + cells_lev[c_nb_ind][H_CI])*cell_n_tot \
                                    + v_ind*cell_n_tot_vi

                                rows_new = c_ind_tot_st + rows_nb_0
                                cols_new = c_nb_ind_tot_st + cols_nb_0
                                data_new = data_nb_0

                                row.extend(rows_new)
                                col.extend(cols_new)
                                data.extend(data_new)

                            else:
                                ## boundary cell
                                mults_bc_F, inds_bc_0 = self.cellsCC[v_ind][lev][n_dir][d_ord][c_ind]
                                #print('inds_bc_0:', inds_bc_0, 'cc_eq_index:', cc_eq_index)
                                assert inds_bc_0==cc_eq_index
                                
                                cc_multip_bc = [m*ccnb_multip for m in mults_bc_F]
                                
                                data_0_bc = [None]*n_el
                                data_nb_0_bc = [None]*n_el_nb
                                ind_0 = 0
                                for i in range(len(arr_I_FL_f0_dir)):
                                    for col_ind in arr_I_FL_f0_dir[i]:
                                        data_0_bc[ind_0] = data_0[ind_0]*cc_multip[i]/cc_multip_bc[i]
                                        ind_0 += 1
                                assert ind_0==n_el
                                ind_0 = 0
                                for i in range(len(arr_I_FL_f0_nb_dir)):
                                    for col_ind in arr_I_FL_f0_nb_dir[i]:
                                        data_nb_0_bc[ind_0] = data_nb_0[ind_0]*cc_multip[i]/cc_multip_bc[i]
                                        ind_0 += 1
                                assert ind_0==n_el_nb
                                
                                data_0_bc = np.array(data_0_bc)
                                data_nb_0_bc = np.array(data_nb_0_bc)

                                c_ind_tot_st = (ind_lev_st + cells_lev[c_ind][H_CI])*cell_n_tot \
                                    + v_ind*cell_n_tot_vi
                                
                                rows_new = c_ind_tot_st + rows_0
                                cols_new = c_ind_tot_st + cols_0
                                data_new = data_0_bc
                                        
                                row.extend(rows_new)
                                col.extend(cols_new)
                                data.extend(data_new)

                                c_nb_ind_tot_st = (ind_lev_st + cells_lev[c_nb_ind][H_CI])*cell_n_tot \
                                    + v_ind*cell_n_tot_vi

                                rows_new = c_ind_tot_st + rows_nb_0
                                cols_new = c_nb_ind_tot_st + cols_nb_0
                                data_new = data_nb_0_bc

                                row.extend(rows_new)
                                col.extend(cols_new)
                                data.extend(data_new)

                    ## 'nl'
                    for lev in range(n_lev):
                        ind_lev_st = N_Indexed_Cumul[lev]
                        ind_levp1_st = None
                        if lev<n_lev-1:
                            ind_levp1_st = N_Indexed_Cumul[lev+1]
                        cells_lev = cellsHier[lev]
                        cells_levp1 = None
                        if lev<n_lev-1:
                            cells_levp1 = cellsHier[lev+1]
                        cells_nb_n_nl_lev_dir = self.cells_nb_dic['n']['nl'][lev][n_dir]
                        if cells_nb_n_nl_lev_dir==None:
                            continue
                        for n_conn in cells_nb_n_nl_lev_dir:
                            arr_inds_dir, arr_inds_dir_nb = self.cc_dic['n']['nl'][lev][n_dir][n_conn][d_ord]
                            arr_I_FL_f0_dir, inds_0_FL_f0_dir = arr_inds_dir
                            arr_I_FL_f0_nb_dir, inds_0_FL_f0_nb_dir = arr_inds_dir_nb
                            
                            cc_multip, cc_eq_index = self.cc_eqindex_and_multiplicity[n_dir][d_ord]
                            assert len(arr_I_FL_f0_dir)==len(cc_eq_index)
                            
                            ccnb_multip = self.cc_nbmul_dic['n']['nl'][n_dir]   ## neigbor multiplicity
                            cc_multip = [m*ccnb_multip for m in cc_multip]
                            
                            n_el = sum([len(el) for el in arr_I_FL_f0_dir])
                            n_el_nb = sum([len(el) for el in arr_I_FL_f0_nb_dir])
                            rows_0, cols_0, data_0 = [None]*n_el, [None]*n_el, [None]*n_el
                            rows_nb_0, cols_nb_0, data_nb_0 = [None]*n_el_nb, [None]*n_el_nb, [None]*n_el_nb
                            ind_0 = 0
                            for i in range(len(arr_I_FL_f0_dir)):
                                for col_ind in arr_I_FL_f0_dir[i]:
                                    rows_0[ind_0] = cc_eq_index[i]
                                    cols_0[ind_0] = col_ind
                                    data_0[ind_0] = arr_I_FL_f0_dir[i][col_ind]/cc_multip[i]
                                    ind_0 += 1
                            assert ind_0==n_el
                            ind_0 = 0
                            for i in range(len(arr_I_FL_f0_nb_dir)):
                                for col_ind in arr_I_FL_f0_nb_dir[i]:
                                    ##rows are set according to cell, not cell_nb
                                    rows_nb_0[ind_0] = cc_eq_index[i] 
                                    cols_nb_0[ind_0] = col_ind
                                    data_nb_0[ind_0] = -arr_I_FL_f0_nb_dir[i][col_ind]/cc_multip[i]
                                    ind_0 += 1
                            assert ind_0==n_el_nb
                        
                            rows_0, cols_0, data_0 = np.array(rows_0), np.array(cols_0), np.array(data_0)
                            rows_nb_0, cols_nb_0, data_nb_0 = np.array(rows_nb_0), np.array(cols_nb_0), np.array(data_nb_0)

                            for i in range(len(cells_nb_n_nl_lev_dir[n_conn])):
                                c_ind, c_nb_ind = cells_nb_n_nl_lev_dir[n_conn][i]
                                assert cells_lev[c_ind][H_CC]==NEXIST and cells_levp1[c_nb_ind][H_CC]==NEXIST
                                assert cells_lev[c_ind][H_CI]>=0 and cells_levp1[c_nb_ind][H_CI]>=0
                                if not cells_bound_marked[lev][c_ind]:
                                    c_ind_tot_st = (ind_lev_st + cells_lev[c_ind][H_CI])*cell_n_tot \
                                        + v_ind*cell_n_tot_vi
                                    
                                    rows_new = c_ind_tot_st + rows_0
                                    cols_new = c_ind_tot_st + cols_0
                                    data_new = data_0
                                            
                                    row.extend(rows_new)
                                    col.extend(cols_new)
                                    data.extend(data_new)

                                    c_nb_ind_tot_st = (ind_levp1_st + cells_levp1[c_nb_ind][H_CI])*cell_n_tot \
                                        + v_ind*cell_n_tot_vi

                                    rows_new = c_ind_tot_st + rows_nb_0
                                    cols_new = c_nb_ind_tot_st + cols_nb_0
                                    data_new = data_nb_0

                                    row.extend(rows_new)
                                    col.extend(cols_new)
                                    data.extend(data_new)
                                else:
                                    ## boundary cell
                                    mults_bc_F, inds_bc_0 = self.cellsCC[v_ind][lev][n_dir][d_ord][c_ind]
                                    #print('inds_bc_0:', inds_bc_0, 'cc_eq_index:', cc_eq_index)
                                    assert inds_bc_0==cc_eq_index
                                    
                                    cc_multip_bc = [m*ccnb_multip for m in mults_bc_F]
                                
                                    data_0_bc = [None]*n_el
                                    data_nb_0_bc = [None]*n_el_nb
                                    ind_0 = 0
                                    for i in range(len(arr_I_FL_f0_dir)):
                                        for col_ind in arr_I_FL_f0_dir[i]:
                                            data_0_bc[ind_0] = data_0[ind_0]*cc_multip[i]/cc_multip_bc[i]
                                            ind_0 += 1
                                    assert ind_0==n_el
                                    ind_0 = 0
                                    for i in range(len(arr_I_FL_f0_nb_dir)):
                                        for col_ind in arr_I_FL_f0_nb_dir[i]:
                                            data_nb_0_bc[ind_0] = data_nb_0[ind_0]*cc_multip[i]/cc_multip_bc[i]
                                            ind_0 += 1
                                    assert ind_0==n_el_nb
                                
                                    data_0_bc = np.array(data_0_bc)
                                    data_nb_0_bc = np.array(data_nb_0_bc)
                    
                                    c_ind_tot_st = (ind_lev_st + cells_lev[c_ind][H_CI])*cell_n_tot \
                                        + v_ind*cell_n_tot_vi
                                    
                                    rows_new = c_ind_tot_st + rows_0
                                    cols_new = c_ind_tot_st + cols_0
                                    data_new = data_0_bc
                                            
                                    row.extend(rows_new)
                                    col.extend(cols_new)
                                    data.extend(data_new)

                                    c_nb_ind_tot_st = (ind_levp1_st + cells_levp1[c_nb_ind][H_CI])*cell_n_tot \
                                        + v_ind*cell_n_tot_vi

                                    rows_new = c_ind_tot_st + rows_nb_0
                                    cols_new = c_nb_ind_tot_st + cols_nb_0
                                    data_new = data_nb_0_bc

                                    row.extend(rows_new)
                                    col.extend(cols_new)
                                    data.extend(data_new)
                    ## 'pl'
                    for lev in range(n_lev):
                        ind_lev_st = N_Indexed_Cumul[lev]
                        cells_lev = cellsHier[lev]
                        cells_levm1 = None
                        if lev>0:
                            cells_levm1 = cellsHier[lev-1]
                        ind_levm1_st = None
                        if lev>0:
                            ind_levm1_st = N_Indexed_Cumul[lev-1]
                        cells_nb_n_pl_lev_dir = self.cells_nb_dic['n']['pl'][lev][n_dir]
                        if cells_nb_n_pl_lev_dir==None:
                            continue
                        for n_conn in cells_nb_n_pl_lev_dir:
                            arr_inds_dir, arr_inds_dir_nb = self.cc_dic['n']['pl'][lev][n_dir][n_conn][d_ord]
                            arr_I_FL_f0_dir, inds_0_FL_f0_dir = arr_inds_dir
                            arr_I_FL_f0_nb_dir, inds_0_FL_f0_nb_dir = arr_inds_dir_nb
                            
                            cc_multip, cc_eq_index = self.cc_eqindex_and_multiplicity[n_dir][d_ord]
                            assert len(arr_I_FL_f0_dir)==len(cc_eq_index)
                            
                            ccnb_multip = self.cc_nbmul_dic['n']['pl'][n_dir]   ## neigbor multiplicity
                            cc_multip = [m*ccnb_multip for m in cc_multip]
                            
                            n_el = sum([len(el) for el in arr_I_FL_f0_dir])
                            n_el_nb = sum([len(el) for el in arr_I_FL_f0_nb_dir])
                            rows_0, cols_0, data_0 = [None]*n_el, [None]*n_el, [None]*n_el
                            rows_nb_0, cols_nb_0, data_nb_0 = [None]*n_el_nb, [None]*n_el_nb, [None]*n_el_nb
                            ind_0 = 0
                            for i in range(len(arr_I_FL_f0_dir)):
                                for col_ind in arr_I_FL_f0_dir[i]:
                                    rows_0[ind_0] = cc_eq_index[i]
                                    cols_0[ind_0] = col_ind
                                    data_0[ind_0] = arr_I_FL_f0_dir[i][col_ind]/cc_multip[i]
                                    ind_0 += 1
                            assert ind_0==n_el
                            ind_0 = 0
                            for i in range(len(arr_I_FL_f0_nb_dir)):
                                for col_ind in arr_I_FL_f0_nb_dir[i]:
                                    ##rows are set according to cell, not cell_nb
                                    rows_nb_0[ind_0] = cc_eq_index[i] 
                                    cols_nb_0[ind_0] = col_ind
                                    data_nb_0[ind_0] = -arr_I_FL_f0_nb_dir[i][col_ind]/cc_multip[i]
                                    ind_0 += 1
                            assert ind_0==n_el_nb
                        
                            rows_0, cols_0, data_0 = np.array(rows_0), np.array(cols_0), np.array(data_0)
                            rows_nb_0, cols_nb_0, data_nb_0 = np.array(rows_nb_0), np.array(cols_nb_0), np.array(data_nb_0)

                            for i in range(len(cells_nb_n_pl_lev_dir[n_conn])):
                                c_ind, c_nb_ind = cells_nb_n_pl_lev_dir[n_conn][i]
                                assert cells_lev[c_ind][H_CI]>=0 and cells_levm1[c_nb_ind][H_CI]>=0
                                if not cells_bound_marked[lev][c_ind]:
                                    c_ind_tot_st = (ind_lev_st + cells_lev[c_ind][H_CI])*cell_n_tot \
                                        + v_ind*cell_n_tot_vi
                                    
                                    rows_new = c_ind_tot_st + rows_0
                                    cols_new = c_ind_tot_st + cols_0
                                    data_new = data_0
                                            
                                    row.extend(rows_new)
                                    col.extend(cols_new)
                                    data.extend(data_new)

                                    c_nb_ind_tot_st = (ind_levm1_st + cells_levm1[c_nb_ind][H_CI])*cell_n_tot \
                                        + v_ind*cell_n_tot_vi

                                    rows_new = c_ind_tot_st + rows_nb_0
                                    cols_new = c_nb_ind_tot_st + cols_nb_0
                                    data_new = data_nb_0

                                    row.extend(rows_new)
                                    col.extend(cols_new)
                                    data.extend(data_new)
                                else:
                                    ## boundary cell
                                    mults_bc_F, inds_bc_0 = self.cellsCC[v_ind][lev][n_dir][d_ord][c_ind]
                                    #print('inds_bc_0:', inds_bc_0, 'cc_eq_index:', cc_eq_index)
                                    assert inds_bc_0==cc_eq_index
                                    
                                    cc_multip_bc = [m*ccnb_multip for m in mults_bc_F]

                                    data_0_bc = [None]*n_el
                                    data_nb_0_bc = [None]*n_el_nb
                                    ind_0 = 0
                                    for i in range(len(arr_I_FL_f0_dir)):
                                        for col_ind in arr_I_FL_f0_dir[i]:
                                            data_0_bc[ind_0] = data_0[ind_0]*cc_multip[i]/cc_multip_bc[i]
                                            ind_0 += 1
                                    assert ind_0==n_el
                                    ind_0 = 0
                                    for i in range(len(arr_I_FL_f0_nb_dir)):
                                        for col_ind in arr_I_FL_f0_nb_dir[i]:
                                            data_nb_0_bc[ind_0] = data_nb_0[ind_0]*cc_multip[i]/cc_multip_bc[i]
                                            ind_0 += 1
                                    assert ind_0==n_el_nb
                                
                                    data_0_bc = np.array(data_0_bc)
                                    data_nb_0_bc = np.array(data_nb_0_bc)

                                    c_ind_tot_st = (ind_lev_st + cells_lev[c_ind][H_CI])*cell_n_tot \
                                        + v_ind*cell_n_tot_vi
                                    
                                    rows_new = c_ind_tot_st + rows_0
                                    cols_new = c_ind_tot_st + cols_0
                                    data_new = data_0_bc
                                            
                                    row.extend(rows_new)
                                    col.extend(cols_new)
                                    data.extend(data_new)

                                    c_nb_ind_tot_st = (ind_levm1_st + cells_levm1[c_nb_ind][H_CI])*cell_n_tot \
                                        + v_ind*cell_n_tot_vi

                                    rows_new = c_ind_tot_st + rows_nb_0
                                    cols_new = c_nb_ind_tot_st + cols_nb_0
                                    data_new = data_nb_0_bc

                                    row.extend(rows_new)
                                    col.extend(cols_new)
                                    data.extend(data_new)
                    
                else:
                    raise NotImplementedError()

        
        ##boundary conditions
        BC_parts, BC_rhs = self.DisintegrateBoundaryConds()
        if not self.bc_dic:
            self.GetCellBoundCondEqs_scaled_masked(poly_orders, polyorders_max=orders_max, der_ord_max=der_orders_max)
        assert self.bc_dic!=None
        assert self.bc_eqindex_and_multiplicity!=None
        if not self.bctoic_map:
            self.MapBCsToICs(self.der_orders_max)
        assert self.bctoic_map!=None
        bctoic_map = self.bctoic_map
        
        for bc_ind, b_cond in enumerate(self.BCs):
            n_dir, expr, face = b_cond['dir'], b_cond['expr'], b_cond['face']
            
            n_dir_eqind, v_ind_eqind, der_ord_eqind = bctoic_map[bc_ind]['dir'],\
                bctoic_map[bc_ind]['v_ind'], bctoic_map[bc_ind]['d_ord']
            
            bc_eq_arr = BC_parts[bc_ind]
            for v_ind in bc_eq_arr:
                eq_v = bc_eq_arr[v_ind]
                #eq_v: the differential equation associated with the v-th variable
                for eq_part in eq_v:
                    #eq_part: an additive part of the differential equation
                    coeff_sym, der_ord, der_vars = eq_part
                    
                    der_orders = self.ChangeDerOrderAndDerVarsToStandardForm(der_ord, der_vars)
                    d_ord = der_orders[n_dir]
                    for i in range(N):
                        if i != n_dir:
                            assert der_orders[i]==0

                    coeff_sym = sympify(coeff_sym)
                    coeff_is_constant = True
                    for par in self.pars_list:
                        if coeff_sym.has(par):
                            coeff_is_constant = False
                            break
                    if coeff_is_constant:
                        ##constant coefficients
                        coeff = complex(coeff_sym)
                        for lev in range(n_lev):
                            ind_lev_st = N_Indexed_Cumul[lev]
                            cells_lev = cellsHier[lev]
                            bc_conns_lev_n_dir = self.boundaryCellConnections[lev]['n'][n_dir]
                            bc_conns_lev_p_dir = self.boundaryCellConnections[lev]['p'][n_dir]
                            
                            if face=='n':
                                arr_I_FL_f0_dir, inds_0_FL_f0_dir = self.bc_dic['n'][n_dir][lev][d_ord]
                                
                                bc_multip, bc_eq_index = self.bc_eqindex_and_multiplicity[n_dir_eqind][der_ord_eqind]
                                assert len(bc_eq_index)==len(inds_0_FL_f0_dir)
                                
                                n_el = sum([len(el) for el in arr_I_FL_f0_dir])
                                rows_0, cols_0, data_0 = [None]*n_el, [None]*n_el, [None]*n_el
                                ind_0 = 0
                                for i in range(len(arr_I_FL_f0_dir)):
                                    for col_ind in arr_I_FL_f0_dir[i]:
                                        rows_0[ind_0] = bc_eq_index[i]
                                        cols_0[ind_0] = col_ind
                                        data_0[ind_0] = arr_I_FL_f0_dir[i][col_ind]/bc_multip[i]
                                        ind_0 += 1
                                assert ind_0==n_el

                                rows_0, cols_0, data_0 = np.array(rows_0), np.array(cols_0), np.array(data_0)
                                #print('lev:', lev, 'n_dir:', n_dir, 'v_ind:', v_ind, 'face:', face, 'der_orders:', der_orders)
                                #print('rows_0:', rows_0, 'cols_0:', cols_0, 'data_0:', data_0, '-'*20, sep='\n')

                                for i in range(len(bc_conns_lev_n_dir)):
                                    c_ind = bc_conns_lev_n_dir[i]
                                    
                                    if cells_lev[c_ind][H_CI]>=0:
                                        mults_F_bc, eqinds_F_bc, inds_0_bc = self.cellsBC[v_ind][lev]['n'][n_dir][d_ord][c_ind]
                                        assert -1 not in eqinds_F_bc
                                        #print('inds_0_bc:', inds_0_bc, 'inds_0_FL_f0_dir:', inds_0_FL_f0_dir, 'bc_eq_index:', bc_eq_index)
                                        assert inds_0_bc==bc_eq_index
                                    
                                        assert len(eqinds_F_bc)==len(arr_I_FL_f0_dir)
                                        rows_0_bc, data_0_bc = [None]*n_el, [None]*n_el
                                        ind_0 = 0
                                        for i in range(len(arr_I_FL_f0_dir)):
                                            for col_ind in arr_I_FL_f0_dir[i]:
                                                rows_0_bc[ind_0] = eqinds_F_bc[i]
                                                data_0_bc[ind_0] = data_0[ind_0]*bc_multip[i]/mults_F_bc[i]
                                                ind_0 += 1
                                        assert ind_0==n_el


                                        #c_ind_tot_st_row = (ind_lev_st + cells_lev[c_ind][H_CI])*cell_n_tot \
                                        #    + v_ind_eqind*cell_n_tot_vi
                                        c_ind_tot_st_col = (ind_lev_st + cells_lev[c_ind][H_CI])*cell_n_tot \
                                            + v_ind*cell_n_tot_vi
                                        
                                        rows_new = rows_0_bc
                                        cols_new = c_ind_tot_st_col + cols_0
                                        data_new = data_0_bc
                                                
                                        row.extend(rows_new)
                                        col.extend(cols_new)
                                        data.extend(data_new)
                            else:
                                assert face=='p'
                                arr_I_FL_f0_dir, inds_0_FL_f0_dir = self.bc_dic['p'][n_dir][lev][d_ord]
                                
                                bc_multip, bc_eq_index = self.bc_eqindex_and_multiplicity[n_dir_eqind][der_ord_eqind]
                                assert len(bc_eq_index)==len(inds_0_FL_f0_dir)

                                n_el = sum([len(el) for el in arr_I_FL_f0_dir])
                                rows_0, cols_0, data_0 = [None]*n_el, [None]*n_el, [None]*n_el
                                ind_0 = 0
                                for i in range(len(arr_I_FL_f0_dir)):
                                    for col_ind in arr_I_FL_f0_dir[i]:
                                        rows_0[ind_0] = bc_eq_index[i]
                                        cols_0[ind_0] = col_ind
                                        data_0[ind_0] = arr_I_FL_f0_dir[i][col_ind]/bc_multip[i]
                                        ind_0 += 1
                                assert ind_0==n_el

                                rows_0, cols_0, data_0 = np.array(rows_0), np.array(cols_0), np.array(data_0)

                                for i in range(len(bc_conns_lev_p_dir)):
                                    c_ind = bc_conns_lev_p_dir[i]
                                    c_nb_ind = bc_conns_lev_n_dir[i]
                                    
                                    if cells_lev[c_ind][H_CI]>=0:
                                        assert cells_lev[c_nb_ind][H_CI]>=0

                                        mults_F_bc, eqinds_F_bc, inds_0_bc = self.cellsBC[v_ind][lev]['p'][n_dir][d_ord][c_ind]
                                        assert -1 not in eqinds_F_bc
                                        #print('n_dir:', n_dir, 'd_ord:', d_ord, 'inds_0_bc:', inds_0_bc, 'inds_0_FL_f0_dir:', \
                                        #    inds_0_FL_f0_dir, 'bc_eq_index:', bc_eq_index)
                                        #print('mults_F_bc:', mults_F_bc, 'eqinds_F_bc:', eqinds_F_bc)
                                        #assert inds_0_bc==bc_eq_index
                                    
                                        assert len(eqinds_F_bc)==len(arr_I_FL_f0_dir)
                                        rows_0_bc, data_0_bc = [None]*n_el, [None]*n_el
                                        ind_0 = 0
                                        for i in range(len(arr_I_FL_f0_dir)):
                                            for col_ind in arr_I_FL_f0_dir[i]:
                                                rows_0_bc[ind_0] = eqinds_F_bc[i]
                                                data_0_bc[ind_0] = data_0[ind_0]*bc_multip[i]/mults_F_bc[i]
                                                ind_0 += 1
                                        assert ind_0==n_el

                                        #c_ind_tot_st_row = (ind_lev_st + cells_lev[c_nb_ind][H_CI])*cell_n_tot \
                                        #    + v_ind_eqind*cell_n_tot_vi
                                        c_ind_tot_st_col = (ind_lev_st + cells_lev[c_ind][H_CI])*cell_n_tot \
                                            + v_ind*cell_n_tot_vi
                                        
                                        rows_new = rows_0_bc
                                        cols_new = c_ind_tot_st_col + cols_0
                                        data_new = data_0_bc
                                                
                                        row.extend(rows_new)
                                        col.extend(cols_new)
                                        data.extend(data_new)

                                    
                    else:
                        raise NotImplementedError()


        ##test
        if self.verbose>0:
            c_startend_inds = [0]
            for lev in range(n_lev):
                c_startend_inds.append(c_startend_inds[lev]+self.indexInfo[lev]['N'])
            ##oor: out of range indices (rows)
            rows_oor = [] 
            for i in range(len(row)):
                if row[i]>=n_total:
                    rows_oor.append(row[i])
            print('rows_oor[{}]:'.format(len(rows_oor)), rows_oor)
            cells_oor = [None]*len(rows_oor)
            for i, r_ind in enumerate(rows_oor):
                c_lev = None
                for lev in range(n_lev):
                    if r_ind>=c_startend_inds[lev]*cell_n_tot and r_ind<c_startend_inds[lev+1]*cell_n_tot:
                        c_lev = lev
                        r_ind_resi = r_ind - c_startend_inds[lev]*cell_n_tot
                        c_ind = int(r_ind_resi/cell_n_tot)
                        cells_oor[i] = (c_lev, c_ind)
                        break
            cells_oor = list(set(cells_oor))
            print('cells_oor[{}]:'.format(len(cells_oor)), cells_oor)
            ##oor: out of range indices (cols)
            rows_oor = [] 
            for i in range(len(row)):
                if col[i]>=n_total:
                    rows_oor.append(row[i])
            print('cols_oor[{}]:'.format(len(rows_oor)), rows_oor)
            cells_oor = [None]*len(rows_oor)
            for i, r_ind in enumerate(rows_oor):
                c_lev = None
                for lev in range(n_lev):
                    if r_ind>=c_startend_inds[lev]*cell_n_tot and r_ind<c_startend_inds[lev+1]*cell_n_tot:
                        c_lev = lev
                        r_ind_resi = r_ind - c_startend_inds[lev]*cell_n_tot
                        c_ind = int(r_ind_resi/cell_n_tot)
                        cells_oor[i] = (c_lev, c_ind)
                        break
            cells_oor = list(set(cells_oor))
            print('cells_oor[{}]:'.format(len(cells_oor)), cells_oor)
            ## unset rows
            print('n_total:', n_total)     
            rows_marked = [0]*n_total
            for i in range(len(row)):
                rows_marked[row[i]] = 1
            rows_unset = [i for i in range(n_total) if rows_marked[i]==0]
            print('rows_unset[{}]:'.format(len(rows_unset)), rows_unset)
            cells_unset = [None]*len(rows_unset)
            for i, r_ind in enumerate(rows_unset):
                c_lev = None
                for lev in range(n_lev):
                    if r_ind>=c_startend_inds[lev]*cell_n_tot and r_ind<c_startend_inds[lev+1]*cell_n_tot:
                        c_lev = lev
                        r_ind_resi = r_ind - c_startend_inds[lev]*cell_n_tot
                        c_ind = int(r_ind_resi/cell_n_tot)
                        cells_unset[i] = (c_lev, c_ind)
                        break
            cells_unset = list(set(cells_unset))
            print('cells_unset[{}]:'.format(len(cells_unset)), cells_unset)
            ## unset columns
            cols_marked = [0]*n_total
            for i in range(len(col)):
                cols_marked[col[i]] = 1
            cols_unset = [i for i in range(n_total) if cols_marked[i]==0]
            print('cols_unset[{}]:'.format(len(cols_unset)), cols_unset)
            cells_unset = [None]*len(cols_unset)
            for i, r_ind in enumerate(cols_unset):
                c_lev = None
                for lev in range(n_lev):
                    if r_ind>=c_startend_inds[lev]*cell_n_tot and r_ind<c_startend_inds[lev+1]*cell_n_tot:
                        c_lev = lev
                        r_ind_resi = r_ind - c_startend_inds[lev]*cell_n_tot
                        c_ind = int(r_ind_resi/cell_n_tot)
                        cells_unset[i] = (c_lev, c_ind)
                        break
            cells_unset = list(set(cells_unset))
            print('cells_unset[{}]:'.format(len(cells_unset)), cells_unset)
            ## rows with zero entries   
            rows_zero = [] 
            for i in range(len(row)):
                if data[i]==0:
                    rows_zero.append(row[i])
            print('rows_zero[{}]:'.format(len(rows_zero)), rows_zero)
            cells_zero = [None]*len(rows_zero)
            for i, r_ind in enumerate(rows_zero):
                c_lev = None
                for lev in range(n_lev):
                    if r_ind>=c_startend_inds[lev]*cell_n_tot and r_ind<c_startend_inds[lev+1]*cell_n_tot:
                        c_lev = lev
                        r_ind_resi = r_ind - c_startend_inds[lev]*cell_n_tot
                        c_ind = int(r_ind_resi/cell_n_tot)
                        cells_zero[i] = (c_lev, c_ind)
                        break
            cells_zero = list(set(cells_zero))
            print('cells_zero[{}]:'.format(len(cells_zero)), cells_zero)


        ##RHS
        data_rhs = np.zeros(n_total, dtype=complex)
        for i_eq in range(len(self.rhs)):
            for i_rp in range(len(self.rhs_parts[i_eq])):
                ##constant coefficients
                if self.rhs_parts[i_eq][i_rp][1]==CoeffType.const:
                    if self.verbose>0:
                        print('CoeffType.const::', 'self.rhs_parts[i_eq][i_rp]', self.rhs_parts[i_eq][i_rp])
                    coeff_sym = self.rhs_parts[i_eq][i_rp][0]
                    coeff = complex(coeff_sym)
                    for lev in range(n_lev):
                        ind_lev_st = N_Indexed_Cumul[lev]
                        cells_lev = self.cellsHier[lev]
                        
                        poly_rhs, inds_0 = self.GetPolyArrayRHSConst_masked_flat(i_eq, poly_orders, orders_max, coeff)
                        for i in range(len(cells_lev)):
                            if cells_lev[i][H_CI]>=0:
                                c_ind_tot_st = (ind_lev_st + cells_lev[i][H_CI])*cell_n_tot \
                                    + i_eq*cell_n_tot_vi    ## i_eq=ind_v
                                rows_new = c_ind_tot_st + inds_0
                                data_new = poly_rhs
                                
                                data_rhs[rows_new] += data_new  ##TODO: check multiplicity
                elif self.rhs_parts[i_eq][i_rp][1]==CoeffType.par_ext:
                    if self.verbose>0:
                        print('CoeffType.par_ext::', 'self.rhs_parts[i_eq][i_rp]', self.rhs_parts[i_eq][i_rp])
                    pass

                elif self.rhs_parts[i_eq][i_rp][1]==CoeffType.par_int or self.rhs_parts[i_eq][i_rp][1]==CoeffType.nonconst_anal:
                    coeff_sym = self.rhs_parts[i_eq][i_rp][0]
                    #coeff_sym_sub = coeff_sym
                    #if self.rhs_parts[i_eq][i_rp][1]==CoeffType.par_int:
                    #    for i_p in range(len(self.pars_list)):
                    #        coeff_sym_sub = coeff_sym_sub.subs(self.pars_list[i_p], self.pars_values[i_p])
                    #        assert coeff_sym_sub.has(self.pars_list[i_p])==False
                    #coeff_sym_taylor = self.SetTaylorCoeffAlllevels(coeff_sym_sub)
                    coeff_sym_taylor = self.SetTayCoAllLev_ParType(coeff_sym)

                    if self.verbose>0:
                        print('CoeffType.par_int or anal::', 'self.rhs_parts[i_eq][i_rp]', self.rhs_parts[i_eq][i_rp])
                        print('coeff_sym_sub:', coeff_sym_sub)
                        
                    for lev in range(n_lev):
                        ind_lev_st = N_Indexed_Cumul[lev]
                        cells_lev = self.cellsHier[lev]
                        
                        poly_I_rhs, inds_0 = self.GetPolyArrayRHSConst_masked_flat_Temp(i_eq, poly_orders, orders_max)
                        #print('poly_I_rhs:', poly_I_rhs, 'inds_0:', inds_0, 'orders_max:', orders_max, sep='\n')
                        for i in range(len(cells_lev)):
                            if cells_lev[i][H_CI]>=0:

                                coeff_vec = coeff_sym_taylor[lev][i]
                                poly_rhs = [None]*len(poly_I_rhs)
                                for j in range(len(poly_I_rhs)):
                                    _ind_, _w_ = poly_I_rhs[j]
                                    poly_rhs[j] = _w_*coeff_vec[self.inds_0_all_F[_ind_][1]]

                                c_ind_tot_st = (ind_lev_st + cells_lev[i][H_CI])*cell_n_tot \
                                    + i_eq*cell_n_tot_vi    ## i_eq=ind_v
                                rows_new = c_ind_tot_st + inds_0
                                data_new = poly_rhs
                                
                                data_rhs[rows_new] += data_new  ##TODO: check multiplicity
                    
                else:
                    raise NotImplementedError()

        ##BC RHS (boundary condition RHS)
        #print('self.BCs: ', self.BCs)
        for bc_ind, bc_cond in enumerate(self.BCs):
            #print('bc_cond:', bc_cond)
            n_dir, expr, face = bc_cond['dir'], bc_cond['expr'], bc_cond['face']
            
            face_eqind, n_dir_eqind, v_ind_eqind, d_ord_eqind = self.bcic_to_bcic_map[bc_ind]['face'], \
                self.bcic_to_bcic_map[bc_ind]['dir'], self.bcic_to_bcic_map[bc_ind]['v_ind'], self.bcic_to_bcic_map[bc_ind]['d_ord']
                                
            #print('bc_ind:', bc_ind, 'n_dir:', n_dir, 'expr:', expr, 'face:', face)
            #print('n_dir_eqind:', n_dir_eqind, 'v_ind_eqind:', v_ind_eqind, 'der_ord_eqind:', der_ord_eqind)
        
            bc_rhs = self.BC_rhs[bc_ind]
            coeff_sym = bc_rhs
            coeff_is_constant = 0  ## 0:constant 1:has par  2:has x,y,z..
            for par in self.pars_list:
                if coeff_sym.has(par):
                    coeff_is_constant = 1
                    break
            if coeff_is_constant==0:
                for x_indep in self.indepVars_list:
                    if coeff_sym.has(x_indep):
                        coeff_is_constant = 2
                        break
                
            #print('coeff_sym:', coeff_sym)
            if coeff_is_constant==0:
                ##constant rhs
                coeff = complex(coeff_sym)
                #print('coeff:', coeff)
                for lev in range(n_lev):
                    ind_lev_st = N_Indexed_Cumul[lev]
                    cells_lev = cellsHier[lev]
                    bc_conns_lev_n_dir = self.boundaryCellConnections[lev]['n'][n_dir]
                    bc_conns_lev_p_dir = self.boundaryCellConnections[lev]['p'][n_dir]
                    
                    if face_eqind=='n':

                        #coeff_vec = np.zeros(len(rows_0), dtype=complex)
                        #coeff_vec[0] = coeff    ##TODO: check consistensy
                        #assert self.inds_0_all_F[0] == [0, tuple([0]*N)]
                        #data_0 = np.ones(len(rows_0))*coeff_vec
                        #print('rows_0:', rows_0, 'data_0:', data_0, sep='\n')
                        
                        #print('bc_conns_lev_n_dir:', bc_conns_lev_n_dir)
                        #print('bc_conns_lev_p_dir:', bc_conns_lev_p_dir)
                        for i in range(len(bc_conns_lev_n_dir)):
                            c_ind = bc_conns_lev_n_dir[i]
                            
                            #print('lev:', lev, 'c_ind:', c_ind)
                            if cells_lev[c_ind][H_CI]>=0:
                                mults_F_bc, eqinds_F_bc, inds_0_bc = self.cellsBC[v_ind_eqind][lev]['n'][n_dir_eqind][d_ord_eqind][c_ind]
                                
                                coeff_vec = np.zeros(len(mults_F_bc), dtype=complex)
                                coeff_vec[0] = coeff
                                assert self.inds_0_all_F[0] == [0, tuple([0]*N)]
                                data_0 = np.ones(len(mults_F_bc))*coeff_vec
                                
                                rows_new = eqinds_F_bc
                                data_new = data_0/np.array(mults_F_bc)
                                        
                                data_rhs[rows_new] += data_new
                                #print('mults_F_bc:', mults_F_bc, 'eqinds_F_bc:',  eqinds_F_bc, 'inds_0_bc:', inds_0_bc)
                    if face_eqind=='p':
                        #coeff_vec = np.zeros(len(rows_0), dtype=complex)
                        #coeff_vec[0] = coeff
                        #data_0 = np.ones(len(rows_0))*coeff_vec
                        #print('rows_0:', rows_0, 'data_0:', data_0, sep='\n')
                        
                        for i in range(len(bc_conns_lev_n_dir)):
                            c_ind = bc_conns_lev_p_dir[i]
                            c_nb_ind = bc_conns_lev_n_dir[i]
                            
                            if cells_lev[c_ind][H_CI]>=0: 
                                assert cells_lev[c_nb_ind][H_CI]>=0
                                mults_F_bc, eqinds_F_bc, inds_0_bc = self.cellsBC[v_ind_eqind][lev]['p'][n_dir_eqind][d_ord_eqind][c_ind]

                                coeff_vec = np.zeros(len(mults_F_bc), dtype=complex)
                                coeff_vec[0] = coeff
                                assert self.inds_0_all_F[0] == [0, tuple([0]*N)]
                                data_0 = np.ones(len(mults_F_bc))*coeff_vec
                                
                                rows_new = eqinds_F_bc
                                data_new = data_0/np.array(mults_F_bc)
                                        
                                data_rhs[rows_new] += data_new

            elif coeff_is_constant==1 or coeff_is_constant==2:
                ##variable rhs
                if self.verbose>0:
                    print('coeff_sym:', coeff_sym)
                coeff_sym_sub = coeff_sym
                if coeff_is_constant==1:
                    for i_p in range(len(self.pars_list)):
                        coeff_sym_sub = coeff_sym_sub.subs(self.pars_list[i_p], self.pars_values[i_p])
                        assert coeff_sym_sub.has(self.pars_list[i_p])==False
                if self.verbose>0:
                    print('coeff_sym_sub:', coeff_sym_sub)
                for lev in range(n_lev):
                    ind_lev_st = N_Indexed_Cumul[lev]
                    cells_lev = cellsHier[lev]
                    bc_conns_lev_n_dir = self.boundaryCellConnections[lev]['n'][n_dir]
                    bc_conns_lev_p_dir = self.boundaryCellConnections[lev]['p'][n_dir]
                    
                    if face_eqind=='n':
                    
                        bcTayco = self.SetTaylorCoeffAlllevels_face(coeff_sym_sub, face_eqind, n_dir, flat=True)
                        
                        #print('bc_conns_lev_n_dir:', bc_conns_lev_n_dir)
                        #print('bc_conns_lev_p_dir:', bc_conns_lev_p_dir)
                        for i in range(len(bc_conns_lev_n_dir)):
                            c_ind = bc_conns_lev_n_dir[i]
                            
                            #print('lev:', lev, 'c_ind:', c_ind)
                            if cells_lev[c_ind][H_CI]>=0:
                                mults_F_bc, eqinds_F_bc, inds_0_bc = self.cellsBC[v_ind_eqind][lev]['n'][n_dir_eqind][d_ord_eqind][c_ind]

                                r_0 = self.nodesPoints[cellsHier[lev][c_ind][H_CN][0]]
                                h = self.dx_levels[lev]
                                #[taylor, inds] = self.GetTaylorCoeffAll_face(f=coeff_sym_sub, r_0=r_0+self.x0, h=h, n_dir=n_dir, flat=True)                        
                                [taylor, _inds_] = bcTayco[lev][i]
                                #assert _inds_==inds_0_bc ##TODO: check relation of _inds_ and inds_0_bc
                                #print('_inds_:', _inds_, 'inds_0_bc:', inds_0_bc, sep='\n')
                                data_0 = taylor/np.array(mults_F_bc)
                                                                
                                rows_new = eqinds_F_bc
                                data_new = data_0
                                        
                                data_rhs[rows_new] += data_new

                                #print('c_ind:', c_ind, 'ind_lev_st:', ind_lev_st, 'c_ind_tot_st_row:', c_ind_tot_st_row)
                                #print('rows_0:', rows_0, 'data_0:', data_0, sep='\n')

                    if face_eqind=='p':
                        
                        bcTayco = self.SetTaylorCoeffAlllevels_face(coeff_sym_sub, face_eqind, n_dir, flat=True)

                        for i in range(len(bc_conns_lev_n_dir)):
                            c_ind = bc_conns_lev_p_dir[i]
                            c_nb_ind = bc_conns_lev_n_dir[i]
                            
                            if cells_lev[c_ind][H_CI]>=0: 
                                mults_F_bc, eqinds_F_bc, inds_0_bc = self.cellsBC[v_ind_eqind][lev]['p'][n_dir_eqind][d_ord_eqind][c_ind]

                                r_0 = self.nodesPoints[cellsHier[lev][c_ind][H_CN][0]]
                                h = self.dx_levels[lev]
                                #[taylor, inds] = self.GetTaylorCoeffAll_face(f=coeff_sym_sub, r_0=r_0+self.x0, h=h, n_dir=n_dir, flat=True)                        
                                [taylor, _inds_] = bcTayco[lev][i]
                                #assert _inds_==inds_0_bc
                                #print('_inds_:', _inds_, 'inds_0_bc:', inds_0_bc, sep='\n')
                                data_0 = taylor/np.array(mults_F_bc)

                                assert cells_lev[c_nb_ind][H_CI]>=0
                                
                                rows_new = eqinds_F_bc
                                data_new = data_0
                                        
                                data_rhs[rows_new] += data_new

                                #print('c_ind:', c_ind, 'c_nb_ind:', c_nb_ind, 'ind_lev_st:', ind_lev_st, 'c_ind_tot_st_row:', c_ind_tot_st_row)
                                #print('rows_0:', rows_0, 'data_0:', data_0, sep='\n')
            else:
                raise NotImplementedError()

        ##construct the sparse matrix
        row = np.array(row)
        col = np.array(col)
        data = np.array(data)
        A_coo = coo_matrix((data, (row, col)), shape=(n_total,n_total), dtype=complex)
        
        ##test
        if self.verbose>0:
            print('min(matrix): ', min(abs(data)), '    max(matrix): ', max(abs(data)))
        #print('data_rhs', data_rhs)
        
        #for i in range(n_total):
        #    print(i, A_coo.getrow(i))
        #    print('-'*20)
        
        """
        rows_marked_same = [None]*n_total
        ind_same = 0
        for i in range(n_total):
            if rows_marked_same[i]!=None:
                continue
            row_i = A_coo.getrow(i).tocoo()
            is_same = False
            for j in range(i, n_total):
                row_j = A_coo.getrow(j).tocoo()
                if row_i.nnz!=row_j.nnz:
                    continue
                if np.all((row_i.col==row_j.col)) and np.all((row_i.data==row_j.data)):
                    is_same==True
                    rows_marked_same[i] = ind_same
                    rows_marked_same[j] = ind_same
            if is_same:
                ind_same += 1
        rows_same = [None]*ind_same
        for i in range(ind_same):
            rows_same[i] = []
            for j in range(n_total):
                if rows_marked_same[j]==ind_same:
                    rows_same[i].append(j)
        print('ind_same:', ind_same)
        for i in range(ind_same):
            print(rows_same[i])
            print(A_coo.getrow(rows_same[i][0]))
            print('-'*20)
        """
        
        
        ##solve the matrix
        if self.verbose>0:
            print('solving matrix..')
        from scipy.sparse import linalg

        #rank = np.linalg.matrix_rank(A_coo)
        #print('Matrix rank = ', rank)
        
        if self.matsolver=='UMF':
            x_res = linalg.spsolve(A_coo.tocsc(), data_rhs)
            self.mat_coo = A_coo
            self.mat_rhs = data_rhs
            self.mat_res = x_res
            return x_res
        elif self.matsolver=='SLU':        
            LU = linalg.splu(A_coo.tocsc())
            x_res = LU.solve(data_rhs)
            self.LU = LU
            self.mat_coo = A_coo
            self.mat_rhs = data_rhs
            self.mat_res = x_res
            return x_res
        else:
            raise NotImplementedError()
        

    def RegisterFieldsAtTheGivenSurface(self, v_ind, n_dir, face, der_order):
        [PolyFaceCoeffs, PolyFaceInds0] = self.SetSolvedPolyCoeffsOnFace(v_ind, face, n_dir, der_order)
    
        val_n_dir = None
        if face=='n':
            val_n_dir = self.x0[n_dir]
        elif face=='p':
            val_n_dir = self.x1[n_dir]
        der_order_tup = tuple(der_order)
            
        if (v_ind, n_dir, face, der_order_tup) not in self.FieldsAtSurf:
            self.FieldsAtSurf[(v_ind, n_dir, face, der_order_tup)] = deque(maxlen=self.que_maxlen)
            
        self.FieldsAtSurf[(v_ind, n_dir, face, der_order_tup)].append([val_n_dir, [PolyFaceCoeffs, PolyFaceInds0]])
        return


    def RegisterICFieldsAtTheGivenSurface(self, v_ind, n_dir, face, der_order):
        [PolyFaceCoeffs, PolyFaceInds0] = self.SetSolvedPolyCoeffsOnFace(v_ind, 'n', n_dir, der_order)
    
        val_n_dir = self.x0[n_dir]
        der_order_tup = tuple(der_order)

        if (v_ind, n_dir, face, der_order_tup) not in self.FieldsAtSurf:
            self.FieldsAtSurf[(v_ind, n_dir, face, der_order_tup)] = deque(maxlen=self.que_maxlen)
            
        self.FieldsAtSurf[(v_ind, n_dir, face, der_order_tup)].append([val_n_dir, [PolyFaceCoeffs, PolyFaceInds0]])
        return
        

    def GetLagrangeCoeffs(self, t, t_0=0.0, dt=None):
        """ t: list of float (sample points)
            dt: scale factor 
            t_0:poly origin, otherwise origin is 0
            convention:   a_n*((t-t_0)/dt)^n + a_(n-1)*((t-t_0)/dt)^n-1 + ... 
        """
        n = len(t)
        p = [None]*n
        for i in range(n):
            p_i = np.poly1d([1.0])
            for j in range(n):
                if j!=i:
                    p_i *= np.poly1d([1, t_0-t[j]])/(t[i]-t[j])
            p[i] = p_i.c.tolist()
            p[i].reverse()
            if dt!=None:
                for j in range(len(p[i])):
                    p[i][j] *= dt**j
        return p
        

    def ExtrapolateRegisteredFields(self, poly_order, v_ind, n_dir, face, der_order):
        assert face=='p'
        N = self.N
        face_op = 'n'
        cellsHier = self.cellsHier
        n_lev = len(cellsHier)
        
        inds_keep = self.inds_keep
        shape = inds_keep.shape
        inds_0_all_F = self.inds_0_all_F
        n_inds = len(inds_0_all_F)
        inds_0_all_F_1 = [inds_0_all_F[i][1] for i in range(n_inds)]
        
        n_p1 = poly_order[n_dir] + 1
        polyssaved = self.FieldsAtSurf[(v_ind, n_dir, face, tuple(der_order))]
        n_p1 = min(n_p1, len(polyssaved))
        
        n_p_sv1 = len(polyssaved)-1
        polyssaved_lagr = [None]*n_p1
        for i in range(n_p1):
            if n_p_sv1-i>=0:
                polyssaved_lagr[i] = polyssaved[n_p_sv1-i]
        assert None not in polyssaved_lagr
        
        t = [None]*n_p1
        dt = self.W[n_dir]
        for i in range(n_p1):
            t[i] = -float(i)*dt
            
        n_lg_max = 4
        c_lagr_arr = [None]*n_lev
        for lev in range(n_lev):
            c_lagr_arr[lev] = [None]*n_lg_max
            dt_lev = self.dx_levels[lev][n_dir]
            for i in range(n_lg_max):
                c_lagr_arr[lev][i] = self.GetLagrangeCoeffs(t, t_0=t[0]+i*dt_lev, dt=dt_lev)
                                
        fieldsExtp = [None]*n_lev
        for lev in range(n_lev):
            fieldsExtp[lev] = [None]*len(cellsHier[lev])
            cells_lev = cellsHier[lev]
            bc_conns_lev_face_dir = self.boundaryCellConnections[lev][face_op][n_dir]

            #dt_lev = self.dx_levels[lev][n_dir]
            #c_lagr = self.GetLagrangeCoeffs(t, t_0=t[0], dt=dt_lev)
            c_lagr = c_lagr_arr[lev][0]

            for i in range(len(bc_conns_lev_face_dir)):
                c_ind = bc_conns_lev_face_dir[i]
                if cells_lev[c_ind][H_CI]>=0:

                    p_extp = np.zeros(shape)
                    for i_ps in range(n_p1):        # _ps: poly saved
                        if  polyssaved_lagr[i_ps]==None:
                            continue
                        [PolyFaceCoeffs, PolyFaceInds0] = polyssaved_lagr[i_ps][1]
                        inds_0_FL_f0_dir = PolyFaceInds0[lev]
                        p_f = np.zeros(shape)
                        
                        n_term_F = len(inds_0_FL_f0_dir)
                        #print('inds_0_all_F:', inds_0_all_F, 'inds_0_FL_f0_dir:', inds_0_FL_f0_dir, sep='\n')
                        for j in range(n_term_F):
                            p_f[inds_0_all_F_1[inds_0_FL_f0_dir[j]]] = PolyFaceCoeffs[lev][i][j]
                        
                        indx_face = [Ellipsis]*N
                        indx_face[n_dir] = 0
                        p_f_face = p_f[indx_face]
                        #print('p_f:', p_f, 'p_f_face:', p_f_face, 'p_extp:', p_extp, sep='\n')
                        for j in range(n_p1):
                            indx = [Ellipsis]*N
                            indx[n_dir] = j
                            p_extp[indx] += c_lagr[i_ps][j]*p_f_face

                    fieldsExtp[lev][c_ind] = np.zeros(n_inds)
                    fieldsExtp_lev_c_ind = fieldsExtp[lev][c_ind]
                    for _i_ in range(n_inds):
                        fieldsExtp_lev_c_ind[_i_] = p_extp[inds_0_all_F_1[_i_]]
                    
                    _c_ind_ = c_ind
                    #_t_0_ = t[0]
                    _ind_lg = 1
                    while True:
                        c_nb_p_dir = self.CellGetFaceConnectedCellsSameLevel(lev, _c_ind_)[1][n_dir]
                        if c_nb_p_dir!=None:
                            ##adjust p_extp TODO: optimize for speed using a template..
                            assert cellsHier[lev][c_ind][H_CI]>=0

                            #_c_lagr_ = self.GetLagrangeCoeffs(t, t_0=_t_0_+dt_lev, dt=dt_lev)
                            _c_lagr_ = c_lagr_arr[lev][_ind_lg]
                            p_extp = np.zeros(shape)
                            for i_ps in range(n_p1):
                                if  polyssaved_lagr[i_ps]==None:
                                    continue
                                [PolyFaceCoeffs, PolyFaceInds0] = polyssaved_lagr[i_ps][1]
                                inds_0_FL_f0_dir = PolyFaceInds0[lev]
                                p_f = np.zeros(shape)
                                
                                n_term_F = len(inds_0_FL_f0_dir)
                                for j in range(n_term_F):
                                    p_f[inds_0_all_F_1[inds_0_FL_f0_dir[j]]] = PolyFaceCoeffs[lev][i][j]
                                
                                indx_face = [Ellipsis]*N
                                indx_face[n_dir] = 0
                                p_f_face = p_f[indx_face]
                                for j in range(n_p1):
                                    indx = [Ellipsis]*N
                                    indx[n_dir] = j
                                    p_extp[indx] += _c_lagr_[i_ps][j]*p_f_face

                            fieldsExtp[lev][c_nb_p_dir] = np.zeros(n_inds)
                            fieldsExtp_lev_c_nb = fieldsExtp[lev][c_nb_p_dir]
                            for _i_ in range(n_inds):
                                fieldsExtp_lev_c_nb[_i_] = p_extp[inds_0_all_F_1[_i_]]
                            _c_ind_ = c_nb_p_dir
                            #_t_0_ += self.dx_levels[n_dir]
                            _ind_lg += 1
                        else:
                            _t_end_ = t[0] + _ind_lg*dt_lev
                            assert (_t_end_ -(t[0]+dt))/(t[0]+dt)<1.0e-5
                            break
                            
        for lev in range(n_lev):
            for c_ind in range(len(cellsHier[lev])):
                if cellsHier[lev][c_ind][H_CI]>=0:
                    assert fieldsExtp[lev][c_ind] != None
                        
        return fieldsExtp   ##flat
        

    def ExtrapolateRegisteredFieldsAll(self, poly_order, n_dir, face, der_order, setAsPolyHier=False):
    
        n_vars = len(self.vars_list)
        cellsHier = self.cellsHier
        n_lev = len(cellsHier)
                
        polyHier = [None]*n_vars
        for v_ind in range(n_vars):
            fieldsExtp = self.ExtrapolateRegisteredFields(poly_order, v_ind, n_dir, face, der_order)
            polyHier[v_ind] = fieldsExtp

        if setAsPolyHier:
            self.polyHier = polyHier
        return polyHier                


    def StorePolyHier(self, polyHier, x_res=None, r_corner=None):
        self.polyStore.append([polyHier, x_res, r_corner])
        return
        
        
    def setupPolyStorageQues(self, maxlen):
        self.polyStore = deque(maxlen=maxlen)
        

    def AdvanceCornerInVarDir(self, var_indep):
        n_dir = self.IndepVarIndex[var_indep]
        dr = np.zeros(self.N)
        dr[n_dir] = self.W[n_dir]
        self.x0 += dr
        self.x1 += dr
        
    ##TODO: if sources are functions of f(r)*g(t), update the taylor coeffs in 
    ## the t direction only
    
    ##TODO: for functions varying only in a small region of space, set the Taylor
    ## coeffs only for non-constant cells

    ##TODO: put all operations in the form of matrix multiplications where possible
    
    ##TODO: vectorize sympy function calcucations (for taylor coeffs)


    def SetTayCoAllLev_ParType(self, expr):
        """ expr: sympy expression
            expr does not contain pars_extern
        """
        N = self.N
        #get all pars inside
        pars_inside = []
        for ind_p, par in enumerate(self.pars_list):
            if expr.has(par):
                pars_inside.append(ind_p)
        
        for par_ext in self.pars_extern:
            if expr.has(par_ext):
                raise NotImplementedError()
        
        if len(pars_inside)==1:
            p_ind = pars_inside[0]
            if self.pars_types[p_ind]==ParType.seperable:
                par = self.pars_list[p_ind]
                eq_tree = symExpr_generate_tree(expr)
                coeff_sym, der_ord, der_vars = sym_getcoeff_setzero(eq_tree, par)
                assert eq_tree[1]==0
                if der_ord==0:
                    f_list = self.pars_values[p_ind]
                    if coeff_sym!=1:
                        f_list[0]*=coeff_sym
                    return self.SetTaylorCoeffAlllevels_SeperableAllDir(f_list)
                else:
                    raise NotImplementedError()
        
        if len(pars_inside)==0:
            expr_isconst = True
            for v_indep in self.indepVars_list:
                if expr.has(v_indep):
                    expr_isconst = False
                    break
            assert not expr_isconst
            return self.SetTaylorCoeffAlllevels_SeperableAllDir(expr)
            
        
        ##general - no optimization possible
        f_pars = [None]*pars_inside
        for i in range(len(pars_inside)): 
            p_ind = pars_inside[i]
            if self.pars_types[p_ind]==ParType.seperable:
                f_pars[i] = Integer(1)
                for n_dir in range(N):
                    f_pars[i] *= self.pars_values[p_ind][n_dir]
            elif self.pars_types[p_ind]==ParType.general:
                f_pars[i] = self.pars_types[p_ind]
            else:
                raise NotImplementedError()
        f = Integer(1)
        for i in range(len(pars_inside)): 
            f *= f_pars[i]
        return self.SetTaylorCoeffAlllevels_SeperableAllDir(f)
        


    def ResetRHSUseCCasIC_Leapfrog(self, bcinds_cc, cc_equiv, n_dir_lf, step_init=False, setMat_res=True):
        """ bcinds_cc: [bc_ind_0, bc_ind_1,..] bc_inds for which to use continuity condition 
        as boundary condition
            cc_equiv: what continituity condition to use for each bc_ind in bcinds_cc
        """            
        N = self.N
        n_total = self.GetNumOfUnknowns()
        cellsHier = self.cellsHier
        n_lev = len(cellsHier)

        poly_orders = self.poly_ords[0]
        orders_max = self.orders_max[0]     ##max poly orders to keep
        der_orders_max = self.der_orders_max[0]
        
        cell_n_tot__vars = self.GetNumOfUnknownsInEachCellForEachVar()
        if self.verbose>0:
            print('cell_n_tot__vars: ', cell_n_tot__vars)
        cell_n_tot = sum(cell_n_tot__vars)
        vars_list = self.vars_list
        n_vars = len(vars_list)
        cell_n_tot_vi = cell_n_tot__vars[0]
        for i in range(n_vars):
            assert cell_n_tot__vars[i]==cell_n_tot_vi

        N_Indexed_Cumul = self.N_indexed_acc

        P0P1_M = self.P0P1_M
        n_F = len(self.inds_0_all_F)
        
        ##RHS
        ##disintegrate self.rhs into coeff, derivative, par and determine 
        ##wether the parameter is in self.pars_list or self.pars_extern
        data_rhs = np.zeros(n_total, dtype=complex)
        for i_eq in range(len(self.rhs)):
            for i_rp in range(len(self.rhs_parts[i_eq])):

                ##constant coefficients
                if self.rhs_parts[i_eq][i_rp][1]==CoeffType.const:
                    coeff_sym = self.rhs_parts[i_eq][i_rp][0]
                    coeff = complex(coeff_sym)
                    for lev in range(n_lev):
                        ind_lev_st = N_Indexed_Cumul[lev]
                        cells_lev = self.cellsHier[lev]
                        
                        poly_rhs, inds_0 = self.GetPolyArrayRHSConst_masked_flat(i_eq, poly_orders, orders_max, coeff)
                        for i in range(len(cells_lev)):
                            if cells_lev[i][H_CI]>=0:
                                c_ind_tot_st = (ind_lev_st + cells_lev[i][H_CI])*cell_n_tot \
                                    + i_eq*cell_n_tot_vi    ## i_eq=ind_v
                                rows_new = c_ind_tot_st + inds_0
                                data_new = poly_rhs
                                
                                data_rhs[rows_new] += data_new  ##TODO: check multiplicity

                elif self.rhs_parts[i_eq][i_rp][1]==CoeffType.par_ext:  
                    if step_init:
                        continue    ##TODO: instead, safer to initialize rgnd_ext.polyHier to 0
                    par_ext, par_type, par_coeff_sym, der_ord, der_vars = self.rhs_parts[i_eq][i_rp]
                    par_coeff_sym = sympify(par_coeff_sym)
                    der_order = self.ChangeDerOrderAndDerVarsToStandardForm(der_ord, der_vars)
                    rgnd_ext = self.pars_extern[par_ext]
                    
                    polyHier_v_ext = None
                    use_direct_diff=False
                    if use_direct_diff:
                        ##TODO: NOT TESTED!!!
                        rgnd_ext.ExtrapolateRegisteredFieldsAll(poly_order=poly_orders, \
                            n_dir=n_dir_lf, face='p', der_order=der_order, setAsPolyHier=True)
                                                        
                        polyHier_v_ext = rgnd_ext.polyHier[rgnd_ext.VarIndex[par_ext]]
                    else:
                        _polyHier_ext_ = rgnd_ext.ExtrapolateRegisteredFieldsAll(poly_order=poly_orders, \
                            n_dir=n_dir_lf, face='p', der_order=tuple(np.zeros(N, dtype=int)), setAsPolyHier=True)
                        polyHier_v_ext = rgnd_ext.TakeDerivativePolyHier(_polyHier_ext_, \
                                v_ind=rgnd_ext.VarIndex[par_ext], der_order=der_order)
                        """
                        for lev in range(n_lev):
                            for c_ind in range(len(rgnd_ext.cellsHier[lev])):
                                if rgnd_ext.cellsHier[lev][c_ind][H_CI]>=0:
                                    print(self.inds_0_allp)
                                    print(self.inds_0_all_F)
                                    print('_polyHier_[lev][c_ind]:', _polyHier_ext_[rgnd_ext.VarIndex[par_ext]][lev][c_ind], sep='\n')
                                    print('_Der Hier_[lev][c_ind]:', polyHier_v_ext[lev][c_ind], sep='\n')
                                    pass
                        """

                    if self.verbose>0:
                        print('CoeffType.par_ext::', 'self.rhs_parts[i_eq][i_rp]:', self.rhs_parts[i_eq][i_rp])
                        print('der_order:', der_order, 'par_coeff_sym:', par_coeff_sym)
                                                    
                    par_coeff_sym_mark = None
                    for _par_ in self.pars_list:
                        if par_coeff_sym.has(_par_):
                            par_coeff_sym_mark = CoeffType.par_int
                            break
                    if par_coeff_sym_mark==None:
                        for _v_ in self.indepVars_list:
                            if par_coeff_sym.has(_v_):
                                par_coeff_sym_mark = CoeffType.nonconst_anal
                                break
                    if par_coeff_sym_mark==None:
                        par_coeff_sym_mark = CoeffType.const
                    
                    
                    par_coeff_sym_taylor = None
                    if par_coeff_sym_mark==CoeffType.par_int or par_coeff_sym_mark==CoeffType.nonconst_anal:
                        #par_coeff_sym_sub = par_coeff_sym
                        #if par_coeff_sym_mark==CoeffType.par_int:
                        #    for i_p in range(len(self.pars_list)):
                        #        par_coeff_sym_sub = par_coeff_sym_sub.subs(self.pars_list[i_p], self.pars_values[i_p])
                        #        assert par_coeff_sym_sub.has(self.pars_list[i_p])==False
                        #par_coeff_sym_taylor = self.SetTaylorCoeffAlllevels(par_coeff_sym_sub)
                        par_coeff_sym_taylor = self.SetTayCoAllLev_ParType(par_coeff_sym)
                    else:
                        par_coeff_sym_taylor = complex(par_coeff_sym)

                    for lev in range(n_lev):
                        ind_lev_st = N_Indexed_Cumul[lev]
                        cells_lev = self.cellsHier[lev]
                        
                        poly_I_rhs, inds_0 = self.GetPolyArrayRHSConst_masked_flat_Temp(i_eq, poly_orders, orders_max)
                        #print('poly_I_rhs:', poly_I_rhs, 'inds_0:', inds_0, sep='\n')
                        for i in range(len(cells_lev)):
                            if cells_lev[i][H_CI]>=0:

                                ##multiply polyHier_v_ext by par_coeff_sym
                                polyHier_v_ext_lev_i = polyHier_v_ext[lev][i]
                                coeff_vec = None
                                if self.rhs_parts[i_eq][i_rp][1]==CoeffType.par_int or self.rhs_parts[i_eq][i_rp][1]==CoeffType.nonconst_anal:
                                    assert False
                                    coeff_vec = np.zeros(n_F)
                                    for i in range(n_F):
                                        for j in range(len(P0P1_M[i])):
                                            ij_01, multip = P0P1_M[i][j]    ##inds (P0, P1), multiplicity
                                            coeff_vec[i] += par_coeff_sym_taylor[ij_01[0]]*polyHier_v_ext_lev_i[ij_01[1]]
                                else:
                                    coeff_vec = par_coeff_sym_taylor*polyHier_v_ext_lev_i
                                
                                if self.verbose>0:
                                    print('polyHier_v_ext_lev_i:', polyHier_v_ext_lev_i)

                                poly_rhs = [None]*len(poly_I_rhs)
                                for j in range(len(poly_I_rhs)):
                                    _ind_, _w_ = poly_I_rhs[j]
                                    poly_rhs[j] = _w_*coeff_vec[_ind_]
                                #print('inds_0:', inds_0, 'coeff_vec:', coeff_vec, 'poly_rhs:', poly_rhs, sep='\n')

                                c_ind_tot_st = (ind_lev_st + cells_lev[i][H_CI])*cell_n_tot \
                                    + i_eq*cell_n_tot_vi    ## i_eq=ind_v
                                rows_new = c_ind_tot_st + inds_0
                                data_new = poly_rhs
                                
                                data_rhs[rows_new] += data_new  ##TODO: check multiplicity
                                    

                elif self.rhs_parts[i_eq][i_rp][1]==CoeffType.par_int or self.rhs_parts[i_eq][i_rp][1]==CoeffType.nonconst_anal:
                    coeff_sym = self.rhs_parts[i_eq][i_rp][0]
                    #coeff_sym_sub = coeff_sym
                    #if self.rhs_parts[i_eq][i_rp][1]==CoeffType.par_int:
                    #    for i_p in range(len(self.pars_list)):
                    #        coeff_sym_sub = coeff_sym_sub.subs(self.pars_list[i_p], self.pars_values[i_p])
                    #        assert coeff_sym_sub.has(self.pars_list[i_p])==False
                    #coeff_sym_taylor = self.SetTaylorCoeffAlllevels(coeff_sym_sub)
                    coeff_sym_taylor = self.SetTayCoAllLev_ParType(coeff_sym)

                    for lev in range(n_lev):
                        ind_lev_st = N_Indexed_Cumul[lev]
                        cells_lev = self.cellsHier[lev]
                        
                        poly_I_rhs, inds_0 = self.GetPolyArrayRHSConst_masked_flat_Temp(i_eq, poly_orders, orders_max)
                        #print('poly_I_rhs:', poly_I_rhs, 'inds_0:', inds_0, sep='\n')
                        for i in range(len(cells_lev)):
                            if cells_lev[i][H_CI]>=0:

                                coeff_vec = coeff_sym_taylor[lev][i]
                                poly_rhs = [None]*len(poly_I_rhs)
                                for j in range(len(poly_I_rhs)):
                                    _ind_, _w_ = poly_I_rhs[j]
                                    poly_rhs[j] = _w_*coeff_vec[self.inds_0_all_F[_ind_][1]]
                                #print('poly_rhs:', poly_rhs)

                                c_ind_tot_st = (ind_lev_st + cells_lev[i][H_CI])*cell_n_tot \
                                    + i_eq*cell_n_tot_vi    ## i_eq=ind_v
                                rows_new = c_ind_tot_st + inds_0
                                data_new = poly_rhs
                                
                                data_rhs[rows_new] += data_new  ##TODO: check multiplicity
                    
                else:
                    raise NotImplementedError()
        
        ##BC_rhs NOT IN bcinds_cc
        for bc_ind, bc_cond in enumerate(self.BCs):
            #print('bc_cond:', bc_cond)
            if bc_ind in bcinds_cc:
                continue
            n_dir, expr, face = bc_cond['dir'], bc_cond['expr'], bc_cond['face']
            
            face_eqind, n_dir_eqind, v_ind_eqind, d_ord_eqind = self.bcic_to_bcic_map[bc_ind]['face'], \
                self.bcic_to_bcic_map[bc_ind]['dir'], self.bcic_to_bcic_map[bc_ind]['v_ind'], self.bcic_to_bcic_map[bc_ind]['d_ord']
                                
            #print('bc_ind:', bc_ind, 'n_dir:', n_dir, 'expr:', expr, 'face:', face)
            #print('n_dir_eqind:', n_dir_eqind, 'v_ind_eqind:', v_ind_eqind, 'der_ord_eqind:', der_ord_eqind)
        
            bc_rhs = self.BC_rhs[bc_ind]
            coeff_sym = bc_rhs
            coeff_is_constant = 0  ## 0:constant 1:has par  2:has x,y,z..
            for par in self.pars_list:
                if coeff_sym.has(par):
                    coeff_is_constant = 1
                    break
            if coeff_is_constant==0:
                for x_indep in self.indepVars_list:
                    if coeff_sym.has(x_indep):
                        coeff_is_constant = 2
                        break
                
            #print('coeff_sym:', coeff_sym)
            if coeff_is_constant==0:
                ##constant rhs
                coeff = complex(coeff_sym)
                #print('coeff:', coeff)
                for lev in range(n_lev):
                    ind_lev_st = N_Indexed_Cumul[lev]
                    cells_lev = cellsHier[lev]
                    bc_conns_lev_n_dir = self.boundaryCellConnections[lev]['n'][n_dir]
                    bc_conns_lev_p_dir = self.boundaryCellConnections[lev]['p'][n_dir]
                    
                    if face_eqind=='n':

                        #coeff_vec = np.zeros(len(rows_0), dtype=complex)
                        #coeff_vec[0] = coeff    ##TODO: check consistensy
                        #assert self.inds_0_all_F[0] == [0, tuple([0]*N)]
                        #data_0 = np.ones(len(rows_0))*coeff_vec
                        #print('rows_0:', rows_0, 'data_0:', data_0, sep='\n')
                        
                        #print('bc_conns_lev_n_dir:', bc_conns_lev_n_dir)
                        #print('bc_conns_lev_p_dir:', bc_conns_lev_p_dir)
                        for i in range(len(bc_conns_lev_n_dir)):
                            c_ind = bc_conns_lev_n_dir[i]
                            
                            #print('lev:', lev, 'c_ind:', c_ind)
                            if cells_lev[c_ind][H_CI]>=0:
                                mults_F_bc, eqinds_F_bc, inds_0_bc = self.cellsBC[v_ind_eqind][lev]['n'][n_dir_eqind][d_ord_eqind][c_ind]
                                
                                coeff_vec = np.zeros(len(mults_F_bc), dtype=complex)
                                coeff_vec[0] = coeff
                                assert self.inds_0_all_F[0] == [0, tuple([0]*N)]
                                data_0 = np.ones(len(mults_F_bc))*coeff_vec
                                
                                rows_new = eqinds_F_bc
                                data_new = data_0/np.array(mults_F_bc)
                                        
                                data_rhs[rows_new] += data_new
                                #print('mults_F_bc:', mults_F_bc, 'eqinds_F_bc:',  eqinds_F_bc, 'inds_0_bc:', inds_0_bc)
                    if face_eqind=='p':
                        #coeff_vec = np.zeros(len(rows_0), dtype=complex)
                        #coeff_vec[0] = coeff
                        #data_0 = np.ones(len(rows_0))*coeff_vec
                        #print('rows_0:', rows_0, 'data_0:', data_0, sep='\n')
                        
                        for i in range(len(bc_conns_lev_n_dir)):
                            c_ind = bc_conns_lev_p_dir[i]
                            c_nb_ind = bc_conns_lev_n_dir[i]
                            
                            if cells_lev[c_ind][H_CI]>=0: 
                                assert cells_lev[c_nb_ind][H_CI]>=0
                                mults_F_bc, eqinds_F_bc, inds_0_bc = self.cellsBC[v_ind_eqind][lev]['p'][n_dir_eqind][d_ord_eqind][c_ind]

                                coeff_vec = np.zeros(len(mults_F_bc), dtype=complex)
                                coeff_vec[0] = coeff
                                assert self.inds_0_all_F[0] == [0, tuple([0]*N)]
                                data_0 = np.ones(len(mults_F_bc))*coeff_vec
                                
                                rows_new = eqinds_F_bc
                                data_new = data_0/np.array(mults_F_bc)
                                        
                                data_rhs[rows_new] += data_new

            elif coeff_is_constant==1 or coeff_is_constant==2:
                ##variable rhs
                if self.verbose>0:
                    print('coeff_sym:', coeff_sym)
                coeff_sym_sub = coeff_sym
                if coeff_is_constant==1:
                    for i_p in range(len(self.pars_list)):
                        coeff_sym_sub = coeff_sym_sub.subs(self.pars_list[i_p], self.pars_values[i_p])
                        assert coeff_sym_sub.has(self.pars_list[i_p])==False
                if self.verbose>0:
                    print('coeff_sym_sub:', coeff_sym_sub)
                for lev in range(n_lev):
                    ind_lev_st = N_Indexed_Cumul[lev]
                    cells_lev = cellsHier[lev]
                    bc_conns_lev_n_dir = self.boundaryCellConnections[lev]['n'][n_dir]
                    bc_conns_lev_p_dir = self.boundaryCellConnections[lev]['p'][n_dir]
                    
                    if face_eqind=='n':
                    
                        bcTayco = self.SetTaylorCoeffAlllevels_face(coeff_sym_sub, face_eqind, n_dir, flat=True)
                        
                        #print('bc_conns_lev_n_dir:', bc_conns_lev_n_dir)
                        #print('bc_conns_lev_p_dir:', bc_conns_lev_p_dir)
                        for i in range(len(bc_conns_lev_n_dir)):
                            c_ind = bc_conns_lev_n_dir[i]
                            
                            #print('lev:', lev, 'c_ind:', c_ind)
                            if cells_lev[c_ind][H_CI]>=0:
                                mults_F_bc, eqinds_F_bc, inds_0_bc = self.cellsBC[v_ind_eqind][lev]['n'][n_dir_eqind][d_ord_eqind][c_ind]

                                r_0 = self.nodesPoints[cellsHier[lev][c_ind][H_CN][0]]
                                h = self.dx_levels[lev]
                                #[taylor, inds] = self.GetTaylorCoeffAll_face(f=coeff_sym_sub, r_0=r_0+self.x0, h=h, n_dir=n_dir, flat=True)                        
                                [taylor, _inds_] = bcTayco[lev][i]
                                #assert _inds_==inds_0_bc ##TODO: check relation of _inds_ and inds_0_bc
                                #print('_inds_:', _inds_, 'inds_0_bc:', inds_0_bc, sep='\n')
                                data_0 = taylor/np.array(mults_F_bc)
                                                                
                                rows_new = eqinds_F_bc
                                data_new = data_0
                                        
                                data_rhs[rows_new] += data_new

                                #print('c_ind:', c_ind, 'ind_lev_st:', ind_lev_st, 'c_ind_tot_st_row:', c_ind_tot_st_row)
                                #print('rows_0:', rows_0, 'data_0:', data_0, sep='\n')

                    if face_eqind=='p':
                        
                        bcTayco = self.SetTaylorCoeffAlllevels_face(coeff_sym_sub, face_eqind, n_dir, flat=True)

                        for i in range(len(bc_conns_lev_n_dir)):
                            c_ind = bc_conns_lev_p_dir[i]
                            c_nb_ind = bc_conns_lev_n_dir[i]
                            
                            if cells_lev[c_ind][H_CI]>=0: 
                                mults_F_bc, eqinds_F_bc, inds_0_bc = self.cellsBC[v_ind_eqind][lev]['p'][n_dir_eqind][d_ord_eqind][c_ind]

                                r_0 = self.nodesPoints[cellsHier[lev][c_ind][H_CN][0]]
                                h = self.dx_levels[lev]
                                #[taylor, inds] = self.GetTaylorCoeffAll_face(f=coeff_sym_sub, r_0=r_0+self.x0, h=h, n_dir=n_dir, flat=True)                        
                                [taylor, _inds_] = bcTayco[lev][i]
                                #assert _inds_==inds_0_bc
                                #print('_inds_:', _inds_, 'inds_0_bc:', inds_0_bc, sep='\n')
                                data_0 = taylor/np.array(mults_F_bc)

                                assert cells_lev[c_nb_ind][H_CI]>=0
                                
                                rows_new = eqinds_F_bc
                                data_new = data_0
                                        
                                data_rhs[rows_new] += data_new

                                #print('c_ind:', c_ind, 'c_nb_ind:', c_nb_ind, 'ind_lev_st:', ind_lev_st, 'c_ind_tot_st_row:', c_ind_tot_st_row)
                                #print('rows_0:', rows_0, 'data_0:', data_0, sep='\n')
            else:
                raise NotImplementedError()
        
        ##BC_rhs IN bcinds_cc
        row, col, data = [], [], []
        for _i_bccc, bc_ind in enumerate(bcinds_cc):
            #print('bc_cond:', bc_cond)
            bc_cond = self.BCs[bc_ind]
            n_dir, expr, face = bc_cond['dir'], bc_cond['expr'], bc_cond['face']
            
            face_eqind, n_dir_eqind, v_ind_eqind, d_ord_eqind = self.bcic_to_bcic_map[bc_ind]['face'], \
                self.bcic_to_bcic_map[bc_ind]['dir'], self.bcic_to_bcic_map[bc_ind]['v_ind'], self.bcic_to_bcic_map[bc_ind]['d_ord']

            assert face=='n'
            assert face_eqind=='n'
                                
            #print('bc_ind:', bc_ind, 'n_dir:', n_dir, 'expr:', expr, 'face:', face)
            #print('n_dir_eqind:', n_dir_eqind, 'v_ind_eqind:', v_ind_eqind, 'der_ord_eqind:', der_ord_eqind)
        

            v_ind, n_dir, d_ord, face = cc_equiv[_i_bccc]['v_ind'], cc_equiv[_i_bccc]['dir'], \
                    cc_equiv[_i_bccc]['der'], cc_equiv[_i_bccc]['face']
            assert face=='n'
            ## 'sl'
            for lev in range(n_lev):
                ind_lev_st = N_Indexed_Cumul[lev]
                cells_lev = cellsHier[lev]
                assert self.cc_dic['n']['sl'][lev][n_dir]!=None
                
                arr_inds_dir, arr_inds_dir_nb = self.cc_dic['n']['sl'][lev][n_dir][d_ord]
                arr_I_FL_f0_dir, inds_0_FL_f0_dir = arr_inds_dir
                arr_I_FL_f0_nb_dir, inds_0_FL_f0_nb_dir = arr_inds_dir_nb
                
                cc_multip, cc_eq_index = self.cc_eqindex_and_multiplicity[n_dir][d_ord]
                assert len(arr_I_FL_f0_dir)==len(cc_eq_index)
                
                ccnb_multip = self.cc_nbmul_dic['n']['sl'][n_dir]   ## neigbor multiplicity
                cc_multip = [m*ccnb_multip for m in cc_multip]
                
                n_el = sum([len(el) for el in arr_I_FL_f0_dir])
                n_el_nb = sum([len(el) for el in arr_I_FL_f0_nb_dir])
                rows_0, cols_0, data_0 = [None]*n_el, [None]*n_el, [None]*n_el
                rows_nb_0, cols_nb_0, data_nb_0 = [None]*n_el_nb, [None]*n_el_nb, [None]*n_el_nb
                ind_0 = 0
                for i in range(len(arr_I_FL_f0_dir)):
                    for col_ind in arr_I_FL_f0_dir[i]:
                        rows_0[ind_0] = cc_eq_index[i]
                        cols_0[ind_0] = col_ind
                        data_0[ind_0] = arr_I_FL_f0_dir[i][col_ind]/cc_multip[i]
                        ind_0 += 1
                assert ind_0==n_el
                ind_0 = 0
                for i in range(len(arr_I_FL_f0_nb_dir)):
                    for col_ind in arr_I_FL_f0_nb_dir[i]:
                        ##rows are set according to cell, not cell_nb
                        rows_nb_0[ind_0] = cc_eq_index[i] 
                        cols_nb_0[ind_0] = col_ind
                        data_nb_0[ind_0] = -arr_I_FL_f0_nb_dir[i][col_ind]/cc_multip[i]
                        ind_0 += 1
                assert ind_0==n_el_nb
                
                rows_0, cols_0, data_0 = np.array(rows_0), np.array(cols_0), np.array(data_0)
                rows_nb_0, cols_nb_0, data_nb_0 = np.array(rows_nb_0), np.array(cols_nb_0), np.array(data_nb_0)
                
                bc_conns_lev_n_dir = self.boundaryCellConnections[lev]['n'][n_dir]
                bc_conns_lev_p_dir = self.boundaryCellConnections[lev]['p'][n_dir]

                for i in range(len(bc_conns_lev_n_dir)):
                    c_ind = bc_conns_lev_n_dir[i]
                    c_nb_ind = bc_conns_lev_p_dir[i]
                    if not cells_lev[c_ind][H_CI]>=0:
                        assert not cells_lev[c_nb_ind][H_CI]>=0
                        continue
                    ## boundary cell
                    mults_F_bc, eqinds_F_bc, inds_0_bc = self.cellsBC[v_ind][lev]['n'][n_dir][d_ord][c_ind]
                    #mults_F_bc, inds_0_bc = self.cellsCC[v_ind][lev][n_dir][d_ord][c_ind]
                    #print('inds_0_bc:', inds_0_bc, 'cc_eq_index:', cc_eq_index, 'eqinds_F_bc:', eqinds_F_bc, sep='\n')
                    
                    cc_multip_bc = [m*ccnb_multip for m in mults_F_bc]
                    
                    
                    #data_0_bc = [None]*n_el
                    #ind_0 = 0
                    #for i in range(len(arr_I_FL_f0_dir)):
                    #    for col_ind in arr_I_FL_f0_dir[i]:
                    #        data_0_bc[ind_0] = data_0[ind_0]*cc_multip[i]/cc_multip_bc[i]
                    #        ind_0 += 1
                    #assert ind_0==n_el
                    
                    #data_0_bc = np.array(data_0_bc)
                    c_ind_tot_st = (ind_lev_st + cells_lev[c_ind][H_CI])*cell_n_tot \
                        + v_ind*cell_n_tot_vi
                    assert np.all(np.array(eqinds_F_bc)==c_ind_tot_st+np.array(cc_eq_index))
                    
                    #rows_new = c_ind_tot_st + rows_0
                    #cols_new = c_ind_tot_st + cols_0
                    #data_new = data_0_bc
                            
                    #row.extend(rows_new)
                    #col.extend(cols_new)
                    #data.extend(data_new)
                    

                    data_nb_0_bc = [None]*n_el_nb
                    ind_0 = 0
                    for i in range(len(arr_I_FL_f0_nb_dir)):
                        for col_ind in arr_I_FL_f0_nb_dir[i]:
                            data_nb_0_bc[ind_0] = data_nb_0[ind_0]*cc_multip[i]/cc_multip_bc[i]
                            ind_0 += 1
                    assert ind_0==n_el_nb

                    data_nb_0_bc = np.array(data_nb_0_bc)
                    c_nb_ind_tot_st = (ind_lev_st + cells_lev[c_nb_ind][H_CI])*cell_n_tot \
                        + v_ind*cell_n_tot_vi

                    rows_new = c_ind_tot_st + rows_nb_0
                    cols_new = c_nb_ind_tot_st + cols_nb_0
                    data_new = data_nb_0_bc

                    row.extend(rows_new)
                    col.extend(cols_new)
                    data.extend(data_new)

        row = np.array(row)
        col = np.array(col)
        data = np.array(data)
        A_ccTobc_csr = coo_matrix((data, (row, col)), shape=(n_total,n_total), dtype=complex).tocsr()
        
        data_rhs -= A_ccTobc_csr.dot(self.mat_res)
        #print('data_rhs:', data_rhs)
        if self.verbose:
            print('np.max(np.abs(data_rhs)): ', np.max(np.abs(data_rhs)))

        ##TODO: for variable coefficients self.mat_coo has to be updated
        from scipy.sparse import linalg
        
        x_res = None
        if self.matsolver=='UMF':
            x_res = linalg.spsolve(self.mat_coo.tocsc(), data_rhs)
        elif self.matsolver=='SLU':
            x_res = self.LU.solve(data_rhs)
        else:
            raise NotImplementedError()
                
        if setMat_res:
            self.mat_res = x_res
        return x_res


    def GetDerivativeTemplate_masked_flat(self, poly_order, der_order, orders_max, level):
        """ template for taking derivative of a flat poly stored in polyHier for example
        """
        N = self.N
        assert len(poly_order)==N and len(der_order)==N and len(orders_max)==N
        dx_lev = self.dx_levels[level]
        poly_order_p1 = tuple([i+1 for i in poly_order])

        #shape = self.inds_keep.shape
        shape_p = self.mask_allp.shape

        arr = np.copy(self.mask_allp).astype(float)    ## arr[i][j][k] --> x^i*y^j*z^k    
        inds = np.copy(self.inds_0_allp)
        inds_0 = np.copy(self.inds_0_allp)
        
        ##take derivative 
        for i in range(N):
            if der_order[i]>0:
                a = np.ones(shape_p[i])
                for j in range(der_order[i]):
                    a *= (np.arange(shape_p[i])-j)/dx_lev[i]
                #print('i:{} a:\n'.format(i), a)
                for j in range(shape_p[i]):
                    indx = [Ellipsis]*arr.ndim
                    indx[i] = j
                    arr[indx] *= a[j]
                    #print('der -- j:{}  arr:\n'.format(j), arr)
        ##roll
        for i in range(N):
            arr = np.roll(arr, -der_order[i], i)
            inds = np.roll(inds, -der_order[i], i)
            #print('rolling -- i:{}  arr:\n'.format(i), arr)

        ##mask zeros
        arr *= self.mask_allp
        inds = inds*self.mask_allp - (self.mask_allp==0)*1
        
        #print('arr:', arr, 'inds:', inds, 'inds_0:', inds_0, sep='\n')
        
        inds_F, inds_0_F, arr_F = [], [], []
        n_term_F = len(self.inds_0_all_F)
        for i in range(n_term_F):
            i_tup = self.inds_0_all_F[i][1]
            if arr[i_tup]!=0.0:
                inds_0_F.append(inds_0[i_tup])
                inds_F.append(inds[i_tup])
                arr_F.append(arr[i_tup])
                
        row = np.array(inds_0_F)
        col = np.array(inds_F)
        data = np.array(arr_F)
        A_csr = coo_matrix((data, (row, col)), shape=(n_term_F,n_term_F)).tocsr()
        
        return A_csr
        

    def SetDerivativeTemplatesForAllLevels(self, der_ords_max):
        N = self.N
        n_lev = len(self.cellsHier)
        n_vars = len(self.vars_list)
        
        der_template = [None]*n_vars
        
        for v_ind in range(n_vars):
            poly_order = self.poly_ords[v_ind]
            orders_max = self.orders_max[v_ind]     ##max poly orders to keep
            
            der_template[v_ind] = {}

            counter = np.zeros(N, dtype=int)
            while True:
                der_order = tuple(counter)
                der_template[v_ind][der_order] = [None]*n_lev
                for lev in range(n_lev):
                    A_csr = self.GetDerivativeTemplate_masked_flat(poly_order, der_order, orders_max, lev)
                    der_template[v_ind][der_order][lev] = A_csr
                if not self._increasePolyCounterIndex(counter, der_ords_max):
                    break
        self.der_template = der_template
        return der_template


    def TakeDerivativePolyHier(self, polyHier, v_ind, der_order):
        cellsHier = self.cellsHier
        n_lev = len(cellsHier)
        der_temp = self.der_template[v_ind][tuple(der_order)]
        der_polyHier_v = [None]*n_lev
        for lev in range(n_lev):
            cellsHier_lev = cellsHier[lev]
            n_cell_lev = len(cellsHier_lev)
            der_polyHier_v[lev] = [None]*n_cell_lev
            der_temp_lev = der_temp[lev]
            der_polyHier_v_lev = der_polyHier_v[lev]
            polyHier_v_lev = polyHier[v_ind][lev]
            for c_ind in range(n_cell_lev):
                if cellsHier_lev[c_ind][H_CI]>=0:
                    der_polyHier_v_lev[c_ind] = der_temp_lev.dot(polyHier_v_lev[c_ind])
             
        return der_polyHier_v

    
    def ResetPolyHier(self, resetMatRes=True):
        n_lev = len(self.cellsHier)
        n_vars = len(self.vars_list)
        n_term_F = len(self.inds_0_all_F)
        for v_ind in range(n_vars):
            for lev in range(n_lev):
                cellsHier_lev = self.cellsHier[lev]
                n_cell_lev = len(cellsHier_lev)
                polyHier_v_lev = self .polyHier[v_ind][lev]
                for c_ind in range(n_cell_lev):
                    if cellsHier_lev[c_ind][H_CI]>=0:
                        polyHier_v_lev[c_ind] = np.zeros(n_term_F)
                    else:
                        polyHier_v_lev[c_ind] = None
        
        if resetMatRes:
            self.PolyHierToMatRes(setSelfMatRes=True)
        return

    def SetupMatrixEQs___ver_0(self):
        ## ** Works for ICs but gives singular matrices for BCs **
        """ take the field continuity on the boundaries, as initial conditions for
        the next cell. i.e. the initial condition for each cell is satisfied through
        the continuity condition at x=0, y=0... and the continuity condition at
        x=1 or y=1.. serves as initial condition for the cell on right or above.
        In case of boundary conditions the conditions at x=x0 and x=x1 serve as 
        filling the requirements of the leftmost (or downmost) cells.. 
        
        choice 1: create the final sparce matrix directly
        
        choice 2: put each sparce row inside the corressponding polyHier hierarchical 
        list, and finally gather them in a matrix
            - the indices can be either local (cell index, poly order index) or 
            final (final index number)
            PROS: may be more efficient for multi-grid methods
        """
        ##we assume all variables have the same polynomial order
        ##TODO: possibility of using unequal orders to be analyzed
        N = self.N
        cellsHier = self.cellsHier
        n_lev = len(cellsHier)
        poly_orders = self.poly_ords[0]
        orders_max = self.orders_max[0]     ##max poly orders to keep
        der_orders_max = self.der_orders_max[0]
        
        cell_n_tot__vars = self.GetNumOfUnknownsInEachCellForEachVar()
        print('cell_n_tot__vars: ', cell_n_tot__vars)
        cell_n_tot = sum(cell_n_tot__vars)
        vars_list = self.vars_list
        n_vars = len(vars_list)
        cell_n_tot_vi = cell_n_tot__vars[0]
        for i in range(n_vars):
            assert cell_n_tot__vars[i]==cell_n_tot_vi
            
        ##cumulative number of indexed cells        
        assert self.N_indexed_acc!=None
        N_Indexed_Cumul = self.N_indexed_acc
        """
        N_Indexed_Cumul = [0]*(n_lev+1)
        for lev in range(n_lev):
            N_Indexed_Cumul[lev+1] = N_Indexed_Cumul[lev]+self.indexInfo[lev]['N']
        self.N_Indexed_Cumul = N_Indexed_Cumul
        """
        
        ##get the total number of unknowns
        n_total = self.GetNumOfUnknowns()
        
        ##TODO:calculate the total number of elements to preallocate arrays
        n_elem = None 

        ##get the sparse matrix components 
        row, col, data = [], [], []

        ##---
        self.SetEqIndicesToKeep(poly_orders, orders_max)

        
        ##set matrix components from interior interactions
        print('self.EQs_parts:', self.EQs_parts)
        for eq_ind, eq_arr in enumerate(self.EQs_parts):
            #eq_arr: a differential eq in the input set of differential equations
            print('eq_ind:', eq_ind, '  eq_arr:', eq_arr)
            for v_ind, eq_v in enumerate(eq_arr):
                #eq_v: the differential equation associated with the v-th variable
                print('v_ind:', v_ind, '  eq_v:', eq_v)
                for eq_part in eq_v:
                    #eq_part: an additive part of the differential equation
                    print('eq_part:', eq_part)
                    coeff_sym, der_ord, der_vars = eq_part
                    
                    display(Math(latex(eq_part)+' : ' + latex(coeff_sym) + \
                    ';' + latex(der_ord) + ';' + latex(der_vars)))
                    
                    der_orders = self.ChangeDerOrderAndDerVarsToStandardForm(der_ord, der_vars)
                    
                    poly_diff = [None]*n_lev
                    inds_diff = [None]*n_lev
                    inds_orig = [None]*n_lev
                    for lev in range(n_lev):
                        poly_diff[lev], inds_diff[lev], inds_orig[lev] = \
                        self.GetPolyArrayDiff_scaled_indices_masked_flat(eq_ind, poly_orders, der_orders, orders_max, lev)
                        #print('lev:', lev, 'poly_diff:', poly_diff[lev], 'inds_diff:',\
                        # inds_diff[lev], 'inds_orig:', inds_orig[lev], sep='\n')
                    
                    coeff_is_constant = True
                    for par in self.pars_list:
                        if coeff_sym.has(par):
                            coeff_is_constant = False
                            break
                    if coeff_is_constant:
                        ##constant coefficients
                        print('coeff_sym: ', coeff_sym)
                        coeff = complex(coeff_sym)
                        for lev in range(n_lev):
                            ind_lev_st = N_Indexed_Cumul[lev]
                            cells_lev = self.cellsHier[lev]
                            
                            for i in range(len(cells_lev)):
                                if cells_lev[i][H_CI]>=0:
                                    c_ind_tot_st_row = (ind_lev_st + cells_lev[i][H_CI])*cell_n_tot \
                                        + eq_ind*cell_n_tot_vi
                                    c_ind_tot_st_col = (ind_lev_st + cells_lev[i][H_CI])*cell_n_tot \
                                        + v_ind*cell_n_tot_vi
                                    rows_new = c_ind_tot_st_row + inds_orig[lev]
                                    cols_new = c_ind_tot_st_col + inds_diff[lev]
                                    data_new = coeff*poly_diff[lev]
                                    
                                    row.extend(rows_new)
                                    col.extend(cols_new)
                                    data.extend(data_new)
                                    
                            assert np.all(inds_orig[lev]<cell_n_tot_vi)
                            assert np.all(inds_diff[lev]<cell_n_tot_vi)
                    else:
                        ##variable coefficients
                        raise NotImplementedError()

        ##TODO: for simplicity assuming all variables are on equal footing in terms
        ## of the derivative orders, to be generalized (poly_orders may need to be
        ## different for each var otherwise)..
        
        ## continuity conditions across cells boundaries
        if not self.cells_nb_dic:
            self.GetNeigborCellTypesInTermsOfConnectedNodes()
        assert self.cells_nb_dic!=None
        if not self.cc_dic:
            self.GetCellContCondEqs_scaled_masked(poly_orders, polyorders_max=orders_max, der_ord_max=der_orders_max)
        assert self.cc_dic!=None
        if not self.cc_nbmul_dic:
            self.SetCCCellMultiplicity()
        assert self.cc_nbmul_dic!=None
        CCs = self.CCs
        for var in CCs:
            v_ind = self.VarIndex[var]    ## variable index (order)
            v_ccs = CCs[var] ## variable's continuity conditions    
            for v_cc in v_ccs:
                n_dir, d_ord, face = v_cc['dir'], v_cc['der'], v_cc['face']
                if face=='n':
                    ## 'sl'
                    for lev in range(n_lev):
                        ind_lev_st = N_Indexed_Cumul[lev]
                        cells_lev = cellsHier[lev]
                        cells_nb_n_sl_lev_dir = self.cells_nb_dic['n']['sl'][lev][n_dir]
                        if cells_nb_n_sl_lev_dir==None:
                            continue
                        if self.cc_dic['n']['sl'][lev][n_dir]==None:
                            continue
                        arr_inds_dir, arr_inds_dir_nb = self.cc_dic['n']['sl'][lev][n_dir][d_ord]
                        arr_I_FL_f0_dir, inds_0_FL_f0_dir = arr_inds_dir
                        arr_I_FL_f0_nb_dir, inds_0_FL_f0_nb_dir = arr_inds_dir_nb
                        
                        cc_multip, cc_eq_index = self.cc_eqindex_and_multiplicity[n_dir][d_ord]
                        assert len(arr_I_FL_f0_dir)==len(cc_eq_index)
                        
                        ccnb_multip = self.cc_nbmul_dic['n']['sl'][n_dir]   ## neigbor multiplicity
                        cc_multip = [m*ccnb_multip for m in cc_multip]
                        
                        n_el = sum([len(el) for el in arr_I_FL_f0_dir])
                        n_el_nb = sum([len(el) for el in arr_I_FL_f0_nb_dir])
                        rows_0, cols_0, data_0 = [None]*n_el, [None]*n_el, [None]*n_el
                        rows_nb_0, cols_nb_0, data_nb_0 = [None]*n_el_nb, [None]*n_el_nb, [None]*n_el_nb
                        ind_0 = 0
                        for i in range(len(arr_I_FL_f0_dir)):
                            for col_ind in arr_I_FL_f0_dir[i]:
                                rows_0[ind_0] = cc_eq_index[i]
                                cols_0[ind_0] = col_ind
                                data_0[ind_0] = arr_I_FL_f0_dir[i][col_ind]/cc_multip[i]
                                ind_0 += 1
                        assert ind_0==n_el
                        ind_0 = 0
                        for i in range(len(arr_I_FL_f0_nb_dir)):
                            for col_ind in arr_I_FL_f0_nb_dir[i]:
                                ##rows are set according to cell, not cell_nb
                                rows_nb_0[ind_0] = cc_eq_index[i] 
                                cols_nb_0[ind_0] = col_ind
                                data_nb_0[ind_0] = -arr_I_FL_f0_nb_dir[i][col_ind]/cc_multip[i]
                                ind_0 += 1
                        assert ind_0==n_el_nb
                        
                        rows_0, cols_0, data_0 = np.array(rows_0), np.array(cols_0), np.array(data_0)
                        rows_nb_0, cols_nb_0, data_nb_0 = np.array(rows_nb_0), np.array(cols_nb_0), np.array(data_nb_0)
                        
                        #print('rows_0:', rows_0)
                        assert np.all(rows_0<cell_n_tot)
                        
                        for i in range(len(cells_nb_n_sl_lev_dir)):
                            c_ind, c_nb_ind = cells_nb_n_sl_lev_dir[i]
                            assert cells_lev[c_ind][H_CC]==NEXIST and cells_lev[c_nb_ind][H_CC]==NEXIST
                            assert cells_lev[c_ind][H_CI]>=0 and cells_lev[c_nb_ind][H_CI]>=0
                            c_ind_tot_st = (ind_lev_st + cells_lev[c_ind][H_CI])*cell_n_tot \
                                + v_ind*cell_n_tot_vi
                            
                            rows_new = c_ind_tot_st + rows_0
                            cols_new = c_ind_tot_st + cols_0
                            data_new = data_0
                                    
                            row.extend(rows_new)
                            col.extend(cols_new)
                            data.extend(data_new)

                            c_nb_ind_tot_st = (ind_lev_st + cells_lev[c_nb_ind][H_CI])*cell_n_tot \
                                + v_ind*cell_n_tot_vi

                            rows_new = c_ind_tot_st + rows_nb_0
                            cols_new = c_nb_ind_tot_st + cols_nb_0
                            data_new = data_nb_0

                            row.extend(rows_new)
                            col.extend(cols_new)
                            data.extend(data_new)

                            assert np.all(rows_0<cell_n_tot_vi)
                            assert np.all(cols_0<cell_n_tot_vi)
                            assert np.all(rows_nb_0<cell_n_tot_vi)
                            assert np.all(cols_nb_0<cell_n_tot_vi)
                    
                    ## 'nl'
                    for lev in range(n_lev):
                        ind_lev_st = N_Indexed_Cumul[lev]
                        ind_levp1_st = None
                        if lev<n_lev-1:
                            ind_levp1_st = N_Indexed_Cumul[lev+1]
                        cells_lev = cellsHier[lev]
                        cells_levp1 = None
                        if lev<n_lev-1:
                            cells_levp1 = cellsHier[lev+1]
                        cells_nb_n_nl_lev_dir = self.cells_nb_dic['n']['nl'][lev][n_dir]
                        if cells_nb_n_nl_lev_dir==None:
                            continue
                        for n_conn in cells_nb_n_nl_lev_dir:
                            arr_inds_dir, arr_inds_dir_nb = self.cc_dic['n']['nl'][lev][n_dir][n_conn][d_ord]
                            arr_I_FL_f0_dir, inds_0_FL_f0_dir = arr_inds_dir
                            arr_I_FL_f0_nb_dir, inds_0_FL_f0_nb_dir = arr_inds_dir_nb
                            
                            cc_multip, cc_eq_index = self.cc_eqindex_and_multiplicity[n_dir][d_ord]
                            assert len(arr_I_FL_f0_dir)==len(cc_eq_index)
                            
                            ccnb_multip = self.cc_nbmul_dic['n']['nl'][n_dir]   ## neigbor multiplicity
                            cc_multip = [m*ccnb_multip for m in cc_multip]
                            
                            n_el = sum([len(el) for el in arr_I_FL_f0_dir])
                            n_el_nb = sum([len(el) for el in arr_I_FL_f0_nb_dir])
                            rows_0, cols_0, data_0 = [None]*n_el, [None]*n_el, [None]*n_el
                            rows_nb_0, cols_nb_0, data_nb_0 = [None]*n_el_nb, [None]*n_el_nb, [None]*n_el_nb
                            ind_0 = 0
                            for i in range(len(arr_I_FL_f0_dir)):
                                for col_ind in arr_I_FL_f0_dir[i]:
                                    rows_0[ind_0] = cc_eq_index[i]
                                    cols_0[ind_0] = col_ind
                                    data_0[ind_0] = arr_I_FL_f0_dir[i][col_ind]/cc_multip[i]
                                    ind_0 += 1
                            assert ind_0==n_el
                            ind_0 = 0
                            for i in range(len(arr_I_FL_f0_nb_dir)):
                                for col_ind in arr_I_FL_f0_nb_dir[i]:
                                    ##rows are set according to cell, not cell_nb
                                    rows_nb_0[ind_0] = cc_eq_index[i] 
                                    cols_nb_0[ind_0] = col_ind
                                    data_nb_0[ind_0] = -arr_I_FL_f0_nb_dir[i][col_ind]/cc_multip[i]
                                    ind_0 += 1
                            assert ind_0==n_el_nb
                        
                            rows_0, cols_0, data_0 = np.array(rows_0), np.array(cols_0), np.array(data_0)
                            rows_nb_0, cols_nb_0, data_nb_0 = np.array(rows_nb_0), np.array(cols_nb_0), np.array(data_nb_0)

                            for i in range(len(cells_nb_n_nl_lev_dir[n_conn])):
                                c_ind, c_nb_ind = cells_nb_n_nl_lev_dir[n_conn][i]
                                assert cells_lev[c_ind][H_CC]==NEXIST and cells_levp1[c_nb_ind][H_CC]==NEXIST
                                assert cells_lev[c_ind][H_CI]>=0 and cells_levp1[c_nb_ind][H_CI]>=0
                                c_ind_tot_st = (ind_lev_st + cells_lev[c_ind][H_CI])*cell_n_tot \
                                    + v_ind*cell_n_tot_vi
                                
                                rows_new = c_ind_tot_st + rows_0
                                cols_new = c_ind_tot_st + cols_0
                                data_new = data_0
                                        
                                row.extend(rows_new)
                                col.extend(cols_new)
                                data.extend(data_new)

                                c_nb_ind_tot_st = (ind_levp1_st + cells_levp1[c_nb_ind][H_CI])*cell_n_tot \
                                    + v_ind*cell_n_tot_vi

                                rows_new = c_ind_tot_st + rows_nb_0
                                cols_new = c_nb_ind_tot_st + cols_nb_0
                                data_new = data_nb_0

                                row.extend(rows_new)
                                col.extend(cols_new)
                                data.extend(data_new)
                    
                    ## 'pl'
                    for lev in range(n_lev):
                        ind_lev_st = N_Indexed_Cumul[lev]
                        cells_lev = cellsHier[lev]
                        cells_levm1 = None
                        if lev>0:
                            cells_levm1 = cellsHier[lev-1]
                        ind_levm1_st = None
                        if lev>0:
                            ind_levm1_st = N_Indexed_Cumul[lev-1]
                        cells_nb_n_pl_lev_dir = self.cells_nb_dic['n']['pl'][lev][n_dir]
                        if cells_nb_n_pl_lev_dir==None:
                            continue
                        for n_conn in cells_nb_n_pl_lev_dir:
                            arr_inds_dir, arr_inds_dir_nb = self.cc_dic['n']['pl'][lev][n_dir][n_conn][d_ord]
                            arr_I_FL_f0_dir, inds_0_FL_f0_dir = arr_inds_dir
                            arr_I_FL_f0_nb_dir, inds_0_FL_f0_nb_dir = arr_inds_dir_nb
                            
                            cc_multip, cc_eq_index = self.cc_eqindex_and_multiplicity[n_dir][d_ord]
                            assert len(arr_I_FL_f0_dir)==len(cc_eq_index)
                            
                            ccnb_multip = self.cc_nbmul_dic['n']['pl'][n_dir]   ## neigbor multiplicity
                            cc_multip = [m*ccnb_multip for m in cc_multip]
                            
                            n_el = sum([len(el) for el in arr_I_FL_f0_dir])
                            n_el_nb = sum([len(el) for el in arr_I_FL_f0_nb_dir])
                            rows_0, cols_0, data_0 = [None]*n_el, [None]*n_el, [None]*n_el
                            rows_nb_0, cols_nb_0, data_nb_0 = [None]*n_el_nb, [None]*n_el_nb, [None]*n_el_nb
                            ind_0 = 0
                            for i in range(len(arr_I_FL_f0_dir)):
                                for col_ind in arr_I_FL_f0_dir[i]:
                                    rows_0[ind_0] = cc_eq_index[i]
                                    cols_0[ind_0] = col_ind
                                    data_0[ind_0] = arr_I_FL_f0_dir[i][col_ind]/cc_multip[i]
                                    ind_0 += 1
                            assert ind_0==n_el
                            ind_0 = 0
                            for i in range(len(arr_I_FL_f0_nb_dir)):
                                for col_ind in arr_I_FL_f0_nb_dir[i]:
                                    ##rows are set according to cell, not cell_nb
                                    rows_nb_0[ind_0] = cc_eq_index[i] 
                                    cols_nb_0[ind_0] = col_ind
                                    data_nb_0[ind_0] = -arr_I_FL_f0_nb_dir[i][col_ind]/cc_multip[i]
                                    ind_0 += 1
                            assert ind_0==n_el_nb
                        
                            rows_0, cols_0, data_0 = np.array(rows_0), np.array(cols_0), np.array(data_0)
                            rows_nb_0, cols_nb_0, data_nb_0 = np.array(rows_nb_0), np.array(cols_nb_0), np.array(data_nb_0)

                            for i in range(len(cells_nb_n_pl_lev_dir[n_conn])):
                                c_ind, c_nb_ind = cells_nb_n_pl_lev_dir[n_conn][i]
                                assert cells_lev[c_ind][H_CI]>=0 and cells_levm1[c_nb_ind][H_CI]>=0
                                c_ind_tot_st = (ind_lev_st + cells_lev[c_ind][H_CI])*cell_n_tot \
                                    + v_ind*cell_n_tot_vi
                                
                                rows_new = c_ind_tot_st + rows_0
                                cols_new = c_ind_tot_st + cols_0
                                data_new = data_0
                                        
                                row.extend(rows_new)
                                col.extend(cols_new)
                                data.extend(data_new)

                                c_nb_ind_tot_st = (ind_levm1_st + cells_levm1[c_nb_ind][H_CI])*cell_n_tot \
                                    + v_ind*cell_n_tot_vi

                                rows_new = c_ind_tot_st + rows_nb_0
                                cols_new = c_nb_ind_tot_st + cols_nb_0
                                data_new = data_nb_0

                                row.extend(rows_new)
                                col.extend(cols_new)
                                data.extend(data_new)
                    
                else:
                    raise NotImplementedError()
                        
        ##boundary conditions
        BC_parts, BC_rhs = self.DisintegrateBoundaryConds()
        if not self.bc_dic:
            self.GetCellBoundCondEqs_scaled_masked(poly_orders, polyorders_max=orders_max, der_ord_max=der_orders_max)
        assert self.bc_dic!=None
        assert self.bc_eqindex_and_multiplicity!=None
        if not self.bctoic_map:
            self.MapBCsToICs(self.der_orders_max)
        assert self.bctoic_map!=None
        bctoic_map = self.bctoic_map
        
        for bc_ind, b_cond in enumerate(self.BCs):
            n_dir, expr, face = b_cond['dir'], b_cond['expr'], b_cond['face']
            
            n_dir_eqind, v_ind_eqind, der_ord_eqind = bctoic_map[bc_ind]['dir'],\
                bctoic_map[bc_ind]['v_ind'], bctoic_map[bc_ind]['d_ord']
            
            bc_eq_arr = BC_parts[bc_ind]
            for v_ind in bc_eq_arr:
                eq_v = bc_eq_arr[v_ind]
                #eq_v: the differential equation associated with the v-th variable
                for eq_part in eq_v:
                    #eq_part: an additive part of the differential equation
                    coeff_sym, der_ord, der_vars = eq_part
                    
                    der_orders = self.ChangeDerOrderAndDerVarsToStandardForm(der_ord, der_vars)
                    d_ord = der_orders[n_dir]
                    for i in range(N):
                        if i != n_dir:
                            assert der_orders[i]==0

                    coeff_is_constant = True
                    for par in self.pars_list:
                        if coeff_sym.has(par):
                            coeff_is_constant = False
                            break
                    if coeff_is_constant:
                        ##constant coefficients
                        coeff = complex(coeff_sym)
                        for lev in range(n_lev):
                            ind_lev_st = N_Indexed_Cumul[lev]
                            cells_lev = cellsHier[lev]
                            bc_conns_lev_n_dir = self.boundaryCellConnections[lev]['n'][n_dir]
                            bc_conns_lev_p_dir = self.boundaryCellConnections[lev]['p'][n_dir]
                            
                            if face=='n':
                                arr_I_FL_f0_dir, inds_0_FL_f0_dir = self.bc_dic['n'][n_dir][lev][d_ord]
                                
                                bc_multip, bc_eq_index = self.bc_eqindex_and_multiplicity[n_dir_eqind][der_ord_eqind]
                                assert len(bc_eq_index)==len(inds_0_FL_f0_dir)
                                
                                n_el = sum([len(el) for el in arr_I_FL_f0_dir])
                                rows_0, cols_0, data_0 = [None]*n_el, [None]*n_el, [None]*n_el
                                ind_0 = 0
                                for i in range(len(arr_I_FL_f0_dir)):
                                    for col_ind in arr_I_FL_f0_dir[i]:
                                        rows_0[ind_0] = bc_eq_index[i]
                                        cols_0[ind_0] = col_ind
                                        data_0[ind_0] = arr_I_FL_f0_dir[i][col_ind]/bc_multip[i]
                                        ind_0 += 1
                                assert ind_0==n_el

                                rows_0, cols_0, data_0 = np.array(rows_0), np.array(cols_0), np.array(data_0)
                                print('lev:', lev, 'n_dir:', n_dir, 'v_ind:', v_ind, 'face:', face, 'der_orders:', der_orders)
                                print('rows_0:', rows_0, 'cols_0:', cols_0, 'data_0:', data_0, '-'*20, sep='\n')

                                for i in range(len(bc_conns_lev_n_dir)):
                                    c_ind = bc_conns_lev_n_dir[i]
                                    
                                    if cells_lev[c_ind][H_CI]>=0:
                                        c_ind_tot_st_row = (ind_lev_st + cells_lev[c_ind][H_CI])*cell_n_tot \
                                            + v_ind_eqind*cell_n_tot_vi
                                        c_ind_tot_st_col = (ind_lev_st + cells_lev[c_ind][H_CI])*cell_n_tot \
                                            + v_ind*cell_n_tot_vi
                                        
                                        rows_new = c_ind_tot_st_row + rows_0
                                        cols_new = c_ind_tot_st_col + cols_0
                                        data_new = data_0
                                                
                                        row.extend(rows_new)
                                        col.extend(cols_new)
                                        data.extend(data_new)
                            else:
                                assert face=='p'
                                arr_I_FL_f0_dir, inds_0_FL_f0_dir = self.bc_dic['p'][n_dir][lev][d_ord]
                                
                                bc_multip, bc_eq_index = self.bc_eqindex_and_multiplicity[n_dir_eqind][der_ord_eqind]
                                assert len(bc_eq_index)==len(inds_0_FL_f0_dir)

                                n_el = sum([len(el) for el in arr_I_FL_f0_dir])
                                rows_0, cols_0, data_0 = [None]*n_el, [None]*n_el, [None]*n_el
                                ind_0 = 0
                                for i in range(len(arr_I_FL_f0_dir)):
                                    for col_ind in arr_I_FL_f0_dir[i]:
                                        rows_0[ind_0] = bc_eq_index[i]
                                        cols_0[ind_0] = col_ind
                                        data_0[ind_0] = arr_I_FL_f0_dir[i][col_ind]/bc_multip[i]
                                        ind_0 += 1
                                assert ind_0==n_el

                                rows_0, cols_0, data_0 = np.array(rows_0), np.array(cols_0), np.array(data_0)

                                for i in range(len(bc_conns_lev_p_dir)):
                                    c_ind = bc_conns_lev_p_dir[i]
                                    c_nb_ind = bc_conns_lev_n_dir[i]
                                    
                                    if cells_lev[c_ind][H_CI]>=0:
                                        assert cells_lev[c_nb_ind][H_CI]>=0
                                        c_ind_tot_st_row = (ind_lev_st + cells_lev[c_nb_ind][H_CI])*cell_n_tot \
                                            + v_ind_eqind*cell_n_tot_vi
                                        c_ind_tot_st_col = (ind_lev_st + cells_lev[c_ind][H_CI])*cell_n_tot \
                                            + v_ind*cell_n_tot_vi
                                        
                                        rows_new = c_ind_tot_st_row + rows_0
                                        cols_new = c_ind_tot_st_col + cols_0
                                        data_new = data_0
                                                
                                        row.extend(rows_new)
                                        col.extend(cols_new)
                                        data.extend(data_new)

                                    
                    else:
                        raise NotImplementedError()

        ##TODO: for 'p' boundary conditions specify redundencies, add redundent 
        ## equations, and fill the remaining rows, with the remaining boundary equations
        

        ##test
        c_startend_inds = [0]
        for lev in range(n_lev):
            c_startend_inds.append(c_startend_inds[lev]+self.indexInfo[lev]['N'])
        ##oor: out of range indices (rows)
        rows_oor = [] 
        for i in range(len(row)):
            if row[i]>=n_total:
                rows_oor.append(row[i])
        print('rows_oor[{}]:'.format(len(rows_oor)), rows_oor)
        cells_oor = [None]*len(rows_oor)
        for i, r_ind in enumerate(rows_oor):
            c_lev = None
            for lev in range(n_lev):
                if r_ind>=c_startend_inds[lev]*cell_n_tot and r_ind<c_startend_inds[lev+1]*cell_n_tot:
                    c_lev = lev
                    r_ind_resi = r_ind - c_startend_inds[lev]*cell_n_tot
                    c_ind = int(r_ind_resi/cell_n_tot)
                    cells_oor[i] = (c_lev, c_ind)
                    break
        cells_oor = list(set(cells_oor))
        print('cells_oor[{}]:'.format(len(cells_oor)), cells_oor)
        ##oor: out of range indices (cols)
        rows_oor = [] 
        for i in range(len(row)):
            if col[i]>=n_total:
                rows_oor.append(row[i])
        print('cols_oor[{}]:'.format(len(rows_oor)), rows_oor)
        cells_oor = [None]*len(rows_oor)
        for i, r_ind in enumerate(rows_oor):
            c_lev = None
            for lev in range(n_lev):
                if r_ind>=c_startend_inds[lev]*cell_n_tot and r_ind<c_startend_inds[lev+1]*cell_n_tot:
                    c_lev = lev
                    r_ind_resi = r_ind - c_startend_inds[lev]*cell_n_tot
                    c_ind = int(r_ind_resi/cell_n_tot)
                    cells_oor[i] = (c_lev, c_ind)
                    break
        cells_oor = list(set(cells_oor))
        print('cells_oor[{}]:'.format(len(cells_oor)), cells_oor)
        ## unset rows
        print('n_total:', n_total)     
        rows_marked = [0]*n_total
        for i in range(len(row)):
            rows_marked[row[i]] = 1
        rows_unset = [i for i in range(n_total) if rows_marked[i]==0]
        print('rows_unset[{}]:'.format(len(rows_unset)), rows_unset)
        cells_unset = [None]*len(rows_unset)
        for i, r_ind in enumerate(rows_unset):
            c_lev = None
            for lev in range(n_lev):
                if r_ind>=c_startend_inds[lev]*cell_n_tot and r_ind<c_startend_inds[lev+1]*cell_n_tot:
                    c_lev = lev
                    r_ind_resi = r_ind - c_startend_inds[lev]*cell_n_tot
                    c_ind = int(r_ind_resi/cell_n_tot)
                    cells_unset[i] = (c_lev, c_ind)
                    break
        cells_unset = list(set(cells_unset))
        print('cells_unset[{}]:'.format(len(cells_unset)), cells_unset)
        ## unset columns
        cols_marked = [0]*n_total
        for i in range(len(col)):
            cols_marked[col[i]] = 1
        cols_unset = [i for i in range(n_total) if cols_marked[i]==0]
        print('cols_unset[{}]:'.format(len(cols_unset)), cols_unset)
        cells_unset = [None]*len(cols_unset)
        for i, r_ind in enumerate(cols_unset):
            c_lev = None
            for lev in range(n_lev):
                if r_ind>=c_startend_inds[lev]*cell_n_tot and r_ind<c_startend_inds[lev+1]*cell_n_tot:
                    c_lev = lev
                    r_ind_resi = r_ind - c_startend_inds[lev]*cell_n_tot
                    c_ind = int(r_ind_resi/cell_n_tot)
                    cells_unset[i] = (c_lev, c_ind)
                    break
        cells_unset = list(set(cells_unset))
        print('cells_unset[{}]:'.format(len(cells_unset)), cells_unset)
        ## rows with zero entries   
        rows_zero = [] 
        for i in range(len(row)):
            if data[i]==0:
                rows_zero.append(row[i])
        print('rows_zero[{}]:'.format(len(rows_zero)), rows_zero)
        cells_zero = [None]*len(rows_zero)
        for i, r_ind in enumerate(rows_zero):
            c_lev = None
            for lev in range(n_lev):
                if r_ind>=c_startend_inds[lev]*cell_n_tot and r_ind<c_startend_inds[lev+1]*cell_n_tot:
                    c_lev = lev
                    r_ind_resi = r_ind - c_startend_inds[lev]*cell_n_tot
                    c_ind = int(r_ind_resi/cell_n_tot)
                    cells_zero[i] = (c_lev, c_ind)
                    break
        cells_zero = list(set(cells_zero))
        print('cells_zero[{}]:'.format(len(cells_zero)), cells_zero)

        ##RHS
        data_rhs = np.zeros(n_total, dtype=complex)
        for i_eq in range(len(self.rhs)):
            coeff_sym = self.rhs[i_eq]
            coeff_is_constant = 0  ## 0:constant 1:has par  2:has x,y,z..
            for par in self.pars_list:
                if coeff_sym.has(par):
                    coeff_is_constant = 1
                    break
            ##constant coefficients
            if coeff_is_constant==0:
                coeff = complex(coeff_sym)
                for lev in range(n_lev):
                    ind_lev_st = N_Indexed_Cumul[lev]
                    cells_lev = self.cellsHier[lev]
                    
                    poly_rhs, inds_0 = self.GetPolyArrayRHSConst_masked_flat(i_eq, poly_orders, orders_max, coeff)
                    for i in range(len(cells_lev)):
                        if cells_lev[i][H_CI]>=0:
                            c_ind_tot_st = (ind_lev_st + cells_lev[i][H_CI])*cell_n_tot \
                                + i_eq*cell_n_tot_vi    ## i_eq=ind_v
                            rows_new = c_ind_tot_st + inds_0
                            data_new = poly_rhs
                            
                            data_rhs[rows_new] += data_new  ##TODO: check multiplicity
            else:
                raise NotImplementedError()
            
        ##BC RHS (boundary condition RHS)
        #print('self.BCs: ', self.BCs)
        for bc_ind, bc_cond in enumerate(self.BCs):
            #print('bc_cond:', bc_cond)
            n_dir, expr, face = bc_cond['dir'], bc_cond['expr'], bc_cond['face']
            
            n_dir_eqind, v_ind_eqind, der_ord_eqind = bctoic_map[bc_ind]['dir'],\
                bctoic_map[bc_ind]['v_ind'], bctoic_map[bc_ind]['d_ord']
                
            #print('bc_ind:', bc_ind, 'n_dir:', n_dir, 'expr:', expr, 'face:', face)
            #print('n_dir_eqind:', n_dir_eqind, 'v_ind_eqind:', v_ind_eqind, 'der_ord_eqind:', der_ord_eqind)
        
            bc_rhs = self.BC_rhs[bc_ind]
            coeff_sym = bc_rhs
            coeff_is_constant = 0  ## 0:constant 1:has par  2:has x,y,z..
            for par in self.pars_list:
                if coeff_sym.has(par):
                    coeff_is_constant = 1
                    break
            if coeff_is_constant==0:
                for x_indep in self.indepVars_list:
                    if coeff_sym.has(x_indep):
                        coeff_is_constant = 2
                        break
                
            #print('coeff_sym:', coeff_sym)
            if coeff_is_constant==0:
                ##constant rhs
                coeff = complex(coeff_sym)
                #print('coeff:', coeff)
                for lev in range(n_lev):
                    ind_lev_st = N_Indexed_Cumul[lev]
                    cells_lev = cellsHier[lev]
                    bc_conns_lev_n_dir = self.boundaryCellConnections[lev]['n'][n_dir]
                    bc_conns_lev_p_dir = self.boundaryCellConnections[lev]['p'][n_dir]
                    
                    if face=='n':
                        bc_multip, bc_eq_index = self.bc_eqindex_and_multiplicity[n_dir_eqind][der_ord_eqind]
                        assert len(bc_eq_index)==len(inds_0_FL_f0_dir)
                        
                        
                        rows_0 = np.array(bc_eq_index)
                        coeff_vec = np.zeros(len(rows_0), dtype=complex)
                        coeff_vec[0] = coeff    ##TODO: check consistensy
                        assert self.inds_0_all_F[0] == [0, tuple([0]*N)]
                        data_0 = np.ones(len(rows_0))*coeff_vec/np.array(bc_multip)
                        #print('rows_0:', rows_0, 'data_0:', data_0, sep='\n')
                        
                        #print('bc_conns_lev_n_dir:', bc_conns_lev_n_dir)
                        #print('bc_conns_lev_p_dir:', bc_conns_lev_p_dir)
                        for i in range(len(bc_conns_lev_n_dir)):
                            c_ind = bc_conns_lev_n_dir[i]
                            
                            #print('lev:', lev, 'c_ind:', c_ind)
                            if cells_lev[c_ind][H_CI]>=0:
                                c_ind_tot_st_row = (ind_lev_st + cells_lev[c_ind][H_CI])*cell_n_tot \
                                    + v_ind_eqind*cell_n_tot_vi
                                
                                rows_new = c_ind_tot_st_row + rows_0
                                data_new = data_0
                                        
                                data_rhs[rows_new] += data_new
                    if face=='p':
                        bc_multip, bc_eq_index = self.bc_eqindex_and_multiplicity[n_dir_eqind][der_ord_eqind]
                        assert len(bc_eq_index)==len(inds_0_FL_f0_dir)
                        
                        
                        rows_0 = np.array(bc_eq_index)
                        coeff_vec = np.zeros(len(rows_0), dtype=complex)
                        coeff_vec[0] = coeff
                        data_0 = np.ones(len(rows_0))*coeff_vec/np.array(bc_multip)
                        #print('rows_0:', rows_0, 'data_0:', data_0, sep='\n')
                        
                        for i in range(len(bc_conns_lev_n_dir)):
                            c_ind = bc_conns_lev_p_dir[i]
                            c_nb_ind = bc_conns_lev_n_dir[i]
                            
                            if cells_lev[c_ind][H_CI]>=0: 
                                assert cells_lev[c_nb_ind][H_CI]>=0
                                c_ind_tot_st_row = (ind_lev_st + cells_lev[c_nb_ind][H_CI])*cell_n_tot \
                                    + v_ind_eqind*cell_n_tot_vi
                                
                                rows_new = c_ind_tot_st_row + rows_0
                                data_new = data_0
                                        
                                data_rhs[rows_new] += data_new

            elif coeff_is_constant==2:
                ##variable rhs
                print('coeff_sym:', coeff_sym)
                for lev in range(n_lev):
                    ind_lev_st = N_Indexed_Cumul[lev]
                    cells_lev = cellsHier[lev]
                    bc_conns_lev_n_dir = self.boundaryCellConnections[lev]['n'][n_dir]
                    bc_conns_lev_p_dir = self.boundaryCellConnections[lev]['p'][n_dir]
                    
                    if face=='n':
                        bc_multip, bc_eq_index = self.bc_eqindex_and_multiplicity[n_dir_eqind][der_ord_eqind]
                        assert len(bc_eq_index)==len(inds_0_FL_f0_dir)
                        
                        rows_0 = np.array(bc_eq_index)
                        
                        #print('bc_conns_lev_n_dir:', bc_conns_lev_n_dir)
                        #print('bc_conns_lev_p_dir:', bc_conns_lev_p_dir)
                        for i in range(len(bc_conns_lev_n_dir)):
                            c_ind = bc_conns_lev_n_dir[i]
                            
                            #print('lev:', lev, 'c_ind:', c_ind)
                            if cells_lev[c_ind][H_CI]>=0:
                                r_0 = self.nodesPoints[cellsHier[lev][c_ind][H_CN][0]]
                                h = self.dx_levels[lev]
                                [taylor, inds] = self.GetTaylorCoeffAll_face(f=coeff_sym, r_0=r_0+self.x0, h=h, n_dir=n_dir, flat=True)                        
                                data_0 = taylor/np.array(bc_multip)
                                
                                c_ind_tot_st_row = (ind_lev_st + cells_lev[c_ind][H_CI])*cell_n_tot \
                                    + v_ind_eqind*cell_n_tot_vi
                                
                                rows_new = c_ind_tot_st_row + rows_0
                                data_new = data_0
                                        
                                data_rhs[rows_new] += data_new

                                print('c_ind:', c_ind, 'ind_lev_st:', ind_lev_st, 'c_ind_tot_st_row:', c_ind_tot_st_row)
                                print('rows_0:', rows_0, 'data_0:', data_0, sep='\n')

                    if face=='p':
                        bc_multip, bc_eq_index = self.bc_eqindex_and_multiplicity[n_dir_eqind][der_ord_eqind]
                        assert len(bc_eq_index)==len(inds_0_FL_f0_dir)

                        rows_0 = np.array(bc_eq_index)
                        
                        for i in range(len(bc_conns_lev_n_dir)):
                            c_ind = bc_conns_lev_p_dir[i]
                            c_nb_ind = bc_conns_lev_n_dir[i]
                            
                            if cells_lev[c_ind][H_CI]>=0: 
                                r_0 = self.nodesPoints[cellsHier[lev][c_ind][H_CN][0]]
                                h = self.dx_levels[lev]
                                [taylor, inds] = self.GetTaylorCoeffAll_face(f=coeff_sym, r_0=r_0+self.x0, h=h, n_dir=n_dir, flat=True)                        
                                data_0 = taylor/np.array(bc_multip)

                                assert cells_lev[c_nb_ind][H_CI]>=0
                                c_ind_tot_st_row = (ind_lev_st + cells_lev[c_nb_ind][H_CI])*cell_n_tot \
                                    + v_ind_eqind*cell_n_tot_vi
                                
                                rows_new = c_ind_tot_st_row + rows_0
                                data_new = data_0
                                        
                                data_rhs[rows_new] += data_new

                                print('c_ind:', c_ind, 'c_nb_ind:', c_nb_ind, 'ind_lev_st:', ind_lev_st, 'c_ind_tot_st_row:', c_ind_tot_st_row)
                                print('rows_0:', rows_0, 'data_0:', data_0, sep='\n')
            else:
                raise NotImplementedError()

        ##construct the sparse matrix
        row = np.array(row)
        col = np.array(col)
        data = np.array(data)
        A_coo = coo_matrix((data, (row, col)), shape=(n_total,n_total), dtype=complex)
        
        ##test
        print('min: ', min(abs(data)), '    max: ', max(abs(data)))
        #print('data_rhs', data_rhs)
        
        #for i in range(n_total):
        #    print(i, A_coo.getrow(i))
        #    print('-'*20)
        
        """
        rows_marked_same = [None]*n_total
        ind_same = 0
        for i in range(n_total):
            if rows_marked_same[i]!=None:
                continue
            row_i = A_coo.getrow(i).tocoo()
            is_same = False
            for j in range(i, n_total):
                row_j = A_coo.getrow(j).tocoo()
                if row_i.nnz!=row_j.nnz:
                    continue
                if np.all((row_i.col==row_j.col)) and np.all((row_i.data==row_j.data)):
                    is_same==True
                    rows_marked_same[i] = ind_same
                    rows_marked_same[j] = ind_same
            if is_same:
                ind_same += 1
        rows_same = [None]*ind_same
        for i in range(ind_same):
            rows_same[i] = []
            for j in range(n_total):
                if rows_marked_same[j]==ind_same:
                    rows_same[i].append(j)
        print('ind_same:', ind_same)
        for i in range(ind_same):
            print(rows_same[i])
            print(A_coo.getrow(rows_same[i][0]))
            print('-'*20)
        """
        
        
        ##solve the matrix
        print('solving matrix..')
        from scipy.sparse import linalg

        #rank = np.linalg.matrix_rank(A_coo)
        #print('Matrix rank = ', rank)
        
        x = linalg.spsolve(A_coo.tocsc(), data_rhs)
        
        self.mat_coo = A_coo
        self.mat_rhs = data_rhs
        self.mat_res = x
        
        return x


    def GetAllCellEquations(self, lev, c_ind):
        ##TODO: probably should be modified for boundary cells
        N = self.N
        cellsHier = self.cellsHier
        n_lev = len(cellsHier)
        
        cell_n_tot__vars = self.GetNumOfUnknownsInEachCellForEachVar()
        cell_n_tot = sum(cell_n_tot__vars)
        vars_list = self.vars_list
        n_vars = len(vars_list)
        cell_n_tot_vi = cell_n_tot__vars[0]
        
        N_Indexed_Cumul = self.N_Indexed_Cumul
        if cellsHier[lev][c_ind][H_CI]>=0:
            ind_lev_st = N_Indexed_Cumul[lev]
            v_ind = 0
            c_ind_tot_st_row = (ind_lev_st + cellsHier[lev][c_ind][H_CI])*cell_n_tot \
                + v_ind*cell_n_tot_vi
            
            A_coo = self.mat_coo
            rhs = self.mat_rhs
            row, col, data = A_coo.row, A_coo.col, A_coo.data
            row_c, col_c, data_c = [], [], []
            for i in range(len(row)):
                if c_ind_tot_st_row<=row[i]<c_ind_tot_st_row+cell_n_tot:
                     row_c.append(row[i])
                     col_c.append(col[i])
                     data_c.append(data[i])
            rhs_c = []
            for i in range(len(rhs)):
                if c_ind_tot_st_row<=i<c_ind_tot_st_row+cell_n_tot:
                     rhs_c.append(rhs[i])
                
            return [row_c, col_c, data_c, rhs_c]
        else:
            print("cell is not indexed!")
            return None 


    def SetPolyBasisCoeffs(self, x_res):
        if self.N_indexed_acc==None:
            self.GetAccumukativeNumOfIndexedCells()
        if self.verbose>0:
            print('self.N_indexed_acc:', self.N_indexed_acc)
    
        n_vars = len(self.vars_list)
        cellsHier = self.cellsHier
        n_lev = len(cellsHier)
        n_pind_cell = sum(self.n_unkn_v)
        
        n_pind_v_acc = [0]*(n_vars+1)
        for v in range(n_vars):
            n_pind_v_acc[v+1] = n_pind_v_acc[v] + self.n_unkn_v[v]
        if self.verbose>0:
            print('n_pind_v_acc:', n_pind_v_acc)
        
        if self.indexing=='children':
            for v in range(n_vars):
                n_pind_v = self.n_unkn_v[v]
                ind_st = n_pind_v_acc[v]
                for lev in range(n_lev):
                    n_cell_lev = len(cellsHier[lev])
                    cellsHier_lev = cellsHier[lev]
                    ## set up self.polyHier
                    self.polyHier[v][lev] = [None]*n_cell_lev
                    polyHier_v_lev = self.polyHier[v][lev]
                    for i in range(n_cell_lev):
                        if cellsHier_lev[i][H_CC]==NEXIST:
                            polyHier_v_lev[i] = x_res[ind_st:ind_st+n_pind_v]
                            ind_st += n_pind_cell
        #print('self.polyHier:', self.polyHier)
        return


    def PolyHierToMatRes(self, polyHier=None, setSelfMatRes=False):
        if self.N_indexed_acc==None:
            self.GetAccumukativeNumOfIndexedCells()
    
        n_vars = len(self.vars_list)
        cellsHier = self.cellsHier
        n_lev = len(cellsHier)
        n_pind_cell = sum(self.n_unkn_v)
        
        n_total = self.GetNumOfUnknowns()
        x_res = np.zeros(n_total)
        
        n_pind_v_acc = [0]*(n_vars+1)
        for v in range(n_vars):
            n_pind_v_acc[v+1] = n_pind_v_acc[v] + self.n_unkn_v[v]
        if self.verbose>0:
            print('n_pind_v_acc:', n_pind_v_acc)
        
        if polyHier==None:
            polyHier = self.polyHier
            
        if self.indexing=='children':
            for v in range(n_vars):
                n_pind_v = self.n_unkn_v[v]
                ind_st = n_pind_v_acc[v]
                for lev in range(n_lev):
                    n_cell_lev = len(cellsHier[lev])
                    cellsHier_lev = cellsHier[lev]
                    polyHier_v_lev = polyHier[v][lev]
                    for i in range(n_cell_lev):
                        if cellsHier_lev[i][H_CC]==NEXIST:
                            x_res[ind_st:ind_st+n_pind_v] = polyHier_v_lev[i]
                            ind_st += n_pind_cell
        if setSelfMatRes:
            self.mat_res = x_res
        return x_res
        

    def SetSolvedPolyCoeffsOnFace(self, v_ind, face, n_dir, der_order=None):
        """ der_order: standard form. for example np.array([2, 0, 1]) for d^3/dx^2dz in 3D
        """
        N = self.N
        cellsHier = self.cellsHier
        n_lev = len(cellsHier)
        n_vars = len(self.vars_list)
        
        der_order_tup = None
        if der_order==None:
            der_order_tup = tuple(np.zeros(N, dtrype=int).tolist())
        elif isinstance(der_order, (list, tuple)):
            der_order_tup = tuple(der_order)
        else:
            assert isinstance(der_order, np.adarray)
            
        PolyFaceCoeffs = [None]*n_lev
        PolyFaceInds0 = [None]*n_lev
        for lev in range(n_lev):
            cells_lev = cellsHier[lev]
            bc_conns_lev_face_dir = self.boundaryCellConnections[lev][face][n_dir]
            #print('bc_conns_lev_face_dir:', bc_conns_lev_face_dir)
            PolyFaceCoeffs[lev] = [None]*len(bc_conns_lev_face_dir)

            arr_I_FL_f0_dir, inds_0_FL_f0_dir = self.pf_dic[face][n_dir][lev][der_order_tup]
            #print('arr_I_FL_f0_dir:', arr_I_FL_f0_dir)
            PolyFaceInds0[lev] = inds_0_FL_f0_dir
            for i in range(len(bc_conns_lev_face_dir)):
                c_ind = bc_conns_lev_face_dir[i]
                if cells_lev[c_ind][H_CI]>=0:
                    poly_F = self.polyHier[v_ind][lev][c_ind]
                    #get poly coeffs on face
                    n_term_F = len(arr_I_FL_f0_dir)
                    #print('poly_F:', poly_F)

                    PolyFaceCoeffs[lev][i] = [None]*n_term_F
                    for j in range(n_term_F):
                        p_coeff_j = 0
                        for p_ind in arr_I_FL_f0_dir[j]:
                            p_coeff_j += arr_I_FL_f0_dir[j][p_ind]*poly_F[p_ind]
                        
                        PolyFaceCoeffs[lev][i][j] = p_coeff_j
                    #print('PolyFaceCoeffs[lev][i]:', PolyFaceCoeffs[lev][i])

        return [PolyFaceCoeffs, PolyFaceInds0]

    
    def SetVarsStartIndexInOutputVec(self):
        if self.N_indexed_acc==None:
            self.GetAccumukativeNumOfIndexedCells()
    
        n_vars = len(self.vars_list)
        cellsHier = self.cellsHier
        n_lev = len(cellsHier)
        n_pind_cell = sum(self.n_unkn_v)
        
        n_pind_v_acc = [0]*(n_vars+1)
        for v in range(n_vars):
            n_pind_v_acc[v+1] = n_pind_v_acc[v] + self.n_unkn_v[v]
        
        if self.indexing=='children':
            VarsStartIndex = [None]*n_vars
            for v in range(n_vars):
                VarsStartIndex[v] = [None]*n_lev
                n_pind_v = self.n_unkn_v[v]
                ind_st = n_pind_v_acc[v]
                for lev in range(n_lev):
                    n_cell_lev = len(cellsHier[lev])
                    cellsHier_lev = cellsHier[lev]
                    VarsStartIndex[v][lev] = [None]*n_cell_lev
                    for i in range(n_cell_lev):
                        if cellsHier_lev[i][H_CC]==NEXIST:
                            VarsStartIndex[v][lev][i] = ind_st
                            ind_st += n_pind_cell
            self.VarsStartIndex = VarsStartIndex
            if self.verbose>0:
                print('self.VarsStartIndex:', self.VarsStartIndex)
        

    def CategorizeFirstLevelCells(self):
        N = self.N
        cellsHier = self.cellsHier
        nx = self.nx
        nx_pow = self.nx_pow

        n_0 = np.zeros(N)
        n_1 = self.W
        
        i__ctr, i__ch, i__cells, i__nodes, i__nodes_in = list(range(5))

        nx_arr = [[np.zeros(N, dtype=int), np.copy(nx)], [], [], [], []]
        def processNxArr(nx_i):
            nx_0, nx_1 = nx_i[i__ctr]
            if np.any( (nx_1-nx_0)>=2 ):
                nx_mid = np.ceil((nx_0+nx_1)/2.0).astype(int)
                b_count = np.zeros(N, dtype=int)
                while True:
                    nx_0_next = nx_0 + b_count*(nx_mid-nx_0)
                    nx_1_next = nx_mid + b_count*(nx_1-nx_mid)
                    if np.all( (nx_1_next-nx_0_next)>=1):
                        nx_i[i__ch].append([[nx_0_next, nx_1_next], [], [], [], []])
                    if not self._binaryCounterIncrease(b_count):
                        break
                for i in range(len(nx_i[i__ch])):
                    processNxArr(nx_i[i__ch][i])

        processNxArr(nx_arr)
        #print('nx_arr:', nx_arr)
                
                
        def insertCind(c_ind, ctr, nx_i):
            nx_0, nx_1 = nx_i[i__ctr]
            if np.all((ctr>=nx_0)*(ctr<nx_1)):
                nx_i[i__cells].append(c_ind)
                for i in range(len(nx_i[i__ch])):
                    insertCind(c_ind, ctr, nx_i[i__ch][i])
                
        counter = np.zeros(N, dtype=int)
        while True:
            cell_ind = self._uniformGridCellInd(counter, nx, nx_pow)
            insertCind(cell_ind, counter, nx_arr)
            if not self._increaseCellCounterIndex(counter, nx):
                break
        
        def setnxarrCornerNodes(nx_i):
            if len(nx_i[i__cells])>0:
                n_min = np.copy(n_1)
                n_max = np.copy(n_0)
                #print('n_min:', n_min, 'n_max:', n_max)
                for c_ind in nx_i[i__cells]:
                    #print(c_ind)
                    nodes = cellsHier[0][c_ind][H_CN]
                    for n in nodes:
                        n_point = self.nodesPoints[n]
                        for n_dir in range(N):
                            if n_point[n_dir]<n_min[n_dir]:
                                n_min[n_dir] = n_point[n_dir]
                            if n_point[n_dir]>n_max[n_dir]:
                                n_max[n_dir] = n_point[n_dir]
                nx_i[i__nodes] = [n_min, n_max]
                #print('n_min:', n_min, 'n_max:', n_max)
            for i in range(len(nx_i[i__ch])):
                setnxarrCornerNodes(nx_i[i__ch][i])
        
        setnxarrCornerNodes(nx_arr)
        
        #print('nx_arr:', nx_arr)
        def printnxarr(nx_i):
            printlist = [nx_i]
            for i in range(len(nx_i[i__ch])):
                printlist.append(nx_i[i__ch][i])
            for nx_ij in printlist:
                print(nx_ij[i__ctr], nx_ij[i__nodes], nx_ij[i__cells], nx_ij[i__nodes_in])
            print('-'*40)
            for i in range(len(nx_i[i__ch])):
                printnxarr(nx_i[i__ch][i])
        
        #printnxarr(nx_arr)
        self.nx_arr = nx_arr
        
        
        
    def GetFirstLevelBoundingCells(self, points):
        N = self.N
        i__ctr, i__ch, i__cells, i__nodes, i__nodes_in = list(range(5))
        n_pts = len(points)
        nx_arr = self.nx_arr
        n_0, n_1 = nx_arr[i__nodes]
        nodes_in = nx_arr[i__nodes_in]
        dx_eps = self.dx_levels[0]*1.0e-5
        for i in range(n_pts):
            pt = points[i]
            if np.all((pt>=n_0)*(pt<=n_1)):
                nodes_in.append(i)
            else:
                pt_adjust_isin = False
                for n_dir in range(N):
                    mask = np.zeros(N)
                    mask[n_dir] = 1
                    pt_epsp = pt + mask*dx_eps
                    pt_epsm = pt - mask*dx_eps
                    if np.all((pt_epsp>=n_0)*(pt_epsp<=n_1)) or np.all((pt_epsm>=n_0)*(pt_epsm<=n_1)):
                        nodes_in.append(i)
                        pt_adjust_isin = True
                        break
                if not pt_adjust_isin:
                    raise NotImplementedError()
        
        def SetNodesInChildren(nx_i):
            if len(nx_i[i__ch])==0:
                return
            nodes_in_p = nx_i[i__nodes_in]
            added = [False]*len(nodes_in_p)
            not_added = []
            for i in range(len(nx_i[i__ch])):
                n_0i, n_1i = nx_i[i__ch][i][i__nodes]
                nodes_in_i = nx_i[i__ch][i][i__nodes_in]
                for j in range(len(nodes_in_p)):
                    if added[j]:
                        continue
                    i_pt = nodes_in_p[j]
                    pt = points[i_pt]
                    if np.all((pt>=n_0i)*(pt<=n_1i)):
                        nodes_in_i.append(i_pt)
                        added[j] = True

            for i in range(len(nodes_in_p)):
                if not added[i]:
                    not_added.append(nodes_in_p[i])

            if len(not_added)>0:
                print('nodes_in_p:', nodes_in_p)
                print('not_added:', not_added)
            
            nodes_in_p = not_added
            added = [False]*len(nodes_in_p)
            not_added = []
            for i in range(len(nx_i[i__ch])):
                n_0i, n_1i = nx_i[i__ch][i][i__nodes]
                nodes_in_i = nx_i[i__ch][i][i__nodes_in]
                for j in range(len(nodes_in_p)):
                    if added[j]:
                        continue
                    i_pt = nodes_in_p[j]
                    pt = points[i_pt]
                    for n_dir in range(N):
                        mask = np.zeros(N)
                        mask[n_dir] = 1
                        pt_epsp = pt + mask*dx_eps
                        pt_epsm = pt - mask*dx_eps
                        if np.all((pt_epsp>=n_0i)*(pt_epsp<=n_1i)) or np.all((pt_epsm>=n_0i)*(pt_epsm<=n_1i)):
                            nodes_in_i.append(i_pt)
                            added[j] = True
                            break
            for i in range(len(nodes_in_p)):
                if not added[i]:
                    not_added.append(nodes_in_p[i])
            assert len(not_added)==0
            if len(not_added)>0:
                print('not_added:', not_added)
                for i in range(len(not_added)):
                    print(points[not_added[i]])
                print('\n')
            
            for i in range(len(nx_i[i__ch])):
                SetNodesInChildren(nx_i[i__ch][i])
        #
        SetNodesInChildren(nx_arr)
        #
        def printnxarr(nx_i):
            printlist = [nx_i]
            for i in range(len(nx_i[i__ch])):
                printlist.append(nx_i[i__ch][i])
            for nx_ij in printlist:
                print(nx_ij[i__ctr], nx_ij[i__nodes], nx_ij[i__cells], nx_ij[i__nodes_in])
            print('-'*40)
            for i in range(len(nx_i[i__ch])):
                printnxarr(nx_i[i__ch][i])
        #
        #printnxarr(nx_arr)
        boundingCells = [None]*n_pts
        def setboundingcells(nx_i):
            if len(nx_i[i__ch])==0:
                assert len(nx_i[i__cells])==1
                c_ind = nx_i[i__cells][0]
                nodes_in = nx_i[i__nodes_in]
                for n in nodes_in:
                    assert boundingCells[n]==None
                    boundingCells[n] = c_ind
            else:
                for i in range(len(nx_i[i__ch])):
                    setboundingcells(nx_i[i__ch][i])
        setboundingcells(nx_arr)
        assert None not in boundingCells
        #print('zip(points, boundingCells):', list(zip(points, boundingCells)))
        return boundingCells
        

    def CleanNxArrPoints(self):
        i__ctr, i__ch, i__cells, i__nodes, i__nodes_in = list(range(5))
        nx_arr = self.nx_arr
        ##
        def cleannxarr(nx_i):
            nx_i[i__nodes_in] = []
            for i in range(len(nx_i[i__ch])):
                cleannxarr(nx_i[i__ch][i])
        cleannxarr(nx_arr)
        return    


    def GetLeafLevelCellsContainingPoints(self, points):
        N = self.N
        _2_pn_1 = 2**N-1
        if self.nx_arr==None:
            self.CategorizeFirstLevelCells()
            
        self.CleanNxArrPoints()
        boundingCells_0 = self.GetFirstLevelBoundingCells(points)
        
        n_pts = len(points)
        boundingCells = [None]*n_pts
        for i in range(n_pts):
            boundingCells[i] = (0, boundingCells_0[i])
        bc_isset = [False]*n_pts
        
        cellsHier = self.cellsHier
        n_lev = len(cellsHier)
        for i_pt in range(n_pts):
            pt = points[i_pt]
            for lev in range(n_lev):
                c_lev, c_ind = boundingCells[i_pt]
                c_ch = cellsHier[c_lev][c_ind][H_CC]
                if c_ch==NEXIST:
                    bc_isset[i_pt] = True
                    continue
                pt_isin = False
                for c in c_ch:
                    n_0 = self.nodesPoints[cellsHier[lev+1][c][H_CN][0]]
                    n_1 = self.nodesPoints[cellsHier[lev+1][c][H_CN][_2_pn_1]]
                    dx = self.dx_levels[lev+1]
                    #print(n_0, n_1, dx)
                    assert np.all((n_1-n_0)>dx/10.0)
                    if np.all((pt>=n_0)*(pt<=n_1)):
                        boundingCells[i_pt] = (lev+1, c)
                        pt_isin = True
                        break
                if not pt_isin:
                    for c in c_ch:
                        n_0 = self.nodesPoints[cellsHier[lev+1][c][H_CN][0]]
                        n_1 = self.nodesPoints[cellsHier[lev+1][c][H_CN][_2_pn_1]]
                        dx = self.dx_levels[lev+1]
                        #assert np.all((n_1-n_0)>dx/10.0)
                        dx_eps = dx*1.0e-5
                        for n_dir in range(N):
                            mask = np.zeros(N)
                            mask[n_dir] = 1.0
                            pt_epsp = pt + mask*dx_eps
                            pt_epsm = pt - mask*dx_eps
                            if np.all((pt_epsp>=n_0)*(pt_epsp<=n_1)) or np.all((pt_epsm>=n_0)*(pt_epsm<=n_1)):
                                boundingCells[i_pt] = (lev+1, c)
                                pt_isin = True
                                break
                        if pt_isin:
                            break
                    assert pt_isin==True
        assert False not in bc_isset
        return boundingCells
        
        
    ##TODO: calculate polynomials efficiently
    """ let the grid select the points, such that the points inside all cells in
    a given level have equal displacements. precalculate the polynomial powers.
    """
    def GetVariableValueAtPointsInsideCell(self, var_ind, lev, c_ind, points, printit=False):
        N = self.N
        dx = self.dx_levels[lev]
        n_pts = len(points)
        assert self.cellsHier[lev][c_ind][H_CC]==NEXIST
        n_0 = self.nodesPoints[self.cellsHier[lev][c_ind][H_CN][0]]
        pts_0 = [None]*n_pts
        for i in range(n_pts):
            pts_0[i] = points[i]-n_0
        
        poly = self.polyHier[var_ind][lev][c_ind]
        poly_vals = [0.0]*n_pts
        inds_0_all_F = self.inds_0_all_F
        n_p = len(inds_0_all_F)
        
        pts_xyz = np.array(pts_0).T
        
        poly_pows = np.zeros((n_p, n_pts))
        for i_p in range(n_p):
            p_ij = inds_0_all_F[i_p][1]
            p_ij_val = np.ones(n_pts)            
            for n in range(N):
                p_ij_val *= (pts_xyz[n,:]/dx[n])**p_ij[n]
            poly_pows[i_p,:] = p_ij_val
            
        poly_vals = poly.dot(poly_pows)
        poly_vals = poly_vals.reshape(n_pts, order='C')
        if printit:
            print('inds_0_all_F:', inds_0_all_F, 'poly:', poly, 'poly_pows.T:', poly_pows.T, \
                'points:', points, 'poly_vals:', poly_vals, sep='\n')
        return poly_vals
                        
    
    def GetVarValuesAtGivenPoints(self, var_ind, points):
        bound_cells = self.GetLeafLevelCellsContainingPoints(points)
        n_pts = len(points)
        var_vals = np.zeros(n_pts)
        for i in range(n_pts):
            lev, c_ind = bound_cells[i]
            pt_i = points[i]
            printit = False
            #if i%int(n_pts/23)==7:
            #    printit=True
            #    print('lev:', lev, 'c_ind:', c_ind)
            var_vals[i] = self.GetVariableValueAtPointsInsideCell(var_ind, lev, c_ind, [pt_i], printit=printit)[0]
        return var_vals


    def GetVarValuesOnMesh(self, var_ind, n_pts_dim, gslice=None):
        """ gslice = [{'dir':0, 'val':0.5}] mesh grid the plane x=0.5
        """
        mesh = self.GetMeshGrid(n_pts_dim, gslice)
        points = self.LayMeshAlong1D(mesh, gslice)
        var_vals = self.GetVarValuesAtGivenPoints(var_ind, points)
        if gslice==None:
            var_vals_mesh = var_vals.reshape(tuple(n_pts_dim), order='C')
            return [mesh, var_vals_mesh]
        else:
            N = self.N
            dim_slice = [d['dir'] for d in gslice]
            val_slice = [d['val'] for d in gslice]
            assert len(set(dim_slice))==len(dim_slice)
            N_R = N-len(gslice)
            assert N_R>=1
            assert len(n_pts_dim)==N
            x_1d = [None]*N_R
            n_pts_dim_R = [n_pts_dim[i] for i in range(N) if i not in dim_slice]

            var_vals_mesh = var_vals.reshape(tuple(n_pts_dim_R), order='C')
            return [mesh, var_vals_mesh]
                

    def GetMeshGrid(self, n_pts_dim, gslice=None):
        if gslice==None:
            N = self.N
            assert len(n_pts_dim)==N
            x_1d = [None]*N
            for i in range(N):
                x_1d[i] = np.linspace(0.0, self.W[i], n_pts_dim[i])
            if N==1:
                return x_1d
            else:
                mesh = np.meshgrid(*tuple(x_1d), indexing='ij')
                return mesh
        else:
            N = self.N
            dim_slice = [d['dir'] for d in gslice]
            val_slice = [d['val'] for d in gslice]
            assert len(set(dim_slice))==len(dim_slice)
            N_R = N-len(gslice)
            assert N_R>=1
            assert len(n_pts_dim)==N
            x_1d = [None]*N_R
            n_pts_dim_R = [n_pts_dim[i] for i in range(N) if i not in dim_slice]
            W_R = [self.W[i] for i in range(N) if i not in dim_slice]
            for i in range(N_R):
                x_1d[i] = np.linspace(0.0, W_R[i], n_pts_dim_R[i])
            if N_R==1:
                return x_1d
            else:
                mesh = np.meshgrid(*tuple(x_1d), indexing='ij')
                return mesh
                
    
    def LayMeshAlong1D(self, mesh, gslice=None):
        if gslice==None:
            N = self.N
            assert len(mesh)==N
            shape = mesh[0].shape
            n_tot = 1
            for i in range(N):
                n_tot *= shape[i]
            
            pts_1d = [None]*n_tot
            mesh_1d = [None]*N
            for i in range(N):
                mesh_1d[i] = mesh[i].reshape(n_tot, order='C')
            
            points = [None]*n_tot
            for i in range(n_tot):
                points[i] = np.zeros(N)
                for n in range(N):
                    points[i][n] = mesh_1d[n][i]
                    
            return points
        else:
            N = self.N
            dim_slice = [d['dir'] for d in gslice]
            val_slice = [d['val'] for d in gslice]
            N_R = N-len(gslice)
            assert len(mesh)==N_R
            shape = mesh[0].shape
            n_tot = 1
            for i in range(N_R):
                n_tot *= shape[i]
            
            pts_1d = [None]*n_tot
            mesh_1d = [None]*N_R
            for i in range(N_R):
                mesh_1d[i] = mesh[i].reshape(n_tot, order='C')
            
            
            dims = list(range(N))
            dims_TF = [True]*N
            dims_tomeshinds = [-1]*N
            m_ind = 0
            for i in range(N):
                if i in dim_slice:
                    dims_TF[i]=False
                else:
                    dims_tomeshinds[i] = m_ind
                    m_ind+=1
            dim_vals = [0]*N
            for i in range(len(dim_slice)):
                dim_vals[dim_slice[i]] = val_slice[i]
                    
            points = [None]*n_tot
            for i in range(n_tot):
                points[i] = np.zeros(N)
                for n in range(N):
                    if dims_TF[n]:
                        points[i][n] = mesh_1d[dims_tomeshinds[n]][i]
                    else:
                        points[i][n] = dim_vals[n]
                        
            return points
    
    
    def SpreadPointsOnMeshgrid(self, points, n_pts_dim, var_vals=None):
        N = self.N
        assert len(n_pts_dim)==N
        n_tot = 1
        for i in range(N):
            n_tot *= n_pts_dim[i]
        assert len(points)==n_tot
        
        mesh_1d = [None]*N
        for n in range(N):
            mesh_1d[n] = np.zeros(n_tot)
            for i in range(n_tot):
                mesh_1d[n][i] = points[i][n]
        
        mesh = [None]*N
        for n in range(N):
            mesh[n] = mesh_1d[n].reshape(tuple(n_pts_dim), order='C')        
        
        if var_vals!=None:
            assert len(var_vals)==n_tot
            var_vals_mesh = var_vals.reshape(tuple(n_pts_dim), order='C')
            return [mesh, var_vals_mesh]
        else:
            return mesh
        
    
    def GetPolyArrayDifferentiated(self, poly_order, der_order):
        """ diff( sum(a_mn*x^m*y*n )
        """
        assert len(poly_order)==len(der_order)
        poly_order_p1 = tuple([i+1 for i in poly_order])
        arr = np.ones(poly_order_p1)
        shape = arr.shape
        for i in range(len(der_order)):
            a = np.ones(shape[i])
            for j in range(der_order[i]):
                a *= np.arange(shape[i])-j
            #print('i:{} a:\n'.format(i), a)
            for j in range(shape[i]):
                indx = [Ellipsis]*arr.ndim
                indx[i] = j
                arr[indx] *= a[j]
                #print('der -- j:{}  arr:\n'.format(j), arr)
        for i in range(len(der_order)):
            arr = np.roll(arr, -der_order[i], i)
            #print('rolling -- i:{}  arr:\n'.format(i), arr)
        return arr
        

    def GetPolyArrayDifferentiated_scaled(self, poly_order, der_order, level):
        """ diff( sum(a_mn*(x/h_x)^m*(y/h_y)*n )
            h_x, h_y... : cell dimensions
        """
        assert len(poly_order)==len(der_order) and len(der_order)==self.N
        dx_lev = self.dx_levels[level]
        poly_order_p1 = tuple([i+1 for i in poly_order])
        arr = np.ones(poly_order_p1)
        shape = arr.shape
        for i in range(len(der_order)):
            a = np.ones(shape[i])
            for j in range(der_order[i]):
                a *= (np.arange(shape[i])-j)/dx_lev[i]
            #print('i:{} a:\n'.format(i), a)
            for j in range(shape[i]):
                indx = [Ellipsis]*arr.ndim
                indx[i] = j
                arr[indx] *= a[j]
                #print('der -- j:{}  arr:\n'.format(j), arr)
        for i in range(len(der_order)):
            arr = np.roll(arr, -der_order[i], i)
            #print('rolling -- i:{}  arr:\n'.format(i), arr)
        return arr
       
    def GetIndsToRemove(self, poly_order, orders_max):
        N = self.N
        poly_order_p1 = tuple([i+1 for i in poly_order])
        inds_to_rem = np.zeros(poly_order_p1, dtype=int)
        if N==1:
            return inds_to_rem
        n_rem = np.array(poly_order) - np.array(orders_max)
        for i in range(N):
            for j in range(n_rem[i]):
                indx = [Ellipsis]*N
                indx[i] = j
                inds_to_rem[indx] = 1
        changed = True
        while changed:
            changed = False
            counter = np.zeros(self.nx.shape, dtype=int)
            while True:
                for i in range(N):
                    if counter[i]>=n_rem[i]:
                        c_i = np.copy(counter)
                        c_i[i] -= n_rem[i]
                        valid_ind = True
                        for j in range(N):
                            c_j = np.copy(c_i)
                            if j!=i:
                                if c_j[j]+n_rem[j]>=poly_order_p1[j]:
                                    valid_ind = False
                                else:
                                    c_j[j] += n_rem[j]
                                    if inds_to_rem[tuple(c_j.tolist())]==0:
                                        valid_ind = False
                        if valid_ind:
                            if inds_to_rem[tuple(counter.tolist())]==0:
                                changed = True
                                inds_to_rem[tuple(counter.tolist())] = 1
                if not self._increasePolyCounterIndex(counter, poly_order):
                    break
        inds_to_rem = (inds_to_rem==0)*1
        #print('inds_to_rem:', inds_to_rem, sep='\n')

        """
        inds_to_rem = np.ones(poly_order_p1, dtype=int)
        for i in range(N):
            for j in range(orders_max[i]+1, poly_order[i]+1):
                indx = [Ellipsis]*N
                indx[i] = j
                inds_to_rem[indx] *= 0
        n_rem = np.array(poly_order) - np.array(orders_max)
        for i in range(N):
            for j in range(n_rem[i]):
                indx = [Ellipsis]*N
                indx[i] = j
                inds_to_rem[indx] += 1
        inds_to_rem = (inds_to_rem==0)*1
        print('inds_to_rem:', inds_to_rem, sep='\n')
        """
        """
        ##remove top-right corner block for 2D or higher
        counter = np.zeros(self.nx.shape, dtype=int)
        inds_to_rem = np.zeros(poly_order_p1)
        rem_porder = np.array(poly_order) - np.array(orders_max) - 1
        #print('rem_porder:', rem_porder)
        rem_ind_init = np.array(orders_max)+1 
        while True:
            inds_to_rem[tuple((rem_ind_init+counter).tolist())] = 1
            if not self._increasePolyCounterIndex(counter, rem_porder):
                break
        print('inds_to_rem: ', inds_to_rem, sep='\n')
        """
        return inds_to_rem
            

    def SetEqIndicesToKeep(self, poly_order, orders_max):
        """
        -consider all derivatives
        -add extra rows and columns containing -1
        -perform derivatives and move indices around for all the derivatives in 
        the equation 
        -mark cells that containg no -1 and at least one x
        -assign resulting poly equations to cells containing x (all rows filled)
        -if there are more equations than xs, assign more poly equations to some 
        xs and average over those equations
        -make sure all xs are present in the final equation sets (are cols are filled)
        """
        N = self.N
        poly_order_p1 = tuple([i+1 for i in poly_order])
        n_rem = np.array(poly_order) - np.array(orders_max)
        inds_to_rem = self.GetIndsToRemove(poly_order, orders_max)
        inds_keep = (inds_to_rem==0)*1
        shape = inds_keep.shape
        shape_p = tuple([shape[i]+n_rem[i] for i in range(N)])
        inds_0_all = -np.ones(shape_p, dtype=int)
        mask_x = np.zeros(shape_p, dtype=int)
        mask_all = np.zeros(shape_p, dtype=int)
        inds_0_all_F = []
        mask_x_F = []
        counter = np.zeros(N, dtype=int)
        ind = 0
        while True:
            counter_tup = tuple((counter).tolist())
            if inds_keep[counter_tup]:
                inds_0_all[counter_tup] = ind
                inds_0_all_F.append([ind, counter_tup])
                mask_all[counter_tup] = 1
                if np.all(counter>=n_rem):
                    mask_x[counter_tup] = 1
                    mask_x_F.append(1)
                else:
                    mask_x_F.append(0)
                ind += 1
            if not self._increasePolyCounterIndex(counter, poly_order):
                break
        n_ind = ind
        assert len(inds_0_all_F)==len(mask_x_F)==n_ind
        for i in range(n_ind):
            assert inds_0_all_F[i][0] == i
        if self.verbose>0:
            print('n_ind:', n_ind, 'inds_keep:', inds_keep, 'inds_0_all:', inds_0_all,\
        'mask_all:', mask_all, 'mask_x:', mask_x, sep='\n')
        if self.verbose>0:
            print('inds_0_all_F:', inds_0_all_F, 'mask_x_F:', mask_x_F, sep='\n')
                
        eqs_inds_all = [None]*len(self.EQs_parts)
                
        for eq_ind, eq_arr in enumerate(self.EQs_parts):
            #eq_arr: a differential equation in the input set of diff eqs
            if self.verbose>0:
                print('eq_ind:', eq_ind, '  eq_arr:', eq_arr)
            eqs_inds_all[eq_ind] = [None]*len(eq_arr)
            for v_ind, eq_v in enumerate(eq_arr):
                #eq_v: the differential equation associated with the v-th variable
                if self.verbose>0:
                    print('v_ind:', v_ind, '  eq_v:', eq_v)
                if eqs_inds_all[eq_ind][v_ind]==None:
                    eqs_inds_all[eq_ind][v_ind] = [None]*n_ind
                    for i in range(n_ind):
                        eqs_inds_all[eq_ind][v_ind][i] = []
                for eq_part in eq_v:
                    #eq_part: an additive part of the differential equation
                    if self.verbose>0:
                        print('eq_part:', eq_part)
                    coeff_sym, der_ord, der_vars = eq_part
                    
                    display(Math(latex(eq_part)+' : ' + latex(coeff_sym) + \
                    ';' + latex(der_ord) + ';' + latex(der_vars)))
                    
                    der_orders = self.ChangeDerOrderAndDerVarsToStandardForm(der_ord, der_vars)
                    
                    a = np.copy(inds_0_all)
                    for i in range(N):
                        a = np.roll(a, -der_orders[i], i)
                    
                    for i in range(n_ind):
                        ind, ind_tup = inds_0_all_F[i]
                        eqs_inds_all[eq_ind][v_ind][ind].append(a[ind_tup])
                                       
        if self.verbose>0:
            print('eqs_inds_all:', eqs_inds_all, sep='\n')
        ##process eqs_inds_all
        n_ind_x = sum(mask_x_F)
        eqs_inds_keep = [None]*len(self.EQs_parts)
        for eq_ind in range(len(self.EQs_parts)):
            eqs_inds_keep[eq_ind] = [None]*len(self.EQs_parts[eq_ind])
            for v_ind in range(len(self.EQs_parts[eq_ind])):
                if eqs_inds_keep[eq_ind][v_ind]==None:
                    eqs_inds_keep[eq_ind][v_ind] = [None]*n_ind

                    for i in range(n_ind):
                        keep = []
                        for ind in eqs_inds_all[eq_ind][v_ind][i]:
                            if ind==-1:
                                keep = []
                                break
                            elif mask_x_F[ind]==1:
                                keep.append(ind)
                        if len(keep)==0:
                            keep = None
                        eqs_inds_keep[eq_ind][v_ind][i] = keep
        if self.verbose>0:
            print('n_ind_x:', n_ind_x, 'eqs_inds_keep:', eqs_inds_keep, sep='\n')
        
        x_assign = [None]*len(self.EQs_parts)
        for eq_ind in range(len(self.EQs_parts)):
            x_assign[eq_ind] = [None]*n_ind
            for i in range(n_ind):
                x_assign[eq_ind][i] = []
                
        n_try, n_try_max = 0, 100
        while True:
            for eq_ind in range(len(self.EQs_parts)):
                for i in range(n_ind):
                    if eqs_inds_keep[eq_ind][eq_ind][i]!=None:
                        for v_ind in range(len(self.EQs_parts[eq_ind])):
                            if v_ind!=eq_ind and eqs_inds_keep[eq_ind][v_ind]==None:
                                raise NotImplementedError()
                        if len(eqs_inds_keep[eq_ind][eq_ind][i])==1:
                            x_assign[eq_ind][eqs_inds_keep[eq_ind][eq_ind][i][0]].append(i)
                            del eqs_inds_keep[eq_ind][eq_ind][i][0]
                        elif len(eqs_inds_keep[eq_ind][eq_ind][i])>1:
                            for j in range(len(eqs_inds_keep[eq_ind][eq_ind][i])-1, -1, -1):
                                if len(x_assign[eq_ind][eqs_inds_keep[eq_ind][eq_ind][i][j]])>0:
                                    del eqs_inds_keep[eq_ind][eq_ind][i][j]
                                    if len(eqs_inds_keep[eq_ind][eq_ind][i])==1:
                                        break
            eqs_inds_keep_empty = True
            for eq_ind in range(len(self.EQs_parts)):
                for i in range(n_ind):
                    if eqs_inds_keep[eq_ind][eq_ind][i]!=None:
                        if len(eqs_inds_keep[eq_ind][eq_ind][i])>0:
                            eqs_inds_keep_empty = False
            n_try += 1
            if eqs_inds_keep_empty:
                break                    
            if n_try>n_try_max:
                print('eqs_inds_keep:', eqs_inds_keep, sep='\n')
                raise NotImplementedError()

        if self.verbose>0:
            print('n_try:', n_try, 'mask_x_F:', mask_x_F, 'x_assign:', x_assign, sep='\n')    
        for eq_ind in range(len(self.EQs_parts)):
            for i in range(n_ind):
                if mask_x_F[i]==1:
                    assert len(x_assign[eq_ind][i])>0
                else:
                    assert len(x_assign[eq_ind][i])==0
                    x_assign[eq_ind][i] = None
        if self.verbose>0:
            print('x_assign:', x_assign, sep='\n')    
        ##indices to keep, associated x index and weight factor
        inds_keep__x_weight = [None]*len(self.EQs_parts)
        for eq_ind in range(len(self.EQs_parts)):
            inds_keep__x_weight[eq_ind] = []
            for i in range(n_ind):
                if x_assign[eq_ind][i]!=None:
                    for j in range(len(x_assign[eq_ind][i])):
                        for k in range(len(inds_keep__x_weight[eq_ind])):
                            assert inds_keep__x_weight[eq_ind][k][0]!=x_assign[eq_ind][i][j]
                        inds_keep__x_weight[eq_ind].append([x_assign[eq_ind][i][j], i, 1.0/len(x_assign[eq_ind][i])])
        if self.verbose>0:
            print('inds_keep__x_weight:', inds_keep__x_weight, sep='\n')

        self.mask_allp = mask_all
        self.mask_xp = mask_x
        self.mask_x_F = mask_x_F
        self.inds_0_allp = inds_0_all
        self.inds_0_all_F = inds_0_all_F
        self.inds_keep = inds_keep
        self.inds_keep__x_weight = inds_keep__x_weight
        
        print('self.mask_allp:', self.mask_allp, 'self.mask_xp:', self.mask_xp, 
        'self.mask_x_F:', self.mask_x_F, 'self.inds_0_allp:', self.inds_0_allp, 
        'self.inds_0_all_F:', self.inds_0_all_F, 'self.inds_keep:', self.inds_keep,
        'self.inds_keep__x_weight:', self.inds_keep__x_weight, '-'*30, sep='\n')



    def GetPolyArrayDiff_scaled_indices_masked_flat(self, eq_ind, poly_order, der_order, orders_max, level):
        """ diff( sum(a_mn*(x/h_x)^m*(y/h_y)*n )
            h_x, h_y... : cell dimensions
            mask polynomial orders greater than or equal to mask_order
        """
        N = self.N
        assert len(poly_order)==N and len(der_order)==N and len(orders_max)==N
        dx_lev = self.dx_levels[level]
        poly_order_p1 = tuple([i+1 for i in poly_order])
        n_rem = np.array(poly_order) - np.array(orders_max)

        #shape = self.inds_keep.shape
        shape_p = self.mask_allp.shape

        arr = np.copy(self.mask_allp).astype(float)    ## arr[i][j][k] --> x^i*y^j*z^k    
        inds = np.copy(self.inds_0_allp)
        inds_0 = np.copy(self.inds_0_allp)
        
        ##take derivative 
        for i in range(N):
            if der_order[i]>0:
                a = np.ones(shape_p[i])
                for j in range(der_order[i]):
                    a *= (np.arange(shape_p[i])-j)/dx_lev[i]
                #print('i:{} a:\n'.format(i), a)
                for j in range(shape_p[i]):
                    indx = [Ellipsis]*arr.ndim
                    indx[i] = j
                    arr[indx] *= a[j]
                    #print('der -- j:{}  arr:\n'.format(j), arr)
        ##roll
        for i in range(N):
            arr = np.roll(arr, -der_order[i], i)
            inds = np.roll(inds, -der_order[i], i)
            #print('rolling -- i:{}  arr:\n'.format(i), arr)

        ##mask zeros
        arr *= self.mask_allp
        inds = inds*self.mask_allp - (self.mask_allp==0)*1
        
        #print('arr:', arr, 'inds:', inds, 'inds_0:', inds_0, sep='\n')
        
        inds_F, inds_0_F, arr_F = [], [], []
        ##inds to keep
        for i in range(len(self.inds_keep__x_weight[eq_ind])):
            i_, i_x_, w_ = self.inds_keep__x_weight[eq_ind][i]
            
            i_tup_ = self.inds_0_all_F[i_][1]      
            i_x_tup_ = self.inds_0_all_F[i_x_][1]  
            
            inds_F.append(inds[i_tup_])
            arr_F.append(w_*arr[i_tup_])

            inds_0_F.append(self.inds_0_allp[i_x_tup_])
            
        if self.verbose>0:
            print('arr_F:', arr_F, 'inds_F:', inds_F, 'inds_0_F:', inds_0_F, '-'*20 ,sep='\n')
        return [np.array(arr_F), np.array(inds_F), np.array(inds_0_F)]        
        

    def GetPolyArrayDiff_scaled_indices_masked_flat_multTemp(self, eq_ind, poly_order, der_order, orders_max, level):
        """ diff( sum(a_mn*(x/h_x)^m*(y/h_y)*n )
            h_x, h_y... : cell dimensions
            mask polynomial orders greater than or equal to mask_order
            multTemp: multiplication template
        """
        N = self.N
        assert len(poly_order)==N and len(der_order)==N and len(orders_max)==N
        dx_lev = self.dx_levels[level]
        poly_order_p1 = tuple([i+1 for i in poly_order])
        n_rem = np.array(poly_order) - np.array(orders_max)

        #shape = self.inds_keep.shape
        shape_p = self.mask_allp.shape

        arr = np.copy(self.mask_allp).astype(float)    ## arr[i][j][k] --> x^i*y^j*z^k    
        inds = np.copy(self.inds_0_allp)
        inds_0 = np.copy(self.inds_0_allp)
        
        ##take derivative 
        for i in range(N):
            if der_order[i]>0:
                a = np.ones(shape_p[i])
                for j in range(der_order[i]):
                    a *= (np.arange(shape_p[i])-j)/dx_lev[i]
                #print('i:{} a:\n'.format(i), a)
                for j in range(shape_p[i]):
                    indx = [Ellipsis]*arr.ndim
                    indx[i] = j
                    arr[indx] *= a[j]
                    #print('der -- j:{}  arr:\n'.format(j), arr)
        ##roll
        for i in range(N):
            arr = np.roll(arr, -der_order[i], i)
            inds = np.roll(inds, -der_order[i], i)
            #print('rolling -- i:{}  arr:\n'.format(i), arr)

        ##mask zeros
        arr *= self.mask_allp
        inds = inds*self.mask_allp - (self.mask_allp==0)*1
        
        #print('arr:', arr, 'inds:', inds, 'inds_0:', inds_0, sep='\n')
        
        ##TODO: multiply arr by another polynomial 
        n_F = len(self.inds_0_all_F)
        inds_F, inds_0_F, arr_F = [None]*n_F, [None]*n_F, [None]*n_F
        for i in range(n_F):
            inds_F[i] = inds[self.inds_0_all_F[i][1]]
            arr_F[i] = arr[self.inds_0_all_F[i][1]]
            inds_0_F[i] = inds_0[self.inds_0_all_F[i][1]]
        
        P0P1_M = self.P0P1_M
        arr_I_mul_F = [None]*n_F
        ## P0*arr
        for i in range(n_F):
            arr_I_mul_F[i] = []
            for j in range(len(P0P1_M[i])):
                ij_01, multip = P0P1_M[i][j]    ##inds (P0, P1), multiplicity
                #print('ij_01:', ij_01, 'multip:', multip)
                if inds_F[ij_01[1]]>=0 and arr_F[ij_01[1]]!=0:
                    arr_I_mul_F[i].append([(ij_01[0], inds_F[ij_01[1]]), multip*arr_F[ij_01[1]]])
                
        #print('inds_0_F:', inds_0_F, 'arr_I_mul_F:', arr_I_mul_F, '-'*20 , sep='\n')

        ##inds to keep
        inds_0_FL, arr_I_mul_FL = [], []
        for i in range(len(self.inds_keep__x_weight[eq_ind])):
            i_, i_x_, w_ = self.inds_keep__x_weight[eq_ind][i]
                        
            el = [None]*len(arr_I_mul_F[i_])
            for j in range(len(arr_I_mul_F[i_])):
                el[j] = [arr_I_mul_F[i_][j][0], w_*arr_I_mul_F[i_][j][1]]
            arr_I_mul_FL.append(el)

            inds_0_FL.append(i_x_)
            
        if self.verbose>0:
            print('inds_0_FL:', inds_0_FL, 'arr_I_mul_FL:', arr_I_mul_FL, '-'*20 , sep='\n')
        return [arr_I_mul_FL, np.array(inds_0_FL)]        



    def GetPolyArrayDiff_scaled_indices_masked_flat___ver_0(self, poly_order, der_order, orders_max, level):
        """ ** defective **
            diff( sum(a_mn*(x/h_x)^m*(y/h_y)*n )
            h_x, h_y... : cell dimensions
            mask polynomial orders greater than or equal to mask_order
        """
        #self.SetEqIndicesToKeep(poly_order, orders_max)
        N = self.N
        assert len(poly_order)==N and len(der_order)==N and len(orders_max)==N
        dx_lev = self.dx_levels[level]
        poly_order_p1 = tuple([i+1 for i in poly_order])
        arr = np.ones(poly_order_p1)    ## arr[i][j][k] --> x^i*y^j*z^k    
        shape = arr.shape
        n_tot = 1 
        for i in range(N):
            n_tot *= shape[i]
        inds = np.arange(n_tot, dtype=int)      ##displaced (by derivatives) indices
        inds = inds.reshape(shape, order='C')
        inds_0 = np.arange(n_tot, dtype=int)    ##original indices
        inds_0 = inds_0.reshape(shape, order='C')
        #print('level:', level, 'dx_lev', dx_lev, 'poly_order:', poly_order, \
        #'der_order:', der_order, 'orders_max:', orders_max, 'inds_0:', inds_0, sep='\n')
        ##take derivative 
        for i in range(N):
            a = np.ones(shape[i])
            for j in range(der_order[i]):
                a *= (np.arange(shape[i])-j)/dx_lev[i]
            #print('i:{} a:\n'.format(i), a)
            for j in range(shape[i]):
                indx = [Ellipsis]*arr.ndim
                indx[i] = j
                arr[indx] *= a[j]
                #print('der -- j:{}  arr:\n'.format(j), arr)
        ##roll
        for i in range(N):
            arr = np.roll(arr, -der_order[i], i)
            inds = np.roll(inds, -der_order[i], i)
            #print('rolling -- i:{}  arr:\n'.format(i), arr)

        ##mask zeros
        for i in range(N):
            for j in range(orders_max[i]+1, poly_order[i]+1):
                arr = np.delete(arr, -1, axis=i)
                inds = np.delete(inds, -1, axis=i)
        for i in range(N):
            for j in range(0, poly_order[i]-orders_max[i]):
                inds_0 = np.delete(inds_0, 0, axis=i)
        ##flatten
        shape_d = arr.shape
        n_tot_d = 1 
        for i in range(N):
            n_tot_d *= shape_d[i]
        arr = arr.reshape(n_tot_d, order='C')
        inds = inds.reshape(n_tot_d, order='C')
        inds_0 = inds_0.reshape(n_tot_d, order='C')
        if self.verbose>0:
            print('arr:', arr, 'inds', inds, 'inds_0', inds_0, sep='\n')
        if self.N>=2:
            ##remove top and right blocks for 2D or higher
            inds_to_rem = self.GetIndsToRemove(poly_order, orders_max)
            for i in range(N):
                for j in range(0, poly_order[i]-orders_max[i]):
                    inds_to_rem = np.delete(inds_to_rem, 0, axis=i)
            inds_to_rem = inds_to_rem.reshape(n_tot_d, order='C')
            #_rm: remove   _fb: forbidden indices        
            arr_rm, inds_rm, inds_0_rm = [], [], []
            inds_fb = []
            for i in range(n_tot_d):
                if inds_to_rem[i]==0:
                    arr_rm.append(arr[i])
                    inds_rm.append(inds[i])
                    inds_0_rm.append(inds_0[i])
                else:
                    inds_fb.append(inds_0[i])
            #verify
            for i_fb in inds_fb:
                assert i_fb not in inds_rm
            if self.verbose>0:
                print('inds_fb', inds_fb, 'arr_rm:', arr_rm, 'inds_rm', inds_rm, 'inds_0_rm', inds_0_rm, sep='\n')
            ##remap indices to 0..ind_final 
            inds0_to_fin_map = np.arange(n_tot, dtype=int)    ##original indices
            inds0_to_fin_map = inds0_to_fin_map.reshape(shape, order='C')
            counter = np.zeros(self.nx.shape, dtype=int)
            ind_map = 0
            while True:
                if inds0_to_fin_map[tuple(counter.tolist())] in inds_fb:
                    inds0_to_fin_map[tuple(counter.tolist())] = -1
                else:
                    inds0_to_fin_map[tuple(counter.tolist())] = ind_map
                    ind_map += 1
                if not self._increasePolyCounterIndex(counter, poly_order):
                    break
            inds0_to_fin_map = inds0_to_fin_map.reshape(n_tot, order='C')
            for i in range(len(inds_0_rm)):
                inds_0_rm[i] = inds0_to_fin_map[inds_0_rm[i]]
                inds_rm[i] = inds0_to_fin_map[inds_rm[i]]
                assert inds_0_rm[i]>=0 and inds_rm[i]>=0
            if self.verbose>0:
                print('inds0_to_fin_map:', inds0_to_fin_map)
                print('mapped', 'arr_rm:', arr_rm, 'inds_rm', inds_rm, 'inds_0_rm', inds_0_rm, sep='\n')
            return [np.array(arr_rm), np.array(inds_rm), np.array(inds_0_rm)]
            
        else:
            return [arr_rm, inds_rm, inds_0_rm]
        

    def GetPolyArrayRHSConst_masked_flat(self, eq_ind, poly_order, orders_max, coeff):
        N = self.N
        assert len(poly_order)==N and len(orders_max)==N
        poly_order_p1 = tuple([i+1 for i in poly_order])

        shape_p = self.mask_allp.shape

        arr = np.zeros(shape_p, dtype=complex)    ## arr[i][j][k] --> x^i*y^j*z^k    
        arr[tuple(np.zeros(N).tolist())] = coeff
        inds_0 = np.copy(self.inds_0_allp)

        #print('arr:', arr, 'inds_0', inds_0, sep='\n')
        ##flatten

        inds_0_F, arr_F = [], []
        ##inds to keep
        for i in range(len(self.inds_keep__x_weight[eq_ind])):
            i_, i_x_, w_ = self.inds_keep__x_weight[eq_ind][i]
            
            i_tup_ = self.inds_0_all_F[i_][1]      
            i_x_tup_ = self.inds_0_all_F[i_x_][1]  
            
            arr_F.append(w_*arr[i_tup_])
            inds_0_F.append(self.inds_0_allp[i_x_tup_])

        return [np.array(arr_F), np.array(inds_0_F)]


    def GetPolyArrayRHSConst_masked_flat_Temp(self, eq_ind, poly_order, orders_max):
        """ _Temp: template (for non constant RHSs)
        """
        N = self.N
        assert len(poly_order)==N and len(orders_max)==N
        poly_order_p1 = tuple([i+1 for i in poly_order])

        shape_p = self.mask_allp.shape

        arr_I = np.copy(self.inds_0_allp)    ##contain coeff indices
        inds_0 = np.copy(self.inds_0_allp)

        #print('arr:', arr, 'inds_0', inds_0, sep='\n')
        ##flatten

        inds_0_F, arr_I_F = [], []
        ##inds to keep
        for i in range(len(self.inds_keep__x_weight[eq_ind])):
            i_, i_x_, w_ = self.inds_keep__x_weight[eq_ind][i]
            
            i_tup_ = self.inds_0_all_F[i_][1]      
            i_x_tup_ = self.inds_0_all_F[i_x_][1]  
            
            arr_I_F.append((arr_I[i_tup_], w_)) ## (ind_coeff, weight)
            inds_0_F.append(self.inds_0_allp[i_x_tup_])

        return [arr_I_F, np.array(inds_0_F)]


    def GetPolyArrayRHS_masked_flat___ver_0(self, poly_order, orders_max, rhs_sym):
        """ **defective**
        """
        N = self.N
        assert len(poly_order)==N and len(orders_max)==N
        poly_order_p1 = tuple([i+1 for i in poly_order])
        coeff_is_constant = True
        for par in self.pars_list:
            if rhs_sym.has(par):
                coeff_is_constant = False
                break
        if coeff_is_constant:
            coeff = complex(rhs_sym)
            arr = np.ones(poly_order_p1)*coeff    ## arr[i][j][k] --> x^i*y^j*z^k    
            shape = arr.shape
            n_tot = 1 
            for i in range(N):
                n_tot *= shape[i]
            inds_0 = np.arange(n_tot, dtype=int)    ##original indices
            inds_0 = inds_0.reshape(shape, order='C')
            #print('arr:', arr, 'inds_0', inds_0, sep='\n')
            ##mask outside orders_max
            for i in range(N):
                for j in range(orders_max[i]+1, poly_order[i]+1):
                    arr = np.delete(arr, -1, axis=i)
            for i in range(N):
                for j in range(0, poly_order[i]-orders_max[i]):
                    inds_0 = np.delete(inds_0, 0, axis=i)
            ##flatten
            shape_d = arr.shape
            n_tot_d = 1 
            for i in range(N):
                n_tot_d *= shape_d[i]
            arr = arr.reshape(n_tot_d, order='C')
            inds_0 = inds_0.reshape(n_tot_d, order='C')
            #print('arr:', arr, 'inds_0', inds_0, sep='\n')
            if self.N>=2:
                inds_to_rem = self.GetIndsToRemove(poly_order, orders_max)
                for i in range(N):
                    for j in range(0, poly_order[i]-orders_max[i]):
                        inds_to_rem = np.delete(inds_to_rem, 0, axis=i)
                inds_to_rem = inds_to_rem.reshape(n_tot_d, order='C')
                #_rm: remove   _fb: forbidden indices        
                arr_rm, inds_0_rm = [], []
                inds_fb = []
                for i in range(n_tot_d):
                    if inds_to_rem[i]==0:
                        arr_rm.append(arr[i])
                        inds_0_rm.append(inds_0[i])
                    else:
                        inds_fb.append(inds_0[i])
                ##remap indices to 0..ind_final 
                inds0_to_fin_map = np.arange(n_tot, dtype=int)    ##original indices
                inds0_to_fin_map = inds0_to_fin_map.reshape(shape, order='C')
                counter = np.zeros(self.nx.shape, dtype=int)
                ind_map = 0
                while True:
                    if inds0_to_fin_map[tuple(counter.tolist())] in inds_fb:
                        inds0_to_fin_map[tuple(counter.tolist())] = -1
                    else:
                        inds0_to_fin_map[tuple(counter.tolist())] = ind_map
                        ind_map += 1
                    if not self._increasePolyCounterIndex(counter, poly_order):
                        break
                inds0_to_fin_map = inds0_to_fin_map.reshape(n_tot, order='C')
                for i in range(len(inds_0_rm)):
                    inds_0_rm[i] = inds0_to_fin_map[inds_0_rm[i]]
                    assert inds_0_rm[i]>=0
                
                return [np.array(arr_rm), np.array(inds_0_rm)]
            else:
                return [arr, inds_0]
        else:
            raise NotImplementedError()

    def ChangeDerOrderAndDerVarsToStandardForm(self, der_ord, der_vars):
        """ der_ord=3, der_vars=[x,x,z] ---> [2, 0, 1] or double derivative with
        respect to x, no derivative with respect to y, single derivative with 
        respect to z
        """
        der_orders = [0]*self.N
        for i in range(self.N):
            n_i = 0
            x_i = self.indepVars_list[i]
            for v in der_vars:
                if v==x_i:
                    n_i += 1
            der_orders[i] = n_i
        assert sum(der_orders)==der_ord
        return der_orders
    

    def GetCellContCondEqs_scaled_masked(self, poly_order, polyorders_max, der_ord_max):
        """ apply the continuity condition for the given derivative order in 
        the given direction. The derivative is assumed in the given direction,
        and the continuity is applied on the face normal to this direction.
        The zero-th derivative fills the available place at index 0
        the first derivative fills the available place at index 1 and so on...
        for a double derivative in the x and y directions the available indices
        are as follows, marked as 0, 1,..
        equations on the shared indices are averaged. 
        |1 1 1 01 11|
        |0 0 0 00 10|
        |x x x 0  1 |
        |x x x 0  1 |
        |x x x 0  1 |
        
        an output is generated for each key in cells_nb_dic
        """
        N = self.N
        cellsHier = self.cellsHier
        nodesPoints = self.nodesPoints
        n_lev = len(cellsHier)
        cc_dic = {'n':{'sl':[None]*n_lev, 'nl':[None]*n_lev, 'pl':[None]*n_lev}}
        cells_nb_dic_neg_sl = self.cells_nb_dic['n']['sl']
        cells_nb_dic_neg_nl = self.cells_nb_dic['n']['nl']
        cells_nb_dic_neg_pl = self.cells_nb_dic['n']['pl']
        
        for lev in range(n_lev):
            cc_dic['n']['sl'][lev] = [None]*N
            cc_dic['n']['nl'][lev] = [None]*N
            cc_dic['n']['pl'][lev] = [None]*N
            cells_nb_dic_neg_sl_lev = cells_nb_dic_neg_sl[lev]
            cells_nb_dic_neg_nl_lev = cells_nb_dic_neg_nl[lev]
            cells_nb_dic_neg_pl_lev = cells_nb_dic_neg_pl[lev]
            dx_lev = np.array(self.dx_levels[lev])
            dx_levp1 = np.array(self.dx_levels[lev+1])
            for n_dir in range(N):
                ## trearing same level neighbors
                cells_nb_dic_neg_sl_lev_dir = cells_nb_dic_neg_sl_lev[n_dir]
                if len(cells_nb_dic_neg_sl_lev_dir)>0:
                    c, c_nb = cells_nb_dic_neg_sl_lev_dir[0]
                    n0 = cellsHier[lev][c][H_CN][0]    ## lower left corner node
                    n0_nb = cellsHier[lev][c_nb][H_CN][0]
                    dr = np.round((nodesPoints[n0]-nodesPoints[n0_nb])/dx_lev)
                    assert np.all((np.abs(dr)==0)+(np.abs(dr)==1))
                    ##set the continuity up to der_ord_max[n_dir]-iem derivative
                    cc_dic['n']['sl'][lev][n_dir] = [None]*der_ord_max[n_dir]
                    for d_ord in range(der_ord_max[n_dir]):
                        der_order = [0]*N
                        der_order[n_dir] = d_ord
                        arr_I_FL_f0_nb, inds_0_FL_f0_nb = self.PolyScaleShiftOrigin(poly_order, 
                            der_order, orders_max=polyorders_max, level=lev, dr=dr, h=None, getface0=True)
                        arr_I_FL_f0, inds_0_FL_f0 = self.PolyScaleShiftOrigin(poly_order, 
                            der_order, orders_max=polyorders_max, level=lev, dr=None, h=None, getface0=True)
                        cc_dic['n']['sl'][lev][n_dir][d_ord] = [[arr_I_FL_f0[n_dir], inds_0_FL_f0[n_dir]], \
                                                    [arr_I_FL_f0_nb[n_dir], inds_0_FL_f0_nb[n_dir]]]
                                                    
                        if self.verbose>0:
                            print('n_dir:', n_dir, 'd_ord:', d_ord, 'lev:', lev, 'dr:', dr)
                            print('arr_I_FL_f0[n_dir]:', arr_I_FL_f0[n_dir], 'inds_0_FL_f0[n_dir]:', inds_0_FL_f0[n_dir], \
                            'arr_I_FL_f0_nb[n_dir]:', arr_I_FL_f0_nb[n_dir], 'inds_0_FL_f0_nb[n_dir]:', inds_0_FL_f0_nb[n_dir], '-'*30, sep='\n')
                                                
                ## treating next level neighbors
                cc_dic['n']['nl'][lev][n_dir] = {}
                cells_nb_dic_neg_nl_lev_dir = cells_nb_dic_neg_nl_lev[n_dir]
                for conn_type in cells_nb_dic_neg_nl_lev_dir:
                    c, c_nb = cells_nb_dic_neg_nl_lev_dir[conn_type][0]
                    n0 = cellsHier[lev][c][H_CN][0]    ## lower left corner node
                    n0_nb = cellsHier[lev+1][c_nb][H_CN][0]
                    dr = np.round((nodesPoints[n0]-nodesPoints[n0_nb])/dx_levp1)
                    assert np.all((np.abs(dr)==0)+(np.abs(dr)==1)+(np.abs(dr)==2))
                    dr /= 2.0
                    h = 2.0*np.ones(N)
                    cc_dic['n']['nl'][lev][n_dir][conn_type] = [None]*der_ord_max[n_dir]
                    ##set the continuity up to der_ord_max[n_dir]-iem derivative
                    for d_ord in range(der_ord_max[n_dir]):
                        der_order = [0]*N
                        der_order[n_dir] = d_ord
                        arr_I_FL_f0_nb, inds_0_FL_f0_nb = self.PolyScaleShiftOrigin(poly_order, 
                            der_order, orders_max=polyorders_max, level=lev+1, dr=dr, h=h, getface0=True)
                        arr_I_FL_f0, inds_0_FL_f0 = self.PolyScaleShiftOrigin(poly_order, 
                            der_order, orders_max=polyorders_max, level=lev, dr=None, h=None, getface0=True)
                        cc_dic['n']['nl'][lev][n_dir][conn_type][d_ord] = [[arr_I_FL_f0[n_dir], inds_0_FL_f0[n_dir]], \
                                                    [arr_I_FL_f0_nb[n_dir], inds_0_FL_f0_nb[n_dir]]]
                ## treating previous level neighbors
                cc_dic['n']['pl'][lev][n_dir] = {}
                cells_nb_dic_neg_pl_lev_dir = cells_nb_dic_neg_pl_lev[n_dir]
                for conn_type in cells_nb_dic_neg_pl_lev_dir:
                    c, c_nb = cells_nb_dic_neg_pl_lev_dir[conn_type][0]
                    n0 = cellsHier[lev][c][H_CN][0]    ## lower left corner node
                    n0_nb = cellsHier[lev-1][c_nb][H_CN][0]
                    dr = np.round((nodesPoints[n0]-nodesPoints[n0_nb])/dx_lev)
                    assert np.all((np.abs(dr)==0)+(np.abs(dr)==1)+(np.abs(dr)==2))
                    h = 0.5*np.ones(N)
                    cc_dic['n']['pl'][lev][n_dir][conn_type] = [None]*der_ord_max[n_dir]
                    ##set the continuity up to der_ord_max[n_dir]-iem derivative
                    for d_ord in range(der_ord_max[n_dir]):
                        der_order = [0]*N
                        der_order[n_dir] = d_ord
                        arr_I_FL_f0_nb, inds_0_FL_f0_nb = self.PolyScaleShiftOrigin(poly_order, 
                            der_order, orders_max=polyorders_max, level=lev-1, dr=dr, h=h, getface0=True)
                        arr_I_FL_f0, inds_0_FL_f0 = self.PolyScaleShiftOrigin(poly_order, 
                            der_order, orders_max=polyorders_max, level=lev, dr=None, h=None, getface0=True)
                        cc_dic['n']['pl'][lev][n_dir][conn_type][d_ord] = [[arr_I_FL_f0[n_dir], inds_0_FL_f0[n_dir]], \
                                                    [arr_I_FL_f0_nb[n_dir], inds_0_FL_f0_nb[n_dir]]]

        self.cc_dic = cc_dic

        ## equation index and multiplicity of each cc 
        cc_eqindex_and_multiplicity = [None]*N
        for n_dir in range(N):
            cc_eqindex_and_multiplicity[n_dir] = [None]*der_ord_max[n_dir]
            for d_ord in range(der_ord_max[n_dir]):
                amul_FL_dir, inds_0_FL_dir = self.SetCCMultiplicityForEachPolyNode(poly_order, \
                    orders_max=polyorders_max, der_ord_max=der_ord_max, n_dir=n_dir, nd_ind=d_ord)
                cc_eqindex_and_multiplicity[n_dir][d_ord] = [amul_FL_dir, inds_0_FL_dir]
        self.cc_eqindex_and_multiplicity = cc_eqindex_and_multiplicity
        if self.verbose>0:
            print('self.cc_eqindex_and_multiplicity[n_dir][d_ord]: ', self.cc_eqindex_and_multiplicity)
        return [cc_dic, cc_eqindex_and_multiplicity]
        


    def PolyScaleShiftOrigin(self, poly_order, der_order, orders_max, level, dr=None, h=None, getface0=False):
        """ h : relative scale factor (relative to cell size)
            dr: relative shift in origin (relative to cell corner)
        """
        ## dr:shift
        N = self.N
        assert len(poly_order)==N and len(der_order)==N
        dx_lev = self.dx_levels[level]

        shape_p = self.mask_allp.shape

        arr = np.copy(self.mask_allp).astype(float)    ## arr[i][j][k] --> x^i*y^j*z^k    
        inds = np.copy(self.inds_0_allp)
        inds_0 = np.copy(self.inds_0_allp)
        
        ##take derivative 
        for i in range(N):
            if der_order[i]>0:
                a = np.ones(shape_p[i])
                for j in range(der_order[i]):
                    a *= (np.arange(shape_p[i])-j)/dx_lev[i]
                #print('i:{} a:\n'.format(i), a)
                for j in range(shape_p[i]):
                    indx = [Ellipsis]*N
                    indx[i] = j
                    arr[indx] *= a[j]
                    #print('der -- j:{}  arr:\n'.format(j), arr)
        ##roll
        for i in range(N):
            arr = np.roll(arr, -der_order[i], i)
            inds = np.roll(inds, -der_order[i], i)
            #print('rolling -- i:{}  arr:\n'.format(i), arr)

        ##mask zeros
        arr *= self.mask_allp
        inds = inds*self.mask_allp - (self.mask_allp==0)*1
        
        #print('arr:', arr, 'inds:', inds, 'inds_0:', inds_0, sep='\n')
        
        ##rescaling
        if h!=None:
            counter = np.zeros(N, dtype=int)
            while True:
                scale_factor = 1.0
                for i in range(N):
                    scale_factor *= h[i]**counter[i]
                arr[tuple(counter.tolist())] = arr[tuple(counter.tolist())]*scale_factor

                if not self._increasePolyCounterIndex(counter, poly_order):
                    break

        #print('h:', h, 'arr:', arr, 'inds:', inds, 'inds_0:', inds_0, sep='\n')
        ##shifting origin
        n_ind = len(self.inds_0_all_F)
        
        inds_F, inds_0_F, arr_F = [None]*n_ind, [None]*n_ind, [None]*n_ind
        arr_I_F = [None]*n_ind
        for i in range(n_ind):
            i_tup = self.inds_0_all_F[i][1]
            inds_F[i] = inds[i_tup]
            inds_0_F[i] = inds_0[i_tup]
            arr_F[i] = arr[i_tup]
            if inds_F[i]>=0 and arr_F[i]!=0:
                arr_I_F[i] = {inds_F[i]: arr_F[i]}
            else:
                arr_I_F[i] = {}
            
        #print('arr_F:', arr_F, 'inds_F:', inds_F, 'inds_0_F:', inds_0_F, sep='\n')
        if dr!=None:
            counter = np.zeros(self.nx.shape, dtype=int)
            while True:
                #counter_val = self.getPolyCounterValue(counter, poly_order)
                for n in range(N):
                    if dr[n]==0:
                        continue
                    c = arr[tuple(counter.tolist())]
                    if c==0:
                        continue
                    b = [binom(counter[n], k)*dr[n]**(counter[n]-k) for k in range(counter[n]+1)]

                    counter_n = np.copy(counter)
                    mask = [i for i in range(N) if i!=n]
                    counter_n[n] = 0
                    #print('b:', b, 'c:', c, sep='\n')
                    a = np.zeros(shape_p, dtype=float)
                    while True:
                        if counter_n[n]>=counter[n]:
                            break
                        a[tuple(counter_n)] = c*b[counter_n[n]]
                                    
                        if not self._increasePolyCounterIndex_Masked(counter_n, poly_order, mask):
                            break

                    a_F = [None]*n_ind
                    for i in range(n_ind):
                        i_tup = self.inds_0_all_F[i][1]
                        a_F[i] = a[i_tup]
                    
                    ind_counter = inds[tuple(counter)]
                    #print(ind_counter, arr_I_F[i])
                    assert 0<=ind_counter<n_ind
                    for i in range(n_ind):
                        if a_F[i]!=0:
                            if ind_counter in arr_I_F[i]:
                                arr_I_F[i][ind_counter] += a_F[i]
                            else:
                                arr_I_F[i][ind_counter] = a_F[i]

                    #print('dr:', dr, 'n:', n, 'counter:', counter, 'b:', b, 'c:', c)
                    #print('a:', a, 'a_F:', a_F, 'inds_F:', inds_F, 'arr_I_F:', arr_I_F, sep='\n')
                if not self._increasePolyCounterIndex(counter, poly_order):
                    break
                
        #print('arr_F:', arr_F, 'inds_F:', inds_F, 'inds_0_F:', inds_0_F, sep='\n')
        ##return the value on a given face
        if getface0:
            ##values at each 0-face (left/down face)
            inds_0_FL_f0 = [None]*N     ## _0f : 0-face
            arr_I_FL_f0 = [None]*N
            for n in range(N):
                inds_0_FL_f0[n] = []
                arr_I_FL_f0[n] = []
                a = np.zeros(shape_p, dtype=int)
                indx = [Ellipsis]*a.ndim
                indx[n] = 0
                a[indx] = 1     ## set a[n=0,:]=1
                a *= self.mask_allp
                 
                a_F = [None]*n_ind
                for i in range(n_ind):
                    i_tup = self.inds_0_all_F[i][1]
                    a_F[i] = a[i_tup]
                
                for i in range(n_ind):
                    if a_F[i]==1:
                        if len(arr_I_F[i])>0:   ##drop empties
                            inds_0_FL_f0[n].append(inds_0_F[i])
                            arr_I_FL_f0[n].append(arr_I_F[i])
                        
            #print('inds_0_FL_f0:', inds_0_FL_f0, 'arr_I_FL_f0:', arr_I_FL_f0, '-'*20, sep='\n')
            return [arr_I_FL_f0, inds_0_FL_f0]
        else:
            ##drop empties
            arr_I_FL, inds_0_FL = [], []
            for i in range(n_ind):
                if len(arr_I_F[i])>0:
                    arr_I_FL.append(arr_I_F[i])
                    inds_0_FL.append(inds_0_F[i])

            return [arr_I_FL, inds_0_FL]


    def PolyScaleShiftOrigin___ver_0(self, poly_order, der_order, orders_max, level, dr=None, h=None, getface0=False):
        """ h : relative scale factor (relative to cell size)
            dr: relative shift in origin (relative to cell corner)
        """
        ## dr:shift
        N = self.N
        assert len(poly_order)==N and len(der_order)==N
        dx_lev = self.dx_levels[level]
        poly_order_p1 = tuple([i+1 for i in poly_order])
        arr = np.ones(poly_order_p1)    ## arr[i][j][k] --> x^i*y^j*z^k    
        shape = arr.shape
        n_tot = 1 
        for i in range(N):
            n_tot *= shape[i]
        inds = np.arange(n_tot, dtype=int) 
        inds = inds.reshape(shape, order='C')
        inds_0 = np.arange(n_tot, dtype=int)  
        inds_0 = inds_0.reshape(shape, order='C')
        
        ##taking the derivatives
        for i in range(N):
            a = np.ones(shape[i])
            for j in range(der_order[i]):
                a *= (np.arange(shape[i])-j)/dx_lev[i]
            #print('i:{} a:\n'.format(i), a)
            for j in range(shape[i]):
                indx = [Ellipsis]*arr.ndim
                indx[i] = j
                arr[indx] *= a[j]
                #print('der -- j:{}  arr:\n'.format(j), arr)
        for i in range(N):
            arr = np.roll(arr, -der_order[i], i)
            inds = np.roll(inds, -der_order[i], i)
        
        #print('arr:', arr, 'inds:', inds, sep='\n')
        
        ##rescaling
        if h!=None:
            counter = np.zeros(self.nx.shape, dtype=int)
            while True:
                scale_factor = 1.0
                for i in range(N):
                    scale_factor *= h[i]**counter[i]
                arr[tuple(counter.tolist())] *= scale_factor

                if not self._increasePolyCounterIndex(counter, poly_order):
                    break

        ##shifting origin
        inds_FL = inds.reshape(n_tot, order='C').tolist()   ## _FL: flat list
        arr_FL = arr.reshape(n_tot, order='C').tolist() 
        inds_0_FL = inds_0.reshape(n_tot, order='C').tolist() 
        ## _I : the list contains indices + array elements
        arr_I_FL = [{inds_FL[i]: arr_FL[i]} for i in range(n_tot)]

        if dr!=None:
            counter = np.zeros(self.nx.shape, dtype=int)
            while True:
                counter_val = self.getPolyCounterValue(counter, poly_order)
                for n in range(N):
                    if dr[n]==0:
                        continue
                    b = [binom(counter[n], k)*dr[n]**(counter[n]-k) for k in range(counter[n]+1)]
                    c = arr[tuple(counter.tolist())]
                    counter_n = np.copy(counter)
                    mask = [i for i in range(N) if i!=n]
                    counter_n[n] = 0
                    #print('b:', b, 'c:', c, sep='\n')
                    while True:
                        if counter_n[n]>=counter[n]:
                            break
                        a = np.zeros(poly_order_p1)
                        a[tuple(counter_n)] = c*b[counter_n[n]]
                        a_F = a.reshape(n_tot, order='C')
                        for i in range(n_tot):
                            if a_F[i]!=0:
                                if inds_FL[counter_val] in arr_I_FL[i]:
                                    arr_I_FL[i][inds_FL[counter_val]] += a_F[i]
                                else:
                                    arr_I_FL[i][inds_FL[counter_val]] = a_F[i]
                                    
                        #print('dr:', dr, 'n:', n, 'counter:', counter, 'counter_n', counter_n, 'b:', b)
                        #print('a:', a, 'a_F:', a_F, 'inds_FL:', inds_FL, 'arr_I_FL:', arr_I_FL, sep='\n')

                        if not self._increasePolyCounterIndex_Masked(counter_n, poly_order, mask):
                            break
                if not self._increasePolyCounterIndex(counter, poly_order):
                    break
        ##drop forbidden indices
        inds0_to_fin_map = None
        if self.N>=2: 
            ##remove top and right blocks for 2D or higher
            inds_to_rem = self.GetIndsToRemove(poly_order, orders_max)
            inds_to_rem = inds_to_rem.reshape(n_tot, order='C')
            #_rm: removed,  _fb=forbidden
            inds_0_FL_rm = []
            arr_I_FL_rm = []
            inds_fb = []
            arr_I_FL__ = [None]*len(inds_0_FL)
            for i in range(n_tot):
                if inds_to_rem[i]==0:
                    inds_0_FL_rm.append(inds_0_FL[i])
                    arr_I_FL_rm.append(arr_I_FL[i])

                    arr_I_FL__[i] = {}
                    for ind in arr_I_FL[i]:
                        if inds_to_rem[ind]==0:
                            arr_I_FL__[i][ind] = arr_I_FL[i][ind]
                else:
                    inds_fb.append(inds_0_FL[i])
                    inds_0_FL[i] = -1
            arr_I_FL = arr_I_FL__

            arr_I_FL_rm_2 = [None]*len(arr_I_FL_rm)
            for i in range(len(arr_I_FL_rm)):
                arr_I_FL_rm_2[i] = {}
                for ind in arr_I_FL_rm[i]:
                    if ind not in inds_fb:            
                        arr_I_FL_rm_2[i][ind] = arr_I_FL_rm[i][ind]
            arr_I_FL_rm = arr_I_FL_rm_2
            if self.verbose>0:
                print('inds_0_FL_rm:', inds_0_FL_rm, 'arr_I_FL_rm:', arr_I_FL_rm, sep='\n')
            ##remap indices to 0..ind_final 
            inds0_to_fin_map = np.arange(n_tot, dtype=int)    ##original indices
            inds0_to_fin_map = inds0_to_fin_map.reshape(shape, order='C')
            counter = np.zeros(self.nx.shape, dtype=int)
            ind_map = 0
            while True:
                if inds0_to_fin_map[tuple(counter.tolist())] in inds_fb:
                    inds0_to_fin_map[tuple(counter.tolist())] = -1
                else:
                    inds0_to_fin_map[tuple(counter.tolist())] = ind_map
                    ind_map += 1
                if not self._increasePolyCounterIndex(counter, poly_order):
                    break
            inds0_to_fin_map = inds0_to_fin_map.reshape(n_tot, order='C')
            #remap
            arr_I_FL_rm2 = [None]*len(arr_I_FL_rm)
            for i in range(len(inds_0_FL_rm)):
                inds_0_FL_rm[i] = inds0_to_fin_map[inds_0_FL_rm[i]]
                arr_I_FL_rm2[i] = {}
                for ind in arr_I_FL_rm[i]:
                    arr_I_FL_rm2[i][inds0_to_fin_map[ind]] = arr_I_FL_rm[i][ind]
                    assert inds0_to_fin_map[ind]>=0
                assert inds_0_FL_rm[i]>=0
            arr_I_FL_rm = arr_I_FL_rm2
            if self.verbose>0:
                print('inds0_to_fin_map:', inds0_to_fin_map, 'inds_0_FL_rm:', inds_0_FL_rm, \
                'arr_I_FL_rm:', arr_I_FL_rm, sep='\n')
            
            arr_I_FL__ = [None]*len(arr_I_FL)
            for i in range(len(inds_0_FL)):
                if inds_0_FL[i]>=0:
                    inds_0_FL[i] = inds0_to_fin_map[inds_0_FL[i]]
                    arr_I_FL__[i] = {}
                    for ind in arr_I_FL[i]:
                        arr_I_FL__[i][inds0_to_fin_map[ind]] = arr_I_FL[i][ind]
                        assert inds0_to_fin_map[ind]>=0
            arr_I_FL = arr_I_FL__
            if self.verbose>0:
                print('inds0_to_fin_map:', inds0_to_fin_map, 'inds_0_FL:', inds_0_FL, \
                'arr_I_FL:', arr_I_FL, sep='\n')
        ##return the value on a given face
        if getface0:
            ##values at each 0-face (left/down face)
            inds_0_FL_f0 = [None]*N     ## _0f : 0-face
            arr_I_FL_f0 = [None]*N
            for n in range(N):
                inds_0_FL_f0[n] = []
                arr_I_FL_f0[n] = []
                a = np.zeros(poly_order_p1, dtype=int)
                indx = [Ellipsis]*a.ndim
                indx[n] = 0
                a[indx] = 1     ## set a[n=0,:]=1 
                a_F = a.reshape(n_tot, order='C')
                for i in range(n_tot):
                    if a_F[i]==1:
                        if inds0_to_fin_map==None:
                            inds_0_FL_f0[n].append(inds_0_FL[i])
                            arr_I_FL_f0[n].append(arr_I_FL[i])
                        elif inds_0_FL[i]>=0:
                            inds_0_FL_f0[n].append(inds_0_FL[i])
                            arr_I_FL_f0[n].append(arr_I_FL[i])
                        
            if self.verbose>0:
                print('inds_0_FL_f0:', inds_0_FL_f0, 'arr_I_FL_f0:', arr_I_FL_f0, sep='\n')
            return [arr_I_FL_f0, inds_0_FL_f0]
        elif self.N>=2:
            return [arr_I_FL_rm, inds_0_FL_rm]
        else:
            assert self.N==1
            return [arr_I_FL, inds_0_FL]



    def SetCCMultiplicityForEachPolyNode(self, poly_order, orders_max, der_ord_max, n_dir=None, nd_ind=None):
        """ in 2D and higher when setting continuity of variable and derivatives
        across cells, the number of equations will be higher than possible unknowns,
        some equations are given a multiplicity number and are averaged according
        to multiplicity number.
        """
        N = self.N
        assert len(poly_order)==N and len(der_ord_max)==N

        inds_0 = np.copy(self.inds_0_allp)
        shape_p = inds_0.shape

        a_mul = np.zeros(shape_p)    ## arr[i][j][k] --> x^i*y^j*z^k    
        
        for n in range(N):
            for d_ord in range(der_ord_max[n]):
                indx = [Ellipsis]*a_mul.ndim
                indx[n] = d_ord
                a_mul[indx] += 1
        a_mul *= self.mask_allp
        
        n_ind = len(self.inds_0_all_F)
        
        inds_0_F, a_mul_F = [], []
        for i in range(n_ind):
            i_tup = self.inds_0_all_F[i][1]
            if a_mul[i_tup]!=0:
                inds_0_F.append(inds_0[i_tup])
                a_mul_F.append(a_mul[i_tup])
                
        if n_dir!=None:
            ## return the result only for the given direction and index 
            ## i.e. for example for n_dir=2, put a_mul[i][j][k=nd_ind] and 
            ## inds_0[i][j][k=nd_ind] in a flat list and return them
            assert nd_ind!=None
            a = np.zeros(shape_p, dtype=int)
            indx = [Ellipsis]*N
            indx[n_dir] = nd_ind
            a[indx] = 1

            a *= self.mask_allp

            a_F, inds_0_F_dir = [], []
            for i in range(n_ind):
                i_tup = self.inds_0_all_F[i][1]
                if a[i_tup]!=0:
                    a_F.append(a_mul[i_tup])
                    inds_0_F_dir.append(inds_0[i_tup])

            return [a_F, inds_0_F_dir]
        else:
            return [a_mul_F, inds_0_F]


    def SetCCMultiplicityForEachPolyNode___ver_0(self, poly_order, orders_max, der_ord_max, n_dir=None, nd_ind=None):
        """ ** defective ** 
        in 2D and higher when setting continuity of variable and derivatives
        across cells, the number of equations will be higher than possible unknowns,
        some equations are given a multiplicity number and are averaged according
        to multiplicity number.
        """
        N = self.N
        assert len(poly_order)==N and len(der_ord_max)==N
        poly_order_p1 = tuple([i+1 for i in poly_order])
        a_mul = np.zeros(poly_order_p1)    ## arr[i][j][k] --> x^i*y^j*z^k    
        shape = a_mul.shape
        n_tot = 1 
        for i in range(N):
            n_tot *= shape[i]
        inds_0 = np.arange(n_tot, dtype=int)  
        inds_0 = inds_0.reshape(shape, order='C')
        
        for n in range(N):
            for d_ord in range(der_ord_max[n]):
                indx = [Ellipsis]*a_mul.ndim
                indx[n] = d_ord
                a_mul[indx] += 1
        
        a_mul_F = a_mul.reshape(n_tot, order='C')
        inds_0_F = inds_0.reshape(n_tot, order='C')
        a_mul_FL, inds_0_FL = [], []
        for i in range(n_tot):
            if a_mul_F[i]!=0:
                a_mul_FL.append(a_mul_F[i])
                inds_0_FL.append(inds_0_F[i])
                
        inds0_to_fin_map = None
        if self.N>=2:
            ##remove top and right blocks for 2D or higher
            inds_to_rem = self.GetIndsToRemove(poly_order, orders_max)
            inds_to_rem = inds_to_rem.reshape(n_tot, order='C')
            #_rm: removed,  _fb=forbidden
            inds_fb = []
            for i in range(n_tot):
                if inds_to_rem[i]==1:
                    inds_fb.append(inds_0_F[i])
            ##remap indices to 0..ind_final 
            inds0_to_fin_map = np.arange(n_tot, dtype=int)    ##original indices
            inds0_to_fin_map = inds0_to_fin_map.reshape(shape, order='C')
            counter = np.zeros(self.nx.shape, dtype=int)
            ind_map = 0
            while True:
                if inds0_to_fin_map[tuple(counter.tolist())] in inds_fb:
                    inds0_to_fin_map[tuple(counter.tolist())] = -1
                else:
                    inds0_to_fin_map[tuple(counter.tolist())] = ind_map
                    ind_map += 1
                if not self._increasePolyCounterIndex(counter, poly_order):
                    break
            inds0_to_fin_map = inds0_to_fin_map.reshape(n_tot, order='C')

            for i in range(len(inds_0_FL)):
                inds_0_FL[i] = inds0_to_fin_map[inds_0_FL[i]]

        if n_dir!=None:
            ## return the result only for the given direction and index 
            ## i.e. for example for n_dir=2, put a_mul[i][j][k=nd_ind] and 
            ## inds_0[i][j][k=nd_ind] in a flat list and return them
            assert nd_ind!=None
            a = np.zeros(poly_order_p1)
            indx = [Ellipsis]*N
            indx[n_dir] = nd_ind
            a[indx] = 1
            a_F = a.reshape(n_tot, order='C')
            a_FL, inds_0_FL_dir = [], []
            for i in range(n_tot):
                if a_F[i]!=0:
                    a_FL.append(a_mul_F[i])
                    if inds0_to_fin_map==None:
                        inds_0_FL_dir.append(inds_0_F[i])
                    else:
                        inds_0_FL_dir.append(inds0_to_fin_map[inds_0_F[i]])
            return [a_FL, inds_0_FL_dir]
        else:
            return [a_mul_FL, inds_0_FL]
        
        
    def SetCCCellMultiplicity(self):
        """ For cells with multiple neighbors on left(down..) the continuity condition
        on its left(down...) face and all neighbor cells is averaged. (to get equal
        number of equations and unknowns)
        """
        N = self.N
        cc_nbmul_dic = {'n':{'sl':[None]*N, 'nl':[None]*N, 'pl':[None]*N}}
        ## nl: left(down...) neighbor cell is smaller
        for n in range(N):
            cc_nbmul_dic['n']['sl'][n] = 1.0
            cc_nbmul_dic['n']['nl'][n] = float(2**(N-1))
            cc_nbmul_dic['n']['pl'][n] = 1.0
        self.cc_nbmul_dic = cc_nbmul_dic
        return cc_nbmul_dic


    def GetCellBoundCondEqs_scaled_masked(self, poly_order, polyorders_max, der_ord_max):
        """ 
        """
        N = self.N
        cellsHier = self.cellsHier
        n_lev = len(cellsHier)
        bc_dic = {'n':[None]*N, 'p':[None]*N}
        for i in range(N):
            bc_dic['n'][i] = [None]*n_lev
            bc_dic['p'][i] = [None]*n_lev
        for n_dir in range(N):
            for lev in range(n_lev):
                bc_dic['n'][n_dir][lev] = [None]*der_ord_max[n_dir]
                for d_ord in range(der_ord_max[n_dir]):
                    der_order = [0]*N
                    der_order[n_dir] = d_ord
                    arr_I_FL_f0, inds_0_FL_f0 = self.PolyScaleShiftOrigin(poly_order, 
                        der_order, orders_max= polyorders_max, level=lev, dr=None, h=None, getface0=True)
                    bc_dic['n'][n_dir][lev][d_ord] = [arr_I_FL_f0[n_dir], inds_0_FL_f0[n_dir]]
                bc_dic['p'][n_dir][lev] = [None]*der_ord_max[n_dir]
                for d_ord in range(der_ord_max[n_dir]):
                    der_order = [0]*N
                    der_order[n_dir] = d_ord
                    dr = np.zeros(N)
                    dr[n_dir] = 1.0
                    arr_I_FL_f0, inds_0_FL_f0 = self.PolyScaleShiftOrigin(poly_order, 
                        der_order, orders_max= polyorders_max, level=lev, dr=dr, h=None, getface0=True)
                    bc_dic['p'][n_dir][lev][d_ord] = [arr_I_FL_f0[n_dir], inds_0_FL_f0[n_dir]]
        self.bc_dic = bc_dic
        if self.verbose>0:
            print('self.bc_dic[n/p][n_dir][lev][d_ord]:', self.bc_dic, sep='\n')
        
        ## equation index and multiplicity of each bc 
        bc_eqindex_and_multiplicity = [None]*N
        for n_dir in range(N):
            bc_eqindex_and_multiplicity[n_dir] = [None]*der_ord_max[n_dir]
            for d_ord in range(der_ord_max[n_dir]):
                amul_FL_dir, inds_0_FL_dir = self.SetCCMultiplicityForEachPolyNode(poly_order, \
                    orders_max=polyorders_max, der_ord_max=der_ord_max, n_dir=n_dir, nd_ind=d_ord)
                bc_eqindex_and_multiplicity[n_dir][d_ord] = [amul_FL_dir, inds_0_FL_dir]
        self.bc_eqindex_and_multiplicity = bc_eqindex_and_multiplicity
        if self.verbose>0:
            print('self.bc_eqindex_and_multiplicity[n_dir][d_ord]:', self.bc_eqindex_and_multiplicity, '-'*30, sep='\n')
        return [bc_dic, bc_eqindex_and_multiplicity]



    def GetCellBoundPolyCoeffOnFaceTemp_scaled_masked(self, poly_order, vp_der_ord_max):
        """ vp_der_ord_max: maximum derivatives for parameters or variables
        """
        N = self.N
        cellsHier = self.cellsHier
        n_lev = len(cellsHier)
        polyorders_max = self.der_orders_max[0]
        pf_dic = {'n':[None]*N, 'p':[None]*N}   ##pf: poly at face
        for i in range(N):
            pf_dic['n'][i] = [None]*n_lev
            pf_dic['p'][i] = [None]*n_lev
        for n_dir in range(N):
            for lev in range(n_lev):
                pf_dic['n'][n_dir][lev] = {}
                pf_dic['p'][n_dir][lev] = {}
                counter = np.zeros(N, dtype=int)
                while True:
                    der_order = counter
                    arr_I_FL_f0, inds_0_FL_f0 = self.PolyScaleShiftOrigin(poly_order, 
                        der_order, orders_max= polyorders_max, level=lev, dr=None, h=None, getface0=True)
                    pf_dic['n'][n_dir][lev][tuple(der_order.tolist())] = [arr_I_FL_f0[n_dir], inds_0_FL_f0[n_dir]]

                    dr = np.zeros(N)
                    dr[n_dir] = 1.0
                    arr_I_FL_f0, inds_0_FL_f0 = self.PolyScaleShiftOrigin(poly_order, 
                        der_order, orders_max= polyorders_max, level=lev, dr=dr, h=None, getface0=True)
                    pf_dic['p'][n_dir][lev][tuple(der_order.tolist())] = [arr_I_FL_f0[n_dir], inds_0_FL_f0[n_dir]]
                    
                    if not self._increasePolyCounterIndex(counter, vp_der_ord_max):
                        break
                    
        self.pf_dic = pf_dic
        if self.verbose>0:
            print('self.pf_dic[n/p][n_dir][lev][d_ord]:', self.pf_dic, sep='\n')
        return pf_dic



    def MapBCsToICs(self, var_der_ords_max):
        """ ** RESULTS IN SINGULAR MATRICES FOR BOUNDARY CONDITIONS IN 2D OR HIGHER **
        map boundary conditions to initial conditions (to set the equation
        index for each BC)
        var_der_ords_max: maximum derivative for each variable
        """
        ##self.BCs: a list [{'dir':dir, 'face':'n'/'p', 'bc':expr}, ...]
        N = self.N
        n_vars = len(self.vars_list)
        assert len(var_der_ords_max)==n_vars
        n_conds = [0]*N    #number of bc/ic conditions
        for v_ind in range(n_vars):
            d_ord_max_v = var_der_ords_max[v_ind]
            for i in range(N):
                n_conds[i] += d_ord_max_v[i]
        n_conds_tot = sum(n_conds)
        n_bcs = len(self.BCs)
        assert n_bcs==n_conds_tot
        
        bc_inds = [None]*N
        for n in range(N):
            bc_inds[n] = [None]*n_vars
            for v_ind in range(n_vars):
                bc_inds[n][v_ind] = [None]*var_der_ords_max[v_ind][n]

        bc_max_ders = [None]*n_bcs
        for i_bc in range(n_bcs):
            bc_max_ders[i_bc] = {}
        
        ##Set max derivatives for each variable
        BC_parts, BC_rhs = self.DisintegrateBoundaryConds()
        for i_bc in range(len(self.BCs)):
            b_cond = self.BCs[i_bc]
            n_dir, expr, face = b_cond['dir'], b_cond['expr'], b_cond['face']
            for v_ind in BC_parts[i_bc]:
                eqs_v = BC_parts[i_bc][v_ind]
                for eq_v_parts in eqs_v:
                    coeff_sym, der_ord, der_vars = eq_v_parts
                    der_orders = self.ChangeDerOrderAndDerVarsToStandardForm(der_ord, der_vars)
                    ##only derivatives normal to the boundary are considered
                    if der_orders[n_dir]!=der_ord:
                        raise NotImplementedError()
                    if 'face' not in bc_max_ders[i_bc]:
                        bc_max_ders[i_bc]['face'] = face 
                    if 'dir' not in bc_max_ders[i_bc]:
                        bc_max_ders[i_bc]['dir'] = n_dir 
                    if v_ind not in bc_max_ders[i_bc]:
                        bc_max_ders[i_bc][v_ind] = der_ord 
                    elif bc_max_ders[i_bc][v_ind]<der_ord:
                        bc_max_ders[i_bc][v_ind] = der_ord

        if self.verbose>0:
            print('bc_max_ders: ', bc_max_ders)
        ##fill face=='n' available indices
        bc_indexed = [0]*len(self.BCs)
        for i_bc in range(n_bcs):
            face = bc_max_ders[i_bc]['face']
            if face=='n':
                n = bc_max_ders[i_bc]['dir']
                v_ind_max, d_max = None, None
                for v_ind in range(n_vars):
                    if v_ind in bc_max_ders[i_bc]:
                        if v_ind_max==None:
                            v_ind_max = v_ind
                            d_max = bc_max_ders[i_bc][v_ind]
                        elif bc_max_ders[i_bc][v_ind]>d_max:
                            v_ind_max = v_ind
                            d_max = bc_max_ders[i_bc][v_ind]
                assert v_ind_max!=None
                if bc_indexed[i_bc]==0 and d_max<var_der_ords_max[v_ind_max][n]:
                    if bc_inds[n][v_ind_max][d_max]==None:
                        bc_inds[n][v_ind_max][d_max] = i_bc
                        bc_indexed[i_bc] = 1
                            
        ##fill the remaining face=='n'/'p' available indices
        for i_bc in range(n_bcs):
            if bc_indexed[i_bc]==0:
                face = bc_max_ders[i_bc]['face']
                n = bc_max_ders[i_bc]['dir']
                for v_ind in range(n_vars):
                    if v_ind in bc_max_ders[i_bc]:
                        for d in range(var_der_ords_max[v_ind][n]):
                            if bc_inds[n][v_ind][d]==None:
                                bc_inds[n][v_ind][d] = i_bc
                                bc_indexed[i_bc] = 1
                        
        ##assign the remaining v_ind s to any available indices in the same direction
        for i_bc in range(n_bcs):
            if bc_indexed[i_bc]==0:
                face = bc_max_ders[i_bc]['face']
                n = bc_max_ders[i_bc]['dir']
                for v_ind in range(n_vars):
                    for d in range(var_der_ords_max[v_ind][n]):
                        if bc_inds[n][v_ind][d]==None:
                            bc_inds[n][v_ind][d] = i_bc
                            bc_indexed[i_bc] = 1
        ##check all bcs are indexed
        assert 0 not in bc_indexed

        ##map each bc index to ---> [n_dir][var][d_ord]
        bctoic_map = [None]*len(self.BCs)
        for n in range(N):
            for v_ind in range(n_vars):
                for d_ord in range(var_der_ords_max[v_ind][n]):
                    i_bc = bc_inds[n][v_ind][d_ord]
                    assert i_bc!=None
                    assert bctoic_map[i_bc]==None
                    bctoic_map[i_bc] = {'dir':n, 'v_ind':v_ind, 'd_ord':d_ord}
        assert None not in bctoic_map

        if self.verbose>0:
            print('self.BCs:', self.BCs , 'bc_inds:', bc_inds, 'bctoic_map:', bctoic_map, sep='\n')
        self.bctoic_map = bctoic_map
        return bctoic_map
        
        
    def MapBCICsToBCICs(self, var_der_ords_max):
        """
        map boundary conditions to possible boundary conditions
        var_der_ords_max: maximum derivative for each variable
        """
        ##self.BCs: a list [{'dir':dir, 'face':'n'/'p', 'bc':expr}, ...]
        N = self.N
        n_vars = len(self.vars_list)
        assert len(var_der_ords_max)==n_vars
        n_conds = [0]*N    #number of bc/ic conditions
        for v_ind in range(n_vars):
            d_ord_max_v = var_der_ords_max[v_ind]
            for i in range(N):
                n_conds[i] += d_ord_max_v[i]
        n_conds_tot = sum(n_conds)
        n_bcs = len(self.BCs)
        assert n_bcs==n_conds_tot
        
        bc_inds = {'n':[None]*N, 'p':[None]*N}
        for n in range(N):
            bc_inds['n'][n] = [None]*n_vars
            bc_inds['p'][n] = [None]*n_vars
            for v_ind in range(n_vars):
                bc_inds['n'][n][v_ind] = [None]*var_der_ords_max[v_ind][n]
                bc_inds['p'][n][v_ind] = [None]*var_der_ords_max[v_ind][n]

        bc_max_ders = [None]*n_bcs
        for i_bc in range(n_bcs):
            bc_max_ders[i_bc] = {}
        
        ##Set max derivatives for each variable
        BC_parts, BC_rhs = self.DisintegrateBoundaryConds()
        for i_bc in range(len(self.BCs)):
            b_cond = self.BCs[i_bc]
            n_dir, expr, face = b_cond['dir'], b_cond['expr'], b_cond['face']
            for v_ind in BC_parts[i_bc]:
                eqs_v = BC_parts[i_bc][v_ind]
                for eq_v_parts in eqs_v:
                    coeff_sym, der_ord, der_vars = eq_v_parts
                    der_orders = self.ChangeDerOrderAndDerVarsToStandardForm(der_ord, der_vars)
                    ##only derivatives normal to the boundary are considered
                    if der_orders[n_dir]!=der_ord:
                        raise NotImplementedError()
                    if 'face' not in bc_max_ders[i_bc]:
                        bc_max_ders[i_bc]['face'] = face 
                    if 'dir' not in bc_max_ders[i_bc]:
                        bc_max_ders[i_bc]['dir'] = n_dir 
                    if v_ind not in bc_max_ders[i_bc]:
                        bc_max_ders[i_bc][v_ind] = der_ord 
                    elif bc_max_ders[i_bc][v_ind]<der_ord:
                        bc_max_ders[i_bc][v_ind] = der_ord
        
        if self.verbose>0:
            print('bc_max_ders: ', bc_max_ders)
        ##fill face=='n'/'p' available indices
        bc_indexed = [0]*len(self.BCs)
        for i_bc in range(n_bcs):
            face = bc_max_ders[i_bc]['face']
            n = bc_max_ders[i_bc]['dir']
            v_ind_max, d_max = None, None
            for v_ind in range(n_vars):
                if v_ind in bc_max_ders[i_bc]:
                    if v_ind_max==None:
                        v_ind_max = v_ind
                        d_max = bc_max_ders[i_bc][v_ind]
                    elif bc_max_ders[i_bc][v_ind]>d_max:
                        v_ind_max = v_ind
                        d_max = bc_max_ders[i_bc][v_ind]
            assert v_ind_max!=None
            if bc_indexed[i_bc]==0 and d_max<var_der_ords_max[v_ind_max][n]:
                if bc_inds[face][n][v_ind_max][d_max]==None:
                    bc_inds[face][n][v_ind_max][d_max] = i_bc
                    bc_indexed[i_bc] = 1
                            
        ##fill the remaining face=='n'/'p' available indices
        for i_bc in range(n_bcs):
            if bc_indexed[i_bc]==0:
                face = bc_max_ders[i_bc]['face']
                n = bc_max_ders[i_bc]['dir']
                for v_ind in range(n_vars):
                    if v_ind in bc_max_ders[i_bc]:
                        for d in range(var_der_ords_max[v_ind][n]):
                            if bc_inds[face][n][v_ind][d]==None:
                                bc_inds[face][n][v_ind][d] = i_bc
                                bc_indexed[i_bc] = 1

        ##check all bcs are indexed
        assert 0 not in bc_indexed

        self.bc_inds = bc_inds
        if self.verbose>0:
            print('self.bc_inds:', self.bc_inds)
        
        bcic_to_bcic_map = [None]*n_bcs
        
        for face in ['n', 'p']:
            for n_dir in range(N):
                for v_ind in range(n_vars):
                    for d_ord in range(var_der_ords_max[v_ind][n_dir]):
                        i_bc = bc_inds[face][n_dir][v_ind][d_ord]
                        if i_bc!=None:
                            assert bcic_to_bcic_map[i_bc]==None
                            bcic_to_bcic_map[i_bc] = {'face':face, 'dir':n_dir, 'v_ind':v_ind, 'd_ord':d_ord}
        self.bcic_to_bcic_map = bcic_to_bcic_map
        if self.verbose>0:
            print('self.bcic_to_bcic_map:', self.bcic_to_bcic_map, sep='\n')
        return        
        


    def MapBCsToICsForEachBoundaryCell(self, poly_order, var_der_ords_max):
        ##TODO: maxwell type systems would be probably not handled.. to check
        ##index ic cells
        N = self.N
        cellsHier = self.cellsHier
        n_lev = len(cellsHier)
        n_vars = len(self.vars_list)
        poly_order_p1 = tuple([i+1 for i in poly_order])
        for v_ind in range(n_vars):
            assert var_der_ords_max[v_ind]==var_der_ords_max[0]

        cell_n_tot__vars = self.GetNumOfUnknownsInEachCellForEachVar()
        cell_n_tot = sum(cell_n_tot__vars)
        cell_n_tot_vi = cell_n_tot__vars[0]
        
        cells_bc = [None]*n_lev
        for lev in range(n_lev):
            cells_bc[lev] = {}
            cells_bc_lev = cells_bc[lev]
            for n_dir in range(N):
                bc_conns_lev_n_dir = self.boundaryCellConnections[lev]['n'][n_dir]        
                cellsHier_lev = cellsHier[lev]
                for c_ind in bc_conns_lev_n_dir:
                    if cellsHier_lev[c_ind][H_CI]>=0:
                        if c_ind in cells_bc_lev:
                            cells_bc_lev[c_ind][0].append((n_dir, 'n'))
                        else:
                            cells_bc_elem = [[(n_dir, 'n')], [None]*n_vars, [None]*n_vars, [None]*n_vars]
                            for v in range(n_vars):
                                cells_bc_elem[1][v] = np.zeros(poly_order_p1, dtype=int)    ##mask
                                cells_bc_elem[2][v] = np.zeros(poly_order_p1, dtype=int)    ##multiplicity
                                cells_bc_elem[3][v] = -np.ones(poly_order_p1, dtype=int)    ##eq inds
                            cells_bc_lev[c_ind] = cells_bc_elem
                bc_conns_lev_p_dir = self.boundaryCellConnections[lev]['p'][n_dir]        
                for c_ind in bc_conns_lev_p_dir:
                    if cellsHier_lev[c_ind][H_CI]>=0:
                        if c_ind in cells_bc_lev:
                            cells_bc_lev[c_ind][0].append((n_dir, 'p'))
                        else:
                            cells_bc_elem = [[(n_dir, 'p')], [None]*n_vars, [None]*n_vars, [None]*n_vars]
                            for v in range(n_vars):
                                cells_bc_elem[1][v] = np.zeros(poly_order_p1, dtype=int)
                                cells_bc_elem[2][v] = np.zeros(poly_order_p1, dtype=int)
                                cells_bc_elem[3][v] = -np.ones(poly_order_p1, dtype=int)
                            cells_bc_lev[c_ind] = cells_bc_elem
                    
        ##set multiplicities
        if not self.bc_inds:
            self.MapBCICsToBCICs(var_der_ords_max)
            
        bc_inds = self.bc_inds  #bc_inds[face][n_dir][v_ind][d_ord] 
                       
        for lev in range(n_lev):
            for c_ind, dir_mat in cells_bc[lev].items():
                dirs, masks, mults, eqinds = dir_mat
                for dir_face in dirs:
                    n_dir, face = dir_face
                    for v_ind in range(n_vars):
                        mul_v = mults[v_ind]
                        msk_v = masks[v_ind]
                        for d_ord in range(var_der_ords_max[v_ind][n_dir]):
                            i_bc = bc_inds[face][n_dir][v_ind][d_ord]
                            if i_bc!=None:
                                if face=='n':
                                    indx = [Ellipsis]*N
                                    indx[n_dir] = d_ord
                                    msk_v[indx] += 1
                                    mul_v[indx] += 1
                                elif face=='p':
                                    indx = [Ellipsis]*N
                                    indx[n_dir] = poly_order[n_dir]-d_ord
                                    msk_v[indx] += 1
                                    mul_v[indx] += 1
                                else:
                                    raise ValueError()

        ##set eq indices for indices intersecting self.inds_0_all_F
        inds_0_all_F = self.inds_0_all_F
        inds_0_F = [inds_0_all_F[i][1] for i in range(len(inds_0_all_F)) if self.mask_x_F[i]==0]
        for lev in range(n_lev):
            for c_ind, dir_mat in cells_bc[lev].items():
                dirs, masks, mults, eqinds = dir_mat
                cells_lev = cellsHier[lev]
                ind_lev_st = self.N_indexed_acc[lev]
                for v_ind in range(n_vars):
                    masks_v = masks[v_ind]
                    eqinds_v = eqinds[v_ind]
                    for ind in inds_0_F:
                        if masks_v[ind]>0:
                            eqinds_v[ind] = (ind_lev_st + cells_lev[c_ind][H_CI])*cell_n_tot \
                                        + v_ind*cell_n_tot_vi + self.inds_0_allp[ind]
                        
        
        ##add cc multiplicities to get total multiplities        
        for lev in range(n_lev):
            for c_ind, dir_mat in cells_bc[lev].items():
                dirs, masks, mults, eqinds = dir_mat
                dir_neg_cc = [True]*N
                for dir_face in dirs:
                    n_dir, face = dir_face
                    if face=='n':
                        dir_neg_cc[n_dir] = False
                for n_dir in range(N):
                    if dir_neg_cc[n_dir]==True:
                        for v_ind in range(n_vars):
                            mul_v = mults[v_ind]
                            for d_ord in range(var_der_ords_max[v_ind][n_dir]):
                                indx = [Ellipsis]*N
                                indx[n_dir] = d_ord
                                mul_v[indx] += 1
        

        ##set empties
        inds_empt = [None]*n_vars
        for v_ind in range(n_vars):
            inds_empt[v_ind] = []
        for lev in range(n_lev):
            for c_ind, dir_mat in cells_bc[lev].items():
                dirs, masks, mults, eqinds = dir_mat
                cells_lev = cellsHier[lev]
                ind_lev_st = self.N_indexed_acc[lev]
                for v_ind in range(n_vars):
                    mults_v = mults[v_ind]
                    eqinds_v = eqinds[v_ind]
                    #print('inds_0_F:', inds_0_F, 'mults_v:', mults_v, sep='\n')
                    for ind in inds_0_F:
                        if mults_v[ind]==0:
                            eqinds_empt = (ind_lev_st + cells_lev[c_ind][H_CI])*cell_n_tot \
                                        + v_ind*cell_n_tot_vi + self.inds_0_allp[ind]
                            inds_empt[v_ind].append(eqinds_empt)
        
        if self.verbose>0:
            print('inds_empt:\n', inds_empt)

        ##fill empties with the rest of bcs
        ##TODO: replace counter in the next block.. for speed.
        inds_empt_last = [0]*n_vars
        for lev in range(n_lev):
            for c_ind, dir_mat in cells_bc[lev].items():
                dirs, masks, mults, eqinds = dir_mat
                ind_lev_st = self.N_indexed_acc[lev]
                for v_ind in range(n_vars):
                    masks_v = masks[v_ind]
                    eqinds_v = eqinds[v_ind]
                    counter = np.zeros(N, dtype=int)
                    while True:
                        ctr_ind = tuple(counter.tolist())
                        #print('masks_v:', masks_v, 'eqinds_v:', eqinds_v, 'ctr_ind:', ctr_ind, sep='\n')
                        if masks_v[ctr_ind]>0 and eqinds_v[ctr_ind]<0:
                            eqinds_v[ctr_ind] = inds_empt[v_ind][inds_empt_last[v_ind]]
                            inds_empt_last[v_ind] += 1
                            
                        if not self._increasePolyCounterIndex(counter, poly_order):
                            break
        
        for v_ind in range(n_vars):
            #print('len(inds_empt[v_ind]):', len(inds_empt[v_ind]), 'inds_empt_last[v_ind]:', inds_empt_last[v_ind])
            assert inds_empt_last[v_ind] == len(inds_empt[v_ind])
            
        ##set cc and bc multiplicity and indices for each cell        
        mask_dir_der = [None]*n_vars
        for v_ind in range(n_vars):
            mask_dir_der[v_ind] = {'n':[None]*N, 'p':[None]*N}
            for n_dir in range(N):
                d_ord_max_v_ndir = var_der_ords_max[v_ind][n_dir]
                mask_dir_der[v_ind]['n'][n_dir] = [None]*d_ord_max_v_ndir
                mask_dir_der[v_ind]['p'][n_dir] = [None]*d_ord_max_v_ndir
                for d_ord in range(d_ord_max_v_ndir):
                    mask_dir_der[v_ind]['n'][n_dir][d_ord] = []
                    mask_dir_der[v_ind]['p'][n_dir][d_ord] = []
                    mask_n = np.zeros(poly_order_p1, dtype=int)
                    indx = [Ellipsis]*N
                    indx[n_dir] = d_ord
                    mask_n[indx] = 1
                    mask_p = np.zeros(poly_order_p1, dtype=int)
                    indx = [Ellipsis]*N
                    indx[n_dir] = poly_order[n_dir] - d_ord
                    mask_p[indx] = 1
                    counter = np.zeros(N, dtype=int)
                    while True:
                        ctr_ind = tuple(counter.tolist())
                        if mask_n[ctr_ind]>0:
                            mask_dir_der[v_ind]['n'][n_dir][d_ord].append(ctr_ind)
                        if mask_p[ctr_ind]>0:
                            mask_dir_der[v_ind]['p'][n_dir][d_ord].append(ctr_ind)
                        if not self._increasePolyCounterIndex(counter, poly_order):
                            break
                
        cellsBC = [None]*n_vars
        for v_ind in range(n_vars):
            cellsBC[v_ind] = [None]*n_lev
            for lev in range(n_lev):
                cellsBC[v_ind][lev] = {'n':[None]*N, 'p':[None]*N}
                cellsBC_vlev = cellsBC[v_ind][lev]
                for n_dir in range(N):
                    cellsBC_vlev['n'][n_dir] = [None]*var_der_ords_max[v_ind][n_dir]
                    cellsBC_vlev['p'][n_dir] = [None]*var_der_ords_max[v_ind][n_dir]
                    for d_ord in range(var_der_ords_max[v_ind][n_dir]):
                        cellsBC_vlev['n'][n_dir][d_ord] = {}
                        cellsBC_vlev['p'][n_dir][d_ord] = {}
                for c_ind, dir_mat in cells_bc[lev].items():
                    dirs, masks, mults, eqinds = dir_mat
                    #print('c_ind', c_ind, 'dirs:', dirs, 'masks:', masks, 'mults:', mults, 'eqinds:', eqinds, sep='\n')
                    for n_dir, face in dirs:
                        for d_ord in range(var_der_ords_max[v_ind][n_dir]):
                            inds = mask_dir_der[v_ind][face][n_dir][d_ord]
                            inds_n = mask_dir_der[v_ind]['n'][n_dir][d_ord] ##neg face 
                            mults_F, eqinds_F, inds_0 = [None]*len(inds), [None]*len(inds), [None]*len(inds)
                            for i in range(len(inds)):
                                mults_F[i] = mults[v_ind][inds[i]]
                                eqinds_F[i] = eqinds[v_ind][inds[i]]
                                inds_0[i] = self.inds_0_allp[inds_n[i]]
                            assert c_ind not in cellsBC_vlev[face][n_dir][d_ord]
                            cellsBC_vlev[face][n_dir][d_ord][c_ind] = [mults_F, eqinds_F, inds_0]
        
        self.cellsBC = cellsBC
        if self.verbose>0:
            print('self.cellsBC:', self.cellsBC)

        cellsCC = [None]*n_vars
        for v_ind in range(n_vars):
            cellsCC[v_ind] = [None]*n_lev
            for lev in range(n_lev):
                cellsCC[v_ind][lev] = [None]*N
                cellsCC_vlev = cellsCC[v_ind][lev]
                for n_dir in range(N):
                    cellsCC_vlev[n_dir] = [None]*var_der_ords_max[v_ind][n_dir]
                    for d_ord in range(var_der_ords_max[v_ind][n_dir]):
                        cellsCC_vlev[n_dir][d_ord] = {}
                    
                for c_ind, dir_mat in cells_bc[lev].items():
                    dirs, masks, mults, eqinds = dir_mat

                    dir_neg_cc = [True]*N
                    for n_dir, face in dirs:
                        if face=='n':
                            dir_neg_cc[n_dir] = False
                    for n_dir in range(N):
                        if dir_neg_cc[n_dir]==True:
                            for d_ord in range(var_der_ords_max[v_ind][n_dir]):
                                inds = mask_dir_der[v_ind]['n'][n_dir][d_ord]
                                mults_F, inds_0 = [None]*len(inds), [None]*len(inds)
                                for i in range(len(inds)):
                                    mults_F[i] = mults[v_ind][inds[i]]
                                    inds_0[i] = self.inds_0_allp[inds[i]]
                                assert c_ind not in cellsCC_vlev[n_dir][d_ord]
                                cellsCC_vlev[n_dir][d_ord][c_ind] = [mults_F, inds_0]
        
        self.cellsCC = cellsCC
        if self.verbose>0:
            print('self.cellsCC:', self.cellsCC)
        return




### --------------------------------------------------------------------------
### --------------------------------------------------------------------------
from sympy import Mul, Add, Derivative, Integer

##TODO: to be extended to non-linear equations

def sym_find_node_with_symbol(node, x):
    """ finds the forst node whose node[1] component is the given symbol x and returns
        it starts from the given node and walks down the tree
    """
    if node[1]==x:
        return node
    else:
        for arg in node[3]:
            subnode_x = sym_find_node_with_symbol(arg, x)
            if subnode_x!=False:
                return subnode_x
    return False


def sym_find_mult_coeff_DerSymb(node, var):
    """ 
    Finds the multiplicative coefficients of a given variable in an expression 
    linear in var.
    var is a Function. Its argumet may contain Symbol(n_str) or Symbol(n_str+'_1')..
    """
    if isinstance(node, list):
        if node[1]==var:
            var_coeff = 1
            var_ind = node[4]
            parent = node[0]
            der_order = 0   #derivative order
            der_vars = []
            while True:
                if parent==None:
                    break
                if not(parent[2]==Mul or parent[2]==Add or parent[2]==Derivative\
                        or parent[2]==Sum):
                    raise NotImplementedError("only Mul, Add, Derivative are treated")
                if parent[2]==Derivative:
                    if der_order==0:
                        der_order = len(parent[3])-1
                        assert der_order>=1
                        if parent[1].args[0]!=var:
                            raise NotImplementedError("function must be the only operand of derivative")
                        der_vars = [parent[3][k][1] for k in range(1, len(parent[3]))]
                    else:
                        raise NotImplementedError("more that one derivation operator")
                if parent[2]==Mul:
                    mul_args_parent = [parent[3][k][1] for k in range(len(parent[3])) if k!=var_ind]
                    var_coeff *= Mul(*tuple(mul_args_parent))
                var_ind = parent[4]
                parent = parent[0]
            return [var_coeff, der_order, der_vars]
        return False
    else:
        expr_tree = symExpr_generate_tree(node)
        node_var = sym_find_node_with_symbol(expr_tree, var)
        if node_var!=False:
            return sym_find_mult_coeff_DerSymb(node_var, var)
        return False


def sym_getcoeff_setzero(node, var):
    """ var is a Symbol.
    """
    node_var = sym_find_node_with_symbol(node, var)
    if node_var!=False:
         var_coeff = sym_find_mult_coeff_DerSymb(node_var, var)
         if node_var[0]!=None:
            if node_var[0][2] == Derivative:
                 node_var[0][1] = Integer(0)
                 symExpr_update_tree(node_var[0])
                 return var_coeff
         node_var[1] = Integer(0)
         symExpr_update_tree(node_var)
         return var_coeff
    return False


             
                
                    
##-----------------------------------------------------------------------
##---------------------- GUI  ------------

from threading import Thread
import queue
from tkinter import *
import time
from random import randint


import wx
import sys
from wx import glcanvas
from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *


class GUIMakerND(Thread):
    
    def __init__(self, W_3D=np.array([1.0, 1.0, 1.0]), que=None):
        self.queue = que
        self.app = None
        self.W_3D = W_3D
        
        Thread.__init__(self)
        #self.daemon = True
        self.start()
        return
        
    def SetCommChannel(self, que):
        self.queue = que
        if self.app!=None:
            self.app.SetQueue(self.queue)
        return

    def run(self):
        self.app = InteractiveGLApp(W_3D=self.W_3D)
        self.app.SetQueue(self.queue)
        #print('GUI queue: ', self.queue)

        self.app.SetCommandLoop()
        self.app.MainLoop()
        self.app.Destroy()

        self.app.deleteGUIItems()
        del self.app


    def attachRectGrid(self, rg):
        if self.app!=None:
            self.app.attachRectGrid(rg)


    def MarkRectGridCell(self, lev, c_ind, color):
        if self.app!=None:
            self.app.MarkRectGridCell(lev, c_ind, color)


class InteractiveGLApp(wx.App):
    
    def __init__(self, W_3D=np.array([1.0, 1.0, 1.0]), que=None):
        wx.App.__init__(self, redirect=False)
        self.W_3D = W_3D        # world dimensions
        self.dW_3D = W_3D/10.0
        self.V_3D = np.array([0.0, 0.0, 0.0])        # View origin
        self.V_phi = 20.0       # view angle phi
        self.V_theta = 30.0     # view angle theta
        
        self.SetFrameDimensions()
        self.queue = que
        self.frame.queue = self.queue
        
        self.modif_state = "select"
        
        self.rectGrid = None
        
        
    def SetQueue(self, que):
        self.queue = que
        self.frame.queue = self.queue
        #print('app queue: ', self.queue)

    def OnInit(self, W_3D=np.array([1.0, 1.0, 1.0])):
        frame = GLFrame(None, -1, "ND Space", pos=(0,0),
                        style=wx.DEFAULT_FRAME_STYLE, name="ND Space Program")
        frame.CreateStatusBar()

        ##------ menubar
        menuBar = wx.MenuBar()
        menu = wx.Menu()
        item = menu.Append(wx.ID_EXIT, "E&xit\tCtrl-Q", "Exit demo")
        self.Bind(wx.EVT_MENU, self.OnExitApp, item)
        menuBar.Append(menu, "&File")
        
        frame.SetMenuBar(menuBar)

        ##----- toolbar
        toolbar = wx.ToolBar(frame, id=-1)
        toolbar.SetToolBitmapSize( (17, 17) ) # Square spacer equired for non-standard size buttons on MSW
        
        tool_save = toolbar.AddTool(wx.ID_ANY, 'Save', wx.Bitmap('./images/icons/save-16x16.png'))
        toolbar.AddSeparator()
        tool_select = toolbar.AddTool(wx.ID_ANY, 'Select', wx.Bitmap('./images/icons/select-16x16.png'), kind=wx.ITEM_RADIO)
        tool_move = toolbar.AddTool(wx.ID_ANY, 'Move', wx.Bitmap('./images/icons/move-16x16.png'), kind=wx.ITEM_RADIO)
        tool_rotate = toolbar.AddTool(wx.ID_ANY, 'Rotate', wx.Bitmap('./images/icons/rotate-16x16.png'), kind=wx.ITEM_RADIO)
        tool_zoom = toolbar.AddTool(wx.ID_ANY, 'Zoom', wx.Bitmap('./images/icons/zoom-16x16.png'), kind=wx.ITEM_RADIO)
        toolbar.AddSeparator()
        tool_toggle_grid = toolbar.AddTool(wx.ID_ANY, 'gridOn', wx.Bitmap('./images/icons/grid-16x16.png'), kind=wx.ITEM_CHECK)
        tool_toggle_axis = toolbar.AddTool(wx.ID_ANY, 'gridOn', wx.Bitmap('./images/icons/axis-16x16.png'), kind=wx.ITEM_CHECK)
        
        
        toolbar.Realize()
        
        frame.SetToolBar(toolbar)
        
        frame.Bind(wx.EVT_TOOL, self.OnToolSave, tool_save)
        frame.Bind(wx.EVT_TOOL, self.OnToolSelect, tool_select)
        frame.Bind(wx.EVT_TOOL, self.OnToolMove, tool_move)
        frame.Bind(wx.EVT_TOOL, self.OnToolRotate, tool_rotate)
        frame.Bind(wx.EVT_TOOL, self.OnToolZoom, tool_zoom)
        frame.Bind(wx.EVT_TOOL, self.OnToolToggleGrid, tool_toggle_grid)
        frame.Bind(wx.EVT_TOOL, self.OnToolToggleAxis, tool_toggle_axis)
               
        
        ##-------
        
        frame.Show(True)
        frame.Bind(wx.EVT_CLOSE, self.OnCloseFrame)

        
        self.SetTopWindow(frame)
        self.frame = frame
        
        self.quit = False
        return True
        
    def OnExitApp(self, evt):
        self.frame.Close(True)
        self.quit = True
        print('OnExitApp')

    def OnCloseFrame(self, evt):
        if hasattr(self, "window") and hasattr(self.window, "ShutdownDemo"):
            self.window.ShutdownDemo()
        evt.Skip()
        self.quit = True
        print('OnCloseFrame')
        
    def OnToolSave(self, evt):
        self.frame.LogText('Save!')
        evt.Skip()
        
    def OnToolSelect(self, evt):
        self.frame.selstate = SelStates.select
        self.frame.LogText('Select! '+repr(self.frame.selstate))
        evt.Skip()

    def OnToolMove(self, evt):
        self.frame.selstate = SelStates.move
        self.frame.LogText('Move! '+repr(self.frame.selstate))
        evt.Skip()

    def OnToolRotate(self, evt):
        self.frame.selstate = SelStates.rotate
        self.frame.LogText('Rotate! '+repr(self.frame.selstate))
        evt.Skip()

    def OnToolZoom(self, evt):
        self.frame.selstate = SelStates.zoom
        self.frame.LogText('Zoom! '+repr(self.frame.selstate))
        evt.Skip()
        
    def OnToolToggleGrid(self, evt):
        self.frame.gridOn = not self.frame.gridOn
        self.frame.UpdateCanvas()
        evt.Skip()

    def OnToolToggleAxis(self, evt):
        self.frame.axisOn = not self.frame.axisOn
        self.frame.UpdateCanvas()
        evt.Skip()

    def SetViewOrigin(self, V_3D):
        self.V_3D = V_3D
        self.SetFrameDimensions()
        
    def SetViewAngles(self, theta, phi):
        self.V_theta = theta
        self.V_phi = phi
        self.SetFrameDimensions()
        
    def SetGuideGridSize(self, dW):
        self.dW_3D = dw
        self.SetFrameDimensions()
        
    def SetFrameDimensions(self):
        self.frame.UpdateDimensions(W_3D=self.W_3D, dW_3D=self.dW_3D, 
            V_3D=self.V_3D, V_theta=self.V_theta, V_phi=self.V_phi)
        
    def deleteGUIItems(self):
        ## deletes window components (frame,...)
        self.quit = True
        del self.frame
        
    def SetCommandLoop(self, time_ms=1000):
        self.time_ms = time_ms
        self.autoCaller = wx.CallLater(time_ms, self.CommandLoop)
        self.autoCaller.Start()
        
        
    def CommandLoop(self):
        if not self.quit:
            if self.queue != None:
                #print('queue: ', self.queue)
                if not self.queue.empty():
                    comm = self.queue.get()
                    self.ProcessCommand(comm)
                    self.queue.task_done()
            self.autoCaller.Restart()
        return

    def ProcessCommand(self, comm):
        if comm[0]=='deleteAllElems':
            self.frame.elements = []
        elif comm[0]=='add elem':
            if comm[1] == 'cube':
                elem = {'name':'cube', 'v_arr': comm[2]}
                self.frame.elements.append(elem)
            elif comm[1] == 'triangle':
                elem = {'name':'triangle', 'v_arr': comm[2]}
                self.frame.elements.append(elem)
            elif comm[1] == 'quad':
                elem = {'name':'quad', 'v_arr': comm[2]}
                self.frame.elements.append(elem)
            self.frame.UpdateCanvas()
        elif comm[0]=="logtext":
            text = comm[1]
            self.frame.LogText(text)
            #print("test command received!")
        return
        
    def addElement(self, elem):
        self.frame.elements.append(elem)

    def attachRectGrid(self, rg):
        self.frame.rectGrid = rg
        #self.frame.attachRectGrid(rg)

    def MarkRectGridCell(self, lev, c_ind, color):
        self.frame.MarkRectGridCell(lev, c_ind, color)


class SelStates(Enum): 
    select = 1
    move = 2
    rotate = 3
    zoom = 4
    
    
class GLFrame(wx.Frame):

    def __init__(self, parent, id, title, pos=wx.DefaultPosition,
                 size=wx.DefaultSize, style=wx.DEFAULT_FRAME_STYLE,
                 name='frame', W_3D=np.array([1.0, 1.0, 1.0])):

        self.W_3D = W_3D          # world dimensions
        self.dW_3D = W_3D/10.0
        self.V_3D = np.array([0.0, 0.0, 0.0])        # View origin
        self.V_phi = 20.0       # view angle phi
        self.V_theta = 30.0     # view angle theta

        self.queue = None       ## can also be used to send commands to the parent window

        self.selstate = SelStates.select
        self.move_started = False
        self.rotate_started = False
        self.zoom_started = False
        
        self.gridOn = True
        self.axisOn = True
        
        self.rectGrid = None
        self.rectGridOn = True
        self.rectGridLevels = None
        self.RGCellsMarked = None

        self.SetInitPerspParams()

        #
        # Forcing a specific style on the window.
        #   Should this include styles passed?
        style = wx.DEFAULT_FRAME_STYLE | wx.NO_FULL_REPAINT_ON_RESIZE

        super(GLFrame, self).__init__(parent, id, title, pos, size, style, name)

        self.GLinitialized = False
        attribList = (glcanvas.WX_GL_RGBA, # RGBA
                      glcanvas.WX_GL_DOUBLEBUFFER, # Double Buffered
                      glcanvas.WX_GL_DEPTH_SIZE, 24) # 24 bit

        #
        # Create the canvas
        self.canvas = glcanvas.GLCanvas(self, attribList=attribList)
        self.context = glcanvas.GLContext(self.canvas)
        
        # create cursors
        self.SetupCursers()
        
        #
        # Set the event handlers.
        self.canvas.Bind(wx.EVT_ERASE_BACKGROUND, self.processEraseBackgroundEvent)
        self.canvas.Bind(wx.EVT_SIZE, self.processSizeEvent)
        self.canvas.Bind(wx.EVT_PAINT, self.processPaintEvent)
        self.canvas.Bind(wx.EVT_LEFT_DOWN, self.OnCanvasLeftMouseDown)
        self.canvas.Bind(wx.EVT_LEFT_UP, self.OnCanvasLeftMouseUp)
        self.canvas.Bind(wx.EVT_MOTION, self.OnCanvasMouseMove)
        
        ## Sizer
        box = wx.BoxSizer(wx.VERTICAL)
        #box.Add((20, 30))

        self.canvas.SetMinSize((500, 400))
        box.Add(self.canvas, 1, wx.ALIGN_CENTER|wx.ALL|wx.EXPAND, 2)

        ##--- textbox
        self.textbox = wx.TextCtrl(self, size=(300, 50), style=wx.TE_MULTILINE)        
        box.Add(self.textbox, 0, wx.ALIGN_CENTER|wx.BOTTOM|wx.EXPAND, 5)

        self.SetAutoLayout(True)
        self.SetSizer(box)
                
                
        ##--- guide grid
        el_guidegrid = {'name':'guidegrid'}
        self.elements = [el_guidegrid]
        
                
    
    def LogText(self, text):
        self.textbox.AppendText(text)
    
    #
    # Canvas Proxy Methods

    def GetGLExtents(self):
        """Get the extents of the OpenGL canvas."""
        return self.canvas.GetClientSize()

    def SwapBuffers(self):
        """Swap the OpenGL buffers."""
        self.canvas.SwapBuffers()

    #
    # wxPython Window Handlers

    def processEraseBackgroundEvent(self, event):
        """Process the erase background event."""
        pass # Do nothing, to avoid flashing on MSWin

    def processSizeEvent(self, event):
        #size_frame = self.GetSize()
        #X_frame, Y_frame = size_frame
        ##---
        #self.textbox.SetWidth(X_frame)
        
        """Process the resize event."""
        size_GL = self.GetGLExtents()
        #self.canvas.SetCurrent(self.context)
        self.OnReshape(size_GL.width, size_GL.height)
        self.canvas.Refresh(False)
        self.canvas.Update()
        event.Skip()

    def processPaintEvent(self, event):
        """Process the drawing event."""
        #self.Show()
        self.canvas.SetCurrent(self.context)

        # This is a 'perfect' time to initialize OpenGL ... only if we need to
        if not self.GLinitialized:
            self.OnInitGL()
            self.GLinitialized = True

        self.OnDraw()
        self.canvas.Refresh(False)
        self.canvas.Update()
        event.Skip()
        
    def OnCanvasLeftMouseDown(self, event):
        size_GL = self.GetGLExtents()
        x = event.GetX()
        y = size_GL.height - event.GetY()
        self.LogText('Canvas: Left mouse down! (x={}, y={})'.format(x, y))
        
        if self.selstate==SelStates.select:
            pass
        elif self.selstate==SelStates.move:
            self.x_sel = x
            self.y_sel = y
            self.move_started = True
            self.ChangeCurser('hand')
        elif self.selstate==SelStates.rotate:
            self.x_sel = x
            self.y_sel = y
            self.rotate_started = True
        elif self.selstate==SelStates.zoom:
            self.x_sel = x
            self.y_sel = y
            self.zoom_started = True
        
        event.Skip()

    def OnCanvasLeftMouseUp(self, event):
        if self.selstate==SelStates.select:
            pass
        elif self.selstate==SelStates.move:
            self.move_started = False
            self.UpdateCanvas()
            self.ChangeCurser('arrow')
        elif self.selstate==SelStates.rotate:
            self.rotate_started = False
            self.UpdateCanvas()
        elif self.selstate==SelStates.zoom:
            self.zoom_started = False
            self.UpdateCanvas()
        
        event.Skip()

    def OnCanvasMouseMove(self, event):
        size_GL = self.GetGLExtents()
        x = event.GetX()
        y = size_GL.height - event.GetY()

        if self.selstate==SelStates.select:
            pass
        elif self.selstate==SelStates.move:
            if self.move_started:
                dx = x - self.x_sel
                dy = y - self.y_sel
                X, Y, Z = self.W_3D[0], self.W_3D[1], self.W_3D[2]
                X_ratio = float(X)/size_GL.width
                Y_ratio = float(Y)/size_GL.height
                self.P_Tx += dx * X_ratio
                self.P_Ty += dy * Y_ratio 
                #self.LogText('dx={}, dy={} '.format(dx, dy))
                #self.queue.put(['logtext', 'dx:{} dy:{} '.format(dx, dy)])
                self.SetPerspective(size_GL.width, size_GL.height)
                self.UpdateCanvas()
                self.x_sel = x
                self.y_sel = y
        elif self.selstate==SelStates.rotate:
            if self.rotate_started:
                dx = x - self.x_sel
                dy = y - self.y_sel
                X, Y, Z = self.W_3D[0], self.W_3D[1], self.W_3D[2]
                X_ratio = float(X)/size_GL.width*25
                Y_ratio = float(Y)/size_GL.height*25
                self.P_Rz += dx * X_ratio
                self.P_Rx -= dy * Y_ratio 
                self.SetPerspective(size_GL.width, size_GL.height)
                self.UpdateCanvas()
                self.x_sel = x
                self.y_sel = y
        elif self.selstate==SelStates.zoom:
            if self.zoom_started:
                dx = x - self.x_sel
                dy = y - self.y_sel
                X, Y, Z = self.W_3D[0], self.W_3D[1], self.W_3D[2]
                #X_ratio = float(X)/size_GL.width*25
                Y_ratio = float(Y)/size_GL.height*15
                self.P_angle -= dy * Y_ratio
                if self.P_angle<10.0:
                    self.P_angle = 10.0
                if self.P_angle>120.0:
                    self.P_angle = 120.0
                self.SetPerspective(size_GL.width, size_GL.height)
                self.UpdateCanvas()
                self.x_sel = x
                self.y_sel = y
        
        event.Skip()

    def SetupCursers(self, model=None):
        if model==None or model=='arrow':
            self.cursor_arrow = wx.Cursor(wx.CURSOR_ARROW)
        if model==None or model=='hand':
            image = wx.Image(r'./images/cursors/hand-16x16.png', wx.BITMAP_TYPE_PNG)
            self.cursor_hand = wx.Cursor(image)

    def ChangeCurser(self, model='arrow'):
        if model=='arrow':
            self.canvas.SetCursor(self.cursor_arrow)
        elif model=='hand':
            self.canvas.SetCursor(self.cursor_hand)

    #
    # GLFrame OpenGL Event Handlers
    
    def UpdateCanvas(self):
        self.canvas.Refresh(True)
        self.canvas.Update()
        

    def OnInitGL(self):
        """Initialize OpenGL for use in the window."""
        glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH)
        glEnable(GL_DEPTH_TEST)
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glClearColor(1, 1, 1, 1)

    def OnReshape(self, width, height):
        """Reshape the OpenGL viewport based on the dimensions of the window."""
        glViewport(0, 0, width, height)

        self.SetPerspective(width, height)
        
        
    def SetInitPerspParams(self):
        X, Y, Z = self.W_3D[0], self.W_3D[1], self.W_3D[2]
        XYZ_Max = max(X, Y, Z)
        XYZ_Min = min(X, Y, Z)

        self.P_angle = 60.0
        self.P_near = 0.1*XYZ_Min
        self.P_far = 6.0*XYZ_Max
        
        self.P_Tx = 0.0
        self.P_Ty = 0.0
        self.P_Tz = -2.0*XYZ_Max
        
        self.P_Rx = -60.0
        self.P_Ry = 0.0
        self.P_Rz = -20.0
                

    def SetPerspective(self, width, height):

        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        gluPerspective(self.P_angle, float(width)/float(height), self.P_near, self.P_far)
        glTranslated(self.P_Tx, self.P_Ty, self.P_Tz)
        glRotated(self.P_Rx, 1.0, 0.0, 0.0)
        glRotated(self.P_Ry, 0.0, 1.0, 0.0)
        glRotated(self.P_Rz, 0.0, 0.0, 1.0)
        
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()

    

    def OnDraw(self, *args, **kwargs):
        "Draw the window."
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        
        for elem in self.elements:
            if elem['name']=='guidegrid':
                self.DrawGuideGrid()
            elif elem['name']=='triangle':
                r0, r1, r2 = elem['v_arr']
                self.DrawTriangle(r0, r1, r2)
            elif elem['name']=='cube':
                r0, r1 = elem['v_arr']
                self.DrawCube(r0, r1)
            elif elem['name']=='quad':
                r0, r1, r2, r3 = elem['v_arr']
                self.DrawQuad(r0, r1, r2, r3)
                
        #print('self.rectGridOn'+str(self.rectGridOn))
        if self.rectGridOn:
            self.DrawRectGrid()
            if self.RGCellsMarked!=None:
                self.DrawMarkedRectGridCells()
            

    def UpdateDimensions(self, W_3D, dW_3D, V_3D, V_theta, V_phi):
        self.W_3D = W_3D
        self.dW_3D = dW_3D
        self.V_3D = V_3D
        self.V_theta = V_theta
        self.V_phi = V_phi
        

    def DrawGuideGrid(self):
        dx = self.dW_3D[0]
        dy = self.dW_3D[1]
        X = self.W_3D[0]
        Y = self.W_3D[1]
        Z = self.W_3D[2]
        
        glLineWidth(1.5)
        glBegin(GL_LINES)

        if self.axisOn:
            glColor3d(1.0, 0.0, 0.0)
            glVertex3d(0.0, 0.0, 0.0)
            glVertex3d(+X/2, 0.0, 0.0)
            
            glColor3d(0.0, 1.0, 0.0)
            glVertex3d(0.0, 0.0, 0.0)
            glVertex3d(0.0, Y/2, 0.0)

            glColor3d(0.0, 0.0, 1.0)
            glVertex3d(0.0, 0.0, 0.0)
            glVertex3d(0.0, 0.0, Z/2)
        glEnd()

        glLineWidth(1.0)
        glBegin(GL_LINES)
        if self.gridOn:
            glColor3d(0.8, 0.8, 0.8)
            for i in range(-5, 6):
                glVertex3d(i*dx, -Y/2, 0.0)
                glVertex3d(i*dx, +Y/2, 0.0)
            for i in range(-5, 6):
                glVertex3d(-X/2, i*dy, 0.0)
                glVertex3d(+X/2, i*dy, 0.0)
        
        glEnd()
        self.SwapBuffers()


    def DrawTriangle(self, r0, r1, r2, color=None):
        # Drawing an example triangle in the middle of the screen
        glBegin(GL_TRIANGLES)
        if color==None:
            glColor3d(random.random(), random.random(), random.random())
        else:
            glColor3d(color[0], color[1], color[2])
        glVertex3d(r0[0], r0[1], r0[2])
        glVertex3d(r1[0], r1[1], r1[2])
        glVertex3d(r2[0], r2[1], r2[2])
        glEnd()

        self.SwapBuffers()
    
    def DrawCube(self, r0, r1, color=[1.0, 0.0, 0.0]):
        x0, y0, z0 = r0
        x1, y1, z1 = r1

        # position viewer
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        #glRotated(0.0, 1.0, 0.0, 0.0)
        #glTranslated(0.0, 0.0, 0.0)

        # draw six faces of a cube
        glBegin(GL_QUADS)
        glNormal3d(0.0, 0.0, 1.0)
        glColor3d(color[0], color[1], color[2])
        glVertex3d(x0, y0, z1)
        glVertex3d(x1, y0, z1)
        glVertex3d(x1, y1, z1)
        glVertex3d(x0, y1, z1)

        glNormal3d( 0.0, 0.0, -1.0)
        glVertex3d(x0, y0, z0)
        glVertex3d(x0, y1, z0)
        glVertex3d(x1, y1, z0)
        glVertex3d(x1, y0, z0)

        glNormal3d( 0.0, 1.0, 0.0)
        glVertex3d(x0, y1, z0)
        glVertex3d(x0, y1, z1)
        glVertex3d(x1, y1, z1)
        glVertex3d(x1, y1, z0)

        glNormal3d( 0.0, -1.0, 0.0)
        glVertex3d(x0, y0, z0)
        glVertex3d(x0, y0, z1)
        glVertex3d(x1, y0, z1)
        glVertex3d(x1, y0, z0)
        glEnd()
        self.SwapBuffers()

    def DrawQuad(self, r0, r1, r2, r3, normal=None, color=None):
        x0, y0, z0 = r0
        x1, y1, z1 = r1
        x2, y2, z2 = r2
        x3, y3, z3 = r3

        # position viewer
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()

        # draw six faces of a cube
        glBegin(GL_QUADS)
        if color==None:
            glColor3d(random.random(), random.random(), random.random())
        else:
            glColor3d(color[0], color[1], color[2])
        if normal==None:
            glNormal3d(0.0, 0.0, 1.0)
        else:
            glNormal3d(normal[0], normal[1], normal[2])
            
        glVertex3d(x0, y0, z0)
        glVertex3d(x1, y1, z1)
        glVertex3d(x2, y2, z2)
        glVertex3d(x3, y3, z3)
        glEnd()
        self.SwapBuffers()

    def DrawRectGrid(self):
        if self.rectGrid!=None:
            #self.LogText('Drawing rect grid.')
            rg = self.rectGrid
            if rg.N==1:
                nodesPoints = rg.nodesPoints
 
                glMatrixMode(GL_MODELVIEW)
                glLoadIdentity()
                
                glPointSize(3.0)

                glBegin(GL_POINTS)
                glColor3d(1.0, 0.0, 0.0)
                glNormal3d(0.0, 0.0, 1.0)
                for r in nodesPoints:
                    glVertex3d(r[0], 0.0, 0.0)
                glEnd()
            
                
                cellsHier = rg.cellsHier
                n_levels = len(cellsHier)
                
                glBegin(GL_LINES)
                glColor3d(0.0, 0.0, 0.0)
                glNormal3d(0.0, 0.0, 1.0)
                for lev in range(n_levels):
                    CH_lev = cellsHier[lev]
                    for c in CH_lev:
                        points = c[H_CN]
                        for p_ind in points:
                            r = nodesPoints[p_ind]
                            glVertex3d(r[0], 0.0, 0.0)
                glEnd()

                self.SwapBuffers()
            elif rg.N==2:
                nodesPoints = rg.nodesPoints
 
                glMatrixMode(GL_MODELVIEW)
                glLoadIdentity()
                
                glPointSize(3.0)

                glBegin(GL_POINTS)
                glColor3d(1.0, 0.0, 0.0)
                glNormal3d(0.0, 0.0, 1.0)
                for r in nodesPoints:
                    glVertex3d(r[0], r[1], 0.0)
                glEnd()
            
                
                cellsHier = rg.cellsHier
                n_levels = len(cellsHier)
                
                glBegin(GL_LINES)
                glColor3d(0.0, 0.0, 0.0)
                glNormal3d(0.0, 0.0, 1.0)
                for lev in range(n_levels):
                    CH_lev = cellsHier[lev]
                    for c in CH_lev:
                        points = c[H_CN]
                        r = nodesPoints[points[0]]
                        glVertex3d(r[0], r[1], 0.0)
                        r = nodesPoints[points[1]]
                        glVertex3d(r[0], r[1], 0.0)
                        glVertex3d(r[0], r[1], 0.0)
                        r = nodesPoints[points[3]]
                        glVertex3d(r[0], r[1], 0.0)
                        glVertex3d(r[0], r[1], 0.0)
                        r = nodesPoints[points[2]]
                        glVertex3d(r[0], r[1], 0.0)
                        glVertex3d(r[0], r[1], 0.0)
                        r = nodesPoints[points[0]]
                        glVertex3d(r[0], r[1], 0.0)
                glEnd()

                self.SwapBuffers()
            elif rg.N==3:
                nodesPoints = rg.nodesPoints
 
                glMatrixMode(GL_MODELVIEW)
                glLoadIdentity()
                
                glPointSize(3.0)

                glBegin(GL_POINTS)
                glColor3d(1.0, 0.0, 0.0)
                glNormal3d(0.0, 0.0, 1.0)
                for r in nodesPoints:
                    glVertex3d(r[0], r[1], r[2])
                glEnd()
            
                
                cellsHier = rg.cellsHier
                n_levels = len(cellsHier)
                
                glBegin(GL_LINES)
                glColor3d(0.0, 0.0, 0.0)
                glNormal3d(0.0, 0.0, 1.0)
                for lev in range(n_levels):
                    CH_lev = cellsHier[lev]
                    for c in CH_lev:
                        points = c[H_CN]
                        for i in range(len(points)):
                            r_i = nodesPoints[points[i]]
                            for j in range(i+1, len(points)):
                                r_j = nodesPoints[points[j]]
                                if np.sum(np.isclose(r_i, r_j))>=2:
                                    glVertex3d(r_i[0], r_i[1], r_i[2])
                                    glVertex3d(r_j[0], r_j[1], r_j[2])
                glEnd()

                self.SwapBuffers()
        return


    def MarkRectGridCell(self, lev, c_ind, color):
        if self.rectGrid!=None:
            rg = self.rectGrid
            cellsHier = rg.cellsHier
            n_levels = len(cellsHier)
            assert lev<n_levels and c_ind<len(cellsHier[lev])
            if self.RGCellsMarked==None:
                self.RGCellsMarked = []    
            self.RGCellsMarked.append([lev, c_ind, color])
            #print(self.RGCellsMarked)

    def DrawMarkedRectGridCells(self):
        if self.rectGrid!=None:
            if self.RGCellsMarked==None:
                return
            if len(self.RGCellsMarked)==0:
                return
                
            #self.LogText('Drawing rect grid.')
            rg = self.rectGrid
            if rg.N==1:
                nodesPoints = rg.nodesPoints
 
                glMatrixMode(GL_MODELVIEW)
                glLoadIdentity()
                
                glPointSize(3.0)
                glLineWidth(2.0)

                glBegin(GL_POINTS)
                glColor3d(1.0, 0.0, 0.0)
                glNormal3d(0.0, 0.0, 1.0)
                for r in nodesPoints:
                    glVertex3d(r[0], 0.0, 0.0)
                glEnd()
                
                cellsHier = rg.cellsHier
                n_levels = len(cellsHier)
                
                glBegin(GL_LINES)
                glNormal3d(0.0, 0.0, 1.0)
                
                for c_marked in self.RGCellsMarked:
                    lev, c_ind, color = c_marked
                    glColor3d(color[0], color[1], color[2])
                    points = cellsHier[lev][c_ind][H_CN]
                    for p_ind in points:
                        r = nodesPoints[p_ind]
                        glVertex3d(r[0], 0.0, 0.0)
                glEnd()
                    
                glLineWidth(1.0)
                self.SwapBuffers()
            elif rg.N==2:
                nodesPoints = rg.nodesPoints
 
                glMatrixMode(GL_MODELVIEW)
                glLoadIdentity()
                
                glPointSize(3.0)
                glLineWidth(2.0)

                glBegin(GL_POINTS)
                glColor3d(1.0, 0.0, 0.0)
                glNormal3d(0.0, 0.0, 1.0)
                for r in nodesPoints:
                    glVertex3d(r[0], r[1], 0.0)
                glEnd()
            
                
                cellsHier = rg.cellsHier
                n_levels = len(cellsHier)
                
                glBegin(GL_LINES)
                glColor3d(0.0, 0.0, 1.0)
                glNormal3d(0.0, 0.0, 1.0)

                for c_marked in self.RGCellsMarked:
                    lev, c_ind, color = c_marked
                    glColor3d(color[0], color[1], color[2])
                    points = cellsHier[lev][c_ind][H_CN]
                    r = nodesPoints[points[0]]
                    glVertex3d(r[0], r[1], 0.0)
                    r = nodesPoints[points[1]]
                    glVertex3d(r[0], r[1], 0.0)
                    glVertex3d(r[0], r[1], 0.0)
                    r = nodesPoints[points[3]]
                    glVertex3d(r[0], r[1], 0.0)
                    glVertex3d(r[0], r[1], 0.0)
                    r = nodesPoints[points[2]]
                    glVertex3d(r[0], r[1], 0.0)
                    glVertex3d(r[0], r[1], 0.0)
                    r = nodesPoints[points[0]]
                    glVertex3d(r[0], r[1], 0.0)
                glEnd()

                glLineWidth(1.0)
                self.SwapBuffers()
            elif rg.N==3:
                nodesPoints = rg.nodesPoints
 
                glMatrixMode(GL_MODELVIEW)
                glLoadIdentity()
                
                glPointSize(3.0)
                glLineWidth(2.0)

                glBegin(GL_POINTS)
                glColor3d(1.0, 0.0, 0.0)
                glNormal3d(0.0, 0.0, 1.0)
                for r in nodesPoints:
                    glVertex3d(r[0], r[1], r[2])
                glEnd()
            
                
                cellsHier = rg.cellsHier
                n_levels = len(cellsHier)
                
                glBegin(GL_LINES)
                glNormal3d(0.0, 0.0, 1.0)

                for c_marked in self.RGCellsMarked:
                    lev, c_ind, color = c_marked
                    glColor3d(color[0], color[1], color[2])
                    points = cellsHier[lev][c_ind][H_CN]
                    for i in range(len(points)):
                        r_i = nodesPoints[points[i]]
                        for j in range(i+1, len(points)):
                            r_j = nodesPoints[points[j]]
                            if np.sum(np.isclose(r_i, r_j))>=2:
                                glVertex3d(r_i[0], r_i[1], r_i[2])
                                glVertex3d(r_j[0], r_j[1], r_j[2])
                glEnd()

                glLineWidth(1.0)
                self.SwapBuffers()
        return




