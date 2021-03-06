## RecGrid2D.py 
## 2D rectangular grid


__all__ = ["RG2D", "GUIMaker2D"]


from Electromagnetics.SymExprTree import *
#from Electromagnetics import Misc
from scipy.sparse import coo_matrix
import numpy as np
import scipy as sp
from sympy import Symbol

from Electromagnetics.Misc import *
from Electromagnetics.Misc import Point2D



SIDE_INSIDE_CELL = -2
SIDE_OUTSIDE_DOMAIN = -1
NUM_NEG_UNUSED = -3     ## negative number not used for any purposes

class RG2D:
    """ 2D structured rectangular grid - each cell can be refined on demand
    - the adjacent cells either are the same size, or one is half the other, if
    this condition does not hold cells are reifined automatically until this 
    condition is satisfied
    """
    def __init__(self, x0, x1, y0, y1, dx, dy, x_sym=Symbol('x'), y_sym=Symbol('y')):
        nx = int((x1 - x0)/dx) + 1
        ny = int((y1 - y0)/dy) + 1
        dx = (x1 - x0)/nx
        dy = (y1 - y0)/ny
        
        self.width = x1 - x0
        self.height = y1 - y0
        
        self.cells = [[]]*(nx*ny)
        self.sides = [[]]*(nx*(ny+1) + (nx+1)*ny)
        self.nodes = [[]]*((nx+1)*(ny+1))
        self.sidesDiv = [[]]        # index zero will remain unused
        self.nodesPoints = [()]*((nx+1)*(ny+1))
                
        for j in range(ny):
            for i in range(nx):
                ## cell :[left side, right side, lower side, upper side]
                self.cells[self._uniformGridCellInd(i, j, nx, ny)] = \
                    [self._uniformGridSideIndY(i, j, nx, ny),\
                     self._uniformGridSideIndY(i+1, j, nx, ny),\
                     self._uniformGridSideIndX(i, j, nx, ny),\
                     self._uniformGridSideIndX(i, j+1, nx, ny),\
                     0]
                
        for j in range(ny+1):
            y = y0 + j*dy
            for i in range(nx+1):
                x = x0 + i*dx
                ## node :[left side, right side, lower side, upper side]
                self.nodes[self._uniformGridNodeInd(i, j, nx, ny)] = \
                [self._uniformGridSideIndX(i-1, j, nx, ny),\
                 self._uniformGridSideIndX(i, j, nx, ny),\
                 self._uniformGridSideIndY(i, j-1, nx, ny),\
                 self._uniformGridSideIndY(i, j, nx, ny)]
                 
                self.nodesPoints[self._uniformGridNodeInd(i, j, nx, ny)] = \
                    Point2D(x,y)
                
        
        for j in range(ny+1):
            for i in range(nx):
                ## side: [left/lower node, right/upper node, 
                ## left/lower cell, right/upper cell, level]
                self.sides[self._uniformGridSideIndX(i, j, nx, ny)] = \
                [self._uniformGridNodeInd(i, j, nx, ny),\
                 self._uniformGridNodeInd(i+1, j, nx, ny),\
                 self._uniformGridCellInd(i, j-1, nx, ny),\
                 self._uniformGridCellInd(i, j, nx, ny), 'X']
                 
        for j in range(ny):
            for i in range(nx+1):
                ## side: [left/lower node, right/upper node, 
                ## left/lower cell, right/upper cell, level]
                self.sides[self._uniformGridSideIndY(i, j, nx, ny)] = \
                [self._uniformGridNodeInd(i, j, nx, ny),\
                 self._uniformGridNodeInd(i, j+1, nx, ny),\
                 self._uniformGridCellInd(i-1, j, nx, ny),\
                 self._uniformGridCellInd(i, j, nx, ny), 'Y']
             
        self.x = x_sym
        self.y = y_sym

        self.dX = dx
        self.dY = dy
        self.SetCellSizesForDifferentLevels()
        self.Vars = {}  ## symbolic variables assigned to cells, nodes, sides...
        self.Regions = {}  ## equations for different regions 
                        ## 'region_name': [[eq0, eq1...eqn], [x0, x1,...xn]]
                        ## each equation is associated with a given variable 
                        ## (the derivatives will be taken with respect to the
                        ## location of that variable)
        self.RegionBoundary = {} ## equations to be used at the boundary of each region

        self.initialMatrix_coo = None
        self.initialMatrix_shared_coo = None
        self.finalMatrix_coo = None
        
        self.ElemFinalToElemReduced = None
        return
                
    def _uniformGridCellInd(self, i, j, nx, ny):
        if 0<=i<nx and 0<=j<ny:
            return j*nx + i
        else:
            return -1   #used to recognize the indices falling outside the domain
    
    def _uniformGridNodeInd(self, i, j, nx, ny):
        if 0<=i<=nx and 0<=j<=ny:
            return j*(nx+1) + i
        else:
            return -1

    def _uniformGridSideIndX(self, i, j, nx, ny):
        if 0<=i<nx and 0<=j<=ny:
            return j*nx + i
        else:
            return -1

    def _uniformGridSideIndY(self, i, j, nx, ny):
        if 0<=i<=nx and 0<=j<ny:
            return nx*(ny+1) + j*(nx+1) + i
        else:
            return -1


    def testGrid(self):
        for i in range(len(self.cells)):
            ## l:left  r:right  d:down  u:up
            sl = self.cells[i][0]   ## side_left
            sr = self.cells[i][1]
            sd = self.cells[i][2]
            su = self.cells[i][3]
            ndl, ndl_, ndr, ndr_, nul, nul_, nur, nur_ = [-1]*8
            if sl>=0:
                ndl = self.sides[sl][0]     ## node_down_left
            else:
                ndl = self.sides[self.sidesDiv[abs(sl)][0]][0]
            if sd>=0:
                ndl_ = self.sides[sd][0]
            else:
                ndl_ = self.sides[self.sidesDiv[abs(sd)][0]][0]
            if sr>=0:
                ndr = self.sides[sr][0]
            else:
                ndr = self.sides[self.sidesDiv[abs(sr)][0]][0]
            if sd>=0:
                ndr_ = self.sides[sd][1]
            else:
                ndr_ = self.sides[self.sidesDiv[abs(sd)][1]][1]
            if sl>=0:
                nul = self.sides[sl][1]
            else:
                nul = self.sides[self.sidesDiv[abs(sl)][1]][1]
            if su>=0:
                nul_ = self.sides[su][0]
            else:
                nul_ = self.sides[self.sidesDiv[abs(su)][0]][0]
            if sr>=0:
                nur = self.sides[sr][1]
            else:
                nur = self.sides[self.sidesDiv[abs(sr)][1]][1]
            if su>=0:
                nur_ = self.sides[su][1]
            else:
                nur_ = self.sides[self.sidesDiv[abs(su)][1]][1]
            
            if ndl!=ndl_ or ndr!=ndr_ or nul!=nul_ or nur!=nur_:
                print('cell: ', i)
                print('ndl:', ndl, ' ndl_:', ndl_)
                print('ndr:', ndr, ' ndr_:', ndr_)
                print('nul:', nul, ' nul_:', nul_)
                print('nur:', nur, ' nur_:', nur_)
                return [False, 'ndl!=ndl_ or ndr!=ndr_ or nul!=nul_ or nur!=nur_']
                
            if sl>=0:
                sl_ = self.nodes[ndl][3]
                sl__ = self.nodes[nul][2]
                if not(sl==sl_==sl__):
                    print('cell= ', i)
                    print('sl=', sl, ' sl_=', sl_, ' sl__=', sl__)
                    return [False, 'not(sl==sl_==sl__)']
            if sr>=0:
                sr_ = self.nodes[ndr][3]
                sr__ = self.nodes[nur][2]
                if not(sr==sr_==sr__):
                    print('cell= ', i)
                    print('sr=', sr, ' sr_=', sr_, ' sr__=', sr__)
                    return [False, 'not(sr==sr_==sr__)']
            if sd>=0:
                sd_ = self.nodes[ndl][1]
                sd__ = self.nodes[ndr][0]
                if not(sd==sd_==sd__):
                    print('cell= ', i)
                    print('sd=', sd, ' sd_=', sd_, ' sd__=', sd__)
                    return [False, 'not(sd==sd_==sd__)']
            if su>=0:
                su_ = self.nodes[nul][1]
                su__ = self.nodes[nur][0]
                if not(su==su_==su__):
                    print('cell= ', i)
                    print('su=', su, ' su_=', su_, ' su__=', su__)
                    return [False, 'not(su==su_==su__)']
            c, c_, c__, c___ = [-1]*4
            if sl>=0:
                c = self.sides[sl][3]
            else:
                c = self.sides[self.sidesDiv[abs(sl)][0]][3]
            if sr>=0:
                c_ = self.sides[sr][2]
            else:
                c_ = self.sides[self.sidesDiv[abs(sr)][0]][2]
            if sd>=0:
                c__ = self.sides[sd][3]
            else:
                c__ = self.sides[self.sidesDiv[abs(sd)][0]][3]
            if su>=0:
                c___ = self.sides[su][2]
            else:
                c___ = self.sides[self.sidesDiv[abs(su)][0]][2]
            
            if not(c==c_==c__==c___==i):
                print('cell= ', i)
                print('sl=', sl, ' sr=', sr, ' sd=', sd, ' su=', su)
                print('c=', c, ' c_=', c_, ' c__=', c__, ' c___=', c___)
                return [False, 'not(c==c_==c__==c___==i)']

        for i in range(len(self.nodes)):
            sl = self.nodes[i][0]
            sr = self.nodes[i][1]
            sd = self.nodes[i][2]
            su = self.nodes[i][3]
            if sl>=0:
                if self.sides[sl][1]!=i:
                    return [False, 'self.sides[sl][1]!=i']
            if sr>=0:
                if self.sides[sr][0]!=i:
                    return [False, 'self.sides[sr][0]!=i']
            if sd>=0:
                if self.sides[sd][1]!=i:
                    return [False, 'self.sides[sd][1]!=i']
            if su>=0:
                if self.sides[su][0]!=i:
                    return [False, 'self.sides[su][0]!=i']
        
        
        for i in range(len(self.sides)):
            nld = self.sides[i][0]
            nru = self.sides[i][1]
            cld = self.sides[i][2]
            cru = self.sides[i][3]
            if self.sides[i][4]=='X':
                if self.nodes[nld][1]!=i:
                    return [False, 'self.nodes[nld][1]!=i']
                if self.nodes[nru][0]!=i:
                    return [False, 'self.nodes[nru][0]!=i']
                if cld>=0:
                    if self.cells[cld][3]>=0:
                        if self.cells[cld][3]!=i:
                            return [False, 'self.cells[cld][3]!=i']
                if cru>=0:
                    if self.cells[cru][2]>=0:
                        if self.cells[cru][2]!=i:
                            return [False, 'self.cells[cru][2]!=i']
            else:
                if self.nodes[nld][3]!=i:
                    print('side = ', i, '   ', self.sides[i])
                    print('nld=', nld, ' nru=', nru, ' cld=', cld, ' cru=', cru)
                    return [False, 'self.nodes[nld][3]!=i']
                if self.nodes[nru][2]!=i:
                    return [False, 'self.nodes[nru][2]!=i']
                if cld>=0:
                    if self.cells[cld][1]>=0:
                        if self.cells[cld][1]!=i:
                            return [False, 'self.cells[cld][1]!=i']
                if cru>=0:
                    if self.cells[cru][0]>=0:
                        if self.cells[cru][0]!=i:
                            return [False, 'self.cells[cru][0]!=i']
        
        return [True]



    def RefineSingleCell(self, cell_ind):
        ## l:left r:right d:down u:up
        n_cells = len(self.cells)
        n_sides = len(self.sides)
        n_nodes = len(self.nodes)
        node_ind_new = n_nodes
        side_ind_new = n_sides
        cells_inside = [cell_ind, n_cells, n_cells+1, n_cells+2] #[ld, lu, rd, ru] 
        cells_outside = [-1]*8  #[ld, lu, rd, ru, dl, dr, ul, ur]
        nodes_corner = [-1]*4   #[ld, lu, rd, ru]
        nodes_sides = [-1]*4    #[l, r, d, u]
        node_center = -1
        sides_sides = [-1]*8    #[ld, lu, rd, ru, dl, dr, ul, ur]
        sides_inside = [-1]*4   #[l, r, d, u]
        sides_corner_outside = [-1]*8   #[ldl, ldd, rdr, rdd, lul, luu, rur, ruu]
        sides_sides_outside = [-1]*4    #[l, r, d, u]
        cell_arr = self.cells[cell_ind]
        cell_level = cell_arr[4]
        nodespoints = []
        #-- nodes_corner, nodes_sides, sides_sides, sides_sides_outside
        ld, lu, rd, ru, dl, dr, ul, ur = range(8) #sides_sides, cells_outside
        l, r, d, u = range(4)   #nodes_sides, sides_sides_outside and self.nodes indices
        ldl, ldd, rdr, rdd, lul, luu, rur, ruu = range(8)
        s_nld, s_nru, s_cld, s_cru = range(4)   #sides::  s:side n:node c:cell
        sl = cell_arr[l]
        if sl>=0:  #left side
            nodes_corner[ld] = self.sides[sl][s_nld]
            nodes_corner[lu] = self.sides[sl][s_nru]
            nodes_sides[l] = node_ind_new
            node_ind_new += 1
            nodespoints.append((self.nodesPoints[nodes_corner[ld]] + \
                self.nodesPoints[nodes_corner[lu]])/2.0)
            sides_sides[ld] = sl
            sides_sides[lu] = side_ind_new
            side_ind_new += 1
            if self.sides[sl][s_cld] >=0:
                sides_sides_outside[l] = SIDE_INSIDE_CELL
            else:
                sides_sides_outside[l] = -1
            cells_outside[ld] = self.sides[sl][s_cld]
            cells_outside[lu] = self.sides[sl][s_cld]
        else:
            nodes_corner[ld] = self.sides[self.sidesDiv[abs(sl)][0]][s_nld]
            nodes_corner[lu] = self.sides[self.sidesDiv[abs(sl)][1]][s_nru]
            nodes_sides[l] = self.sides[self.sidesDiv[abs(sl)][0]][s_nru]
            sides_sides[ld] = self.sidesDiv[abs(sl)][0]
            sides_sides[lu] = self.sidesDiv[abs(sl)][1]
            sides_sides_outside[l] = self.nodes[nodes_sides[l]][l]
            cells_outside[ld] = self.sides[self.sidesDiv[abs(sl)][0]][s_cld]
            cells_outside[lu] = self.sides[self.sidesDiv[abs(sl)][1]][s_cld]
        sr = cell_arr[r]
        if sr>=0:  #right side
            nodes_corner[rd] = self.sides[sr][s_nld]
            nodes_corner[ru] = self.sides[sr][s_nru]
            nodes_sides[r] = node_ind_new
            node_ind_new += 1
            nodespoints.append((self.nodesPoints[nodes_corner[rd]] + \
                self.nodesPoints[nodes_corner[ru]])/2.0)
            sides_sides[rd] = sr
            sides_sides[ru] = side_ind_new
            side_ind_new += 1
            if self.sides[sr][s_cru] >=0:
                sides_sides_outside[r] = SIDE_INSIDE_CELL
            else:
                sides_sides_outside[r] = -1
            cells_outside[rd] = self.sides[sr][s_cru]
            cells_outside[ru] = self.sides[sr][s_cru]
        else:
            nodes_corner[rd] = self.sides[self.sidesDiv[abs(sr)][0]][s_nld]
            nodes_corner[ru] = self.sides[self.sidesDiv[abs(sr)][1]][s_nru]
            nodes_sides[r] = self.sides[self.sidesDiv[abs(sr)][0]][s_nru]
            sides_sides[rd] = self.sidesDiv[abs(sr)][0]
            sides_sides[ru] = self.sidesDiv[abs(sr)][1]
            sides_sides_outside[r] = self.nodes[nodes_sides[r]][r]
            cells_outside[rd] = self.sides[self.sidesDiv[abs(sr)][0]][s_cru]
            cells_outside[ru] = self.sides[self.sidesDiv[abs(sr)][1]][s_cru]
        sd = cell_arr[d]
        if sd>=0:  #down side
            nodes_sides[d] = node_ind_new
            node_ind_new += 1
            nodespoints.append((self.nodesPoints[nodes_corner[ld]] + \
                self.nodesPoints[nodes_corner[rd]])/2.0)
            sides_sides[dl] = sd
            sides_sides[dr] = side_ind_new
            side_ind_new += 1
            if self.sides[sd][s_cld] >=0:
                sides_sides_outside[d] = SIDE_INSIDE_CELL
            else:
                sides_sides_outside[d] = -1
            cells_outside[dl] = self.sides[sd][s_cld]
            cells_outside[dr] = self.sides[sd][s_cld]
        else:
            nodes_sides[d] = self.sides[self.sidesDiv[abs(sd)][0]][s_nru]
            sides_sides[dl] = self.sidesDiv[abs(sd)][0]
            sides_sides[dr] = self.sidesDiv[abs(sd)][1]
            sides_sides_outside[d] = self.nodes[nodes_sides[d]][d]
            cells_outside[dl] = self.sides[self.sidesDiv[abs(sd)][0]][s_cld]
            cells_outside[dr] = self.sides[self.sidesDiv[abs(sd)][1]][s_cld]
        su = cell_arr[u]
        if su>=0:  #up side
            nodes_sides[u] = node_ind_new
            node_ind_new += 1
            nodespoints.append((self.nodesPoints[nodes_corner[lu]] + \
                self.nodesPoints[nodes_corner[ru]])/2.0)
            sides_sides[ul] = su
            sides_sides[ur] = side_ind_new
            side_ind_new += 1
            if self.sides[su][s_cru] >=0:
                sides_sides_outside[u] = SIDE_INSIDE_CELL
            else:
                sides_sides_outside[u] = -1
            cells_outside[ul] = self.sides[su][s_cru]
            cells_outside[ur] = self.sides[su][s_cru]
        else:
            nodes_sides[u] = self.sides[self.sidesDiv[abs(su)][0]][s_nru]
            sides_sides[ul] = self.sidesDiv[abs(su)][0]
            sides_sides[ur] = self.sidesDiv[abs(su)][1]
            sides_sides_outside[u] = self.nodes[nodes_sides[u]][u]
            cells_outside[ul] = self.sides[self.sidesDiv[abs(su)][0]][s_cru]
            cells_outside[ur] = self.sides[self.sidesDiv[abs(su)][1]][s_cru]
        node_center = node_ind_new
        node_ind_new += 1
        nodespoints.append((
    self.nodesPoints[nodes_corner[lu]] + self.nodesPoints[nodes_corner[ld]] +\
    self.nodesPoints[nodes_corner[ru]] + self.nodesPoints[nodes_corner[rd]] )/4.0)
        #-- sides_inside
        l, r, d, u = range(4)
        sides_inside[l] = side_ind_new
        side_ind_new += 1
        sides_inside[r] = side_ind_new
        side_ind_new += 1
        sides_inside[d] = side_ind_new
        side_ind_new += 1
        sides_inside[u] = side_ind_new
        side_ind_new += 1
        #-- sides_corner_outside
        sides_corner_outside[ldl] = self.nodes[nodes_corner[ld]][l]   # ldl
        sides_corner_outside[ldd] = self.nodes[nodes_corner[ld]][d]   # ldd
        sides_corner_outside[rdr] = self.nodes[nodes_corner[rd]][r]   # rdr
        sides_corner_outside[rdd] = self.nodes[nodes_corner[rd]][d]   # rdd
        sides_corner_outside[lul] = self.nodes[nodes_corner[lu]][l]   # lul
        sides_corner_outside[luu] = self.nodes[nodes_corner[lu]][u]   # luu
        sides_corner_outside[rur] = self.nodes[nodes_corner[ru]][r]   # rur
        sides_corner_outside[ruu] = self.nodes[nodes_corner[ru]][u]   # ruu
        ##---- set new list components
        #-- cells
        cells_inside_new = [
        [sides_sides[ld], sides_inside[d], sides_sides[dl], sides_inside[l], cell_level+1],  #ld
        [sides_sides[lu], sides_inside[u], sides_inside[l], sides_sides[ul], cell_level+1],  #lu
        [sides_inside[d], sides_sides[rd], sides_sides[dr], sides_inside[r], cell_level+1],  #rd
        [sides_inside[u], sides_sides[ru], sides_inside[r], sides_sides[ur], cell_level+1]]  #ru
        #-- sides
        sides_sides_new = [
        [nodes_corner[ld], nodes_sides[l], cells_outside[ld], cells_inside[ld], 'Y'],    #ld
        [nodes_sides[l], nodes_corner[lu], cells_outside[lu], cells_inside[lu], 'Y'],    #lu
        [nodes_corner[rd], nodes_sides[r], cells_inside[rd], cells_outside[rd], 'Y'],    #rd
        [nodes_sides[r], nodes_corner[ru], cells_inside[ru], cells_outside[ru], 'Y'],    #ru
        [nodes_corner[ld], nodes_sides[d], cells_outside[dl], cells_inside[ld], 'X'],    #dl
        [nodes_sides[d], nodes_corner[rd], cells_outside[dr], cells_inside[rd], 'X'],    #dr
        [nodes_corner[lu], nodes_sides[u], cells_inside[lu], cells_outside[ul], 'X'],    #ul
        [nodes_sides[u], nodes_corner[ru], cells_inside[ru], cells_outside[ur], 'X']]    #ur
        sides_inside_new = [
        [nodes_sides[l], node_center, cells_inside[ld], cells_inside[lu], 'X'], #l
        [node_center, nodes_sides[r], cells_inside[rd], cells_inside[ru], 'X'], #r
        [nodes_sides[d], node_center, cells_inside[ld], cells_inside[rd], 'Y'], #d
        [node_center, nodes_sides[u], cells_inside[lu], cells_inside[ru], 'Y']] #u
        #-- nodes
        nodes_corner_new = [
        [sides_corner_outside[ldl], sides_sides[dl], sides_corner_outside[ldd], sides_sides[ld]], #ld
        [sides_corner_outside[lul], sides_sides[ul], sides_sides[lu], sides_corner_outside[luu]], #lu
        [sides_sides[dr], sides_corner_outside[rdr], sides_corner_outside[rdd], sides_sides[rd]], #rd
        [sides_sides[ur], sides_corner_outside[rur], sides_sides[ru], sides_corner_outside[ruu]]] #ru
        nodes_sides_new = [
        [sides_sides_outside[l], sides_inside[l], sides_sides[ld], sides_sides[lu]], #l
        [sides_inside[r], sides_sides_outside[r], sides_sides[rd], sides_sides[ru]], #r
        [sides_sides[dl], sides_sides[dr], sides_sides_outside[d], sides_inside[d]], #d
        [sides_sides[ul], sides_sides[ur], sides_inside[u], sides_sides_outside[u]]] #u
        node_center_new = [sides_inside[l], sides_inside[r], sides_inside[d], sides_inside[u]]
        ##---- extend lists
        n_cell_to_add = 3
        n_side_to_add = 4
        for s in sides_sides:
            if s >= n_sides:
                n_side_to_add += 1
        n_node_to_add = 1
        for n in nodes_sides:
            if n >= n_nodes:
                n_node_to_add += 1
        self.cells.extend([[]]*n_cell_to_add)
        self.sides.extend([[]]*n_side_to_add)
        self.nodes.extend([[]]*n_node_to_add)
        ##--- adjusting outside cells (devided sides.. if any)
        sl = cell_arr[l]
        if sl>=0:
            cl = self.sides[sl][s_cld]
            if cl>=0:
                self.cells[cl][r] = -len(self.sidesDiv)
                self.sidesDiv.append([sides_sides[ld], sides_sides[lu]])
        sr = cell_arr[r]
        if sr>=0:
            cr = self.sides[sr][s_cru]
            if cr>=0:
                self.cells[cr][l] = -len(self.sidesDiv)
                self.sidesDiv.append([sides_sides[rd], sides_sides[ru]])
        sd = cell_arr[d]
        if sd>=0:
            cd = self.sides[sd][s_cld]
            if cd>=0:
                self.cells[cd][u] = -len(self.sidesDiv)
                self.sidesDiv.append([sides_sides[dl], sides_sides[dr]])
        su = cell_arr[u]
        if su>=0:
            cu = self.sides[su][s_cru]
            if cu>=0:
                self.cells[cu][d] = -len(self.sidesDiv)
                self.sidesDiv.append([sides_sides[ul], sides_sides[ur]])
        ##--- copying final lists
        for i, c in enumerate(cells_inside):
            self.cells[c] = cells_inside_new[i]
        for i, n in enumerate(nodes_corner):
            self.nodes[n] = nodes_corner_new[i]
        for i, n in enumerate(nodes_sides):
            self.nodes[n] = nodes_sides_new[i]
        self.nodes[node_center] = node_center_new
        for i, s in enumerate(sides_sides):
            self.sides[s] = sides_sides_new[i]
        for i, s in enumerate(sides_inside):
            self.sides[s] = sides_inside_new[i]
        self.nodesPoints.extend(nodespoints)
        return
        
    def RefineCells(self, cells_in):
        cells_in_extended = self.MarkCellsToRefine(cells_in)
        print('Number of cells to refine: ', len(cells_in_extended))
        for c in cells_in_extended:
            self.RefineSingleCell(c)
        return

    def MarkCellsToRefine(self, cells_in):
        """It takes the cells_in input cells and marks cells in a way that after 
        refinement for each cell, the difference in level between the neighboring cells
        is at most +1 or -1
        """
        marked_cells = self.getCellsMarked(cells_in)
        cells_iter = cells_in
        mark_next = 2
        cells_to_add = []
        while(True):
            n_marked = 0
            for c in cells_iter:
                level = self.cells[c][4]
                cells_nb = self.Cell_GetNeighborCells(c)
                for c_nb in cells_nb:
                    if c_nb>=0:
                        if self.cells[c_nb][4]<level:
                            if marked_cells[c_nb] == 0:
                                marked_cells[c_nb] = mark_next
                                n_marked += 1
            if n_marked==0:
                break
            cells_marked_next = self.PackMarkedIndices(marked_cells, mark_next)
            mark_next += 1
            cells_iter = cells_marked_next
            cells_to_add.extend(cells_marked_next)
        return self.SortCellsWRTLevels(cells_in + cells_to_add)
        
    
    def Cell_GetNeighborCells(self, cell):
        ld, lu, rd, ru, dl, dr, ul, ur = range(8) #sides_sides, cells_outside
        s_nld, s_nru, s_cld, s_cru = range(4)   #sides::  s:side n:node c:cell
        cell_nb = [-1]*8
        sl, sr, sd, su, _ = self.cells[cell]
        if sl>=0:
            cell_nb[ld] = self.sides[sl][s_cld]
        else:
            cell_nb[ld] = self.sides[self.sidesDiv[abs(sl)][0]][s_cld]
            cell_nb[lu] = self.sides[self.sidesDiv[abs(sl)][1]][s_cld]
        if sr>=0:
            cell_nb[rd] = self.sides[sr][s_cru]
        else:
            cell_nb[rd] = self.sides[self.sidesDiv[abs(sr)][0]][s_cru]
            cell_nb[ru] = self.sides[self.sidesDiv[abs(sr)][1]][s_cru]
        if sd>=0:
            cell_nb[dl] = self.sides[sd][s_cld]
        else:
            cell_nb[dl] = self.sides[self.sidesDiv[abs(sd)][0]][s_cld]
            cell_nb[dr] = self.sides[self.sidesDiv[abs(sd)][1]][s_cld]
        if su>=0:
            cell_nb[ul] = self.sides[su][s_cru]
        else:
            cell_nb[ul] = self.sides[self.sidesDiv[abs(su)][0]][s_cru]
            cell_nb[ur] = self.sides[self.sidesDiv[abs(su)][1]][s_cru]
        return cell_nb
        
    def Node_GetNeighboringCells(self, node):
        s_nld, s_nru, s_cld, s_cru = range(4)   #sides::  s:side n:node c:cell
        ld, lu, rd, ru = range(4)
        cell_nb = [NUM_NEG_UNUSED]*4
        sl, sr, sd, su = self.nodes[node]
        if sl>=0:
            cell_nb[ld] = self.sides[sl][s_cld]
            cell_nb[lu] = self.sides[sl][s_cru]
        if sr>=0:
            cell_nb[rd] = self.sides[sr][s_cld]
            cell_nb[ru] = self.sides[sr][s_cru]
        if sd>=0:
            cell_nb[ld] = self.sides[sd][s_cld]
            cell_nb[rd] = self.sides[sd][s_cru]
        if su>=0:
            cell_nb[lu] = self.sides[su][s_cld]
            cell_nb[ru] = self.sides[su][s_cru]
        return cell_nb
        
        
    def getCellsMarked(self, cell_arr):
        ## cells_marked[cell_arr[i]] = 1
        cells_marked = [0]*len(self.cells)
        for c in cell_arr:
            cells_marked[c] = 1
        return cells_marked
        

    def PackMarkedIndices(self, marked_indices, mark):
        n_mark = 0
        for c in marked_indices:
            if c==mark:
                n_mark += 1
        packed_arr = [-1]*n_mark
        n_mark = 0
        for i, c in enumerate(marked_indices):
            if c==mark:
                packed_arr[n_mark] = i
                n_mark += 1
        return packed_arr
        
    def SortCellsWRTLevels(self, cells_in):
        levels = [-1]*len(cells_in)
        for i, c in enumerate(cells_in):
            levels[i] = self.cells[c][4]
        z = dict(zip(cells_in, levels))
        zs = sorted(z, key=z.get)
        return zs
        
    def GetBoundaryNodesAndSides(self):
        ## The algorithm assumes a rectangular domain
        l, r, d, u = range(4)
        s_nld, s_nru, s_cld, s_cru = range(4)   #sides::  s:side n:node c:cell
        n_ll = self.GetDomainNodeLeftDown()
        n_b = [n_ll]    #nodes boundary
        s_b = []    #sides boundary
        n_last = n_ll
        while True:
            sr = self.nodes[n_last][r]
            if sr>=0:
                s_b.append(sr)
                n_last = self.sides[sr][s_nru]
                n_b.append(n_last)
            else:
                break
        while True:
            su = self.nodes[n_last][u]
            if su>=0:
                s_b.append(su)
                n_last = self.sides[su][s_nru]
                n_b.append(n_last)
            else:
                break
        while True:
            sl = self.nodes[n_last][l]
            if sl>=0:
                s_b.append(sl)
                n_last = self.sides[sl][s_nld]
                n_b.append(n_last)
            else:
                break
        while True:
            sd = self.nodes[n_last][d]
            if sd>=0:
                s_b.append(sd)
                n_last = self.sides[sd][s_nld]
                n_b.append(n_last)
            else:
                break
        assert n_b[-1]==n_ll
        n_b.pop()
        return [n_b, s_b]
        
    def GetDomainNodeLeftDown(self):
        ## Find the lower left node of the rectangular domain
        n_ll = 0
        nll_x = self.nodesPoints[0].x
        nll_y = self.nodesPoints[0].y
        for i, n in enumerate(self.nodesPoints):
            if n.x<nll_x or n.y<nll_y:
                n_ll = i
                nll_x = n.x
                nll_y = n.y
        return n_ll
        
    def SetCellSizesForDifferentLevels(self, n_levels=50):
        self.sizeX = [0.0]*n_levels
        self.sizeY = [0.0]*n_levels
        self.sizeX[0] = self.dX
        self.sizeY[0] = self.dY
        for i in range(1, n_levels):
            self.sizeX[i] = self.sizeX[i-1]/2.0
            self.sizeY[i] = self.sizeY[i-1]/2.0
        return
        
    def GetCellSizeX(self,level):
        return self.sizeX[level]
    def GetCellSizeY(self, level):
        return self.sizeY[level]
        
    def GetSideLevel(self, side):
        s_nld, s_nru, s_cld, s_cru = range(4)   #sides::  s:side n:node c:cell
        lev = 4
        cld = self.sides[side][s_cld]
        cru = self.sides[side][s_cru]
        cld_level = -1
        cru_level = -1
        if cld>=0:
            cld_level = self.cells[cld][lev]
        if cru>=0:
            cru_level = self.cells[cru][lev]
        side_level = max(cld_level, cru_level)
        return side_level
        
    def GetSideSize(self, side):
        XYind = 4
        side_XY = self.sides[side][XYind]
        side_level = self.GetSideLevel(side)
        if side_XY=='X':
            return self.GetCellSizeX(side_level)
        else:
            assert side_XY=='Y'
            return self.GetCellSizeY(side_level)
            
    def GetCellSize(self, cell):
        lev = 4
        cell_arr = self.cells[cell]
        level = cell_arr[lev]
        return [self.sizeX[level], self.sizeY[level]]
        
    def Cell_dxSY(self, cell):
        ## d(Side_Y)/dx at the center of the cell
        sl, sr, sd, su, level = self.cells[cell]
        dxSY_ind = []
        dxSY = []
        dx = self.GetCellSizeX(level)
        if sl>=0:
            dxSY_ind.append(sl)
            dxSY.append(-1.0/dx)
        else:
            dxSY_ind.append(self.sidesDiv[abs(sl)][0])
            dxSY_ind.append(self.sidesDiv[abs(sl)][1])
            dxSY.append(-0.5/dx)
            dxSY.append(-0.5/dx)
        if sr>=0:
            dxSY_ind.append(sr)
            dxSY.append(1.0/dx)
        else:
            dxSY_ind.append(self.sidesDiv[abs(sr)][0])
            dxSY_ind.append(self.sidesDiv[abs(sr)][1])
            dxSY.append(0.5/dx)
            dxSY.append(0.5/dx)
        return [dxSY_ind, dxSY]
                
    def Cell_dySX(self, cell):
        ## d(Side_X)/dy at the center of the cell
        sl, sr, sd, su, level = self.cells[cell]
        dySX_ind = []
        dySX = []
        dy = self.GetCellSizeY(level)
        if sd>=0:
            dySX_ind.append(sd)
            dySX.append(-1.0/dy)
        else:
            dySX_ind.append(self.sidesDiv[abs(sd)][0])
            dySX_ind.append(self.sidesDiv[abs(sd)][1])
            dySX.append(-0.5/dy)
            dySX.append(-0.5/dy)
        if su>=0:
            dySX_ind.append(su)
            dySX.append(1.0/dy)
        else:
            dySX_ind.append(self.sidesDiv[abs(su)][0])
            dySX_ind.append(self.sidesDiv[abs(su)][1])
            dySX.append(0.5/dy)
            dySX.append(0.5/dy)
        return [dySX_ind, dySX]
        
    def Node_dxSX(self, node, onesided=False):
        # SX: side along X
        # onesided = 'ld' returns the portion of the FD related to the left/down section
        # onesided = 'ru' returns the portion of the FD related to the right/up section
        l, r, d, u = range(4)
        s_nld, s_nru, s_cld, s_cru = range(4)   #sides::  s:side n:node c:cell
        sl, sr, sd, su = self.nodes[node]
        sl_level, sr_level = -1, -1
        dxSX_ind = []
        dxSX = []
        if sl>=0:
            sl_level = self.GetSideLevel(sl)
        if sr>=0:
            sr_level = self.GetSideLevel(sr)
        if sl_level==sr_level:
            dx = self.GetCellSizeX(sl_level)
            if not onesided:
                dxSX_ind.append(sl)
                dxSX.append(-1.0/dx)
                dxSX_ind.append(sr)
                dxSX.append(1.0/dx)
            elif onesided=='ld':
                dxSX_ind.append(sl)
                dxSX.append(-1.0/dx)
            else:
                assert onesided=='ru'
                dxSX_ind.append(sr)
                dxSX.append(1.0/dx)
        elif sl_level>=0 and sr_level>=0:
            if sl_level>sr_level:
                dx = self.GetCellSizeX(sr_level)
                sll = self.nodes[self.sides[sl][s_nld]][l]
                assert sll>=0
                if not onesided:
                    dxSX_ind.append(sl)
                    dxSX.append(-0.5/dx)
                    dxSX_ind.append(sll)
                    dxSX.append(-0.5/dx)
                    dxSX_ind.append(sr)
                    dxSX.append(1.0/dx)
                elif onesided=='ld':
                    dxSX_ind.append(sl)
                    dxSX.append(-0.5/dx)
                    dxSX_ind.append(sll)
                    dxSX.append(-0.5/dx)
                else:
                    assert onesided=='ru'
                    dxSX_ind.append(sr)
                    dxSX.append(1.0/dx)
            else:
                dx = self.GetCellSizeX(sl_level)
                srr = self.nodes[self.sides[sr][s_nru]][r]
                assert srr>=0
                if not onesided:
                    dxSX_ind.append(sl)
                    dxSX.append(-1.0/dx)
                    dxSX_ind.append(sr)
                    dxSX.append(0.5/dx)
                    dxSX_ind.append(srr)
                    dxSX.append(0.5/dx)
                elif onesided=='ld':
                    dxSX_ind.append(sl)
                    dxSX.append(-1.0/dx)
                else:
                    assert onesided=='ru'
                    dxSX_ind.append(sr)
                    dxSX.append(0.5/dx)
                    dxSX_ind.append(srr)
                    dxSX.append(0.5/dx)
        elif sl_level>=0 and sr==SIDE_INSIDE_CELL:
            ##TODO: use higher order interpolation for the right side
            assert su>=0 and sd>=0
            dx = self.GetCellSizeX(sl_level)
            cr = self.sides[su][s_cru]  #right cell
            s_ru = self.cells[cr][u]
            s_rd = self.cells[cr][d]
            if not onesided:
                dxSX_ind.append(sl)
                dxSX.append(-2.0/3.0/dx)
                if s_ru>=0:
                    dxSX_ind.append(s_ru)
                    dxSX.append(+1.0/3.0/dx)
                else:
                    s_ru__0 = self.sidesDiv[abs(s_ru)][0]
                    s_ru__1 = self.sidesDiv[abs(s_ru)][1]
                    dxSX_ind.append(s_ru__0)
                    dxSX.append(+1.0/6.0/dx)
                    dxSX_ind.append(s_ru__1)
                    dxSX.append(+1.0/6.0/dx)
                if s_rd>=0:
                    dxSX_ind.append(s_rd)
                    dxSX.append(+1.0/3.0/dx)
                else:
                    s_rd__0 = self.sidesDiv[abs(s_rd)][0]
                    s_rd__1 = self.sidesDiv[abs(s_rd)][1]
                    dxSX_ind.append(s_rd__0)
                    dxSX.append(+1.0/6.0/dx)
                    dxSX_ind.append(s_rd__1)
                    dxSX.append(+1.0/6.0/dx)
            elif onesided=='ld':
                dxSX_ind.append(sl)
                dxSX.append(-2.0/3.0/dx)
            else:
                assert onesided=='ru'
                if s_ru>=0:
                    dxSX_ind.append(s_ru)
                    dxSX.append(1.0/3.0/dx)
                else:
                    s_ru__0 = self.sidesDiv[abs(s_ru)][0]
                    s_ru__1 = self.sidesDiv[abs(s_ru)][1]
                    dxSX_ind.append(s_ru__0)
                    dxSX.append(+1.0/6.0/dx)
                    dxSX_ind.append(s_ru__1)
                    dxSX.append(+1.0/6.0/dx)
                if s_rd>=0:
                    dxSX_ind.append(s_rd)
                    dxSX.append(+1.0/3.0/dx)
                else:
                    s_rd__0 = self.sidesDiv[abs(s_rd)][0]
                    s_rd__1 = self.sidesDiv[abs(s_rd)][1]
                    dxSX_ind.append(s_rd__0)
                    dxSX.append(+1.0/6.0/dx)
                    dxSX_ind.append(s_rd__1)
                    dxSX.append(+1.0/6.0/dx)
        elif sr_level>=0 and sl==SIDE_INSIDE_CELL:
            ##TODO: use higher order interpolation for the left side
            assert su>=0 and sd>=0
            dx = self.GetCellSizeX(sr_level)
            cl = self.sides[su][s_cld]  #left cell
            s_lu = self.cells[cl][u]
            s_ld = self.cells[cl][d]
            if not onesided:
                if s_lu>=0:
                    dxSX_ind.append(s_lu)
                    dxSX.append(-1.0/3.0/dx)
                else:
                    s_lu__0 = self.sidesDiv[abs(s_lu)][0]
                    s_lu__1 = self.sidesDiv[abs(s_lu)][1]
                    dxSX_ind.append(s_lu__0)
                    dxSX.append(-1.0/6.0/dx)
                    dxSX_ind.append(s_lu__1)
                    dxSX.append(-1.0/6.0/dx)
                if s_ld>=0:
                    dxSX_ind.append(s_ld)
                    dxSX.append(-1.0/3.0/dx)
                else:
                    s_ld__0 = self.sidesDiv[abs(s_ld)][0]
                    s_ld__1 = self.sidesDiv[abs(s_ld)][1]
                    dxSX_ind.append(s_ld__0)
                    dxSX.append(-1.0/6.0/dx)
                    dxSX_ind.append(s_ld__1)
                    dxSX.append(-1.0/6.0/dx)
                dxSX_ind.append(sr)
                dxSX.append(2.0/3.0/dx)
            elif onesided=='ld':
                if s_lu>=0:
                    dxSX_ind.append(s_lu)
                    dxSX.append(-1.0/3.0/dx)
                else:
                    s_lu__0 = self.sidesDiv[abs(s_lu)][0]
                    s_lu__1 = self.sidesDiv[abs(s_lu)][1]
                    dxSX_ind.append(s_lu__0)
                    dxSX.append(-1.0/6.0/dx)
                    dxSX_ind.append(s_lu__1)
                    dxSX.append(-1.0/6.0/dx)
                if s_ld>=0:
                    dxSX_ind.append(s_ld)
                    dxSX.append(-1.0/3.0/dx)
                else:
                    s_ld__0 = self.sidesDiv[abs(s_ld)][0]
                    s_ld__1 = self.sidesDiv[abs(s_ld)][1]
                    dxSX_ind.append(s_ld__0)
                    dxSX.append(-1.0/6.0/dx)
                    dxSX_ind.append(s_ld__1)
                    dxSX.append(-1.0/6.0/dx)
            else:
                assert onesided=='ru'
                dxSX_ind.append(sr)
                dxSX.append(2.0/3.0/dx)
        elif sl_level>=0 and sr==-1:
            ## assuming the side value outside the grid region is zero
            ##TODO: use higher order backward difference
            assert onesided==False or onesided=='ld'
            dx = self.GetCellSizeX(sl_level)
            dxSX_ind.append(sl)
            dxSX.append(-1.0/dx)
        else:
            ## assuming the side value outside the grid region is zero
            ##TODO: use higher order forward difference
            assert sr_level>=0 and sl==-1
            assert onesided==False or onesided=='ru'
            dx = self.GetCellSizeX(sr_level)
            dxSX_ind.append(sr)
            dxSX.append(1.0/dx)
        return [dxSX_ind, dxSX]

    def Node_dySY(self, node, onesided=False):
        # SY: side along Y
        l, r, d, u = range(4)
        s_nld, s_nru, s_cld, s_cru = range(4)   #sides::  s:side n:node c:cell
        sl, sr, sd, su = self.nodes[node]
        sd_level, su_level = -1, -1
        dySY_ind = []
        dySY = []
        if sd>=0:
            sd_level = self.GetSideLevel(sd)
        if su>=0:
            su_level = self.GetSideLevel(su)
        if sd_level==su_level:
            dy = self.GetCellSizeY(sd_level)
            if not onesided:
                dySY_ind.append(sd)
                dySY.append(-1.0/dy)
                dySY_ind.append(su)
                dySY.append(1.0/dy)
            elif onesided=='ld':
                dySY_ind.append(sd)
                dySY.append(-1.0/dy)
            else:
                assert onesided=='ru'
                dySY_ind.append(su)
                dySY.append(1.0/dy)
        elif sd_level>=0 and su_level>=0:
            if sd_level>su_level:
                dy = self.GetCellSizeY(su_level)
                sdd = self.nodes[self.sides[sd][s_nld]][d]
                if not onesided:
                    dySY_ind.append(sd)
                    dySY.append(-0.5/dy)
                    dySY_ind.append(sdd)
                    dySY.append(-0.5/dy)
                    dySY_ind.append(su)
                    dySY.append(1.0/dy)
                elif onesided=='ld':
                    dySY_ind.append(sd)
                    dySY.append(-0.5/dy)
                    dySY_ind.append(sdd)
                    dySY.append(-0.5/dy)
                else:
                    assert onesided=='ru'
                    dySY_ind.append(su)
                    dySY.append(1.0/dy)
            else:
                dy = self.GetCellSizeY(sd_level)
                suu = self.nodes[self.sides[su][s_nru]][u]
                if not onesided:
                    dySY_ind.append(sd)
                    dySY.append(-1.0/dy)
                    dySY_ind.append(su)
                    dySY.append(0.5/dy)
                    dySY_ind.append(suu)
                    dySY.append(0.5/dy)
                elif onesided=='ld':
                    dySY_ind.append(sd)
                    dySY.append(-1.0/dy)
                else:
                    assert onesided=='ru'
                    dySY_ind.append(su)
                    dySY.append(0.5/dy)
                    dySY_ind.append(suu)
                    dySY.append(0.5/dy)
        elif sd_level>=0 and su==SIDE_INSIDE_CELL:
            ##TODO: use higher order interpolation for the right side
            assert sr>=0 and sl>=0
            dy = self.GetCellSizeY(sd_level)
            cu = self.sides[sl][s_cru]  #up cell
            s_ur = self.cells[cu][r]
            s_ul = self.cells[cu][l]
            if not onesided:
                dySY_ind.append(sd)
                dySY.append(-2.0/3.0/dy)
                if s_ur>=0:
                    dySY_ind.append(s_ur)
                    dySY.append(1.0/3.0/dy)
                else:
                    s_ur__0 = self.sidesDiv[abs(s_ur)][0]
                    s_ur__1 = self.sidesDiv[abs(s_ur)][1]
                    dySY_ind.append(s_ur__0)
                    dySY.append(1.0/6.0/dy)
                    dySY_ind.append(s_ur__1)
                    dySY.append(1.0/6.0/dy)
                if s_ul>=0:
                    dySY_ind.append(s_ul)
                    dySY.append(1.0/3.0/dy)
                else:
                    s_ul__0 = self.sidesDiv[abs(s_ul)][0]
                    s_ul__1 = self.sidesDiv[abs(s_ul)][1]
                    dySY_ind.append(s_ul__0)
                    dySY.append(1.0/6.0/dy)
                    dySY_ind.append(s_ul__1)
                    dySY.append(1.0/6.0/dy)
            elif onesided=='ld':
                dySY_ind.append(sd)
                dySY.append(-2.0/3.0/dy)
            else:
                assert onesided=='ru'
                if s_ur>=0:
                    dySY_ind.append(s_ur)
                    dySY.append(1.0/3.0/dy)
                else:
                    s_ur__0 = self.sidesDiv[abs(s_ur)][0]
                    s_ur__1 = self.sidesDiv[abs(s_ur)][1]
                    dySY_ind.append(s_ur__0)
                    dySY.append(1.0/6.0/dy)
                    dySY_ind.append(s_ur__1)
                    dySY.append(1.0/6.0/dy)
                if s_ul>=0:
                    dySY_ind.append(s_ul)
                    dySY.append(1.0/3.0/dy)
                else:
                    s_ul__0 = self.sidesDiv[abs(s_ul)][0]
                    s_ul__1 = self.sidesDiv[abs(s_ul)][1]
                    dySY_ind.append(s_ul__0)
                    dySY.append(1.0/6.0/dy)
                    dySY_ind.append(s_ul__1)
                    dySY.append(1.0/6.0/dy)
        elif su_level>=0 and sd==SIDE_INSIDE_CELL:
            ##TODO: use higher order interpolation for the left side
            assert sr>=0 and sl>=0
            dy = self.GetCellSizeY(su_level)
            cd = self.sides[sr][s_cld]  #down cell
            s_dl = self.cells[cd][l]
            s_dr = self.cells[cd][r]
            if not onesided:
                if s_dl>=0:
                    dySY_ind.append(s_dl)
                    dySY.append(-1.0/3.0/dy)
                else:
                    s_dl__0 = self.sidesDiv[abs(s_dl)][0]
                    s_dl__1 = self.sidesDiv[abs(s_dl)][1]
                    dySY_ind.append(s_dl__0)
                    dySY.append(-1.0/6.0/dy)
                    dySY_ind.append(s_dl__1)
                    dySY.append(-1.0/6.0/dy)
                if s_dr>=0:
                    dySY_ind.append(s_dr)
                    dySY.append(-1.0/3.0/dy)
                else:
                    s_dr__0 = self.sidesDiv[abs(s_dr)][0]
                    s_dr__1 = self.sidesDiv[abs(s_dr)][1]
                    dySY_ind.append(s_dr__0)
                    dySY.append(-1.0/6.0/dy)
                    dySY_ind.append(s_dr__1)
                    dySY.append(-1.0/6.0/dy)
                dySY_ind.append(su)
                dySY.append(2.0/3.0/dy)
            elif onesided=='ld':
                if s_dl>=0:
                    dySY_ind.append(s_dl)
                    dySY.append(-1.0/3.0/dy)
                else:
                    s_dl__0 = self.sidesDiv[abs(s_dl)][0]
                    s_dl__1 = self.sidesDiv[abs(s_dl)][1]
                    dySY_ind.append(s_dl__0)
                    dySY.append(-1.0/6.0/dy)
                    dySY_ind.append(s_dl__1)
                    dySY.append(-1.0/6.0/dy)
                if s_dr>=0:
                    dySY_ind.append(s_dr)
                    dySY.append(-1.0/3.0/dy)
                else:
                    s_dr__0 = self.sidesDiv[abs(s_dr)][0]
                    s_dr__1 = self.sidesDiv[abs(s_dr)][1]
                    dySY_ind.append(s_dr__0)
                    dySY.append(-1.0/6.0/dy)
                    dySY_ind.append(s_dr__1)
                    dySY.append(-1.0/6.0/dy)
            else:
                assert onesided=='ru'
                dySY_ind.append(su)
                dySY.append(2.0/3.0/dy)
        elif sd_level>=0 and su==-1:
            ## assuming the side value outside the grid region is zero
            ##TODO: use higher order backward difference
            assert onesided==False or onesided=='ld'
            dy = self.GetCellSizeY(sd_level)
            dySY_ind.append(sd)
            dySY.append(-1.0/dy)
        else:
            ## assuming the side value outside the grid region is zero
            ##TODO: use higher order forward difference
            assert su_level>=0 and sd==-1
            assert onesided==False or onesided=='ru'
            dy = self.GetCellSizeY(su_level)
            dySY_ind.append(su)
            dySY.append(1.0/dy)
        return [dySY_ind, dySY]
        
    def Side_dxyN(self, side):
        ## dNode/dx or dNode/dy depending on the orientation of the side
        nld, nru, cld, cru, side_XY = self.sides[side]
        dxyN_ind = []
        dxyN = []
        dxy = self.GetSideSize(side)
        dxyN_ind.append(nld)
        dxyN.append(-1.0/dxy)
        dxyN_ind.append(nru)
        dxyN.append(1.0/dxy)
        return [dxyN_ind, dxyN]
        
    def Side_dxyC(self, side, onesided=False):
        ## dCell/dx or dCell/dy depending on the orientation of the side
        ## onesided = 'ld'  --> only the coefficient of the left/down cell is returned
        ## onesided = 'ru'  --> only the coefficient of the right/up cell is returned
        l, r, d, u = range(4)
        s_nld, s_nru, s_cld, s_cru = range(4)   #sides::  s:side n:node c:cell
        nld, nru, cld, cru, side_XY = self.sides[side]
        lev = 4
        dxyC_ind = []
        dxyC = []
        side_level = self.GetSideLevel(side)
        dxy = 0.0
        if side_XY=='X':
            dxy = self.GetCellSizeY(side_level)
        else:
            dxy = self.GetCellSizeX(side_level)
        cld_level, cru_level = -1, -1
        if cld>=0:
            cld_level = self.cells[cld][lev]
        if cru>=0:
            cru_level = self.cells[cru][lev]
        if cld_level==cru_level:
            if not onesided:
                dxyC_ind.append(cld)
                dxyC.append(-1.0/dxy)
                dxyC_ind.append(cru)
                dxyC.append(1.0/dxy)
            elif onesided=='ld':
                dxyC_ind.append(cld)
                dxyC.append(-1.0/dxy)
            else:
                assert onesided=='ru'
                dxyC_ind.append(cru)
                dxyC.append(1.0/dxy)
        elif cld_level<cru_level:
            if cld>=0:
                p_ld = self.Cell_GetCenterPoint(cld)
                if side_XY=='X':
                    cru_ = -1
                    assert self.cells[cld][u]<0
                    if self.sidesDiv[abs(self.cells[cld][u])][0]==side:
                        cru_ = self.sides[self.sidesDiv[abs(self.cells[cld][u])][1]][s_cru]
                    else:
                        assert self.sidesDiv[abs(self.cells[cld][u])][1]==side
                        cru_ = self.sides[self.sidesDiv[abs(self.cells[cld][u])][0]][s_cru]
                    p_ru = self.Cell_GetCenterPoint(cru)
                    p_ru_ = self.Cell_GetCenterPoint(cru_)
                    p = p_ru - Point2D(0.0, dxy)
                    lagr_coeffs = self.InterpolateCoeffs([p_ld, p_ru, p_ru_], p)
                    if not onesided:
                        dxyC_ind.append(cld)
                        dxyC.append(-lagr_coeffs[0]/dxy)
                        dxyC_ind.append(cru)
                        dxyC.append(-lagr_coeffs[1]/dxy)
                        dxyC_ind.append(cru_)
                        dxyC.append(-lagr_coeffs[2]/dxy)
                        dxyC_ind.append(cru)
                        dxyC.append(1.0/dxy)
                    elif onesided=='ld':
                        dxyC_ind.append(cld)
                        dxyC.append(-lagr_coeffs[0]/dxy)
                        dxyC_ind.append(cru)
                        dxyC.append(-lagr_coeffs[1]/dxy)
                        dxyC_ind.append(cru_)
                        dxyC.append(-lagr_coeffs[2]/dxy)
                    else:
                        assert onesided=='ru'
                        dxyC_ind.append(cru)
                        dxyC.append(1.0/dxy)
                else:
                    assert side_XY=='Y'
                    cru_ = -1
                    assert self.cells[cld][r]<0
                    if self.sidesDiv[abs(self.cells[cld][r])][0]==side:
                        cru_ = self.sides[self.sidesDiv[abs(self.cells[cld][r])][1]][s_cru]
                    else:
                        assert self.sidesDiv[abs(self.cells[cld][r])][1]==side
                        cru_ = self.sides[self.sidesDiv[abs(self.cells[cld][r])][0]][s_cru]
                    p_ru = self.Cell_GetCenterPoint(cru)
                    p_ru_ = self.Cell_GetCenterPoint(cru_)
                    p = p_ru - Point2D(dxy, 0.0)
                    lagr_coeffs = self.InterpolateCoeffs([p_ld, p_ru, p_ru_], p)
                    if not onesided:
                        dxyC_ind.append(cld)
                        dxyC.append(-lagr_coeffs[0]/dxy)
                        dxyC_ind.append(cru)
                        dxyC.append(-lagr_coeffs[1]/dxy)
                        dxyC_ind.append(cru_)
                        dxyC.append(-lagr_coeffs[2]/dxy)
                        dxyC_ind.append(cru)
                        dxyC.append(1.0/dxy)
                    elif onesided=='ld':
                        dxyC_ind.append(cld)
                        dxyC.append(-lagr_coeffs[0]/dxy)
                        dxyC_ind.append(cru)
                        dxyC.append(-lagr_coeffs[1]/dxy)
                        dxyC_ind.append(cru_)
                        dxyC.append(-lagr_coeffs[2]/dxy)
                    else:
                        assert onesided=='ru'
                        dxyC_ind.append(cru)
                        dxyC.append(1.0/dxy)
            else:
                ##TODO:cells outside are considered 0. bachward difference may be used
                assert onesided==False or onesided=='ru'
                dxyC_ind.append(cru)
                dxyC.append(1.0/dxy)
        else:
            if cru>=0:
                p_ru = self.Cell_GetCenterPoint(cru)
                if side_XY=='X':
                    cld_ = -1
                    assert self.cells[cru][d]<0
                    if self.sidesDiv[abs(self.cells[cru][d])][0]==side:
                        cld_ = self.sides[self.sidesDiv[abs(self.cells[cru][d])][1]][s_cld]
                    else:
                        assert self.sidesDiv[abs(self.cells[cru][d])][1]==side
                        cld_ = self.sides[self.sidesDiv[abs(self.cells[cru][d])][0]][s_cld]
                    p_ld = self.Cell_GetCenterPoint(cld)
                    p_ld_ = self.Cell_GetCenterPoint(cld_)
                    p = p_ld + Point2D(0.0, dxy)
                    lagr_coeffs = self.InterpolateCoeffs([p_ld, p_ld_, p_ru], p)
                    if not onesided:
                        dxyC_ind.append(cld)
                        dxyC.append(-1.0/dxy)
                        dxyC_ind.append(cld)
                        dxyC.append(lagr_coeffs[0]/dxy)
                        dxyC_ind.append(cld_)
                        dxyC.append(lagr_coeffs[1]/dxy)
                        dxyC_ind.append(cru)
                        dxyC.append(lagr_coeffs[2]/dxy)
                    elif onesided=='ld':
                        dxyC_ind.append(cld)
                        dxyC.append(-1.0/dxy)
                    else:
                        assert onesided=='ru'
                        dxyC_ind.append(cld)
                        dxyC.append(lagr_coeffs[0]/dxy)
                        dxyC_ind.append(cld_)
                        dxyC.append(lagr_coeffs[1]/dxy)
                        dxyC_ind.append(cru)
                        dxyC.append(lagr_coeffs[2]/dxy)
                else:
                    assert side_XY=='Y'
                    cld_ = -1
                    assert self.cells[cru][l]<0
                    if self.sidesDiv[abs(self.cells[cru][l])][0]==side:
                        cld_ = self.sides[self.sidesDiv[abs(self.cells[cru][l])][1]][s_cld]
                    else:
                        assert self.sidesDiv[abs(self.cells[cru][l])][1]==side
                        cld_ = self.sides[self.sidesDiv[abs(self.cells[cru][l])][0]][s_cld]
                    p_ld = self.Cell_GetCenterPoint(cld)
                    p_ld_ = self.Cell_GetCenterPoint(cld_)
                    p = p_ld + Point2D(dxy, 0.0)
                    lagr_coeffs = self.InterpolateCoeffs([p_ld, p_ld_, p_ru], p)
                    if not onesided:
                        dxyC_ind.append(cld)
                        dxyC.append(-1.0/dxy)
                        dxyC_ind.append(cld)
                        dxyC.append(lagr_coeffs[0]/dxy)
                        dxyC_ind.append(cld_)
                        dxyC.append(lagr_coeffs[1]/dxy)
                        dxyC_ind.append(cru)
                        dxyC.append(lagr_coeffs[2]/dxy)
                    elif onesided=='ld':
                        dxyC_ind.append(cld)
                        dxyC.append(-1.0/dxy)
                    else:
                        assert onesided=='ru'
                        dxyC_ind.append(cld)
                        dxyC.append(lagr_coeffs[0]/dxy)
                        dxyC_ind.append(cld_)
                        dxyC.append(lagr_coeffs[1]/dxy)
                        dxyC_ind.append(cru)
                        dxyC.append(lagr_coeffs[2]/dxy)
            else:
                ##TODO:cells outside are considered 0. bachward difference may be used
                assert onesided==False or onesided=='ld'
                dxyC_ind.append(cld)
                dxyC.append(-1.0/dxy)
        return [dxyC_ind, dxyC]
        
    def SideX_SX(self, side):
        ## Side SX (type SX) at location of side SX
        SX_ind = [side]
        SX = [1.0]
        return [SX_ind, SX]
        
    def SideY_SY(self, side):
        ## Side SY at location of side SY
        SY_ind = [side]
        SY = [1.0]
        return [SY_ind, SY]
        
    def SideX_N(self, side):
        s_nld, s_nru, s_cld, s_cru = range(4)   #sides::  s:side n:node c:cell
        nl, nr = self.sides[side][s_nld], self.sides[side][s_nru]
        N_ind = [nl, nr]
        N = [0.5, 0.5]
        return [N_ind, N]
        
    def SideY_N(self, side):
        s_nld, s_nru, s_cld, s_cru = range(4)   #sides::  s:side n:node c:cell
        nd, nu = self.sides[side][s_nld], self.sides[side][s_nru]
        N_ind = [nd, nu]
        N = [0.5, 0.5]
        return [N_ind, N]

    def NodeN_N(self, node):
        N_ind = [node]
        N = [1.0]
        return [N_ind, N]
        
    def NodeN_SX(self, node):
        l, r, d, u = range(4)
        s_nld, s_nru, s_cld, s_cru = range(4)   #sides::  s:side n:node c:cell
        SX_ind = []
        SX = []
        sl, sr, sd, su = self.nodes[node]
        sl_lev, sr_lev = self.GetSideLevel(sl), self.GetSideLevel(sr)
        if sl>=0 and sr>=0:
            if sl_lev==sr_lev:
                SX_ind = [sl, sr]
                SX = [0.5, 0.5]
            elif sl_lev>sr_lev:
                SX_ind = [sl, sr]
                SX = [2.0/3.0, 1.0/3.0]
            else:
                assert sl_lev<sr_lev
                SX_ind = [sl, sr]
                SX = [1.0/3.0, 2.0/3.0]
        elif sl<0:
            if sl==SIDE_INSIDE_CELL:
                assert sd>=0 and su>=0
                cl = self.sides[sd][s_cld]
                sld = self.cells[cl][d]
                slu = self.cells[cl][u]
                if sld>=0 and slu>=0:
                    SX_ind = [sld, slu, sr]
                    SX = [1.0/6.0, 1.0/6.0, 2.0/3.0]
                elif sld<0 and slu>=0:
                    sld0 = self.sidesDiv[abs(sld)][0]
                    sld1 = self.sidesDiv[abs(sld)][1]
                    SX_ind = [sld0, sld1, slu, sr]
                    SX = [1.0/12.0, 1.0/12.0, 1.0/6.0, 2.0/3.0]
                elif sld>=0 and slu<0:
                    slu0 = self.sidesDiv[abs(slu)][0]
                    slu1 = self.sidesDiv[abs(slu)][1]
                    SX_ind = [sld, slu0, slu1, sr]
                    SX = [1.0/6.0, 1.0/12.0, 1.0/12.0, 2.0/3.0]
                else:
                    assert sld<0 and slu<0
                    sld0 = self.sidesDiv[abs(sld)][0]
                    sld1 = self.sidesDiv[abs(sld)][1]
                    slu0 = self.sidesDiv[abs(slu)][0]
                    slu1 = self.sidesDiv[abs(slu)][1]
                    SX_ind = [sld0, sld1, slu0, slu1, sr]
                    SX = [1.0/12.0, 1.0/12.0, 1.0/12.0, 1.0/12.0, 2.0/3.0]
            else:
                assert sl==SIDE_OUTSIDE_DOMAIN
                SX_ind = [sr]
                SX = [1.0]
        else:
            assert sr<0
            if sr==SIDE_INSIDE_CELL:
                assert sd>=0 and su>=0
                cr = self.sides[sd][s_cru]
                srd = self.cells[cr][d]
                sru = self.cells[cr][u]
                if srd>=0 and sru>=0:
                    SX_ind = [sl, srd, sru]
                    SX = [2.0/3.0, 1.0/6.0, 1.0/6.0]
                elif srd<0 and sru>=0:
                    srd0 = self.sidesDiv[abs(srd)][0]
                    srd1 = self.sidesDiv[abs(srd)][1]
                    SX_ind = [sl, srd0, srd1, sru]
                    SX = [2.0/3.0, 1.0/12.0, 1.0/12.0, 1.0/6.0]
                elif srd>=0 and sru<0:
                    sru0 = self.sidesDiv[abs(sru)][0]
                    sru1 = self.sidesDiv[abs(sru)][1]
                    SX_ind = [sl, srd, sru0, sru1]
                    SX = [2.0/3.0, 1.0/6.0, 1.0/12.0, 1.0/12.0]
                else:
                    assert srd<0 and sru<0
                    srd0 = self.sidesDiv[abs(srd)][0]
                    srd1 = self.sidesDiv[abs(srd)][1]
                    sru0 = self.sidesDiv[abs(sru)][0]
                    sru1 = self.sidesDiv[abs(sru)][1]
                    SX_ind = [sl, srd0, srd1, sru0, sru1]
                    SX = [2.0/3.0, 1.0/12.0, 1.0/12.0, 1.0/12.0, 1.0/12.0]
            else:
                assert sr==SIDE_OUTSIDE_DOMAIN
                SX_ind = [sl]
                SX = [1.0]
        return [SX_ind, SX]

    def NodeN_SY(self, node):
        l, r, d, u = range(4)
        s_nld, s_nru, s_cld, s_cru = range(4)   #sides::  s:side n:node c:cell
        SY_ind = []
        SY = []
        sl, sr, sd, su = self.nodes[node]
        sd_lev, su_lev = self.GetSideLevel(sd), self.GetSideLevel(su)
        if sd>=0 and su>=0:
            if sd_lev==su_lev:
                SY_ind = [sd, su]
                SY = [0.5, 0.5]
            elif sd_lev>su_lev:
                SY_ind = [sd, su]
                SY = [2.0/3.0, 1.0/3.0]
            else:
                assert sd_lev<su_lev
                SY_ind = [sd, su]
                SY = [1.0/3.0, 2.0/3.0]
        elif sd<0:
            if sd==SIDE_INSIDE_CELL:
                assert sl>=0 and sr>=0
                cd = self.sides[sl][s_cld]
                sdl = self.cells[cd][l]
                sdr = self.cells[cd][r]
                if sdl>=0 and sdr>=0:
                    SY_ind = [sdl, sdr, su]
                    SY = [1.0/6.0, 1.0/6.0, 2.0/3.0]
                elif sdl<0 and sdr>=0:
                    sdl0 = self.sidesDiv[abs(sdl)][0]
                    sdl1 = self.sidesDiv[abs(sdl)][1]
                    SY_ind = [sdl0, sdl1, sdr, su]
                    SY = [1.0/12.0, 1.0/12.0, 1.0/6.0, 2.0/3.0]
                elif sdl>=0 and sdr<0:
                    sdr0 = self.sidesDiv[abs(sdr)][0]
                    sdr1 = self.sidesDiv[abs(sdr)][1]
                    SY_ind = [sdl, sdr0, sdr1, su]
                    SY = [1.0/6.0, 1.0/12.0, 1.0/12.0, 2.0/3.0]
                else:
                    assert sdl<0 and sdr<0
                    sdl0 = self.sidesDiv[abs(sdl)][0]
                    sdl1 = self.sidesDiv[abs(sdl)][1]
                    sdr0 = self.sidesDiv[abs(sdr)][0]
                    sdr1 = self.sidesDiv[abs(sdr)][1]
                    SY_ind = [sdl0, sdl1, sdr0, sdr1, su]
                    SY = [1.0/12.0, 1.0/12.0, 1.0/12.0, 1.0/12.0, 2.0/3.0]
            else:
                assert sd==SIDE_OUTSIDE_DOMAIN
                SY_ind = [su]
                SY = [1.0]
        else:
            assert su<0
            if su==SIDE_INSIDE_CELL:
                assert sl>=0 and sr>=0
                cu = self.sides[sl][s_cru]
                sul = self.cells[cu][l]
                sur = self.cells[cu][r]
                if sul>=0 and sur>=0:
                    SY_ind = [sd, sul, sur]
                    SY = [2.0/3.0, 1.0/6.0, 1.0/6.0]
                elif sul<0 and sur>=0:
                    sul0 = self.sidesDiv[abs(sul)][0]
                    sul1 = self.sidesDiv[abs(sul)][1]
                    SY_ind = [sd, sul0, sul1, sur]
                    SY = [2.0/3.0, 1.0/12.0, 1.0/12.0, 1.0/6.0]
                elif sul>=0 and sur<0:
                    sur0 = self.sidesDiv[abs(sur)][0]
                    sur1 = self.sidesDiv[abs(sur)][1]
                    SY_ind = [sd, sul, sur0, sur1]
                    SY = [2.0/3.0, 1.0/6.0, 1.0/12.0, 1.0/12.0]
                else:
                    assert sul<0 and sur<0
                    sul0 = self.sidesDiv[abs(sul)][0]
                    sul1 = self.sidesDiv[abs(sul)][1]
                    sur0 = self.sidesDiv[abs(sur)][0]
                    sur1 = self.sidesDiv[abs(sur)][1]
                    SY_ind = [sd, sul0, sul1, sur0, sur1]
                    SY = [2.0/3.0, 1.0/12.0, 1.0/12.0, 1.0/12.0, 1.0/12.0]
            else:
                assert su==SIDE_OUTSIDE_DOMAIN
                SY_ind = [sd]
                SY = [1.0]
        return [SY_ind, SY]

    
    def CellC_C(self, cell):
        C_ind = [cell]
        C = [1.0]
        return [C_ind, C]
        
    def Cell_IntplNodes(self, cell):
        ## interpolates the center of the cell with the corner nodes
        l, r, d, u = range(4)
        s_nld, s_nru, s_cld, s_cru = range(4)   #sides::  s:side n:node c:cell
        nld, nlu, nrd, nru = [-1]*4
        sl = self.cells[cell][l]
        if sl>=0:
            nld = self.sides[sl][s_nld]
            nlu = self.sides[sl][s_nru]
        else:
            nld = self.sides[self.sidesDiv[abs(sl)][0]][s_nld]
            nlu = self.sides[self.sidesDiv[abs(sl)][1]][s_nru]
        sr = self.cells[cell][r]
        if sr>=0:
            nrd = self.sides[sr][s_nld]
            nru = self.sides[sr][s_nru]
        else:
            nrd = self.sides[self.sidesDiv[abs(sr)][0]][s_nld]
            nru = self.sides[self.sidesDiv[abs(sr)][1]][s_nru]
        C_ind = [nld, nlu, nrd, nru]    ## corner points
        C = [0.25, 0.25, 0.25, 0.25]
        return [C_ind, C]

    def Cell_GetCornerNodes(self, cell):
        ## interpolates the center of the cell with the corner nodes
        l, r, d, u = range(4)
        s_nld, s_nru, s_cld, s_cru = range(4)   #sides::  s:side n:node c:cell
        nld, nlu, nrd, nru = [-1]*4
        sl = self.cells[cell][l]
        if sl>=0:
            nld = self.sides[sl][s_nld]
            nlu = self.sides[sl][s_nru]
        else:
            nld = self.sides[self.sidesDiv[abs(sl)][0]][s_nld]
            nlu = self.sides[self.sidesDiv[abs(sl)][1]][s_nru]
        sr = self.cells[cell][r]
        if sr>=0:
            nrd = self.sides[sr][s_nld]
            nru = self.sides[sr][s_nru]
        else:
            nrd = self.sides[self.sidesDiv[abs(sr)][0]][s_nld]
            nru = self.sides[self.sidesDiv[abs(sr)][1]][s_nru]
        C_ind = [nld, nlu, nrd, nru]    ## corner points
        return C_ind


    def Cell_GetCenterPoint(self, cell):
        ## calculates the coordinates of the center of the cell
        C_ind, _ = self.Cell_IntplNodes(cell)
        p = Point2D(0.0, 0.0)
        for i in range(4):
            p += self.nodesPoints[C_ind[i]]
        return p/4.0
        
    def InterpolateCoeffs(self, ri, r):
        assert len(ri)==3
        return barycentricInterpolation_coeffs(ri, r)

    def AddExpression(self, region_key, expr_lst, var_lst):
        """ adds {region_key: expr} to self.Regions dictionary, region_key is a string
        describing the equation
        expr is a list of sympy expressions for a system of equations describing
        the same region of the computational domain
        var is a list of sympy variables
        'region_key': [[eq0, eq1...eqn], [x0, x1,...xn]]
        each equation is associated with a given variable (the derivatives are 
        taken with respect to the center of that variable)
        """
        self.Regions[region_key] = [expr_lst, var_lst]
        return
        
    def AssignVariableToCells(self, var, cells, region_key):
        ## var=sympy.Symbol    cells: list of cells indices
        if region_key in self.Vars:
            self.Vars[region_key].append([var, cells, 'C'])
        else:
            self.Vars[region_key] = [[var, cells, 'C']]
        return
        
    def AssignVariableToNodes(self, var, nodes, region_key):
        ## var=sympy.Symbol    nodes: list of nodes indices
        if region_key in self.Vars:
            self.Vars[region_key].append([var, nodes, 'N'])
        else:
            self.Vars[region_key] = [[var, nodes, 'N']]
        return
        
    def AssignVariableToSides(self, var, sides, region_key, orient='X'):
        ## for orient=='X' variable is assigned only to x directed sides inside
        ## sides list. similar for 'Y'
        if orient=='X':
            if region_key in self.Vars:
                self.Vars[region_key].append([var, self.GetSidesXDirected(sides), 'SX'])
            else:
                self.Vars[region_key] = [[var, self.GetSidesXDirected(sides), 'SX']]
        elif orient=='Y':
            if region_key in self.Vars:
                self.Vars[region_key].append([var, self.GetSidesYDirected(sides), 'SY'])
            else:
                self.Vars[region_key] = [[var, self.GetSidesYDirected(sides), 'SY']]
        return
        
    def _GetVariablesAndAssociatedElements(self):
        """ var_type = {[var:['C'/'N'/'S', n]]} where n is its index
        the first 'C' has indice 0, the sencond 1 ... if more than one variable
        is associated with cells for example
        """
        ##TODO: The assignment of indices to types (C/N/S) with more than one variable
        ## might be different on a second call (to be checked)
        var_type = {}
        for reg in self.Regions:
            vars_ = self.Vars[reg]
            for v in vars_:
                var_type[v[0]] = [v[2], 0]
        n_vars = len(var_type)
        n_C, n_N, n_SX, n_SY = [0]*4
        for v_, t_ in var_type.items():
            if t_[0]=='C':
                n_C += 1
                t_[1] = n_C - 1
            if t_[0]=='N':
                n_N += 1
                t_[1] = n_N - 1
            if t_[0]=='SX':
                n_SX += 1
                t_[1] = n_SX - 1
            if t_[0]=='SY':
                n_SY += 1
                t_[1] = n_SY - 1
        n_vars_dic = dict(zip(['C', 'N', 'SX', 'SY'], [n_C, n_N, n_SX, n_SY]))
        self.VarTypeDics = [var_type, n_vars_dic]
        return self.VarTypeDics

    def SetTotalIndicesForVariables(self, SameBoundaryNumber=True):
        """ Arranges all variables and gives them unique index numbers to be used
        in a matrix representation
        It also identifies the shared cells, sides or nodes between different regions
        SameBoundaryNumber==True --> the shared elements will be given the same 
        index, otherwise different
        """
        var_type_dic, n_vars_dic = self._GetVariablesAndAssociatedElements()
        #print('var_type_dic: ', var_type_dic)
        self.CellsToIndexTotal = [[None]*len(self.cells) for i in range(n_vars_dic['C'])]
        self.NodesToIndexTotal = [[None]*len(self.nodes) for i in range(n_vars_dic['N'])]
        self.SidesToIndexTotal = [[None]*len(self.sides) for i in range(max(n_vars_dic['SX'],n_vars_dic['SY']))]
        self.CellsToIndexTotalShared = [-1]
        self.NodesToIndexTotalShared = [-1]
        self.SidesToIndexTotalShared = [-1]
        ind_tot = 0
        regions_endIndex = []
        regions_names = []
        for reg in self.Regions:
            regions_names.append(reg)
            regions_endIndex.append(ind_tot)
            vars_ind = self.Vars[reg]
            for var in vars_ind:
                v = var[0]
                v_inds = var[1]
                v_type = var[2]
                v_cns_ = var_type_dic[v][1]
                if v_type == 'C':                
                    for c in v_inds:
                        if self.CellsToIndexTotal[v_cns_][c]==None:
                            self.CellsToIndexTotal[v_cns_][c] = ind_tot
                            ind_tot += 1
                        elif self.CellsToIndexTotal[v_cns_][c]>=0:
                            ind_tot_shared = self.CellsToIndexTotal[v_cns_][c]
                            reg_shared = None
                            for i in range(len(regions_names)-1):
                                if ind_tot_shared>=regions_endIndex[i] and ind_tot_shared<regions_endIndex[i+1]:
                                    reg_shared = regions_names[i]
                                    break
                            assert reg_shared!=None
                            if SameBoundaryNumber:
                                self.CellsToIndexTotal[v_cns_][c] = -len(self.CellsToIndexTotalShared)
                                self.CellsToIndexTotalShared.append([c, (reg_shared, ind_tot_shared), (reg, ind_tot_shared)])
                            else:
                                self.CellsToIndexTotal[v_cns_][c] = -len(self.CellsToIndexTotalShared)
                                self.CellsToIndexTotalShared.append([c, (reg_shared, ind_tot_shared), (reg, ind_tot)])
                                ind_tot += 1
                        else:
                            assert self.CellsToIndexTotal[v_cns_][c]<0
                            if SameBoundaryNumber:
                                ind_tot_shared = self.CellsToIndexTotalShared[abs(self.CellsToIndexTotal[v_cns_][c])][1][1]
                                self.CellsToIndexTotalShared[abs(self.CellsToIndexTotal[v_cns_][c])].append((reg, ind_tot_shared))
                            else:
                                self.CellsToIndexTotalShared[abs(self.CellsToIndexTotal[v_cns_][c])].append((reg, ind_tot))
                                ind_tot += 1
                elif v_type == 'SX' or v_type == 'SY':
                    for s in v_inds:
                        if self.SidesToIndexTotal[v_cns_][s]==None:
                            self.SidesToIndexTotal[v_cns_][s] = ind_tot
                            ind_tot += 1
                        elif self.SidesToIndexTotal[v_cns_][s]>=0:
                            ind_tot_shared = self.SidesToIndexTotal[v_cns_][s]
                            reg_shared = None
                            for i in range(len(regions_names)-1):
                                if ind_tot_shared>=regions_endIndex[i] and ind_tot_shared<regions_endIndex[i+1]:
                                    reg_shared = regions_names[i]
                                    break
                            assert reg_shared!=None
                            if SameBoundaryNumber:
                                self.SidesToIndexTotal[v_cns_][s] = -len(self.SidesToIndexTotalShared)
                                self.SidesToIndexTotalShared.append([s, (reg_shared, ind_tot_shared), (reg, ind_tot_shared)])
                            else:
                                self.SidesToIndexTotal[v_cns_][s] = -len(self.SidesToIndexTotalShared)
                                self.SidesToIndexTotalShared.append([s, (reg_shared, ind_tot_shared), (reg, ind_tot)])
                                ind_tot += 1
                        else:
                            assert self.SidesToIndexTotal[v_cns_][s]<0
                            if SameBoundaryNumber:
                                ind_tot_shared = self.SidesToIndexTotalShared[abs(self.SidesToIndexTotal[v_cns_][s])][1][1]
                                self.SidesToIndexTotalShared[abs(self.SidesToIndexTotal[v_cns_][s])].append((reg, ind_tot_shared))
                            else:
                                self.SidesToIndexTotalShared[abs(self.SidesToIndexTotal[v_cns_][s])].append((reg, ind_tot))
                                ind_tot += 1
                else:
                    assert v_type == 'N'
                    for n in v_inds:
                        if self.NodesToIndexTotal[v_cns_][n]==None:
                            self.NodesToIndexTotal[v_cns_][n] = ind_tot
                            ind_tot += 1
                        elif self.NodesToIndexTotal[v_cns_][n]>=0:
                            ind_tot_shared = self.NodesToIndexTotal[v_cns_][n]
                            reg_shared = None
                            for i in range(len(regions_names)-1):
                                if ind_tot_shared>=regions_endIndex[i] and ind_tot_shared<regions_endIndex[i+1]:
                                    reg_shared = regions_names[i]
                                    break
                            assert reg_shared!=None
                            if SameBoundaryNumber:
                                self.NodesToIndexTotal[v_cns_][n] = -len(self.NodesToIndexTotalShared)
                                self.NodesToIndexTotalShared.append([n, (reg_shared, ind_tot_shared), (reg, ind_tot_shared)])
                            else:
                                self.NodesToIndexTotal[v_cns_][n] = -len(self.NodesToIndexTotalShared)
                                self.NodesToIndexTotalShared.append([n, (reg_shared, ind_tot_shared), (reg, ind_tot)])
                                ind_tot += 1
                        else:
                            assert self.NodesToIndexTotal[v_cns_][n]<0
                            if SameBoundaryNumber:
                                ind_tot_shared = self.NodesToIndexTotalShared[abs(self.NodesToIndexTotal[v_cns_][n])][1][1]
                                self.NodesToIndexTotalShared[abs(self.NodesToIndexTotal[v_cns_][n])].append((reg, ind_tot_shared))
                            else:
                                self.NodesToIndexTotalShared[abs(self.NodesToIndexTotal[v_cns_][n])].append((reg, ind_tot))
                                ind_tot += 1
        return
        
        
    def GetCellIndexToIndexTotal(self, var, cell, region_key=None):
        var_type_dic, n_vars_dic = self.VarTypeDics
        ind_tot = self.CellsToIndexTotal[var_type_dic[var][1]][cell]
        if ind_tot>=0:
            return ind_tot
        elif region_key==None:
            return self.CellsToIndexTotalShared[abs(ind_tot)]
        else:
            ind_tot_shared = self.CellsToIndexTotalShared[abs(ind_tot)]
            for i in range(1, len(ind_tot_shared)):
                if ind_tot_shared[i][0]==region_key:
                    return ind_tot_shared[i][1]
        print(self.CellsToIndexTotalShared[abs(ind_tot)])
        raise ValueError("index total not found for the given region_key")
        
    def GetNodeIndexToIndexTotal(self, var, node, region_key=None):
        var_type_dic, n_vars_dic = self.VarTypeDics
        ind_tot = self.NodesToIndexTotal[var_type_dic[var][1]][node]
        if ind_tot>=0:
            return ind_tot
        elif region_key==None:
            return self.NodesToIndexTotalShared[abs(ind_tot)]
        else:
            ind_tot_shared = self.NodesToIndexTotalShared[abs(ind_tot)]
            for i in range(1, len(ind_tot_shared)):
                if ind_tot_shared[i][0]==region_key:
                    return ind_tot_shared[i][1]
        print('region_key: ', region_key)
        print(self.NodesToIndexTotalShared[abs(ind_tot)])
        raise ValueError("index total not found for the given region_key")

    def GetSideIndexToIndexTotal(self, var, side, region_key=None):
        var_type_dic, n_vars_dic = self.VarTypeDics
        ind_tot = self.SidesToIndexTotal[var_type_dic[var][1]][side]
        if ind_tot>=0:
            return ind_tot
        elif region_key==None:
            return self.SidesToIndexTotalShared[abs(ind_tot)]
        else:
            ind_tot_shared = self.SidesToIndexTotalShared[abs(ind_tot)]
            for i in range(1, len(ind_tot_shared)):
                if ind_tot_shared[i][0]==region_key:
                    return ind_tot_shared[i][1]
        print(self.SidesToIndexTotalShared[abs(ind_tot)])
        raise ValueError("index total not found for the given region_key")
            
    def _GetVariableIndexAssociatedWithRegionsVarsAndType(self, region_name):
        """ looks at self.Vars[region_name] variables and finds the variable
        index associated with each variable in self.Regions[region_name][1]
        plus it returns its type also 'C'/'N'/..
        """
        reg_vars = self.Regions[region_name][1]
        var_varsList = self.Vars[region_name]
        assocInd = [-1]*len(reg_vars)
        for i in range(len(reg_vars)):
            v_i = reg_vars[i]
            for j in range(len(var_varsList)):
                v_j = var_varsList[j][0]
                if v_i==v_j:
                    assert assocInd[i]==-1
                    assocInd[i] = [j, var_varsList[j][2]]   ## [ind, type]
        for i in assocInd:
            assert i!=-1
        return assocInd
        
    def _GetTotalNumberOfVariables(self):
        n_total = 0
        for reg_name, eqs_vars_list in self.Regions.items():
            RegVars = self.Vars[reg_name]
            for i_v in range(len(RegVars)):
                var_inds = RegVars[i_v][1]
                n_total += len(var_inds)
        return n_total
        
    def GetMarkedSharedCells(self):
        cells_marked = [0]*len(self.cells)
        for i in range(1, len(self.CellsToIndexTotalShared)):
            cells_marked[self.CellsToIndexTotalShared[i][0]] = 1
        return cells_marked

    def GetMarkedSharedSides(self):
        sides_marked = [0]*len(self.sides)
        for i in range(1, len(self.SidesToIndexTotalShared)):
            sides_marked[self.SidesToIndexTotalShared[i][0]] = 1
        return sides_marked

    def GetMarkedSharedNodes(self):
        nodes_marked = [0]*len(self.nodes)
        for i in range(1, len(self.NodesToIndexTotalShared)):
            nodes_marked[self.NodesToIndexTotalShared[i][0]] = 1
        return nodes_marked
        
    def GetMarkedRegionCells(self, region_key):
        cells_marked = [0]*len(self.cells)
        for var_ind_type in self.Vars[region_key]:   #([var, cells, 'C'])
            _, inds, v_type = var_ind_type
            if v_type=='C':
                for c in inds:
                    cells_marked[c] = 1
        return cells_marked
                    
    def GetMarkedRegionSides(self, region_key):
        sides_marked = [0]*len(self.sides)
        for var_ind_type in self.Vars[region_key]:
            _, inds, v_type = var_ind_type
            if v_type=='SX' or v_type=='SY':
                for s in inds:
                    sides_marked[s] = 1
        return sides_marked

    def GetMarkedRegionNodes(self, region_key):
        nodes_marked = [0]*len(self.nodes)
        for var_ind_type in self.Vars[region_key]:
            _, inds, v_type = var_ind_type
            if v_type=='N':
                for n in inds:
                    nodes_marked[n] = 1
        return nodes_marked

    def ConstructInitialMatrix(self, save=False):
        cells_shared_marked = self.GetMarkedSharedCells()
        sides_shared_marked = self.GetMarkedSharedSides()
        nodes_shared_marked = self.GetMarkedSharedNodes()
        row = []
        col = []
        data = []
        row_ind = 0
        for reg_name, eqs_vars_list in self.Regions.items():
            EQs = eqs_vars_list[0]
            EqVars = eqs_vars_list[1]
            RegVars = self.Vars[reg_name]
            selfRegVar_types = self._GetVariableIndexAssociatedWithRegionsVarsAndType(reg_name)
            #print(selfRegVar_types)
            for i_eq in range(len(EQs)):
                eq = EQs[i_eq]
                vcen = EqVars[i_eq]     #center variable
                vcen_type = selfRegVar_types[i_eq][1]
                vcen_ind = selfRegVar_types[i_eq][0]    #index inside self.Vars[reg_name]
                vcen_varInds = RegVars[vcen_ind][1]     #elements (node,cell..) associated with the center variable
                for i_v in range(len(RegVars)):
                    var = RegVars[i_v][0]
                    var_inds = RegVars[i_v][1]
                    var_type = RegVars[i_v][2]
                    if eq.has(var):
                        eq_tree = symExpr_generate_tree(eq)
                        while eq_tree[1].has(var):
                            coeff_sym, der_ord, der_vars = sym_getcoeff_setzero(eq_tree, var)
                            coeff = complex(coeff_sym)
                            if der_ord==0:
                                if vcen_type=='SX' and var_type=='SX':
                                    for s in vcen_varInds:
                                        if sides_shared_marked[s]==0:
                                            row_ind = self.GetSideIndexToIndexTotal(vcen, s, reg_name)
                                            ind_init, vals = self.SideX_SX(s)
                                            for i in range(len(ind_init)):
                                                ind_tot = self.GetSideIndexToIndexTotal(var, ind_init[i], reg_name)
                                                row.append(row_ind)
                                                col.append(ind_tot)
                                                data.append(coeff*vals[i])
                                elif vcen_type=='SY' and var_type=='SY':
                                    for s in vcen_varInds:
                                        if sides_shared_marked[s]==0:
                                            row_ind = self.GetSideIndexToIndexTotal(vcen, s, reg_name)
                                            ind_init, vals = self.SideY_SY(s)
                                            for i in range(len(ind_init)):
                                                ind_tot = self.GetSideIndexToIndexTotal(var, ind_init[i], reg_name)
                                                row.append(row_ind)
                                                col.append(ind_tot)
                                                data.append(coeff*vals[i])
                                elif vcen_type=='N' and var_type=='N':
                                    for n in vcen_varInds:
                                        if nodes_shared_marked[n]==0:
                                            row_ind = self.GetNodeIndexToIndexTotal(vcen, n, reg_name)
                                            ind_init, vals = self.NodeN_N(n)
                                            for i in range(len(ind_init)):
                                                ind_tot = self.GetNodeIndexToIndexTotal(var, ind_init[i], reg_name)
                                                row.append(row_ind)
                                                col.append(ind_tot)
                                                data.append(coeff*vals[i])
                                elif vcen_type=='C' and var_type=='C':
                                    for c in vcen_varInds:
                                        if cells_shared_marked[c]==0:
                                            row_ind = self.GetCellIndexToIndexTotal(vcen, c, reg_name)
                                            ind_init, vals = self.CellC_C(c)
                                            for i in range(len(ind_init)):
                                                ind_tot = self.GetCellIndexToIndexTotal(var, ind_init[i], reg_name)
                                                row.append(row_ind)
                                                col.append(ind_tot)
                                                data.append(coeff*vals[i])
                                elif vcen_type=='SX' and var_type=='N':
                                    for s in vcen_varInds:
                                        if sides_shared_marked[s]==0:
                                            row_ind = self.GetSideIndexToIndexTotal(vcen, s, reg_name)
                                            ind_init, vals = self.SideX_N(s)
                                            for i in range(len(ind_init)):
                                                ind_tot = self.GetNodeIndexToIndexTotal(var, ind_init[i], reg_name)
                                                row.append(row_ind)
                                                col.append(ind_tot)
                                                data.append(coeff*vals[i])
                                elif vcen_type=='SY' and var_type=='N':
                                    for s in vcen_varInds:
                                        if sides_shared_marked[s]==0:
                                            row_ind = self.GetSideIndexToIndexTotal(vcen, s, reg_name)
                                            ind_init, vals = self.SideY_N(s)
                                            for i in range(len(ind_init)):
                                                ind_tot = self.GetNodeIndexToIndexTotal(var, ind_init[i], reg_name)
                                                row.append(row_ind)
                                                col.append(ind_tot)
                                                data.append(coeff*vals[i])
                                elif vcen_type=='N' and var_type=='SX':
                                    for n in vcen_varInds:
                                        if nodes_shared_marked[n]==0:
                                            row_ind = self.GetNodeIndexToIndexTotal(vcen, n, reg_name)
                                            ind_init, vals = self.NodeN_SX(n)
                                            for i in range(len(ind_init)):
                                                ind_tot = self.GetSideIndexToIndexTotal(var, ind_init[i], reg_name)
                                                row.append(row_ind)
                                                col.append(ind_tot)
                                                data.append(coeff*vals[i])
                                elif vcen_type=='N' and var_type=='SY':
                                    for n in vcen_varInds:
                                        if nodes_shared_marked[n]==0:
                                            row_ind = self.GetNodeIndexToIndexTotal(vcen, n, reg_name)
                                            ind_init, vals = self.NodeN_SY(n)
                                            for i in range(len(ind_init)):
                                                ind_tot = self.GetSideIndexToIndexTotal(var, ind_init[i], reg_name)
                                                row.append(row_ind)
                                                col.append(ind_tot)
                                                data.append(coeff*vals[i])
                                else:
                                    print('vcen_type:', vcen_type, 'var_type:', var_type)
                                    raise NotImplementedError("Not implemented")
                            elif der_ord==1:
                                if vcen_type=='C' and var_type=='SX':
                                    assert der_vars==[self.y]
                                    for c in vcen_varInds:
                                        if cells_shared_marked[c]==0:
                                            row_ind = self.GetCellIndexToIndexTotal(vcen, c, reg_name)
                                            ind_init, vals = self.Cell_dySX(c)
                                            for i in range(len(ind_init)):
                                                ind_tot = self.GetSideIndexToIndexTotal(var, ind_init[i], reg_name)
                                                row.append(row_ind)
                                                col.append(ind_tot)
                                                data.append(coeff*vals[i])
                                elif vcen_type=='C' and var_type=='SY':
                                    assert der_vars==[self.x]
                                    for c in vcen_varInds:
                                        if cells_shared_marked[c]==0:
                                            row_ind = self.GetCellIndexToIndexTotal(vcen, c, reg_name)
                                            ind_init, vals = self.Cell_dxSY(c)
                                            for i in range(len(ind_init)):
                                                ind_tot = self.GetSideIndexToIndexTotal(var, ind_init[i], reg_name)
                                                row.append(row_ind)
                                                col.append(ind_tot)
                                                data.append(coeff*vals[i])
                                elif vcen_type=='N' and var_type=='SX':
                                    assert der_vars==[self.x]
                                    for n in vcen_varInds:
                                        if nodes_shared_marked[n]==0:
                                            row_ind = self.GetNodeIndexToIndexTotal(vcen, n, reg_name)
                                            ind_init, vals = self.Node_dxSX(n)
                                            for i in range(len(ind_init)):
                                                ind_tot = self.GetSideIndexToIndexTotal(var, ind_init[i], reg_name)
                                                row.append(row_ind)
                                                col.append(ind_tot)
                                                data.append(coeff*vals[i])
                                elif vcen_type=='N' and var_type=='SY':
                                    assert der_vars==[self.y]
                                    for n in vcen_varInds:
                                        if nodes_shared_marked[n]==0:
                                            row_ind = self.GetNodeIndexToIndexTotal(vcen, n, reg_name)
                                            ind_init, vals = self.Node_dySY(n)
                                            for i in range(len(ind_init)):
                                                ind_tot = self.GetSideIndexToIndexTotal(var, ind_init[i], reg_name)
                                                row.append(row_ind)
                                                col.append(ind_tot)
                                                data.append(coeff*vals[i])
                                elif vcen_type=='SX' and var_type=='N':
                                    assert der_vars==[self.x]
                                    for s in vcen_varInds:
                                        if sides_shared_marked[s]==0:
                                            row_ind = self.GetSideIndexToIndexTotal(vcen, s, reg_name)
                                            ind_init, vals = self.Side_dxyN(s)
                                            for i in range(len(ind_init)):
                                                ind_tot = self.GetNodeIndexToIndexTotal(var, ind_init[i], reg_name)
                                                row.append(row_ind)
                                                col.append(ind_tot)
                                                data.append(coeff*vals[i])
                                elif vcen_type=='SY' and var_type=='N':
                                    assert der_vars==[self.y]
                                    for s in vcen_varInds:
                                        if sides_shared_marked[s]==0:
                                            row_ind = self.GetSideIndexToIndexTotal(vcen, s, reg_name)
                                            ind_init, vals = self.Side_dxyN(s)
                                            for i in range(len(ind_init)):
                                                ind_tot = self.GetNodeIndexToIndexTotal(var, ind_init[i], reg_name)
                                                row.append(row_ind)
                                                col.append(ind_tot)
                                                data.append(coeff*vals[i])
                                elif vcen_type=='SX' and var_type=='C':
                                    assert der_vars==[self.y]
                                    for s in vcen_varInds:
                                        if sides_shared_marked[s]==0:
                                            row_ind = self.GetSideIndexToIndexTotal(vcen, s, reg_name)
                                            ind_init, vals = self.Side_dxyC(s)
                                            for i in range(len(ind_init)):
                                                ind_tot = self.GetCellIndexToIndexTotal(var, ind_init[i], reg_name)
                                                row.append(row_ind)
                                                col.append(ind_tot)
                                                data.append(coeff*vals[i])
                                elif vcen_type=='SY' and var_type=='C':
                                    assert der_vars==[self.x]
                                    for s in vcen_varInds:
                                        if sides_shared_marked[s]==0:
                                            row_ind = self.GetSideIndexToIndexTotal(vcen, s, reg_name)
                                            ind_init, vals = self.Side_dxyC(s)
                                            for i in range(len(ind_init)):
                                                ind_tot = self.GetCellIndexToIndexTotal(var, ind_init[i], reg_name)
                                                row.append(row_ind)
                                                col.append(ind_tot)
                                                data.append(coeff*vals[i])
                                else:
                                    raise NotImplementedError("Not implemented")
                            else:
                                raise NotImplementedError("Not implemented")
        row = np.array(row)
        col = np.array(col)
        data = np.array(data)
        n_total = self._GetTotalNumberOfVariables()
        A_coo = coo_matrix((data, (row, col)), shape=(n_total,n_total), dtype=complex)
        if save==True:
            self.initialMatrix_coo = A_coo      #@
        return A_coo
        
    def InitRegionBoundaryEquations(self):
        self.RegionBoundary = {}        #@
        return
        
    def SetRegionBoundaryEquations(self, region_key, eq_list, var_list, var_type):
        """ the eq_list element is a dictionary
        for var_type=='SX'/'SY':
        eq_list[i] = {'zero': eq_i_0, 'one':eq_i_1, 'two':eq_i_2}
        eq_i_1 to be used if the side covers 1 cell
        eq_i_2 to be used if the side covers 2 cells (the two cells on its sides
        belong to the same region)
        for var_type=='N' :
        eq_list[i] = {'zero': eq_i_0,'one':eq_i_1, 'two':eq_i_2, 'three':eq_i_3, 'four':eq_i_4}
        eq_i_1 to be used if the node is on corner and covers 1 cell
        eq_i_2 to be used if the node is on the side of region and covers 2 cell
        eq_i_3 and eq_i_4 if the node covers 3 or 4 cells (the 4 cells around it 
        belong to the same region)
        for var_type=='C' :
        nomally not shared (exceptions to be handled in the future versions)
        """
        self.RegionBoundary[region_key] = [eq_list, var_list, var_type]
        return

    def _Side_GetOneSidedness(self, s, reg_cells_marked):
        s_nld, s_nru, s_cld, s_cru = range(4)   #sides::  s:side n:node c:cell
        cld = self.sides[s][s_cld]
        cru = self.sides[s][s_cru]
        onesided_ld, onesided_ru = False, False
        if cld>=0 and reg_cells_marked[cld]==1:
            onesided_ld = True
        if cru>=0 and reg_cells_marked[cru]==1:
            onesided_ru = True
        onesided_ = None
        n_cell_cov = 0
        if onesided_ld==True and onesided_ru==True:
            onesided_ = False
            n_cell_cov = 2
        elif onesided_ld==True:
            onesided_ = 'ld'
            n_cell_cov = 1
        elif onesided_ru==True:
            onesided_ = 'ru'
            n_cell_cov = 1
        else:
            onesided_ = False
            n_cell_cov = 0
            #raise ValueError("onesided_ld==False and onesided_ru==False")
        return [onesided_, n_cell_cov]
    
    def _Node_GetOneSidedness(self, n, reg_sides_marked, reg_cells_marked):
        l, r, d, u = range(4)   #nodes_sides, sides_sides_outside and self.nodes indices
        s_nld, s_nru, s_cld, s_cru = range(4)   #sides::  s:side n:node c:cell
        sl, sr, sd, su = self.nodes[n]
        onesided_l, onesided_r, onesided_d, onesided_u = False, False, False, False
        if (sl>=0 and reg_sides_marked[sl]==1):
            onesided_l = True
        elif sl==SIDE_INSIDE_CELL:
            if reg_cells_marked[self.sides[sd][s_cld]]==1:
                onesided_l = True
        if (sr>=0 and reg_sides_marked[sr]==1):
            onesided_r = True
        elif sr==SIDE_INSIDE_CELL:
            if reg_cells_marked[self.sides[sd][s_cru]]==1:
                onesided_r = True
        if (sd>=0 and reg_sides_marked[sd]==1):
            onesided_d = True
        elif sd==SIDE_INSIDE_CELL:
            if reg_cells_marked[self.sides[sl][s_cld]]==1:
                onesided_d = True
        if (su>=0 and reg_sides_marked[su]==1):
            onesided_u = True
        elif su==SIDE_INSIDE_CELL:
            if reg_cells_marked[self.sides[sl][s_cru]]==1:
                onesided_u = True
        onesided_lr, onesided_du = None, None
        if onesided_l==True and onesided_r==True:
            onesided_lr = False
        elif onesided_l==True:
            onesided_lr = 'ld'
        elif onesided_r==True:
            onesided_lr = 'ru'
        else:
            raise ValueError("onesided_ld==False and onesided_ru==False")
        if onesided_d==True and onesided_u==True:
            onesided_du = False
        elif onesided_d==True:
            onesided_du = 'ld'
        elif onesided_u==True:
            onesided_du = 'ru'
        else:
            raise ValueError("onesided_ld==False and onesided_ru==False")
        return [onesided_lr, onesided_du]
        
    def Node_GetNumRegionCellsCovered(self, node, region_cells_marked):
        c_nb = self.Node_GetNeighboringCells(node)
        n_c = 0
        for c in c_nb:
            if c>=0:
                if region_cells_marked[c]==1:
                    n_c += 1
        return n_c

    def ConstructInitialMatrix_Shared(self, A_coo, save=False):
        cells_shared_marked = self.GetMarkedSharedCells()
        sides_shared_marked = self.GetMarkedSharedSides()
        nodes_shared_marked = self.GetMarkedSharedNodes()
        row, col, data = list(A_coo.row), list(A_coo.col), list(A_coo.data)
        n_total = A_coo.shape[0]
        row_ind = len(row)
        n_cell_cov_to_str_dic = {'zero':0, 'one':1, 'two':2, 'three':3, 'four':4, 'two_x':2, 'two_y':2}
        for reg_name, eqs_vars_list in self.RegionBoundary.items():
            EQs = eqs_vars_list[0]
            EqVars = eqs_vars_list[1]
            EqVars_types = eqs_vars_list[2]
            RegVars = self.Vars[reg_name]
            selfRegVar_types = self._GetVariableIndexAssociatedWithRegionsVarsAndType(reg_name)
            reg_cells_marked = self.GetMarkedRegionCells(reg_name)
            reg_sides_marked = self.GetMarkedRegionSides(reg_name)
            for i_eq in range(len(EQs)):
                eq_dic = EQs[i_eq]
                vcen = EqVars[i_eq]     #center variable
                vcen_type = EqVars_types[i_eq]
                vcen_type__ = selfRegVar_types[i_eq][1]
                assert vcen_type==vcen_type__
                vcen_ind = selfRegVar_types[i_eq][0]    #index inside self.Vars[reg_name]
                vcen_varInds = RegVars[vcen_ind][1]     #elements (node,cell..) associated with the center variable
                for i_v in range(len(RegVars)):
                    var = RegVars[i_v][0]
                    var_inds = RegVars[i_v][1]
                    var_type = RegVars[i_v][2]
                    for eq_key, eq in eq_dic.items():
                        if eq.has(var):
                            eq_tree = symExpr_generate_tree(eq)
                            while eq_tree[1].has(var):
                                coeff_sym, der_ord, der_vars = sym_getcoeff_setzero(eq_tree, var)
                                coeff = complex(coeff_sym)
                                if der_ord==0:
                                    if vcen_type=='SX' and var_type=='SX':
                                        for s in vcen_varInds:
                                            if sides_shared_marked[s]==1:
                                                onesided_, n_cell_cov = self._Side_GetOneSidedness(s, reg_cells_marked)
                                                eq_key_num = n_cell_cov_to_str_dic[eq_key]
                                                if n_cell_cov!=eq_key_num:
                                                    continue
                                                row_ind = self.GetSideIndexToIndexTotal(vcen, s, reg_name)
                                                ind_init, vals = self.SideX_SX(s)
                                                for i in range(len(ind_init)):
                                                    ind_tot = self.GetSideIndexToIndexTotal(var, ind_init[i], reg_name)
                                                    row.append(row_ind)
                                                    col.append(ind_tot)
                                                    data.append(coeff*vals[i])
                                    elif vcen_type=='SY' and var_type=='SY':
                                        for s in vcen_varInds:
                                            if sides_shared_marked[s]==1:
                                                onesided_, n_cell_cov = self._Side_GetOneSidedness(s, reg_cells_marked)
                                                eq_key_num = n_cell_cov_to_str_dic[eq_key]
                                                if n_cell_cov!=eq_key_num:
                                                    continue
                                                row_ind = self.GetSideIndexToIndexTotal(vcen, s, reg_name)
                                                ind_init, vals = self.SideY_SY(s)
                                                for i in range(len(ind_init)):
                                                    ind_tot = self.GetSideIndexToIndexTotal(var, ind_init[i], reg_name)
                                                    row.append(row_ind)
                                                    col.append(ind_tot)
                                                    data.append(coeff*vals[i])
                                    elif vcen_type=='N' and var_type=='N':
                                        for n in vcen_varInds:
                                            if nodes_shared_marked[n]==1:
                                                n_cell_cov = self.Node_GetNumRegionCellsCovered(n, reg_cells_marked)
                                                eq_key_num = n_cell_cov_to_str_dic[eq_key]
                                                if n_cell_cov!=eq_key_num:
                                                    continue
                                                if eq_key=='two_x' or eq_key=='two_y':
                                                    continue
                                                row_ind = self.GetNodeIndexToIndexTotal(vcen, n, reg_name)
                                                ind_init, vals = self.NodeN_N(n)
                                                for i in range(len(ind_init)):
                                                    ind_tot = self.GetNodeIndexToIndexTotal(var, ind_init[i], reg_name)
                                                    row.append(row_ind)
                                                    col.append(ind_tot)
                                                    data.append(coeff*vals[i])
                                    elif vcen_type=='C' and var_type=='C':
                                        for c in vcen_varInds:
                                            if cells_shared_marked[c]==1:
                                                raise NotImplementedError("Not implemented")
                                                row_ind = self.GetCellIndexToIndexTotal(vcen, c, reg_name)
                                                ind_init, vals = self.CellC_C(c)
                                                for i in range(len(ind_init)):
                                                    ind_tot = self.GetCellIndexToIndexTotal(var, ind_init[i], reg_name)
                                                    row.append(row_ind)
                                                    col.append(ind_tot)
                                                    data.append(coeff*vals[i])
                                    elif vcen_type=='SX' and var_type=='N':
                                        for s in vcen_varInds:
                                            if sides_shared_marked[s]==1:
                                                onesided_, n_cell_cov = self._Side_GetOneSidedness(s, reg_cells_marked)
                                                eq_key_num = n_cell_cov_to_str_dic[eq_key]
                                                if n_cell_cov!=eq_key_num:
                                                    continue
                                                row_ind = self.GetSideIndexToIndexTotal(vcen, s, reg_name)
                                                ind_init, vals = self.SideX_N(s)
                                                for i in range(len(ind_init)):
                                                    ind_tot = self.GetNodeIndexToIndexTotal(var, ind_init[i], reg_name)
                                                    row.append(row_ind)
                                                    col.append(ind_tot)
                                                    data.append(coeff*vals[i])
                                    elif vcen_type=='SY' and var_type=='N':
                                        for s in vcen_varInds:
                                            if sides_shared_marked[s]==1:
                                                onesided_, n_cell_cov = self._Side_GetOneSidedness(s, reg_cells_marked)
                                                eq_key_num = n_cell_cov_to_str_dic[eq_key]
                                                if n_cell_cov!=eq_key_num:
                                                    continue
                                                row_ind = self.GetSideIndexToIndexTotal(vcen, s, reg_name)
                                                ind_init, vals = self.SideY_N(s)
                                                for i in range(len(ind_init)):
                                                    ind_tot = self.GetNodeIndexToIndexTotal(var, ind_init[i], reg_name)
                                                    row.append(row_ind)
                                                    col.append(ind_tot)
                                                    data.append(coeff*vals[i])
                                    elif vcen_type=='N' and var_type=='SX':
                                        for n in vcen_varInds:
                                            if nodes_shared_marked[n]==1:
                                                n_cell_cov = self.Node_GetNumRegionCellsCovered(n, reg_cells_marked)
                                                eq_key_num = n_cell_cov_to_str_dic[eq_key]
                                                if n_cell_cov!=eq_key_num:
                                                    continue
                                                if eq_key=='two_x' or eq_key=='two_y':
                                                    continue
                                                row_ind = self.GetNodeIndexToIndexTotal(vcen, n, reg_name)
                                                ind_init, vals = self.NodeN_SX(n)
                                                for i in range(len(ind_init)):
                                                    ind_tot = self.GetSideIndexToIndexTotal(var, ind_init[i], reg_name)
                                                    row.append(row_ind)
                                                    col.append(ind_tot)
                                                    data.append(coeff*vals[i])
                                    elif vcen_type=='N' and var_type=='SY':
                                        for n in vcen_varInds:
                                            if nodes_shared_marked[n]==1:
                                                n_cell_cov = self.Node_GetNumRegionCellsCovered(n, reg_cells_marked)
                                                eq_key_num = n_cell_cov_to_str_dic[eq_key]
                                                if n_cell_cov!=eq_key_num:
                                                    continue
                                                if eq_key=='two_x' or eq_key=='two_y':
                                                    continue
                                                row_ind = self.GetNodeIndexToIndexTotal(vcen, n, reg_name)
                                                ind_init, vals = self.NodeN_SY(n)
                                                for i in range(len(ind_init)):
                                                    ind_tot = self.GetSideIndexToIndexTotal(var, ind_init[i], reg_name)
                                                    row.append(row_ind)
                                                    col.append(ind_tot)
                                                    data.append(coeff*vals[i])
                                    else:
                                        print('vcen_type: ', vcen_type, ' var_type: ', var_type)
                                        raise NotImplementedError("Not implemented")
                                elif der_ord==1:
                                    if vcen_type=='C' and var_type=='SX':
                                        assert der_vars==[self.y]
                                        for c in vcen_varInds:
                                            if cells_shared_marked[c]==1:
                                                raise NotImplementedError("Not implemented")
                                    elif vcen_type=='C' and var_type=='SY':
                                        assert der_vars==[self.x]
                                        for c in vcen_varInds:
                                            if cells_shared_marked[c]==1:
                                                raise NotImplementedError("Not implemented")
                                    elif vcen_type=='N' and var_type=='SX':
                                        assert der_vars==[self.x]
                                        for n in vcen_varInds:
                                            if nodes_shared_marked[n]==1:
                                                n_cell_cov = self.Node_GetNumRegionCellsCovered(n, reg_cells_marked)
                                                eq_key_num = n_cell_cov_to_str_dic[eq_key]
                                                if n_cell_cov!=eq_key_num:
                                                    continue
                                                if eq_key=='two':
                                                    continue
                                                onesided_lr, onesided_du = self._Node_GetOneSidedness(n, reg_sides_marked, reg_cells_marked)
                                                if onesided_lr==False and eq_key=='two_y':  # x directed boundary
                                                    continue
                                                if onesided_du==False and eq_key=='two_x':  # y directed boundary
                                                    continue
                                                ind_init, vals = self.Node_dxSX(n, onesided=onesided_lr)
                                                row_ind = self.GetNodeIndexToIndexTotal(vcen, n, reg_name)
                                                for i in range(len(ind_init)):
                                                    ind_tot = self.GetSideIndexToIndexTotal(var, ind_init[i], reg_name)
                                                    row.append(row_ind)
                                                    col.append(ind_tot)
                                                    data.append(coeff*vals[i])
                                    elif vcen_type=='N' and var_type=='SY':
                                        assert der_vars==[self.y]
                                        for n in vcen_varInds:
                                            if nodes_shared_marked[n]==1:
                                                n_cell_cov = self.Node_GetNumRegionCellsCovered(n, reg_cells_marked)
                                                eq_key_num = n_cell_cov_to_str_dic[eq_key]
                                                if n_cell_cov!=eq_key_num:
                                                    continue
                                                if eq_key=='two':
                                                    continue
                                                onesided_lr, onesided_du = self._Node_GetOneSidedness(n, reg_sides_marked, reg_cells_marked)
                                                if onesided_lr==False and eq_key=='two_y':
                                                    continue
                                                if onesided_du==False and eq_key=='two_x':
                                                    continue
                                                ind_init, vals = self.Node_dySY(n, onesided=onesided_du)
                                                row_ind = self.GetNodeIndexToIndexTotal(vcen, n, reg_name)
                                                for i in range(len(ind_init)):
                                                    ind_tot = self.GetSideIndexToIndexTotal(var, ind_init[i], reg_name)
                                                    row.append(row_ind)
                                                    col.append(ind_tot)
                                                    data.append(coeff*vals[i])
                                    elif vcen_type=='SX' and var_type=='N':
                                        assert der_vars==[self.x]
                                        for s in vcen_varInds:
                                            if sides_shared_marked[s]==1:
                                                onesided_, n_cell_cov = self._Side_GetOneSidedness(s, reg_cells_marked)
                                                eq_key_num = n_cell_cov_to_str_dic[eq_key]
                                                if n_cell_cov!=eq_key_num:
                                                    continue
                                                row_ind = self.GetSideIndexToIndexTotal(vcen, s, reg_name)
                                                ind_init, vals = self.Side_dxyN(s)
                                                for i in range(len(ind_init)):
                                                    ind_tot = self.GetNodeIndexToIndexTotal(var, ind_init[i], reg_name)
                                                    row.append(row_ind)
                                                    col.append(ind_tot)
                                                    data.append(coeff*vals[i])
                                    elif vcen_type=='SY' and var_type=='N':
                                        assert der_vars==[self.y]
                                        for s in vcen_varInds:
                                            if sides_shared_marked[s]==1:
                                                onesided_, n_cell_cov = self._Side_GetOneSidedness(s, reg_cells_marked)
                                                eq_key_num = n_cell_cov_to_str_dic[eq_key]
                                                if n_cell_cov!=eq_key_num:
                                                    continue
                                                row_ind = self.GetSideIndexToIndexTotal(vcen, s, reg_name)
                                                ind_init, vals = self.Side_dxyN(s)
                                                for i in range(len(ind_init)):
                                                    ind_tot = self.GetNodeIndexToIndexTotal(var, ind_init[i], reg_name)
                                                    row.append(row_ind)
                                                    col.append(ind_tot)
                                                    data.append(coeff*vals[i])
                                    elif vcen_type=='SX' and var_type=='C':
                                        assert der_vars==[self.y]
                                        for s in vcen_varInds:
                                            if sides_shared_marked[s]==1:
                                                onesided_, n_cell_cov = self._Side_GetOneSidedness(s, reg_cells_marked)
                                                eq_key_num = n_cell_cov_to_str_dic[eq_key]
                                                if n_cell_cov!=eq_key_num:
                                                    continue
                                                ind_init, vals = self.Side_dxyC(s, onesided=onesided_)
                                                row_ind = self.GetSideIndexToIndexTotal(vcen, s, reg_name)
                                                for i in range(len(ind_init)):
                                                    ind_tot = self.GetCellIndexToIndexTotal(var, ind_init[i], reg_name)
                                                    row.append(row_ind)
                                                    col.append(ind_tot)
                                                    data.append(coeff*vals[i])
                                    elif vcen_type=='SY' and var_type=='C':
                                        assert der_vars==[self.x]
                                        for s in vcen_varInds:
                                            if sides_shared_marked[s]==1:
                                                onesided_, n_cell_cov = self._Side_GetOneSidedness(s, reg_cells_marked)
                                                eq_key_num = n_cell_cov_to_str_dic[eq_key]
                                                if n_cell_cov!=eq_key_num:
                                                    continue
                                                ind_init, vals = self.Side_dxyC(s, onesided=onesided_)
                                                row_ind = self.GetSideIndexToIndexTotal(vcen, s, reg_name)
                                                for i in range(len(ind_init)):
                                                    ind_tot = self.GetCellIndexToIndexTotal(var, ind_init[i], reg_name)
                                                    row.append(row_ind)
                                                    col.append(ind_tot)
                                                    data.append(coeff*vals[i])
                                    else:
                                        raise NotImplementedError("Not implemented")
                                else:
                                    raise NotImplementedError("Not implemented")
        row = np.array(row)
        col = np.array(col)
        data = np.array(data)
        A_coo_nx = coo_matrix((data, (row, col)), shape=(n_total,n_total), dtype=complex)
        if save==True:
            self.initialMatrix_shared_coo = A_coo_nx      #@
        return A_coo_nx
        
    '''
    def ApplyContinuityConditionBetweenRegions(self, A_coo, conds):
        """ conds = [[[reg1, var1], [reg2, var2], [a0, a1, a2]], ... ]
            a0*var1 + a1*var*2 + a2 = 0
        """
        row, col, data = list(A_coo.row), list(A_coo.col), list(A_coo.data)
        n_total = A_coo.shape[0]
        ind_row = len(row)
        sh_list_list = [self.CellsToIndexTotalShared, self.NodesToIndexTotalShared,\
                 self.SidesToIndexTotalShared]
        for sh_list in sh_list_list:        ## any of three shared lists
            for i in range(1, len(sh_list)):
                elems_sh = [e[0] for e in sh_list[i][1:len(sh_list[i])]]
                n_sh = len(elems_sh)
                assert n_sh>=2
                for j in range(n_sh-1):
                    row.append(ind_row)
                    col.append(elems_sh[j])     ## to be fixed (replace with total indices)
                    data.append(1.0)
                    row.append(ind_row)
                    col.append(elems_sh[j+1])
                    data.append(-1.0)
                    ind_row += 1
        assert ind_row==n_total
        row = np.array(row)
        col = np.array(col)
        data = np.array(data)
        A_coo = coo_matrix((data, (row, col)), shape=(n_total, n_total), dtype=complex)
        return A_coo
    '''
    
    def ApplyDirichletBoundaryCondition(self, A_coo, save=False):
        """Assumes zero boundary condition"""
        row, col, data = list(A_coo.row), list(A_coo.col), list(A_coo.data)
        n_total = A_coo.shape[0]
        n_b, s_b = self.GetBoundaryNodesAndSides()  # node boundary, side boundary
        n_b_vr = self.nodes_GetAssociatedVariablesAndRegions(n_b)   #n_b vars, regions
        s_b_vr = self.sides_GetAssociatedVariablesAndRegions(s_b)   #s_b vars, regions
        indtot_mkd = [0]*n_total
        for i in range(len(n_b)):
            n = n_b[i]
            for vr in n_b_vr[i]:
                ind_tot = self.GetNodeIndexToIndexTotal(vr[0], n, vr[1])
                indtot_mkd[ind_tot] = 1
        for i in range(len(s_b)):
            s = s_b[i]
            for vr in s_b_vr[i]:
                ind_tot = self.GetSideIndexToIndexTotal(vr[0], s, vr[1])
                indtot_mkd[ind_tot] = 1
        elems_to_del = [0]*len(row)
        n_elem_next = 0
        for i in range(len(row)):
            if indtot_mkd[row[i]]==1 or indtot_mkd[col[i]]==1:
                elems_to_del[i] = 1
            else:
                n_elem_next += 1
        self.ElemTotalToElemFinal = [-1]*n_total      #@    :new data member
        ind_next = 0
        for i in range(n_total):
            if indtot_mkd[i]==0:
                self.ElemTotalToElemFinal[i] = ind_next
                ind_next += 1
        n_total_nx = ind_next
        row_nx, col_nx, data_nx = [-1]*n_elem_next, [-1]*n_elem_next, [0j]*n_elem_next
        elem_nx = 0
        for i in range(len(row)):
            if elems_to_del[i]==0:
                row_nx[elem_nx] = self.ElemTotalToElemFinal[row[i]]
                col_nx[elem_nx] = self.ElemTotalToElemFinal[col[i]]
                data_nx[elem_nx] = data[i]
                elem_nx += 1
        assert elem_nx==n_elem_next
        n_total = n_total_nx
        row_nx = np.array(row_nx)
        col_nx = np.array(col_nx)
        data_nx = np.array(data_nx)
        A_coo_nx = coo_matrix((data_nx, (row_nx, col_nx)), shape=(n_total, n_total), dtype=complex)
        if save==True:
            self.finalMatrix_coo = A_coo_nx
        return A_coo_nx
        
    def EliminateVarible_Direct(self, A_coo, var_list):
        """ It eliminates the requested variables by direct substitution
        """
        row, col, data = A_coo.row, A_coo.col, A_coo.data
        n_total = A_coo.shape[0]
        indfinal_mkd = [0]*n_total
        var_type_dic, _ = self._GetVariablesAndAssociatedElements()
        for var in var_list:
            v_type = var_type_dic[var][0]
            for reg_name, eqs_vars_list in self.Regions.items():
                RegVars = self.Vars[reg_name]
                selfRegVar_types = self._GetVariableIndexAssociatedWithRegionsVarsAndType(reg_name)
                EQs = eqs_vars_list[0]
                n_EQs = len(EQs)
                for i_eq in range(n_EQs):
                    vcen_type = selfRegVar_types[i_eq][1]
                    vcen_ind = selfRegVar_types[i_eq][0]    #index inside self.Vars[reg_name]
                    vcen_varInds = RegVars[vcen_ind][1]     #elements (node,cell..) associated with the center variable
                    if (v_type=='SX' or v_type=='SY') and vcen_type==v_type:
                        for s in vcen_varInds:
                            ind_tot = self.GetSideIndexToIndexTotal(var, s, reg_name)
                            ind_final = self.ElemTotalToElemFinal[ind_tot]
                            if ind_final>=0:
                                indfinal_mkd[ind_final] = 1
                    elif v_type=='C' and vcen_type==v_type:
                        for c in vcen_varInds:
                            ind_tot = self.GetCellIndexToIndexTotal(var, c, reg_name)
                            ind_final = self.ElemTotalToElemFinal[ind_tot]
                            if ind_final>=0:
                                indfinal_mkd[ind_final] = 1
                    elif v_type=='N' and vcen_type==v_type:
                        for n in vcen_varInds:
                            ind_tot = self.GetNodeIndexToIndexTotal(var, n, reg_name)
                            ind_final = self.ElemTotalToElemFinal[ind_tot]
                            if ind_final>=0:
                                indfinal_mkd[ind_final] = 1
        n_ind_fin_mkd = 0
        for i in range(n_total):
            if indfinal_mkd[i]==1:
                n_ind_fin_mkd += 1
        self.ElemFinalToElemReduced = [-1]*n_total          #@ new data member
        n_ind_nonmkd = n_total - n_ind_fin_mkd
        ind_mkd = 0
        ind_nonmkd = 0
        for i in range(n_total):
            if indfinal_mkd[i]==0:
                self.ElemFinalToElemReduced[i] = ind_nonmkd
                ind_nonmkd += 1
            else:
                assert indfinal_mkd[i]==1
                self.ElemFinalToElemReduced[i] = n_ind_nonmkd + ind_mkd
                ind_mkd += 1
        for i in range(n_total):
            assert self.ElemFinalToElemReduced[i]>=0
            
        row_perm, col_perm, data_perm = [-1]*n_total, list(range(n_total)), [1.0]*n_total
        for i in range(n_total):
            row_perm[i] = self.ElemFinalToElemReduced[i]
        mat_perm = coo_matrix((data_perm, (row_perm, col_perm)), shape=(n_total, n_total), dtype=float).tocsr()
        A_p_csr = mat_perm*A_coo.tocsr()*mat_perm.T     #permutated matrix
        N0 = n_ind_nonmkd
        A_00 = A_p_csr[0:N0, 0:N0]
        A_01 = A_p_csr[0:N0, N0:n_total]
        A_10 = A_p_csr[N0:n_total, 0:N0]
        A_11 = A_p_csr[N0:n_total, N0:n_total]
        print('A_11.shape: ', A_11.shape, 'A_11.nnz: ', A_11.nnz)
        #A_red = A_00 - A_01*sp.sparse.linalg.inv(A_11.tocsc())*A_10         #reduced matrix
        #TODO: handle sp.sparse.linalg.inv manually for efficiency, sp.sparse.linalg.inv acts unreasonably slow
        A_red = A_00 - A_01*self.InvertMatrixComponents_csr(A_11)*A_10         #reduced matrix

        #A_red = A_00 - A_01*sp.sparse.linalg.spsolve(A_11.tocsc(), sp.sparse.eye(A_11.shape[0], \
        #    A_11.shape[1], format='csr'), permc_spec='COLAMD', use_umfpack=False)*A_10         #reduced matrix
        return A_red
        
    def InvertMatrixComponents_csr(self, A_csr):
        print("warning: this inversion mechanism only works for diagonal matrices")
        assert A_csr.nnz == A_csr.shape[0]
        data, indices, indptr = A_csr.data, A_csr.indices, A_csr.indptr
        data_inv = np.zeros(data.shape, dtype=data.dtype)
        for i in range(len(data)):
            data_inv[i] = 1.0/data[i]
        return sp.sparse.csr_matrix((data_inv, indices, indptr), dtype=data.dtype)
        
    def SolveEigenvalues(self, A_coo, n_eig=1, lambda_0=None, vec_0=None, tol=0, maxiter=None):
        A_csr = A_coo.tocsr()
        vals, vecs = sp.sparse.linalg.eigs(A_csr, k=n_eig, sigma=lambda_0, v0=vec_0, maxiter=maxiter, tol=tol)
        return [vals, vecs]

    def GetVecReducedToVecFinal(self, vec_red, A_coo_final):
        """ vec_red : column vector, reduced vector
            A_coo_final : A_coo final from which the reduced matrix is constructed
            (the one after the application of the boundary condition)
        """
        assert self.ElemFinalToElemReduced!=None
        n_total = A_coo_final.shape[0]
        row_perm, col_perm, data_perm = [-1]*n_total, list(range(n_total)), [1.0]*n_total
        for i in range(n_total):
            row_perm[i] = self.ElemFinalToElemReduced[i]
        mat_perm = coo_matrix((data_perm, (row_perm, col_perm)), shape=(n_total, n_total), dtype=float).tocsr()
        A_p_csr = mat_perm*A_coo_final.tocsr()*mat_perm.T     #permutated matrix
        N0 = vec_red.shape[0]
        A_00 = A_p_csr[0:N0, 0:N0]
        A_01 = A_p_csr[0:N0, N0:n_total]
        A_10 = A_p_csr[N0:n_total, 0:N0]
        A_11 = A_p_csr[N0:n_total, N0:n_total]
        #vec_extra = -sp.sparse.linalg.inv(A_11.tocsc())*A_10*vec_red
        #TODO: the following inversion only works for diagonal matrices
        vec_extra = -self.InvertMatrixComponents_csr(A_11)*A_10*vec_red
        
        #print('vec_red.shape:', vec_red.shape)
        #print('vec_extra.shape:', vec_extra.shape)
        vec_fin = np.concatenate((vec_red, vec_extra))
        vec_fin = mat_perm.T*vec_fin
        return vec_fin
        
        
    def GetFromSolutionIndexToInitialIndex(self, vec, A_coo_final=None):
        """ Takes the solution vectors and assigns respective values to each 
            element
        """ 
        #TODO: shared elements that don't have the same index for the sharing regions,
        # are not treated (it is assumed that shared elements have the same index)
        vec_fin = vec
        if self.ElemFinalToElemReduced!=None:
            assert A_coo_final!=None
            vec_fin = self.GetVecReducedToVecFinal(vec, A_coo_final)
        var_type_dic, n_vars_dic = self._GetVariablesAndAssociatedElements()
        #print('var_type_dic: ', var_type_dic)
        #print('n_vars_dic:', n_vars_dic)
        self.CellsIndexInitFinalValue = [[0j]*len(self.cells) for i in range(n_vars_dic['C'])]
        self.NodesIndexInitFinalValue = [[0j]*len(self.nodes) for i in range(n_vars_dic['N'])]
        self.SidesIndexInitFinalValue = [[0j]*len(self.sides) for i in range(max(n_vars_dic['SX'],n_vars_dic['SY']))]
        for var, v_type__Nv in var_type_dic.items():
            v_type, v_ind = var_type_dic[var]
            for reg_name, eqs_vars_list in self.Regions.items():
                RegVars = self.Vars[reg_name]
                selfRegVar_types = self._GetVariableIndexAssociatedWithRegionsVarsAndType(reg_name)
                EQs = eqs_vars_list[0]
                n_EQs = len(EQs)
                for i_eq in range(n_EQs):
                    vcen_type = selfRegVar_types[i_eq][1]
                    vcen_ind = selfRegVar_types[i_eq][0]    #index inside self.Vars[reg_name]
                    vcen_varInds = RegVars[vcen_ind][1]     #elements (node,cell..) associated with the center variable
                    if (v_type=='SX' or v_type=='SY') and vcen_type==v_type:
                        for s in vcen_varInds:
                            ind_tot = self.GetSideIndexToIndexTotal(var, s, reg_name)
                            ind_final = self.ElemTotalToElemFinal[ind_tot]
                            if ind_final>=0:
                                self.SidesIndexInitFinalValue[v_ind][s] = vec_fin[ind_final]
                    if v_type=='C' and vcen_type==v_type:
                        for c in vcen_varInds:
                            ind_tot = self.GetCellIndexToIndexTotal(var, c, reg_name)
                            ind_final = self.ElemTotalToElemFinal[ind_tot]
                            if ind_final>=0:
                                self.CellsIndexInitFinalValue[v_ind][c] = vec_fin[ind_final]
                    if v_type=='N' and vcen_type==v_type:
                        for n in vcen_varInds:
                            ind_tot = self.GetNodeIndexToIndexTotal(var, n, reg_name)
                            ind_final = self.ElemTotalToElemFinal[ind_tot]
                            if ind_final>=0:
                                self.NodesIndexInitFinalValue[v_ind][n] = vec_fin[ind_final]
        return
           
    def InterpolateVariablesValuesToCellCenters(self, var_list):
        """ It eliminates the requested variables by direct substitution
        """
        l, r, d, u = range(4)
        vals_interp = [[0j]*len(self.cells) for i in range(len(var_list))]
        var_type_dic, _ = self._GetVariablesAndAssociatedElements()
        for i_v, var in enumerate(var_list):
            v_type, v_ind = var_type_dic[var]
            if v_type=='SX':
                for i_c in range(len(self.cells)):
                    cell_arr = self.cells[i_c]
                    sd = cell_arr[d]
                    s_d_val = None
                    if sd>=0:
                        s_d_val = self.SidesIndexInitFinalValue[v_ind][sd]
                    else:
                        s_d_val = (self.SidesIndexInitFinalValue[v_ind][self.sidesDiv[abs(sd)][0]] +\
                            self.SidesIndexInitFinalValue[v_ind][self.sidesDiv[abs(sd)][1]])/2.0
                    su = cell_arr[u]
                    s_u_val = None
                    if su>=0:
                        s_u_val = self.SidesIndexInitFinalValue[v_ind][su]
                    else:
                        s_u_val = (self.SidesIndexInitFinalValue[v_ind][self.sidesDiv[abs(su)][0]] +\
                            self.SidesIndexInitFinalValue[v_ind][self.sidesDiv[abs(su)][1]])/2.0
                    vals_interp[i_v][i_c] = (s_d_val + s_u_val)/2.0
            elif v_type=='SY':
                for i_c in range(len(self.cells)):
                    cell_arr = self.cells[i_c]
                    sl = cell_arr[l]
                    s_l_val = None
                    if sl>=0:
                        s_l_val = self.SidesIndexInitFinalValue[v_ind][sl]
                    else:
                        s_l_val = (self.SidesIndexInitFinalValue[v_ind][self.sidesDiv[abs(sl)][0]] +\
                            self.SidesIndexInitFinalValue[v_ind][self.sidesDiv[abs(sl)][1]])/2.0
                    sr = cell_arr[r]
                    s_r_val = None
                    if sr>=0:
                        s_r_val = self.SidesIndexInitFinalValue[v_ind][sr]
                    else:
                        s_r_val = (self.SidesIndexInitFinalValue[v_ind][self.sidesDiv[abs(sr)][0]] +\
                            self.SidesIndexInitFinalValue[v_ind][self.sidesDiv[abs(sr)][1]])/2.0
                    vals_interp[i_v][i_c] = (s_l_val + s_r_val)/2.0
            elif v_type=='C':
                raise NotImplementedError()
            elif v_type=='N':
                raise NotImplementedError()
        return vals_interp
    
    def FindCellIndicesWithHighestValues(self, var_list, percent=10.0):
        """ It finds the elements with n percent higher values
        """
        l, r, d, u = range(4)
        vals_interp = [0.0]*len(self.cells)
        var_type_dic, _ = self._GetVariablesAndAssociatedElements()
        for i_v, var in enumerate(var_list):
            v_type, v_ind = var_type_dic[var]
            if v_type=='SX':
                for i_c in range(len(self.cells)):
                    cell_arr = self.cells[i_c]
                    sd = cell_arr[d]
                    s_d_val = None
                    if sd>=0:
                        s_d_val = self.SidesIndexInitFinalValue[v_ind][sd]
                    else:
                        s_d_val = (self.SidesIndexInitFinalValue[v_ind][self.sidesDiv[abs(sd)][0]] +\
                            self.SidesIndexInitFinalValue[v_ind][self.sidesDiv[abs(sd)][1]])/2.0
                    su = cell_arr[u]
                    s_u_val = None
                    if su>=0:
                        s_u_val = self.SidesIndexInitFinalValue[v_ind][su]
                    else:
                        s_u_val = (self.SidesIndexInitFinalValue[v_ind][self.sidesDiv[abs(su)][0]] +\
                            self.SidesIndexInitFinalValue[v_ind][self.sidesDiv[abs(su)][1]])/2.0
                    interp_val = abs(s_d_val + s_u_val)/2.0
                    if vals_interp[i_c]<interp_val:
                        vals_interp[i_c] = interp_val
            elif v_type=='SY':
                for i_c in range(len(self.cells)):
                    cell_arr = self.cells[i_c]
                    sl = cell_arr[l]
                    s_l_val = None
                    if sl>=0:
                        s_l_val = self.SidesIndexInitFinalValue[v_ind][sl]
                    else:
                        s_l_val = (self.SidesIndexInitFinalValue[v_ind][self.sidesDiv[abs(sl)][0]] +\
                            self.SidesIndexInitFinalValue[v_ind][self.sidesDiv[abs(sl)][1]])/2.0
                    sr = cell_arr[r]
                    s_r_val = None
                    if sr>=0:
                        s_r_val = self.SidesIndexInitFinalValue[v_ind][sr]
                    else:
                        s_r_val = (self.SidesIndexInitFinalValue[v_ind][self.sidesDiv[abs(sr)][0]] +\
                            self.SidesIndexInitFinalValue[v_ind][self.sidesDiv[abs(sr)][1]])/2.0
                    interp_val = abs(s_l_val + s_r_val)/2.0
                    if vals_interp[i_c]<interp_val:
                        vals_interp[i_c] = interp_val
            elif v_type=='C':
                raise NotImplementedError()
            elif v_type=='N':
                raise NotImplementedError()
        return self.getNpercentHigherstValues(np.array(vals_interp), percent)
        
    def getNpercentHigherstValues(self, vals, percent):
        ratio = percent/100.0
        n_ratio = int(ratio*len(vals))
        val_max = np.amax(vals)
        print('n_ratio: ', n_ratio, '  val_max: ', val_max)
        val_th = val_max/2
        n_high = np.sum(vals>val_th)
        if n_high<n_ratio:
            while True:
                n_high = np.sum(vals>val_th)
                if n_high>n_ratio:
                    break
                val_th = val_th/2.0
        """else:
            while True:
                n_high = np.sum(vals>val_th)
                if n_high<n_ratio:
                    break
                val_th += (val_max-val_th)/2.0"""
        inds_mkd = (vals>val_th)
        n_high = np.sum(inds_mkd)
        print('n_high: ', n_high)
        inds_high = [-1]*n_high
        i_next = 0
        for i in range(len(inds_mkd)):
            if inds_mkd[i]==True:
                inds_high[i_next] = i
                i_next += 1
        return inds_high
        
    def CalculateSurfaceIntegral(self, vec1, vec2):
        """ it estimates the integral(var1*var2)ds on the grid
            vec1 and vec2 are the values at the center of each cell (numpy vectors)
        """
        Nc = len(self.cells)
        assert len(vec1)==Nc
        assert len(vec2)==Nc
        ds = np.zeros(Nc)
        for c in range(len(self.cells)):
            dx, dy = self.GetCellSize(c)
            ds[c] = dx*dy
        I_s = (vec1*vec2).dot(ds)
        return I_s
                
    def cells_GetAssociatedVariablesAndRegions(self, cells):
        cells_marked = [-1]*len(self.cells)
        varReg = [[] for i in range(len(cells))]
        for i in range(len(cells)):
            cells_marked[cells[i]] = i
        for reg_name in self.Vars:
            var_ind_type_list = self.Vars[reg_name]
            for var_ind_type in var_ind_type_list: 
                if var_ind_type[2]=='C':
                    var = var_ind_type[0]
                    inds = var_ind_type[1]
                    for ind in inds:
                        if cells_marked[ind]>=0:
                            varReg[cells_marked[ind]].append([var, reg_name])
        if self.Vars!={}:
            for vr in varReg:
                assert len(vr)>0
        return varReg

    def nodes_GetAssociatedVariablesAndRegions(self, nodes):
        nodes_marked = [-1]*len(self.nodes)
        varReg = [[] for i in range(len(nodes))]
        for i in range(len(nodes)):
            nodes_marked[nodes[i]] = i
        for reg_name in self.Vars:
            var_ind_type_list = self.Vars[reg_name]
            for var_ind_type in var_ind_type_list: 
                if var_ind_type[2]=='N':
                    var = var_ind_type[0]
                    inds = var_ind_type[1]
                    for ind in inds:
                        if nodes_marked[ind]>=0:
                            varReg[nodes_marked[ind]].append([var, reg_name])
        if self.Vars!={}:
            for vr in varReg:
                assert len(vr)>0
        return varReg

    def sides_GetAssociatedVariablesAndRegions(self, sides):
        sides_marked = [-1]*len(self.sides)
        varReg = [[] for i in range(len(sides))]
        for i in range(len(sides)):
            sides_marked[sides[i]] = i
        for reg_name in self.Vars:
            var_ind_type_list = self.Vars[reg_name]
            for var_ind_type in var_ind_type_list: 
                if var_ind_type[2]=='SX' or var_ind_type[2]=='SY':
                    var = var_ind_type[0]
                    inds = var_ind_type[1]
                    for ind in inds:
                        if sides_marked[ind]>=0:
                            varReg[sides_marked[ind]].append([var, reg_name])
        if self.Vars!={}:
            for vr in varReg:
                assert len(vr)>0
        return varReg

    def GetSidesXDirected(self, sides):
        ## sides: list of side indices
        XYind = 4
        n_x = 0
        for s in sides:
            if self.sides[s][XYind] == 'X':
                n_x += 1
        sides_x = [-1]*n_x
        n_x = 0
        for s in sides:
            if self.sides[s][XYind] == 'X':
                sides_x[n_x] = s
                n_x += 1
        return sides_x

    def GetSidesYDirected(self, sides):
        ## sides: list of side indices
        XYind = 4
        n_y = 0
        for s in sides:
            if self.sides[s][XYind] == 'Y':
                n_y += 1
        sides_y = [-1]*n_y
        n_y = 0
        for s in sides:
            if self.sides[s][XYind] == 'Y':
                sides_y[n_y] = s
                n_y += 1
        return sides_y
        
    def FindNodesInsideRectangle(self, P00, P11):
        """ P00: lower left corner
            P11: upper right corner
        """
        nodesPoints = self.nodesPoints
        N = len(nodesPoints)
        x0, y0 = P00.x, P00.y
        x1, y1 = P11.x, P11.y
        nodes_inside = []
        for i in range(N):
            if x0<= nodesPoints[i].x <=x1 and y0<= nodesPoints[i].y <=y1:
                nodes_inside.append(i)
        return nodes_inside
        
    def FindCellsInsideRectangle(self, P00, P11):
        """ P00: lower left corner
            P11: upper right corner
        """
        nodesPoints = self.nodesPoints
        N = len(self.cells)
        x0, y0 = P00.x, P00.y
        x1, y1 = P11.x, P11.y
        cells_inside = []
        for i in range(N):
            n_corners = self.Cell_GetCornerNodes(i)
            cell_in = True
            for n in n_corners:
                if not (x0<= nodesPoints[n].x <=x1 and y0<= nodesPoints[n].y <=y1):
                    cell_in = False
                    break
            if cell_in:
                cells_inside.append(i)
        return cells_inside
        
    def FindCellsWithNodesInsideRectangle(self, P00, P11):
        """ if a cell has one or more nodes inside the rectangle it is selected
            P00: lower left corner
            P11: upper right corner
        """
        nodesPoints = self.nodesPoints
        N = len(self.cells)
        x0, y0 = P00.x, P00.y
        x1, y1 = P11.x, P11.y
        cells_inside = []
        for i in range(N):
            n_corners = self.Cell_GetCornerNodes(i)
            cell_in = False
            for n in n_corners:
                if (x0<= nodesPoints[n].x <=x1 and y0<= nodesPoints[n].y <=y1):
                    cell_in = True
                    break
            if cell_in:
                cells_inside.append(i)
        return cells_inside

    def Cells_GetNodesSide(self, cells, XY=True):
        # returns nodes and sides connected to the input cell
        l, r, d, u = range(4)
        s_nld, s_nru, s_cld, s_cru, XYind = range(5)   #sides::  s:side n:node c:cell
        sides_marked = [0]*len(self.sides)
        nodes_marked = [0]*len(self.nodes)
        for c in cells:
            cell_arr = self.cells[c]
            for i in range(4):
                s = cell_arr[i]
                if s>=0:
                    sides_marked[s] = 1
                    nodes_marked[self.sides[s][s_nld]] = 1
                    nodes_marked[self.sides[s][s_nru]] = 1
                else:
                    s0 = self.sidesDiv[abs(s)][0]
                    sides_marked[s0] = 1
                    nodes_marked[self.sides[s0][s_nld]] = 1
                    nodes_marked[self.sides[s0][s_nru]] = 1
                    s1 = self.sidesDiv[abs(s)][1]
                    sides_marked[s1] = 1
                    nodes_marked[self.sides[s1][s_nld]] = 1
                    nodes_marked[self.sides[s1][s_nru]] = 1
        n_node_mk = 0
        for i in range(len(nodes_marked)):
            if nodes_marked[i]==1:
                n_node_mk += 1
        nodes = [-1]*n_node_mk
        n_node_mk = 0
        for i in range(len(nodes_marked)):
            if nodes_marked[i]==1:
                nodes[n_node_mk] = i
                n_node_mk += 1
        n_side_mk = 0
        for i in range(len(sides_marked)):
            if sides_marked[i]==1:
                n_side_mk += 1
        sides = [-1]*n_side_mk
        n_side_mk = 0
        for i in range(len(sides_marked)):
            if sides_marked[i]==1:
                sides[n_side_mk] = i
                n_side_mk += 1
        if XY == True:
            return [nodes, self.GetSidesXDirected(sides), self.GetSidesYDirected(sides)]
        else:
            return [nodes, sides]


    def FindElementsMarkedByRectangle(self, P00, P11):
        cells = self.FindCellsInsideRectangle(P00, P11)
        nodes, sidesX, sidesY = self.Cells_GetNodesSide(cells, XY=True)
        return [cells, nodes, sidesX, sidesY]

    def FindCellsContactingRectangleBorder(self, P00, P11):
        """ P00: lower left corner
            P11: upper right corner
        """
        nodesPoints = self.nodesPoints
        N = len(self.cells)
        x0, y0 = P00.x, P00.y
        x1, y1 = P11.x, P11.y
        cells_border = []
        for i in range(N):
            nld, nlu, nrd, nru = self.Cell_GetCornerNodes(i)
            cell_bord = False
            if nodesPoints[nld].x<= x0 <=nodesPoints[nrd].x and (y0<= nodesPoints[nld].y <=y1 or y0<= nodesPoints[nlu].y <=y1):
                cell_bord = True
            elif nodesPoints[nld].x<= x1 <=nodesPoints[nrd].x and (y0<= nodesPoints[nld].y <=y1 or y0<= nodesPoints[nlu].y <=y1):
                cell_bord = True
            elif nodesPoints[nld].y<= y0 <=nodesPoints[nlu].y and (x0<= nodesPoints[nld].x <=x1 or x0<= nodesPoints[nrd].x <=x1):
                cell_bord = True
            elif nodesPoints[nld].y<= y1 <=nodesPoints[nlu].y and (x0<= nodesPoints[nld].x <=x1 or x0<= nodesPoints[nrd].x <=x1):
                cell_bord = True
            if cell_bord:
                cells_border.append(i)
        return cells_border
        
        
    def FindCellsSpanningOverThinRectangle(self, P00, P11):
        """ P00: lower left corner
            P11: upper right corner
        """
        nodesPoints = self.nodesPoints
        N = len(self.cells)
        x0, y0 = P00.x, P00.y
        x1, y1 = P11.x, P11.y
        assert x0<x1 and y0<y1
        cells_marked = [0]*N
        for i in range(N):
            nld, nlu, nrd, nru = self.Cell_GetCornerNodes(i)
            if nodesPoints[nld].x <= x0 and nodesPoints[nrd].x >= x1 \
                    and (y0<= nodesPoints[nld].y <=y1 or y0<= nodesPoints[nlu].y <=y1):
                cells_marked[i] = 1
            elif nodesPoints[nld].y <= y0 and nodesPoints[nru].y >= y1 \
                    and (x0<= nodesPoints[nld].x <=x1 or x0<= nodesPoints[nrd].x <=x1):
                cells_marked[i] = 1
        
        ## marking nodes inside rectangle whose connected cells are all outside the rectangle
        del_x_2 = (x1 - x0)/2.01
        del_y_2 = (y1 - y0)/2.01
        print('cecking nodes inside')
        nodes_inside = self.FindNodesInsideRectangle(P00, P11)
        for n in nodes_inside:
            cells_nb = self.Node_GetNeighboringCells(n)
            mark_all = True     #mark all neigboring cells
            for c in cells_nb:
                if c >=0:
                    if cells_marked[c]==1:
                        mark_all = False 
                        continue
            if mark_all:
                for c in cells_nb:
                    if c >=0:
                        level = self.cells[c][4]
                        if self.GetCellSizeX(level)>=del_x_2 or self.GetCellSizeY(level)>=del_y_2:
                            #print(self.GetCellSizeX(level), self.GetCellSizeY(level))
                            cells_marked[c] = 1

        n_marked = 0
        for cm in cells_marked:
            if cm==1:
                n_marked += 1
        
        cells_span = [-1]*n_marked
        n_marked = 0
        for i in range(N):
            if cells_marked[i]==1:
                cells_span[n_marked] = i
                n_marked += 1
        return cells_span
        

    def RefineRectangleBorders(self, rect):
        P00, P11 = rect
        cells_border = self.FindCellsContactingRectangleBorder(P00, P11)
        self.RefineCells(cells_border)
        return

    def RefineRectangle(self, rect):
        P00, P11 = rect
        cells_inside = self.FindCellsWithNodesInsideRectangle(P00, P11)
        self.RefineCells(cells_inside)
        return
        
    def RefineThinRectangle(self, rect):
        ##TODO: the meshed object might become several disconnected objects
        P00, P11 = rect
        while True:
            cells_in = self.FindCellsSpanningOverThinRectangle(P00, P11)
            if len(cells_in)>0:
                self.RefineCells(cells_in)
            else:
                break
        return
        
    def RegionsRemoveSharedCells(self, cells_list):
        """ cells_list : [cells0, cells1, ...]
            it removes the shared cells between cells0, cells1, ...
            cells0<cells1<cells2... in terms of priority for receiving the 
            shared cell (cells0: lowest priority)
        """
        N = len(self.cells)
        cells_all = [-1]*N
        for i in range(len(cells_list)):
            cells_i = cells_list[i]
            for j in range(len(cells_i)):
                cells_all[cells_i[j]] = i
        n_cell_i = [0]*len(cells_list)
        for i in range(N):
            if cells_all[i]>=0:
                n_cell_i[cells_all[i]] += 1
        cells_list_new = [None]*len(cells_list)
        for i in range(len(cells_list)):
            cells_list_new[i] = [-1]*n_cell_i[i]
        n_cell_i = [0]*len(cells_list)
        for i in range(N):
            li = cells_all[i]
            if li>=0:
                cells_list_new[li][n_cell_i[li]] = i
                n_cell_i[li] += 1
        return cells_list_new
        
        
    def GetDebugInfo(self, elem_str, elem_ind, info_type='connections'):
        """ elem_str : 'cell', 'side', 'node'
            info_type : 'connections', 'connections_connections'
        """
        def cell_getAllAssociatedIndices(elem_ind, connections=True):
            l, r, d, u = range(4)
            cell_arr = self.cells[elem_ind]
            cell_size = self.GetCellSize(elem_ind)
            inds = []
            varReg = self.cells_GetAssociatedVariablesAndRegions([elem_ind])
            for var_reg_name in varReg[0]:
                if len(var_reg_name)>0:
                    var, reg_name = var_reg_name
                    ind_tot = self.GetCellIndexToIndexTotal(var, elem_ind, reg_name)
                    ind_final = self.ElemTotalToElemFinal[ind_tot]
                    inds.append(['cell', elem_ind, reg_name, var, ind_tot, ind_final, cell_size])
            if connections:
                if cell_arr[l]>=0:
                    inds = inds + side_getAllAssociatedIndices(cell_arr[l], connections=False)
                else:
                    inds = inds + side_getAllAssociatedIndices(self.sidesDiv[abs(cell_arr[l])][0], connections=False)
                    inds = inds + side_getAllAssociatedIndices(self.sidesDiv[abs(cell_arr[l])][1], connections=False)
                if cell_arr[r]>=0:
                    inds = inds + side_getAllAssociatedIndices(cell_arr[r], connections=False)
                else:
                    inds = inds + side_getAllAssociatedIndices(self.sidesDiv[abs(cell_arr[r])][0], connections=False)
                    inds = inds + side_getAllAssociatedIndices(self.sidesDiv[abs(cell_arr[r])][1], connections=False)
                if cell_arr[d]>=0:
                    inds = inds + side_getAllAssociatedIndices(cell_arr[d], connections=False)
                else:
                    inds = inds + side_getAllAssociatedIndices(self.sidesDiv[abs(cell_arr[d])][0], connections=False)
                    inds = inds + side_getAllAssociatedIndices(self.sidesDiv[abs(cell_arr[d])][1], connections=False)
                if cell_arr[u]>=0:
                    inds = inds + side_getAllAssociatedIndices(cell_arr[u], connections=False)
                else:
                    inds = inds + side_getAllAssociatedIndices(self.sidesDiv[abs(cell_arr[u])][0], connections=False)
                    inds = inds + side_getAllAssociatedIndices(self.sidesDiv[abs(cell_arr[u])][1], connections=False)
            return inds
        def side_getAllAssociatedIndices(elem_ind, connections=True):
            s_nld, s_nru, s_cld, s_cru, XYind = range(5)   #sides::  s:side n:node c:cell
            side_arr = self.sides[elem_ind]
            side_size = self.GetSideSize(elem_ind)
            inds = []
            varReg = self.sides_GetAssociatedVariablesAndRegions([elem_ind])
            for var_reg_name in varReg[0]:
                if len(var_reg_name)>0:
                    var, reg_name = var_reg_name
                    ind_tot = self.GetSideIndexToIndexTotal(var, elem_ind, reg_name)
                    ind_final = self.ElemTotalToElemFinal[ind_tot]
                    inds.append(['side', elem_ind, reg_name, var, ind_tot, ind_final, side_size])
            if connections:
                inds = inds + node_getAllAssociatedIndices(side_arr[s_nld], connections=False)
                inds = inds + node_getAllAssociatedIndices(side_arr[s_nru], connections=False)
                if side_arr[s_cld]>=0:
                    inds = inds + cell_getAllAssociatedIndices(side_arr[s_cld], connections=False)
                if side_arr[s_cru]>=0:
                    inds = inds + cell_getAllAssociatedIndices(side_arr[s_cru], connections=False)
            return inds
        
        def node_getAllAssociatedIndices(elem_ind, connections=True):
            l, r, d, u = range(4)
            node_arr = self.nodes[elem_ind]
            inds = []
            varReg = self.nodes_GetAssociatedVariablesAndRegions([elem_ind])
            for var_reg_name in varReg[0]:
                if len(var_reg_name)>0:
                    var, reg_name = var_reg_name
                    ind_tot = self.GetNodeIndexToIndexTotal(var, elem_ind, reg_name)
                    ind_final = self.ElemTotalToElemFinal[ind_tot]
                    inds.append(['node', elem_ind, reg_name, var, ind_tot, ind_final])
            if connections:
                if node_arr[l]>=0:
                    inds = inds + side_getAllAssociatedIndices(node_arr[l], connections=False)
                if node_arr[r]>=0:
                    inds = inds + side_getAllAssociatedIndices(node_arr[r], connections=False)
                if node_arr[d]>=0:
                    inds = inds + side_getAllAssociatedIndices(node_arr[d], connections=False)
                if node_arr[u]>=0:
                    inds = inds + side_getAllAssociatedIndices(node_arr[u], connections=False)
            return inds

        if info_type=='connections':
            str_list = []
            if elem_str=='cell':
                l, r, d, u = range(4)
                cell_arr = self.cells[elem_ind]
                str_list.append("cell_ind: {ind}".format(ind=elem_ind))
                str_list.append("sides left: {left}  right: {right}  down: {down}  up: {up}"\
                .format(left=cell_arr[l], right=cell_arr[r], down=cell_arr[d], up=cell_arr[u]))
                str_list.append("associated indices: \n{}".format(cell_getAllAssociatedIndices(elem_ind)))
                varReg = self.cells_GetAssociatedVariablesAndRegions([elem_ind])
                for var_reg_name in varReg[0]:
                    if len(var_reg_name)>0:
                        var, reg_name = var_reg_name
                        ind_tot = self.GetCellIndexToIndexTotal(var, elem_ind, reg_name)
                        str_list.append("region: {}, var: {}".format(reg_name, var))
                        str_list.append("cell ind: {} total: {}".format(elem_ind, ind_tot))
                        if self.initialMatrix_coo!=None:
                            initRow = self.initialMatrix_coo.getrow(ind_tot)
                            str_list.append("Matrix_Init Row: {} \n{}".format(ind_tot, initRow))
                        if self.initialMatrix_shared_coo!=None:
                            initSharedRow = self.initialMatrix_shared_coo.getrow(ind_tot)
                            str_list.append("Matrix_Init_shared Row: {} \n{}".format(ind_tot, initSharedRow))
                        if self.finalMatrix_coo!=None:
                            finalRow = self.finalMatrix_coo.getrow(self.ElemTotalToElemFinal[ind_tot])
                            str_list.append("Matrix_final Row: {} \n{}".format(self.ElemTotalToElemFinal[ind_tot], finalRow))
            elif elem_str=='side':
                s_nld, s_nru, s_cld, s_cru, XYind = range(5)   #sides::  s:side n:node c:cell
                side_arr = self.sides[elem_ind]
                str_list.append("side_ind: {ind}".format(ind=elem_ind))
                if side_arr[XYind]=='X':
                    str_list.append("Direction : X")
                    str_list.append("nodes left: {n0}  right: {n1}  cells down: {c0}  up: {c1}"\
                    .format(n0=side_arr[s_nld], n1=side_arr[s_nru], c0=side_arr[s_cld], c1=side_arr[s_cru]))
                else:
                    assert side_arr[XYind]=='Y'
                    str_list.append("Direction : Y")
                    str_list.append("nodes down: {n0}  up: {n1}  cells left: {c0}  right: {c1}"\
                    .format(n0=side_arr[s_nld], n1=side_arr[s_nru], c0=side_arr[s_cld], c1=side_arr[s_cru]))
                str_list.append("associated indices: \n{}".format(side_getAllAssociatedIndices(elem_ind)))
                varReg = self.sides_GetAssociatedVariablesAndRegions([elem_ind])
                for var_reg_name in varReg[0]:
                    if len(var_reg_name)>0:
                        var, reg_name = var_reg_name
                        ind_tot = self.GetSideIndexToIndexTotal(var, elem_ind, reg_name)
                        str_list.append("region: {}, var: {}".format(reg_name, var))
                        str_list.append("side ind: {} total: {}".format(elem_ind, ind_tot))
                        if self.initialMatrix_coo!=None:
                            initRow = self.initialMatrix_coo.getrow(ind_tot)
                            str_list.append("Matrix_Init Row: {} \n{}".format(ind_tot, initRow))
                        if self.initialMatrix_shared_coo!=None:
                            initSharedRow = self.initialMatrix_shared_coo.getrow(ind_tot)
                            str_list.append("Matrix_Init_shared Row: {} \n{}".format(ind_tot, initSharedRow))
                        if self.finalMatrix_coo!=None:
                            finalRow = self.finalMatrix_coo.getrow(self.ElemTotalToElemFinal[ind_tot])
                            str_list.append("Matrix_final Row: {} \n{}".format(self.ElemTotalToElemFinal[ind_tot], finalRow))
            elif elem_str=='node':
                l, r, d, u = range(4)
                node_arr = self.nodes[elem_ind]
                node_pt = self.nodesPoints[elem_ind]
                str_list.append("node_ind: {ind}  position({pt_x}, {pt_y})".format(ind=elem_ind, pt_x=node_pt.x, pt_y=node_pt.y))
                str_list.append("sides left: {left}  right: {right}  down: {down}  up: {up}"\
                .format(left=node_arr[l], right=node_arr[r], down=node_arr[d], up=node_arr[u]))
                str_list.append("associated indices: \n{}".format(node_getAllAssociatedIndices(elem_ind)))
                varReg = self.nodes_GetAssociatedVariablesAndRegions([elem_ind])
                for var_reg_name in varReg[0]:
                    if len(var_reg_name)>0:
                        var, reg_name = var_reg_name
                        ind_tot = self.GetNodeIndexToIndexTotal(var, elem_ind, reg_name)
                        str_list.append("region: {}, var: {}".format(reg_name, var))
                        str_list.append("node ind: {} total: {}".format(elem_ind, ind_tot))
                        if self.initialMatrix_coo!=None:
                            initRow = self.initialMatrix_coo.getrow(ind_tot)
                            str_list.append("Matrix_Init Row: {} \n{}".format(ind_tot, initRow))
                        if self.initialMatrix_shared_coo!=None:
                            initSharedRow = self.initialMatrix_shared_coo.getrow(ind_tot)
                            str_list.append("Matrix_Init_shared Row: {} \n{}".format(ind_tot, initSharedRow))
                        if self.finalMatrix_coo!=None:
                            finalRow = self.finalMatrix_coo.getrow(self.ElemTotalToElemFinal[ind_tot])
                            str_list.append("Matrix_final Row: {} \n{}".format(self.ElemTotalToElemFinal[ind_tot], finalRow))
            str_out = ""
            for i in range(len(str_list)):
                str_out += str_list[i] + '\n'
            return str_out
        elif info_type=='connections_connections':
            str_list = []
            if elem_str=='cell':
                str_list.append(self.GetDebugInfo('cell', elem_ind, 'connections'))
                l, r, d, u = range(4)
                cell_arr = self.cells[elem_ind]
                str_list.append('cell ===> left side:')                
                if cell_arr[l]>=0:
                    str_list.append(self.GetDebugInfo('side', cell_arr[l], 'connections'))
                else:
                    str_list.append(self.GetDebugInfo('side', self.sidesDiv[abs(cell_arr[l])][0], 'connections'))
                    str_list.append(self.GetDebugInfo('side', self.sidesDiv[abs(cell_arr[l])][1], 'connections'))
                str_list.append('cell ===> right side:')                
                if cell_arr[r]>=0:
                    str_list.append(self.GetDebugInfo('side', cell_arr[r], 'connections'))
                else:
                    str_list.append(self.GetDebugInfo('side', self.sidesDiv[abs(cell_arr[r])][0], 'connections'))
                    str_list.append(self.GetDebugInfo('side', self.sidesDiv[abs(cell_arr[r])][1], 'connections'))
                str_list.append('cell ===> down side:')                
                if cell_arr[d]>=0:
                    str_list.append(self.GetDebugInfo('side', cell_arr[d], 'connections'))
                else:
                    str_list.append(self.GetDebugInfo('side', self.sidesDiv[abs(cell_arr[d])][0], 'connections'))
                    str_list.append(self.GetDebugInfo('side', self.sidesDiv[abs(cell_arr[d])][1], 'connections'))
                str_list.append('cell ===> up side:')                
                if cell_arr[u]>=0:
                    str_list.append(self.GetDebugInfo('side', cell_arr[u], 'connections'))
                else:
                    str_list.append(self.GetDebugInfo('side', self.sidesDiv[abs(cell_arr[u])][0], 'connections'))
                    str_list.append(self.GetDebugInfo('side', self.sidesDiv[abs(cell_arr[u])][1], 'connections'))
            elif elem_str=='side':
                str_list.append(self.GetDebugInfo('side', elem_ind, 'connections'))
                s_nld, s_nru, s_cld, s_cru, XYind = range(5)   #sides::  s:side n:node c:cell
                side_arr = self.sides[elem_ind]
                str_list.append('side ===> left/down node:')                
                str_list.append(self.GetDebugInfo('node', side_arr[s_nld], 'connections'))
                str_list.append('side ===> right/up node:')                
                str_list.append(self.GetDebugInfo('node', side_arr[s_nru], 'connections'))
                if side_arr[s_cld]>=0:
                    str_list.append('side ===> left/down cell:')                
                    str_list.append(self.GetDebugInfo('cell', side_arr[s_cld], 'connections'))
                if side_arr[s_cru]>=0:
                    str_list.append('side ===> right/up cell:')                
                    str_list.append(self.GetDebugInfo('cell', side_arr[s_cru], 'connections'))
            elif elem_str=='node':
                str_list.append(self.GetDebugInfo('node', elem_ind, 'connections'))
                l, r, d, u = range(4)
                node_arr = self.nodes[elem_ind]
                if node_arr[l]>=0:
                    str_list.append('node ===> left side:')                
                    str_list.append(self.GetDebugInfo('side', node_arr[l], 'connections'))
                if node_arr[r]>=0:
                    str_list.append('node ===> right side:')                
                    str_list.append(self.GetDebugInfo('side', node_arr[r], 'connections'))
                if node_arr[d]>=0:
                    str_list.append('node ===> down side:')                
                    str_list.append(self.GetDebugInfo('side', node_arr[d], 'connections'))
                if node_arr[u]>=0:
                    str_list.append('node ===> up side:')                
                    str_list.append(self.GetDebugInfo('side', node_arr[u], 'connections'))
            str_out = ""
            for i in range(len(str_list)):
                str_out += str_list[i] + '\n'
            return str_out
        
       
### --------------------------------------------------------------------------
### --------------------------------------------------------------------------
from sympy import Mul, Add, Derivative, Integer

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

        

### --------------------------------------------------------------------------
### --------------------------------------------------------------------------
        
       
##------------------------   GUI -------------------

#from multiprocessing import Process, Queue
from threading import Thread
import queue
from tkinter import *
import time
from random import randint

class GUIMaker2D(Thread):
    
    def __init__(self, wx=1.0, wy=1.0, wxp=200, wyp=200, que=None):
        self.wxp = wxp  #p:pixel
        self.wyp = wyp
        self.wx = wx    #clipping window width
        self.wy = wy
        self.wx0 = 0.0  #clipping window lower left corner
        self.wy0 = 0.0
        self.x0p_offset = 10
        self.y0p_offset = 10
        self.queue = que
        self.recGrid2d = None
        self.showGrid = False
        self.testMode = False
        self.quit = False
        self.elements = []
        Thread.__init__(self)
        #self.daemon = True
        self.start()
        return
        
    def SetCommChannel(self, que):
        self.queue = que
        return
        
    def RealWorldToCanvasPixel(self, x, y):
        c_x = self.x0p_offset + int((x - self.wx0)/(self.wx)*(self.wxp - 2*self.x0p_offset))
        c_y = self.y0p_offset + int((y - self.wy0)/(self.wy)*(self.wyp - 2*self.y0p_offset))
        c_y = self.wyp - c_y
        return [c_x, c_y]
        
    def CanvasPixelToRealWorld(self, xp, yp):
        yp = self.wyp - yp
        x = self.wx0 + (xp - self.x0p_offset)/(self.wxp - 2*self.x0p_offset)*self.wx
        y = self.wy0 + (yp - self.y0p_offset)/(self.wyp - 2*self.y0p_offset)*self.wy
        return [x, y]

    def run(self):
        self.root = Tk()
        self.canvas = Canvas(self.root, bg='white', width=self.wxp, height=self.wyp)
        self.canvas.pack(fill='both', expand=YES)
        self.DrawAxis()

        """
        self.text_box = Text(self.root, height=5, width=int(self.wxp/10))
        self.text_box.pack(fill=BOTH)
        """
        self.scrol_text = Scrollbar(self.root, width=16)
        self.text_box = Text(self.root, height=10, width=int(self.wxp/8))
        self.scrol_text.pack(side=RIGHT, fill=Y)
        self.text_box.pack(side=LEFT, fill=X, expand=True)
        self.scrol_text.config(command=self.text_box.yview)
        self.text_box.config(yscrollcommand=self.scrol_text.set)
        
        self.root.event_add("<<DebugInfoEv>>", "<Shift-Control-Button-1>")
        self.root.event_add("<<PanToEv>>", "<Shift-Button-1>")
        self.root.event_add("<<ZoomInEv>>", "<Control-Button-1>")
        self.root.event_add("<<ZoomOutEv>>", "<Control-Button-3>")
        #self.canvas.bind("<<ProcCommEv>>", self.DoCommand)
        self.canvas.bind("<Configure>", self.adjustsize)
        self.canvas.bind("<<DebugInfoEv>>", self.clickSide)
        self.canvas.bind("<<PanToEv>>", self.panOnClick)
        self.canvas.bind("<<ZoomInEv>>", self.zoomInClick)
        self.canvas.bind("<<ZoomOutEv>>", self.zoomOutClick)

        self.root.protocol("WM_DELETE_WINDOW", self.TerminateRoot)
        
        self.root.after(1000, self.CommandLoop)
        self.root.mainloop()
        print('thread exiting..!')
        return
        
    def TerminateRoot(self):
        print('stopping thread!')
        self.quit = True
        self.root.quit()
        self.root.destroy()
        del self.text_box
        del self.scrol_text
        del self.canvas
        del self.root
        return

    def CommandLoop(self):
        if not self.quit:
            if self.queue != None:
                if not self.queue.empty():                
                    comm = self.queue.get()
                    self.ProcessCommand(comm)
                    self.queue.task_done()
            self.canvas.update_idletasks()
            self.root.after(1000, self.CommandLoop)
        return

    def DoCommand(self, event):
        if self.queue != None:
            if not self.queue.empty():                
                comm = self.queue.get()
                self.ProcessCommand(comm)
                self.queue.task_done()
        return

    def adjustsize(self, event):
        self.canvas.delete("all")
        self.wxp, self.wyp = event.width, event.height
        self.ClearCanvas()
        if self.showGrid:
            self.UpdateGrid()
        self.text_box.config(width=int(self.wxp/8))
        self.UpdateElementsOnCanvas()
        return
        
    def panOnClick(self, event):
        ##TODO: update canvas elements instead of deleting and adding them
        xp, yp = event.x, event.y
        x, y = self.CanvasPixelToRealWorld(xp, yp)
        self.SetWindowCenter(x, y)
        self.ClearCanvas()
        self.UpdateGrid()

    def zoomInClick(self, event):
        xp, yp = event.x, event.y
        x, y = self.CanvasPixelToRealWorld(xp, yp)
        self.ZoomInOnPoint(x, y)
        self.ClearCanvas()
        self.UpdateGrid()
                        
    def zoomOutClick(self, event):
        xp, yp = event.x, event.y
        x, y = self.CanvasPixelToRealWorld(xp, yp)
        self.ZoomOutOnPoint(x, y)
        self.ClearCanvas()
        self.UpdateGrid()

    def ProcessCommand(self, comm):
        if comm[0]=="updategrid":
            self.UpdateGrid()
        elif comm[0]=="clearcanvas":
            self.ClearCanvas()
        elif comm[0]=="attachgrid":
            self.AttachGrid(comm[1])
        elif comm[0]=="drawline":
            x0, y0, x1, y1, fill = comm[1]
            self.DrawLine(x0, y0, x1, y1, fill)
        elif comm[0]=="terminate":
            self.root.quit()
            self.root.destroy()
        elif comm[0]=="toggleTestMode":
            if self.testMode == False:
                self.testMode = True
            else:
                self.testMode = False
        elif comm[0]=="clearlog":
            self.clearText()
        elif comm[0]=="logtext":
            text = comm[1]
            self.logText(text)
        elif comm[0]=="addelement":
            elname, eltype, elparam = comm[1], comm[2], comm[3]
            if isinstance(elparam, list):
                elparam = tuple(elparam)
            self.AddElement(elname, eltype, elparam)
        elif comm[0]=="toggleelemvisibility":
            elname = comm[1]
            self.ToggleElementVisibility(elname)
        elif comm[0]=="MarkGridSides":
            sides, color = comm[1], comm[2]
            self.MarkGridSides(sides, color)
        elif comm[0]=="ColorCells":
            clrs = comm[1]
            self.ColorCells(clrs)
        return
        
    def AttachGrid(self, rg2d):
        self.recGrid2d = rg2d
        self.wx = rg2d.width
        self.wy = rg2d.height
        self.wx0 = 0.0
        self.wy0 = 0.0
        return
        
    def UpdateGrid(self):
        self.showGrid = True
        sides = self.recGrid2d.sides
        nodes = self.recGrid2d.nodes
        nodesPoints = self.recGrid2d.nodesPoints
        s_nld, s_nru, s_cld, s_cru = range(4)   #sides::  s:side n:node c:cell
        fill = 'black'
        for i in range(len(sides)):
            s = sides[i]
            p_ld = nodesPoints[s[s_nld]]
            p_ru = nodesPoints[s[s_nru]]
            x0, y0, x1, y1 = p_ld.x, p_ld.y, p_ru.x, p_ru.y
            self.DrawLine(x0, y0, x1, y1, fill=fill, tags='side_'+str(i))
            
    def MarkGridSides(self, sides, color):
        for s in sides:
            items = self.canvas.find_withtag('side_'+str(s))
            for it in items:
                self.canvas.itemconfig(it, fill=color)
        return

    def DrawLine(self, x0, y0, x1, y1, fill='red', tags="notag"):
        cx_0, cy_0 = self.RealWorldToCanvasPixel(x0, y0)
        cx_1, cy_1 = self.RealWorldToCanvasPixel(x1, y1)
        coord = cx_0, cy_0, cx_1, cy_1
        self.canvas.create_line(coord, fill=fill, tags=tags)
        #self.canvas.update_idletasks()
        return        
        
    def DrawRectangle(self, x0, y0, x1, y1, fill='red', outline='red', tags="notag"):
        cx_0, cy_0 = self.RealWorldToCanvasPixel(x0, y0)
        cx_1, cy_1 = self.RealWorldToCanvasPixel(x1, y1)
        coord = cx_0, cy_0, cx_1, cy_1
        self.canvas.create_rectangle(coord, fill=fill, outline=outline, tags=tags)
        #self.canvas.update_idletasks()
        return        

    def DrawAxis(self):
        fill = 'red'
        self.DrawLine(0.0, 0.0, self.wx, 0.0, fill=fill, tags="axis_x")
        self.DrawLine(0.0, 0.0, 0.0, self.wy, fill=fill, tags="axis_y")
        return

    def ClearCanvas(self):
        self.canvas.delete(ALL)
        self.DrawAxis()
        return
        
    def EnterTestMode(self):
        self.testMode = True
        return
        
    def clickSide(self, event):
        if self.testMode:
            if self.canvas.find_withtag(CURRENT):
                tags = self.canvas.gettags(CURRENT)[0]
                assert tags != "current"
                #self.logText(tags)
                if tags[0:5] == 'side_':
                    s = int(tags[5:])
                    self.logText(self.recGrid2d.GetDebugInfo('side', s, info_type='connections_connections'))
                    self.canvas.itemconfig(CURRENT, fill="yellow")
                    self.canvas.update_idletasks()
                    #self.canvas.after(500)
                    #self.canvas.itemconfig(CURRENT, fill="black")

    def logText(self, s):
        self.text_box.insert(END, s + '\n')

    def clearText(self):
        self.text_box.delete(1.0, END)

    def AddElement(self, elname, eltype, elparam, visible=True):
        """ elname: str
            eltype: 'rectangle', 'line'
            elparam: tuple - depends on eltype
            visible: True/False
        """
        elem = {'name':elname, 'type':eltype, 'param':elparam, 'visible':False}
        self.elements.append(elem)
        if visible:
            self.DrawElement(elname)
        
    def DeleteElement(self, name):
        for i in range(len(self.elements)):
            el = self.elements[i]
            if el['name']==name:
                if el['visible']==True:
                    self.RemoveFromCanvas('elem_'+name)
                del self.elements[i]
                break
        
    def DrawElement(self, name):
        for el in self.elements:
            if el['name']==name and el['visible']==False:
                if el['type']=='line':
                    self.DrawLine(*el['param'], tags='elem_{}'.format(name))
                    el['visible'] = True
                    break
                if el['type']=='rectangle':
                    self.DrawRectangle(*el['param'], tags='elem_{}'.format(name))
                    el['visible'] = True
                    break
    
    def RemoveFromCanvas(self, tag):
        items = self.canvas.find_withtag(tag)
        for it in items:
            self.canvas.delete(it)

    def ToggleElementVisibility(self, name):
        for el in self.elements:
            if el['name'] == name:
                if el['visible']==False:
                    self.DrawElement(name)
                else:
                    self.RemoveFromCanvas('elem_'+name)
                    el['visible']=False
    
    def UpdateElementsOnCanvas(self):
        for i in range(len(self.elements)):
            el = self.elements[i]
            if el['visible']==True:
                el['visible']=False
                self.DrawElement(el['name'])
        
    def ColorCells(self, clrs):
        cells = self.recGrid2d.cells
        nodesPoints = self.recGrid2d.nodesPoints
        for i in range(len(cells)):
            nld, nlu, nrd, nru = self.recGrid2d.Cell_GetCornerNodes(i)
            x0 = nodesPoints[nld].x
            y0 = nodesPoints[nld].y
            x1 = nodesPoints[nru].x
            y1 = nodesPoints[nru].y
            self.DrawRectangle(x0, y0, x1, y1, fill=clrs[i], outline=clrs[i])
            

    def SetWindowCorner(self, x0, y0):
        self.wx0 = x0
        self.wy0 = y0
        
    def SetWindowCenter(self, x0, y0):
        self.wx0 = x0 - self.wx/2.0
        self.wy0 = y0 - self.wy/2.0
        
    def ChangeWindowSize(self, r_x, r_y):
        ## zooms in or out by the ration r_x, r_y
        assert r_x>0.0 and r_y>0.0
        self.wx *= r_x
        self.wy *= r_y
        
    def ZoomInOnPoint(self, x, y):
        self.ChangeWindowSize(0.5, 0.5)
        self.SetWindowCenter(x, y)

    def ZoomOutOnPoint(self, x, y):
        self.ChangeWindowSize(2.0, 2.0)
        self.SetWindowCenter(x, y)

    def ZoomOnCenter(self):
        x, y = self.wx0 + self.wx/2.0, self.wy0 + self.wy/2.0
        self.ZoomInOnPoint(x, y)
        
    def ZoomOutOnCenter(self):
        x, y = self.wx0 + self.wx/2.0, self.wy0 + self.wy/2.0
        self.ZoomOutOnPoint(x, y)
        
        
        
        
        

