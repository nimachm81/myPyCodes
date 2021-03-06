## FourierBlochND.py
## N-Dimensional Fourier Analysis


__all__ = ["PDEFourierSeriesND"]


import sympy
from sympy import zeros, abc, Matrix, eye, Symbol, Sum, I, oo, Function,\
             expand, Mul, Add, Tuple, Integer, lambdify, KroneckerDelta, latex


import numpy as np
from scipy.sparse import csr_matrix, coo_matrix
from scipy import linalg, fftpack
#from cmath import *
import cmath

from IPython.display import display, Math, Latex

from Electromagnetics.SymExprTree import *
from Electromagnetics.Misc import tic, toc


class PDEFourierSeriesND:
    """ N-Dimensional Fourier series analysis of differential equations
    """
    #TODO: not tested for more than one eigen coefficient
    def __init__(self, expr, vars_fourier, x_str='x', X_str='X', n_str='n', n_dim=1, 
            mark='\\tilde', harmonic=None, harmonicMul=None, symmSpectrum=False, usexyz=True, vecRecip=None):
        assert n_dim>=0
        self.expr = expr
        self.varsfourier = vars_fourier
        self.x_str = x_str
        self.n_str = n_str
        self.N_Dim = n_dim
        self.harmonicMul = harmonicMul
        
        self.xs = [None]*self.N_Dim
        self.Xs = [None]*self.N_Dim
        self.ns = [None]*self.N_Dim
        self.dummys_str = ['m', 'p', 's']

        for i in range(self.N_Dim):
            self.xs[i] = Symbol(x_str + '_' + str(i))
            self.Xs[i] = Symbol(X_str + '_' + str(i))
            self.ns[i] = Symbol(n_str + '_' + str(i))
            
        if usexyz:
            if self.N_Dim==1:
                self.xs[0] = Symbol('x')
                self.Xs[0] = Symbol('X')
                if vecRecip==None:
                    self.ns[0] = Symbol('n_x')
            elif self.N_Dim==2:
                self.xs[0] = Symbol('x')
                self.xs[1] = Symbol('y')
                self.Xs[0] = Symbol('X')
                self.Xs[1] = Symbol('Y')
                if vecRecip==None:
                    self.ns[0] = Symbol('n_x')
                    self.ns[1] = Symbol('n_y')
            elif self.N_Dim==3:
                self.xs[0] = Symbol('x')
                self.xs[1] = Symbol('y')
                self.xs[2] = Symbol('z')
                self.Xs[0] = Symbol('X')
                self.Xs[1] = Symbol('Y')
                self.Xs[2] = Symbol('Z')
                if vecRecip==None:
                    self.ns[0] = Symbol('n_x')
                    self.ns[1] = Symbol('n_y')
                    self.ns[2] = Symbol('n_z')
        
        if harmonic:
            self.harmonic = harmonic
        else:
            self.harmonic = Integer(1)
            if vecRecip==None:
                for i in range(self.N_Dim):
                    self.harmonic = self.harmonic*sympy.exp(I*self.ns[i]*2*sympy.pi/self.Xs[i]*self.xs[i])
            else:
                assert len(vecRecip)==self.N_Dim
                for i in range(self.N_Dim):
                    b_i_dot_r = vecRecip[i][0]*self.xs[0]
                    for j in range(1, self.N_Dim):
                        b_i_dot_r += vecRecip[i][j]*self.xs[j]
                    self.harmonic = self.harmonic*sympy.exp(I*self.ns[i]*b_i_dot_r)
                            
        if harmonicMul!=None:
            ##Trouble: only vars take Bloch coefficients, pars do not..
            assert False
            self.harmonic = self.harmonic*harmonicMul
            
        self.harmonic = (self.harmonic).simplify()
        
        print('harmonic: ')
        display(Math(latex(self.harmonic)))
        
        self.harmonic_alis = Symbol('HARM')
        
        self.mark = mark
        self.setup_symmetric_spectrum(symmSpectrum)
        
        self.convProds = []
        self.n_pow_dic = None
        self.save_n_as_powers = True
        self.save_m_as_powers = True
        self.useFFTforConvs = True
        self.calculateDenseMatrices = False
        
        
    def putSums(self):
        #TODO should support x dependant functions (instead of Symbols) 
        # otherwise the time derivatives have to be performed after 
        # calling this function not before
        var_subs = [None]*len(self.varsfourier)
        for i in range(len(self.varsfourier)):
            var_subs[i] = Sum(Symbol(self.mark+'{'+self.varsfourier[i].name+'}')(*tuple(self.ns))
                *self.harmonic,(self.ns[0],-oo,+oo))
            for j in range(1, self.N_Dim):
                var_subs[i] = Sum(var_subs[i], (self.ns[j],-oo,+oo))
            
        
        #the following (commented) simple substitution does not work in the new version of sympy
        #expr_new = self.expr.subs([(self.varsfourier[i], var_subs[i]) for i in range(len(self.varsfourier))])
        
        expr_new = self.expr
        for i in range(len(self.varsfourier)):
            expr_new = symExp_replaceSymbol(expr_new, self.varsfourier[i], var_subs[i])
        

        self.varsHarm = [None]*len(self.varsfourier)
        for i in range(len(self.varsfourier)):
            self.varsHarm[i] = Function(self.mark+'{'+self.varsfourier[i].name+'}')
        return expr_new
        
    def createVarsHarm(self):
        self.varsHarm = [None]*len(self.varsfourier)
        for i in range(len(self.varsfourier)):
            self.varsHarm[i] = Function(self.mark+'{'+self.varsfourier[i].name+'}')
        
    def isFourierSeries(self, expr):
        if expr.func != Sum:
            return False
        args = expr.args
        if len(args)!=self.N_Dim+1:
            return False
        #if not args[0].has(self.harmonic):
        #    ##TODO: not safe - sometimes self.harmonic not recognized
        #    return False
        for i in range(self.N_Dim):
            if args[i+1]!=(self.ns[i], -oo, +oo):
                return False
        return True
        
        
    def getCoefficientPartOfSeries(self, expr):
        if expr.func != Sum:
            raise ValueError('expression is not a Fourier series!')
        args = expr.args
        if len(args)!=self.N_Dim+1:
            raise ValueError('expression is not a Fourier series!')
        #if not args[0].has(self.harmonic): ##TODO: to be fixed
        #    raise ValueError('expression is not a Fourier series!')
        for i in range(self.N_Dim):
            if args[i+1]!=(self.ns[i], -oo, +oo):
                raise ValueError('expression is not a Fourier series!')
        expr_coeff = (args[0]/self.harmonic).simplify()         ##TODO: time consuming
        return expr_coeff
        
        
    def MultiplySeries_Convulve(self, expr1, expr2):
        expr1_coeff = self.getCoefficientPartOfSeries(expr1)
        expr2_coeff = self.getCoefficientPartOfSeries(expr2)
        dummys = [None]*self.N_Dim
        for i in range(self.N_Dim):
            dummys[i] = Symbol(self.dummys_str[0] + '_' + str(i))
        expr1_coeff = expr1_coeff.subs([(self.ns[i], dummys[i]) for i in range(self.N_Dim)])
        expr2_coeff = expr2_coeff.subs([(self.ns[i], self.ns[i] - dummys[i]) for i in range(self.N_Dim)])
        
        expr_new_coeff = Sum(expr1_coeff*expr2_coeff, (dummys[0],-oo,oo))
        for j in range(1, self.N_Dim):
            expr_new_coeff = Sum(expr_new_coeff, (dummys[j],-oo,+oo))
        
        expr_new = Sum(expr_new_coeff*self.harmonic, (self.ns[0], -oo, oo))
        for j in range(1, self.N_Dim):
            expr_new = Sum(expr_new, (self.ns[j],-oo,+oo))

        return expr_new


    def applyConvolutionToNode(self, node):
        """Apply convolution to node if it contains multiplication of series
        """
        if node[2] == Mul:
            for i, arg_i in  enumerate(node[3]):
                if self.isFourierSeries(arg_i[1]):
                    for j, arg_j in enumerate(node[3]):
                        if i<j:
                            if self.isFourierSeries(arg_j[1]):
                                expr_conv = self.MultiplySeries_Convulve(arg_i[1], arg_j[1])
                                node[3][i][1] = expr_conv
                                del node[3][j]
                                node[1] = node[2](*tuple(node[3][k][1] for k in range(len(node[3]))))
                                symExpr_update_tree(node)
                                return True
        return False
        
        
    def applyConvolutions_walk_tree_INACTIVE(self, node):
        self.applyConvolutionToNode(node)
        for arg_i in node[3]:
            if len(arg_i)>0:
                self.applyConvolutions_walk_tree(arg_i)
        return

    def applyConvolutions_walk_tree(self, node):
        nodes_list = [node]
        ind_next = 0
        while True:
            self.applyConvolutionToNode(nodes_list[ind_next])
            node_next = nodes_list[ind_next]
            for arg in node_next[3]:
                if len(arg)>0:
                    nodes_list.append(arg)
            ind_next += 1
            if ind_next>=len(nodes_list):
                break
        return

    def applyConvolutions(self, expr):
        """Apply convolutions
        the equations should be expanded first in order to properly recognize
        the multiplication of 2 series
        """
        expr_expanded = expand(expr)
        expr_tree = symExpr_generate_tree(expr_expanded)
        self.applyConvolutions_walk_tree(expr_tree)
        return expr_tree[1]

    def aSumb_to_Sumab(self, node):
        """ If the expression at a given node is a*Sum(b*harmonic)
        it will be modified to Sum(a*b*harmonic)
        """
        if node[2]==Mul:
            for i, arg_i in enumerate(node[3]):
                if self.isFourierSeries(arg_i[1]):
                    coeff_in = self.getCoefficientPartOfSeries(arg_i[1])
                    mul_args_except_sum = [node[3][k][1] for k in range(len(node[3]))]
                    del mul_args_except_sum[i]
                    coeff_out = Mul(*tuple(mul_args_except_sum))
                    
                    node[1] = Sum(coeff_out*coeff_in*self.harmonic, (self.ns[0],-oo,oo))
                    for j in range(1, self.N_Dim):
                        node[1] = Sum(node[1], (self.ns[j],-oo,+oo))
                    symExpr_update_tree(node)
                    return True
        return False
        
    
    def aSumb_to_Sumab_walk_tree_INACTIVE(self, node):
        self.aSumb_to_Sumab(node)
        for arg_i in node[3]:
            if len(arg_i)>0:
                self.aSumb_to_Sumab_walk_tree(arg_i)
        return
        
    def aSumb_to_Sumab_walk_tree(self, node):
        nodes_list = [node]
        ind_next = 0
        while True:
            self.aSumb_to_Sumab(nodes_list[ind_next])
            node_next = nodes_list[ind_next]
            for arg in node_next[3]:
                if len(arg)>0:
                    nodes_list.append(arg)
            ind_next += 1
            if ind_next>=len(nodes_list):
                break
        return

    def SumaPlusSumb_to_SumaPlusb(self, node):
        """ If the expression at a given node is Sum(a*harmonic) + Sum(b*harmonic)
        it will be modified to Sum((a+b)*harmonic)
        """
        if node[2]==Add:
            for i, arg_i in enumerate(node[3]):
                if self.isFourierSeries(arg_i[1]):
                    for j, arg_j in enumerate(node[3]):
                        if i<j:
                            if self.isFourierSeries(arg_j[1]):
                                ##TODO: the simplify operation inside this function becomes expensive 
                                ## as the series argument grows
                                coeff_i = self.getCoefficientPartOfSeries(arg_i[1])     
                                coeff_j = self.getCoefficientPartOfSeries(arg_j[1])
                                ##
                                #print('len(node[3]):', len(node[3]), end=' ')
                                if len(node[3])==2:
                                    node[1] = Sum((coeff_i+coeff_j)*self.harmonic, (self.ns[0],-oo,oo))
                                    for k in range(1, self.N_Dim):
                                        node[1] = Sum(node[1], (self.ns[k],-oo,+oo))
                                    symExpr_update_tree(node)
                                    return 1
                                else:
                                    node_1_arg = [node[3][k][1] for k in range(len(node[3]))]
                                    del node_1_arg[j]
                                    del node_1_arg[i]
                                    sum_i_j = Sum((coeff_i+coeff_j)*self.harmonic, (self.ns[0],-oo,oo))
                                    for k in range(1, self.N_Dim):
                                        sum_i_j = Sum(sum_i_j, (self.ns[k],-oo,+oo))
                                    node_1_arg.append(sum_i_j)
                                    node[1] = node[2](*tuple(node_1_arg))
                                    symExpr_update_tree(node)
                                    ## the following line should be uncommented for 
                                    ## SumaPlusSumb_to_SumaPlusb_walk_tree_INACTIVE
                                    #self.SumaPlusSumb_to_SumaPlusb(node)
                                    return 2
        return 0
        

    def SumaPlusSumb_to_SumaPlusb_walk_tree_INACTIVE(self, node):
        ## modification in SumaPlusSumb_to_SumaPlusb is required for 
        ## proper operation
        self.SumaPlusSumb_to_SumaPlusb(node)
        for arg_i in node[3]:
            if len(arg_i)>0:
                self.SumaPlusSumb_to_SumaPlusb_walk_tree(arg_i)
        return


    def SumaPlusSumb_to_SumaPlusb_walk_tree(self, node):
        nodes_list = [node]
        ind_next = 0
        #tic()
        while True:
            res = self.SumaPlusSumb_to_SumaPlusb(nodes_list[ind_next])
            #print(toc(), end='  ')
            if res != 2:
                node_next = nodes_list[ind_next]
                for arg in node_next[3]:
                    if len(arg)>0:
                        nodes_list.append(arg)
                ind_next += 1
                if ind_next>=len(nodes_list):
                    break
        return


    def addSeriesManually(self, expr):
        expr_expanded = expand(expr)
        expr_tree = symExpr_generate_tree(expr_expanded)
        self.aSumb_to_Sumab_walk_tree(expr_tree)
        self.SumaPlusSumb_to_SumaPlusb_walk_tree(expr_tree)
        return expr_tree[1]


    def applyOrthogonalities(self, expr):
        expr = self.addSeriesManually(expr)
        if expr.func == Sum:
            return self.getCoefficientPartOfSeries(expr)
        elif expr == 0:
            return 0
        else:
            display(Math('\\text{expr} = ' + latex(expr)))
            print(expr)
            raise ValueError("Orthogonality was not assured..")
            return None



    def aSumb_to_Sumab_truncated(self, node):
        """ If the expression at a given node is a*Sum(b, limits)
        it will be modified to Sum(a*b, limits)
        truncated: limits are not infinity
        """
        dummys = [None]*self.N_Dim
        for i in range(self.N_Dim):
            dummys[i] = Symbol(self.dummys_str[0] + '_' + str(i))

        if node[2]==Mul:
            for i, arg_i in enumerate(node[3]):
                if arg_i[2]==Sum:
                    assert len(arg_i[3])==self.N_Dim+1
                    coeff_in = arg_i[3][0][1]
                    for j in range(self.N_Dim):
                        sum_limits_arg = arg_i[3][j+1][1]
                        assert sum_limits_arg.func==Tuple
                        dummy = sum_limits_arg.args[0].name    # 'n_1' or 'n_2' ...
                        #print(dummy, dummys)
                        assert Symbol(dummy) in dummys
                        #dummy_nodigits = ''.join([k for k in dummy if not (k.isdigit() or k=='_')])
                    mul_args_except_sum = [node[3][k][1] for k in range(len(node[3]))]
                    del mul_args_except_sum[i]
                    coeff_out = Mul(*tuple(mul_args_except_sum))
                    
                    for j in range(self.N_Dim):
                        sum_limits_arg = arg_i[3][j+1][1]
                        assert not coeff_out.has(sum_limits_arg.args[0])
                            
                    sum_limits_arg = arg_i[3][1][1]
                    node[1] = Sum(coeff_out*coeff_in, sum_limits_arg)
                    for j in range(1, self.N_Dim):
                        sum_limits_arg = arg_i[3][j+1][1]
                        node[1] = Sum(node[1], sum_limits_arg)
                    symExpr_update_tree(node)
                    return True
        return False

    def aSumb_to_Sumab_truncated_walk_tree_INACTIVE(self, node):
        self.aSumb_to_Sumab_truncated(node)
        for arg_i in node[3]:
            if len(arg_i)>0:
                self.aSumb_to_Sumab_truncated_walk_tree(arg_i)
        return
        
    def aSumb_to_Sumab_truncated_walk_tree(self, node):
        nodes_list = [node]
        ind_next = 0
        while True:
            self.aSumb_to_Sumab_truncated(nodes_list[ind_next])
            node_next = nodes_list[ind_next]
            for arg in node_next[3]:
                if len(arg)>0:
                    nodes_list.append(arg)
            ind_next += 1
            if ind_next>=len(nodes_list):
                break
        return

    def find_node_with_function(self, node, f):
        """ finds the first node whose node[2] component is the given function f and returns
            it starts from the given node and walks down the tree
        """
        if node[2]==f:
            return node
        else:
            for arg in node[3]:
                subnode_f = self.find_node_with_function(arg, f)
                if subnode_f!=False:
                    return subnode_f
        return False

    def find_node_with_symbol(self, node, x):
        """ finds the first node whose node[1] component is the given symbol x and returns
            it starts from the given node and walks down the tree
        """
        if node[1]==x:
            return node
        else:
            for arg in node[3]:
                subnode_x = self.find_node_with_symbol(arg, x)
                if subnode_x!=False:
                    return subnode_x
        return False


    def find_mult_coeff_of_func(self, node, var):
        """ 
        Finds the multiplicative coefficients of a given variable in an expression 
        linear in var.
        var is a Function. Its argumet may contain Symbol(n_str) or Symbol(n_str+'_1')..
        """
        if isinstance(node, list):
            if node[2]==var:
                var_coeff = 1
                var_arg = node[1].args
                var_ind = node[4]
                parent = node[0]
                while parent!=None:
                    if not(parent[2]==Mul or parent[2]==Add or parent[2]==Sum):
                        print(parent[2])
                        raise NotImplementedError("only Mul, Add are treated")
                    if parent[2]==Mul:
                        mul_args_parent = [parent[3][k][1] for k in range(len(parent[3])) if k!=var_ind]
                        var_coeff *= Mul(*tuple(mul_args_parent))
                    var_ind = parent[4]     ##++
                    parent = parent[0]
                coeff_arg = [var_coeff, var_arg]
                return coeff_arg
            return False
        else:
            expr_tree = symExpr_generate_tree(node)
            node_var = find_node_with_function(expr_tree, var)
            if node_var!=False:
                return find_mult_coeff_of_func(node_var, var)
            return False


    def getcoeff_setzero(self, node, var):
        """ var is a Function. Its argumet may contain Symbol(n_str) or Symbol(n_str+'_1')..
        """
        node_var = self.find_node_with_function(node, var)
        if node_var!=False:
             var_coeff_and_arg = self.find_mult_coeff_of_func(node_var, var)
             node_var[1] = Integer(0)
             symExpr_update_tree(node_var)
             return var_coeff_and_arg    # [var_coeff, var_arg]
        return False

    def remove_sums_of_zero_arg(self, node):
        """
        removes Sum functions whose argumet is zero
        """
        if node[2]==Sum:
            if node[1].args[0]==0:
                node[1] = Integer(0)    
                symExpr_update_tree(node)
        for arg in node[3]:
            if len(arg)>0:
                self.remove_sums_of_zero_arg(arg)
        return


    #-- numerical part
    
    def setup_symmetric_spectrum(self, sym):
        if sym==True:
            self.symmetric_spec = True
        elif sym==False:
            self.symmetric_spec = False
        else:
            raise ValueError('takes boolean input')
        return

    def setupNumericalParameters(self, expr_list, Ns, vars, pars, pars_vecs, eig_vars):
        """
           expr_list: sympy expression list (the given system of equations: expr_list[i]=0)
           N: number of Fourier harmonics (-N...N-1) or (-N...N)
           vars: list containing unknown fourier variables
           pars: list containing known Fourier parameters
           pars_vecs: list containig numpy arrays holding the coefficients of each variable
                  in pars
           eig_vars: list containing unknown scalars (usually the eigenvalue)
        """
        self.expr_list = expr_list
        self.Ns = Ns
        self.vars = vars
        self.pars = pars
        self.pars_vecs = pars_vecs
        self.eig_vars = eig_vars
        
        N_mp = 1
        for i in range(self.N_Dim):
            N = Ns[i]
            NF_0 = -N
            NF_1 = N
            if self.symmetric_spec==True:
                NF_1 = N + 1
            N_mp *= NF_1 - NF_0  #2*N  # Fourier indexes _N...N-1: 2*N (for speed)
        
        self.N_mp = N_mp
        self.n_total = N_mp*len(self.vars)
        if len(self.expr_list)!=len(self.vars):
            raise ValueError('number of equations should be equal to the number of \
                variables')
                        


    def orthogonalToNumpyMatrix(self):
        """It takes the orthogonal Fourier expression and generates a numpy matrix
            to solve.
        """
        expr_list = self.expr_list
        Ns = self.Ns
        pars = self.pars
        pars_vecs = self.pars_vecs
        eig_vars = self.eig_vars
        vars = self.vars
        
        self.convProds = []
         
        n_total = self.n_total

        ### matrices
        A_lin_sp = csr_matrix( (n_total,n_total), dtype=complex)
        A_eig_sp_list = []
        A_eig_dense_list = []
        A_convs = None
        if self.calculateDenseMatrices or not self.useFFTforConvs:
            A_convs = np.zeros((n_total, n_total), dtype=complex)
        b_rhs = np.zeros((n_total, 1))
        for eq_ind, expr in enumerate(expr_list):
            expr_expanded = expand(expr)
            expr_tree = symExpr_generate_tree(expr_expanded)
            self.aSumb_to_Sumab_truncated_walk_tree(expr_tree)
            for ind_v, var in enumerate(vars):
                while True:
                    coeff_arg = self.getcoeff_setzero(expr_tree, var)
                    #print('coeff_arg :')
                    #display(Math(latex(coeff_arg)))
                    if coeff_arg==False:
                        break
                    else:
                        coeff, arg = coeff_arg
                        for var_ in vars:
                            if coeff.has(var_):
                                raise ValueError('Equation is not linear.')
                        n_par = 0    # number of periodic parameters involved
                        for par in pars:
                            if coeff.has(par):
                                n_par += 1
                        n_eg = 0    # number of eigen variables involved
                        for eg in eig_vars:
                            if coeff.has(eg):
                                n_eg += 1
                        if n_par==0:
                            if n_eg==0:
                                a_lin_sp_coo = self.orth_ToNumpyMatrix_const_coeff(coeff,
                                       arg, ind_v, eq_ind*self.N_mp)
                                A_lin_sp = A_lin_sp + a_lin_sp_coo.tocsr()
                            else:
                                A_eig, a_eig_sp_coo = self.orth_toNumpyMatrix_eig_coeff(coeff,
                                      arg, ind_v, eq_ind*self.N_mp)
                                a_eig_sp_csr = a_eig_sp_coo.tocsr()
                                A_eig_sp_list.append([A_eig, a_eig_sp_csr])
                        elif n_par==1:
                            if self.calculateDenseMatrices or not self.useFFTforConvs:
                                A_eig_coeff, A_dense = self.orth_ToNumpyMatrix_convolutions(coeff, arg, 
                                     ind_v, eq_ind*self.N_mp)
                                if A_eig_coeff==None:
                                    A_convs += A_dense
                                else:
                                    A_eig_dense_list.append([A_eig_coeff, A_dense])
                            if self.useFFTforConvs: 
                                self.setupConvsForMatVec(coeff, arg, ind_v, eq_ind_start=eq_ind*self.N_mp)
                        else:
                            raise NotImplementedError('more than one periodic coefficient..')
            self.remove_sums_of_zero_arg(expr_tree)
            if expr_tree[1]!=0:
                # get the rhs
                print(expr_tree[1])
                display(Math(latex(expr_tree[1])))
                raise NotImplementedError('non-zero rhs not implemented..')
        return [A_lin_sp, A_convs, A_eig_sp_list, A_eig_dense_list, b_rhs]

    def increase_i_n_INACTIVE(self, i_n):
        Ns = self.Ns
        if self.symmetric_spec:
            for i in range(self.N_Dim+1):
                if i_n[i]<Ns[i]:
                    i_n[i] += 1
                    return True
                else:
                    i_n[i] = -Ns[i]
        else:
            for i in range(self.N_Dim):
                if i_n[i]<Ns[i]-1:
                    i_n[i] += 1
                    return True
                else:
                    i_n[i] = -Ns[i]
        return False

    def increase_i_n(self, i_n):
        Ns = self.Ns
        if self.symmetric_spec:
            for i in range(self.N_Dim-1, -1, -1):
                if i_n[i]<Ns[i]:
                    i_n[i] += 1
                    return True
                else:
                    i_n[i] = -Ns[i]
        else:
            for i in range(self.N_Dim-1, -1, -1):
                if i_n[i]<Ns[i]-1:
                    i_n[i] += 1
                    return True
                else:
                    i_n[i] = -Ns[i]
        return False

    def initialize_i_n(self):
        i_n = [0]*self.N_Dim
        for i in range(self.N_Dim):
            i_n[i] = -self.Ns[i]
        return i_n
    
    def eq_index_INACTIVE(self, i_n, eq_ind_start):
        Ns = self.Ns
        eqi_ind = i_n[0] + Ns[0]
        if self.symmetric_spec==True:
            for i in range(1, self.N_Dim):
                eqi_ind += (i_n[i] + Ns[i])*(2*Ns[i-1]+1)
        else:
            for i in range(1, self.N_Dim):
                eqi_ind += (i_n[i] + Ns[i])*2*Ns[i-1]
        return eq_ind_start + eqi_ind
    
    def eq_index(self, i_n, eq_ind_start):
        Ns = self.Ns
        eqi_ind = i_n[self.N_Dim-1] + Ns[self.N_Dim-1]
        if self.symmetric_spec==True:
            for i in range(self.N_Dim-2, -1, -1):
                eqi_ind += (i_n[i] + Ns[i])*(2*Ns[i+1]+1)
        else:
            for i in range(self.N_Dim-2, -1, -1):
                eqi_ind += (i_n[i] + Ns[i])*2*Ns[i+1]
        return eq_ind_start + eqi_ind
    
    def var_index_checkrange_INACTIVE(self, i_n, v):
        # v: ind_v
        Ns = self.Ns
        for i in range(self.N_Dim):
            NF_0 = -Ns[i]
            NF_1 = Ns[i]
            if self.symmetric_spec:
                NF_1 = Ns[i]+1
            if not NF_0<=i_n[i]<NF_1:
                return -1
        vi_ind = i_n[0] + Ns[0]
        if self.symmetric_spec==True:
            for i in range(1, self.N_Dim):
                vi_ind += (i_n[i] + Ns[i])*(2*Ns[i-1]+1)
        else:
            for i in range(1, self.N_Dim):
                vi_ind += (i_n[i] + Ns[i])*2*Ns[i-1]
        return v*self.N_mp + vi_ind

    def var_index_INACTIVE(self, i_n, v):
        # v: ind_v
        Ns = self.Ns
        vi_ind = i_n[0] + Ns[0]
        if self.symmetric_spec==True:
            for i in range(1, self.N_Dim):
                vi_ind += (i_n[i] + Ns[i])*(2*Ns[i-1]+1)
        else:
            for i in range(1, self.N_Dim):
                vi_ind += (i_n[i] + Ns[i])*2*Ns[i-1]
        return v*self.N_mp + vi_ind

    def var_index_checkrange(self, i_n, v):
        # v: ind_v
        Ns = self.Ns
        for i in range(self.N_Dim):
            NF_0 = -Ns[i]
            NF_1 = Ns[i]
            if self.symmetric_spec:
                NF_1 = Ns[i] + 1
            if not NF_0<=i_n[i]<NF_1:
                return -1
        vi_ind = i_n[self.N_Dim-1] + Ns[self.N_Dim-1]
        if self.symmetric_spec==True:
            for i in range(self.N_Dim-2, -1, -1):
                vi_ind += (i_n[i] + Ns[i])*(2*Ns[i+1]+1)
        else:
            for i in range(self.N_Dim-2, -1, -1):
                vi_ind += (i_n[i] + Ns[i])*2*Ns[i+1]
        return v*self.N_mp + vi_ind

    def var_index(self, i_n, v):
        # v: ind_v
        Ns = self.Ns
        vi_ind = i_n[self.N_Dim-1] + Ns[self.N_Dim-1]
        if self.symmetric_spec==True:
            for i in range(self.N_Dim-2, -1, -1):
                vi_ind += (i_n[i] + Ns[i])*(2*Ns[i+1]+1)
        else:
            for i in range(self.N_Dim-2, -1, -1):
                vi_ind += (i_n[i] + Ns[i])*2*Ns[i+1]
        return v*self.N_mp + vi_ind


    def orth_ToNumpyMatrix_const_coeff(self, coeff, arg, ind_v, eq_ind_start):
        """
        it takes a coefficient and argument (coeff*var(arg)) depennding only on 
        n_str and returns a sparse numpy matrix
        coeff*vars[var_ind](arg)
        coeff: coefficient expression
        arg: argument expression
        ind_v: index of the variable inside vars
        eq_ind_start: index of the first equation (row indices - 0-based)
        """
        Ns = self.Ns
            
        var_index = lambda i_n: self.var_index_checkrange(i_n, ind_v)
        
        #var_index_np = np.frompyfunc(var_index,1,1)

        eq_index = lambda i_n: self.eq_index(i_n, eq_ind_start)
            
        #eq_index_np = np.frompyfunc(eq_index,1,1) 
        increase_i_n = self.increase_i_n
            
        f_coeff = lambdify(tuple(self.ns), coeff)
        f_arg = lambdify(tuple(self.ns), arg)
        
        i_n = self.initialize_i_n()
            
        col = [-1]*self.N_mp
        row = [-1]*self.N_mp
        data = [-1]*self.N_mp
        
        ind = 0
        while True:
            col[ind] = var_index(list(f_arg(*tuple(i_n))))
            row[ind] = eq_index(i_n)
            data[ind] = f_coeff(*tuple(i_n))
            ind += 1
            if not increase_i_n(i_n):
                break
        assert ind==self.N_mp
        
        ind_to_remove = []
        for i in range(len(col)-1, -1, -1):
            if col[i]<0:
                del row[i]
                del col[i]
                del data[i]
        
        col = np.array(col)
        row = np.array(row)
        data = np.array(data)
        mat_sp = coo_matrix((data, (row,col)), shape=(self.n_total, self.n_total), dtype=complex)
        return mat_sp


    def split_expr(self, expr, symbs_1, symbs_2):
        """
        splits an expression into A(symbs_1)+B(symbs_2) or A(symbs_1)*B(symbs_2)
        symbs_1: list of symbols
        symbs_2: list of symbols
        expr: sympy expression
        """
        #TODO it works only if the top level functions is Add or Mul, 
        # the rest to be done
        if expr.func == Add or expr.func == Mul:
            A, B = [], []
            for arg in list(expr.args):
                toA, toB = True, True
                for x1 in symbs_1:
                    for x2 in symbs_2:
                        if arg.has(x1)==True:
                            toB = False
                        if arg.has(x2)==True:
                            toA = False
                if toA==True:
                    A.append(arg)
                elif toB==True:
                    B.append(arg)
                else:
                    return [False]
            return [True, expr.func, A, B]
        else:
            A_only = True
            for x in symbs_1:
                if expr.has(x):
                    A_only = False
                    break
            if A_only==True:
                return [True, Mul, [expr], [Integer(1)]]
            A_only = True
            for x in symbs_2:
                if expr.has(x):
                    A_only = False
                    break
            if A_only==True:
                return [True, Mul, [expr], [Integer(1)]]
                    
            raise NotImplementedError('Functions other than Add and Mul not imple\
            mented')
        return [False]


    def orth_toNumpyMatrix_eig_coeff(self, coeff, arg, ind_v, eq_ind_start):
        """
        it takes a coefficient and argument (coeff(eig_vars)*var(arg)), where 
        coeff is a function of eig_vars and n_str 
        and returns a sparse numpy matrix multiplied by an expression depending
        only on eig_vars as symbols
        coeff*vars[var_ind](arg)
        coeff: coefficient expression
        arg: argument expression
        n_str: fourier index (usually n)
        N: number of Fourier harmonics to keep
        n_total: total number of variables
        ind_v: index of the variable inside vars
        eq_ind_start: index of the first equation (row indices)
        eig_vars: eigenvalue parameters
        """
        eig_vars = self.eig_vars
        
        expr_split = self.split_expr(coeff.expand(), eig_vars, self.ns)
        if expr_split[0]==False:
            raise NotImplementedError('not implemented')
        if expr_split[1]==Add:
            raise NotImplementedError('not implemented: expression should have \
            been already expanded')
        ## coeff = A(eigs)*B(n)
        A_eigs = expr_split[1](*tuple(expr_split[2]))
        B_n = expr_split[1](*tuple(expr_split[3]))
        #print('A_eigs :', A_eigs)
        #print('B_n :', B_n)
        mat_sp = self.orth_ToNumpyMatrix_const_coeff(B_n, arg, ind_v, eq_ind_start)
        return [A_eigs, mat_sp] ## A_eigs*mat_sp



    def get_func_arg(self, expr, func):
        """
        finds the first instance of the func and returns its argument
        """
        expr_tree = symExpr_generate_tree(expr)
        node_func = self.find_node_with_function(expr_tree, func)
        if node_func!=False:
            return  node_func[1].args
        return False


    def orth_ToNumpyMatrix_convolutions(self, coeff, arg, ind_v, eq_ind_start):
        """
        it takes a coefficient and argument (coeff*var(arg)) depennding only on 
        n_str or n_str_i (a convolution sum) and returns a sparse numpy matrix
        coeff*vars[var_ind](arg)
        coeff: coefficient expression involving functions taking n, n_1 as arguments
        arg: argument expression (involving n and/or n_1)
        n_str: fourier index (usually n)
        N: number of Fourier harmonics to keep
        n_total: total number of variables
        ind_v: index of the variable inside vars
        eq_ind_start: index of the first equation (row indices)
        pars: list of periodic parameters
        pars_vecs: list of Fourier coefficients of periodic parameters
        """
        pars = self.pars 
        pars_vecs = self.pars_vecs
        Ns = self.Ns
        
        var_index = lambda i_n: self.var_index_checkrange(i_n, ind_v)
        
        #var_index_np = np.frompyfunc(var_index,1,1)

        eq_index = lambda i_n: self.eq_index(i_n, eq_ind_start)
            
        #eq_index_np = np.frompyfunc(eq_index,1,1) 
        
        increase_i_n = self.increase_i_n

        # get argumet and coeff of the periodic parameter
        p_coeff = None
        p_arg = None
        p_vec = None
        for i, p in enumerate(pars):
            if coeff.has(p):
                p_arg_list = self.get_func_arg(coeff, p)
                if p_arg_list==False or len(p_arg_list)==0:
                    print("p_arg_list: ", p_arg_list)
                    raise ValueError('par should be a periodic function')
                p_arg = p_arg_list
                p_coeff = symExp_replaceFunction(coeff, p, Symbol('AAAZZ'))\
                     .coeff(Symbol('AAAZZ'))
                p_vec = pars_vecs[i]
                break
        for p in pars:
            if p_coeff.has(p):
                raise NotImplementedError('Only one periodic parameter can be handled')

        # dummys                
        dummys = [None]*self.N_Dim
        for i in range(self.N_Dim):
            dummys[i] = Symbol(self.dummys_str[0] + '_' + str(i))

        # eig dependance
        has_eig = False
        for eig_var in self.eig_vars:
            if p_coeff.has(eig_var):
                has_eig = True
                break
                #raise NotImplementedError('the convolution term contains an eigen variable')
        p_coeff_eigs, p_coeff_nm = None, None
        if has_eig:
            expr_split = self.split_expr(p_coeff.expand(), self.eig_vars, self.ns+dummys)
            if expr_split[0]==False:
                raise NotImplementedError('not implemented')
            if expr_split[1]==Add:
                raise NotImplementedError('not implemented: expression should have been already expanded')
            ## p_coeff = A(eigs)*B(n,m)
            p_coeff_eigs = expr_split[1](*tuple(expr_split[2]))
            p_coeff_nm = expr_split[1](*tuple(expr_split[3]))

        # constructing functions
        f_p_coeff = None
        if has_eig==False:
            f_p_coeff = lambdify(tuple(self.ns + dummys), p_coeff)
        else:
            f_p_coeff = lambdify(tuple(self.ns + dummys), p_coeff_nm)
        f_p_arg = lambdify(tuple(self.ns + dummys), p_arg)
        f_arg = lambdify(tuple(self.ns + dummys), arg)
        
        #print('arg:', arg, '\np_arg:', p_arg, '\np_coeff:', p_coeff)
        
        #print('arg = ', arg)
        # constructing the dense coefficien2t matrix
        A_mat = np.zeros((self.n_total, self.n_total), dtype=complex)

        i_n = self.initialize_i_n()

        while True:
            row_ind = eq_index(i_n)
            
            i_m = self.initialize_i_n()
            while True:
                __c = f_p_coeff(*tuple(i_n+i_m))
                _ind_v_ = f_arg(*tuple(i_n+i_m))
                _ind_p_ = f_p_arg(*tuple(i_n+i_m))
                #print('i_n: ', i_n, '   i_m:', i_m, '   _ind_v_:', _ind_v_, '   _ind_p_:', _ind_p_, '   ind_v: ', ind_v)
                
                is_in_range = True
                for i in range(self.N_Dim):
                    NF_0 = -Ns[i]
                    NF_1 = Ns[i]
                    if self.symmetric_spec:
                        NF_1 = Ns[i]+1
                    if not(NF_0<=_ind_v_[i]<NF_1 and NF_0<=_ind_p_[i]<NF_1):
                        is_in_range = False
                        break
                if is_in_range:
                    col_ind = var_index(list(_ind_v_))
                    _ind_p_vec = [0]*self.N_Dim
                    for i in range(self.N_Dim):
                        _ind_p_vec[i] = _ind_p_[i] + Ns[i]
                    #print('[row, col]: ', [row_ind, col_ind], '  _ind_p_vec:', _ind_p_vec)
                    A_mat[row_ind, col_ind] += __c*p_vec[tuple(_ind_p_vec)]
                    #print('A: ', A_mat[row_ind, col_ind], '  p_vec: ', p_vec[tuple(_ind_p_vec)])
                if not increase_i_n(i_m):
                    break

            if not increase_i_n(i_n):
                break
        return [p_coeff_eigs, A_mat]
    
    
    def setupConvsForMatVec(self, coeff, arg, ind_v, eq_ind_start):
        """ setup convolutions for matrix vector product
            it takes the variable and its coefficient (including a aperiodic parameter)
            and defines the appropriate variables for matrix vector product using FFTs
        """
        pars = self.pars 
        pars_vecs = self.pars_vecs
        Ns = self.Ns

        # get argumet and coeff of the periodic parameter
        p_coeff = None
        p_arg = None
        p_vec = None
        for i, p in enumerate(pars):
            if coeff.has(p):
                p_arg_list = self.get_func_arg(coeff, p)
                if p_arg_list==False or len(p_arg_list)==0:
                    print("p_arg_list: ", p_arg_list)
                    raise ValueError('par should be a periodic function')
                p_arg = p_arg_list
                p_coeff = symExp_replaceFunction(coeff, p, Symbol('AAAZZ'))\
                     .coeff(Symbol('AAAZZ'))
                p_vec = pars_vecs[i]
                break
        for p in pars:
            if p_coeff.has(p):
                raise NotImplementedError('Only one periodic parameter can be handled')

        # dummys                
        dummys = [None]*self.N_Dim
        for i in range(self.N_Dim):
            dummys[i] = Symbol(self.dummys_str[0] + '_' + str(i))

        # eig dependance
        has_eig = False
        for eig_var in self.eig_vars:
            if p_coeff.has(eig_var):
                has_eig = True
                break
                #raise NotImplementedError('the convolution term contains an eigen variable')
        p_coeff_eigs, p_coeff_nm = None, None
        if has_eig:
            expr_split = self.split_expr(p_coeff.expand(), self.eig_vars, self.ns+dummys)
            if expr_split[0]==False:
                raise NotImplementedError('not implemented')
            if expr_split[1]==Add:
                raise NotImplementedError('not implemented: expression should have been already expanded')
            ## p_coeff = A(eigs)*B(n,m)
            p_coeff_eigs = expr_split[1](*tuple(expr_split[2]))
            p_coeff_nm = expr_split[1](*tuple(expr_split[3]))

        if p_coeff_nm==None:
            p_coeff_nm = p_coeff
        
        # n and m dependance
        has_n, has_m = False, False
        for i in range(self.N_Dim):
            if p_coeff_nm.has(self.ns[i]):
                has_n = True
                break
        for i in range(self.N_Dim):
            if p_coeff_nm.has(dummys[i]):
                has_m = True
                break
        
        p_coeff_n, p_coeff_m, p_coeff_const = None, None, None
        if has_n and has_m:
            expr_split = self.split_expr(p_coeff_nm.expand(), self.ns, dummys)
            if expr_split[0]==False:
                raise NotImplementedError('not implemented')
            if expr_split[1]==Add:
                raise NotImplementedError('not implemented: expression should have been already expanded')
            ## p_coeff_nm = A(n)*B(m)
            p_coeff_n = expr_split[1](*tuple(expr_split[2]))
            p_coeff_m = expr_split[1](*tuple(expr_split[3]))
        elif has_n:
            p_coeff_n = p_coeff_nm
        elif has_m:
            p_coeff_m = p_coeff_nm
        else:
            p_coeff_const = p_coeff_nm
            
        # get the corressponding sparse (diagonal) matrix multiplier if p_coeff_n != None
        eq_index = lambda i_n: self.eq_index(i_n, eq_ind_start)

        increase_i_n = self.increase_i_n
        
        ## p_coeff_n
        POW_MAX_ = 5
        mat_sp_coeff_n = None
        n_is_power_coeff = False
        n_powers = [0]*self.N_Dim
        if p_coeff_n!=None:
            for n_ind in range(self.N_Dim):
                n = self.ns[n_ind]
                if p_coeff_n.has(n):
                    n_pow = 1
                    coeff_div_n = (p_coeff_n/n).simplify()
                    is_pow_form = True
                    for i in range(POW_MAX_):
                        if coeff_div_n.has(n):
                            n_pow += 1
                            coeff_div_n = (coeff_div_n/n).simplify()
                            is_pow_form = False
                        else:
                            is_pow_form = True
                            break
                    if is_pow_form:
                        n_powers[n_ind] = n_pow
                        n_is_power_coeff = True
                    else:
                        display(Math(latex(p_coeff_n)))
                        n_is_power_coeff = False
            if n_is_power_coeff:
                ns_ = Integer(1)
                for i in range(self.N_Dim):
                    ns_ = ns_ * self.ns[i]**n_powers[i]
                p_coeff_n_const = (p_coeff_n/ns_).simplify()
                for i in range(self.N_Dim):
                    assert p_coeff_n_const.has(self.ns[i])==False

                p_coeff_n = ns_
                if p_coeff_const==None:
                    p_coeff_const = p_coeff_n_const
                else:
                    p_coeff_const = p_coeff_const*p_coeff_n_const

            if has_n and (n_is_power_coeff==False or self.save_n_as_powers==False):
                #display(Math(latex(p_coeff_n)))
                f_p_coeff_n = lambdify(tuple(self.ns), p_coeff_n)
                
                i_n = self.initialize_i_n()

                rows = [-1]*self.N_mp
                cols = [-1]*self.N_mp
                data = [-1]*self.N_mp
                ind_ = 0
                while True:
                    row_ind = eq_index(i_n)
                    
                    rows[ind_] = row_ind
                    cols[ind_] = row_ind
                    data[ind_] = f_p_coeff_n(*tuple(i_n))
                    ind_ += 1
                    
                    if not increase_i_n(i_n):
                        break
                mat_sp_coeff_n = coo_matrix((data, (rows,cols)), shape=(self.n_total, self.n_total), dtype=complex)

        ## p_coeff_m
        m_is_power_coeff = False
        m_powers = [0]*self.N_Dim
        arr_coeff_m = None
        if p_coeff_m!=None:
            for m_ind in range(self.N_Dim):
                m = dummys[m_ind]
                if p_coeff_m.has(m):
                    m_pow = 1
                    coeff_div_m = (p_coeff_m/m).simplify()
                    is_pow_form = True
                    for i in range(POW_MAX_):
                        if coeff_div_m.has(m):
                            m_pow += 1
                            coeff_div_m = (coeff_div_m/m).simplify()
                            is_pow_form = False
                        else:
                            is_pow_form = True
                            break
                    if is_pow_form:
                        m_powers[m_ind] = m_pow
                        m_is_power_coeff = True
                    else:
                        display(Math(latex(p_coeff_m)))
                        m_is_power_coeff = False
            if m_is_power_coeff:
                ms = Integer(1)
                for i in range(self.N_Dim):
                    ms = ms * dummys[i]**m_powers[i]
                p_coeff_m_const = (p_coeff_m/ms).simplify()
                for i in range(self.N_Dim):
                    assert p_coeff_m_const.has(dummys[i])==False
                
                p_coeff_m = ms
                if p_coeff_const==None:
                    p_coeff_const = p_coeff_m_const
                else:
                    p_coeff_const = p_coeff_const*p_coeff_m_const

            if has_m and (m_is_power_coeff==False or self.save_m_as_powers==False):
                #display(Math(latex(p_coeff_m)))
                f_p_coeff_m = lambdify(tuple(dummys), p_coeff_m)
                
                i_m = self.initialize_i_n()

                arr_shape = [-1]*self.N_Dim
                for i in range(self.N_Dim):
                    if self.symmetric_spec:
                        arr_shape[i] = 2*self.Ns[i]+1
                    else:
                        arr_shape[i] = 2*self.Ns[i]
                arr_shape = tuple(arr_shape)

                arr_coeff_m = np.zeros(arr_shape, dtype=complex)
                _ind_m_arr = [0]*self.N_Dim
                N_Dim = self.N_Dim
                while True:
                    for i in range(N_Dim):
                        _ind_m_arr[i] = i_m[i] + Ns[i]
                    arr_coeff_m[tuple(_ind_m_arr)] = f_p_coeff_m(*tuple(i_m))
                    if not increase_i_n(i_m):
                        break
                        
        ## p_coeff_eigs
        if p_coeff_eigs!=None and len(self.eig_vars)==1:
            if (p_coeff_eigs/self.eig_vars[0]).simplify().has(self.eig_vars[0]) == False:
                p_coeff_eigs = self.eig_vars[0]
                if p_coeff_const==None:
                    p_coeff_const = (p_coeff_eigs/self.eig_vars[0]).simplify()
                else:
                    p_coeff_const = p_coeff_const*(p_coeff_eigs/self.eig_vars[0]).simplify()

        ## arranging different terms                
        n_coeff_term = None
        if n_is_power_coeff and self.save_n_as_powers:
            n_coeff_term = [n_is_power_coeff, n_powers]
        else:
            #assert mat_sp_coeff_n!=None
            n_coeff_term = [False, mat_sp_coeff_n]
        
        m_coeff_term = None
        if m_is_power_coeff and self.save_m_as_powers:
            m_coeff_term = [m_is_power_coeff, m_powers]
        else:
            #assert arr_coeff_m!=None
            #print('arr_coeff_m: \n', arr_coeff_m)
            m_coeff_term = [False, arr_coeff_m]
        
        if p_coeff_const!=None:
            #display(Math('\\text{p_coeff_const}:'+latex(p_coeff_const)))
            p_coeff_const = complex(p_coeff_const)
        
        self.convProds.append([[p_coeff_eigs, n_coeff_term, m_coeff_term, p_coeff_const],
                [arg, p_arg], [ind_v, eq_ind_start], p_vec])
        
        
    def set_powers_of_n(self):
        n_pow_dic = {}
        for conv in self.convProds:
            n_coeff_term = conv[0][1]
            eq_ind_start = conv[2][1]
            if n_coeff_term[0]==False:
                continue
            n_pows = n_coeff_term[1]
            if eq_ind_start in n_pow_dic:
                pows_list = n_pow_dic[eq_ind_start]
                for i in range(self.N_Dim):
                    if n_pows[i]!=0:
                        pow_i_exists = False
                        for j in range(len(pows_list[i])):
                            if pows_list[i][j][0]==n_pows[i]:
                                pow_i_exists = True
                                break
                        if pow_i_exists==False:
                            pows_list[i].append([n_pows[i], None])
            else:
                n_pow_dic[eq_ind_start] = []
                for i in range(self.N_Dim):
                    n_pow_dic[eq_ind_start].append([])
                pows_list = n_pow_dic[eq_ind_start]
                for i in range(self.N_Dim):
                    if n_pows[i]!=0:
                        pow_i_exists = False
                        for j in range(len(pows_list[i])):
                            if pows_list[i][j][0]==n_pows[i]:
                                pow_i_exists = True
                                break
                        if pow_i_exists==False:
                            pows_list[i].append([n_pows[i], None])
        self.n_pow_dic = n_pow_dic             
            

    def set_powers_of_n_mats(self):
        for eq_ind_start in self.n_pow_dic:
            pows_list = self.n_pow_dic[eq_ind_start]
            for i in range(self.N_Dim):
                n_powers = [0]*self.N_Dim
                for j in range(len(pows_list[i])):
                    assert pows_list[i][j][0]!=0
                    n_powers[i] = pows_list[i][j][0]
                    assert pows_list[i][j][1]==None
                    mat_sp_coeff_n = self.get_n_coeff_Mat(n_powers, eq_ind_start)
                    assert mat_sp_coeff_n!=None
                    pows_list[i][j][1] = mat_sp_coeff_n
                    

    def get_n_coeff_Mat(self, n_powers, eq_ind_start):
        """ get the matrix corressponding to multiplication by powers of self.ns
        """
        Ns = self.Ns

        eq_index = lambda i_n: self.eq_index(i_n, eq_ind_start)

        increase_i_n = self.increase_i_n

        if n_powers==[0]*self.N_Dim:
            return None        
        
        ns_ = Integer(1)
        for i in range(self.N_Dim):
            ns_ = ns_ * self.ns[i]**n_powers[i]

        f_pow_n = lambdify(tuple(self.ns), ns_)
        
        i_n = self.initialize_i_n()

        rows = [-1]*self.N_mp
        cols = [-1]*self.N_mp
        data = [-1]*self.N_mp
        ind_ = 0
        while True:
            row_ind = eq_index(i_n)
            
            rows[ind_] = row_ind
            cols[ind_] = row_ind
            data[ind_] = f_pow_n(*tuple(i_n))
            ind_ += 1
            
            if not increase_i_n(i_n):
                break
        mat_sp_coeff_n = coo_matrix((data, (rows,cols)), shape=(self.n_total, self.n_total), dtype=complex)
        return mat_sp_coeff_n


    def multiplyArrayByPowersOfm(self, arr, m_powers):
        """ multipy a self.N_Dim dimensional array by powers of m appearing inside
            the convolution
        """
        if self.N_Dim==1:
            if m_powers[0]!=0:
                if self.symmetric_spec:
                    return arr*np.power(np.arange(-self.Ns[0], self.Ns[0]+1), m_powers[0])
                else:
                    return arr*np.power(np.arange(-self.Ns[0], self.Ns[0]), m_powers[0])
            else:
                return arr
        elif self.N_Dim==2:
            A = arr.copy()
            if m_powers[0]!=0:
                if self.symmetric_spec:
                    for i in range(2*self.Ns[1]+1):
                        A[:,i] = A[:,i]*np.power(np.arange(-self.Ns[0], self.Ns[0]+1), m_powers[0])
                else:
                    for i in range(2*self.Ns[1]):
                        A[:,i] = A[:,i]*np.power(np.arange(-self.Ns[0], self.Ns[0]), m_powers[0])
            if m_powers[1]!=0:
                if self.symmetric_spec:
                    for i in range(2*self.Ns[0]+1):
                        A[i,:] = A[i,:]*np.power(np.arange(-self.Ns[1], self.Ns[1]+1), m_powers[1])
                else:
                    for i in range(2*self.Ns[0]):
                        A[i,:] = A[i,:]*np.power(np.arange(-self.Ns[1], self.Ns[1]), m_powers[1])
            return A
        elif self.N_Dim==3:
            A = arr.copy()
            if m_powers[0]!=0:
                if self.symmetric_spec:
                    for i in range(2*self.Ns[1]+1):
                        for j in range(2*self.Ns[2]+1):
                            A[:,i,j] = A[:,i,j]*np.power(np.arange(-self.Ns[0], self.Ns[0]+1), m_powers[0])
                else:
                    for i in range(2*self.Ns[1]):
                        for j in range(2*self.Ns[2]):
                            A[:,i,j] = A[:,i,j]*np.power(np.arange(-self.Ns[0], self.Ns[0]), m_powers[0])
            if m_powers[1]!=0:
                if self.symmetric_spec:
                    for i in range(2*self.Ns[0]+1):
                        for j in range(2*self.Ns[2]+1):
                            A[i,:,j] = A[i,:,j]*np.power(np.arange(-self.Ns[1], self.Ns[1]+1), m_powers[1])
                else:
                    for i in range(2*self.Ns[0]):
                        for j in range(2*self.Ns[2]):
                            A[i,:,j] = A[i,:,j]*np.power(np.arange(-self.Ns[1], self.Ns[1]), m_powers[1])
            if m_powers[2]!=0:
                if self.symmetric_spec:
                    for i in range(2*self.Ns[0]+1):
                        for j in range(2*self.Ns[1]+1):
                            A[i,j,:] = A[i,j,:]*np.power(np.arange(-self.Ns[2], self.Ns[2]+1), m_powers[2])
                else:
                    for i in range(2*self.Ns[0]):
                        for j in range(2*self.Ns[1]):
                            A[i,j,:] = A[i,j,:]*np.power(np.arange(-self.Ns[2], self.Ns[2]), m_powers[2])
            return A
        else:
            raise NotImplementedError()
            

    def resetConvProds(self):
        self.convProds = []
        

    def wrapVariableIndices_INACTIVE(self, x_vec):
        """ it takes the unwraped final variables vector and wraps it to 
            the corressponding multidimensional array for each variable
        """
        Ns = self.Ns
        N_Dim = self.N_Dim
        
        N_t = [0]*N_Dim
        for i in range(N_Dim):
            if self.symmetric_spec:
                N_t[i] = 2*Ns[i] + 1
            else:
                N_t[i] = 2*Ns[i]

        n_vars = len(self.vars)
        vars_ND_array = [None]*n_vars
        for i in range(n_vars):
            vars_ND_array[i] = np.zeros(tuple(N_t), dtype=complex)
        
        var_index = self.var_index
        
        increase_i_n = self.increase_i_n

        if len(x_vec.shape)==1:
            for v in range(n_vars):
                v_arr = vars_ND_array[v]
            
                i_n = self.initialize_i_n()
                
                while True:
                    i_n_v = [0]*self.N_Dim
                    for i in range(self.N_Dim):
                        i_n_v[i] = i_n[i] + Ns[i]
                    v_arr[tuple(i_n_v)] = x_vec[var_index(i_n, v)]
                    if not increase_i_n(i_n):
                        break
        elif len(x_vec.shape)==2:
            assert x_vec.shape[1]==1
            for v in range(n_vars):
                v_arr = vars_ND_array[v]
            
                i_n = self.initialize_i_n()
                
                while True:
                    i_n_v = [0]*self.N_Dim
                    for i in range(self.N_Dim):
                        i_n_v[i] = i_n[i] + Ns[i]
                    v_arr[tuple(i_n_v)] = x_vec[var_index(i_n, v), 0]
                    if not increase_i_n(i_n):
                        break
        else:
            raise ValueError()
            
        self.vars_ND_array = vars_ND_array

    def wrapVariableIndices(self, x_vec):
        """ it takes the unwraped final variables vector and wraps it to 
            the corressponding multidimensional array for each variable
        """
        ## this wrapping is correct only if indices are treated like C arrays
        ## in var_index and eq_index methods
        x_vec_arr_shape = [len(self.vars)]
        for i in range(self.N_Dim):
            if self.symmetric_spec:
                x_vec_arr_shape.append(2*self.Ns[i]+1)
            else:
                x_vec_arr_shape.append(2*self.Ns[i])
        self.vars_ND_array = x_vec.reshape(tuple(x_vec_arr_shape))


    def getConvProductsFFT(self, x_vec, eig_var_vals=None):
        """ uses FFT for calculating matrix-vector products for convolution terms
        """
        N_Dim = self.N_Dim
        Ns = self.Ns
        
        eq_index = self.eq_index
        
        increase_i_n = self.increase_i_n
            

        if self.n_pow_dic==None:
            self.set_powers_of_n()
            self.set_powers_of_n_mats()
        
        self.wrapVariableIndices(x_vec)        
        
        dummys = [None]*self.N_Dim
        for i in range(self.N_Dim):
            dummys[i] = Symbol(self.dummys_str[0] + '_' + str(i))

        arg_ms, arg_ns_ms = [0]*self.N_Dim, [0]*self.N_Dim
        for i in range(self.N_Dim):
            arg_ms[i] = dummys[i]
            arg_ns_ms[i] = self.ns[i] - dummys[i]
        arg_ms = tuple(arg_ms)
        arg_ns_ms = tuple(arg_ns_ms)
        
        assert x_vec.shape[0]==self.n_total
        y_tot = np.zeros(x_vec.shape, dtype=complex)
        y_tot_eigs = []
        
        for conv_term in self.convProds:
            coeffs__, args__, inds__, p_vec = conv_term
            p_coeff_eigs, n_coeff_term, m_coeff_term, p_coeff_const = coeffs__
            arg, p_arg = args__ 
            ind_v, eq_ind_start = inds__
        
            y_i = np.zeros(x_vec.shape, dtype=complex)
        
            arr_coeff = self.vars_ND_array[ind_v]
            arr_par = p_vec
            if arg==arg_ms and p_arg==arg_ns_ms:
                if m_coeff_term[0]==True:
                    m_powers = m_coeff_term[1]
                    if m_powers!=[0]*self.N_Dim:
                        arr_coeff = self.multiplyArrayByPowersOfm(arr_coeff, m_powers)
                else:
                    arr_coeff_m = m_coeff_term[1]
                    if arr_coeff_m!=None:
                        arr_coeff = arr_coeff * arr_coeff_m
            elif p_arg==arg_ms and arg==arg_ns_ms:
                if m_coeff_term[0]==True:
                    m_powers = m_coeff_term[1]
                    if m_powers!=[0]*self.N_Dim:
                        arr_par = self.multiplyArrayByPowersOfm(arr_par, m_powers)
                else:
                    arr_coeff_m = m_coeff_term[1]
                    if arr_coeff_m!=None:
                        arr_par = arr_par * arr_coeff_m
            else:
                display(Math('arg: ' + latex(arg)))
                display(Math('p_arg: ' + latex(p_arg)))
                display(Math('arg_ms: ' + latex(arg_ms)))
                display(Math('arg_ns_ms: ' + latex(arg_ns_ms)))
                raise NotImplementedError()

            arr_mul = self.getInverseFourierCoeffsNoShift(arr_coeff) * self.getInverseFourierCoeffsNoShift(arr_par)
            #arr_mul = self.getInverseFourierCoeffs(arr_coeff) * self.getInverseFourierCoeffs(arr_par)
            arr_mul = self.getFourierCoeffs_(arr_mul)
            
            if p_coeff_const!=None:
                arr_mul *= p_coeff_const
            
            
            if len(x_vec.shape)==1:
                y_i[eq_ind_start:eq_ind_start+self.N_mp] = arr_mul.reshape((self.N_mp,))
            elif len(x_vec.shape)==2:
                assert x_vec.shape[1]==1
                y_i[eq_ind_start:eq_ind_start+self.N_mp, 0] = arr_mul.reshape((self.N_mp,))
            else:
                raise ValueError()
            
            """
            i_n = self.initialize_i_n()
            i_n_arr = [0]*N_Dim
            if len(x_vec.shape)==1:
                while True:
                    for i in range(N_Dim):
                        i_n_arr[i] = i_n[i] + Ns[i]
                    y_i[eq_index(i_n, eq_ind_start)] = arr_mul[tuple(i_n_arr)]
                    if not increase_i_n(i_n):
                        break
            elif len(x_vec.shape)==2:
                assert x_vec.shape[1]==1
                while True:
                    for i in range(N_Dim):
                        i_n_arr[i] = i_n[i] + Ns[i]
                    y_i[eq_index(i_n, eq_ind_start), 0] = arr_mul[tuple(i_n_arr)]
                    if not increase_i_n(i_n):
                        break
            else:
                raise ValueError()
            """
            
            if n_coeff_term[0]==True:
                n_powers = n_coeff_term[1]
                if n_powers!=[0]*self.N_Dim:
                    for i in range(self.N_Dim):
                        p_i = n_powers[i]
                        ps_i = self.n_pow_dic[eq_ind_start][i]  # precalculated powers
                        mat_n = None
                        for j in range(len(ps_i)):
                            if ps_i[j][0]==p_i:
                                mat_n = ps_i[j][1]
                                break
                        assert mat_n!=None
                        y_i = mat_n.dot(y_i)
            else:
                mat_sp_coeff_n = n_coeff_term[1]
                if mat_sp_coeff_n != None:
                    y_i = mat_sp_coeff_n.dot(y_i)
            
            p_coeff_eigs_val = None
            if p_coeff_eigs!=None:
                f_p_coeff_eigs = lambdify(tuple(self.eig_vars), p_coeff_eigs)
                if eig_var_vals!=None:
                    p_coeff_eigs_val = f_p_coeff_eigs(*tuple(eig_var_vals))

            if p_coeff_eigs_val!=None:
                y_tot += p_coeff_eigs_val*y_i
            elif p_coeff_eigs==None:
                y_tot += y_i
            else:
                y_tot_eigs.append(y_i, p_coeff_eigs)
            
        return [y_tot, y_tot_eigs]
            
            

    """
    def getMatrix(self, eigvar_vals, A_eqs):
        eig_vars = self.eig_vars
        A_lin_sp, A_convs, A_eig_sp_list, b_rhs = A_eqs
        A_tot = A_convs.copy()
        if A_lin_sp.getnnz()>0:
            A_tot += A_lin_sp
        for i in range(len( A_eig_sp_list)):
            C = complex(A_eig_sp_list[i][0].evalf(subs=dict(zip(eig_vars, eigvar_vals))))
            A_tot +=  C*A_eig_sp_list[i][1]
        return A_tot
    """

    def getMatrices(self, eigvar_vals, A_eqs_list):
        if not self.calculateDenseMatrices:
            raise ValueError('set self.calculateDenseMatrices = True !')
        eig_vars = self.eig_vars
        mats = []
        for A_eqs in A_eqs_list:
            A_lin_sp, A_convs, A_eig_sp_list, A_eig_dense_list, b_rhs = A_eqs
            A_tot = A_convs.copy()
            if A_lin_sp.getnnz()>0:
                A_tot += A_lin_sp
            for i in range(len( A_eig_sp_list)):
                C = complex(A_eig_sp_list[i][0].evalf(subs=dict(zip(eig_vars, eigvar_vals))))
                A_tot +=  C*A_eig_sp_list[i][1]
            for i in range(len( A_eig_dense_list)):
                C = complex(A_eig_dense_list[i][0].evalf(subs=dict(zip(eig_vars, eigvar_vals))))
                A_tot +=  C*A_eig_dense_list[i][1]
            mats.append(A_tot)
        return mats


    def getDeterminant(self, eigvar_vals, A_eqs_list):
        """
        find the determinant of the eigenvalue system where the eigenvalues are
        given in eigvar_vals
        """
        if not self.calculateDenseMatrices:
            raise ValueError('set self.calculateDenseMatrices = True !')
        eig_vars = self.eig_vars
        dets = []
        for A_eqs in A_eqs_list:
            A_lin_sp, A_convs, A_eig_sp_list, A_eig_dense_list, b_rhs = A_eqs
            A_tot = A_convs.copy()
            if A_lin_sp.getnnz()>0:
                A_tot += A_lin_sp
            for i in range(len( A_eig_sp_list)):
                C = complex(A_eig_sp_list[i][0].evalf(subs=dict(zip(eig_vars, eigvar_vals))))
                A_tot +=  C*A_eig_sp_list[i][1]
            for i in range(len( A_eig_dense_list)):
                C = complex(A_eig_dense_list[i][0].evalf(subs=dict(zip(eig_vars, eigvar_vals))))
                A_tot +=  C*A_eig_dense_list[i][1]
            dets.append(linalg.det(A_tot))
        return dets
        


    def getLogDeterminant(self, eigvar_vals, A_eqs_list):
        """
        find the logarithm of the determinant of the eigenvalue system where 
        the eigenvalues are given in eigvar_vals
        """
        if not self.calculateDenseMatrices:
            raise ValueError('set self.calculateDenseMatrices = True !')
        eig_vars = self.eig_vars
        dets = []
        for A_eqs in A_eqs_list:
            A_lin_sp, A_convs, A_eig_sp_list, A_eig_dense_list, b_rhs = A_eqs
            A_tot = A_convs.copy()
            if A_lin_sp.getnnz()>0:
                A_tot += A_lin_sp
            for i in range(len( A_eig_sp_list)):
                C = complex(A_eig_sp_list[i][0].evalf(subs=dict(zip(eig_vars, eigvar_vals))))
                A_tot +=  C*A_eig_sp_list[i][1]
            for i in range(len( A_eig_dense_list)):
                C = complex(A_eig_dense_list[i][0].evalf(subs=dict(zip(eig_vars, eigvar_vals))))
                A_tot +=  C*A_eig_dense_list[i][1]
            P, L, U = plu = linalg.lu(A_tot)
            log_det = 0.0
            for i in range(A_tot.shape[0]):
                if U[i,i]==0.0:
                    print('U[i,i]==0.0')
                    #print(A_tot)
                    log_det = -300.0    ## TODO: handle log(0)
                    break
                log_det += cmath.log((U[i,i]))
            log_det += cmath.log(linalg.det(P))
            dets.append(log_det)
        return dets

        
        
    def solveDeterminant(self, eigvar_vals_0, A_eqs_list,
            solver='lm', handle_overflow=False, roots_prev=None, 
            tol=1.0e-10, tol_relative=True, maxiter=1000, ftol=1.e-5, maxfev=1000,
            getMats=False):
        """
        solves the determinant for eigenvalues
        x_indep: independant sympy variables (list), 
        x_dep: dependant sympy variables (expressible in terms of x_indep), 
        f_x_dep: python function f describing x_dep = f(x_indep) 
        x_indep_0: x_indep guess value (starting point)
        """
        eig_vars = self.eig_vars
        if len(A_eqs_list)!=len(eig_vars):
            raise ValueError("Error: len(A_eqs_list)!=len(eig_vars")
        if solver=='muller':
            from Electromagnetics.Misc import solveMuller
            if len(eig_vars)==1:
                log_det_0 = 0.0
                if handle_overflow == True:
                    log_det_0 = self.getLogDeterminant(eigvar_vals_0, A_eqs_list)[0]
                def f(x):
                    eigvar_vals = [x]
                    if handle_overflow == False:
                        det_ = self.getDeterminant(eigvar_vals, A_eqs_list)[0]
                        return det_
                    else:
                        log_det = self.getLogDeterminant(eigvar_vals, A_eqs_list)[0]
                        det_ = np.exp(log_det - log_det_0)
                        if np.real(log_det - log_det_0) < -700:
                            det_ = np.exp(-700)
                        if np.real(log_det - log_det_0) > +700:
                            det_ = np.exp(+700)
                        return det_
                x_0 = eigvar_vals_0[0]
                x_1 = 1.1*x_0
                x_2 = 1.2*x_0
                if x_0==0.0:
                    x_1 = 1.0
                    x_2 = 2.0
                x_muller = solveMuller(f, x_0, x_1, x_2)
                if x_muller[3]==None:
                    print('number of tries: ', x_muller[2])
                    return [[x_muller[0]], [x_muller[1]]]
                else:
                    raise ValueError('Solver did not converge: ' + x_muller[3])
            else:
                raise NotImplementedError('Muller handles only one independant \
                    eigenvalue variable')
        else:
            methods_list = ['hybr', 'lm', 'broyden1', 'broyden2', 'anderson',
                            'linearmixing', 'diagbroyden', 'excitingmixing',
                            'krylov']
            if solver in methods_list==False:
                raise ValueError('Available methods: {0}'.format(methods_list))
            det_0 = self.getDeterminant(eigvar_vals_0, A_eqs_list)
            log_det_0 = self.getLogDeterminant(eigvar_vals_0, A_eqs_list)
            #print("xxxxxx det_0", det_0,  "  log_det_0: ", log_det_0)
            def f(x_ri):
                eigvar_vals = []
                for i in range(int(len(x_ri)/2)):
                    eigvar_vals.append(x_ri[2*i] + 1.0j*x_ri[2*i+1])
                if handle_overflow == False:
                    dets_ = self.getDeterminant(eigvar_vals, A_eqs_list)
                    dets_ri = []
                    for i in range(len(dets_)):
                        det_ = dets_[i]
                        if tol_relative:
                            det_ /= det_0[i]
                        if roots_prev!=None:
                            for j in range(len(roots_prev)):
                                for k in range(len(eigvar_vals)):
                                    if eigvar_vals[k] != roots_prev[j][k]:
                                        det_ /= (eigvar_vals[k] - roots_prev[j][k])
                        dets_ri.append(det_.real)
                        dets_ri.append(det_.imag)
                    return dets_ri
                else:
                    log_dets_ = self.getLogDeterminant(eigvar_vals, A_eqs_list)
                    dets_ri = []
                    for i in range(len(log_dets_)):
                        det_ = np.exp(log_dets_[i])
                        if np.real(log_dets_[i]) < -700:
                            det_ = np.exp(-700)
                        if np.real(log_dets_[i]) > +700:
                            det_ = np.exp(+700)
                        if tol_relative:
                            det_ = np.exp(log_dets_[i] - log_det_0[i])
                            if np.real(log_dets_[i] - log_det_0[i]) < -700:
                                det_ = np.exp(-700)
                            if np.real(log_dets_[i] - log_det_0[i]) > +700:
                                det_ = np.exp(+700)
                        ## TODO: depending on log_det_0 might return a math range error
                        if roots_prev!=None:
                            for j in range(len(roots_prev)):
                                for k in range(len(eigvar_vals)):
                                    det_ /= (eigvar_vals[k] - roots_prev[j][k])
                        dets_ri.append(det_.real)
                        dets_ri.append(det_.imag)
                    #print(x_ri, ' '*20, dets_ri)
                    return dets_ri
            from scipy import optimize
            #print('optimize.root started')
            x_0 = []
            for i in range(len(eigvar_vals_0)):
                x_0.append(eigvar_vals_0[i].real)
                x_0.append(eigvar_vals_0[i].imag)
            sol = optimize.root(f, x_0, method=solver, tol=tol, 
                options={'maxfev':maxfev, 'maxiter':maxiter, 'ftol': ftol})
            #print('tol: ', tol, '  maxfev: ', maxfev, '  maxiter:', maxiter)
            sol_x = sol.x
            sol_f = sol.fun
            x_ = []
            for i in range(int(len(sol_x)/2)):
                x_.append(sol_x[2*i] + 1j*sol_x[2*i+1])
            f_ = []
            for i in range(int(len(sol_f)/2)):
                f_.append(sol_f[2*i] + 1j*sol_f[2*i+1])
            if not getMats:
                return [x_, f_, sol]
            else:
                mats = self.getMatrices(x_, A_eqs_list)
                return [x_, f_, sol, mats]
                
        return
        
        
    def getVarsFields(self, x_vec):
        print('x_vec.shape:', x_vec.shape)
        self.wrapVariableIndices(x_vec)
        n_vars = len(self.vars)
        vars_fields = [None]*n_vars
        for i in range(n_vars):
            vars_fields[i] = self.getInverseFourierCoeffs(self.vars_ND_array[i])
        return vars_fields
        

    def getFourierCoeffs(self, f, x_0, x_1, N):
        """
        returns the multi-dimensional Fourier coeffs of a function that takes a numpy array as input
        harmonics from NF_0...NF_1-1
        x_0, x_1: limit points of the periodic domain
        N: maximum harmonic index 
        """
        n_dim = len(x_0)
        N_t = [0]*n_dim
        for i in range(n_dim):
            if self.symmetric_spec:
                N_t[i] = 2*N[i] + 1
            else:
                N_t[i] = 2*N[i]
        x = [None]*n_dim
        N_tot = 1
        for i in range(n_dim):
            x[i] = np.linspace(x_0[i], x_1[i], N_t[i], endpoint=False)
            N_tot *= N_t[i]
            
        x_mesh = np.meshgrid(*tuple(x))
        f_x = f(x_mesh)
        #print('f_x:\n', f_x)
        #print(np.min(f_x), np.max(f_x))
        #print('x_mesh:', x_mesh)
        fft_x = fftpack.fftshift(fftpack.fftn(f_x))/N_tot
        #print('fft_x:\n', np.abs(fft_x))
        return fft_x
        
        
    def getInverseFourierCoeffs(self, fft_x):
        """
        returns the function samples given its 1d Fourier coeffs
        """
        shape = fft_x.shape
        N_tot = 1
        for i in range(len(shape)):
            N_tot *= shape[i]
        
        f_x = fftpack.ifftn(fftpack.ifftshift(fft_x))*N_tot
        return f_x
        
    def getInverseFourierCoeffsNoShift(self, fft_x):
        """
        returns the function samples given its 1d Fourier coeffs
        """
        shape = fft_x.shape
        N_tot = 1
        for i in range(len(shape)):
            N_tot *= shape[i]
        
        f_x = fftpack.ifftn(fft_x)*N_tot
        return f_x

    def getFourierCoeffs_(self, f_x):
        """
        returns the function samples given its 1d Fourier coeffs
        """
        shape = f_x.shape
        N_tot = 1
        for i in range(len(shape)):
            N_tot *= shape[i]
        
        fft_x = fftpack.fftshift(fftpack.fftn(f_x))/N_tot
        return fft_x





        
