## FourierBloch.py
## Fourier Bloch analysis


from sympy import zeros, abc, Matrix, eye, Symbol, Sum, exp, I, oo, Function,\
             expand, Mul, Add, Tuple, Integer, lambdify

from Electromagnetics.SymExprTree import *

__all__ = ["d1_putSums", "d1_applyConvolutions", "d1_applyOrthogonalities",
           "d1_setup_symmetric_spectrum",
           "d1_orthogonalToNumpyMatrix", "d1_getMatrix", "d1_getDeterminant", 
           "d1_solveDeterminant",
           "d1_getFourierCoeffs", "d1_getInverseFourierCoeffs"]



###--------  Fourier Bloch Analysis

#TODO check if expression is linear in a given argument

### TODO: carefull with simplify when using this package as the sum of more 
### than 2 series in sympy returns the wrong response (issue #8596)
### print simplify(Sum(x,(n,-oo,oo))+Sum(y,(n,-oo,oo))+Sum(z,(n,-oo,oo))+Sum(t,(n,-oo,oo)))
### ---> Sum(x + y, (n, -oo, oo))


def d1_putSums(expr, vars, t=Symbol('t'), harmonic=exp(I*Symbol('n')*Symbol('\\omega_0')*Symbol('t')), 
        n_str='n', mark='\\tilde'):
    """Replaces the vars with their periodic Fourier sums
    """
    #TODO should support time dependant functions (instead of Symbols) 
    # otherwise the time derivatives have to be performed after 
    # calling this function not before
    n = Symbol(n_str)
    expr_new = expr.subs([(v, Sum(Symbol(mark+'{'+v.name+'}')(n)
        *harmonic,(n,-oo,+oo))) for v in vars])
    for i in range(len(vars)):
        vars[i] = Function(mark+'{'+vars[i].name+'}')
    return expr_new


def d1_isFourierSeries(expr, t=Symbol('t'), harmonic=exp(I*Symbol('n')*Symbol('\\omega_0')*Symbol('t')), 
        n_str='n', mark='\\tilde'):
    n = Symbol(n_str)
    #print('expr: ', expr)
    if expr.func == Sum:
        args = expr.args
        if args[0].has(harmonic) and args[1]==(n, -oo, +oo):
            return True
    return False


def d1_getCoefficientPartOfSeries(expr, t=Symbol('t'), harmonic=exp(I*Symbol('n')*Symbol('\\omega_0')*Symbol('t')), 
        n_str='n', mark='\\tilde'):
    n = Symbol(n_str)
    #print('expr: ', expr)
    if expr.func == Sum:
        args = expr.args
        if args[0].has(harmonic) and args[1]==(n, -oo, +oo):
            expr_coeff = args[0]/harmonic
            #if expr_coeff.has(t):
            #    raise ValueError('expression is not a proper Fourier series!')
            #else:
            #    return expr_coeff
            return expr_coeff
    raise ValueError('expression is not a Fourier series!')
    return None


def d1_MultiplySeries_Convulve(expr1, expr2, t=Symbol('t'), harmonic=exp(I*Symbol('n')*Symbol('\\omega_0')*Symbol('t')), 
        n_str='n', mark='\\tilde'):
    n = Symbol(n_str)
    expr1_coeff = d1_getCoefficientPartOfSeries(expr1, t, harmonic, n_str, mark)
    expr2_coeff = d1_getCoefficientPartOfSeries(expr2, t, harmonic, n_str, mark)
    i = 1;
    n_dummy = Symbol(n_str + '_' + str(i))
    while expr1_coeff.has(n_dummy) or expr2_coeff.has(n_dummy):
        i += 1
        n_dummy = Symbol(n_str + '_' + str(i))
    expr1_coeff = expr1_coeff.subs(n, n_dummy)
    expr2_coeff = expr2_coeff.subs(n, n - n_dummy)
    expr_new_coeff = Sum(expr1_coeff*expr2_coeff, (n_dummy,-oo,oo))
    expr_new = Sum(expr_new_coeff*harmonic, (n, -oo, oo))
    return expr_new


def d1_applyConvolutionToNode(node, t=Symbol('t'), harmonic=exp(I*Symbol('n')*Symbol('\\omega_0')*Symbol('t')), 
        n_str='n', mark='\\tilde'):
    """Apply convolution to node if it contains multiplication of series
    """
    if node[2] == Mul:
        for i, arg_i in  enumerate(node[3]):
            if d1_isFourierSeries(arg_i[1], t, harmonic, n_str, mark):
                for j, arg_j in enumerate(node[3]):
                    if i<j:
                        if d1_isFourierSeries(arg_j[1], t, harmonic, n_str, mark):
                            expr_conv = d1_MultiplySeries_Convulve(arg_i[1], arg_j[1], t, harmonic, n_str, mark)
                            node[3][i][1] = expr_conv
                            del node[3][j]
                            node[1] = node[2](*tuple(node[3][k][1] for k in range(len(node[3]))))
                            symExpr_update_tree(node)
                            return True
    return False



def d1_applyConvolutions_walk_tree(node, t=Symbol('t'), harmonic=exp(I*Symbol('n')*Symbol('\\omega_0')*Symbol('t')), 
        n_str='n', mark='\\tilde'):
    d1_applyConvolutionToNode(node, t, harmonic, n_str, mark)
    for arg_i in node[3]:
        if len(arg_i)>0:
            d1_applyConvolutions_walk_tree(arg_i, t, harmonic, n_str, mark)
    return



def d1_applyConvolutions(expr, t=Symbol('t'), harmonic=exp(I*Symbol('n')*Symbol('\\omega_0')*Symbol('t')), 
        n_str='n', mark='\\tilde'):
    """Apply convolutions
    the equations should be expanded first in order to properly recognize
    the multiplication of 2 series
    """
    expr_expanded = expand(expr)
    expr_tree = symExpr_generate_tree(expr_expanded)
    d1_applyConvolutions_walk_tree(expr_tree, t, harmonic, n_str, mark)
    return expr_tree[1]

def d1_aSumb_to_Sumab(node, t=Symbol('t'), harmonic=exp(I*Symbol('n')*Symbol('\\omega_0')*Symbol('t')), 
        n_str='n', mark='\\tilde'):
    """ If the expression at a given node is a*Sum(b*harmonic)
    it will be modified to Sum(a*b*harmonic)
    """
    """###implementation using wild matches - not used for safety 
    ####(the outcome is not exactly predictable)
    AA = Wild('AA')
    BB = Wild('BB')
    n = Symbol(n_str)
    match_pattern = node[1].match(AA*Sum(BB*harmonic, (n,-oo,oo)))
    if match_pattern != None:
        node[1] = Sum(AA*BB*harmonic, (n,-oo,oo))
        symExpr_update_tree(node)
    else:
        return
    """
    n = Symbol(n_str)
    if node[2]==Mul:
        for i, arg_i in enumerate(node[3]):
            if d1_isFourierSeries(arg_i[1], t, harmonic, n_str, mark):
                coeff_in = d1_getCoefficientPartOfSeries(arg_i[1], t, harmonic, n_str, mark)
                mul_args_except_sum = [node[3][k][1] for k in range(len(node[3]))]
                del mul_args_except_sum[i]
                coeff_out = Mul(*tuple(mul_args_except_sum))
                node[1] = Sum(coeff_out*coeff_in*harmonic, (n,-oo,oo))
                symExpr_update_tree(node)
                return True
    return False

def d1_aSumb_to_Sumab_walk_tree(node, t=Symbol('t'), harmonic=exp(I*Symbol('n')*Symbol('\\omega_0')*Symbol('t')), 
        n_str='n', mark='\\tilde'):
    d1_aSumb_to_Sumab(node, t, harmonic, n_str, mark)
    for arg_i in node[3]:
        if len(arg_i)>0:
            d1_aSumb_to_Sumab_walk_tree(arg_i, t, harmonic, n_str, mark)
    return

def d1_SumaPlusSumb_to_SumaPlusb(node, t=Symbol('t'), harmonic=exp(I*Symbol('n')*Symbol('\\omega_0')*Symbol('t')), 
        n_str='n', mark='\\tilde'):
    """ If the expression at a given node is Sum(a*harmonic) + Sum(b*harmonic)
    it will be modified to Sum((a+b)*harmonic)
    """
    n = Symbol(n_str)
    if node[2]==Add:
        for i, arg_i in enumerate(node[3]):
            if d1_isFourierSeries(arg_i[1], t, harmonic, n_str, mark):
                for j, arg_j in enumerate(node[3]):
                    if i<j:
                        if d1_isFourierSeries(arg_j[1], t, harmonic, n_str, mark):
                            coeff_i = d1_getCoefficientPartOfSeries(arg_i[1], t, harmonic, n_str, mark)
                            coeff_j = d1_getCoefficientPartOfSeries(arg_j[1], t, harmonic, n_str, mark)
                            if len(node[3])==2:
                                node[1] = Sum((coeff_i+coeff_j)*harmonic, (n,-oo,oo))
                                symExpr_update_tree(node)
                                return True
                            else:
                                node_1_arg = [node[3][k][1] for k in range(len(node[3]))]
                                del node_1_arg[j]
                                del node_1_arg[i]
                                sum_i_j = Sum((coeff_i+coeff_j)*harmonic, (n,-oo,oo))
                                node_1_arg.append(sum_i_j)
                                node[1] = node[2](*tuple(node_1_arg))
                                symExpr_update_tree(node)
                                d1_SumaPlusSumb_to_SumaPlusb(node, t, harmonic, n_str, mark)
                                return True
    return False


def d1_SumaPlusSumb_to_SumaPlusb_walk_tree(node, t=Symbol('t'), harmonic=exp(I*Symbol('n')*Symbol('\\omega_0')*Symbol('t')), 
        n_str='n', mark='\\tilde'):
    d1_SumaPlusSumb_to_SumaPlusb(node, t, harmonic, n_str, mark)
    for arg_i in node[3]:
        if len(arg_i)>0:
            d1_SumaPlusSumb_to_SumaPlusb_walk_tree(arg_i, t, harmonic, n_str, mark)
    return


def d1_addSeriesManually(expr, t=Symbol('t'), harmonic=exp(I*Symbol('n')*Symbol('\\omega_0')*Symbol('t')), 
        n_str='n', mark='\\tilde'):
    expr_expanded = expand(expr)
    expr_tree = symExpr_generate_tree(expr_expanded)
    d1_aSumb_to_Sumab_walk_tree(expr_tree, t, harmonic, n_str, mark)
    d1_SumaPlusSumb_to_SumaPlusb_walk_tree(expr_tree, t, harmonic, n_str, mark)
    return expr_tree[1]

def d1_applyOrthogonalities(expr, t=Symbol('t'), harmonic=exp(I*Symbol('n')*Symbol('\\omega_0')*Symbol('t')), 
        n_str='n', mark='\\tilde'):
    ##from sympy import Wild
    ##AA = Wild('AA')
    ##BB = Wild('BB')
    ##n = Symbol(n_str)
    ## TODO : redo after this sympy bug is fixed
    ##match_pattern = (1*simplify(expr)).match(AA*Sum(BB*harmonic, (n,-oo,oo)))
    ##
    ## this code should work, however simplify(addition of series) returns the wrong
    ## result if there are more than 2 series.. (BUG in sympy)
    ## the series are added manually in the following
    ##
    n = Symbol(n_str)
    expr = d1_addSeriesManually(expr, t, harmonic, n_str, mark)
    if expr.func == Sum:
        return d1_getCoefficientPartOfSeries(expr, t, harmonic, n_str, mark)
    elif expr == 0:
        return 0
    else:
        from IPython.display import display, Math, Latex
        display(Math('expr = ' + latex(expr)))
        print(expr)
        raise ValueError("Orthogonality was not assured..")
        return None


def aSumb_to_Sumab_1d(node, n_str='n'):
    """ If the expression at a given node is a*Sum(b, limits)
    it will be modified to Sum(a*b, limits)
    n_srt: dummy indices of Sums are n_1, n_2, ...
           n_str is the index of some variables f(n) (Function('f')(Symbol(n_str)))
    """
    n = Symbol(n_str)
    if node[2]==Mul:
        for i, arg_i in enumerate(node[3]):
            if arg_i[2]==Sum:
                coeff_in = arg_i[3][0][1]
                sum_limits_arg = arg_i[3][1][1]
                if sum_limits_arg.func!=Tuple:
                    raise NotImplementedError('It is assumed Sum has a Tuple as the second argument.')
                dummy = sum_limits_arg.args[0].name    # 'n_1' or 'n_2' ...
                dummy_nodigits = ''.join([k for k in dummy if not (k.isdigit() or k=='_')])
                if dummy_nodigits!=n_str:
                    continue
                mul_args_except_sum = [node[3][k][1] for k in range(len(node[3]))]
                del mul_args_except_sum[i]
                coeff_out = Mul(*tuple(mul_args_except_sum))
                if coeff_out.has(sum_limits_arg.args[0]):
                    raise NotImplementedError('It is assumed the coefficient outside \
                        the sum does not include the dummy index of the sum')
                node[1] = Sum(coeff_out*coeff_in, sum_limits_arg)
                symExpr_update_tree(node)
                return True
    return False

def aSumb_to_Sumab_1d_walk_tree(node, n_str='n'):
    aSumb_to_Sumab_1d(node, n_str)
    for arg_i in node[3]:
        if len(arg_i)>0:
            aSumb_to_Sumab_1d_walk_tree(arg_i, n_str)
    return


def find_node_with_function(node, f):
    """ finds the forst node whose node[2] component is the given function f and returns
        it starts from the given node and walks down the tree
    """
    if node[2]==f:
        return node
    else:
        for arg in node[3]:
            subnode_f = find_node_with_function(arg, f)
            if subnode_f!=False:
                return subnode_f
    return False

def find_node_with_symbol(node, x):
    """ finds the forst node whose node[1] component is the given symbol x and returns
        it starts from the given node and walks down the tree
    """
    if node[1]==x:
        return node
    else:
        for arg in node[3]:
            subnode_x = find_node_with_symbol(arg, x)
            if subnode_x!=False:
                return subnode_x
    return False


def d1_find_mult_coeff_of_func(node, var):
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
            return d1_find_mult_coeff_of_func(node_var, var)
        return False


def d1_getcoeff_setzero(node, var):
    """ var is a Function. Its argumet may contain Symbol(n_str) or Symbol(n_str+'_1')..
    """
    node_var = find_node_with_function(node, var)
    if node_var!=False:
         var_coeff_and_arg = d1_find_mult_coeff_of_func(node_var, var)
         node_var[1] = Integer(0)
         symExpr_update_tree(node_var)
         return var_coeff_and_arg    # [var_coeff, var_arg]
    return False

def d1_remove_sums_of_zero_arg(node):
    """
    removes Sum functions whose argumet is zero
    """
    if node[2]==Sum:
        if node[1].args[0]==0:
            node[1] = Integer(0)    
            symExpr_update_tree(node)
    for arg in node[3]:
        if len(arg)>0:
            d1_remove_sums_of_zero_arg(arg)
    return


##----------------------  Numerical part (Numpy, Scipy)

import numpy as np
from scipy.sparse import csr_matrix, coo_matrix
from scipy import linalg, fftpack
from cmath import *
import cmath

d1_symmetric_spec = False       # symmetric spectrum -N...N

def d1_setup_symmetric_spectrum(sym):
    global d1_symmetric_spec
    if sym==True:
        d1_symmetric_spec = True
    elif sym==False:
        d1_symmetric_spec = False
    else:
        raise ValueError('takes boolean input')
    return

def d1_orthogonalToNumpyMatrix(expr_list, n_str, N, vars, pars, pars_vecs, eig_vars):
    """It takes the orthogonal Fourier expression and generates a numpy matrix
        to solve.
       expr_list: sympy expression list (the given system of equations: expr_list[i]=0)
       N: number of Fourier harmonics (-N...N-1)
       vars: list containing unknown fourier variables
       pars: list containing known Fourier parameters
       pars_vecs: list containig numpy arrays holding the coefficients of each variable
              in pars
       eig_vars: list containing unknown scalars (usually the eigenvalue)
    """
    NF_0 = -N
    NF_1 = N
    if d1_symmetric_spec==True:
        NF_1 = N + 1
    N_mp = NF_1 - NF_0  #2*N  # Fourier indexes _N...N-1: 2*N (for speed)
    def var_index(v, i):
        return v*N_mp+i+N
    n_total = N_mp*len(vars)
    if len(expr_list)!=len(vars):
        raise ValueError('number of equations should be equal to the number of \
            variables')
    ### matrices
    A_lin_sp = csr_matrix( (n_total,n_total), dtype=complex)
    A_eig_sp_list = []
    A_convs = np.zeros((n_total, n_total), dtype=complex)
    b_rhs = np.zeros((n_total, 1))
    for eq_ind, expr in enumerate(expr_list):
        expr_expanded = expand(expr)
        expr_tree = symExpr_generate_tree(expr_expanded)
        aSumb_to_Sumab_1d_walk_tree(expr_tree, n_str)
        for ind_v, var in enumerate(vars):
            while True:
                coeff_arg = d1_getcoeff_setzero(expr_tree, var)
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
                            a_lin_sp_coo = d1_orth_ToNumpyMatrix_const_coeff(coeff,
                                   arg, n_str, N, n_total, ind_v, eq_ind*N_mp)
                            A_lin_sp = A_lin_sp + a_lin_sp_coo.tocsr()
                        else:
                            A_eig, a_eig_sp_coo = d1_orth_toNumpyMatrix_eig_coeff(coeff,
                                  arg, n_str, N, n_total, ind_v, eq_ind*N_mp, eig_vars)
                            a_eig_sp_csr = a_eig_sp_coo.tocsr()
                            A_eig_sp_list.append([A_eig, a_eig_sp_csr])
                    elif n_par==1:
                        A_dense = d1_orth_ToNumpyMatrix_convolutions(coeff, arg, n_str, N, n_total, 
                         ind_v, eq_ind*N_mp, pars, pars_vecs)
                        A_convs += A_dense
                    else:
                        raise NotImplementedError('more than one periodic coefficient..')
        d1_remove_sums_of_zero_arg(expr_tree)
        if expr_tree[1]!=0:
            # get the rhs
            from IPython.display import display, Math
            print(expr_tree[1])
            display(Math(latex(expr_tree[1])))
            raise NotImplementedError('non-zero rhs not implemented..')
    return [A_lin_sp, A_convs, A_eig_sp_list, b_rhs]



def d1_orth_ToNumpyMatrix_const_coeff(coeff, arg, n_str, N, n_total, 
    ind_v, eq_ind_start):
    """
    it takes a coefficient and argument (coeff*var(arg)) depennding only on 
    n_str and returns a sparse numpy matrix
    coeff*vars[var_ind](arg)
    coeff: coefficient expression
    arg: argument expression
    n_str: fourier index (default n)
    N: number of Fourier harmonics to keep
    n_total: total number of variables
    ind_v: index of the variable inside vars
    eq_ind_start: index of the first equation (row indices)
    """
    NF_0 = -N
    NF_1 = N
    if d1_symmetric_spec==True:
        NF_1 = N + 1
    N_mp = NF_1 - NF_0  #2*N  # Fourier indexes -N...N-1: 2*N (for speed)
    def var_index(i):
        if i>=NF_0 and i<NF_1:
            return ind_v*N_mp+i+N
        else:
            return -1
    var_index_np = np.frompyfunc(var_index,1,1)
    def eq_index(i):
        return eq_ind_start+i+N
    eq_index_np = np.frompyfunc(eq_index,1,1) 
    n = Symbol(n_str)
    f_coeff = lambdify(n, coeff, "numpy")
    f_arg = lambdify(n, arg[0], "numpy")
    n_arr = np.array(range(NF_0,NF_1))
    col = var_index_np(f_arg(n_arr))
    row = eq_index_np(n_arr)
    data = None
    if coeff.has(n):
        data = f_coeff(n_arr)
    else:
        data = np.array([coeff.evalf()]*len(n_arr))
    ind_to_remove = []
    for i in range(len(col)-1, -1, -1):
        if col[i]<0:
            del row[i]
            del col[i]
            del data[i]
    mat_sp = coo_matrix((data, (row,col)), shape=(n_total,n_total), dtype=complex)
    return mat_sp


def d1_split_expr(expr, symbs_1, symbs_2):
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


def d1_orth_toNumpyMatrix_eig_coeff(coeff, arg, n_str, N, n_total, ind_v, 
    eq_ind_start, eig_vars):
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
    expr_split = d1_split_expr(coeff.expand(), eig_vars, [Symbol(n_str)])
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
    mat_sp = d1_orth_ToNumpyMatrix_const_coeff(B_n, arg, n_str, N, n_total, 
        ind_v, eq_ind_start)
    return [A_eigs, mat_sp] ## A_eigs*mat_sp


def d1_get_func_arg(expr, func):
    """
    finds the first instance of the func and returns its argument
    """
    expr_tree = symExpr_generate_tree(expr)
    node_func = find_node_with_function(expr_tree, func)
    if node_func!=False:
        return  node_func[1].args
    return False


def d1_orth_ToNumpyMatrix_convolutions(coeff, arg, n_str, N, n_total, 
    ind_v, eq_ind_start, pars, pars_vecs):
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
    NF_0 = -N
    NF_1 = N
    if d1_symmetric_spec==True:
        NF_1 = N + 1
    N_mp = NF_1 - NF_0  #2*N  # Fourier indexes _N...N-1: 2*N (for speed)
    def var_index(i):
        return ind_v*N_mp+i+N
    def eq_index(i):
        return eq_ind_start+i+N
    eq_index_np = np.frompyfunc(eq_index,1,1) 
    # get argumet and coeff of the periodic parameter
    p_coeff = None
    p_arg = None
    p_vec = None
    for i, p in enumerate(pars):
        if coeff.has(p):
            p_arg_list = d1_get_func_arg(coeff, p)
            if p_arg_list==False or len(p_arg_list)==0:
                print("p_arg_list: ", p_arg_list)
                raise ValueError('par should be a periodic function')
            p_arg = p_arg_list[0]
            p_coeff = symExp_replaceFunction(coeff, p, Symbol('AAAZZ'))\
                 .coeff(Symbol('AAAZZ'))
            p_vec = pars_vecs[i]
            break
    for p in pars:
        if p_coeff.has(p):
            raise NotImplementedError('Only one periodic parameter can be handled')
    # constructing functions
    n = Symbol(n_str)
    n_1 = Symbol(n_str+'_1')
    f_p_coeff = lambdify(n, p_coeff)
    f_p_arg = lambdify((n, n_1), p_arg)
    f_arg = lambdify((n, n_1), arg[0])
    #print('arg = ', arg)
    # constructing the dense coefficien2t matrix
    A_mat = np.zeros((n_total, n_total), dtype=complex)
    n_arr = range(NF_0,NF_1)
    for i in n_arr:
        row_ind = eq_index(i)
        __c = f_p_coeff(i)
        for j in n_arr:
            _ind_v_ = f_arg(i, j)
            _ind_p_ = f_p_arg(i, j)
            if NF_0<=_ind_v_<NF_1 and NF_0<=_ind_p_<NF_1:
                col_ind = var_index(_ind_v_)
                A_mat[row_ind, col_ind] += __c*p_vec[_ind_p_+N]
    return A_mat


def d1_getMatrix(eig_vars, eigvar_vals, A_eqs):
    A_lin_sp, A_convs, A_eig_sp_list, b_rhs = A_eqs
    A_tot = A_convs.copy()
    if A_lin_sp.getnnz()>0:
        A_tot += A_lin_sp
    for i in range(len( A_eig_sp_list)):
        C = complex(A_eig_sp_list[i][0].evalf(subs=dict(zip(eig_vars, eigvar_vals))))
        A_tot +=  C*A_eig_sp_list[i][1]
    return A_tot
    

def d1_getDeterminant(eig_vars, eigvar_vals, A_eqs_list):
    """
    find the determinant of the eigenvalue system where the eigenvalues are
    given in eigvar_vals
    """
    dets = []
    for A_eqs in A_eqs_list:
        A_lin_sp, A_convs, A_eig_sp_list, b_rhs = A_eqs
        A_tot = A_convs.copy()
        if A_lin_sp.getnnz()>0:
            A_tot += A_lin_sp
        for i in range(len( A_eig_sp_list)):
            C = complex(A_eig_sp_list[i][0].evalf(subs=dict(zip(eig_vars, eigvar_vals))))
            A_tot +=  C*A_eig_sp_list[i][1]
        dets.append(linalg.det(A_tot))
    return dets
    

def d1_getLogDeterminant(eig_vars, eigvar_vals, A_eqs_list):
    """
    find the logarithm of the determinant of the eigenvalue system where 
    the eigenvalues are given in eigvar_vals
    """
    dets = []
    for A_eqs in A_eqs_list:
        A_lin_sp, A_convs, A_eig_sp_list, b_rhs = A_eqs
        A_tot = A_convs.copy()
        if A_lin_sp.getnnz()>0:
            A_tot += A_lin_sp
        for i in range(len( A_eig_sp_list)):
            C = complex(A_eig_sp_list[i][0].evalf(subs=dict(zip(eig_vars, eigvar_vals))))
            A_tot +=  C*A_eig_sp_list[i][1]
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

    
def d1_solveDeterminant(eig_vars, eigvar_vals_0, A_eqs_list,
        solver='lm', handle_overflow=False, roots_prev=None):
    """
    solves the determinant for eigenvalues
    x_indep: independant sympy variables (list), 
    x_dep: dependant sympy variables (expressible in terms of x_indep), 
    f_x_dep: python function f describing x_dep = f(x_indep) 
    x_indep_0: x_indep guess value (starting point)
    """
    if len(A_eqs_list)!=len(eig_vars):
        raise ValueError("Error: len(A_eqs_list)!=len(eig_vars")
    if solver=='muller':
        from Misc import solveMuller
        if len(eig_vars)==1:
            log_det_0 = 0.0
            if handle_overflow == True:
                log_det_0 = d1_getLogDeterminant(eig_vars, eigvar_vals_0, 
                            A_eqs_list)[0]
            def f(x):
                eigvar_vals = [x]
                if handle_overflow == False:
                    det_ = d1_getDeterminant(eig_vars, eigvar_vals, 
                           A_eqs_list)[0]
                    return det_
                else:
                    log_det = d1_getLogDeterminant(eig_vars, eigvar_vals, 
                            A_eqs_list)[0]
                    det_ = cmath.exp(log_det - log_det_0)
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
        log_det_0 = d1_getLogDeterminant(eig_vars, eigvar_vals_0, A_eqs_list)
        det_0 = d1_getDeterminant(eig_vars, eigvar_vals_0, A_eqs_list)
        #print("xxxxxx det_0", det_0,  "  log_det_0: ", log_det_0)
        def f(x_ri):
            eigvar_vals = []
            for i in range(int(len(x_ri)/2)):
                eigvar_vals.append(x_ri[2*i] + 1.0j*x_ri[2*i+1])
            if handle_overflow == False:
                dets_ = d1_getDeterminant(eig_vars, eigvar_vals, 
                       A_eqs_list)
                dets_ri = []
                for i in range(len(dets_)):
                    det_ = dets_[i]/det_0[i]
                    if roots_prev!=None:
                        for j in range(len(roots_prev)):
                            for k in range(len(eigvar_vals)):
                                det_ /= (eigvar_vals[k] - roots_prev[j][k])
                    dets_ri.append(det_.real)
                    dets_ri.append(det_.imag)
                return dets_ri
            else:
                log_dets_ = d1_getLogDeterminant(eig_vars, eigvar_vals, 
                        A_eqs_list)
                dets_ri = []
                for i in range(len(log_dets_)):
                    det_ = cmath.exp(log_dets_[i] - log_det_0[i])
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
        tol, maxfev, maxiter = 1.0e-10, 1000, 1000
        sol = optimize.root(f, x_0, method=solver, tol=tol, 
            options={'maxfev':maxfev, 'maxiter':maxiter, 'ftol': 1.0e-5})
        #print('tol: ', tol, '  maxfev: ', maxfev, '  maxiter:', maxiter)
        sol_x = sol.x
        sol_f = sol.fun
        x_ = []
        for i in range(int(len(sol_x)/2)):
            x_.append(sol_x[2*i] + 1j*sol_x[2*i+1])
        f_ = []
        for i in range(int(len(sol_f)/2)):
            f_.append(sol_f[2*i] + 1j*sol_f[2*i+1])
        return [x_, f_, sol]
        
    return

'''
def d1_getLogDeterminant_plot__(A_lin_sp, A_convs, A_eig_sp_list,
        x_indep, x_dep, f_x_dep, x_indep_vals):
    """
    helper function for plotting the determinant
    """
    n_x_dep = len(x_dep)
    x_dep_vals = []
    for i in range(n_x_dep):
        x_dep_vals.append(f_x_dep[i](*tuple(x_indep_vals)))
    eig_vars = x_indep + x_dep
    eigvar_vals = x_indep_vals + x_dep_vals
    log_det = d1_getLogDeterminant(eig_vars, eigvar_vals, A_lin_sp, 
        A_convs, A_eig_sp_list)
    return log_det
    
def d1_getLogDeterminant_plot(A_lin_sp, A_convs, A_eig_sp_list,
        x_indep, x_dep, f_x_dep, x_indep_vals, ind_x_indep,
        x_r_limits, x_i_limits, n_x_r=100, n_x_i=100):
    """
    creating a 2D pcolor at the given range for the independent variable
    x_indep[ind_x_indep]
    the values for other independant variable are given at x_indep_vals,
    x_indep_vals[ind_x_indep] is unused
    """
    x_r = np.linspace(x_r_limits[0], x_r_limits[1], n_x_r)   
    x_i = np.linspace(x_i_limits[0], x_i_limits[1], n_x_i)
    x_r, x_i = np.meshgrid(x_r, x_i)
    z = np.zeros(x_r.shape, dtype=complex)
    for i in range(z.shape[0]):
        for j in range(z.shape[1]):
            x_indep_vals[ind_x_indep] = x_r[i,j] + 1.0j*x_i[i,j]
            z[i,j] = d1_getLogDeterminant_plot__(A_lin_sp, A_convs, A_eig_sp_list,
                x_indep, x_dep, f_x_dep, x_indep_vals)
    return [x_r, x_i, z]
'''

def d1_getFourierCoeffs(f, x_0, x_1, N):
    """
    returns the 1d Fourier coeffs of a function that takes a numpy array as input
    harmonics from NF_0...NF_1-1
    x_0, x_1: limit points of the periodic domain
    N: maximum harmonic index 
    """
    NF_0 = -N
    NF_1 = N
    if d1_symmetric_spec==True:
        NF_1 = N + 1
    N_t = NF_1 - NF_0
    x = np.linspace(x_0, x_1, N_t, endpoint=False)
    f_x = f(x)
    fft_x = fftpack.fftshift(fftpack.fft(f_x))/N_t
    return fft_x
    

def d1_getInverseFourierCoeffs(fft_x):
    """
    returns the function samples given its 1d Fourier coeffs
    """
    f_x = fftpack.ifft(fftpack.ifftshift(fft_x))*len(fft_x)
    return f_x   
    
    

