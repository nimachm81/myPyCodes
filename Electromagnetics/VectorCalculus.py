## VectorCalculus.py   Symbolic




__all__ = ["gradient_r", "divergence_r", "curl_r", "dotproduct", "crossproduct",
           "del_square_sc_r", "del_square_vec_r",
           "gradient_cy", "divergence_cy", "curl_cy",
           "getJacobianMatrix", "pdeChangeOfVariables"]


from sympy import Derivative, diff, zeros, abc, Symbol

from Electromagnetics.SymExprTree import *

def dotproduct(vec1, vec2):
    dot = vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]
    return dot

def crossproduct(vec1, vec2):
    cross = zeros(1,3)
    cross[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1]
    cross[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2]
    cross[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0]
    return cross

##rectangular
def gradient_r(scalar_field, x=abc.x, y=abc.y, z=abc.z):
    grad = zeros(1,3)
    grad[0] = Derivative(scalar_field, x)
    grad[1] = Derivative(scalar_field, y)
    grad[2] = Derivative(scalar_field, z)
    return grad

def divergence_r(vector_field, x=abc.x, y=abc.y, z=abc.z):
    div = Derivative(vector_field[0], x)
    div += Derivative(vector_field[1], y)
    div += Derivative(vector_field[2], z)
    return div

def curl_r(vector_field, x=abc.x, y=abc.y, z=abc.z):
    curl = zeros(1,3)
    curl[0] = Derivative(vector_field[2], y) - Derivative(vector_field[1], z)
    curl[1] = Derivative(vector_field[0], z) - Derivative(vector_field[2], x)
    curl[2] = Derivative(vector_field[1], x) - Derivative(vector_field[0], y)
    return curl

def del_square_sc_r(scalar_field, x=abc.x, y=abc.y, z=abc.z):
    grad = gradient_r(scalar_field, x, y, z)
    del_grad = divergence_r(grad, x, y, z)
    return del_grad

def del_square_vec_r(vector_field, x=abc.x, y=abc.y, z=abc.z):
    d2 = zeros(1,3)
    d2[0] = del_square_sc_r(vector_field[0], x, y, z)
    d2[1] = del_square_sc_r(vector_field[1], x, y, z)
    d2[2] = del_square_sc_r(vector_field[2], x, y, z)
    return d2

##cylindrical
def gradient_cy(scalar_field, rho=Symbol(r'\rho'), phi=Symbol(r'\phi'), z=abc.z):
    grad = zeros(1,3)
    grad[0] = Derivative(scalar_field, rho)
    grad[1] = Derivative(scalar_field, phi)/rho
    grad[2] = Derivative(scalar_field, z)
    return grad

def divergence_cy(vector_field, rho=Symbol(r'\rho'), phi=Symbol(r'\phi'), z=abc.z):
    div = Derivative(vector_field[0]*rho, rho)/rho
    div += Derivative(vector_field[1], phi)/rho
    div += Derivative(vector_field[2], z)
    return div

def curl_cy(vector_field, rho=Symbol(r'\rho'), phi=Symbol(r'\phi'), z=abc.z):
    curl = zeros(1,3)
    curl[0] = Derivative(vector_field[2], phi)/rho - Derivative(vector_field[1], z)
    curl[1] = Derivative(vector_field[0], z) - Derivative(vector_field[2], rho)
    curl[2] = Derivative(vector_field[1]*rho, rho)/rho - Derivative(vector_field[0], phi)/rho
    return curl





##------ PDE change of variables --------

def getJacobianMatrix(vars_old, vars_new, funcs):
    N = len(vars_old)
    assert len(vars_new)==N and len(funcs)==N
    
    Jac = zeros(N,N)
    for i in range(N):
        for j in range(N):
            Jac[i, j] = Derivative(funcs[i], vars_old[j])
            
    return Jac


def pdeChangeOfVariables(pde, vars_old, vars_new, funcs, funcs_inv):
    jac = getJacobianMatrix(vars_old, vars_new, funcs)
    
    expr_tree = symExpr_generate_tree(pde)
    symExp_replaceDerivative_walk_tree(expr_tree, vars_old, vars_new, jac, funcs_inv)
    return expr_tree[1]
    
    

def symExp_replaceDerivative_node(node, vars_old, vars_new, jac, funcs_inv):
    if node[1].func==Derivative:
        N_vars = len(vars_old)
        args = node[1].args
        assert len(args)>0
        f = args[0]
        args_ind_to_remove = []
        var_old_exist = False
        for i in range(1, len(args)):
            for j in range(N_vars):
                if args[i]==vars_old[j]:
                    var_old_exist = True
                    f_ = 0
                    for k in range(N_vars):
                        df_dnew = Derivative(f, vars_new[k])
                        f_ += jac[k, j].doit().subs([(vars_old[i_], funcs_inv[i_]) \
                            for i_ in range(N_vars)])*df_dnew
                    f = f_
                    args_ind_to_remove.append(i)
        args = [args[i] for i in range(len(args)) if not i in args_ind_to_remove]
        args[0] = f
        #print(args)
        if len(args)>1:
            node[1] = Derivative(*tuple(args))
        else:
            node[1] = args[0]
        symExpr_update_tree(node)
        if var_old_exist:
            return True
        else:
            return False
    else:
        return False


def symExp_replaceDerivative_walk_tree(node, vars_old, vars_new, jac, funcs_inv):
    nodes_list = [node]
    ind_next = 0
    while True:
        replaced = symExp_replaceDerivative_node(nodes_list[ind_next], vars_old, vars_new, jac, funcs_inv)
        if replaced==False:
            node_next = nodes_list[ind_next]
            for arg in node_next[3]:
                if len(arg)>0:
                    nodes_list.append(arg)
        ind_next += 1
        if ind_next>=len(nodes_list):
            break
    return 












