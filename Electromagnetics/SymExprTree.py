## SymExprTree.py
## Manipulating sympy expressions through expression tree
## utilising expr.func expr.args srepr

__all__ = ["symExpr_generate_tree", "symExpr_update_tree", "symExp_replaceSymbol", "symExp_replaceFunction",
            "symExp_replaceFunction_withArgs", "symExp_CompareFunctionArgs",
            "symExp_FindReplaceFunctionArgs"]


from sympy import Function, srepr, Symbol, latex
from IPython.display import Math, display

###---------  sympy expression tree manipulation

def symExpr_get_tree_root(expr):
    """ node[0]: parent node
        node[1]: node expression
        node[2]: node function
        node[3]: branches (list, each element is a branch node)
        node[4]: index of the node expression inside the parent's branch list
    """
    root = [None, expr, expr.func, list(expr.args), -1]
    return root

def symExpr_add_branches(node):
    for i, arg_i in enumerate(node[3]):
        node[3][i] = [node, node[3][i], arg_i.func, list(arg_i.args), i]
        symExpr_add_branches(node[3][i])
    return

def symExpr_generate_tree(expr):
    root = symExpr_get_tree_root(expr)
    symExpr_add_branches(root)
    return root


def symExpr_update_parents(node):
    parent = node[0]
    if parent != None:
        parent[1] = parent[2](*tuple(parent[3][k][1] for k in range(len(parent[3]))))
        symExpr_update_parents(parent)
    return
    
def symExpr_update_tree(node):
    node[2] = node[1].func
    node[3] = list(node[1].args)
    symExpr_add_branches(node)
    symExpr_update_parents(node)
    return

##------------ replace symbol (useful when symbol is inside a derivative operator)

def symExp_replaceSymbol_node(node, var, var_new):
    #display(Math(latex(node[1])))
    if node[1]==var:
        node[1] = var_new
        #print("var replaced!")
        symExpr_update_tree(node)
        return True
    return False

def symExp_replaceSymbol_walk_tree(node, var, var_new):
    nodes_list = [node]
    ind_next = 0
    while True:
        replaced = symExp_replaceSymbol_node(nodes_list[ind_next], var, var_new)
        node_next = nodes_list[ind_next]
        if not replaced:
            for arg in node_next[3]:
                if len(arg)>0:
                    nodes_list.append(arg)
        ind_next += 1
        if ind_next>=len(nodes_list):
            break
    return 


def symExp_replaceSymbol(expr, var, var_expr_new):
    """It searches the expr for symbol var and replaces it with the new expression var_new
    """
    expr_tree = symExpr_generate_tree(expr)
    symExp_replaceSymbol_walk_tree(expr_tree, var, var_expr_new)
    return expr_tree[1]
    


##------------ replace function (ignore its arguments)

def symExp_replaceFunction_walk_tree_INACTIVE(node, f, expr_new):
    if node[1].func==f:
        node[1] = expr_new
        symExpr_update_tree(node)
        #print('found sample: ' + srepr(f))
    for arg in node[3]:
        if len(arg)>0:
            symExp_replaceFunction_walk_tree(arg, f, expr_new)
    return 

def symExp_replaceFunction_node(node, f, expr_new):
    if node[1].func==f:
        node[1] = expr_new
        symExpr_update_tree(node)
    return 

def symExp_replaceFunction_walk_tree(node, f, expr_new):
    nodes_list = [node]
    ind_next = 0
    while True:
        symExp_replaceFunction_node(nodes_list[ind_next], f, expr_new)
        node_next = nodes_list[ind_next]
        for arg in node_next[3]:
            if len(arg)>0:
                nodes_list.append(arg)
        ind_next += 1
        if ind_next>=len(nodes_list):
            break
    return 


def symExp_replaceFunction(expr, f, f_expr_new):
    """It searches the expr for function f and replaces it with the new expression expr_new
        (it ignores the argument of the f i.e. f(x) is relaced with expr_new no matter
        what is (x))
    """
    expr_tree = symExpr_generate_tree(expr)
    symExp_replaceFunction_walk_tree(expr_tree, f, f_expr_new)
    return expr_tree[1]
    
##-------------------- Replace function with arguments

def symExp_replaceFunction_withArgs_node(node, f, f_args, expr_new):
    if node[1].func==f:
        args = node[1].args
        assert len(args)==len(f_args)
        node_new = expr_new
        for i in range(len(args)):
            node_new = node_new.subs(f_args[i], args[i])
            
        node[1] = node_new
        symExpr_update_tree(node)
    return 


def symExp_replaceFunction_withArgs_walk_tree(node, f, f_args, expr_new):
    nodes_list = [node]
    ind_next = 0
    while True:
        symExp_replaceFunction_withArgs_node(nodes_list[ind_next], f, f_args, expr_new)
        node_next = nodes_list[ind_next]
        for arg in node_next[3]:
            if len(arg)>0:
                nodes_list.append(arg)
        ind_next += 1
        if ind_next>=len(nodes_list):
            break
    return 

def symExp_replaceFunction_withArgs(expr, f, f_args, f_expr_new):
    """It searches the expr for function f and replaces it with the new expression expr_new
        
    """
    expr_tree = symExpr_generate_tree(expr)
    symExp_replaceFunction_withArgs_walk_tree(expr_tree, f, f_args, f_expr_new)
    return expr_tree[1]


##-------------------------------------

def symExp_changeFunctionName_walk_tree_INACTIVE(node, f, f_new_name):
    # f_new_name: string
    #print(node[1])
    if node[1].func==f:
        node[1] = Function(f_new_name)(*node[1].args)
        symExpr_update_tree(node)
    for arg in node[3]:
        if len(arg)>0:
            symExp_changeFunctionName_walk_tree(arg, f, f_new_name)
    return 

def symExp_changeFunctionName_node(node, f, f_new_name):
    # f_new_name: string
    #print(node[1])
    if node[1].func==f:
        node[1] = Function(f_new_name)(*node[1].args)
        symExpr_update_tree(node)
    return

def symExp_changeFunctionName_walk_tree(node, f, f_new_name):
    # f_new_name: string
    #print(node[1])
    nodes_list = [node]
    ind_next = 0
    while True:
        symExp_changeFunctionName_node(nodes_list[ind_next], f, expr_new)
        node_next = nodes_list[ind_next]
        for arg in node_next[3]:
            if len(arg)>0:
                nodes_list.append(arg)
        ind_next += 1
        if ind_next>=len(nodes_list):
            break
    return 

def symExp_changeFunctionName(expr, f, f_new_name):
    """It searches the expr for function f and replaces it with the new name f_new_name
       f_new_name : string
    """
    expr_tree = symExpr_generate_tree(expr)
    symExp_changeFunctionName_walk_tree(expr_tree, f, f_new_name)
    return expr_tree[1]


##------------------
## find functions with a specific argument

def symExp_CompareFunctionArgs_node(node, f_args, nodes_match):
    """ nodes_match: output list containing nodes whose arguments are f_args
    """
    if node[1].args==f_args:
        nodes_match.append(node)
    return


def symExp_CompareFunctionArgs_walk_tree(node, f_args):
    nodes_list = [node]
    nodes_match = []
    ind_next = 0
    while True:
        symExp_CompareFunctionArgs_node(nodes_list[ind_next], f_args, nodes_match)
        node_next = nodes_list[ind_next]
        for arg in node_next[3]:
            if len(arg)>0:
                nodes_list.append(arg)
        ind_next += 1
        if ind_next>=len(nodes_list):
            break
    return nodes_match

def symExp_CompareFunctionArgs(expr, f_args):
    """It searches the expr for functions with arguments f_arg
    """
    expr_tree = symExpr_generate_tree(expr)
    nodes_match = symExp_CompareFunctionArgs_walk_tree(expr_tree, f_args)
    return [nodes_match[i][1] for i in range(len(nodes_match))]


##------------------
## find functions with a specific argument and replace argument

def symExp_FindReplaceFunctionArgs_node(node, f_args, f_args_new, found_list, replaced_list):
    """ found_list: output list containing nodes whose arguments are f_args
    """
    if node[1].args==f_args:
        found_list.append(node[1])
        #print(node[1], node[1].func)
        node[1] = (node[1].func)(*f_args_new)
        symExpr_update_tree(node)
        replaced_list.append(node[1])
    return


def symExp_FindReplaceFunctionArgs_walk_tree(node, f_args, f_args_new):
    nodes_list = [node]
    found_list, replaced_list = [], []
    ind_next = 0
    while True:
        symExp_FindReplaceFunctionArgs_node(nodes_list[ind_next], f_args, f_args_new, found_list, replaced_list)
        node_next = nodes_list[ind_next]
        for arg in node_next[3]:
            if len(arg)>0:
                nodes_list.append(arg)
        ind_next += 1
        if ind_next>=len(nodes_list):
            break
    return found_list, replaced_list
    

def symExp_FindReplaceFunctionArgs(expr, f_args, f_args_new):
    """It searches the expr for functions with arguments f_arg
    """
    expr_tree = symExpr_generate_tree(expr)
    found_list, replaced_list = symExp_FindReplaceFunctionArgs_walk_tree(expr_tree, f_args, f_args_new)
    return expr_tree[1], found_list, replaced_list



###TODO add code for adding sympy series




##-----  change of variables in differenitial equations  ---







