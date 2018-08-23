## solidstate_sym.py


from sympy import Symbol, symbols, Matrix, I, latex
from sympy.utilities.lambdify import lambdastr
from IPython.display import display, Math, Latex
from Electromagnetics.VectorCalculus import *
import sympy
from Electromagnetics import Misc


__all__ = ["FermiDistribution_sym",
            "FerriteModel_sym"]


def FermiDistribution_sym(E=Symbol('E'), mu_c=Symbol('\mu_c'), T=Symbol('T'), k_B=Symbol('k_B')):
    return 1/(1 + exp(E-mu_c)/k_B*T)





##-----------------  Ferrites  -------------------

def FerriteModel_sym():
    """ Derivation of magnetically biased ferrite permeability model
    """
    x, y, z = symbols('x y z')
    Mx = Symbol('M_{x}')
    My = Symbol('M_{y}')
    Mz = Symbol('M_{z}')
    Hx = Symbol('H_{x}')
    Hy = Symbol('H_{y}')
    Hz = Symbol('H_{z}')
    Ms = Symbol('M_{s}')
    B0 = Symbol('B_{0}')
    alpha = Symbol('\\alpha', real=True)
    gamma = Symbol('\\gamma', real=True)
    omega = Symbol('\\omega', real=True)
    omega_0 = Symbol('\\omega_0', real=True)
    omega_m = Symbol('\\omega_m', real=True)
    mu0, mu_r = symbols('\\mu_0 \\mu_r')

    M = Matrix([[Mx, My, Mz]])
    H = Matrix([[Hx, Hy, Hz]])
    a_z = Matrix([[0, 0, 1]])

    LHS = I*omega*M + gamma*crossproduct(M, B0*a_z) - I*omega*alpha*crossproduct(a_z, M)
    RHS = -gamma*mu0*crossproduct(Ms*a_z, H)

    display(Math('LHS = ' + latex(LHS)))
    display(Math('RHS = ' + latex(RHS)))

    A = Matrix([[LHS[0].subs([(Mx, 1), (My, 0)]), LHS[0].subs([(Mx, 0), (My, 1)])], 
                [LHS[1].subs([(Mx, 1), (My, 0)]), LHS[1].subs([(Mx, 0), (My, 1)])]])
    b = Matrix([RHS[0], RHS[1]])

    display(Math('A = ' + latex(A)))
    display(Math('b = ' + latex(b)))

    x = A.inv()*b
    display(Math('x = ' + latex(x)))

    chi_xx = x[0].subs([(Hx, 1), (Hy, 0)]).simplify()
    chi_xy = x[0].subs([(Hx, 0), (Hy, 1)]).simplify()
    chi_yx = x[1].subs([(Hx, 1), (Hy, 0)]).simplify()
    chi_yy = x[1].subs([(Hx, 0), (Hy, 1)]).simplify()

    display(Math('\\chi_{xx} = ' + latex(chi_xx)))
    display(Math('\\chi_{xy} = ' + latex(chi_xy)))
    display(Math('\\chi_{yx} = ' + latex(chi_yx)))
    display(Math('\\chi_{yy} = ' + latex(chi_yy)))

    chi_xx = chi_xx.subs([(Ms, omega_m/(gamma*mu0)), (B0, omega_0/gamma)])
    chi_xy = chi_xy.subs([(Ms, omega_m/(gamma*mu0)), (B0, omega_0/gamma)])
    chi_yx = chi_yx.subs([(Ms, omega_m/(gamma*mu0)), (B0, omega_0/gamma)])
    chi_yy = chi_yy.subs([(Ms, omega_m/(gamma*mu0)), (B0, omega_0/gamma)])

    display(Math(latex(omega_0) + ' = ' + latex(gamma*B0)))
    display(Math(latex(omega_m) + ' = ' + latex(gamma*mu0*Ms)))

    print('-'*20)
    display(Math('\\chi_{xx} = ' + latex(chi_xx)))
    display(Math('\\chi_{xy} = ' + latex(chi_xy)))
    display(Math('\\chi_{yx} = ' + latex(chi_yx)))
    display(Math('\\chi_{yy} = ' + latex(chi_yy)))

    chi_xx_str = lambdastr((omega, omega_0, omega_m, alpha), chi_xx).replace('\\', '')
    chi_xy_str = lambdastr((omega, omega_0, omega_m, alpha), chi_xy).replace('\\', '')

    print('-'*30)
    print('\nchi_xx_str: \n', chi_xx_str)
    print('\nchi_xy_str: \n', chi_xy_str)
    
    chi_xx_str = Misc.replace_whole_word(chi_xx_str, 'I', '1j')
    chi_xy_str = Misc.replace_whole_word(chi_xy_str, 'I', '1j')
    
    print('-'*10, 'replacing special characters')
    print('\nchi_xx_str: \n', chi_xx_str)
    print('\nchi_xy_str: \n', chi_xy_str)

    print('-'*20)
    chi_xx_r = sympy.re(chi_xx).simplify()
    chi_xx_i = sympy.im(chi_xx).simplify()
    display(Math('\\Re{\\chi_{xx}} = ' + latex(chi_xx_r)))
    display(Math('\\Im{\\chi_{xx}} = ' + latex(chi_xx_i)))

    chi_xy_r = sympy.re(chi_xy).simplify()
    chi_xy_i = sympy.im(chi_xy).simplify()
    display(Math('\\Re{\\chi_{xy}} = ' + latex(chi_xy_r)))
    display(Math('\\Im{\\chi_{xy}} = ' + latex(chi_xy_i)))
    
    chi_xx_r_str = lambdastr((omega, omega_0, omega_m, alpha), chi_xx_r).replace('\\', '')
    chi_xx_i_str = lambdastr((omega, omega_0, omega_m, alpha), chi_xx_i).replace('\\', '')
    chi_xy_r_str = lambdastr((omega, omega_0, omega_m, alpha), chi_xy_r).replace('\\', '')
    chi_xy_i_str = lambdastr((omega, omega_0, omega_m, alpha), chi_xy_i).replace('\\', '')

    print('\nchi_xx_r_str: \n', chi_xx_r_str)
    print('\nchi_xx_i_str: \n', chi_xx_i_str)
    print('\nchi_xy_r_str: \n', chi_xy_r_str)
    print('\nchi_xy_i_str: \n', chi_xy_i_str)
    
    return
    
    
##------------------ Drude Model ----------



    
    
