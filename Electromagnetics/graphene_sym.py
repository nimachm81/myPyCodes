## graphene_sym.py

from sympy import *

__all__ = ["condKuboLorentzian", "densityOfStates", "carrierDensity",
           "InfPdicFrStGrapheneMP_Sym"]

## Kubo conductivity

def condKuboLorentzian(mu_c=Symbol('\mu_c'), B_0=Symbol('B_0'), Gamma=Symbol('\Gamma'), 
    omega=Symbol('\omega'), f_d=Function('f_d'), hbar=Symbol('hbar'), Delta=Symbol('\Delta'), 
    v_F=Symbol('v_F'), e=Symbol('e'), M=Function('M'), n=Symbol('n')):
    """
    returns graphene Lorentzian Kubo conductivity in symbolic form
    time convention: exp(+i*omega*t) is assumed
    """
    sigma_d_coeff = e**2*v_F**2*abs(e*B_0)*(omega-I*2*Gamma)*hbar/(-I*pi)
    
    sigma_d_intra_arg = (f_d(M(n)) - f_d(M(n+1)) + f_d(-M(n+1)) - f_d(-M(n)))\
    /((M(n+1)-M(n))**2 - (omega-I*2*Gamma)**2*hbar**2) \
    *(1-Delta**2/(M(n)*M(n+1)))*(1/(M(n+1)-M(n)))
    sigma_d_inter_arg = (f_d(-M(n)) - f_d(M(n+1)) + f_d(-M(n+1)) - f_d(M(n)))\
    /((M(n+1)+M(n))**2 - (omega-I*2*Gamma)**2*hbar**2) \
    *(1+Delta**2/(M(n)*M(n+1)))*(1/(M(n+1)+M(n)))
    
    sigma_o_coeff = -e**2*v_F**2*e*B_0/pi
    
    sigma_o_intra_arg = (f_d(M(n)) - f_d(M(n+1)) - f_d(-M(n+1)) + f_d(-M(n)))\
    /((M(n+1)-M(n))**2 - (omega-I*2*Gamma)**2*hbar**2) \
    *(1-Delta**2/(M(n)*M(n+1)))
    
    sigma_o_inter_arg = (f_d(M(n)) - f_d(M(n+1)) - f_d(-M(n+1)) + f_d(-M(n)))\
    /((M(n+1)+M(n))**2 - (omega-I*2*Gamma)**2*hbar**2) \
    *(1+Delta**2/(M(n)*M(n+1)))
    
    sigma_d = Sum(sigma_d_coeff*(sigma_d_intra_arg+sigma_d_inter_arg), (n,0,oo))
    sigma_o = Sum(sigma_o_coeff*(sigma_o_intra_arg+sigma_o_inter_arg), (n,0,oo))
    
    M_n = sqrt(Delta**2 + 2*n*v_F**2*abs(e*B_0)*hbar)
    
    return [sigma_d, sigma_o, M_n]
    
    
    
def densityOfStates(E=Symbol('E'), hbar=Symbol('hbar'), v_F=Symbol('v_F')):
    g_v = Symbol('g_v')
    g_s = Symbol('g_s')
    DOS = g_s*g_v/(2*pi*(hbar*v_F)**2)*abs(E)
    return DOS


def carrierDensity(E=Symbol('E'), T=Symbol('T'), f_d=Function('f_d'),
    hbar=Symbol('hbar'), v_F=Symbol('v_F')):
    DOS = densityOfStates(E=E)
    n = Integral(DOS*f_d(E), (E, 0, oo))
    p = Integral(DOS*f_d(E), (E, -oo, 0))
    return [n, p]


from Electromagnetics.VectorCalculus import *
from Electromagnetics.FourierBloch import *
from Electromagnetics.SymExprTree import *
from IPython.display import display, Math, Latex
from scipy import constants
import math

class InfPdicFrStGrapheneMP_Sym:
    """
    Infinite periodic free standing graphene magnetoplasmons
    """
    ## TODO: test eps_r!=1 or mu_r!=1
    def __init__(self, vb=False):
        self.x, self.y, self.z = symbols('x y z')
        self.Ex1 = Symbol('E_{x1}')
        self.Ey1 = Symbol('E_{y1}')
        self.Ez1 = Symbol('E_{z1}')
        self.alpha = Symbol('alpha')
        self.k, self.k0 = symbols('k k_0')
        self.omega = Symbol('\\omega')
        self.eps, self.mu, self.eta = symbols('\\epsilon \\mu \\eta')
        self.sigma_d, self.sigma_o = symbols('\\sigma_d, \\sigma_o')
        self.K_0z = Symbol('K_{0z}') 
        self.n_str = 'n'
        self.vbose = vb
        self.setBoundaryConds()
        return

    def getFields(self):
        E1 = Matrix([[self.Ex1, self.Ey1, self.Ez1]])\
            *exp(-I*self.k*self.z-self.alpha*self.y)       # upper region
        E2 = Matrix([[self.Ex1, -self.Ey1, self.Ez1]])\
            *exp(-I*self.k*self.z+self.alpha*self.y)       # upper region
        if self.vbose:
            display(Math('E_1 = ' + latex(E1)))
            display(Math('E_2 = ' + latex(E2)))
        H1 = -1/(I*self.omega*self.mu)*curl_r(E1)
        H2 = -1/(I*self.omega*self.mu)*curl_r(E2)
        if self.vbose:
            display(Math('H_1 = ' + latex(H1)))
            display(Math('H_2 = ' + latex(H2)))
        return [E1, E2, H1, H2]

    def getSigmaTensor(self):
        sigma = Matrix([[self.sigma_d, 0, -self.sigma_o], 
                        [0, 0, 0], 
                        [self.sigma_o, 0, self.sigma_d]])
        if self.vbose:
            display(Math('\\sigma = ' + latex(sigma)))
        return sigma
        
    def setBoundaryConds(self):
        E1, E2, H1, H2 = self.getFields()
        sigma = self.getSigmaTensor()
        a_y = Matrix([[0, 1, 0]])

        #boundary condition
        BC = crossproduct(a_y, (H1-H2)) - (sigma*E1.T).T    
        if self.vbose:
            print('Boundary condition:')
            display(Math('BC = ' + latex(BC.T)))
                
        n = Symbol(self.n_str)
        self.harmonic = exp(I*n*self.K_0z*self.z)
        vars_fourier = [self.Ex1, self.Ey1, self.Ez1, self.sigma_d, self.sigma_o]
        BC = d1_putSums(BC, vars_fourier, self.z, self.harmonic)
        if self.vbose:
            display(Math('BC = ' + latex(BC.T)))
        
        self.Ex1_tilde, self.Ey1_tilde, self.Ez1_tilde, self.sigma_d_tilde,\
            self.sigma_o_tilde = vars_fourier
        if self.vbose:
            print('Fourier Variables: ', vars_fourier)

        BC = Matrix([[BC[i].doit() for i in range(BC.cols)]])
        if self.vbose:
            display(Math('BC = ' + latex(BC.T)))
        BC = BC.subs(self.y, 0)
        if self.vbose:
            print('setting y=0')
            display(Math('BC = ' + latex(BC.T)))

        BC = d1_applyConvolutions(BC, self.z, self.harmonic)
        if self.vbose:
            print('applying convolutions')
            display(Math('BC = ' + latex(BC.T)))

        BC = Matrix([[d1_applyOrthogonalities(BC[i], self.z, self.harmonic)\
                for i in range(BC.cols)]])
        if self.vbose:
            display(Math('BC = ' + latex(BC.T)))

        # alpha k relation
        AK = (del_square_vec_r(E1)+ self.k0**2*E1)     
        if self.vbose:
            print('alpha-k relation:')
            display(Math('AK = ' + latex(AK.T)))

        vars_fourier = [self.Ex1, self.Ey1, self.Ez1, self.sigma_d, self.sigma_o]
        AK = d1_putSums(AK, vars_fourier, self.z, self.harmonic)
        if self.vbose:
            print('alpha-k relation:')
            display(Math('AK = ' + latex(AK.T)))

        AK = Matrix([[AK[i].doit() for i in range(AK.cols)]])
        if self.vbose:
            print('alpha-k relation:')
            display(Math('AK = ' + latex(AK.T)))
        AK = Matrix([[(d1_applyOrthogonalities(AK[i], self.z, self.harmonic)*\
        exp(I*self.k*self.z+self.alpha*self.y)).expand() for i in range(AK.cols)]])
        if self.vbose:
            display(Math('AK = ' + latex(AK.T)))
        # normal condition
        NC = divergence_r(E1)    
        if self.vbose:
            print('getting normal component')
            display(Math('NC = ' + latex(NC)))
        NC = d1_putSums(NC, [self.Ex1, self.Ey1, self.Ez1, self.sigma_d, self.sigma_o],
            self.z, self.harmonic)
        if self.vbose:
            display(Math('NC = ' + latex(NC)))
        NC = NC.doit()
        if self.vbose:
            display(Math('NC = ' + latex(NC)))
        #Ey1_rep = solve(NC, Ey1)[0]
        if self.vbose:
            print('applying convolutions and orthogonalities')
        NC = d1_applyConvolutions(NC, self.z, self.harmonic)
        NC = d1_applyOrthogonalities(NC, self.z, self.harmonic)

        NC = simplify(NC.subs(self.y, 0)*exp(I*self.k*self.z))
        if self.vbose:
            display(Math('NC = ' + latex(NC)))

        self.Ey1_rep = solve(NC, self.Ey1_tilde(n))[0]
        Ey1_rep = self.Ey1_rep      ## same for upper and lower regions
        if self.vbose:
            display(Math(latex(self.Ey1_tilde) + ' = ' + latex(Ey1_rep)))

        BC = Matrix([[symExp_replaceFunction(BC[i], self.Ey1_tilde, Ey1_rep)\
            for i in range(BC.cols)]])*exp(I*self.k*self.z)
        if self.vbose:
            display(Math('BC = ' + latex(BC.T)))
            display(Math('BC = ' + latex(BC.T.subs(self.z,0))))
        
        BC = Matrix([[BC[i].simplify() for i in range(BC.cols)]])
        if self.vbose:
            display(Math('BC = ' + latex(BC.T)))

        self.BC = [BC[0].expand(), BC[2].expand()]
        self.AK = [AK[0].expand()]

        return [self.BC, self.AK]

    
    def getBoundaryConds(self):
        return [self.BC, self.AK]
        
    def subsParams(self, freq, Z_period, eps_r=1.0, mu_r=1.0):
        k_0 = 2.0*math.pi*freq*math.sqrt(constants.mu_0*mu_r*constants.epsilon_0*eps_r)
        
        BC_fin = [self.BC[i].subs([(self.omega, Float(2.0*math.pi*freq, 20)), (self.mu, Float(constants.mu_0*mu_r, 20)), 
            (self.K_0z, Float(2.0*math.pi/Z_period, 20)), (I, 1j)]) for i in range(len(self.BC))]
        AK_fin = [self.AK[i].subs([(self.k0, Float(k_0, 20)), (self.K_0z, Float(2.0*math.pi/Z_period, 20)),
                (I, 1j)])  for i in range(len(self.AK))]

        if self.vbose:
            display(Math('BC = ' + latex(Matrix(BC_fin))))
            display(Math('BC = ' + latex(expand(Matrix(BC_fin)))))
            display(Math('AK = ' + latex(Matrix(AK_fin))))

        eq_list_bc = [BC_fin[0].expand(), BC_fin[1].expand()]
        eq_list_ak = [AK_fin[0].expand()]
        return [eq_list_bc, eq_list_ak]

    def getVarsPars(self):
        vars_bc = [self.Ex1_tilde, self.Ez1_tilde]
        pars_bc = [self.sigma_d_tilde, self.sigma_o_tilde]

        vars_ak = [self.Ex1_tilde]
        pars_ak = []
        return [vars_bc, pars_bc, vars_ak, pars_ak]
        
    def getEigVars(self):
        return [self.k, self.alpha]
        
    def getElectricFieldHarmonics(self):
        E1, E2, H1, H2 = self.getFields()
        a_y = Matrix([[0, 1, 0]])
        n = Symbol(self.n_str)
        vars_fourier = [self.Ex1, self.Ey1, self.Ez1, self.sigma_d, self.sigma_o]
        E1 = d1_putSums(E1, vars_fourier, self.z, self.harmonic)
        vars_fourier = [self.Ex1, self.Ey1, self.Ez1, self.sigma_d, self.sigma_o]
        E2 = d1_putSums(E2, vars_fourier, self.z, self.harmonic)
        if self.vbose:
            display(Math('E_1 = ' + latex(E1.T)))
            display(Math('E_2 = ' + latex(E2.T)))
        E1 = Matrix([[E1[i].doit() for i in range(E1.cols)]])
        E2 = Matrix([[E2[i].doit() for i in range(E2.cols)]])
        
        E1 = d1_applyConvolutions(E1, self.z, self.harmonic)
        E2 = d1_applyConvolutions(E2, self.z, self.harmonic)
        if self.vbose:
            display(Math('E_1 = ' + latex(E1.T)))
            display(Math('E_2 = ' + latex(E2.T)))
        
        E1 = Matrix([[d1_applyOrthogonalities(E1[i], self.z, self.harmonic)\
                for i in range(E1.cols)]])
        E2 = Matrix([[d1_applyOrthogonalities(E2[i], self.z, self.harmonic)\
                for i in range(E2.cols)]])
        if self.vbose:
            display(Math('\\tilde{E}_1 = ' + latex(E1.T)))
            display(Math('\\tilde{E}_2 = ' + latex(E2.T)))
        E1 = Matrix([[simplify(symExp_replaceFunction(E1[i], self.Ey1_tilde, self.Ey1_rep)\
            /exp(-I*self.k*self.z-self.alpha*self.y)) for i in range(H1.cols)]])
        E2 = Matrix([[simplify(symExp_replaceFunction(E2[i], self.Ey1_tilde, self.Ey1_rep)\
            /exp(-I*self.k*self.z+self.alpha*self.y)) for i in range(H2.cols)]])
        if self.vbose:
            display(Math('\\tilde{E}_1 = ' + latex(E1.T)))
            display(Math('\\tilde{E}_2 = ' + latex(E2.T)))
        return [E1, E2]

    def getElectricFieldHarmonics_subsParams(self, freq, Z_period, k, alpha, 
            eps_r=1.0, mu_r=1.0):
        E1, E2 = self.getElectricFieldHarmonics()

        k_0 = 2.0*math.pi*freq*math.sqrt(constants.mu_0*mu_r*constants.epsilon_0*eps_r)
        E1 = E1.subs([(self.omega, Float(2.0*math.pi*freq, 20)), (self.mu, Float(constants.mu_0*mu_r, 20)),
            (self.eps, Float(constants.epsilon_0*eps_r, 20)), (self.k, Float(k.real, 20)+I*Float(k.imag, 20)), 
            (self.alpha, Float(alpha.real, 20)+I*Float(alpha.imag, 20)), 
            (self.K_0z, Float(2.0*math.pi/Z_period, 20)), (I, 1j)])
        E2 = E2.subs([(self.omega, Float(2.0*math.pi*freq, 20)), (self.mu, Float(constants.mu_0*mu_r, 20)),
            (self.eps, Float(constants.epsilon_0*eps_r, 20)), (self.k, Float(k.real, 20)+I*Float(k.imag, 20)), 
            (self.alpha, Float(alpha.real, 20)+I*Float(alpha.imag, 20)), 
            (self.K_0z, Float(2.0*math.pi/Z_period, 20)), (I, 1j)])
        if self.vbose:
            display(Math('\\tilde{E}_1 = ' + latex(E1.T)))
            display(Math('\\tilde{E}_2 = ' + latex(E2.T)))
        return [E1, E2]   
    
    def getMagneticFieldHarmonics(self):
        E1, E2, H1, H2 = self.getFields()
        a_y = Matrix([[0, 1, 0]])
        n = Symbol(self.n_str)
        vars_fourier = [self.Ex1, self.Ey1, self.Ez1, self.sigma_d, self.sigma_o]
        H1 = d1_putSums(H1, vars_fourier, self.z, self.harmonic)
        vars_fourier = [self.Ex1, self.Ey1, self.Ez1, self.sigma_d, self.sigma_o]
        H2 = d1_putSums(H2, vars_fourier, self.z, self.harmonic)
        if self.vbose:
            display(Math('H_1 = ' + latex(H1.T)))
            display(Math('H_2 = ' + latex(H2.T)))
        H1 = Matrix([[H1[i].doit() for i in range(H1.cols)]])
        H2 = Matrix([[H2[i].doit() for i in range(H2.cols)]])

        H1 = d1_applyConvolutions(H1, self.z, self.harmonic)
        H2 = d1_applyConvolutions(H2, self.z, self.harmonic)
        if self.vbose:
            display(Math('H_1 = ' + latex(H1.T)))
            display(Math('H_2 = ' + latex(H2.T)))

        H1 = Matrix([[d1_applyOrthogonalities(H1[i], self.z, self.harmonic)\
                for i in range(H1.cols)]])
        H2 = Matrix([[d1_applyOrthogonalities(H2[i], self.z, self.harmonic)\
                for i in range(H2.cols)]])
        if self.vbose:
            display(Math('\\tilde{H}_1 = ' + latex(H1.T)))
            display(Math('\\tilde{H}_2 = ' + latex(H2.T)))
        H1 = Matrix([[simplify(symExp_replaceFunction(H1[i], self.Ey1_tilde, self.Ey1_rep)\
            /exp(-I*self.k*self.z-self.alpha*self.y)) for i in range(H1.cols)]])
        H2 = Matrix([[simplify(symExp_replaceFunction(H2[i], self.Ey1_tilde, self.Ey1_rep)\
            /exp(-I*self.k*self.z+self.alpha*self.y)) for i in range(H2.cols)]])
        if self.vbose:
            display(Math('\\tilde{H}_1 = ' + latex(H1.T)))
            display(Math('\\tilde{H}_2 = ' + latex(H2.T)))
        return [H1, H2]
        
    def getMagneticFieldHarmonics_subsParams(self, freq, Z_period, k, alpha, 
            eps_r=1.0, mu_r=1.0):
        H1, H2 = self.getMagneticFieldHarmonics()

        k_0 = 2.0*math.pi*freq*math.sqrt(constants.mu_0*mu_r*constants.epsilon_0*eps_r)
        H1 = H1.subs([(self.omega, Float(2.0*math.pi*freq, 20)), (self.mu, Float(constants.mu_0*mu_r, 20)),
            (self.eps, Float(constants.epsilon_0*eps_r, 20)), (self.k, Float(k.real, 20)+I*Float(k.imag, 20)), 
            (self.alpha, Float(alpha.real, 20)+I*Float(alpha.imag, 20)), 
            (self.K_0z, Float(2.0*math.pi/Z_period, 20)), (I, 1j)])
        H2 = H2.subs([(self.omega, Float(2.0*math.pi*freq, 20)), (self.mu, Float(constants.mu_0*mu_r, 20)),
            (self.eps, Float(constants.epsilon_0*eps_r, 20)), (self.k, Float(k.real, 20)+I*Float(k.imag, 20)), 
            (self.alpha, Float(alpha.real, 20)+I*Float(alpha.imag, 20)), 
            (self.K_0z, Float(2.0*math.pi/Z_period, 20)), (I, 1j)])
        if self.vbose:
            display(Math('\\tilde{H}_1 = ' + latex(H1.T)))
            display(Math('\\tilde{H}_2 = ' + latex(H2.T)))
        return [H1, H2]   
    
    
    def getCurrentHarmonics(self):
        E1, E2, H1, H2 = self.getFields()
        sigma = self.getSigmaTensor()

        #boundary condition
        J = (sigma*E1.T).T    
        if self.vbose:
            print('Current density:')
            display(Math('J = ' + latex(J.T)))
        n = Symbol(self.n_str)
        vars_fourier = [self.Ex1, self.Ey1, self.Ez1, self.sigma_d, self.sigma_o]
        J = d1_putSums(J, vars_fourier, self.z, self.harmonic)
        if self.vbose:
            display(Math('J = ' + latex(J.T)))
        J = Matrix([[J[i].doit() for i in range(J.cols)]])

        J = d1_applyConvolutions(J, self.z, self.harmonic)
        if self.vbose:
            display(Math('J = ' + latex(J.T)))

        J = Matrix([[d1_applyOrthogonalities(J[i], self.z, self.harmonic)\
                for i in range(J.cols)]])
        if self.vbose:
            display(Math('\\tilde{J} = ' + latex(J.T)))
        J = Matrix([[simplify(symExp_replaceFunction(J[i], self.Ey1_tilde, self.Ey1_rep)\
            /exp(-I*self.k*self.z-self.alpha*self.y)) for i in range(J.cols)]])
        if self.vbose:
            display(Math('\\tilde{J} = ' + latex(J.T)))
        for i in range(J.cols):
            if J[i].func==Sum:
                J[i] = J[i].args[0]
        if self.vbose:
            print("removing sums in front of convolutions")
            display(Math('\\tilde{J} = ' + latex(J.T)))
        return J
                
    


