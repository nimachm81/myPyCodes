## graphene.py


from scipy import constants, integrate, special
from numpy import *
from scipy import *
from Electromagnetics.solidstate import *
import math
import cmath
import numpy as np
import sympy
from sympy import lambdify, Symbol
from Electromagnetics.FourierBloch import *
from Electromagnetics.Misc import null, solveMuller, RootsMultipleComplex
from Electromagnetics.SymExprTree import symExp_changeFunctionName
from scipy import linalg, fftpack
from scipy.optimize import newton_krylov



__all__ = ["condKuboLorentzian", "condKuboNoMagneticField", "condDrude",
            "LandauLevel", "densityOfStates", "carrierDensity",
            "magnetoplasmonDispersion_cond_normalized",
            "mobilityToScatteringTime", 
            "magnetoplasmonDispersion_normalized", "InfPdicFrStGrapheneMP", "plasmonDispersion",
            "plasmonDispersion_normalized", "plasmonDispersion_normalized_cond", 
            "plasmonDispersion_TE_normalized", "plasmonDispersion_TE_normalized_cond",
            "GrapheneOnFerriteSubstrate", "GrapheneOnFerriteSubstrateSuperstrate"]


fermi_velocity = 1.0e6


## Kubo conductivity

def condKuboLorentzian(mu_c, B_0, tau, omega, T, Delta=0.0, 
    seperate_intra_inter=False):
    """
    returns graphene Lorentzian Kubo conductivity
    time convention: exp(+i*omega*t) is assumed
    """
    #TODO: ensure convergence
    if B_0==0.0:
        if seperate_intra_inter==False:
            sigma_d = condKuboNoMagneticField(mu_c, tau, omega, T)
            sigma_o = np.array([0.0+0j]*len(omega))
            return [sigma_d, sigma_o]
        else:
            sigma_d_intra_inter = condKuboNoMagneticField(mu_c, tau, omega,
                 T, seperate_intra_inter)
            sigma_o_intra_inter = [np.array([0.0+0j]*len(omega)), 
                    np.array([0.0+0j]*len(omega))] 
            return [sigma_d_intra_inter, sigma_o_intra_inter]
        
        
    v_F = fermi_velocity
    e = -constants.e
    hbar = constants.hbar
    Gamma = 1.0/(2.0*tau)
    if mu_c<0:
        e *= -1
    
    def M(n):
        return sqrt(Delta**2 + 2*n*v_F**2*abs(e*B_0)*hbar)
        
    def f_d(E):
        return FermiDistribution(E, mu_c, T)
        
    decay_min = 1.0e-20 # minimum decay in fermi distribution (for n_max)
    E_max = abs(mu_c) - constants.k*T*log(decay_min)
    #E_max = max(E_max, 10.0*hbar*max(omega))
    #print('E_max: ', E_max/constants.eV , 'eV')
    

    sigma_d_coeff_intra = e**2*(omega-1j*2.0*Gamma)/(-2j*pi)
    sigma_d_coeff_inter = e**2*v_F**2*abs(e*B_0)*(omega-1j*2.0*Gamma)*hbar/(-1j*pi)
    sigma_o_coeff = -e**2*v_F**2*e*B_0/pi
    
    alpha = 2.0*v_F**2*abs(e*B_0)*hbar
    n_max = 10 + int(math.ceil(abs(E_max**2 - Delta**2)/alpha))
    #print('n_max = ', n_max)
    E_max = sqrt(2*n_max*v_F**2*abs(e*B_0)*hbar)
    
    len_w = len(omega)
    sigma_d_intra_arg, sigma_d_inter_arg = [0.0+0j]*len_w, [0.0+0j]*len_w
    sigma_o_intra_arg, sigma_o_inter_arg = [0.0+0j]*len_w, [0.0+0j]*len_w

    hbar_sq = hbar*hbar
    two_j_gamma = 1j*2*Gamma
    alpha_sq = alpha*alpha
    for n in range(n_max):
        M_n, M_np1 = M(n), M(n+1)
        f_d_Mn, f_d_mMn, f_d_Mnp1, f_d_mMnp1 =\
            f_d(M_n), f_d(-M_n), f_d(M_np1), f_d(-M_np1)
        A = f_d_Mn - f_d_Mnp1 + f_d_mMnp1 - f_d_mMn
        B = f_d_mMn - f_d_Mnp1 + f_d_mMnp1 - f_d_Mn
        C = f_d_Mn - f_d_Mnp1 - f_d_mMnp1 + f_d_mMn
        M_np1_p_M_n = M_np1 + M_n
        M_np1_p_M_n_sq = M_np1_p_M_n*M_np1_p_M_n
        Dsq_MnMnp1 = Delta**2
        if Dsq_MnMnp1!=0.0:
            Dsq_MnMnp1 /= (M_n*M_np1)
        
        w_m_2jg_hb_sq = (omega - two_j_gamma)*(omega - two_j_gamma)*hbar_sq
        sigma_d_intra_arg += A /(alpha_sq/M_np1_p_M_n_sq - w_m_2jg_hb_sq) \
        *M_np1_p_M_n*(1.0-Dsq_MnMnp1)
        
        sigma_d_inter_arg += B /(M_np1_p_M_n_sq - w_m_2jg_hb_sq) \
        /M_np1_p_M_n*(1.0+Dsq_MnMnp1)
        
        sigma_o_intra_arg += C / (alpha_sq/M_np1_p_M_n_sq - w_m_2jg_hb_sq)\
        *(1.0-Dsq_MnMnp1) 
        
        sigma_o_inter_arg += C /(M_np1_p_M_n_sq - w_m_2jg_hb_sq)\
        *(1.0+Dsq_MnMnp1)

    ## correction for the interband term  TODO: to be checked for B_0 very large
    ## incorrect for Delta!=0.0 
    ## TODO: similar correction to be added for the off diagonala interband term
    a__ = (omega - 1j*2.0*Gamma)*hbar
    sigma_d_inter_arg += 1.0/(4.0*a__)*log((2.0*E_max+a__)/(2.0*E_max-a__))/(v_F**2*abs(e*B_0)*hbar)


    if seperate_intra_inter==False:
        sigma_d = sigma_d_coeff_intra*sigma_d_intra_arg +\
            sigma_d_coeff_inter*sigma_d_inter_arg
        sigma_o = sigma_o_coeff*(sigma_o_intra_arg + sigma_o_inter_arg)
            
        return [sigma_d, sigma_o]
    else:
        sigma_d_intra = sigma_d_coeff_intra*sigma_d_intra_arg
        sigma_d_inter = sigma_d_coeff_inter*sigma_d_inter_arg
        sigma_o_intra = sigma_o_coeff*sigma_o_intra_arg
        sigma_o_inter = sigma_o_coeff*sigma_o_inter_arg
        
        return [[sigma_d_intra, sigma_d_inter], 
               [sigma_o_intra, sigma_o_inter]]


def condKuboNoMagneticField(mu_c, tau, omega, T, seperate_intra_inter=False):
    """
    diagonal conductivity for B_0==0
    time convention: exp(+i*omega*t) is assumed
    """
    ##TODO: handle exceptions in the integral method
    v_F = fermi_velocity
    e = -constants.e
    hbar = constants.hbar
    Gamma = 1.0/(2.0*tau)
    if mu_c<0:
        e *= -1.0

    f_d = FermiDistribution
    #import cmath

    sigma_d_intra = -1j*e**2*constants.k*T/(pi*hbar**2*(omega-1j*2.0*Gamma))\
    *(mu_c/(constants.k*T) + 2.0*math.log(math.exp(-mu_c/(constants.k*T)) + 1.0))
    
    decay_min = 1.0e-20 # minimum decay in fermi distribution (for n_max)
    E_max = abs(mu_c) - constants.k*T*log(decay_min)
    #print('E_max: ', E_max/constants.eV , 'eV')
    #E_max = np.inf

    sigma_d_inter = np.array([0.0+0j]*len(omega))
    for i in range(len(omega)):
        def arg_inter_re(E):
            return np.real((f_d(E, mu_c, T) - f_d(-E, mu_c, T))/(((omega[i] - 1j*2.0*Gamma)*hbar)**2\
             - 4.0*E**2))
        def arg_inter_im(E):
            return np.imag((f_d(E, mu_c, T) - f_d(-E, mu_c, T))/(((omega[i] - 1j*2.0*Gamma)*hbar)**2\
             - 4.0*E**2))
        #E_end = max(E_max, 10.0*constants.hbar*omega[i])
        integ_res_r = integrate.quad(arg_inter_re, 0, E_max, points=[abs(mu_c)])
        integ_res_i = integrate.quad(arg_inter_im, 0, E_max, points=[abs(mu_c)])
        #print('integration error r: ' , integ_res_r[1] , '+j ', integ_res_i[1])
        integ_res = integ_res_r[0] + 1j*integ_res_i[0]
        a__ = (omega[i] - 1j*2.0*Gamma)*hbar
        integ_res += 1.0/(4.0*a__)*cmath.log((2.0*E_max+a__)/(2.0*E_max-a__))
        sigma_d_inter[i] = integ_res*1j*e**2*(omega[i] - 1j*2.0*Gamma)/(pi)
    
    if seperate_intra_inter==False:
        return sigma_d_intra+sigma_d_inter
    else:
        return [sigma_d_intra, sigma_d_inter]
    



def condDrude(mu_c, B_0, tau, omega, T):
    ##TODO: handle T=0 and mu_c=0
    v_F = fermi_velocity
    e = -constants.e
    hbar = constants.hbar
    Gamma = 1.0/(2.0*tau)
    if mu_c<0:
        e *= -1.0
        
    omega_c = e*B_0*v_F**2/mu_c
    
    sigma_0 = tau*e**2*constants.k*T/(pi*hbar**2)\
    *(mu_c/(constants.k*T) + 2.0*math.log(math.exp(-mu_c/(constants.k*T)) + 1.0))
    
    A = (1.0+1j*omega*tau)**2 + (omega_c*tau)**2
    
    sigma_d = sigma_0*(1.0+1j*omega*tau)/A
    
    sigma_o = -sigma_0*omega_c*tau/A
    
    return [sigma_d, sigma_o]


def LandauLevel(B_0, n, Delta=0.0):
    """
    Landau level n=0,1,...
    """
    v_F = fermi_velocity
    e = -constants.e
    hbar = constants.hbar
    return  sqrt(Delta**2 + 2*n*v_F**2*abs(e*B_0)*hbar)


def densityOfStates(E):
    v_F = fermi_velocity
    hbar = constants.hbar
    g_v = 2
    g_s = 2
    DOS = g_s*g_v/(2.0*pi*(hbar*v_F)**2)*abs(E)
    return DOS


def carrierDensity(mu_c, T):
    ##TODO: handle big exponential arguments
    v_F = fermi_velocity
    hbar = constants.hbar
    if T>0:
        k_B = constants.k
        kT= k_B*T
        C = 2.0/pi*(kT/(hbar*v_F))**2
        n = -C*special.spence(1.0+exp(mu_c/kT))
        p = -C*special.spence(1.0+exp(-mu_c/kT))
        return [n, p]
    elif T==0.0:
        n_p = mu_c**2/(pi*(hbar*v_F)**2)
        n = n_p * ((np.sign(mu_c)+1)/2)
        p = -n_p * ((np.sign(mu_c)-1)/2)
        return [n, p]
    else:
        raise ValueError('negative temperature.')
        

def mobilityToScatteringTime(mu, mu_c):
    """ mu: mobility in SI units 
    mu_c : chemical potential in SI
    """
    v_F = fermi_velocity
    m = mu_c/v_F**2
    tau = mu*m/constants.e
    return tau    


def magnetoplasmonDispersion_cond_normalized(omega, sigma_L, sigma_H):
    """
    gives the dispersion of the magnetoplasmon fields for free standing graphene
    inside vaccum
    two complex arrays are returned. The one with real part greater than 1.0
    would be a surface wave
    normalization : k/k0
    omega: numpy array
    """
    etha_0 = math.sqrt(constants.mu_0/constants.epsilon_0)
    
    s_L = etha_0*sigma_L/2.0
    s_H = etha_0*sigma_H/2.0

    A = 1.0 + s_L**2 + s_H**2
    a_y_k0_p = -1.0j*(A + sqrt(A**2-4.0*s_L**2))/(2.0*s_L)
    a_y_k0_m = -1.0j*(A - sqrt(A**2-4.0*s_L**2))/(2.0*s_L)
        
    k_z_k0_p = sqrt(1.0 + a_y_k0_p**2)
    k_z_k0_m = sqrt(1.0 + a_y_k0_m**2)
    return [k_z_k0_p, k_z_k0_m]
    

def magnetoplasmonDispersion_normalized(mu_c, B_0, omega, tau, T, cond='Kubo'):
    """
    gives the dispersion of the magnetoplasmon fields for free standing graphene
    inside vaccum
    two complex arrays are returned. The one with real part greater than 1.0
    would be a surface wave
    normalization : k/k0
    omega: numpy array
    """
    if cond=='Kubo':
        sigma_d, sigma_o = condKuboLorentzian(mu_c, B_0, tau, omega, T)
        return magnetoplasmonDispersion_cond_normalized(omega, sigma_d, sigma_o)
    elif cond=='Drude':
        sigma_d, sigma_o = condDrude(mu_c, B_0, tau, omega, T)
        return magnetoplasmonDispersion_cond_normalized(omega, sigma_d, sigma_o)
    else:
        raise ValueError('cond: either "Kubo" or "Drude"')
    return
    

def plasmonDispersion_cond(omega, sigma, eps_r1=1.0, mu_r1=1.0, 
        eps_r2=1.0, mu_r2=1.0):
    """ plasmon dispersion for a 2DEG surrounded by 2 media
        eps_r1, mu_r1: relative parameters of the upper medium
        eps_r2, mu_r2: relative parameters of the lower medium
        omega: numpy array
    """
    ## TODO: add non-locality
    eps1 = eps_r1*constants.epsilon_0
    eps2 = eps_r2*constants.epsilon_0
    mu1 = mu_r1*constants.mu_0
    mu2 = mu_r2*constants.mu_0
    if eps_r1==eps_r2 and mu_r1==mu_r2:
        k_TM = sqrt(omega**2*mu1*eps1 + (2.0j*omega*eps1/sigma)**2)
        k_TM = k_TM*(np.imag(sigma)<0.0)
        k_TE = sqrt(omega**2*mu1*eps1 + (sigma*mu1*omega/2.0j)**2)
        k_TE = k_TE*(np.imag(sigma)>0.0)
        k = k_TM + k_TE
        return k
    else:
        raise NotImplementedError() ## TODO: to be done
        
def plasmonDispersion_normalized_cond(omega, sigma, eps_r1=1.0, mu_r1=1.0, 
        eps_r2=1.0, mu_r2=1.0):
    """ plasmon dispersion for a 2DEG surrounded by 2 media
        eps_r1, mu_r1: relative parameters of the upper medium
        eps_r2, mu_r2: relative parameters of the lower medium
        omega: numpy array
    """
    ## TODO: add non-locality
    mu_1 = mu_r1
    epsilon_1 = eps_r1
    eta_0 = math.sqrt(constants.mu_0/constants.epsilon_0)
    sigma = sigma*eta_0
    if eps_r1==eps_r2 and mu_r1==mu_r2:
        k_k0_TM = np.sqrt(epsilon_1*(-4.0*epsilon_1/sigma**2 + mu_1))
        k_k0_TM = k_k0_TM*(np.imag(sigma)<0.0)
        k_k0_TE = np.sqrt(mu_1*(4.0*epsilon_1 - mu_1*sigma**2))/2.0    
        k_k0_TE = k_k0_TE*(np.imag(sigma)>0.0)
        k_k0 = k_k0_TM + k_k0_TE
        return k_k0
    else:
        raise NotImplementedError() ## TODO: to be done

def plasmonDispersion(mu_c, omega, tau, T, eps_r1=1.0, mu_r1=1.0, 
        eps_r2=1.0, mu_r2=1.0, cond='Kubo'):
    """
    gives the dispersion of the oplasmons in free standing graphene
    omega: numpy array
    """
    if cond=='Kubo':
        B_0 = 0.0
        sigma_d, sigma_o = condKuboLorentzian(mu_c, B_0, tau, omega, T)
        return plasmonDispersion_cond(omega, sigma_d, eps_r1=eps_r1, mu_r1=mu_r1, 
            eps_r2=eps_r2, mu_r2=mu_r2)
    elif cond=='Drude':
        B_0 = 0.0
        sigma_d, sigma_o = condDrude(mu_c, B_0, tau, omega, T)
        return plasmonDispersion_cond(omega, sigma_d, eps_r1=eps_r1, mu_r1=mu_r1, 
            eps_r2=eps_r2, mu_r2=mu_r2)
    else:
        raise ValueError('cond: either "Kubo" or "Drude"')
    return

def plasmonDispersion_normalized(mu_c, omega, tau, T, eps_r1=1.0, mu_r1=1.0, 
        eps_r2=1.0, mu_r2=1.0, cond='Kubo'):
    """
    gives the dispersion of the oplasmons in free standing graphene
    omega: numpy array
    """
    if cond=='Kubo':
        B_0 = 0.0
        sigma_d, sigma_o = condKuboLorentzian(mu_c, B_0, tau, omega, T)
        return plasmonDispersion_normalized_cond(omega, sigma_d, eps_r1=eps_r1, mu_r1=mu_r1, 
            eps_r2=eps_r2, mu_r2=mu_r2)
    elif cond=='Drude':
        B_0 = 0.0
        sigma_d, sigma_o = condDrude(mu_c, B_0, tau, omega, T)
        return plasmonDispersion_normalized_cond(omega, sigma_d, eps_r1=eps_r1, mu_r1=mu_r1, 
            eps_r2=eps_r2, mu_r2=mu_r2)
    else:
        raise ValueError('cond: either "Kubo" or "Drude"')
    return

def plasmonDispersion_TE_normalized_cond__(sigma, eps_r1=1.0, mu_r1=1.0, 
        eps_r2=1.0, mu_r2=1.0):
    """ k/k_0 dispersion for TE surface plasmons
        eps_r1, mu_r1: relative parameters of the upper medium
        eps_r2, mu_r2: relative parameters of the lower medium
        calculates k directly
    """
    ##TODO: proper and improper reults to be checked
    if eps_r1==eps_r2 and mu_r1==mu_r2:
        mu_1 = mu_r1
        epsilon_1 = eps_r1
        eta_0 = math.sqrt(constants.mu_0/constants.epsilon_0)
        sigma = sigma*eta_0
        k_k0_TE = np.sqrt(mu_1*(4.0*epsilon_1 - mu_1*sigma**2))/2.0    
        return k_k0_TE
    elif mu_r1==mu_r2:
        mu_1 = mu_r1
        epsilon_1 = eps_r1
        epsilon_2 = eps_r2
        eta_0 = math.sqrt(constants.mu_0/constants.epsilon_0)
        sigma = sigma*eta_0
        
        k_k0_TE = np.sqrt((-epsilon_1**2 + 2.0*epsilon_1*epsilon_2 + 2.0*epsilon_1*mu_1*sigma**2 - epsilon_2**2 + 2.0*epsilon_2*mu_1*sigma**2 - mu_1**2*sigma**4)/sigma**2)/2.0
        return k_k0_TE
    else:
        mu_1 = mu_r1
        epsilon_1 = eps_r1
        mu_2 = mu_r2
        epsilon_2 = eps_r2
        eta_0 = math.sqrt(constants.mu_0/constants.epsilon_0)
        sigma = sigma*eta_0
        #k_k0_TE = (np.sqrt((-epsilon_1*mu_1**3*mu_2**2 + epsilon_1*mu_1*mu_2**4 + epsilon_2*mu_1**4*mu_2 - epsilon_2*mu_1**2*mu_2**3 - mu_1**4*mu_2**2*sigma**2 - mu_1**2*mu_2**4*sigma**2 - 2*np.sqrt(mu_1**4*mu_2**4*sigma**2*(epsilon_1*mu_1**3 - epsilon_1*mu_1*mu_2**2 - epsilon_2*mu_1**2*mu_2 + epsilon_2*mu_2**3 + mu_1**2*mu_2**2*sigma**2)))/(mu_1**4 - 2*mu_1**2*mu_2**2 + mu_2**4)))
        k_k0_TE = (np.sqrt((-epsilon_1*mu_1**3*mu_2**2 + epsilon_1*mu_1*mu_2**4 + epsilon_2*mu_1**4*mu_2 - epsilon_2*mu_1**2*mu_2**3 - mu_1**4*mu_2**2*sigma**2 - mu_1**2*mu_2**4*sigma**2 + 2*np.sqrt(mu_1**4*mu_2**4*sigma**2*(epsilon_1*mu_1**3 - epsilon_1*mu_1*mu_2**2 - epsilon_2*mu_1**2*mu_2 + epsilon_2*mu_2**3 + mu_1**2*mu_2**2*sigma**2)))/(mu_1**4 - 2*mu_1**2*mu_2**2 + mu_2**4)))
        return k_k0_TE

def plasmonDispersion_TE_normalized_cond(sigma, eps_r1=1.0, mu_r1=1.0, 
        eps_r2=1.0, mu_r2=1.0):
    """ k/k_0 dispersion for TE surface plasmons
        eps_r1, mu_r1: relative parameters of the upper medium
        eps_r2, mu_r2: relative parameters of the lower medium
        
        calculates alpha_1, alpha_2 and if both are positive it calculates k
    """
    ##TODO: proper and improper reults to be checked
    if eps_r1==eps_r2 and mu_r1==mu_r2:
        mu_1 = mu_r1
        epsilon_1 = eps_r1
        eta_0 = math.sqrt(constants.mu_0/constants.epsilon_0)
        sigma = sigma*eta_0
        k_k0_TE = np.sqrt(mu_1*(4.0*epsilon_1 - mu_1*sigma**2))/2.0    
        k_k0_TE *= (np.imag(sigma)>0.0)
        return k_k0_TE
    elif mu_r1==mu_r2:
        mu_1 = mu_r1
        epsilon_1 = eps_r1
        epsilon_2 = eps_r2
        eta_0 = math.sqrt(constants.mu_0/constants.epsilon_0)
        sigma = sigma*eta_0
        
        alpha_1 = 1j*(-epsilon_1 + epsilon_2 - mu_1*sigma**2)/(2*sigma)
        alpha_2 = 1j*(epsilon_1 - epsilon_2 - mu_1*sigma**2)/(2*sigma)
        
        k_k0_TE = np.sqrt(alpha_1**2 + 1.0)
        k_k0_TE *= (np.real(alpha_1)>0.0)*(np.real(alpha_2)>0.0)
        return k_k0_TE
    else:
        mu_1 = mu_r1
        epsilon_1 = eps_r1
        mu_2 = mu_r2
        epsilon_2 = eps_r2
        eta_0 = math.sqrt(constants.mu_0/constants.epsilon_0)
        sigma = sigma*eta_0
        
        ## -/+ sign before sqrt
        alpha_1 = (mu_1*(1j*mu_2**3*sigma - np.sqrt(mu_2**2*(-epsilon_1*mu_1**3 + epsilon_1*mu_1*mu_2**2 + epsilon_2*mu_1**2*mu_2 - epsilon_2*mu_2**3 - mu_1**2*mu_2**2*sigma**2)))/(mu_2*(mu_1**2 - mu_2**2)))
        alpha_2 = ((-1j*mu_1**2*mu_2*sigma + np.sqrt(-mu_2**2*(mu_1**4*sigma**2 - (mu_1**2 - mu_2**2)*(-epsilon_1*mu_1 + epsilon_2*mu_2 + mu_1**2*sigma**2))))/(mu_1**2 - mu_2**2))
        
        k_k0_TE_1 = np.sqrt(alpha_1**2 + 1.0)
        k_k0_TE_1 *= (np.real(alpha_1)>0.0)*(np.real(alpha_2)>0.0)
        
        ## +/- sign before sqrt
        alpha_1 = (mu_1*(1j*mu_2**3*sigma + np.sqrt(mu_2**2*(-epsilon_1*mu_1**3 + epsilon_1*mu_1*mu_2**2 + epsilon_2*mu_1**2*mu_2 - epsilon_2*mu_2**3 - mu_1**2*mu_2**2*sigma**2)))/(mu_2*(mu_1**2 - mu_2**2)))
        alpha_2 = (-(1j*mu_1**2*mu_2*sigma + np.sqrt(-mu_2**2*(epsilon_1*mu_1**3 - epsilon_1*mu_1*mu_2**2 - epsilon_2*mu_1**2*mu_2 + epsilon_2*mu_2**3 + mu_1**2*mu_2**2*sigma**2)))/(mu_1**2 - mu_2**2))
        
        k_k0_TE_2 = np.sqrt(alpha_1**2 + 1.0)
        k_k0_TE_2 *= (np.real(alpha_1)>0.0)*(np.real(alpha_2)>0.0)

        k_k0_TE = k_k0_TE_1*(k_k0_TE_2==0.0) + k_k0_TE_2*(k_k0_TE_1==0.0)

        return k_k0_TE

        

def plasmonDispersion_TE_normalized(mu_c, omega, tau, T, eps_r1=1.0, mu_r1=1.0, 
        eps_r2=1.0, mu_r2=1.0, cond='Kubo'):
    """
    returns k/k_0
    gives the dispersion of the oplasmons in free standing graphene
    omega: numpy array
    """
    if cond=='Kubo':
        B_0 = 0.0
        sigma_d, sigma_o = condKuboLorentzian(mu_c, B_0, tau, omega, T)
        return plasmonDispersion_TE_normalized_cond(sigma_d, eps_r1=eps_r1, mu_r1=mu_r1, 
            eps_r2=eps_r2, mu_r2=mu_r2)
    elif cond=='Drude':
        B_0 = 0.0
        sigma_d, sigma_o = condDrude(mu_c, B_0, tau, omega, T)
        return plasmonDispersion_TE_normalized_cond(sigma_d, eps_r1=eps_r1, mu_r1=mu_r1, 
            eps_r2=eps_r2, mu_r2=mu_r2)
    else:
        raise ValueError('cond: either "Kubo" or "Drude"')
    return


        

class InfPdicFrStGrapheneMP:
    """
    infinite periodic free standing graphene magnetoplasmons
    using Fourier Bloch analysis
    """    
    
    def __init__(self, N=16, vb=False, solver='lm', symm=False,
            eps_r=1.0, mu_r=1.0, handle_overflow=True):
        ## function generating a periodic conductivity 
        ## [sig_d, sig_o] = f_cond_per(z, f), 
        self.f_cond_per = None  
        self.Z_period = None
        self.N = N
        self.symmetric = symm   ## True:-N...N   False:-N...N-1
        self.freq = None
        self.vbose = vb
        self.solver = solver
        self.eps_r = eps_r
        self.mu_r = mu_r
        self.handle_overflow = handle_overflow
        return
        
    def set_f_step_B(self, mu_c, tau, T, B_0, B_1, dz_0, dz_1):
        self.Z_period = dz_0 + dz_1
        self.param_f_step_B = [mu_c, tau, T, B_0, B_1, dz_0, dz_1]
        self.f_cond_per = self.f_step_B
        return
        
    def set_f_step_mu_c(self, mu_c0, mu_c1, tau, T, B_0, dz_0, dz_1):
        self.Z_period = dz_0 + dz_1
        self.param_f_step_mu_c = [mu_c0, mu_c1, tau, T, B_0, dz_0, dz_1]
        self.f_cond_per = self.f_step_mu_c
        return
        
    def set_f_step_mu_c_tau(self, mu_c0, mu_c1, tau0, tau1, T, B_0, dz_0, dz_1):
        self.Z_period = dz_0 + dz_1
        self.param_f_step_mu_c_tau = [mu_c0, mu_c1, tau0, tau1, T, B_0, dz_0, dz_1]
        self.f_cond_per = self.f_step_mu_c_tau
        return

    def set_f_general_B(self, mu_c, tau, T, B_0, dz, fun):
        """ f_general_B = sigma_0 * f(z/Z_period)
        """
        self.Z_period = dz
        self.param_f_general_B = [mu_c, tau, T, B_0, dz, fun]
        self.f_cond_per = self.f_general_B
        return

    def f_step_B(self, z, f):
        mu_c, tau, T, B_0, B_1, dz_0, dz_1 = self.param_f_step_B
        omega_0 = np.array([2.0*math.pi*f])
        sigma_d_0, sigma_o_0 = condKuboLorentzian(mu_c, B_0, tau, omega_0, T)
        sigma_d_1, sigma_o_1 = condKuboLorentzian(mu_c, B_1, tau, omega_0, T)
        sig_d = (z<=dz_0)*(sigma_d_0[0]) + (z>dz_0)*(sigma_d_1[0])
        sig_o = (z<=dz_0)*(sigma_o_0[0]) + (z>dz_0)*(sigma_o_1[0])
        return [sig_d, sig_o]
    
    def f_step_mu_c(self, z, f):
        mu_c0, mu_c1, tau, T, B_0, dz_0, dz_1 = self.param_f_step_mu_c
        omega_0 = np.array([2.0*math.pi*f])
        sigma_d_0, sigma_o_0 = condKuboLorentzian(mu_c0, B_0, tau, omega_0, T)
        sigma_d_1, sigma_o_1 = condKuboLorentzian(mu_c1, B_0, tau, omega_0, T)
        sig_d = (z<=dz_0)*(sigma_d_0[0]) + (z>dz_0)*(sigma_d_1[0])
        sig_o = (z<=dz_0)*(sigma_o_0[0]) + (z>dz_0)*(sigma_o_1[0])
        return [sig_d, sig_o]
        
    def f_step_mu_c_tau(self, z, f):
        mu_c0, mu_c1, tau0, tau1, T, B_0, dz_0, dz_1 = self.param_f_step_mu_c_tau
        omega_0 = np.array([2.0*math.pi*f])
        sigma_d_0, sigma_o_0 = condKuboLorentzian(mu_c0, B_0, tau0, omega_0, T)
        sigma_d_1, sigma_o_1 = condKuboLorentzian(mu_c1, B_0, tau1, omega_0, T)
        sig_d = (z<=dz_0)*(sigma_d_0[0]) + (z>dz_0)*(sigma_d_1[0])
        sig_o = (z<=dz_0)*(sigma_o_0[0]) + (z>dz_0)*(sigma_o_1[0])
        return [sig_d, sig_o]

    def f_general_B(self, z, f):
        mu_c, tau, T, B_0, dz, fun = self.param_f_general_B
        omega_0 = np.array([2.0*math.pi*f])
        sigma_d_0, sigma_o_0 = condKuboLorentzian(mu_c, B_0, tau, omega_0, T)
        sig_d = sigma_d_0 * fun(z/self.Z_period)
        sig_o = sigma_o_0 * fun(z/self.Z_period)
        return [sig_d, sig_o]
    

    def getCondFCs(self):
        """ get conductivity Fourier Coeffs 
        calculates and returns the Fourier Coefficients of the conductivity
        at frequency f
        """
        f = self.freq
        def f_sigma_d_np(z):
            [sig_d, sig_o] = self.f_cond_per(z, f)
            return sig_d
        sigma_d_tilde_vec = d1_getFourierCoeffs(f_sigma_d_np, 0.0, 
                            self.Z_period, self.N)
        def f_sigma_o_np(z):
            [sig_d, sig_o] = self.f_cond_per(z, f)
            return sig_o
        sigma_o_tilde_vec = d1_getFourierCoeffs(f_sigma_o_np, 0.0, 
                            self.Z_period, self.N)
        return [sigma_d_tilde_vec, sigma_o_tilde_vec]

    
    def getMatrices(self, infPdicGraphMP_sym):
        eq_list_bc, eq_list_ak = infPdicGraphMP_sym.subsParams(self.freq, 
            self.Z_period, eps_r=self.eps_r, mu_r=self.mu_r)
        
        vars_bc, pars_bc, vars_ak, pars_ak = infPdicGraphMP_sym.getVarsPars()
        eig_vars = infPdicGraphMP_sym.getEigVars()
        
        pars_vecs_bc = self.getCondFCs()    ## Fourier Coeffs
        pars_vecs_ak = []

        if self.vbose:
            sigma_d_tilde_vec, sigma_o_tilde_vec = pars_vecs_bc
            print('sigma_d_tilde_vec: \n', sigma_d_tilde_vec)
            print('sigma_o_tilde_vec: \n', sigma_o_tilde_vec)

        n_str = infPdicGraphMP_sym.n_str
        self.A_eq_bc = d1_orthogonalToNumpyMatrix(eq_list_bc, n_str, self.N, vars_bc, 
                pars_bc, pars_vecs_bc, eig_vars)

        self.A_eq_ak = d1_orthogonalToNumpyMatrix(eq_list_ak, n_str, self.N, vars_ak, 
                pars_ak, pars_vecs_ak, eig_vars)
        return
        
    def solveDeterminant(self, infPdicGraphMP_sym, eigvar_vals_0, n_roots=1):
        """ n_roots : number of roots to look for (it tries looking for new 
        roots a maximum number of n_roots times)
        """
        if n_roots == 1:
            k_0 = 2.0*math.pi*self.freq*math.sqrt(constants.mu_0*constants.epsilon_0)
            if self.vbose:
                print('k_0 (free space): ', k_0)
                print('lambda_0:', 2.0*math.pi/k_0)
                print('k_nz/k_0: (n=1)', (2.0*math.pi/self.Z_period)/k_0)
                
            eig_vars = infPdicGraphMP_sym.getEigVars()
            A_eqs_list = [self.A_eq_bc, self.A_eq_ak]
            res = d1_solveDeterminant(eig_vars, eigvar_vals_0, A_eqs_list, 
                    solver=self.solver, handle_overflow=self.handle_overflow)

            if self.vbose:
                print('res:', res)
                print('k/k_0 = ', res[0][0]/k_0)
                print('alpha/k_0 =  ', res[0][1]/k_0, '   sqrt(k/k_0**2 -1): ', \
                        cmath.sqrt((res[0][0]/k_0)**2 - 1.0))
                print('det = ', res[1][0])

            if res[2].success==False:
                return [False, res[2].message]
                                
            det_norm = np.linalg.norm(res[1][0])
            if det_norm>1.0e-3:
                message = "No convergence, det = "+str(det_norm)
                state = [False, message]
                print(message)
                return [False, message]

            eigvar_vals = res[0]
            A_tot_bc = d1_getMatrix(eig_vars, eigvar_vals, self.A_eq_bc)

            null_rtol = 1.0e-7
            rank, x_null, sing_vals = null(A_tot_bc, rtol=null_rtol)

            if self.vbose:
                print('null space relative smallness : ', null_rtol)
                print('sing_vals : ', sing_vals)
                print('rank: ', rank)
                print('x_null-shape: ', x_null.shape)
                print('null-space : \n', x_null)

            return [[eigvar_vals, x_null]]
        else:
            if n_roots<1:
                raise ValueError('n_roots should be a positive integer')
            eig_results = []
            state = [True]
            for n_r in range(n_roots):
                k_0 = 2.0*math.pi*self.freq*math.sqrt(constants.mu_0*constants.epsilon_0)
                if self.vbose:
                    print('k_0 (free space): ', k_0)
                    print('lambda_0:', 2.0*math.pi/k_0)
                    print('k_nz/k_0: (n=1)', (2.0*math.pi/self.Z_period)/k_0)
                    
                eig_vars = infPdicGraphMP_sym.getEigVars()
                A_eqs_list = [self.A_eq_bc, self.A_eq_ak]
                roots_ = [eig_results[i][0] for i in range(len(eig_results))]
                res = d1_solveDeterminant(eig_vars, eigvar_vals_0, A_eqs_list, 
                        solver=self.solver, handle_overflow=self.handle_overflow, 
                        roots_prev = roots_)

                if self.vbose:
                    print('res:', res)
                    print('k/k_0 = ', res[0][0]/k_0)
                    print('alpha/k_0 =  ', res[0][1]/k_0, '   sqrt(k/k_0**2 -1): ', \
                            cmath.sqrt((res[0][0]/k_0)**2 - 1.0))
                    print('det = ', res[1][0])

                if res[2].success==False:
                    state = [False, res[2].message]
                    break
                    
                det_norm = np.linalg.norm(res[1][0])
                if det_norm>1.0e-4:
                    message = "No convergence, det = "+str(det_norm)
                    state = [False, message]
                    print(message)
                    break
                
                eigvar_vals = res[0]
                A_tot_bc = d1_getMatrix(eig_vars, eigvar_vals, self.A_eq_bc)

                null_rtol = 1.0e-7      ## TODO: make it a data member
                rank, x_null, sing_vals = null(A_tot_bc, rtol=null_rtol)

                if self.vbose:
                    print('null space relative smallness : ', null_rtol)
                    print('sing_vals : ', sing_vals)
                    print('rank: ', rank)
                    print('null-space : \n', x_null)
                eig_results.append([eigvar_vals, x_null])
            if len(eig_results)>0:
                return eig_results
            else:
                return state
        
    def sweepFrequency(self, infPdicGraphMP_sym, eigvar_vals_0, x_null_0, f_start, f_end, f_step, record_data=False):
        freq_sweep = ParamSweeper()
        freq_sweep.setLimitPoints(f_start, f_end)
        freq_sweep.setMaxStepSize(f_step)
        freq_sweep.addNext()
        eigvar_vals_start = eigvar_vals_0
        freq_sweep.setLastY(eigvar_vals_start)
        if self.vbose:
            freq_sweep.printXY()
        
        x_null_prev = x_null_0[:,0]
        x_null_prev /= np.linalg.norm(x_null_prev)
        
        x_null_sweep = []
        if record_data:
            x_null_sweep.append(x_null_prev)
            
        n_diverge = 0
        n_diverg_stat = []
        while(freq_sweep.addNext()):
            self.freq = freq_sweep.getLastX()
            k_init, alpha_init = freq_sweep.extrapolateLast()
            eigvar_vals_0 = [k_init, alpha_init]
            self.getCondFCs()
            self.getMatrices(infPdicGraphMP_sym)
            solu_res = self.solveDeterminant(infPdicGraphMP_sym, eigvar_vals_0)
            eigvar_vals, x_null = None, None
            if solu_res[0]!=False:
                eigvar_vals, x_null = solu_res[0]
            else:
                eigvar_vals, x_null = solu_res
            
            if eigvar_vals!=False:
                print('frequency : ', self.freq, '   k_z/k_0: ', eigvar_vals[0].real/\
            (2.0*math.pi*self.freq*math.sqrt(constants.mu_0*constants.epsilon_0)),\
            '  lambda_z/z_period: ', (2.0*math.pi/eigvar_vals[0].real)/self.Z_period)

            if eigvar_vals==False or x_null==[] or x_null.shape[1]==0:
                freq_sweep.refineLast()
                n_diverge += 1
                if n_diverge>3:
                    freq_sweep.decreaseMaxStepSize()
                    n_diverge = 0
                if self.vbose:
                    if eigvar_vals==False:
                        print("eigvar_vals==False")
                    if x_null==[] or x_null.shape[1]==0:
                        print("x_null==[] or x_null.shape[1]==0")
                    print('making refinement...')
                    print('=-'*50)
                print('making refinement...')
                continue
            
            x_null = x_null[:,0]
            x_null /= np.linalg.norm(x_null)
            

            x_correlation = x_null.dot(np.conj(x_null_prev))
            
            if self.vbose:
                print('x_correlation: ', x_correlation, '    abs(x_correlation): ', abs(x_correlation))
            k_correlation_r = np.linalg.norm(np.real(np.array(eigvar_vals) - np.array(eigvar_vals_0)))/\
                max(np.linalg.norm(np.real(np.array(eigvar_vals))), np.linalg.norm(np.real(np.array(eigvar_vals_0))))
            k_correlation_i = np.linalg.norm(np.imag(np.array(eigvar_vals) - np.array(eigvar_vals_0)))/\
                max(np.linalg.norm(np.imag(np.array(eigvar_vals))), np.linalg.norm(np.imag(np.array(eigvar_vals_0))))
            if self.vbose:
                print('k_correlation_r: ', k_correlation_r, '    k_correlation_i: ', k_correlation_i)
            if np.linalg.norm(np.real(np.array(eigvar_vals)))/np.linalg.norm(np.imag(np.array(eigvar_vals)))>1.0e4:
                k_correlation_i = 0.0
                if self.vbose:
                    print('**** k_correlation_i set to zero ****')
            
            n_diverg_stat.append(n_diverge)
            if abs(x_correlation)<0.9 or abs(k_correlation_r)>0.1 or abs(k_correlation_i)>0.1:
                freq_sweep.refineLast()
                n_diverge += 1
                if n_diverge>3:
                    freq_sweep.decreaseMaxStepSize()
                    n_diverge = 0
                if self.vbose:
                    print('making refinement...')
                    print('=-'*50)
                print("abs(x_correlation)=", abs(x_correlation), "abs(k_correlation_r)=", abs(k_correlation_r)\
                     , "abs(k_correlation_i)=", abs(k_correlation_i))
                print('making refinement...')
                continue
            
            x_null_prev = x_null
            freq_sweep.setLastY(eigvar_vals)
            if self.vbose:
                freq_sweep.printXY()
                print('--'*50)

            if record_data:
                x_null_sweep.append(x_null)

        print('divergence stats: ', n_diverg_stat)
        return [freq_sweep, x_null_sweep]
        
    def GetFieldsT_Exz(self, x_null):
        N_p = 0 #number of Fourier components
        if self.symmetric:
            N_p = 2*self.N + 1
        else:
            N_p = 2*self.N
        Ex_field_inv__ = x_null[range(0, N_p)]
        Ez_field_inv__ = x_null[range(N_p, 2*N_p)]
        Ex_field_inv = np.array([Ex_field_inv__[j,0] for j in range(N_p)])
        Ez_field_inv = np.array([Ez_field_inv__[j,0] for j in range(N_p)])
        Ex_field = d1_getInverseFourierCoeffs(Ex_field_inv)
        Ez_field = d1_getInverseFourierCoeffs(Ez_field_inv)
        z_axis = np.linspace(0.0, self.Z_period, N_p, endpoint=False)
        return [Ex_field, Ez_field, z_axis]

    def GetFieldsFourierCoeffsT_Exz(self, x_null):
        N_p = 0 #number of Fourier components
        if self.symmetric:
            N_p = 2*self.N + 1
        else:
            N_p = 2*self.N
        Ex_field_inv__ = x_null[range(0, N_p)]
        Ez_field_inv__ = x_null[range(N_p, 2*N_p)]
        Ex_field_inv = np.array([Ex_field_inv__[j,0] for j in range(N_p)])
        Ez_field_inv = np.array([Ez_field_inv__[j,0] for j in range(N_p)])
        z_axis = np.array(range(N_p))
        z_axis = z_axis - self.N
        return [Ex_field_inv, Ez_field_inv, z_axis]
        
    def GetFieldHarmonics_Exyz(self, infPdicGraphMP_sym, freq, k, alpha, x_null):
        Ex_coeffs, Ez_coeffs, z_axis = self.GetFieldsFourierCoeffsT_Exz(x_null)
        E1_sym, E2_sym = infPdicGraphMP_sym.getElectricFieldHarmonics_subsParams(freq, self.Z_period, k, alpha, 
            eps_r=self.eps_r, mu_r=self.mu_r)
        NF0, NF1 = -self.N, self.N
        if self.symmetric:
            NF1 = self.N + 1
        def f_Ex_n(n):
            if NF0<=n<NF1:
                return Ex_coeffs[n-NF0]
            else:
                return 0.0
        def f_Ez_n(n):
            if NF0<=n<NF1:
                return Ez_coeffs[n-NF0]
            else:
                return 0.0
        Ex1_tilde_new_name = 'tildeEx1'
        Ez1_tilde_new_name = 'tildeEz1'
        E1_sym = sympy.Matrix([[symExp_changeFunctionName(E1_sym[i], infPdicGraphMP_sym.Ex1_tilde, Ex1_tilde_new_name)\
         for i in range(E1_sym.cols)]])
        E1_sym = sympy.Matrix([[symExp_changeFunctionName(E1_sym[i], infPdicGraphMP_sym.Ez1_tilde, Ez1_tilde_new_name)\
         for i in range(E1_sym.cols)]])
        E2_sym = sympy.Matrix([[symExp_changeFunctionName(E2_sym[i], infPdicGraphMP_sym.Ex1_tilde, Ex1_tilde_new_name)\
         for i in range(E2_sym.cols)]])
        E2_sym = sympy.Matrix([[symExp_changeFunctionName(E2_sym[i], infPdicGraphMP_sym.Ez1_tilde, Ez1_tilde_new_name)\
         for i in range(E2_sym.cols)]])
        myfuncs = {Ex1_tilde_new_name:f_Ex_n, Ez1_tilde_new_name:f_Ez_n, 'ImmutableMatrix':np.array, 'I': 1j}
        n = Symbol(infPdicGraphMP_sym.n_str)
        E1_fun = lambdify(n, E1_sym, myfuncs)
        E2_fun = lambdify(n, E2_sym, myfuncs)
        nF = NF1 - NF0
        assert nF==len(Ex_coeffs)
        E1_coeffs = np.zeros((3,nF), dtype=complex)
        E2_coeffs = np.zeros((3,nF), dtype=complex)
        for i in range(NF0, NF1):
            E1_i = E1_fun(i)
            E1_coeffs[0, i-NF0] = E1_i[0][0]
            E1_coeffs[1, i-NF0] = E1_i[0][1]
            E1_coeffs[2, i-NF0] = E1_i[0][2]
            E2_i = E2_fun(i)
            E2_coeffs[0, i-NF0] = E2_i[0][0]
            E2_coeffs[1, i-NF0] = E2_i[0][1]
            E2_coeffs[2, i-NF0] = E2_i[0][2]
        return [E1_coeffs, E2_coeffs]
    
    def GetFieldHarmonics_Hxyz(self, infPdicGraphMP_sym, freq, k, alpha, x_null):
        Ex_coeffs, Ez_coeffs, z_axis = self.GetFieldsFourierCoeffsT_Exz(x_null)
        H1_sym, H2_sym = infPdicGraphMP_sym.getMagneticFieldHarmonics_subsParams(freq, self.Z_period, k, alpha, 
            eps_r=self.eps_r, mu_r=self.mu_r)
        NF0, NF1 = -self.N, self.N
        if self.symmetric: 
            NF1 = self.N + 1
        def f_Ex_n(n):
            if NF0<=n<NF1:
                return Ex_coeffs[n-NF0]
            else:
                return 0.0
        def f_Ez_n(n):
            if NF0<=n<NF1:
                return Ez_coeffs[n-NF0]
            else:
                return 0.0
        Ex1_tilde_new_name = 'tildeEx1'
        Ez1_tilde_new_name = 'tildeEz1'
        H1_sym = sympy.Matrix([[symExp_changeFunctionName(H1_sym[i], infPdicGraphMP_sym.Ex1_tilde, Ex1_tilde_new_name)\
         for i in range(H1_sym.cols)]])
        H1_sym = sympy.Matrix([[symExp_changeFunctionName(H1_sym[i], infPdicGraphMP_sym.Ez1_tilde, Ez1_tilde_new_name)\
         for i in range(H1_sym.cols)]])
        H2_sym = sympy.Matrix([[symExp_changeFunctionName(H2_sym[i], infPdicGraphMP_sym.Ex1_tilde, Ex1_tilde_new_name)\
         for i in range(H2_sym.cols)]])
        H2_sym = sympy.Matrix([[symExp_changeFunctionName(H2_sym[i], infPdicGraphMP_sym.Ez1_tilde, Ez1_tilde_new_name)\
         for i in range(H2_sym.cols)]])
        myfuncs = {Ex1_tilde_new_name:f_Ex_n, Ez1_tilde_new_name:f_Ez_n, 'ImmutableMatrix':np.array, 'I': 1j}
        n = Symbol(infPdicGraphMP_sym.n_str)
        H1_fun = lambdify(n, H1_sym, myfuncs)
        H2_fun = lambdify(n, H2_sym, myfuncs)
        nF = NF1 - NF0
        assert nF==len(Ex_coeffs)
        H1_coeffs = np.zeros((3,nF), dtype=complex)
        H2_coeffs = np.zeros((3,nF), dtype=complex)
        for i in range(NF0, NF1):
            H1_i = H1_fun(i)
            H1_coeffs[0, i-NF0] = H1_i[0][0]
            H1_coeffs[1, i-NF0] = H1_i[0][1]
            H1_coeffs[2, i-NF0] = H1_i[0][2]
            H2_i = H2_fun(i)
            H2_coeffs[0, i-NF0] = H2_i[0][0]
            H2_coeffs[1, i-NF0] = H2_i[0][1]
            H2_coeffs[2, i-NF0] = H2_i[0][2]
        return [H1_coeffs, H2_coeffs]
        
        
    def GetFieldsT_Hxz(self, infPdicGraphMP_sym, freq, k, alpha, x_null):
        H1_coeffs, H2_coeffs = self.GetFieldHarmonics_Hxyz(infPdicGraphMP_sym, freq,
                 k, alpha, x_null)
        N_p = 0 #number of Fourier components
        if self.symmetric:
            N_p = 2*self.N + 1
        else:
            N_p = 2*self.N
        Hx1_inv = H1_coeffs[0,:]
        Hz1_inv = H1_coeffs[2,:]
        Hx2_inv = H2_coeffs[0,:]
        Hz2_inv = H2_coeffs[2,:]
        Hx1_field = d1_getInverseFourierCoeffs(Hx1_inv)
        Hz1_field = d1_getInverseFourierCoeffs(Hz1_inv)
        Hx2_field = d1_getInverseFourierCoeffs(Hx2_inv)
        Hz2_field = d1_getInverseFourierCoeffs(Hz2_inv)
        z_axis = np.linspace(0.0, self.Z_period, N_p, endpoint=False)
        return [Hx1_field, Hz1_field, Hx2_field, Hz2_field, z_axis]
        
    def GetFieldsFourierCoeffsT_Hxz(self, infPdicGraphMP_sym, freq, k, alpha, x_null):
        N_p = 0 #number of Fourier components
        if self.symmetric:
            N_p = 2*self.N + 1
        else:
            N_p = 2*self.N
        H1_coeffs, H2_coeffs = self.GetFieldHarmonics_Hxyz(infPdicGraphMP_sym, freq,
                 k, alpha, x_null)
        Hx1_inv = H1_coeffs[0,:]
        Hz1_inv = H1_coeffs[2,:]
        Hx2_inv = H2_coeffs[0,:]
        Hz2_inv = H2_coeffs[2,:]
        z_axis = np.array(range(N_p))
        z_axis = z_axis - self.N
        return [Hx1_inv, Hz1_inv, Hx2_inv, Hz2_inv, z_axis]
        
    def GetCurrentT_Jxz(self, freq, x_null):
        Ex_field, Ez_field, z_axis = self.GetFieldsT_Exz(x_null)
        N_p = len(z_axis)
        Jx = np.zeros(N_p, dtype=complex) 
        Jz = np.zeros(N_p, dtype=complex) 
        for i in range(N_p):
            z = z_axis[i]
            sig_d, sig_o = self.f_cond_per(z, freq)
            ## TODO: warning: +/- signs depend on convention used in infPdicGraphMP_sym
            Jx[i] = sig_d*Ex_field[i] - sig_o*Ez_field[i]   
            Jz[i] = sig_d*Ez_field[i] + sig_o*Ex_field[i]
        return [Jx, Jz, z_axis]
        
    def GetFieldsFourierCoeffsT_Jxz(self, freq, x_null):
        Jx, Jz, _ = self.GetCurrentT_Jxz(freq, x_null)
        N_p = len(Jz)
        Jx_inv = fftpack.fftshift(fftpack.fft(Jx))/N_p
        Jz_inv = fftpack.fftshift(fftpack.fft(Jz))/N_p
        z_axis = np.array(range(N_p))
        z_axis = z_axis - self.N
        return [Jx_inv, Jz_inv, z_axis]
        
        
    def GetCurrentHarmonics_Jxyz(self, infPdicGraphMP_sym, freq, x_null):
        ## uses direct convolutions : J_bar = conv(sigma_bar, E_bar)
        Ex_coeffs, Ez_coeffs, z_axis = self.GetFieldsFourierCoeffsT_Exz(x_null)
        self.freq = freq
        sigma_d_coeffs, sigma_o_coeffs = self.getCondFCs()
        J_sym = infPdicGraphMP_sym.getCurrentHarmonics()
        NF0, NF1 = -self.N, self.N
        if self.symmetric: 
            NF1 = self.N + 1
        def f_Ex_n(n):
            if NF0<=n<NF1:
                return Ex_coeffs[n-NF0]
            else:
                return 0.0
        def f_Ez_n(n):
            if NF0<=n<NF1:
                return Ez_coeffs[n-NF0]
            else:
                return 0.0
        def f_Sig_d_n(n):
            if NF0<=n<NF1:
                return sigma_d_coeffs[n-NF0]
            else:
                return 0.0
        def f_Sig_o_n(n):
            if NF0<=n<NF1:
                return sigma_o_coeffs[n-NF0]
            else:
                return 0.0
        Ex1_tilde_new_name = 'tildeEx1'
        Ez1_tilde_new_name = 'tildeEz1'
        sigma_d_tilde_new_name = 'tildeSigma_d'
        sigma_o_tilde_new_name = 'tildeSigma_o'
        J_sym = sympy.Matrix([[symExp_changeFunctionName(J_sym[i], infPdicGraphMP_sym.Ex1_tilde, Ex1_tilde_new_name)\
         for i in range(J_sym.cols)]])
        J_sym = sympy.Matrix([[symExp_changeFunctionName(J_sym[i], infPdicGraphMP_sym.Ez1_tilde, Ez1_tilde_new_name)\
         for i in range(J_sym.cols)]])
        J_sym = sympy.Matrix([[symExp_changeFunctionName(J_sym[i], infPdicGraphMP_sym.sigma_d_tilde, sigma_d_tilde_new_name)\
         for i in range(J_sym.cols)]])
        J_sym = sympy.Matrix([[symExp_changeFunctionName(J_sym[i], infPdicGraphMP_sym.sigma_o_tilde, sigma_o_tilde_new_name)\
         for i in range(J_sym.cols)]])
        myfuncs = {Ex1_tilde_new_name:f_Ex_n, Ez1_tilde_new_name:f_Ez_n, 
        sigma_d_tilde_new_name:f_Sig_d_n, sigma_o_tilde_new_name:f_Sig_o_n, 'ImmutableMatrix':np.array, 'I': 1j}
        n = Symbol(infPdicGraphMP_sym.n_str)
        n1 = Symbol(infPdicGraphMP_sym.n_str+'_1')
        J_fun = lambdify((n, n1), J_sym, myfuncs)
        nF = NF1 - NF0
        assert nF==len(Ex_coeffs)
        J_coeffs = np.zeros((3,nF), dtype=complex)
        for n_ in range(NF0, NF1):
            for n1_ in range(NF0, NF1):
                J_nn1 = J_fun(n_, n1_)
                J_coeffs[0, n_-NF0] += J_nn1[0][0]
                J_coeffs[1, n_-NF0] += J_nn1[0][1]
                J_coeffs[2, n_-NF0] += J_nn1[0][2]
        return J_coeffs
        
    def GetCurrentT_Jxz_Conv(self, infPdicGraphMP_sym, freq, x_null):
        ## uses convolution and inverse Fourier transform to calculate current
        J_xyz_inv = self.GetCurrentHarmonics_Jxyz(infPdicGraphMP_sym, freq, x_null)
        Jx = d1_getInverseFourierCoeffs(J_xyz_inv[0,:])
        Jz = d1_getInverseFourierCoeffs(J_xyz_inv[2,:])
        N_p = len(Jz)
        z_axis = np.linspace(0.0, self.Z_period, N_p, endpoint=False)
        return [Jx, Jz, z_axis]
        
        
    
##--------------   Graphene sheet on Ferrite substrate ---    
    
class GrapheneOnFerriteSubstrate:
    def __init__(self, ferrite):
        self.ferrite = ferrite
        self.vbose = False
        return
        
    def setFrequency(self, f):
        self.freq = f
        return
    
    def setTopMaterialParameters(self, eps_r1, mu_r1):
        self.eps_r1 = eps_r1
        self.mu_r1 = mu_r1
        return

    def setButtomMaterialParameters(self, eps_r2):
        self.eps_r2 = eps_r2
        return

    def setGrapheneSheetParameters(self, mu_c, tau, T):
        self.mu_c = mu_c
        self.tau = tau
        self.T = T
        return
        
    def findTEmodeOverFrequencyRange(self, omega_0, omega_1, N):
        """ calculates the conductivity of graphene over the given frequency range
            and gives a range of frequencies where imag(sigma)>0.0
        """
        omega = np.linspace(omega_0, omega_1, N)
        B_0 = 0.0
        sigma, _ = condKuboLorentzian(self.mu_c, B_0, self.tau, omega, self.T)
        omega = omega*(np.imag(sigma)>0.0)
        omega = np.ma.masked_equal(omega, 0.0).compressed()
        if self.vbose:
            mu_rd2, mu_ro2 = self.ferrite.getRelativePermeabilities(omega)
            sigma, _ = condKuboLorentzian(self.mu_c, 0.0, self.tau, omega, self.T)
            for i in range(len(omega)):
                print(i, '  omega: ', omega[i], '  sig_i: ', sigma[i], '\n'+' '*10, '  mu_rd: ', mu_rd2[i], '  mu_ro: ', mu_ro2[i])
        return omega
        
    def func_dispersion(self, K):
        ## function to solve for finding dispersion curves
        ## K = k/k_0
        omega = 2.0*math.pi*self.freq
        B_0 = 0.0
        sigma, _ = condKuboLorentzian(self.mu_c, B_0, self.tau, np.array([omega]), self.T)
        sigma = sigma[0]
        self.mu_rd2, self.mu_ro2 = self.ferrite.getRelativePermeabilities(omega)
        #print('mu_rd2: ', self.mu_rd2, '  mu_ro2: ', self.mu_ro2)
        mu_0, epsilon_0 = constants.mu_0, constants.epsilon_0
        #eq = (-mu_0*self.mu_r1*sigma*(self.mu_rd2**2 + self.mu_ro2**2) + sqrt(epsilon_0*mu_0)*(self.mu_r1*(K*self.mu_ro2 + 1j*self.mu_rd2*sqrt((-self.eps_r2*self.mu_ro2**2 + self.mu_rd2*(K**2 - self.eps_r2*self.mu_rd2))/self.mu_rd2)) + 1j*sqrt(K**2 - self.eps_r1*self.mu_r1)*(self.mu_rd2**2 + self.mu_ro2**2))/(mu_0*self.mu_r1*(self.mu_rd2**2 + self.mu_ro2**2)))
        
        #better conditionned
        eq = -sigma + sqrt(epsilon_0/mu_0)*(self.mu_r1*(K*self.mu_ro2 + 1j*self.mu_rd2*sqrt((-self.eps_r2*self.mu_ro2**2 + self.mu_rd2*(K**2 - self.eps_r2*self.mu_rd2))/self.mu_rd2)) + 1j*sqrt(K**2 - self.eps_r1*self.mu_r1)*(self.mu_rd2**2 + self.mu_ro2**2))/(self.mu_r1*(self.mu_rd2**2 + self.mu_ro2**2))

        return eq
        
    def solveDispersion(self, freq, k_k0_init=None, n_roots=4, n_max_iter=1000, 
            tol_y_abs=1.0e-3, tol_y_rel=None, solver='lm'):
        self.freq = freq
        k_i = 2.0*math.pi*freq*cmath.sqrt(self.mu_r1*self.eps_r1*constants.mu_0*constants.epsilon_0)
        if k_k0_init!=None:
            k_i = k_k0_init
        if solver=='muller':
            res = solveMuller(self.func_dispersion, k_i, 1.01*k_i, 1.02*k_i, 
                n_roots=n_roots, toly_rel=tol_y_rel, tolx_abs=None, toly_abs=tol_y_abs, n_max=n_max_iter)
            return res
        else:
            res = RootsMultipleComplex(self.func_dispersion, k_i, n_roots=n_roots, 
                solver=solver, tol=tol_y_abs, ftol=tol_y_abs, maxfev=n_max_iter, maxiter=n_max_iter)
            return res
        
    def func_dispersion_alpha_k(self, x):
        ## function to solve for finding dispersion parameters alpha_1, alpha_2, k
        ## x = [Re(alpha_1/k0), Im(alpha_1/k0), Re(alpha_2/k0), Im(alpha_2/k0), Re(k/k0), Im(k/k0)]
        omega = 2.0*math.pi*self.freq
        B_0 = 0.0
        sigma, _ = condKuboLorentzian(self.mu_c, B_0, self.tau, np.array([omega]), self.T)
        sigma = sigma[0]
        self.mu_rd2, self.mu_ro2 = self.ferrite.getRelativePermeabilities(omega)
        #print('mu_rd2: ', self.mu_rd2, '  mu_ro2: ', self.mu_ro2)
        mu_0, epsilon_0 = constants.mu_0, constants.epsilon_0
        
        alpha_1 = x[0] + 1j*x[1]
        alpha_2 = x[2] + 1j*x[3]
        k = x[4] + 1j*x[5]
        
        epsilon_r1, mu_r1 = self.eps_r1, self.mu_r1
        epsilon_r2, mu_rd2, mu_ro2 = self.eps_r2, self.mu_rd2, self.mu_ro2
        sigma = sigma*np.sqrt(mu_0/epsilon_0)
        
        eq_1 = (1j*alpha_1*(mu_rd2**2 + mu_ro2**2) - mu_r1*sigma*(mu_rd2**2 + mu_ro2**2) + mu_r1*(1j*alpha_2*mu_rd2 + mu_ro2*k))
        eq_2 = (alpha_1**2 + epsilon_r1*mu_r1 - k**2)
        eq_3 = (alpha_2**2*mu_rd2 + epsilon_r2*mu_rd2**2 + epsilon_r2*mu_ro2**2 - mu_rd2*k**2)
        
        return [np.real(eq_1), np.imag(eq_1), np.real(eq_2), np.imag(eq_2), np.real(eq_3), np.imag(eq_3)]


    def solveDispersion_alpha_k(self, freq, k_k0_init=None, n_roots=4, n_max_iter=1000, 
            tol_y_abs=1.0e-3, tol_y_rel=None, solver='lm'):
        self.freq = freq
        omega = 2.0*math.pi*self.freq
        self.mu_rd2, self.mu_ro2 = self.ferrite.getRelativePermeabilities(omega)

        k_i = np.sqrt((self.mu_r1+self.mu_rd2)/2.0*(self.eps_r1+self.eps_r2)/2.0)
        if k_k0_init!=None:
            k_i = k_k0_init
        alpha_1_i = np.sqrt(k_i**2 - self.eps_r1*self.mu_r1)
        alpha_2_i = np.sqrt(k_i**2 - self.eps_r2*self.mu_rd2 - self.eps_r2*self.mu_ro2**2/self.mu_rd2)
        
        x_0 = np.array([np.real(alpha_1_i), np.imag(alpha_1_i), \
            np.real(alpha_2_i), np.imag(alpha_2_i), np.real(k_i), np.imag(k_i)])
        
        res = newton_krylov(self.func_dispersion_alpha_k, x_0, f_tol=tol_y_abs, maxiter=n_max_iter)
        return res



class GrapheneOnFerriteSubstrateSuperstrate:
    def __init__(self, ferrite):
        self.ferrite = ferrite
        self.vbose = False
        return
        
    def setFrequency(self, f):
        self.freq = f
        return
    
    def setTopMaterialParameters(self, eps_r1, mu_r1):
        self.eps_r1 = eps_r1
        self.mu_r1 = mu_r1
        return

    def setButtomMaterialParameters(self, eps_r2):
        self.eps_r2 = eps_r2
        return

    def setGrapheneSheetParameters(self, mu_c, tau, T):
        self.mu_c = mu_c
        self.tau = tau
        self.T = T
        return
        
    def findTEmodeOverFrequencyRange(self, omega_0, omega_1, N):
        """ calculates the conductivity of graphene over the given frequency range
            and gives a range of frequencies where imag(sigma)>0.0
        """
        omega = np.linspace(omega_0, omega_1, N)
        B_0 = 0.0
        sigma, _ = condKuboLorentzian(self.mu_c, B_0, self.tau, omega, self.T)
        omega = omega*(np.imag(sigma)>0.0)
        omega = np.ma.masked_equal(omega, 0.0).compressed()
        if self.vbose:
            mu_rd2, mu_ro2 = self.ferrite.getRelativePermeabilities(omega)
            sigma, _ = condKuboLorentzian(self.mu_c, 0.0, self.tau, omega, self.T)
            for i in range(len(omega)):
                print(i, '  omega: ', omega[i], '  sig_i: ', sigma[i], '\n'+' '*10, '  mu_rd: ', mu_rd2[i], '  mu_ro: ', mu_ro2[i])
        return omega
                
    def func_dispersion_alpha_k(self, x):
        ## function to solve for finding dispersion parameters alpha_1, alpha_2, k
        ## x = [Re(alpha_1/k0), Im(alpha_1/k0), Re(alpha_2/k0), Im(alpha_2/k0), Re(k/k0), Im(k/k0)]
        omega = 2.0*math.pi*self.freq
        B_0 = 0.0
        sigma, _ = condKuboLorentzian(self.mu_c, B_0, self.tau, np.array([omega]), self.T)
        sigma = sigma[0]
        self.mu_rd2, self.mu_ro2 = self.ferrite.getRelativePermeabilities(omega)
        #print('mu_rd2: ', self.mu_rd2, '  mu_ro2: ', self.mu_ro2)
        mu_0, epsilon_0 = constants.mu_0, constants.epsilon_0
        
        alpha_1 = x[0] + 1j*x[1]
        alpha_2 = x[2] + 1j*x[3]
        k = x[4] + 1j*x[5]
        
        epsilon_r1, mu_rd1, mu_ro1 = self.eps_r1, self.mu_rd2, -self.mu_ro2
        epsilon_r2, mu_rd2, mu_ro2 = self.eps_r2, self.mu_rd2, self.mu_ro2
        sigma = sigma*np.sqrt(mu_0/epsilon_0)
        
        eq_1 = (-sigma*(mu_rd1**2 + mu_ro1**2)*(mu_rd2**2 + mu_ro2**2) + (mu_rd1**2 + mu_ro1**2)*(1j*alpha_2*mu_rd2 + mu_ro2*k) + (mu_rd2**2 + mu_ro2**2)*(1j*alpha_1*mu_rd1 - mu_ro1*k))
        eq_2 = (alpha_1**2*mu_rd1 + epsilon_r1*mu_rd1**2 + epsilon_r1*mu_ro1**2 - mu_rd1*k**2)
        eq_3 = (alpha_2**2*mu_rd2 + epsilon_r2*mu_rd2**2 + epsilon_r2*mu_ro2**2 - mu_rd2*k**2)
        
        return [np.real(eq_1), np.imag(eq_1), np.real(eq_2), np.imag(eq_2), np.real(eq_3), np.imag(eq_3)]


    def solveDispersion_alpha_k(self, freq, k_k0_init=None, n_roots=4, n_max_iter=1000, 
            tol_y_abs=1.0e-3, tol_y_rel=None, solver='lm'):
        self.freq = freq
        omega = 2.0*math.pi*self.freq
        self.mu_rd2, self.mu_ro2 = self.ferrite.getRelativePermeabilities(omega)

        k_i = np.sqrt((self.mu_r1+self.mu_rd2)/2.0*(self.eps_r1+self.eps_r2)/2.0)
        if k_k0_init!=None:
            k_i = k_k0_init
        alpha_1_i = np.sqrt(k_i**2 - self.eps_r1*self.mu_rd2 - self.eps_r1*self.mu_ro2**2/self.mu_rd2)
        alpha_2_i = np.sqrt(k_i**2 - self.eps_r2*self.mu_rd2 - self.eps_r2*self.mu_ro2**2/self.mu_rd2)
        
        x_0 = np.array([np.real(alpha_1_i), np.imag(alpha_1_i), \
            np.real(alpha_2_i), np.imag(alpha_2_i), np.real(k_i), np.imag(k_i)])
        
        res = newton_krylov(self.func_dispersion_alpha_k, x_0, f_tol=tol_y_abs, maxiter=n_max_iter)
        return res







