## EMEQs_sym.py
## Electromagnetic equations in symbolic form

__all__ = ["MaxwellEqs_r"]

from sympy import Symbol, symbols, Matrix, I, exp
from Electromagnetics.VectorCalculus import curl_r, divergence_r

class MaxwellEqs_r:
    """ Maxwell equations in rectangular coordinate system
    timeConv='+': e^+j*omega*t time convention
    """
    def __init__(self, mediumInd=-1, eps_type='scalar', mu_type='scalar', sigma_type='scalar'):
        """ types: 'scalar' or 'tensor'
        """
        self.eps_type = eps_type
        self.mu_type = mu_type
        self.sigma_type = sigma_type
        Ex_str = 'E_{x}'
        Ey_str = 'E_{y}'
        Ez_str = 'E_{z}'
        Hx_str = 'H_{x}'
        Hy_str = 'H_{y}'
        Hz_str = 'H_{z}'
        k_str = 'k'
        eps_str = '\\epsilon'
        mu_str = '\\mu'
        eta_str = '\\eta'
        sigma_str = '\\sigma'
        sigma_xx_str = '\\sigma_{xx}'
        sigma_xy_str = '\\sigma_{xy}'
        sigma_xz_str = '\\sigma_{xz}'
        sigma_yx_str = '\\sigma_{yx}'
        sigma_yy_str = '\\sigma_{yy}'
        sigma_yz_str = '\\sigma_{yz}'
        sigma_zx_str = '\\sigma_{zx}'
        sigma_zy_str = '\\sigma_{zy}'
        sigma_zz_str = '\\sigma_{zz}'
        if mediumInd>=0:
            Ex_str = 'E_{x' + str(mediumInd) + '}'
            Ey_str = 'E_{y' + str(mediumInd) + '}'
            Ez_str = 'E_{z' + str(mediumInd) + '}'
            Hx_str = 'H_{x' + str(mediumInd) + '}'
            Hy_str = 'H_{y' + str(mediumInd) + '}'
            Hz_str = 'H_{z' + str(mediumInd) + '}'
            k_str = 'k_' + str(mediumInd)
            eps_str = '\\epsilon_' + str(mediumInd)
            mu_str = '\\mu_' + str(mediumInd)
            eta_str = '\\eta_' + str(mediumInd)
            sigma_str = '\\sigma_' + str(mediumInd)
            sigma_xx_str = '\\sigma_{xx' + str(mediumInd) + '}'
            sigma_xy_str = '\\sigma_{xy' + str(mediumInd) + '}'
            sigma_xz_str = '\\sigma_{xz' + str(mediumInd) + '}'
            sigma_yx_str = '\\sigma_{yx' + str(mediumInd) + '}'
            sigma_yy_str = '\\sigma_{yy' + str(mediumInd) + '}'
            sigma_yz_str = '\\sigma_{yz' + str(mediumInd) + '}'
            sigma_zx_str = '\\sigma_{zx' + str(mediumInd) + '}'
            sigma_zy_str = '\\sigma_{zy' + str(mediumInd) + '}'
            sigma_zz_str = '\\sigma_{zz' + str(mediumInd) + '}'
            
            
        self.x, self.y, self.z = symbols('x y z')
        self.Ex = Symbol(Ex_str)
        self.Ey = Symbol(Ey_str)
        self.Ez = Symbol(Ez_str)
        self.Hx = Symbol(Hx_str)
        self.Hy = Symbol(Hy_str)
        self.Hz = Symbol(Hz_str)
        self.k = Symbol(k_str)
        self.omega = Symbol('\\omega')
        self.eps = Symbol(eps_str)
        self.mu = Symbol(mu_str)
        self.eta = Symbol(eta_str)
        self.sigma = Symbol(sigma_str)
        self.sigma_xx = Symbol(sigma_xx_str)
        self.sigma_xy = Symbol(sigma_xy_str)
        self.sigma_xz = Symbol(sigma_xz_str)
        self.sigma_yx = Symbol(sigma_yx_str)
        self.sigma_yy = Symbol(sigma_yy_str)
        self.sigma_yz = Symbol(sigma_yz_str)
        self.sigma_zx = Symbol(sigma_zx_str)
        self.sigma_zy = Symbol(sigma_zy_str)
        self.sigma_zz = Symbol(sigma_zz_str)
        self.sigma_tensor = Matrix([[self.sigma_xx, self.sigma_xy, self.sigma_xz],
                                    [self.sigma_yx, self.sigma_yy, self.sigma_yz],
                                    [self.sigma_zx, self.sigma_zy, self.sigma_zz]])
        return
    
    
    def getMaxwellFaraday(self):
        if self.mu_type=='scalar':
            E = Matrix([[self.Ex, self.Ey, self.Ez]])
            H = Matrix([[self.Hx, self.Hy, self.Hz]])
            lhs = curl_r(E)
            rhs = -I*self.omega*self.mu*H
            return [lhs, rhs]
        else:
            raise NotImplementedError()
        
    def getMaxwellAmper(self, hasSigma=False):
        if self.eps_type=='scalar' and self.sigma_type=='scalar':
            E = Matrix([[self.Ex, self.Ey, self.Ez]])
            H = Matrix([[self.Hx, self.Hy, self.Hz]])
            lhs = curl_r(H)
            rhs = I*self.omega*self.eps*E
            if hasSigma:
                rhs = (I*self.omega*self.eps + self.sigma)*E
            return [lhs, rhs]
        elif self.eps_type=='scalar' and self.sigma_type=='tensor':
            E = Matrix([[self.Ex, self.Ey, self.Ez]])
            H = Matrix([[self.Hx, self.Hy, self.Hz]])
            lhs = curl_r(H)
            rhs = I*self.omega*self.eps*E
            if hasSigma:
                rhs = I*self.omega*self.eps*E + (self.sigma_tensor*E.T).T
            return [lhs, rhs]
        else:
            raise NotImplementedError()
        
        
    def getDivE(self):
        if self.eps_type=='scalar':
            E = Matrix([[self.Ex, self.Ey, self.Ez]])
            lhs = divergence_r(E)
            rhs = 0
            return [lhs, rhs]
        else:
            raise NotImplementedError()
        
    def getDivH(self):
        if self.mu_type=='scalar':
            H = Matrix([[self.Hx, self.Hy, self.Hz]])
            lhs = divergence_r(H)
            rhs = 0
            return [lhs, rhs]
        else:
            raise NotImplementedError()
        
    def getZPropEqs(self, expr):
        """ replaces field components by A*e^-j*k*z
        """
        harm = exp(-I*self.k*self.z)
        expr_z =  expr.subs([(self.Ex, self.Ex(self.x, self.y)*harm), (self.Ey, self.Ey(self.x, self.y)*harm),
            (self.Ez, self.Ez(self.x, self.y)*harm), (self.Hx, self.Hx(self.x, self.y)*harm), 
            (self.Hy, self.Hy(self.x, self.y)*harm), (self.Hz, self.Hz(self.x, self.y)*harm)])
        return expr_z
        
    def ZPropEqs_ReplFunBySym(self, expr):
        """ replaces E_x(x, y) by E_x ...
        """
        expr_z =  expr.subs([(self.Ex(self.x, self.y), self.Ex), (self.Ey(self.x, self.y), self.Ey),
            (self.Ez(self.x, self.y), self.Ez), (self.Hx(self.x, self.y), self.Hx), 
            (self.Hy(self.x, self.y), self.Hy), (self.Hz(self.x, self.y), self.Hz)])
        return expr_z
    
    
    
    
    
    
    
