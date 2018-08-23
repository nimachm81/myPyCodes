## TransElec.py
## Transformation Electromagnetics


__all__ = ["TransElec"]


from Electromagnetics.VectorCalculus import getJacobianMatrix

from Electromagnetics.Misc import SymMatrixdoit


class TransElec:
    
    def __init__(self, vars_old=None, vars_new=None, funcs=None):
        self.vars_old = vars_old
        self.vars_new = vars_new
        self.funcs = funcs
        
    def SetVarsOld(self, vars_old, metric_old):
        """ symbolic variables (list) of the old coordinate system and the metric tensor
        as a Matrix
        """
        self.vars_old = vars_old
        self.metric_old = metric_old
        
    def SetVarsNew(self, vars_new, metric_new):
        """ symbolic variables (list) of the new coordinate system and the metric tensor
        as a Matrix
        """
        self.vars_new = vars_new
        self.metric_new = metric_new
        
    def SetFuncs(self, funcs):
        """ funcs: list of symbolic functions
            vars_old[i] = funcs[i](vars_new[0], vars_new[1]...)
        """
        self.funcs = funcs
        
    def GetJacobianMatrix(self):
        jac = getJacobianMatrix(self.vars_old, self.vars_new, self.funcs)  
        SymMatrixdoit(jac)
        return jac
        
        
        
        
        
        
