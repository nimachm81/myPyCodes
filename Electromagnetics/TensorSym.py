## TensorSym.py

""" symbolic tensor calculus
"""


__all__ = ["v__unset___arg", "IndexType", "TensorIndex", "Tensor"]


import numpy as np
from sympy import Symbol

from Electromagnetics.SymExprTree import symExp_FindReplaceFunctionArgs

from enum import Enum
class IndexType(Enum):
    covariant = 1
    contravariant = 2


v__unset___arg = Symbol('.')


class TensorIndex:
    def __init__(self, name=None, type=None, range=None):
        """ name: abitrary name, just for convenience
            typr: contravariant or covariant
            range: 2-tuple
        """
        self.name = name
        self.type = type
        self.range = range
        if range is not None:
            assert range[1]>range[0]
            


class Tensor:
    def __init__(self, label, indices, indexLabels=None):
        self.label = label
        self.indices = indices
        self.indexLabels = indexLabels
        if indexLabels==None:
            self.indexLabels = [indices[i].name for i in range(len(indices))]
        assert len(self.indexLabels)==len(self.indexLabels)
        self.getShape()
        if len(self.shape)>0:
            self.items = np.ndarray(self.shape, dtype=object)
        else:
            self.items = None
        return
                
    def getShape(self):
        order = len(self.indices)
        shape = [None]*order
        for i in range(order):
            range_i = self.indices[i].range
            shape[i] = range_i[1] - range_i[0]
        self.shape = shape
        return shape
        
    def __getitem__(self, key):
        if len(self.shape)>0:
            return self.items[key]
        else:
            assert key is None
            return self.items

    def __setitem__(self, key, value):
        if len(self.shape)>0:
            self.items[key] = value
        else:
            assert key is None
            self.items = value

    def __str__(self):
        out_str = self.label
        order = len(self.indices)
        for i in range(order):
            ind_i = self.indices[i]
            ind_label_i = self.indexLabels[i]
            if ind_i.type==IndexType.covariant:
                out_str += r'\vphantom1_{'+ind_label_i+'}'
            elif ind_i.type==IndexType.contravariant:
                out_str += r'\vphantom1^{'+ind_label_i+'}'
            else:
                raise NotImplementedError()
        return out_str
        
        
    def resetIndexLabel(self, ind, label_new):
        """ ind: tensor index position inside self.indices
            label_new: the new tensor index label
        """
        self.indexLabels[ind] = label_new
        
        
    @staticmethod
    def TensorProd(T1, T2):
        global v__unset___arg
        
        label_1 = T1.label
        indices_1 = T1.indices 
        indexLabels_1 = T1.indexLabels
        shape_1 = T1.shape
        items_1 = T1.items

        label_2 = T2.label
        indices_2 = T2.indices 
        indexLabels_2 = T2.indexLabels
        shape_2 = T2.shape
        items_2 = T2.items

        label = "({}.{})".format(label_1, label_2)
        indices = indices_1 + indices_2
        indexLabels = indexLabels_1 + indexLabels_2
        shape = shape_1 + shape_2

        items = None        
        if len(shape_1)==0 or len(shape_2)==0:
            items = items_1*items_2
        else:
            items = np.ndarray(shape, dtype=object)
            
            inds_1_cycled = np.zeros(len(shape_1), dtype=int)    #cycles through indices_1
            
            def increaseIndsCyclicLeft(inds, shape):
                N = len(shape)
                assert len(inds)==N
                end = True
                for i in range(N):
                    inds[i] += 1
                    if inds[i]==shape[i]:
                        inds[i] = 0
                        end = False
                    if end:
                        break
                if np.all(inds==0):
                    return False
                else:
                    return True
            
            while True:
                f_1 = items_1[tuple(inds_1_cycled)]
                                
                if type(f_1) not in [int, float, complex] and f_1.has(v__unset___arg):
                    ##iterating explicitely
                    inds_2_cycled = np.zeros(len(shape_2), dtype=int)
                    while True:
                        f_2 = items_2[tuple(inds_2_cycled)]
                        inds_all = inds_1_cycled.tolist() + inds_2_cycled.tolist()

                        items[tuple(inds_all)] = symExp_FindReplaceFunctionArgs(f_1, (v__unset___arg,), (f_2,))[0]

                        if not increaseIndsCyclicLeft(inds_2_cycled, shape_2):
                            break
                else:
                    ##use numpy slicing
                    inds_all = inds_1_cycled.tolist() + [None]*len(shape_2)
                    items[tuple(inds_all)] = f_1*items_2

                if not increaseIndsCyclicLeft(inds_1_cycled, shape_1):
                    break
                

            
        T = Tensor(label, indices, indexLabels)
        T.items = items
        return T
        

    def Contract(self, _ind_0, _ind_1):
        ind_0 = min(_ind_0, _ind_1)
        ind_1 = max(_ind_0, _ind_1)
        
        assert ind_0!=ind_1
        assert self.shape[ind_0]==self.shape[ind_1]

        inds_new = list(range(len(self.shape)))
        del inds_new[ind_1]
        del inds_new[ind_0]
        
        
        shape_new = (np.array(self.shape)[inds_new]).tolist()
        indices_new = (np.array(self.indices)[inds_new]).tolist()
        indexLabels_new = (np.array(self.indexLabels)[inds_new]).tolist()
        
        items_new = np.trace(self.items, axis1=ind_0, axis2=ind_1, dtype=object)
        
        self.items = items_new
        self.shape = shape_new
        self.indices = indices_new
        self.indexLabels = indexLabels_new
        
        
    def raiseIndices(self, ind, metric_inv):
        
        return
        
        
    def lowerIndices(self, ind, metric):
    
        return
        
    def transformCoordinates(self, A):
        """ x' = A*x
            To support custom (non-tensor) coordinate transformations such as
            Christophel symbols this class can be extended and this member can
            then be overrriden
        """
        
        return


class TensorField:
    """ components:
            - coordinate variables    
            - metric tensor
            - basis vectors
            - member tenror fields
            
        supports:
            - lowering and raising of indices
            - coordinate transformations
            - expressing member tensor fields in a different coordinate system
        
    """

    
    def __init__(self):
        return        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        



