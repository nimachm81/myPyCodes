## Misc.py


__all__ = ["svd_symbolic", "solve_linear_system_nontrivial",
            "SymMatrixSimplify", "SymMatrixdiff", "SymMatrixdoit", "SymMatrixsubs", 
            "SymFracSimplify",
            "SimplifyDiffEquationCoeffs", "GetDiffEquationCoeffs",
            "solveMuller", "RootsMultipleComplex",
            "lagrange1d_coeffs", "lagrange2d_coeffs", "barycentricInterpolation_coeffs",
            "Point2D", "Point3D", "GetReiprocalVecs",
            "null",
            "tic", "toc",
            "replace_whole_word",
            "NormalizeArrayTo0n1", "GetMaxComplexVec2D",
            "ColorMaps", "GetColorMap",
            "Fourier1D", "InvFourier1D", "Fourier2D", "InvFourier2D",
            "ParamSweeper",
            "testReload",
            "TextInputParameterReader"]



##---

def testReload():
    print('6')
    return


###--------   Some Matrix Operations
import sympy
from sympy import Matrix, eye, zeros, Symbol, diff, fraction, simplify,\
                  Derivative, Integer
import operator

def svd_symbolic(A):
    """
    decomposition U.inv()*S*V.inv() = A or U*A*V = S
    U and V are not unitary (to be checked)
    """
    if not isinstance(A, Matrix):
        raise TypeError("A should be a sympy Matrix")
    M, N = A.rows, A.cols
    S, U, V = A, eye(M), eye(N)
    #lower triangular part
    for i in range(0,min(M,N)):
        if S[i,i]==0:
            for j in range(i+1,M):
                if S[j,i]!=0:
                    S.row_swap(i,j)
                    U.row_swap(i,j)
                    break
        for j in range(i+1,M):
            if S[j,i]!=0:
                P = eye(M)
                P[j,i] = S[j,i]
                P[j,j] = -S[i,i]
                U = P*U
                S = P*S
                #S[j*N] = S[j,i]*S.row(i) - S[i,i]*S.row(j) #TODO should replace S=P*S
    #upper triangular part
    for i in range(min(M,N)-1,0,-1):
        if S[i,i]!=0:
            for j in range(0,i):
                if S[j,i]!=0:
                    P = eye(M)
                    P[j,i] = S[j,i]
                    P[j,j] = -S[i,i]
                    U = P*U
                    #S = P*S
                    S[j*N] = S[j,i]*S.row(i) - S[i,i]*S.row(j)
    #right triangular part
    for i in range(min(M,N)-1,-1,-1):
        if S[i,i]==0:
            for j in range(i+1,N):
                if S[i,j]!=0:
                    S.col_swap(i,j)
                    V.col_swap(i,j)
                    break
        for j in range(i+1,N):
            if S[i,j]!=0:
                P = eye(N)
                P[i,j] = S[i,j]
                P[j,j] = -S[i,i]
                V = V*P
                #S = S*P
                S[j] = S[i,j]*S.col(i) - S[i,i]*S.col(j)
    #sort 0s
    for i in range(min(M,N)-1,-1,-1):
        for j in range(0,i):
            if S[j,j]==0:
                U.row_swap(j,j+1)
                V.col_swap(j,j+1)
                #S[j,j], S[j+1,j+1] = S[j+1,j+1], S[j,j]
                S.row_swap(j,j+1)
                S.col_swap(j,j+1)
    return [S, U, V]


def solve_linear_system_nontrivial(A0, b0, dummy='x'):
    """Finds the non-trivial solution of singular linear systems of equations
    """
    if not isinstance(A0, Matrix):
        raise TypeError("A should be a sympy Matrix")
    if not isinstance(b0, Matrix):
        raise TypeError("b should be a sympy Matrix")
    if A0.rows != b0.rows:
        raise ShapeError("Matrix size mismatch")
    M, N = A0.rows, A0.cols
    A, b = A0.copy(), b0.copy()
    #lower triangular part
    for i in range(0,min(M,N)):
        if A[i,i]==0:
            for j in range(i+1,M):
                if A[j,i]!=0:
                    A.row_swap(i,j)
                    b.row_swap(i,j)
                    break
        for j in range(i+1,M):
            if A[j,i]!=0:
                P = eye(M)
                P[j,i] = A[j,i]
                P[j,j] = -A[i,i]
                A = P*A
                b = P*b
                #A[j*N] = A[i,i]*A.row(j)-A[j,i]*A.row(i)     #TODO should replace A=P*A
                #b[j] = A[i,i]*b[j]-A[j,i]*b[i]
    #upper triangular part
    for i in range(min(M,N)-1,0,-1):
        if A[i,i]!=0:
            for j in range(0,i):
                if A[j,i]!=0:
                    P = eye(M)
                    P[j,i] = A[j,i]
                    P[j,j] = -A[i,i]
                    A = P*A
                    b = P*b
                    #A[j*N] = A[i,i]*A.row(j)-A[j,i]*A.row(i)    #TODO should replace A=P*A
                    #b[j] = A[i,i]*b[j]-A[j,i]*b[i]
    #zero diagonals
    for i in range(min(M,N)-2,-1,-1):
        if A[i,i]==0:
            for j in range(i+1,min(M,N)):
                if A[j,j]==0:
                    if A[i,j]!=0:
                        A.row_swap(i,j)
                        b.row_swap(i,j)
                        for k in range(0,j):
                            if A[k,j]!=0:
                                A[k*N] = A[j,j]*A.row(k)-A[k,j]*A.row(j)
                                b[k] = A[j,j]*b[k]-A[k,j]*b[j]
                        break
    #solve
    n_dummies = 0
    Conditions = []    #conditions expr_i==0
    x = zeros(M, 1)    #solution
    nnz = [0]*M    #number of nonzeros in each row
    for i in range(0,M):
        for j in range(0, N):
            if A[i,j]!=0:
                nnz[i]+=1
    nnz_dict = dict(zip(range(0, M), nnz))
    variable_solved = dict(zip(range(0, N), [[False, 0]]*N))
    for (i, v) in sorted(nnz_dict.items(), key=operator.itemgetter(1)):
        if v==0:
            Conditions.append(b[i])
        elif v==1:
            for j in range(0, N):
                if A[i,j]!=0:
                    variable_solved[j] = [True, b[i]/A[i,j]]
                    x[j] = b[i]/A[i,j]
                    break
        else:
            x_not_set = []    #variables not solved yet
            for j in range(0, N):
                if A[i,j]!=0 and variable_solved[j][0]==False:
                    x_not_set.append(j)
            for k in range(1, len(x_not_set)):
                x_dummy_next = dummy + str(n_dummies)    #next dummy variable
                n_dummies += 1
                variable_solved[x_not_set[k]] = [True, Symbol(x_dummy_next)]
                x[x_not_set[k]] = Symbol(x_dummy_next)
            if len(x_not_set)>0:
                x_ind = x_not_set[0]    #the index to solve (x_not_set[i>0] are set to dummies)
                x_val = b[x_ind]    #x_val:the value of the index to solve
                for j in range(0, N):
                    if A[i,j]!=0 and variable_solved[j][0]==True:
                        x_val -= A[i,j]*variable_solved[j][1]
                x_val /= A[i, x_ind]
                variable_solved[x_ind] = [True, x_val]
                x[x_ind] = x_val
    return [A, b, x, Conditions]


def SymFracSimplify(A):
    N, D = fraction(A)
    return simplify(N)/simplify(D)

def SymMatrixSimplify(A, frac=False):
    B = zeros(*A.shape)
    for i in range(A.rows):
        for j in range(A.cols):
            if frac==False:
                B[i, j] = A[i, j].simplify()
            else:
                B[i, j] = SymFracSimplify(A[i, j])
    return B

def SymMatrixdiff(A, x, n=1):
    B = zeros(A.rows, A.cols)
    for i in range(A.rows):
        for j in range(A.cols):
            B[i, j] = diff(A[i, j], x, n)
    return B

def SymMatrixdoit(A):
    for i in range(A.rows):
        for j in range(A.cols):
            A[i, j] = A[i, j].doit()
    return A

def SymMatrixsubs(A, a_a_sub):
    ## a_a_sub: list of tuples
    for i in range(A.rows):
        for j in range(A.cols):
            A[i, j] = A[i, j].subs(a_a_sub)
    return A
    
def SymEqSysLinearToMat(Eqs, vars):
    """ System of equation to Matrix
        Eqs: column matrix, each row representing 1 equation
        vars: list of variables
    """
    assert Eqs.cols == 1
    mat_coeff = zeros(Eqs.rows, len(vars))
    for i in range(Eqs.rows):
        for j in range(len(vars)):
            a_ij = Eqs[i]
            for k in range(len(vars)):
                if j==k:
                    a_ij = a_ij.subs(vars[k], 1)
                else:
                    a_ij = a_ij.subs(vars[k], 0)
            mat_coeff[i, j] = a_ij
    return mat_coeff
    

def SimplifyDiffEquationCoeffs(EQ, f, x, D_max):
    """ EQ: differential equation in f
        x: independant variable
        D_max: maximum derivative present in the equation
    """
    EQ_new = 0
    EQ_inhomog_part = EQ.subs(f, 0)     ## inhomogeneous part
    #EQ_inhomog_part = EQ_inhomog_part.subs(Derivative(Integer(0), x), 0)
    EQ = EQ - EQ_inhomog_part
    for i in range(D_max, -1, -1):
        coeff_i = EQ
        for j in range(D_max, -1, -1):
            if j!=i:
                coeff_i = coeff_i.subs(Derivative(f, x, j), 0)
            else:
                coeff_i = coeff_i.subs(Derivative(f, x, j), 1)
                
        coeff_i = coeff_i.doit().simplify()
        
        EQ_new = EQ_new + coeff_i*Derivative(f, x, i)
    return EQ_new + EQ_inhomog_part


def GetDiffEquationCoeffs(EQ, f, x, D_max):
    """ returns the coefficients of each derivative from 0 to D_max
    """
    EQ_inhomog_part = EQ.subs(f, 0)     ## inhomogeneous part
    #EQ_inhomog_part = EQ_inhomog_part.subs(Derivative(Integer(0), x), 0)
    EQ = EQ - EQ_inhomog_part       
    D_coeffs = [None]*(D_max+1)
    for i in range(D_max, -1, -1):
        coeff_i = EQ
        for j in range(D_max, -1, -1):
            if j!=i:
                coeff_i = coeff_i.subs(Derivative(f, x, j), 0)
            else:
                coeff_i = coeff_i.subs(Derivative(f, x, j), 1)
                
        coeff_i = coeff_i.doit().simplify()
        D_coeffs[i] = coeff_i
    return [D_coeffs, EQ_inhomog_part]


##---------------------  Muller's Method
from cmath import sqrt

def solveMuller(f, x_0, x_1, x_2, tolx_abs=None, toly_abs=1.0e-6,
        toly_rel=1.0e-6, n_max=1000, n_roots=1):
    """
    Muller's method for finding complex roots of function f
    """
    if n_roots==1:
        y_0 = f(x_0)
        y_1 = f(x_1)
        y_2 = f(x_2)

        #print('-:    x_i: ',  x_0, x_1, x_2, '\n   y_i: ', y_0, y_1, y_2)
        
        y_init = y_0
        message = None
        n_try = 0
        for i in range(n_max):
            y_01 = (y_1 - y_0)/(x_1 - x_0)
            y_12 = (y_2 - y_1)/(x_2 - x_1)
            y_02 = (y_2 - y_0)/(x_2 - x_0)

            y_012 = (y_12 - y_01)/(x_2 - x_0)

            w = y_01 + y_02 - y_12

            denom_p = w + sqrt(w**2 - 4*y_0*y_012)
            denom_m = w - sqrt(w**2 - 4*y_0*y_012)

            if abs(denom_p)>=abs(denom_m):
                x_next = x_0 - (2*y_0)/denom_p
            else:
                x_next = x_0 - (2*y_0)/denom_m

            y_next = f(x_next)
            break_loop = False
            if tolx_abs!=None:
                if abs(x_0-x_next)<tolx_abs:
                    break_loop = True 
            if toly_rel!=None:
                if abs(y_next)<toly_rel*abs(y_init):
                    break_loop = True
            if toly_abs!=None:
                if abs(y_next)<toly_abs:
                    break_loop = True
            if break_loop:
                x_0 = x_next
                y_0 = y_next
                n_try = i+1
                break
            
            x_2 = x_1;  y_2 = y_1
            x_1 = x_0;  y_1 = y_0
            x_0 = x_next
            y_0 = y_next
            
            #print('i: ', i, '   x_i: ',  x_0, '   y_i: ', y_0)
            if i==n_max-1:
                n_try = i+1
                message = 'Maximum number of iterations reached'
          
        return [x_0, y_0, n_try, message]
    else:
        assert n_roots>1
        roots = []
        resids = []
        n_tries = []
        message = None
        res = solveMuller(f, x_0, x_1, x_2, tolx_abs, toly_abs, toly_rel, n_max, n_roots=1)
        #print('solveMuller res: ', res)
        while res[3]==None and len(roots)<n_roots:
            roots.append(res[0])
            resids.append(res[1])
            n_tries.append(res[2])
            def f_new(x):
                denom = 1.0
                for x_i in roots:
                    denom *= (x - x_i)
                return f(x)/denom
            res = solveMuller(f_new, x_0, x_1, x_2, tolx_abs, toly_abs, toly_rel, n_max, n_roots=1)
            #print('solveMuller res: ', res)
            message = res[3]
        return [roots, resids, n_tries, message]
            
            
from scipy.optimize import root 

def RootsMultipleComplex(f, x_0=0.0, n_roots=2, solver='lm', tol=1.0e-10, ftol=1.0e-5,
        maxfev=1000, maxiter=1000):
    roots = []
    resids = []
    n_tries = []
    message = None
    while True:
        def f_ri(x_ri):
            x = x_ri[0] + 1j*x_ri[1]
            fx = f(x)
            for x_p in roots:
                fx /= (x - x_p)
            return [fx.real, fx.imag]
        x_0ri = [x_0.real, x_0.imag]
        sol = root(f_ri, x_0ri, method=solver, tol=tol, 
            options={'maxfev':maxfev, 'maxiter':maxiter, 'ftol': ftol})
        #print('sol: \n', sol)
        sol_x = sol.x
        sol_f = sol.fun
        if sol.success==True:
            roots.append(sol_x[0]+1j*sol_x[1])
            resids.append(sol_f[0]+1j*sol_f[1])
            n_tries.append(sol.nfev)
        else:
            message = sol.message
            break
        if len(roots)>=n_roots:
            break
    return [roots, resids, n_tries, message]


def lagrange1d_coeffs(x_s, x):
    # x_s: sample points    x : target point
    # f(x) = a[0]*f(x_s[0]) + a[1]*f(x_s[1]) + a[2]*f(x_s[2]) + ...
    N = len(x_s)
    a = [0.0]*N
    for i in range(N):
        a[i] = 1.0
        for j in range(N):
            if i != j:
                a[i] *= (x - x_s[j])/(x_s[i] - x_s[j])
    return a
                

def lagrange2d_coeffs(r_i, r):
    """ r_i = [(x0, y0), (x1, y1), (x2, y2), ...]  : sample points
        r = (x, y) : target point
        returns: [a0, a1, a2,...]
        f(r) = a0*f(r0) + a1*f(r1) + a2*f(r2) + ...
    """
    ## TODO: not tested
    N = len(r_i)
    a = [0.0]*N
    for i in range(N):
        a[i] = 1.0
        for j in range(N):
            if i != j:
                if r_i[i].x!=r_i[j].x and r_i[i].y!=r_i[j].y:
                    a[i] *= (r.x-r_i[j].x)*(r.y-r_i[j].y)/((r_i[i].x-r_i[j].x)*(r_i[i].y-r_i[j].y))
                elif r_i[i].x==r_i[j].x:
                    a[i] *= (r.y-r_i[j].y)/(r_i[i].y-r_i[j].y)
                else:
                    a[i] *= (r.x-r_i[j].x)/(r_i[i].x-r_i[j].x)
    return a

def barycentricInterpolation_coeffs(ri, r):
    """ r_i = [p0, p1, p2]  : sample points p0=Point2D(x0, y0)
        r = p : target point Point2D
        returns: [a0, a1, a2]
        f(r) = a0*f(r0) + a1*f(r1) + a2*f(r2)
    """
    A = (ri[1]-ri[0])^(ri[2]-ri[0])
    A0 = (ri[1]-r)^(ri[2]-r)
    A1 = (ri[2]-r)^(ri[0]-r)
    A2 = (ri[0]-r)^(ri[1]-r)    
    return [A0/A, A1/A, A2/A]



#-----------------  nullspace of a matrix
import numpy as np
import scipy as sp

def null(a, rtol=1.0e-6, forceRank=-1):
    ## TODO: the ratio of the smallest singular value to average is used for 
    ## determining the zero singular values.. The selection criteria can be improved
    u, s, v_h = sp.linalg.svd(a)
    s_avg = (np.abs(s)).sum()/len(s)
    #print(' s_avg: ', s_avg)
    rank = (np.abs(s) > abs(rtol*s_avg)).sum()
    if forceRank>0:
        rank = forceRank
    return [rank, v_h[rank:].T.conjugate(), s]


#----------------  timing
import time
time_start = 0
def tic():
    global time_start
    time_start = time.perf_counter()
    print('time_start: ', time_start)
    return time_start
    
def toc(t_start=None):
    if t_start==None:
        t_start = time_start
    time_end = time.perf_counter()
    print('time_end: ', time_end)
    seconds = time_end - t_start
    minute, seconds = divmod(seconds, 60)
    hour, minute = divmod(minute, 60)
    return "%d:%02d:%02d" % (hour, minute, seconds)
        


#---------------  strings -------

ch_spec = ['', ' ', '+', '-', '*', '.', '(', ')', '[', ']', '{', '}', '/']
def replace_whole_word(s, w, w_new):
    i_st = 0
    while True:
        n = s.find(w, i_st)
        if n<0:
            break
        i_st = n+1
        ch_prev = ''
        ch_next = ''
        if n>0:
            ch_prev = s[n-1]
        if n+len(w)<len(s):
            ch_next = s[n+len(w)]
        if (ch_prev in ch_spec) and (ch_next in ch_spec):
            a = ch_prev + w + ch_next
            b = ch_prev + w_new + ch_next
            s = s[0:n] + w_new + s[n+len(w):]
            i_st += len(w_new)-1
    return s


#---------------  point2D -------
import math

class Point2D:
    def __init__(self, x, y):
        self.x = x
        self.y = y
        return
        
    def __add__(self, p2):
        x = self.x + p2.x
        y = self.y + p2.y
        return Point2D(x, y)
        
    def __sub__(self, p2):
        x = self.x - p2.x
        y = self.y - p2.y
        return Point2D(x, y)
        
    def __neg__(self):
        x = -self.x
        y = -self.y
        return Point2D(x, y)
        
    def __iadd__(self, p2):
        self.x += p2.x
        self.y += p2.y
        return self
    
    def __mul__(self, d):
        x = self.x * d
        y = self.y * d
        return Point2D(x, y)

    def __truediv__(self, d):
        x = self.x / d
        y = self.y / d
        return Point2D(x, y)
        
    def __and__(self, p2):
        return self.x*p2.x + self.y*p2.y
        
    def __xor__(self, p2):
        return self.x*p2.y - self.y*p2.x
        
    def __str__(self):
        return "({}, {})".format(self.x, self.y)

    def __getitem__(self, i):
        if i==0:
            return self.x
        elif i==1:
            return self.y
        else:
            raise ValueError('Point2D: Index out of bounds')
 
    def __len__(self):
        return 2

    def norm(self):
        return math.sqrt(abs(self.x)**2 + abs(self.y)**2)
        
    def MagnitudeAngle(self):
        mag = self.norm()
        angle = math.atan2(self.y, self.x)
        if abs(angle)<1.0e-15:
            angle = 0.0
        return [mag, angle]
    


class Point3D:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
        return
        
    def __add__(self, p2):
        x = self.x + p2.x
        y = self.y + p2.y
        z = self.z + p2.z
        return Point3D(x, y, z)
        
    def __sub__(self, p2):
        x = self.x - p2.x
        y = self.y - p2.y
        z = self.z - p2.z
        return Point3D(x, y, z)
        
    def __neg__(self):
        x = -self.x
        y = -self.y
        z = -self.z
        return Point3D(x, y, z)
        
    def __iadd__(self, p2):
        self.x += p2.x
        self.y += p2.y
        self.z += p2.z
        return self
    
    def __mul__(self, d):
        x = self.x * d
        y = self.y * d
        z = self.z * d
        return Point3D(x, y, z)

    def __truediv__(self, d):
        x = self.x / d
        y = self.y / d
        z = self.z / d
        return Point3D(x, y, z)
        
    def __and__(self, p2):
        return self.x*p2.x + self.y*p2.y + self.z*p2.z
        
    def __xor__(self, p2):
        x = self.y*p2.z - self.z*p2.y
        y = self.z*p2.x - self.x*p2.z
        z = self.x*p2.y - self.y*p2.x
        return Point3D(x, y, z)
        
    def __str__(self):
        return "({}, {}, {})".format(self.x, self.y, self.z)
        
    def __getitem__(self, i):
        if i==0:
            return self.x
        elif i==1:
            return self.y
        elif i==2:
            return self.z
        else:
            raise ValueError('Point3D: Index out of bounds')
            
    def __len__(self):
        return 3
 
    def norm(self):
        return math.sqrt(abs(self.x)**2 + abs(self.y)**2 + abs(self.z)**2)
    
##------------------- crystallography

def GetReiprocalVecs(vecs):
    N = len(vecs)
    assert 1<=N<=3
    if N==1:
        d = vecs[0].norm()
        a_0 = vecs[0]/d
        b_0 = a_0*(2.0*math.pi/d)
        return [b_0]
    elif N==2:
        a_n = vecs[0]^vecs[1]
        a_n = a_n/a_n.norm()
        b_0 = (a_n^vecs[1]/(vecs[0]^vecs[1]).norm())*(-2.0*math.pi)
        b_1 = (a_n^vecs[0]/(vecs[0]^vecs[1]).norm())*(2.0*math.pi)
        return [b_0, b_1]
    elif N==3:
        b_0 = (vecs[1]^vecs[2])/(vecs[0]&(vecs[1]^vecs[2]))*(2.0*math.pi)
        b_1 = (vecs[2]^vecs[0])/(vecs[1]&(vecs[2]^vecs[0]))*(2.0*math.pi)
        b_2 = (vecs[0]^vecs[1])/(vecs[2]&(vecs[0]^vecs[1]))*(2.0*math.pi)
        return [b_0, b_1, b_2]
        
##---------------------------------

import numpy as np

def GetMaxComplexVec2D(F_x, F_y):
    """ Returns the maximum possible value of a time harmonic complex vector
    """
    f_x = np.abs(F_x)
    phi_x = np.angle(F_x)    #cmath.phase(F_x)
    f_y = np.abs(F_y)
    phi_y = np.angle(F_y)   #cmath.phase(F_y)

    A = f_x * f_x * np.cos(2.0 * phi_x) + f_y * f_y * np.cos(2.0 * phi_y)
    B = f_x * f_x * np.sin(2.0 * phi_x) + f_y * f_y * np.sin(2.0 * phi_y)
    norm = np.sqrt(A * A + B * B)
    A = A/norm
    B = B/norm

    res_sq = f_x * f_x * (1.0 + A * np.cos(2.0 * phi_x) + B * np.sin(2.0 * phi_x))\
            + f_y * f_y * (1.0 + A * np.cos(2.0 * phi_y) + B * np.sin(2.0 * phi_y))
    res_sq /= 2.0
    return np.sqrt(res_sq)

def NormalizeArrayTo0n1(A, logscale=False, logrange=5.0):
    """ normalizes numpy array to [0, 1]
    """
    if logscale==False:
        n = len(A)
        A_max = np.amax(A)
        A_min = np.amin(A)
        A_norm = np.zeros(n)
        denom = A_max - A_min
        if (denom == 0):
            return A_norm
        A_norm = (A - A_min) / denom
        return A_norm
    else:
        n = len(A)
        A_max = np.amax(A)
        A_min = np.amin(A)
        A_p = NormalizeArrayTo0n1(A, logscale=False)
        A_p = A_p*(A_p >= math.pow(10.0, -logrange)) + math.pow(10.0, -logrange)*(A_p < math.pow(10.0, -logrange))
        A_plog = np.log10(A_p)
        A_plog_01 = NormalizeArrayTo0n1(A_plog, logscale=False)
        return A_plog_01
    

from enum import Enum
class ColorMaps(Enum):
    BlackAndWhite = 1
    RedAndBlue = 2
    RedGreenBlue = 3
    
def GetColorMap(value, cmap=ColorMaps.RedGreenBlue):
    """ gets a value array between 0 and 1.0
        returns a color
    """
    if cmap==ColorMaps.RedAndBlue:
        color = [0]*len(value)
        for i in range(len(value)):
            value_01 = value[i]
            if (value_01 < 0.0): 
                value_01 = 0.0
            if (value_01 > 1.0): 
                value_01 = 1.0

            red = int(value_01 * 255)
            blue = int(255 - red)
            alpha = 255
            green = 0

            color[i] = '#{:02x}{:02x}{:02x}'.format(red, green, blue)
        return color
    elif (cmap == ColorMaps.RedGreenBlue):
        color = [0]*len(value)
        for i in range(len(value)):
            value_01 = value[i]
            if (value_01 < 0.0): 
                value_01 = 0.0
            if (value_01 > 1.0): 
                value_01 = 1.0

            red = 0
            blue = 0
            green = 0
            alpha = 255

            if value_01 <= 1.0 and value_01 >= 0.5:
                red = int((value_01 - 0.5) / 0.5 * 255)
                green = int(255 - red)
                blue = 0
            else:
                assert (value_01 >= 0.0 and value_01 <= 0.5)
                green = int(value_01 / 0.5 * 255)
                blue = int(255 - green)
                red = 0

            color[i] = '#{:02x}{:02x}{:02x}'.format(red, green, blue)
        return color
    else:
        raise NotImplementedError('cmap not implemented')



##----------------------  Fourier --------------------------------

def Fourier1D(x, fx):
    N = len(x)
    dx = x[1]-x[0]
    _Dx = 1/dx
    _x = np.linspace(-_Dx/2, _Dx/2, N, endpoint=False)
    _fx = np.fft.fftshift(np.fft.fft(fx))/N*np.exp(-1j*2.0*np.pi*x[0]*_x)
    return _x, _fx

def InvFourier1D(x, _x, _fx):
    N = len(x)
    fx = np.fft.ifft(np.fft.ifftshift(_fx*np.exp(+1j*2.0*np.pi*x[0]*_x)))*N
    return fx

def Fourier2D(x, y, fxy):
    Nx, Ny = x.shape
    dx = x[1,0]-x[0,0]
    dy = y[0,1]-y[0,0]
    _Dx = 1.0/dx
    _Dy = 1.0/dy
    _x = np.linspace(-_Dx/2, _Dx/2, Nx, endpoint=False)
    _y = np.linspace(-_Dy/2, _Dy/2, Ny, endpoint=False)
    _X, _Y = np.meshgrid(_x, _y, indexing='ij')
    _fxy = np.fft.fftshift(np.fft.fft2(fxy))/(Nx*Ny)*np.exp(-1j*2.0*np.pi*x[0,0]*_X)*np.exp(-1j*2.0*np.pi*y[0,0]*_Y)
    return _X, _Y, _fxy

def InvFourier2D(x, y, _x, _y, _fxy):
    Nx, Ny = x.shape
    fxy = np.fft.ifft2(np.fft.ifftshift(_fxy*np.exp(+1j*2.0*np.pi*x[0,0]*_x)*np.exp(+1j*2.0*np.pi*y[0,0]*_y)))*Nx*Ny
    return fxy




    
##----------------------  manage frequency points   -----------

from scipy.interpolate import lagrange

class ParamSweeper:
    """
    it adds the frequency points, decreases or increases the frequency steps 
    on demand, extrapolates the last frequency point
    """
    def __init__(self):
        self.X = []
        self.Y = []
        return
        
    def setMaxStepSize(self, dx, dx_0_ratio=8.0):
        self.dx = dx/dx_0_ratio
        self.dx_max = dx
        return
        
    def decreaseMaxStepSize(self):
        self.dx_max /= 2.0
        return

    def increaseMaxStepSize(self):
        self.dx_max *= 2.0
        return

    def setLimitPoints(self, x_0, x_N):
        self.x_0 = x_0
        self.x_N = x_N
        return
        
    def getLastX(self):
        return self.X[-1]
        
    def setLastY(self, y):
        self.Y[-1] = y
        return
        
    def RemoveLast(self):
        if len(self.X)>0:
            assert len(self.Y)==len(self.X)
            del self.X[-1]
            del self.Y[-1]
        
    def RemoveNones(self):
        N = len(self.X)
        assert len(self.Y)==N
        for i in range(N-1, -1, -1):
            if self.Y[i]==None:
                del self.X[i]
                del self.Y[i]
        
    def printXY(self):
        print('X : ', self.X)
        print('Y : ', self.Y)

    def addNext(self):        
        if len(self.X)==0:
            self.X.append(self.x_0)
            self.Y.append(None)
        else:
            if abs(self.dx) < abs(self.dx_max):
                self.dx *= 2.0
                if abs(self.dx) > abs(self.dx_max):
                    self.dx = self.dx_max
            x_next = self.X[len(self.X)-1] + self.dx
            if (self.dx>0.0 and x_next>self.x_N) or \
                    (self.dx<0.0 and x_next<self.x_N):
                return False
            if self.dx==0.0:
                return False
            self.X.append(x_next)
            self.Y.append(None)
        return True
    
    def refineLast(self):
        del self.X[-1]
        del self.Y[-1]
        self.dx /= 4.0
        #self.addNext()
        return

    def getXY(self):
        return [self.X[:], self.Y[:]]

    def extrapolateLast(self, n=3, seprateRI=False):
        if seprateRI==False:
            return self.extrapolateLastComplexY(n=n)
        else:
            y_ext_r = self.extrapolateLastRealY(n=n)
            y_ext_i = self.extrapolateLastImagY(n=n)
            if y_ext_r==None or y_ext_i==None:
                return None
            else:
                ny = len(y_ext_r)
                assert len(y_ext_i)==ny
                y_ext = [None]*ny
                for i in range(ny):
                    y_ext[i] = y_ext_r[i] + y_ext_i[i]*1j
                return y_ext

    def extrapolateLastComplexY(self, n):
        N_xy = len(self.X)
        X = [None]*N_xy
        Y = [None]*N_xy
        for i in range(N_xy):
            X[i] = self.X[i]
            Y[i] = self.Y[i]
            
        if len(X)>n:
            res = []
            for i in range(len(Y[0])):
                Y_i = [Y[j][i] for j in range(len(Y)-1)]
                P_lag = lagrange(X[-n:-1], Y_i[-n+1:])
                res.append(np.polyval(P_lag, X[-1]))
            return res
        elif len(X)>2:
            m = len(X) - 1
            res = []
            for i in range(len(Y[0])):
                Y_i = [Y[j][i] for j in range(len(Y)-1)]
                P_lag = lagrange(X[-m:-1], Y_i[-m+1:])
                res.append(np.polyval(P_lag, X[-1]))
            return res
        elif len(X)==2:
            return Y[0]
        else:
            return None
        return None

    def extrapolateLastRealY(self, n):
        N_xy = len(self.X)
        X = [None]*N_xy
        Y = [None]*N_xy
        for i in range(N_xy):
            X[i] = self.X[i]
            if self.Y[i]==None:
                Y[i] = None
            else:
                N_yi = len(self.Y[i])
                Y[i] = [None]*N_yi
                for j in range(N_yi):
                    Y[i][j] = self.Y[i][j].real
            
        if len(X)>n:
            res = []
            for i in range(len(Y[0])):
                Y_i = [Y[j][i] for j in range(len(Y)-1)]
                P_lag = lagrange(X[-n:-1], Y_i[-n+1:])
                res.append(np.polyval(P_lag, X[-1]))
            return res
        elif len(X)>2:
            m = len(X) - 1
            res = []
            for i in range(len(Y[0])):
                Y_i = [Y[j][i] for j in range(len(Y)-1)]
                P_lag = lagrange(X[-m:-1], Y_i[-m+1:])
                res.append(np.polyval(P_lag, X[-1]))
            return res
        elif len(X)==2:
            return Y[0]
        else:
            return None
        return None
        
    def extrapolateLastImagY(self, n):
        N_xy = len(self.X)
        X = [None]*N_xy
        Y = [None]*N_xy
        for i in range(N_xy):
            X[i] = self.X[i]
            if self.Y[i]==None:
                Y[i] = None
            else:
                N_yi = len(self.Y[i])
                Y[i] = [None]*N_yi
                for j in range(N_yi):
                    Y[i][j] = self.Y[i][j].imag
            
        if len(X)>n:
            res = []
            for i in range(len(Y[0])):
                Y_i = [Y[j][i] for j in range(len(Y)-1)]
                P_lag = lagrange(X[-n:-1], Y_i[-n+1:])
                res.append(np.polyval(P_lag, X[-1]))
            return res
        elif len(X)>2:
            m = len(X) - 1
            res = []
            for i in range(len(Y[0])):
                Y_i = [Y[j][i] for j in range(len(Y)-1)]
                P_lag = lagrange(X[-m:-1], Y_i[-m+1:])
                res.append(np.polyval(P_lag, X[-1]))
            return res
        elif len(X)==2:
            return Y[0]
        else:
            return None
        return None




##--------------------  Parameter reader: reads input parameters from file

class TextInputParameterReader:
    
    def __init__(self, file_name):
        self.file_name = file_name
        myfile = open(file_name, 'r')
        file_data = myfile.read()
        myfile.close()
        self.file_data = file_data
        
        
    def find_block(self, data, param):
        i_param = data.find(param)
        if i_param>=0:
            i_beg = data.find(r'{', i_param+len(param))
            i_end = data.find(r'}', i_param+len(param))
            if i_beg>=0 and i_end>=0 and i_end>i_beg:
                param = data[i_beg+1:i_end]
                return param
        return None

    def find_param(self, data, param):
        i_param = data.find(param)
        if i_param>=0:
            i_beg = data.find(r'(', i_param+len(param))
            i_end = data.find(r')', i_param+len(param))
            if i_beg>=0 and i_end>=0 and i_end>i_beg:
                param = data[i_beg+1:i_end]
                return param.split(',')
        return None
            
    def liststr_to_floats(self, list_str):
        list_float = [float(list_str[i]) for i in range(len(list_str))]
        return list_float
            




