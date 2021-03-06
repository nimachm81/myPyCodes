## GSTCTD.py time domain generalized sheet transition condition (GSTC)

__all__ = ["GSTCTD1D", "GSTCTD1Dcond"]


import numpy as np
from scipy.integrate import romb, quad, simps
from sympy import lambdify
from scipy.optimize import brentq, bisect

class GSTCTD1D:
    
    def __init__(self):
        self.eps_0 = 1.0
        self.mu_0 = 1.0
        self.eta_0 = np.sqrt(self.mu_0/self.eps_0)
        self.c = 1.0/np.sqrt(self.mu_0*self.eps_0)
        return
        
        
    def SetFieldsSympy(self, t, z, E_i, E_r, E_t):
        self.t_sym = t
        self.z_sym = z
        self.E_i_sym = E_i
        self.E_r_sym = E_r
        self.E_t_sym = E_t
        self.E_i = lambdify((t, z), self.E_i_sym, 'numpy')
        self.E_r = lambdify((t, z), self.E_r_sym, 'numpy')
        self.E_t = lambdify((t, z), self.E_t_sym, 'numpy')
        return
        
        
    def GetPMXeeXmm(self, t):
        """ t: numpt vector
        """
        P_t = np.zeros(len(t))
        M_t = np.zeros(len(t))
        E_i, E_r, E_t = self.E_i, self.E_r, self.E_t
        eps_0, mu_0, eta_0 = self.eps_0, self.mu_0, self.eta_0
        for i in range(len(t)):
            P_t[i] = quad(lambda t: -(E_t(t, 0.0)-E_i(t, 0.0)+E_r(t, 0.0))/eta_0, t[0], t[i])[0]
            M_t[i] = quad(lambda t: 1.0/mu_0*(E_t(t, 0.0)-E_i(t, 0.0)-E_r(t, 0.0)), t[0], t[i])[0]

        ##TODO: prevent devide by zeros
        X_ee = (2.0*P_t)/(E_t(t, 0.0)+E_i(t, 0.0)+E_r(t, 0.0))/eps_0
        X_mm = (2.0*eta_0*M_t)/(E_t(t, 0.0)+E_i(t, 0.0)-E_r(t, 0.0))

        self.t = t
        self.P_t = P_t
        self.M_t = M_t
        self.X_ee = X_ee
        self.X_mm = X_mm
        return



    def TimeWindowSusceptibilities(self, wnd_ee, wnd_mm):
        self.wnd_ee = wnd_ee
        self.wnd_mm = wnd_mm
        t = self.t
        wnd = wnd_ee
        self.X_ee = self.X_ee*np.logical_and(t>wnd[0], t<wnd[1])
        wnd = wnd_mm
        self.X_mm = self.X_mm*np.logical_and(t>wnd[0], t<wnd[1])



    def ClipSusceptibilities(self, X_ee_max, X_mm_max):
        assert X_ee_max>0 and X_mm_max>0
        t = self.t
        X_ee, X_mm = self.X_ee, self.X_mm
        X_ee[np.logical_not(np.isfinite(X_ee))] = X_ee_max 
        X_mm[np.logical_not(np.isfinite(X_mm))] = X_mm_max
        
        X_min, X_max = -X_ee_max, X_ee_max
        self.X_ee =  X_ee*(np.logical_and(X_ee<=X_max,X_ee>=X_min)) + (X_ee>X_max)*X_max + (X_ee<X_min)*X_min
        X_min, X_max = -X_mm_max, X_mm_max
        self.X_mm =  X_mm*(np.logical_and(X_mm<=X_max,X_mm>=X_min)) + (X_mm>X_max)*X_max + (X_mm<X_min)*X_min



    def FindSusceptZeros(self):
        t = self.t
        X_ee = self.X_ee
        X_mm = self.X_mm
        X_ee_intervals = []
        X_mm_intervals = []
        
        X_ee_max = np.max(X_ee)
        X_mm_max = np.max(X_mm)
        
        print('X_ee_max :', X_ee_max)
        
        X_ee_small = 0.5
        X_mm_small = 0.5
        
        N = len(t)
        for i in range(N-1):
            if X_ee[i]*X_ee[i+1]<=0.0:
                if X_ee[i+1]!=0:
                    if np.abs(X_ee[i+1]-X_ee[i])<X_ee_small:
                        X_ee_intervals.append([t[i], t[i+1]])
                if i>0 and X_ee[i]==0 and X_ee[i+1]==0 and X_ee[i-1]!=0:
                    if np.abs(X_ee[i-1]-X_ee[i])<X_ee_small:
                        X_ee_intervals.append([t[i-1], t[i]])
            if X_mm[i]*X_mm[i+1]<=0.0:
                if X_mm[i+1]!=0:
                    if np.abs(X_mm[i+1]-X_mm[i])<X_mm_small:
                        X_mm_intervals.append([t[i], t[i+1]])
                if i>0 and X_mm[i]==0 and X_mm[i+1]==0 and X_mm[i-1]!=0:
                    if np.abs(X_mm[i-1]-X_mm[i])<X_mm_small:
                        X_mm_intervals.append([t[i-1], t[i]])
        if X_ee[N-1]==0.0 and X_ee[N-2]!=0.0:
            if np.abs(X_ee[i+1]-X_ee[i])<X_ee_max/10.0:
                X_ee_intervals.append([t[N-1], t[N-2]])
        if X_mm[N-1]==0.0 and X_mm[N-2]!=0.0:
            if np.abs(X_mm[i+1]-X_mm[i])<X_mm_max/10.0:
                X_mm_intervals.append([t[N-1], t[N-2]])

        self.X_ee_intervals = X_ee_intervals
        self.X_mm_intervals = X_mm_intervals
        print('self.X_ee_intervals:', self.X_ee_intervals)
        print('self.X_mm_intervals:', self.X_mm_intervals)
        
        
        X_ee_roots = []
        X_mm_roots = []
        
        E_i, E_r, E_t = self.E_i, self.E_r, self.E_t
        eps_0, mu_0, eta_0 = self.eps_0, self.mu_0, self.eta_0
        t_0 = self.t[0]
        def f_X_ee(t):
            P_t = quad(lambda t: -(E_t(t, 0.0)-E_i(t, 0.0)+E_r(t, 0.0))/eta_0, t_0, t, limit=300)[0]
            X_ee_t = (2.0*P_t)/(E_t(t, 0.0)+E_i(t, 0.0)+E_r(t, 0.0))/eps_0
            return X_ee_t
            
        def f_X_mm(t):
            M_t = quad(lambda t: 1.0/mu_0*(E_t(t, 0.0)-E_i(t, 0.0)-E_r(t, 0.0)), t_0, t, limit=300)[0]
            X_mm_t = (2.0*eta_0*M_t)/(E_t(t, 0.0)+E_i(t, 0.0)-E_r(t, 0.0))
            return X_mm_t

        for i in range(len(X_ee_intervals)):
            a, b = X_ee_intervals[i]
            if f_X_ee(a)*f_X_ee(b)<=0.0:
                t_root = bisect(f_X_ee, a, b)
                X_ee_roots.append(t_root)

        for i in range(len(X_mm_intervals)):
            a, b = X_mm_intervals[i]
            t_root = None
            if f_X_mm(a)*f_X_mm(b)<=0.0:
                t_root = bisect(f_X_mm, a, b)
                X_mm_roots.append(t_root)

        self.X_ee_roots = X_ee_roots
        self.X_mm_roots = X_mm_roots
        print('self.X_ee_roots : ', self.X_ee_roots)


    def AdjustWindowToZeros(self):
        wnd_ee, wnd_mm = None, None
        if len(self.X_ee_roots)>=2:
            wnd_ee = [self.X_ee_roots[0], self.X_ee_roots[-1]]
        if len(self.X_mm_roots)>=2:
            wnd_mm = [self.X_mm_roots[0], self.X_mm_roots[-1]]
        if wnd_ee is not None and wnd_mm is not None:
            self.TimeWindowSusceptibilities(wnd_ee, wnd_mm)
            
    
    def GetSubintervalsWithSingleZero(self):
        ## it is assumed the windows start and end with zeros (of susceptibilities)
        intervals_ee = []    
        intervals_mm = []
        
        intervals_ee.append(self.wnd_ee[0])
        dw_ee = self.wnd_ee[1] - self.wnd_ee[0]
        for i in range(len(self.X_ee_roots)):
            if self.X_ee_roots[i]>intervals_ee[-1]: #abs(self.X_ee_roots[i]-intervals_ee[-1])/dw_ee>1.0e-5:
                intervals_ee.append(self.X_ee_roots[i])
        if self.wnd_ee[1]>intervals_ee[-1]:
            intervals_ee.append(self.wnd_ee[1])
        
        intervals_mm.append(self.wnd_mm[0])
        dw_mm = self.wnd_mm[1] - self.wnd_mm[0]
        for i in range(len(self.X_mm_roots)):
            if self.X_mm_roots[i]>intervals_mm[-1]: #abs(self.X_mm_roots[i]-intervals_mm[-1])/dw_mm>1.0e-5:
                intervals_mm.append(self.X_mm_roots[i])
        if self.wnd_mm[1]>intervals_mm[-1]:
            intervals_mm.append(self.wnd_mm[1])
        
        N = len(intervals_ee)
        subintervals_ee = [None]*(N+1)
        subinterval_zeros_ee = [None]*N
        subintervals_ee[0] = intervals_ee[0]
        subinterval_zeros_ee[0] = intervals_ee[0]   ## assuming window starts with a zero
        for i in range(1, N):
            subintervals_ee[i] = (intervals_ee[i-1]+intervals_ee[i])/2.0
            subinterval_zeros_ee[i] = intervals_ee[i]
        subintervals_ee[N] = intervals_ee[N-1]

        N = len(intervals_mm)
        subintervals_mm = [None]*(N+1)
        subinterval_zeros_mm = [None]*N
        subintervals_mm[0] = intervals_mm[0]
        subinterval_zeros_mm[0] = intervals_mm[0]   ## assuming window starts with a zero
        for i in range(1, N):
            subintervals_mm[i] = (intervals_mm[i-1]+intervals_mm[i])/2.0
            subinterval_zeros_mm[i] = intervals_mm[i]
        subintervals_mm[N] = intervals_mm[N-1]

        self.subintervals_ee = subintervals_ee
        self.subintervals_mm = subintervals_mm
        self.subinterval_zeros_ee = subinterval_zeros_ee
        self.subinterval_zeros_mm = subinterval_zeros_mm
        

    def GetSubintervalsWithSingleZero__(self):
        ## it is assumed the windows start and end with zeros (of susceptibilities)
        intervals_ee = []    
        intervals_mm = []
        
        intervals_ee.append(self.wnd_ee[0])
        dw_ee = self.wnd_ee[1] - self.wnd_ee[0]
        for i in range(len(self.X_ee_roots)):
            if self.X_ee_roots[i]>intervals_ee[-1]: #abs(self.X_ee_roots[i]-intervals_ee[-1])/dw_ee>1.0e-5:
                intervals_ee.append(self.X_ee_roots[i])
        if self.wnd_ee[1]>intervals_ee[-1]:
            intervals_ee.append(self.wnd_ee[1])
        
        intervals_mm.append(self.wnd_mm[0])
        dw_mm = self.wnd_mm[1] - self.wnd_mm[0]
        for i in range(len(self.X_mm_roots)):
            if self.X_mm_roots[i]>intervals_mm[-1]: #abs(self.X_mm_roots[i]-intervals_mm[-1])/dw_mm>1.0e-5:
                intervals_mm.append(self.X_mm_roots[i])
        if self.wnd_mm[1]>intervals_mm[-1]:
            intervals_mm.append(self.wnd_mm[1])
        
        N = len(intervals_ee)
        subintervals_ee = [None]*(2*N-1)
        subinterval_zeros_ee = [None]*(2*N-2)
        subintervals_ee[0] = intervals_ee[0]
        subinterval_zeros_ee[0] = intervals_ee[0]   ## assuming window starts with a zero
        for i in range(1, N):
            subintervals_ee[2*i-1] = (intervals_ee[i-1]+intervals_ee[i])/2.0
            subintervals_ee[2*i  ] = intervals_ee[i]
            subinterval_zeros_ee[2*i-1] = intervals_ee[i]
            if i<N-1:
                subinterval_zeros_ee[2*i  ] = intervals_ee[i]

        N = len(intervals_mm)
        subintervals_mm = [None]*(2*N-1)
        subinterval_zeros_mm = [None]*(2*N-2)
        subintervals_mm[0] = intervals_mm[0]
        subinterval_zeros_mm[0] = intervals_mm[0]   ## assuming window starts with a zero
        for i in range(1, N):
            subintervals_mm[2*i-1] = (intervals_mm[i-1]+intervals_mm[i])/2.0
            subintervals_mm[2*i  ] = intervals_mm[i]
            subinterval_zeros_mm[2*i-1] = intervals_mm[i]
            if i<N-1:
                subinterval_zeros_mm[2*i  ] = intervals_mm[i]

        self.subintervals_ee = subintervals_ee
        self.subintervals_mm = subintervals_mm
        self.subinterval_zeros_ee = subinterval_zeros_ee
        self.subinterval_zeros_mm = subinterval_zeros_mm
        
        
    def SolvePinSubinterval(self, P_0, interval, t_p, dt_min, n_t_min=10):
        ## a = 2/(eps_0*eta_0)*X_ee_inv*(t-t_p)
        ## g = 2/eta_0*E_i
        t_0, t_1 = interval
        n_t = max(int((t_1 - t_0)/dt_min+1), n_t_min)
        t = np.linspace(t_0, t_1, n_t)
        assert len(t)==n_t
        
        ## get X_ee_inv*(t-t_p)
        E_i, E_r, E_t = self.E_i, self.E_r, self.E_t
        eps_0, mu_0, eta_0 = self.eps_0, self.mu_0, self.eta_0
        X_ee_inv_t_tp = np.zeros(n_t)
        for i in range(n_t):
            t_i = t[i]
            P_t = quad(lambda t: -(E_t(t, 0.0)-E_i(t, 0.0)+E_r(t, 0.0))/eta_0, self.t[0], t_i, limit=300)[0]
            if np.abs(t[i] - t_p)>1.0e-6:
                X_ee_inv_t_tp[i] = (E_t(t_i, 0.0)+E_i(t_i, 0.0)+E_r(t_i, 0.0))*eps_0/(2.0*P_t)*(t_i - t_p)
            else:
                t_i = t_p+1.0e-6
                P_t = quad(lambda t: -(E_t(t, 0.0)-E_i(t, 0.0)+E_r(t, 0.0))/eta_0, self.t[0], t_i)[0]
                X_ee_inv_t_tp[i] = (E_t(t_i, 0.0)+E_i(t_i, 0.0)+E_r(t_i, 0.0))*eps_0/(2.0*P_t)*(t_i - t_p)
        
        a = 2.0/(eps_0*eta_0)*X_ee_inv_t_tp
        g = 2.0/eta_0*E_i(t, 0.0)
        
        #print('X_ee_inv_t_tp:', X_ee_inv_t_tp)
        
        P = self.SolvePde1D(t, a, g, t_p, f_0=P_0)
        return [t, P]
        

    def SolveMinSubinterval(self, M_0, interval, t_p, dt_min, n_t_min=10):
        ## a = 2*eta_0/(mu_0)*X_mm_inv*(t-t_p)
        ## g = 2/mu_0*E_i
        t_0, t_1 = interval
        n_t = max(int((t_1 - t_0)/dt_min+1), n_t_min)
        t = np.linspace(t_0, t_1, n_t)
        assert len(t)==n_t
        
        ## get X_mm_inv*(t-t_p)
        E_i, E_r, E_t = self.E_i, self.E_r, self.E_t
        eps_0, mu_0, eta_0 = self.eps_0, self.mu_0, self.eta_0
        X_mm_inv_t_tp = np.zeros(n_t)
        for i in range(n_t):
            t_i = t[i]
            M_t = quad(lambda t: 1.0/mu_0*(E_t(t, 0.0)-E_i(t, 0.0)-E_r(t, 0.0)), self.t[0], t_i, limit=300)[0]
            if np.abs(t[i] - t_p)>1.0e-6:
                X_mm_inv_t_tp[i] = (E_t(t_i, 0.0)+E_i(t_i, 0.0)-E_r(t_i, 0.0))/(2.0*eta_0*M_t)*(t_i - t_p)
            else:
                t_i = t_p+1.0e-6
                M_t = quad(lambda t: 1.0/mu_0*(E_t(t, 0.0)-E_i(t, 0.0)-E_r(t, 0.0)), self.t[0], t_i, limit=300)[0]
                X_mm_inv_t_tp[i] = 0#(E_t(t_i, 0.0)+E_i(t_i, 0.0)-E_r(t_i, 0.0))/(2.0*eta_0*M_t)*(t_i - t_p)
        
        a = 2.0*eta_0/(mu_0)*X_mm_inv_t_tp
        g = 2.0/mu_0*E_i(t, 0.0)
        
        #print('X_mm_inv_t_tp:', X_mm_inv_t_tp)
        #print('a:', a)
        #print('t_p ind:', np.argmax(t>=t_p))
        #print('ind total:', len(t))
        
        M = self.SolvePde1D(t, a, g, t_p, f_0=M_0)
        return [t, M]


    def SolveP(self, dt_min, n_t_min=10):
        N = len(self.subintervals_ee)-1
        P_0 = 0.0
        t_P_list = []
        for i in range(N):
            interval = [self.subintervals_ee[i], self.subintervals_ee[i+1]]
            t_p = self.subinterval_zeros_ee[i]
            t, P = self.SolvePinSubinterval(P_0, interval, t_p, dt_min, n_t_min)
            P_0 = P[-1]
            t_P_list.append([t, P])
            
        return t_P_list


    def SolveM(self, dt_min, n_t_min=10):
        N = len(self.subintervals_mm)-1
        M_0 = 0.0
        t_M_list = []
        for i in range(N):
            interval = [self.subintervals_mm[i], self.subintervals_mm[i+1]]
            t_p = self.subinterval_zeros_mm[i]
            t, M = self.SolveMinSubinterval(M_0, interval, t_p, dt_min, n_t_min)
            M_0 = M[-1]
            t_M_list.append([t, M])
            
        return t_M_list


    def SolvePde1D(self, t, a, g, t_p, f_0):
        """ df/dt + a(t)*f/(t-t_p) = g
        """
        dt = t[1]-t[0]
        assert t.shape==a.shape==g.shape
        N = len(t)

        def GetTheta(a, t_p):
            i__p1 = np.argmax(t>t_p)+1
            if t[-1]<=t_p:
                i__p1 = N
            i__m1 = np.argmax(t>t_p)-2
            if t[0]>=t_p:
                i__m1 = -1
            #da/dt
            a_p = np.zeros(N)
            a_p[0:N-1] = (a[1:N] - a[0:N-1])/dt
            a_p[N-1] = (a[N-1]-a[N-2])/dt
            #
            def LnInteg(x):
                if x==0.0:
                    return 0.0
                return x*np.log(x)-x
            
            theta = np.zeros(N)
            for i in range(N):
                if t[i]>=t_p:
                    b = 0.0
                    if i>i__p1:
                        I_ln_dt = LnInteg(np.abs(t[i__p1]-t_p))*a_p[i__p1-1]
                        I_arg = a_p[i__p1:i+1]*np.log(np.abs(t[i__p1:i+1]-t_p))
                        I_ln = simps(I_arg, dx=dt)
                        b = np.exp(I_ln_dt + I_ln)
                    else:
                        I_ln_dt = LnInteg(np.abs(t[i]-t_p))*a_p[i]
                        b = np.exp(+I_ln_dt)
                    theta[i] = np.abs(t[i]-t_p)**(a[i])/b
                    if t[i]==t_p and a[i]<0.0:
                        theta[i] = 1.0e100
                    #print('b: ', b, theta[i])
                else:
                    b = 0.0
                    if i<i__m1:
                        I_ln_dt = LnInteg(np.abs(t[i__m1]-t_p))*a_p[i__p1-1]
                        I_arg = a_p[i:i__m1+1]*np.log(np.abs(t[i:i__m1+1]-t_p))
                        I_ln = simps(I_arg, dx=dt)
                        b = np.exp(-I_ln_dt - I_ln)
                    else:
                        I_ln_dt = LnInteg(np.abs(t[i]-t_p))*a_p[i]
                        b = np.exp(-I_ln_dt)
                    theta[i] = np.abs(t[i]-t_p)**a[i]/b
                    #print('b: ', b, theta[i])
            return theta
        
        theta = GetTheta(a, t_p)
                
        i__tp = np.argmax(t>=t_p)
        if t[-1]<t_p:
            i__tp = None
        if t[0]>t_p:
            i__tp = None

        #assert i__tp!=None

        
        if i__tp!=None and a[i__tp]>=0.0:
            int_theta_g = np.zeros(N)
            for i in range(N):
                if i==i__tp:
                    int_theta_g[i] = 0.0
                elif i<i__tp:
                    I_arg = theta[i:i__tp+1]*g[i:i__tp+1]
                    int_theta_g[i] = -simps(I_arg, dx=dt)
                else:
                    I_arg = theta[i__tp:i+1]*g[i__tp:i+1]
                    int_theta_g[i] = simps(I_arg, dx=dt)
            
            c = 0.0
            
            is_close_tp = np.abs(t-t_p)<=2.0*dt
            f = (int_theta_g + c)/(theta)*np.logical_not(is_close_tp) + \
                (t - t_p)*g/a*is_close_tp
            return f
        elif False and i__tp!=None and a[i__tp]<0.0:
            t_p_rev = t[0]+(t[-1]-t_p)
            g_rev = -g[::-1]
            a_rev = -a[::-1]
            theta_rev = GetTheta(a_rev, t_p_rev)
        
            int_theta_g = np.zeros(N)
            for i in range(N):
                if i==i__tp:
                    int_theta_g[i] = 0.0
                elif i<i__tp:
                    I_arg = theta_rev[i:i__tp+1]*g_rev[i:i__tp+1]
                    int_theta_g[i] = -simps(I_arg, dx=dt)
                else:
                    I_arg = theta_rev[i__tp:i+1]*g_rev[i__tp:i+1]
                    int_theta_g[i] = simps(I_arg, dx=dt)
            
            c = 0.0
            
            is_close_tp = np.abs(t-t_p_rev)<=2.0*dt
            f = (int_theta_g + c)/(theta_rev)*np.logical_not(is_close_tp) + \
                (t - t_p_rev)*g_rev/a_rev*is_close_tp
            return f[::-1]
        elif i__tp!=None and a[i__tp]<0.0:
            int_theta_g = np.zeros(N)
            
            #n_di = 3
            #i_l = max(0, i__tp-n_di)
            #i_u = min(N-1, i__tp+n_di)
            #for i in range(i_l, i_u+1):
            for i in range(0, N):
                if i==i__tp:
                    int_theta_g[i] = 0.0
                elif i<i__tp:
                    I_arg = theta[0:i+1]*g[0:i+1]
                    int_theta_g[i] = simps(I_arg, dx=dt)
                elif i<len(g)-1:
                    I_arg = theta[i:N]*g[i:N]
                    int_theta_g[i] = -simps(I_arg, dx=dt)
            
            """
            c_mid = 0.0
            f_l = (int_theta_g[i_l] + c_mid)/(theta[i_l])
            f_u = (int_theta_g[i_u] + c_mid)/(theta[i_u])

            for i in range(i_l):
                I_arg = theta[i:i_l+1]*g[i:i_l+1]
                int_theta_g[i] = -simps(I_arg, dx=dt)
            for i in range(i_u+1, N):
                I_arg = theta[i_u:i+1]*g[i_u:i+1]
                int_theta_g[i] = simps(I_arg, dx=dt)
            
            
            c_l = theta[i_l-1]*f_l
            c_u = theta[i_u+1]*f_u
                        
            is_close_tp = np.abs(t-t_p)<=2.0*dt
            is_in_l = (np.arange(N)<i_l)*np.logical_not(is_close_tp)
            is_in_u = (np.arange(N)>i_u)*np.logical_not(is_close_tp)
            is_in_mid = (np.arange(N)>=i_l)*(np.arange(N)<=i_u)*np.logical_not(is_close_tp)

            f = (int_theta_g + c_l)/(theta)*is_in_l + (int_theta_g + c_u)/(theta)*is_in_u +\
                (int_theta_g + c_mid)/(theta)*is_in_mid + (t - t_p)*g/a*is_close_tp
            """
            
            c = 0.0
            
            is_close_tp = np.abs(t-t_p)<=2.0*dt
            f = (int_theta_g + c)/(theta)*np.logical_not(is_close_tp) + \
                (t - t_p)*g/a*is_close_tp
            return f
        elif False and i__tp!=None and a[i__tp]<0.0:
            int_theta_g = np.zeros(N)
            
            for i in range(0, N):
                if i==i__tp:
                    int_theta_g[i] = 0.0
                elif i<i__tp:
                    I_arg = theta[0:i+1]*g[0:i+1]
                    int_theta_g[i] = simps(I_arg, dx=dt)
                elif i<len(g)-1:
                    I_arg = theta[i:N]*g[i:N]
                    int_theta_g[i] = -simps(I_arg, dx=dt)
            
            
            i_l = i__tp-4
            i_u = i__tp+4
            c_l = -int_theta_g[i_l] + (t[i_l] - t_p)*g[i_l]/a[i_l]
            c_u = -int_theta_g[i_u] + (t[i_u] - t_p)*g[i_u]/a[i_u]
            print('c_l:', c_l, '  c_u:', c_u)
                        
            is_close_tp = np.abs(t-t_p)<=2.0*dt
            is_in_l = (np.arange(N)<i_l)*np.logical_not(is_close_tp)
            is_in_u = (np.arange(N)>i_u)*np.logical_not(is_close_tp)
            is_in_mid = (np.arange(N)>=i_l)*(np.arange(N)<=i_u)*np.logical_not(is_close_tp)

            f = (int_theta_g + c_l)/(theta)*is_in_l + (int_theta_g + c_u)/(theta)*is_in_u +\
                (t - t_p)*g/a*is_close_tp
            
            return f
        else:
            int_theta_g = np.zeros(N)
            for i in range(1, N):
                I_arg = theta[0:i+1]*g[0:i+1]
                int_theta_g[i] = simps(I_arg, dx=dt)

            c = theta[0]*f_0

            f = (int_theta_g + c)/(theta)
            return f



class GSTCTD1Dcond:
    
    def __init__(self):
        self.eps_0 = 1.0
        self.mu_0 = 1.0
        self.eta_0 = np.sqrt(self.mu_0/self.eps_0)
        self.c = 1.0/np.sqrt(self.mu_0*self.eps_0)
        return
        
        
    def SetFieldsSympy(self, t, z, E_i, E_r, E_t):
        self.t_sym = t
        self.z_sym = z
        self.E_i_sym = E_i
        self.E_r_sym = E_r
        self.E_t_sym = E_t
        self.E_i = lambdify((t, z), self.E_i_sym, 'numpy')
        self.E_r = lambdify((t, z), self.E_r_sym, 'numpy')
        self.E_t = lambdify((t, z), self.E_t_sym, 'numpy')
        return
        
        
    def GetJemSigeeSigmm(self, t):
        """ t: numpt vector
        """
        E_i, E_r, E_t = self.E_i, self.E_r, self.E_t
        eps_0, mu_0, eta_0 = self.eps_0, self.mu_0, self.eta_0

        Je_t = -(E_t(t, 0.0)-E_i(t, 0.0)+E_r(t, 0.0))/eta_0
        Jm_t = (-E_t(t, 0.0)+E_i(t, 0.0)+E_r(t, 0.0))

        ##TODO: prevent devide by zeros
        Sig_ee = (2.0*Je_t)/(E_t(t, 0.0)+E_i(t, 0.0)+E_r(t, 0.0))
        Sig_mm = (2.0*eta_0*Jm_t)/(E_t(t, 0.0)+E_i(t, 0.0)-E_r(t, 0.0))

        self.t = t
        self.Je_t = Je_t
        self.Jm_t = Jm_t
        self.Sig_ee = Sig_ee
        self.Sig_mm = Sig_mm
        return



    def TimeWindowConductivities(self, wnd_ee, wnd_mm):
        self.wnd_ee = wnd_ee
        self.wnd_mm = wnd_mm
        t = self.t
        wnd = wnd_ee
        self.Sig_ee = self.Sig_ee*np.logical_and(t>wnd[0], t<wnd[1])
        wnd = wnd_mm
        self.Sig_mm = self.Sig_mm*np.logical_and(t>wnd[0], t<wnd[1])



    def ClipConductivities(self, Sig_ee_max, Sig_mm_max):
        assert Sig_ee_max>0 and Sig_mm_max>0
        t = self.t
        Sig_ee, Sig_mm = self.Sig_ee, self.Sig_mm
        Sig_ee[np.logical_not(np.isfinite(Sig_ee))] = Sig_ee_max 
        Sig_mm[np.logical_not(np.isfinite(Sig_mm))] = Sig_mm_max
        
        Sig_min, Sig_max = -Sig_ee_max, Sig_ee_max
        self.Sig_ee =  Sig_ee*(np.logical_and(Sig_ee<=Sig_max,Sig_ee>=Sig_min)) +\
                         (Sig_ee>Sig_max)*Sig_max + (Sig_ee<Sig_min)*Sig_min
        Sig_min, Sig_max = -Sig_mm_max, Sig_mm_max
        self.Sig_mm =  Sig_mm*(np.logical_and(Sig_mm<=Sig_max,Sig_mm>=Sig_min)) +\
                         (Sig_mm>Sig_max)*Sig_max + (Sig_mm<Sig_min)*Sig_min



    def FindCondZeros(self):
        t = self.t
        Sig_ee = self.Sig_ee
        Sig_mm = self.Sig_mm
        Sig_ee_intervals = []
        Sig_mm_intervals = []
        
        Sig_ee_max = np.max(Sig_ee)
        Sig_mm_max = np.max(Sig_mm)
        
        N = len(t)
        for i in range(N-1):
            if Sig_ee[i]*Sig_ee[i+1]<=0.0:
                if Sig_ee[i+1]!=0:
                    if np.abs(Sig_ee[i+1]-Sig_ee[i])<Sig_ee_max/10.0:
                        Sig_ee_intervals.append([t[i], t[i+1]])
                if i>0 and Sig_ee[i]==0 and Sig_ee[i+1]==0 and Sig_ee[i-1]!=0:
                    if np.abs(Sig_ee[i-1]-Sig_ee[i])<Sig_ee_max/10.0:
                        Sig_ee_intervals.append([t[i-1], t[i]])
            if Sig_mm[i]*Sig_mm[i+1]<=0.0:
                if Sig_mm[i+1]!=0:
                    if np.abs(Sig_mm[i+1]-Sig_mm[i])<Sig_mm_max/10.0:
                        Sig_mm_intervals.append([t[i], t[i+1]])
                if i>0 and Sig_mm[i]==0 and Sig_mm[i+1]==0 and Sig_mm[i-1]!=0:
                    if np.abs(Sig_mm[i-1]-Sig_mm[i])<Sig_mm_max/10.0:
                        Sig_mm_intervals.append([t[i-1], t[i]])
        if Sig_ee[N-1]==0.0 and Sig_ee[N-2]!=0.0:
            if np.abs(Sig_ee[i+1]-Sig_ee[i])<Sig_ee_max/10.0:
                Sig_ee_intervals.append([t[N-1], t[N-2]])
        if Sig_mm[N-1]==0.0 and Sig_mm[N-2]!=0.0:
            if np.abs(Sig_mm[i+1]-Sig_mm[i])<Sig_mm_max/10.0:
                Sig_mm_intervals.append([t[N-1], t[N-2]])

        self.Sig_ee_intervals = Sig_ee_intervals
        self.Sig_mm_intervals = Sig_mm_intervals
        print('self.Sig_ee_intervals:', self.Sig_ee_intervals)
        print('self.Sig_mm_intervals:', self.Sig_mm_intervals)
        
        
        Sig_ee_roots = []
        Sig_mm_roots = []
        
        E_i, E_r, E_t = self.E_i, self.E_r, self.E_t
        eps_0, mu_0, eta_0 = self.eps_0, self.mu_0, self.eta_0
        t_0 = self.t[0]
        def f_Sig_ee(t):
            Je_t = -(E_t(t, 0.0)-E_i(t, 0.0)+E_r(t, 0.0))/eta_0
            Sig_ee_t = (2.0*Je_t)/(E_t(t, 0.0)+E_i(t, 0.0)+E_r(t, 0.0))
            return Sig_ee_t
            
        def f_Sig_mm(t):
            Jm_t = (-E_t(t, 0.0)+E_i(t, 0.0)+E_r(t, 0.0))
            Sig_mm_t = (2.0*eta_0*Jm_t)/(E_t(t, 0.0)+E_i(t, 0.0)-E_r(t, 0.0))
            return Sig_mm_t

        for i in range(len(Sig_ee_intervals)):
            a, b = Sig_ee_intervals[i]
            if f_Sig_ee(a)*f_Sig_ee(b)<=0.0:
                t_root = bisect(f_Sig_ee, a, b)
                Sig_ee_roots.append(t_root)

        for i in range(len(Sig_mm_intervals)):
            a, b = Sig_mm_intervals[i]
            t_root = None
            if f_Sig_mm(a)*f_Sig_mm(b)<=0.0:
                t_root = bisect(f_Sig_mm, a, b)
                Sig_mm_roots.append(t_root)

        self.Sig_ee_roots = Sig_ee_roots
        self.Sig_mm_roots = Sig_mm_roots


    def AdjustWindowToZeros(self):
        wnd_ee, wnd_mm = None, None
        if len(self.Sig_ee_roots)>=2:
            wnd_ee = [self.Sig_ee_roots[0], self.Sig_ee_roots[-1]]
        if len(self.Sig_mm_roots)>=2:
            wnd_mm = [self.Sig_mm_roots[0], self.Sig_mm_roots[-1]]
        self.TimeWindowConductivities(wnd_ee, wnd_mm)


    def SolveErEt(self):
        eps_0, mu_0, eta_0 = self.eps_0, self.mu_0, self.eta_0
        t = self.t
        E_i = self.E_i(t, 0.0)
        
        sig_ee = self.Sig_ee
        sig_mm = self.Sig_mm

        E_r = (-2.0*(eta_0**2*sig_ee - sig_mm)*E_i/(2.0*eta_0**2*sig_ee + \
            eta_0*sig_ee*sig_mm + 4.0*eta_0 + 2.0*sig_mm))
        
        E_t = (-eta_0*(sig_ee*sig_mm - 4.0)*E_i/(2.0*eta_0**2*sig_ee + \
            eta_0*sig_ee*sig_mm + 4.0*eta_0 + 2.0*sig_mm))        

        return [E_r, E_t]


        
        
        
        
        

