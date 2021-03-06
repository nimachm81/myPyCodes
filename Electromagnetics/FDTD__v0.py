## FDTD 1D 2D 3D

__all__ = ["FDTDYEE"]

import numpy as np
import inspect


class FDTDYEE():
    
    def __init__(self, r0, r1, dr, dt, N_is_even=None):
        assert len(r0)==len(r1)==len(dr)
        self.n_dim=len(r0)
    
        self.r0 = np.copy(r0)   #3 vectors
        self.r1 = np.copy(r1)
        self.W = self.r1-self.r0
        self.N = (self.W/dr).astype(int)
        if N_is_even!=None:
            assert len(N_is_even)==self.n_dim
            for i in range(self.n_dim):
                if N_is_even[i]==True:
                    if self.N[i]%2==1:
                        self.N[i] += 1
                else:
                    if self.N[i]%2==0:
                        self.N[i] += 1
                
        self.dr = self.W/self.N
        self.dt = dt

        self.dtype = float
        
        if self.n_dim==3:
            self.Nsx = self.N+np.array([0, 1, 1])
            self.Nsy = self.N+np.array([1, 0, 1])
            self.Nsz = self.N+np.array([1, 1, 0])

            self.Nfx = self.N+np.array([1, 0, 0])
            self.Nfy = self.N+np.array([0, 1, 0])
            self.Nfz = self.N+np.array([0, 0, 1])
        elif self.n_dim==2:
            self.Nsx = self.N+np.array([0, 1])
            self.Nsy = self.N+np.array([1, 0])
            self.Nsz = self.N+np.array([1, 1])

            self.Nfx = self.N+np.array([1, 0])
            self.Nfy = self.N+np.array([0, 1])
            self.Nfz = self.N+np.array([0, 0])
        elif self.n_dim==1:
            self.Nsx = self.N+np.array([0])
            self.Nsy = self.N+np.array([1])
            self.Nsz = self.N+np.array([1])

            self.Nfx = self.N+np.array([1])
            self.Nfy = self.N+np.array([0])
            self.Nfz = self.N+np.array([0])
        
        
        self.ViewPlanes = []

        self.cart_ind = None
        self.cart_ind_tot = None
        return


    def SetCartesianIndex(self, cart_ind=None, cart_ind_tot=None):
        """ cart_ind: np.array(n_dim)
        """
        self.cart_ind = cart_ind
        self.cart_ind_tot = cart_ind_tot


    def AllocateSideArr_list(self, str_list):
        for A_str in str_list:
            self.allocateSideArrays(A_str)

    def AllocateFaceArr_list(self, str_list):
        for A_str in str_list:
            self.allocateFaceArrays(A_str)

    def allocateSideArrays(self, A_str):
        cmd_x = "%sx = np.zeros(tuple(self.Nsx), dtype=self.dtype)"%A_str
        cmd_y = "%sy = np.zeros(tuple(self.Nsy), dtype=self.dtype)"%A_str
        cmd_z = "%sz = np.zeros(tuple(self.Nsz), dtype=self.dtype)"%A_str
        exec(cmd_x)
        exec(cmd_y)
        exec(cmd_z)
        
        cmd = "self.%s = [%sx, %sy, %sz]"%(A_str, A_str, A_str, A_str)
        exec(cmd)
        exec("print('self.%s[0].shape:', self.%s[0].shape)"%(A_str, A_str))

    def allocateFaceArrays(self, A_str):
        cmd_x = "%sx = np.zeros(tuple(self.Nfx), dtype=self.dtype)"%A_str
        cmd_y = "%sy = np.zeros(tuple(self.Nfy), dtype=self.dtype)"%A_str
        cmd_z = "%sz = np.zeros(tuple(self.Nfz), dtype=self.dtype)"%A_str
        exec(cmd_x)
        exec(cmd_y)
        exec(cmd_z)

        cmd = "self.%s = [%sx, %sy, %sz]"%(A_str, A_str, A_str, A_str)
        exec(cmd)
        exec("print('self.%s[0].shape:', self.%s[0].shape)"%(A_str, A_str))


    def UpdateSide_dAdt_CurlB_C(self, A, B, C_list=None, a=None, b=None, c_list=None):
        """ Wall indices for A are not updated, it will be later forced to 0
            a*d/dt A = b*curl(B) + c*C
            a, A, b, c, C are side elems
            B is face elem
            it starts with D then C then B
            A is on the left hand side.. the rest on the sight hand side
        """
        if a==None:
            a = np.ones(3)
        if b==None:
            b = np.ones(3)
        Nfx, Nfy, Nfz = self.Nfx, self.Nfy, self.Nfz
        
        ax, ay, az = a
        bx, by, bz = b

        Ax, Ay, Az = A
        Bx, By, Bz = B

        #print('side', '-'*50)
        #print(np.max(np.abs(Ax)), np.max(np.abs(Ay)), np.max(np.abs(Az)))
        #print(np.max(np.abs(A[0])), np.max(np.abs(A[1])), np.max(np.abs(A[2])))


        dx, dy, dz = [None]*3
        if self.n_dim==3:
            dx, dy, dz = self.dr[0], self.dr[1], self.dr[2]
        elif self.n_dim==2:
            dx, dy, dz = self.dr[0], self.dr[1], None
        elif self.n_dim==1:
            dx, dy, dz = self.dr[0], None, None
        dt = self.dt
        
        if C_list!=None:
            if c_list==None:
                c_list = [None]*len(C_list)
            for i in range(len(C_list)):
                C = C_list[i]
                c = c_list[i]
                if c==None:
                    c = np.ones(3)
                    
                Cx, Cy, Cz = C
                cx, cy, cz = c
                Ax += Cx*(cx*dt/ax)
                Ay += Cy*(cy*dt/ay)
                Az += Cz*(cz*dt/az)

        ax_, ay_, az_ = ax, ay, az
        bx_, by_, bz_ = bx, by, bz

        if self.n_dim==3:
            if type(ax)==np.ndarray:
                ax_ = ax[:, 1:-1, 1:-1]
            if type(bx)==np.ndarray:
                bx_ = bx[:, 1:-1, 1:-1]

            if type(ay)==np.ndarray:
                ay_ = ay[1:-1, :, 1:-1]
            if type(by)==np.ndarray:
                by_ = by[1:-1, :, 1:-1]

            if type(az)==np.ndarray:
                az_ = az[1:-1, 1:-1, :]
            if type(bz)==np.ndarray:
                bz_ = bz[1:-1, 1:-1, :]
        elif self.n_dim==2:
            if type(ax)==np.ndarray:
                ax_ = ax[:, 1:-1]
            if type(bx)==np.ndarray:
                bx_ = bx[:, 1:-1]

            if type(ay)==np.ndarray:
                ay_ = ay[1:-1, :]
            if type(by)==np.ndarray:
                by_ = by[1:-1, :]

            if type(az)==np.ndarray:
                az_ = az[1:-1, 1:-1]
            if type(bz)==np.ndarray:
                bz_ = bz[1:-1, 1:-1]
        elif self.n_dim==1:
            if type(ax)==np.ndarray:
                ax_ = ax[:]
            if type(bx)==np.ndarray:
                bx_ = bx[:]

            if type(ay)==np.ndarray:
                ay_ = ay[1:-1]
            if type(by)==np.ndarray:
                by_ = by[1:-1]

            if type(az)==np.ndarray:
                az_ = az[1:-1]
            if type(bz)==np.ndarray:
                bz_ = bz[1:-1]

        if self.n_dim==3:
            Ax[:, 1:-1, 1:-1] -= (By[:, 1:-1, 1:Nfy[2]] - By[:, 1:-1, 0:Nfy[2]-1])*(bx_*dt/(dz*ax_)) 
            Ax[:, 1:-1, 1:-1] += (Bz[:, 1:Nfz[1], 1:-1] - Bz[:, 0:Nfz[1]-1, 1:-1])*(bx_*dt/(dy*ax_)) 

            Ay[1:-1, :, 1:-1] -= (Bz[1:Nfz[0], :, 1:-1] - Bz[0:Nfz[0]-1, :, 1:-1])*(by_*dt/(dx*ay_)) 
            Ay[1:-1, :, 1:-1] += (Bx[1:-1, :, 1:Nfx[2]] - Bx[1:-1, :, 0:Nfx[2]-1])*(by_*dt/(dz*ay_)) 

            Az[1:-1, 1:-1, :] -= (Bx[1:-1, 1:Nfx[1], :] - Bx[1:-1, 0:Nfx[1]-1, :])*(bz_*dt/(dy*az_)) 
            Az[1:-1, 1:-1, :] += (By[1:Nfy[0], 1:-1, :] - By[0:Nfy[0]-1, 1:-1, :])*(bz_*dt/(dx*az_)) 
        elif self.n_dim==2:
            #Ax[:, 1:-1] -= (By[:, 1:-1    ] - By[:, 1:-1      ])*(bx_*dt/(dz*ax_)) 
            Ax[:, 1:-1] += (Bz[:, 1:Nfz[1]] - Bz[:, 0:Nfz[1]-1])*(bx_*dt/(dy*ax_)) 

            Ay[1:-1, :] -= (Bz[1:Nfz[0], :] - Bz[0:Nfz[0]-1, :])*(by_*dt/(dx*ay_)) 
            #Ay[1:-1, :] += (Bx[1:-1, :    ] - Bx[1:-1, :      ])*(by_*dt/(dz*ay_)) 

            Az[1:-1, 1:-1] -= (Bx[1:-1, 1:Nfx[1]] - Bx[1:-1, 0:Nfx[1]-1])*(bz_*dt/(dy*az_)) 
            Az[1:-1, 1:-1] += (By[1:Nfy[0], 1:-1] - By[0:Nfy[0]-1, 1:-1])*(bz_*dt/(dx*az_)) 
        elif self.n_dim==1:
            #Ax[:] -= (By[:] - By[:])*(bx_*dt/(dz*ax_)) 
            #Ax[:] += (Bz[:] - Bz[:])*(bx_*dt/(dy*ax_)) 

            Ay[1:-1] -= (Bz[1:Nfz[0]] - Bz[0:Nfz[0]-1])*(by_*dt/(dx*ay_)) 
            #Ay[1:-1] += (Bx[1:-1    ] - Bx[1:-1      ])*(by_*dt/(dz*ay_)) 

            #Az[1:-1] -= (Bx[1:-1    ] - Bx[1:-1      ])*(bz_*dt/(dy*az_)) 
            Az[1:-1] += (By[1:Nfy[0]] - By[0:Nfy[0]-1])*(bz_*dt/(dx*az_)) 
                

    
    def GetCurlFace(self, B, b=None):
        """ b*curl(B) B is Face element
            returns side element
        """
        if b==None:
            b = np.ones(3)

        Nfx, Nfy, Nfz = self.Nfx, self.Nfy, self.Nfz
        Bx, By, Bz = B
        bx, by, bz = b

        bx_, by_, bz_ = bx, by, bz

        if self.n_dim==3:
            if type(bx)==np.ndarray:
                bx_ = bx[:, 1:-1, 1:-1]
            if type(by)==np.ndarray:
                by_ = by[1:-1, :, 1:-1]
            if type(bz)==np.ndarray:
                bz_ = bz[1:-1, 1:-1, :]
        elif self.n_dim==2:
            if type(bx)==np.ndarray:
                bx_ = bx[:, 1:-1]
            if type(by)==np.ndarray:
                by_ = by[1:-1, :]
            if type(bz)==np.ndarray:
                bz_ = bz[1:-1, 1:-1]
        elif self.n_dim==1:
            if type(bx)==np.ndarray:
                bx_ = bx[:, 1:-1]
            if type(by)==np.ndarray:
                by_ = by[1:-1, :]
            if type(bz)==np.ndarray:
                bz_ = bz[1:-1, 1:-1]

        dx, dy, dz = [None]*3
        if self.n_dim==3:
            dx, dy, dz = self.dr[0], self.dr[1], self.dr[2]
        elif self.n_dim==2:
            dx, dy, dz = self.dr[0], self.dr[1], None
        elif self.n_dim==1:
            dx, dy, dz = self.dr[0], None, None

        Ax = np.zeros(self.Nsx)
        Ay = np.zeros(self.Nsy)
        Az = np.zeros(self.Nsz)
        
        if self.n_dim==3:
            Ax[:, 1:-1, 1:-1] -= (By[:, 1:-1, 1:Nfy[2]] - By[:, 1:-1, 0:Nfy[2]-1])*(bx_/dz) 
            Ax[:, 1:-1, 1:-1] += (Bz[:, 1:Nfz[1], 1:-1] - Bz[:, 0:Nfz[1]-1, 1:-1])*(bx_/dy) 

            Ay[1:-1, :, 1:-1] -= (Bz[1:Nfz[0], :, 1:-1] - Bz[0:Nfz[0]-1, :, 1:-1])*(by_/dx) 
            Ay[1:-1, :, 1:-1] += (Bx[1:-1, :, 1:Nfx[2]] - Bx[1:-1, :, 0:Nfx[2]-1])*(by_/dz) 

            Az[1:-1, 1:-1, :] -= (Bx[1:-1, 1:Nfx[1], :] - Bx[1:-1, 0:Nfx[1]-1, :])*(bz_/dy) 
            Az[1:-1, 1:-1, :] += (By[1:Nfy[0], 1:-1, :] - By[0:Nfy[0]-1, 1:-1, :])*(bz_/dx)
        elif self.n_dim==2:
            #Ax[:, 1:-1] -= (By[:, 1:-1    ] - By[:, 1:-1      ])*(bx_/dz) 
            Ax[:, 1:-1] += (Bz[:, 1:Nfz[1]] - Bz[:, 0:Nfz[1]-1])*(bx_/dy) 

            Ay[1:-1, :] -= (Bz[1:Nfz[0], :] - Bz[0:Nfz[0]-1, :])*(by_/dx) 
            #Ay[1:-1, :] += (Bx[1:-1, :    ] - Bx[1:-1, :      ])*(by_/dz) 

            Az[1:-1, 1:-1] -= (Bx[1:-1, 1:Nfx[1]] - Bx[1:-1, 0:Nfx[1]-1])*(bz_/dy) 
            Az[1:-1, 1:-1] += (By[1:Nfy[0], 1:-1] - By[0:Nfy[0]-1, 1:-1])*(bz_/dx)
        elif self.n_dim==1:
            #Ax[:] -= (By[:] - By[:])*(bx_/dz) 
            #Ax[:] += (Bz[:] - Bz[:])*(bx_/dy) 

            Ay[1:-1] -= (Bz[1:Nfz[0]] - Bz[0:Nfz[0]-1])*(by_/dx) 
            #Ay[1:-1] += (Bx[1:-1    ] - Bx[1:-1      ])*(by_/dz) 

            #Az[1:-1] -= (Bx[1:-1    ] - Bx[1:-1      ])*(bz_/dy) 
            Az[1:-1] += (By[1:Nfy[0]] - By[0:Nfy[0]-1])*(bz_/dx)
        
        return [Ax, Ay, Az]

        
    def UpdateFace_dAdt_CurlB_C(self, A, B, C_list=None, a=None, b=None, c_list=None):
        """ Wall indices for A are not updated, it will be later forced to 0
            a*d/dt A = b*curl(B) + c*C
            a, A, b, c, C are face elems
            B is side elem
        """
        if a==None:
            a = np.ones(3)
        if b==None:
            b = np.ones(3)

        Nsx, Nsy, Nsz = self.Nsx, self.Nsy, self.Nsz
        ax, ay, az = a
        bx, by, bz = b

        Ax, Ay, Az = A
        Bx, By, Bz = B
        dx, dy, dz = [None]*3
        if self.n_dim==3:
            dx, dy, dz = self.dr[0], self.dr[1], self.dr[2]
        elif self.n_dim==2:
            dx, dy, dz = self.dr[0], self.dr[1], None
        elif self.n_dim==1:
            dx, dy, dz = self.dr[0], None, None
        dt = self.dt

        if C_list!=None:
            if c_list==None:
                c_list = [None]*len(C_list)
            for i in range(len(C_list)):
                C = C_list[i]
                c = c_list[i]
                if c==None:
                    c = np.ones(3)
                    
                Cx, Cy, Cz = C
                cx, cy, cz = c
                Ax += Cx*(cx*dt/ax)
                Ay += Cy*(cy*dt/ay)
                Az += Cz*(cz*dt/az)

        if self.n_dim==3:
            Ax -= (By[:, :, 1:Nsy[2]] - By[:, :, 0:Nsy[2]-1])*(bx*dt/(dz*ax)) 
            Ax += (Bz[:, 1:Nsz[1], :] - Bz[:, 0:Nsz[1]-1, :])*(bx*dt/(dy*ax)) 

            Ay -= (Bz[1:Nsz[0], :, :] - Bz[0:Nsz[0]-1, :, :])*(by*dt/(dx*ay)) 
            Ay += (Bx[:, :, 1:Nsx[2]] - Bx[:, :, 0:Nsx[2]-1])*(by*dt/(dz*ay)) 

            Az -= (Bx[:, 1:Nsx[1], :] - Bx[:, 0:Nsx[1]-1, :])*(bz*dt/(dy*az)) 
            Az += (By[1:Nsy[0], :, :] - By[0:Nsy[0]-1, :, :])*(bz*dt/(dx*az)) 
        elif self.n_dim==2:
            #Ax -= (By[:, :       ] - By[:, :         ])*(bx*dt/(dz*ax)) 
            Ax += (Bz[:, 1:Nsz[1]] - Bz[:, 0:Nsz[1]-1])*(bx*dt/(dy*ax)) 

            Ay -= (Bz[1:Nsz[0], :] - Bz[0:Nsz[0]-1, :])*(by*dt/(dx*ay)) 
            #Ay += (Bx[:, :       ] - Bx[:, :         ])*(by*dt/(dz*ay)) 

            Az -= (Bx[:, 1:Nsx[1]] - Bx[:, 0:Nsx[1]-1])*(bz*dt/(dy*az)) 
            Az += (By[1:Nsy[0], :] - By[0:Nsy[0]-1, :])*(bz*dt/(dx*az)) 
        elif self.n_dim==1:
            #Ax -= (By[:] - By[:])*(bx*dt/(dz*ax)) 
            #Ax += (Bz[:] - Bz[:])*(bx*dt/(dy*ax)) 

            Ay -= (Bz[1:Nsz[0]] - Bz[0:Nsz[0]-1])*(by*dt/(dx*ay)) 
            #Ay += (Bx[:       ] - Bx[:         ])*(by*dt/(dz*ay)) 

            #Az -= (Bx[:       ] - Bx[:         ])*(bz*dt/(dy*az)) 
            Az += (By[1:Nsy[0]] - By[0:Nsy[0]-1])*(bz*dt/(dx*az)) 


    def GetCurlSide(self, B, b=None):
        """ b*curl(B) B is side element
            returns face element
        """
        if b==None:
            b = np.ones(3)
            
        Nsx, Nsy, Nsz = self.Nsx, self.Nsy, self.Nsz
        Bx, By, Bz = B
        bx, by, bz = b

        dx, dy, dz = [None]*3
        if self.n_dim==3:
            dx, dy, dz = self.dr[0], self.dr[1], self.dr[2]
        elif self.n_dim==2:
            dx, dy, dz = self.dr[0], self.dr[1], None
        elif self.n_dim==1:
            dx, dy, dz = self.dr[0], None, None
        dt = self.dt

        Ax, Ay, Az = np.zeros(self.Nfx), np.zeros(self.Nfy), np.zeros(self.Nfz)
        if self.n_dim==3:
            Ax -= (By[:, :, 1:Nsy[2]] - By[:, :, 0:Nsy[2]-1])*(bx/dz) 
            Ax += (Bz[:, 1:Nsz[1], :] - Bz[:, 0:Nsz[1]-1, :])*(bx/dy) 

            Ay -= (Bz[1:Nsz[0], :, :] - Bz[0:Nsz[0]-1, :, :])*(by/dx) 
            Ay += (Bx[:, :, 1:Nsx[2]] - Bx[:, :, 0:Nsx[2]-1])*(by/dz) 

            Az -= (Bx[:, 1:Nsx[1], :] - Bx[:, 0:Nsx[1]-1, :])*(bz/dy) 
            Az += (By[1:Nsy[0], :, :] - By[0:Nsy[0]-1, :, :])*(bz/dx) 
        if self.n_dim==2:
            #Ax -= (By[:, :       ] - By[:, :         ])*(bx/dz) 
            Ax += (Bz[:, 1:Nsz[1]] - Bz[:, 0:Nsz[1]-1])*(bx/dy) 

            Ay -= (Bz[1:Nsz[0], :] - Bz[0:Nsz[0]-1, :])*(by/dx) 
            #Ay += (Bx[:, :       ] - Bx[:, :         ])*(by/dz) 

            Az -= (Bx[:, 1:Nsx[1]] - Bx[:, 0:Nsx[1]-1])*(bz/dy) 
            Az += (By[1:Nsy[0], :] - By[0:Nsy[0]-1, :])*(bz/dx) 
        if self.n_dim==1:
            #Ax -= (By[:] - By[:])*(bx/dz) 
            #Ax += (Bz[:] - Bz[:])*(bx/dy) 

            Ay -= (Bz[1:Nsz[0]] - Bz[0:Nsz[0]-1])*(by/dx) 
            #Ay += (Bx[:       ] - Bx[:         ])*(by/dz) 

            #Az -= (Bx[:       ] - Bx[:         ])*(bz/dy) 
            Az += (By[1:Nsy[0]] - By[0:Nsy[0]-1])*(bz/dx) 
        return [Ax, Ay, Az]


    def Update_adAdt_bB(self, A, B, C_list=None, a=None, b=None, c_list=None):
        """ a*d/dt A = b*B
            all are side/face elems
        """
        if a==None:
            a = np.ones(3)
        if b==None:
            b = np.ones(3)

        ax, ay, az = a
        bx, by, bz = b

        Ax, Ay, Az = A
        Bx, By, Bz = B
        
        dt = self.dt
        
        if C_list!=None:
            if c_list==None:
                c_list = [None]*len(C_list)
            for i in range(len(C_list)):
                C = C_list[i]
                c = c_list[i]
                if c==None:
                    c = np.ones(3)
                    
                Cx, Cy, Cz = C
                cx, cy, cz = c
                Ax += Cx*(cx*dt/ax)
                Ay += Cy*(cy*dt/ay)
                Az += Cz*(cz*dt/az)

        Ax += Bx*(bx*dt/ax)
        Ay += By*(by*dt/ay)
        Az += Bz*(bz*dt/az)
        

    def Update_aA_bB(self, A, B, C_list=None, a=None, b=None, c_list=None):
        """ a*A = b*B
            all are side/face elems
        """
        if a==None:
            a = np.ones(3)
        if b==None:
            b = np.ones(3)

        ax, ay, az = a
        bx, by, bz = b

        Bx, By, Bz = B
                
        A[0] = Bx*(bx/ax)
        A[1] = By*(by/ay)
        A[2] = Bz*(bz/az)

        Ax, Ay, Az = A

        if C_list!=None:
            if c_list==None:
                c_list = [None]*len(C_list)
            for i in range(len(C_list)):
                C = C_list[i]
                c = c_list[i]
                if c==None:
                    c = np.ones(3)
                    
                Cx, Cy, Cz = C
                cx, cy, cz = c
                Ax += Cx*(cx/ax)
                Ay += Cy*(cy/ay)
                Az += Cz*(cz/az)


    def ResetSideWalls(self, A):
        """ A is side element
        """
        Nsx, Nsy, Nsz = self.Nsx, self.Nsy, self.Nsz
        Ax, Ay, Az = A
        
        if self.n_dim==3:
            Ax[:       , 0       , :       ] = 0.0
            Ax[:       , Nsx[1]-1, :       ] = 0.0
            Ax[:       , :       , 0       ] = 0.0
            Ax[:       , :       , Nsx[2]-1] = 0.0
            
            Ay[0       , :       , :       ] = 0.0
            Ay[Nsy[0]-1, :       , :       ] = 0.0
            Ay[:       , :       , 0       ] = 0.0
            Ay[:       , :       , Nsy[2]-1] = 0.0

            Az[0       , :       , :       ] = 0.0
            Az[Nsz[0]-1, :       , :       ] = 0.0
            Az[:       , 0       , :       ] = 0.0
            Az[:       , Nsz[1]-1, :       ] = 0.0
        elif self.n_dim==2:
            Ax[:       , 0       ] = 0.0
            Ax[:       , Nsx[1]-1] = 0.0
            
            Ay[0       , :       ] = 0.0
            Ay[Nsy[0]-1, :       ] = 0.0

            Az[0       , :       ] = 0.0
            Az[Nsz[0]-1, :       ] = 0.0
            Az[:       , 0       ] = 0.0
            Az[:       , Nsz[1]-1] = 0.0
        elif self.n_dim==1:            
            Ay[0       ] = 0.0
            Ay[Nsy[0]-1] = 0.0

            Az[0       ] = 0.0
            Az[Nsz[0]-1] = 0.0


    def SetSpatialPointsSide(self):     
        if self.n_dim==3:   
            dx, dy, dz = self.dr[0], self.dr[1], self.dr[2]
            x0, y0, z0 = self.r0[0], self.r0[1], self.r0[2]
            x1, y1, z1 = self.r1[0], self.r1[1], self.r1[2]
            nx, ny, nz = self.N[0],  self.N[1],  self.N[2]
            
            x_1d = np.linspace(x0+dx/2, x1-dx/2, nx, endpoint=True)
            y_1d = np.linspace(y0+dy/2, y1-dy/2, ny, endpoint=True)
            z_1d = np.linspace(z0+dz/2, z1-dz/2, nz, endpoint=True)

            x_1dp = np.linspace(x0, x1, nx+1, endpoint=True)
            y_1dp = np.linspace(y0, y1, ny+1, endpoint=True)
            z_1dp = np.linspace(z0, z1, nz+1, endpoint=True)
             
            rx_x, rx_y, rx_z = np.meshgrid(x_1d, y_1dp, z_1dp, indexing='ij')
            self.rsx = [rx_x, rx_y, rx_z]

            ry_x, ry_y, ry_z = np.meshgrid(x_1dp, y_1d, z_1dp, indexing='ij')
            self.rsy = [ry_x, ry_y, ry_z]

            rz_x, rz_y, rz_z = np.meshgrid(x_1dp, y_1dp, z_1d, indexing='ij')
            self.rsz = [rz_x, rz_y, rz_z]
        elif self.n_dim==2:   
            dx, dy, dz = self.dr[0], self.dr[1], None
            x0, y0, z0 = self.r0[0], self.r0[1], None
            x1, y1, z1 = self.r1[0], self.r1[1], None
            nx, ny, nz = self.N[0],  self.N[1],  None
            
            x_1d = np.linspace(x0+dx/2, x1-dx/2, nx, endpoint=True)
            y_1d = np.linspace(y0+dy/2, y1-dy/2, ny, endpoint=True)

            x_1dp = np.linspace(x0, x1, nx+1, endpoint=True)
            y_1dp = np.linspace(y0, y1, ny+1, endpoint=True)
             
            rx_x, rx_y = np.meshgrid(x_1d, y_1dp, indexing='ij')
            self.rsx = [rx_x, rx_y]

            ry_x, ry_y = np.meshgrid(x_1dp, y_1d, indexing='ij')
            self.rsy = [ry_x, ry_y]

            rz_x, rz_y = np.meshgrid(x_1dp, y_1dp, indexing='ij')
            self.rsz = [rz_x, rz_y]
        elif self.n_dim==1:   
            dx, dy, dz = self.dr[0], None, None
            x0, y0, z0 = self.r0[0], None, None
            x1, y1, z1 = self.r1[0], None, None
            nx, ny, nz = self.N[0],  None,  None
            
            x_1d = np.linspace(x0+dx/2, x1-dx/2, nx, endpoint=True)

            x_1dp = np.linspace(x0, x1, nx+1, endpoint=True)
             
            rx_x = x_1d
            self.rsx = [rx_x]

            ry_x = x_1dp
            self.rsy = [ry_x]

            rz_x = x_1dp
            self.rsz = [rz_x]


    def SetSpatialPointsFace(self):        
        if self.n_dim==3:
            dx, dy, dz = self.dr[0], self.dr[1], self.dr[2]
            x0, y0, z0 = self.r0[0], self.r0[1], self.r0[2]
            x1, y1, z1 = self.r1[0], self.r1[1], self.r1[2]
            nx, ny, nz = self.N[0],  self.N[1],  self.N[2]
            
            x_1d = np.linspace(x0+dx/2, x1-dx/2, nx, endpoint=True)
            y_1d = np.linspace(y0+dy/2, y1-dy/2, ny, endpoint=True)
            z_1d = np.linspace(z0+dz/2, z1-dz/2, nz, endpoint=True)

            x_1dp = np.linspace(x0, x1, nx+1, endpoint=True)
            y_1dp = np.linspace(y0, y1, ny+1, endpoint=True)
            z_1dp = np.linspace(z0, z1, nz+1, endpoint=True)
            
            rx_x, rx_y, rx_z = np.meshgrid(x_1dp, y_1d, z_1d, indexing='ij')
            self.rfx = [rx_x, rx_y, rx_z]
             
            ry_x, ry_y, ry_z = np.meshgrid(x_1d, y_1dp, z_1d, indexing='ij')
            self.rfy = [ry_x, ry_y, ry_z]

            rz_x, rz_y, rz_z = np.meshgrid(x_1d, y_1d, z_1dp, indexing='ij')
            self.rfz = [rz_x, rz_y, rz_z]
        elif self.n_dim==2:
            dx, dy, dz = self.dr[0], self.dr[1], None
            x0, y0, z0 = self.r0[0], self.r0[1], None
            x1, y1, z1 = self.r1[0], self.r1[1], None
            nx, ny, nz = self.N[0],  self.N[1],  None
            
            x_1d = np.linspace(x0+dx/2, x1-dx/2, nx, endpoint=True)
            y_1d = np.linspace(y0+dy/2, y1-dy/2, ny, endpoint=True)

            x_1dp = np.linspace(x0, x1, nx+1, endpoint=True)
            y_1dp = np.linspace(y0, y1, ny+1, endpoint=True)
            
            rx_x, rx_y = np.meshgrid(x_1dp, y_1d, indexing='ij')
            self.rfx = [rx_x, rx_y]
             
            ry_x, ry_y = np.meshgrid(x_1d, y_1dp, indexing='ij')
            self.rfy = [ry_x, ry_y]

            rz_x, rz_y = np.meshgrid(x_1d, y_1d, indexing='ij')
            self.rfz = [rz_x, rz_y]
        elif self.n_dim==1:
            dx, dy, dz = self.dr[0], None, None
            x0, y0, z0 = self.r0[0], None, None
            x1, y1, z1 = self.r1[0], None, None
            nx, ny, nz = self.N[0],  None,  None
            
            x_1d = np.linspace(x0+dx/2, x1-dx/2, nx, endpoint=True)

            x_1dp = np.linspace(x0, x1, nx+1, endpoint=True)
            
            rx_x = x_1dp
            self.rfx = [rx_x]
             
            ry_x = x_1d
            self.rfy = [ry_x]

            rz_x = x_1d
            self.rfz = [rz_x]
    

    def UpdateFuncSide_space(self, A, f):
        fx, fy, fz = f
        rsx, rsy, rsz = self.rsx, self.rsy, self.rsz
        if self.n_dim==1:
            rsx, rsy, rsz = self.rsx[0], self.rsy[0], self.rsz[0]
        
        if fx!=None:                
            A[0] = fx(rsx)
        if fy!=None:                
            A[1] = fy(rsy)
        if fz!=None:                
            A[2] = fz(rsz)
        A_copy = [np.copy(A[0]), np.copy(A[1]), np.copy(A[2])]
        return A_copy
            

    def UpdateFuncSide_spacetime(self, A, f, t):
        fx, fy, fz = f
        rsx, rsy, rsz = self.rsx, self.rsy, self.rsz
        if self.n_dim==1:
            rsx, rsy, rsz = self.rsx[0], self.rsy[0], self.rsz[0]
        
        if fx!=None:                
            A[0] = fx(rsx, t)
        if fy!=None:                
            A[1] = fy(rsy, t)
        if fz!=None:                
            A[2] = fz(rsz, t)
        A_copy = [np.copy(A[0]), np.copy(A[1]), np.copy(A[2])]
        return A_copy

    
    def UpdateFuncFace_space(self, A, f):
        fx, fy, fz = f
        rfx, rfy, rfz = self.rfx, self.rfy, self.rfz
        if self.n_dim==1:
            rfx, rfy, rfz = self.rfx[0], self.rfy[0], self.rfz[0]

        if fx!=None:                
            A[0] = fx(rfx)
        if fy!=None:                
            A[1] = fy(rfy)
        if fz!=None:                
            A[2] = fz(rfz)
        A_copy = [np.copy(A[0]), np.copy(A[1]), np.copy(A[2])]
        return A_copy
        

    def UpdateFuncFace_spacetime(self, A, f, t):
        fx, fy, fz = f
        rfx, rfy, rfz = self.rfx, self.rfy, self.rfz
        if self.n_dim==1:
            rfx, rfy, rfz = self.rfx[0], self.rfy[0], self.rfz[0]

        if fx!=None:                
            A[0] = fx(rfx, t)
        if fy!=None:                
            A[1] = fy(rfy, t)
        if fz!=None:                
            A[2] = fz(rfz, t)
        A_copy = [np.copy(A[0]), np.copy(A[1]), np.copy(A[2])]
        return A_copy


    def UpdateSeperableFunc_Time(self, A, A_0, f, t):
        A[0] = f(t)*A_0[0]
        A[1] = f(t)*A_0[1]
        A[2] = f(t)*A_0[2]

        
    def ComplexEpsToEpsSigmaE(self, eps_c, omega):
        eps = np.real(eps_c)
        sig = -np.imag(eps_c)*omega
        return [eps, sig]


    def ComplexMuToMuSigmaM(self, mu_c, omega):
        mu = np.real(mu_c)
        sig = -np.imag(mu_c)*omega
        return [eps, sig]


    def GetWallsAllDic__ver_0(self, dw, s):
        if self.n_dim==3:
            dx, dy, dz = dw[0], dw[1], dw[2]
            sx, sy, sz = s
            walls = {('x', 'n'):[sx, dx], ('x', 'p'):[sx, dx],
                     ('y', 'n'):[sy, dy], ('y', 'p'):[sy, dy],
                     ('z', 'n'):[sz, dz], ('z', 'p'):[sz, dz]}
            return walls
        elif self.n_dim==2:
            dx, dy = dw[0], dw[1]
            sx, sy = s
            walls = {('x', 'n'):[sx, dx], ('x', 'p'):[sx, dx],
                     ('y', 'n'):[sy, dy], ('y', 'p'):[sy, dy]}
            return walls
        elif self.n_dim==1:
            dx = dw[0]
            sx = s
            walls = {('x', 'n'):[sx, dx], ('x', 'p'):[sx, dx]}
            return walls


    def GetWallsAllDic(self, dw, s):
        if self.n_dim==3:
            dx, dy, dz = dw[0], dw[1], dw[2]
            sx, sy, sz = s
            walls = {}
            if dx>0.0:
                walls[('x', 'n')] = [sx, dx]
                walls[('x', 'p')] = [sx, dx]
            if dy>0.0:
                walls[('y', 'n')] = [sy, dy]
                walls[('y', 'p')] = [sy, dy]
            if dz>0.0:
                walls[('z', 'n')] = [sz, dz]
                walls[('z', 'p')] = [sz, dz]
            return walls
        elif self.n_dim==2:
            dx, dy = dw[0], dw[1]
            sx, sy = s
            walls = {}
            if dx>0.0:
                walls[('x', 'n')] = [sx, dx]
                walls[('x', 'p')] = [sx, dx]
            if dy>0.0:
                walls[('y', 'n')] = [sy, dy]
                walls[('y', 'p')] = [sy, dy]
            return walls
        elif self.n_dim==1:
            dx = dw[0]
            sx = s
            walls = {}
            if dx>0.0:
                walls[('x', 'n')] = [sx, dx]
                walls[('x', 'p')] = [sx, dx]
            return walls

    def GetUPMLFactor(self, walls, eps_or_s='s', side_or_face='side'):
        """ walls:{('x', 'n'):[sx, dx], ('y', 'p'):[sy, dy]..}
            ('x', 'n'):[sx, dx] --> plane normal to x direction, negative side, 
                stretch factor sx, thickness dx
            eps_or_s='eps'/'s'
            returns boxes with complex upml permittivity/permeability factors
            if eps_or_s=='s' returns where sx/sy/sz are non-zero
        """
        Nx, Ny, Nz = [None]*3
        if side_or_face=='side':
            Nx, Ny, Nz = self.Nsx, self.Nsy, self.Nsz
        else:
            assert side_or_face=='face'
            Nx, Ny, Nz = self.Nfx, self.Nfy, self.Nfz
            
        if eps_or_s=='eps':
            eps_xx = np.ones(Nx)
            eps_yy = np.ones(Ny)
            eps_zz = np.ones(Nz)
            for wall_dir in walls:
                w, n_dir = wall_dir
                f = None
                if w=='x':
                    r1 = np.copy(self.r1)
                    sx, dx = walls[wall_dir]
                    if n_dir=='n':
                        r1[0] = self.r0[0] + dx
                        r0 = np.copy(self.r0)
                        f_xx = self.FunctionBox(r0, r1, a=1.0/sx, b=1.0)    
                        f_yy = self.FunctionBox(r0, r1, a=sx, b=1.0)    
                        f_zz = self.FunctionBox(r0, r1, a=sx, b=1.0)    
                        f = [f_xx, f_yy, f_zz]
                    else:
                        assert n_dir=='p'
                        r0 = np.copy(self.r0)
                        r0[0] = self.r1[0] - dx
                        f_xx = self.FunctionBox(r0, r1, a=1.0/sx, b=1.0)    
                        f_yy = self.FunctionBox(r0, r1, a=sx, b=1.0)    
                        f_zz = self.FunctionBox(r0, r1, a=sx, b=1.0)    
                        f = [f_xx, f_yy, f_zz]
                elif w=='y':
                    assert self.n_dim>=2
                    r1 = np.copy(self.r1)
                    sy, dy = walls[wall_dir]
                    if n_dir=='n':
                        r1[1] = self.r0[1] + dy
                        r0 = np.copy(self.r0)
                        f_xx = self.FunctionBox(r0, r1, a=sy, b=1.0)    
                        f_yy = self.FunctionBox(r0, r1, a=1.0/sy, b=1.0)    
                        f_zz = self.FunctionBox(r0, r1, a=sy, b=1.0)    
                        f = [f_xx, f_yy, f_zz]
                    else:
                        assert n_dir=='p'
                        r0 = np.copy(self.r0)
                        r0[1] = self.r1[1] - dy
                        f_xx = self.FunctionBox(r0, r1, a=sy, b=1.0)    
                        f_yy = self.FunctionBox(r0, r1, a=1.0/sy, b=1.0)    
                        f_zz = self.FunctionBox(r0, r1, a=sy, b=1.0)    
                        f = [f_xx, f_yy, f_zz]
                else:
                    assert w=='z'
                    assert self.n_dim==3
                    r1 = np.copy(self.r1)
                    sz, dz = walls[wall_dir]
                    if n_dir=='n':
                        r1[2] = self.r0[2] + dz
                        r0 = np.copy(self.r0)
                        f_xx = self.FunctionBox(r0, r1, a=sz, b=1.0)    
                        f_yy = self.FunctionBox(r0, r1, a=sz, b=1.0)    
                        f_zz = self.FunctionBox(r0, r1, a=1.0/sz, b=1.0)    
                        f = [f_xx, f_yy, f_zz]
                    else:
                        assert n_dir=='p'
                        r0 = np.copy(self.r0)
                        r0[2] = self.r1[2] - dz
                        f_xx = self.FunctionBox(r0, r1, a=sz, b=1.0)    
                        f_yy = self.FunctionBox(r0, r1, a=sz, b=1.0)    
                        f_zz = self.FunctionBox(r0, r1, a=1.0/sz, b=1.0)    
                        f = [f_xx, f_yy, f_zz]
                
                A = [None]*3
                if side_or_face=='side':
                    self.UpdateFuncSide_space(A, f)
                else:
                    self.UpdateFuncFace_space(A, f)
                eps_xx *= A[0]
                eps_yy *= A[1]
                eps_zz *= A[2]
            return [eps_xx, eps_yy, eps_zz]
        else:
            assert eps_or_s == 's'
            sx_arr, sy_arr, sz_arr = [None]*3, [None]*3, [None]*3
            sx_arr[0] = np.ones(Nx, dtype=complex)   ## s_x for sides x
            sx_arr[1] = np.ones(Ny, dtype=complex)   ## s_x for sides y
            sx_arr[2] = np.ones(Nz, dtype=complex)

            sy_arr[0] = np.ones(Nx, dtype=complex)   ## s_y for sides x
            sy_arr[1] = np.ones(Ny, dtype=complex)   ## s_y for sides y
            sy_arr[2] = np.ones(Nz, dtype=complex)

            sz_arr[0] = np.ones(Nx, dtype=complex)   ## s_z for sides x
            sz_arr[1] = np.ones(Ny, dtype=complex)   ## s_z for sides y
            sz_arr[2] = np.ones(Nz, dtype=complex)

            for wall_dir in walls:
                w, n_dir = wall_dir
                f = None
                if w=='x':
                    r1 = np.copy(self.r1)
                    sx, dx = walls[wall_dir]
                    if n_dir=='n':
                        r1[0] = self.r0[0] + dx
                        r0 = np.copy(self.r0)
                        f = self.FunctionBox(r0, r1, a=sx-1.0, b=0.0)    
                    else:
                        assert n_dir=='p'
                        r0 = np.copy(self.r0)
                        r0[0] = self.r1[0] - dx
                        f = self.FunctionBox(r0, r1, a=sx-1.0, b=0.0)    
                    A = [None]*3
                    if side_or_face=='side':
                        self.UpdateFuncSide_space(A, [f, f, f])
                    else:
                        self.UpdateFuncFace_space(A, [f, f, f])
                    sx_arr[0] += A[0]
                    sx_arr[1] += A[1]
                    sx_arr[2] += A[2]
                elif w=='y':
                    assert self.n_dim>=2
                    r1 = np.copy(self.r1)
                    sy, dy = walls[wall_dir]
                    if n_dir=='n':
                        r1[1] = self.r0[1] + dy
                        r0 = np.copy(self.r0)
                        f = self.FunctionBox(r0, r1, a=sy-1.0, b=0.0)    
                    else:
                        assert n_dir=='p'
                        r0 = np.copy(self.r0)
                        r0[1] = self.r1[1] - dy
                        f = self.FunctionBox(r0, r1, a=sy-1.0, b=0.0)    
                    A = [None]*3
                    if side_or_face=='side':
                        self.UpdateFuncSide_space(A, [f, f, f])
                    else:
                        self.UpdateFuncFace_space(A, [f, f, f])
                    sy_arr[0] += A[0]
                    sy_arr[1] += A[1]
                    sy_arr[2] += A[2]
                else:
                    assert w=='z'
                    assert self.n_dim==3
                    r1 = np.copy(self.r1)
                    sz, dz = walls[wall_dir]
                    if n_dir=='n':
                        r1[2] = self.r0[2] + dz
                        r0 = np.copy(self.r0)
                        f = self.FunctionBox(r0, r1, a=sz-1.0, b=0.0)    
                    else:
                        assert n_dir=='p'
                        r0 = np.copy(self.r0)
                        r0[2] = self.r1[2] - dz
                        f = self.FunctionBox(r0, r1, a=sz-1.0, b=0.0)    
                    A = [None]*3
                    if side_or_face=='side':
                        self.UpdateFuncSide_space(A, [f, f, f])
                    else:
                        self.UpdateFuncFace_space(A, [f, f, f])
                    sz_arr[0] += A[0]
                    sz_arr[1] += A[1]
                    sz_arr[2] += A[2]
            return [sx_arr, sy_arr, sz_arr]
                    


    def FunctionSphere(self, r0, rad, a, b=0.0):
        """ r0: center
            rad: radius
            a: amplitude on shpere
            b: amplitude outside sphere
        """
        if self.n_dim==3:
            x0, y0, z0 = r0[0], r0[1], r0[2]        
            def f(r):
                x, y, z = r
                R2 = (x-x0)**2 + (y-y0)**2 + (z-z0)**2
                isin = (R2 <= rad**2)
                isout = np.logical_not(isin)
                res = None
                if inspect.isfunction(a):
                    res = a(r)*isin
                else:
                    res = a*isin
                if inspect.isfunction(b):
                    res += b(r)*isout
                else:
                    res += b*isout
                return res
            return f
        if self.n_dim==2:
            x0, y0 = r0[0], r0[1]        
            def f(r):
                x, y = r
                R2 = (x-x0)**2 + (y-y0)**2
                isin = (R2 <= rad**2)
                isout = np.logical_not(isin)
                res = None
                if inspect.isfunction(a):
                    res = a(r)*isin
                else:
                    res = a*isin
                if inspect.isfunction(b):
                    res += b(r)*isout
                else:
                    res += b*isout
                return res
            return f
        if self.n_dim==1:
            x0 = r0[0]        
            def f(r):
                x = r   ##TODO: make consistent --> x=r[0]
                R2 = (x-x0)**2
                isin = (R2 <= rad**2)
                isout = np.logical_not(isin)
                res = None
                if inspect.isfunction(a):
                    res = a(r)*isin
                else:
                    res = a*isin
                if inspect.isfunction(b):
                    res += b(r)*isout
                else:
                    res += b*isout
                return res
            return f

    def FunctionBox(self, r0, r1, a, b=0.0):
        """ r0: lowe left corner
            r1: upper right corner
            a: amplitude
        """
        if self.n_dim==3:
            x0, y0, z0 = r0[0], r0[1], r0[2]        
            x1, y1, z1 = r1[0], r1[1], r1[2]        
            assert x0<x1 and y0<y1 and z0<z1
            def f(r, t=0):
                x, y, z = r
                isin = (x>=x0) & (x<=x1) & (y>=y0) & (y<=y1) & (z>=z0) & (z<=z1)
                isout = np.logical_not(isin)
                res = None
                if inspect.isfunction(a):
                    a_args = inspect.getargspec(a)[0]
                    if len(a_args)==2:
                        res = a(r, t)*isin
                    else:
                        assert len(a_args)==1
                        res = a(r)*isin
                else:
                    res = a*isin
                if inspect.isfunction(b):
                    b_args = inspect.getargspec(b)[0]
                    if len(b_args)==2:
                        res += b(r, t)*isout
                    else:
                        assert len(b_args)==1
                        res += b(r)*isout
                else:
                    res += b*isout
                return res
            return f
        elif self.n_dim==2:
            x0, y0 = r0[0], r0[1]        
            x1, y1 = r1[0], r1[1]        
            assert x0<x1 and y0<y1
            def f(r, t=0):
                x, y = r
                isin = (x>=x0) & (x<=x1) & (y>=y0) & (y<=y1)
                isout = np.logical_not(isin)
                res = None
                if inspect.isfunction(a):
                    a_args = inspect.getargspec(a)[0]
                    if len(a_args)==2:
                        res = a(r, t)*isin
                    else:
                        assert len(a_args)==1
                        res = a(r)*isin
                else:
                    res = a*isin
                if inspect.isfunction(b):
                    b_args = inspect.getargspec(b)[0]
                    if len(b_args)==2:
                        res += b(r, t)*isout
                    else:
                        assert len(b_args)==1
                        res += b(r)*isout
                else:
                    res += b*isout
                return res
            return f
        elif self.n_dim==1:
            x0 = r0[0]        
            x1 = r1[0]        
            assert x0<x1
            def f(r, t=0):
                x = r
                isin = (x>=x0) & (x<=x1)
                isout = np.logical_not(isin)
                res = None
                if inspect.isfunction(a):
                    a_args = inspect.getargspec(a)[0]
                    if len(a_args)==2:
                        res = a(r, t)*isin
                    else:
                        assert len(a_args)==1
                        res = a(r)*isin
                else:
                    res = a*isin
                if inspect.isfunction(b):
                    b_args = inspect.getargspec(b)[0]
                    if len(b_args)==2:
                        res += b(r, t)*isout
                    else:
                        assert len(b_args)==1
                        res += b(r)*isout
                else:
                    res += b*isout
                return res
            return f
        
            

    def GetPointSourceFunc(self, r_0, mag=1.0, src_dir='x', em_type='e'):
        if em_type=='e':
            i_xyz_min, r_xyz_min, d_xyz_min = self.FindClosestSides(r_0)
            f_je = None
            if src_dir=='x':
                i_j, r_j, d_j = i_xyz_min[0], r_xyz_min[0], d_xyz_min[0]

                f_jex = self.FunctionSphere(r_j, d_j/2.0, mag)
                f_jey = lambda r : 0.0
                f_jez = lambda r : 0.0
                f_je = [f_jex, f_jey, f_jez]
            elif src_dir=='y':
                i_j, r_j, d_j = i_xyz_min[1], r_xyz_min[1], d_xyz_min[1]

                f_jex = lambda r : 0.0
                f_jey = self.FunctionSphere(r_j, d_j/2.0, mag)
                f_jez = lambda r : 0.0
                f_je = [f_jex, f_jey, f_jez]
            elif src_dir=='z':
                i_j, r_j, d_j = i_xyz_min[2], r_xyz_min[2], d_xyz_min[2]

                f_jex = lambda r : 0.0
                f_jey = lambda r : 0.0
                f_jez = self.FunctionSphere(r_j, d_j/2.0, mag)
                f_je = [f_jex, f_jey, f_jez]
            else:
                raise ValueError()
            
            return f_je
        elif em_type=='m' or em_type=='h':
            i_xyz_min, r_xyz_min, d_xyz_min = self.FindClosestFaces(r_0)
            f_jm = None
            if src_dir=='x':
                i_j, r_j, d_j = i_xyz_min[0], r_xyz_min[0], d_xyz_min[0]

                f_jmx = self.FunctionSphere(r_j, d_j/2.0, mag)
                f_jmy = lambda r : 0.0
                f_jmz = lambda r : 0.0
                f_jm = [f_jmx, f_jmy, f_jmz]
            elif src_dir=='y':
                i_j, r_j, d_j = i_xyz_min[1], r_xyz_min[1], d_xyz_min[1]

                f_jmx = lambda r : 0.0
                f_jmy = self.FunctionSphere(r_j, d_j/2.0, mag)
                f_jmz = lambda r : 0.0
                f_jm = [f_jmx, f_jmy, f_jmz]
            elif src_dir=='z':
                i_j, r_j, d_j = i_xyz_min[2], r_xyz_min[2], d_xyz_min[2]

                f_jmx = lambda r : 0.0
                f_jmy = lambda r : 0.0
                f_jmz = self.FunctionSphere(r_j, d_j/2.0, mag)
                f_jm = [f_jmx, f_jmy, f_jmz]
            else:
                raise ValueError()

            return f_jm
        else:
            raise ValueError()


    def GetSheetSourceFunc(self, r_0, mag=1.0, src_dir='x', norm_dir='x', em_type='e'):
        assert src_dir in ['x', 'y', 'z']
        assert norm_dir in ['x', 'y', 'z']
        assert em_type in ['e', 'm', 'h']
        if em_type=='e':
            i_xyz_min, r_xyz_min, d_xyz_min = self.FindClosestSides(r_0)
            f_je = None
            if src_dir=='x':
                i_j, r_j, d_j = i_xyz_min[0], r_xyz_min[0], d_xyz_min[0]
                d_min = min(self.dr)

                r0 = np.copy(self.r0)
                r1 = np.copy(self.r1)
                if norm_dir=='x':
                    r0[0] = r_j[0]-d_min/3
                    r1[0] = r_j[0]+d_min/3
                elif norm_dir=='y':
                    assert self.n_dim>=2
                    r0[1] = r_j[1]-d_min/3
                    r1[1] = r_j[1]+d_min/3
                elif norm_dir=='z':
                    assert self.n_dim==3
                    r0[2] = r_j[2]-d_min/3
                    r1[2] = r_j[2]+d_min/3
                
                f_jex = self.FunctionBox(r0, r1, mag)
                f_jey = lambda r : 0.0
                f_jez = lambda r : 0.0
                f_je = [f_jex, f_jey, f_jez]
            elif src_dir=='y':
                i_j, r_j, d_j = i_xyz_min[1], r_xyz_min[1], d_xyz_min[1]
                d_min = min(self.dr)

                r0 = np.copy(self.r0)
                r1 = np.copy(self.r1)
                if norm_dir=='x':
                    r0[0] = r_j[0]-d_min/3
                    r1[0] = r_j[0]+d_min/3
                elif norm_dir=='y':
                    assert self.n_dim>=2
                    r0[1] = r_j[1]-d_min/3
                    r1[1] = r_j[1]+d_min/3
                elif norm_dir=='z':
                    assert self.n_dim==3
                    r0[2] = r_j[2]-d_min/3
                    r1[2] = r_j[2]+d_min/3
                

                f_jex = lambda r : 0.0
                f_jey = self.FunctionBox(r0, r1, mag)
                f_jez = lambda r : 0.0
                f_je = [f_jex, f_jey, f_jez]
            elif src_dir=='z':
                i_j, r_j, d_j = i_xyz_min[2], r_xyz_min[2], d_xyz_min[2]
                d_min = min(self.dr)

                r0 = np.copy(self.r0)
                r1 = np.copy(self.r1)
                if norm_dir=='x':
                    r0[0] = r_j[0]-d_min/3
                    r1[0] = r_j[0]+d_min/3
                elif norm_dir=='y':
                    assert self.n_dim>=2
                    r0[1] = r_j[1]-d_min/3
                    r1[1] = r_j[1]+d_min/3
                elif norm_dir=='z':
                    assert self.n_dim==3
                    r0[2] = r_j[2]-d_min/3
                    r1[2] = r_j[2]+d_min/3

                f_jex = lambda r : 0.0
                f_jey = lambda r : 0.0
                f_jez = self.FunctionBox(r0, r1, mag)
                f_je = [f_jex, f_jey, f_jez]
            else:
                raise ValueError()
            
            return f_je
        elif em_type=='m' or em_type=='h':
            i_xyz_min, r_xyz_min, d_xyz_min = self.FindClosestFaces(r_0)
            f_jm = None
            if src_dir=='x':
                i_j, r_j, d_j = i_xyz_min[0], r_xyz_min[0], d_xyz_min[0]
                d_min = min(self.dr)

                r0 = np.copy(self.r0)
                r1 = np.copy(self.r1)
                if norm_dir=='x':
                    r0[0] = r_j[0]-d_min/3
                    r1[0] = r_j[0]+d_min/3
                elif norm_dir=='y':
                    assert self.n_dim>=2
                    r0[1] = r_j[1]-d_min/3
                    r1[1] = r_j[1]+d_min/3
                elif norm_dir=='z':
                    assert self.n_dim==3
                    r0[2] = r_j[2]-d_min/3
                    r1[2] = r_j[2]+d_min/3

                f_jmx = self.FunctionBox(r0, r1, mag)
                f_jmy = lambda r : 0.0
                f_jmz = lambda r : 0.0
                f_jm = [f_jmx, f_jmy, f_jmz]
            elif src_dir=='y':
                i_j, r_j, d_j = i_xyz_min[1], r_xyz_min[1], d_xyz_min[1]
                d_min = min(self.dr)

                r0 = np.copy(self.r0)
                r1 = np.copy(self.r1)
                if norm_dir=='x':
                    r0[0] = r_j[0]-d_min/3
                    r1[0] = r_j[0]+d_min/3
                elif norm_dir=='y':
                    assert self.n_dim>=2
                    r0[1] = r_j[1]-d_min/3
                    r1[1] = r_j[1]+d_min/3
                elif norm_dir=='z':
                    assert self.n_dim==3
                    r0[2] = r_j[2]-d_min/3
                    r1[2] = r_j[2]+d_min/3
                
                f_jmx = lambda r : 0.0
                f_jmy = self.FunctionBox(r0, r1, mag)
                f_jmz = lambda r : 0.0
                f_jm = [f_jmx, f_jmy, f_jmz]
            elif src_dir=='z':
                i_j, r_j, d_j = i_xyz_min[2], r_xyz_min[2], d_xyz_min[2]
                d_min = min(self.dr)

                r0 = np.copy(self.r0)
                r1 = np.copy(self.r1)
                if norm_dir=='x':
                    r0[0] = r_j[0]-d_min/3
                    r1[0] = r_j[0]+d_min/3
                elif norm_dir=='y':
                    assert self.n_dim>=2
                    r0[1] = r_j[1]-d_min/3
                    r1[1] = r_j[1]+d_min/3
                elif norm_dir=='z':
                    assert self.n_dim==3
                    r0[2] = r_j[2]-d_min/3
                    r1[2] = r_j[2]+d_min/3

                f_jmx = lambda r : 0.0
                f_jmy = lambda r : 0.0
                f_jmz = self.FunctionBox(r0, r1, mag)
                f_jm = [f_jmx, f_jmy, f_jmz]
            else:
                raise ValueError()

            return f_jm
        else:
            raise ValueError()


            

    def FindClosestSides(self, r0):
        """ r0:point
        """
        drx2=None
        if self.n_dim==3:
            drx2 = (self.rsx[0]-r0[0])**2 + (self.rsx[1]-r0[1])**2 + (self.rsx[2]-r0[2])**2
        elif self.n_dim==2:
            drx2 = (self.rsx[0]-r0[0])**2 + (self.rsx[1]-r0[1])**2
        elif self.n_dim==1:
            drx2 = (self.rsx[0]-r0[0])**2
        inds_minx_F = np.argmin(drx2)
        inds_min_x = np.unravel_index(inds_minx_F, self.Nsx)
        x_min = np.array([self.rsx[i][inds_min_x] for i in range(self.n_dim)])
        dx_min = np.sqrt(np.min(drx2))
        print('inds_min_x:', inds_min_x, 'd:', dx_min)

        dry2 = None
        if self.n_dim==3:
            dry2 = (self.rsy[0]-r0[0])**2 + (self.rsy[1]-r0[1])**2 + (self.rsy[2]-r0[2])**2
        elif self.n_dim==2:
            dry2 = (self.rsy[0]-r0[0])**2 + (self.rsy[1]-r0[1])**2
        elif self.n_dim==1:
            dry2 = (self.rsy[0]-r0[0])**2
        inds_miny_F = np.argmin(dry2)
        inds_min_y = np.unravel_index(inds_miny_F, self.Nsy)
        y_min = np.array([self.rsy[i][inds_min_y] for  i in range(self.n_dim)])
        dy_min = np.sqrt(np.min(dry2))
        print('inds_min_y:', inds_min_y, 'd:', np.sqrt(np.min(dry2)))

        drz2 = None
        if self.n_dim==3:
            drz2 = (self.rsz[0]-r0[0])**2 + (self.rsz[1]-r0[1])**2 + (self.rsz[2]-r0[2])**2
        elif self.n_dim==2:
            drz2 = (self.rsz[0]-r0[0])**2 + (self.rsz[1]-r0[1])**2 
        elif self.n_dim==1:
            drz2 = (self.rsz[0]-r0[0])**2 
        inds_minz_F = np.argmin(drz2)
        inds_min_z = np.unravel_index(inds_minz_F, self.Nsz)
        z_min = np.array([self.rsz[i][inds_min_z] for i in range(self.n_dim)])
        dz_min = np.sqrt(np.min(drz2))
        print('inds_min_z:', inds_min_z, 'd:', np.sqrt(np.min(drz2)))
        
        return [inds_min_x, inds_min_y, inds_min_z], [x_min, y_min, z_min],\
            np.array([dx_min, dy_min, dz_min])


    def FindClosestFaces(self, r0):
        """ r0:point
        """
        drx2=None
        if self.n_dim==3:
            drx2 = (self.rfx[0]-r0[0])**2 + (self.rfx[1]-r0[1])**2 + (self.rfx[2]-r0[2])**2
        elif self.n_dim==2:
            drx2 = (self.rfx[0]-r0[0])**2 + (self.rfx[1]-r0[1])**2
        elif self.n_dim==1:
            drx2 = (self.rfx[0]-r0[0])**2
        inds_minx_F = np.argmin(drx2)
        inds_min_x = np.unravel_index(inds_minx_F, self.Nfx)
        x_min = np.array([self.rfx[i][inds_min_x] for i in range(self.n_dim)])
        dx_min = np.sqrt(np.min(drx2))
        print('inds_min_x:', inds_min_x, 'd:', dx_min)

        dry2 = None
        if self.n_dim==3:
            dry2 = (self.rfy[0]-r0[0])**2 + (self.rfy[1]-r0[1])**2 + (self.rfy[2]-r0[2])**2
        elif self.n_dim==2:
            dry2 = (self.rfy[0]-r0[0])**2 + (self.rfy[1]-r0[1])**2
        elif self.n_dim==1:
            dry2 = (self.rfy[0]-r0[0])**2
        inds_miny_F = np.argmin(dry2)
        inds_min_y = np.unravel_index(inds_miny_F, self.Nfy)
        y_min = np.array([self.rfy[i][inds_min_y] for  i in range(self.n_dim)])
        dy_min = np.sqrt(np.min(dry2))
        print('inds_min_y:', inds_min_y, 'd:', np.sqrt(np.min(dry2)))

        drz2 = None
        if self.n_dim==3:
            drz2 = (self.rfz[0]-r0[0])**2 + (self.rfz[1]-r0[1])**2 + (self.rfz[2]-r0[2])**2
        elif self.n_dim==2:
            drz2 = (self.rfz[0]-r0[0])**2 + (self.rfz[1]-r0[1])**2 
        elif self.n_dim==1:
            drz2 = (self.rfz[0]-r0[0])**2 
        inds_minz_F = np.argmin(drz2)
        inds_min_z = np.unravel_index(inds_minz_F, self.Nfz)
        z_min = np.array([self.rfz[i][inds_min_z] for i in range(self.n_dim)])
        dz_min = np.sqrt(np.min(drz2))
        print('inds_min_z:', inds_min_z, 'd:', np.sqrt(np.min(drz2)))
        
        return [inds_min_x, inds_min_y, inds_min_z], [x_min, y_min, z_min],\
            np.array([dx_min, dy_min, dz_min])


    def SetViewPlane_Side(self, r0, A_dir_list):
        """ 3 normal planes passing through r0 will be returned 
            A: the fields to output for example [{'A':self.Ex, 'A_dir:''x', 
            'O_dir':'y', 'v_out':v), ..]
            var_out: save output to this variable (a list) when demanded
        """
        inds = self.FindClosestSides(r0)[0]
        ind_x, ind_y, ind_z = inds
        
        Ax_r = [None]*self.n_dim
        ix_p = [None]*self.n_dim
        for n in range(self.n_dim):
            ix_p[n] = [Ellipsis]*self.n_dim    ##plane normal to n
            ix_p[n][n] = ind_x[n]
            Ax_r[n] = [self.rsx[i][ix_p[n]] for i in range(self.n_dim)]
        
        ix_all = [Ellipsis]*self.n_dim   
        Ax_r_all = self.rsx
        if self.n_dim==0:
            Ax_r_all = self.rsx[0]


        Ay_r = [None]*self.n_dim
        iy_p = [None]*self.n_dim
        for n in range(self.n_dim):
            iy_p[n] = [Ellipsis]*self.n_dim    ##plane normal to n
            iy_p[n][n] = ind_y[n]
            Ay_r[n] = [self.rsy[i][iy_p[n]] for i in range(self.n_dim)]
        
        iy_all = [Ellipsis]*self.n_dim   
        Ay_r_all = self.rsy
        if self.n_dim==0:
            Ay_r_all = self.rsy[0]


        Az_r = [None]*self.n_dim
        iz_p = [None]*self.n_dim
        for n in range(self.n_dim):
            iz_p[n] = [Ellipsis]*self.n_dim    ##plane normal to n
            iz_p[n][n] = ind_z[n]
            Az_r[n] = [self.rsz[i][iz_p[n]] for i in range(self.n_dim)]
        
        iz_all = [Ellipsis]*self.n_dim   
        Az_r_all = self.rsz
        if self.n_dim==0:
            Az_r_all = self.rsz[0]


        for view_dic in A_dir_list:
            A = view_dic['A']           ##field
            A_dir = view_dic['A_dir']   ##field direction
            O_dir = view_dic['O_dir']   ##output plane is normal to this direction
            v_out = view_dic['v_out']   ##output list
            
            if A_dir=='x':
                n_dir=0
                if O_dir=='x':
                    self.ViewPlanes.append( {'v_out':v_out, 'A':A, 'A_dir':n_dir, 'ind':ix_p[0], 'r':Ax_r[0]})
                elif O_dir=='y':
                    assert self.n_dim>=2
                    self.ViewPlanes.append( {'v_out':v_out, 'A':A, 'A_dir':n_dir, 'ind':ix_p[1], 'r':Ax_r[1]})
                elif O_dir=='z':
                    assert self.n_dim==3
                    self.ViewPlanes.append( {'v_out':v_out, 'A':A, 'A_dir':n_dir, 'ind':ix_p[2], 'r':Ax_r[2]})
                elif O_dir==None:
                    self.ViewPlanes.append( {'v_out':v_out, 'A':A, 'A_dir':n_dir, 'ind':ix_all, 'r':Ax_r_all})
                else:
                    raise ValueError()
            elif A_dir=='y':
                n_dir=1
                if O_dir=='x':
                    self.ViewPlanes.append( {'v_out':v_out, 'A':A, 'A_dir':n_dir, 'ind':iy_p[0], 'r':Ay_r[0]})
                elif O_dir=='y':
                    assert self.n_dim>=2
                    self.ViewPlanes.append( {'v_out':v_out, 'A':A, 'A_dir':n_dir, 'ind':iy_p[1], 'r':Ay_r[1]})
                elif O_dir=='z':
                    assert self.n_dim==3
                    self.ViewPlanes.append( {'v_out':v_out, 'A':A, 'A_dir':n_dir, 'ind':iy_p[2], 'r':Ay_r[2]})
                elif O_dir==None:
                    self.ViewPlanes.append( {'v_out':v_out, 'A':A, 'A_dir':n_dir, 'ind':iy_all, 'r':Ay_r_all})
                else:
                    raise ValueError()
            elif A_dir=='z':
                n_dir=2
                if O_dir=='x':
                    self.ViewPlanes.append( {'v_out':v_out, 'A':A, 'A_dir':n_dir, 'ind':iz_p[0], 'r':Az_r[0]})
                elif O_dir=='y':
                    assert self.n_dim>=2
                    self.ViewPlanes.append( {'v_out':v_out, 'A':A, 'A_dir':n_dir, 'ind':iz_p[1], 'r':Az_r[1]})
                elif O_dir=='z':
                    assert self.n_dim==3
                    self.ViewPlanes.append( {'v_out':v_out, 'A':A, 'A_dir':n_dir, 'ind':iz_p[2], 'r':Az_r[2]})
                elif O_dir==None:
                    self.ViewPlanes.append( {'v_out':v_out, 'A':A, 'A_dir':n_dir, 'ind':iz_all, 'r':Az_r_all})
                else:
                    raise ValueError()
            else:
                raise NotImplementedError()
                    
    
    def SaveOutputs(self):
        for vp in self.ViewPlanes:
            v_out = vp['v_out']
            A = vp['A']
            n_dir = vp['A_dir']
            ind = vp['ind']
            r = vp['r']
            
            v_out.append(np.copy(A[n_dir][ind]))
        
    def GetOutputs(self):
        outs = []
        for vp in self.ViewPlanes:
            v_out = vp['v_out']
            A = vp['A']
            n_dir = vp['A_dir']
            ind = vp['ind']
            r = vp['r']
            
            outs.append([r, v_out])
        return outs
        
        


