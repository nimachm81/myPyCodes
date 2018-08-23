## plane wave expansion method
## given the transverse fields on a plane finds the fields everywhere on the right
## or left half-space


__all__ = ["FourierReciprocalSpace", "Fourier1D", "InvFourier1D", "Fourier2D", "InvFourier2D",
           "PWE2D", "PWE3D"]



import numpy as np
from scipy import constants



""" Fourier transforms:
    The input direct space range is arbitrary, the output reciprocal space range is
    chosen symmetrically 
"""

def FourierReciprocalSpace(x):
    N = len(x)
    dx = x[1]-x[0]
    _Dx = 1/dx
    _x = np.linspace(-_Dx/2, _Dx/2, N, endpoint=False)
    return _x

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





class PWE2D:
    """ 2D problem in the yz plane
        uniform along x
        fields are specified at z=0
        the wave propagates towards +z or -z
        exp( -j*wt + j*k*r) convention
    """
    def __init__(self, f, y, Ex, Ey, epsilon_r=1, mu_r=1):
        self.mu = constants.mu_0*mu_r
        self.epsilon = constants.epsilon_0*epsilon_r
        self.w = 2.0*np.pi*f
        self.k0 = self.w*np.sqrt(self.mu*self.epsilon)
        self.y = y
        self.Ex = Ex
        self.Ey = Ey
        return
        
    
    def setTransverseField(self, y, Ex, Ey):
        self.y = y
        self.Ex = Ex
        self.Ey = Ey
        
        
    def setFrequency(self, f):
        self.w = 2.0*np.pi*f
        self.k0 = self.w*np.sqrt(self.mu*self.epsilon)


    def GetK(self):
        _y, _Ex = Fourier1D(self.y, self.Ex)
        _y, _Ey = Fourier1D(self.y, self.Ey)
        kx = 0.0
        ky = 2.0*np.pi*_y
        self._y = _y
        self.ky = ky
        self.Ex_f = _Ex
        self.Ey_f = _Ey
        
        kz = np.sqrt(self.k0**2 - ky**2 + 0j)
        self.kz = kz

        _Ez = -ky*_Ey/kz
        self.Ez_f = _Ez
        
        ## jw H = jk x E
        self.Hx_f = (ky*self.Ez_f - kz*self.Ey_f)/(self.w*self.mu)
        self.Hy_f = (kz*self.Ex_f - kx*self.Ez_f)/(self.w*self.mu)
        self.Hz_f = (kx*self.Ey_f - ky*self.Ex_f)/(self.w*self.mu)
        
        

    
    def GetFields(self, z):
        """ finds the fields in the [z0, z1] range
        """
        y = self.y
        _y = self._y
        kz = self.kz
        
        Z, Y = np.meshgrid(z, y, indexing='ij')
        
        Nz, Ny = len(z), len(y)
        Ex_yz_f = np.zeros((Nz, Ny), dtype=complex)
        Ey_yz_f = np.zeros((Nz, Ny), dtype=complex)
        Ez_yz_f = np.zeros((Nz, Ny), dtype=complex)
        Hx_yz_f = np.zeros((Nz, Ny), dtype=complex)
        Hy_yz_f = np.zeros((Nz, Ny), dtype=complex)
        Hz_yz_f = np.zeros((Nz, Ny), dtype=complex)
        for i in range(Nz):
            Ex_yz_f[i,:] = self.Ex_f*np.exp(1j*kz*z[i])
            Ey_yz_f[i,:] = self.Ey_f*np.exp(1j*kz*z[i])
            Ez_yz_f[i,:] = self.Ez_f*np.exp(1j*kz*z[i])
            Hx_yz_f[i,:] = self.Hx_f*np.exp(1j*kz*z[i])
            Hy_yz_f[i,:] = self.Hy_f*np.exp(1j*kz*z[i])
            Hz_yz_f[i,:] = self.Hz_f*np.exp(1j*kz*z[i])

        Ex_yz = np.zeros((Nz, Ny), dtype=complex)
        Ey_yz = np.zeros((Nz, Ny), dtype=complex)
        Ez_yz = np.zeros((Nz, Ny), dtype=complex)
        Hx_yz = np.zeros((Nz, Ny), dtype=complex)
        Hy_yz = np.zeros((Nz, Ny), dtype=complex)
        Hz_yz = np.zeros((Nz, Ny), dtype=complex)
        for i in range(Nz):
            Ex_yz[i, :] = InvFourier1D(y, _y, Ex_yz_f[i,:])
            Ey_yz[i, :] = InvFourier1D(y, _y, Ey_yz_f[i,:])
            Ez_yz[i, :] = InvFourier1D(y, _y, Ez_yz_f[i,:])
            Hx_yz[i, :] = InvFourier1D(y, _y, Hx_yz_f[i,:])
            Hy_yz[i, :] = InvFourier1D(y, _y, Hy_yz_f[i,:])
            Hz_yz[i, :] = InvFourier1D(y, _y, Hz_yz_f[i,:])
        
        
        return Z, Y, [Ex_yz, Ey_yz, Ez_yz], [Hx_yz, Hy_yz, Hz_yz]


    def BackPropagateWithGain(self, Ex_z0, Ey_z0, z0):
        """ Ex and Ey are provided at z=z0
            The fields are then back propagated to z=0 plane while the evanescent
            fields are amplified during back propagation
        """
        y = self.y
        _y = self._y
        kz = self.kz

        Ny = len(y)
        
        _y, Ex_z0_f = Fourier1D(y, Ex_z0)
        _y, Ey_z0_f = Fourier1D(y, Ey_z0)
                
        Ex_0_f = Ex_z0_f*np.exp(-1j*kz*z0)
        Ey_0_f = Ey_z0_f*np.exp(-1j*kz*z0)
        
        Ex_0 = InvFourier1D(y, _y, Ex_z0_f)
        Ey_0 = InvFourier1D(y, _y, Ey_z0_f)
    
        return Ex_0, Ey_0



class PWE3D:
    """ 3D problem in the xyz plane
        fields are specified at z=0
        the wave propagates towards +z or -z
        exp( -j*wt + j*k*r) convention
    """
    def __init__(self, f, X, Y, Ex, Ey, epsilon_r=1, mu_r=1):
        """ X, Y should be arranged with ij indexing in the meshgrid
        """
        self.mu = constants.mu_0*mu_r
        self.epsilon = constants.epsilon_0*epsilon_r
        self.w = 2.0*np.pi*f
        self.k0 = self.w*np.sqrt(self.mu*self.epsilon)
        self.x = X[:,0]
        self.y = Y[0,:]
        self._x = FourierReciprocalSpace(self.x)
        self._y = FourierReciprocalSpace(self.y)
        self.X = X
        self.Y = Y
        self.Ex = Ex
        self.Ey = Ey
        return
        
    
    def setTransverseField(self, x, y, Ex, Ey):
        self.x = X[:,0]
        self.y = Y[0,:]
        self._x = FourierReciprocalSpace(self.x)
        self._y = FourierReciprocalSpace(self.y)
        self.X = X
        self.Y = Y
        self.Ex = Ex
        self.Ey = Ey
        
        
    def setFrequency(self, f):
        self.w = 2.0*np.pi*f
        self.k0 = self.w*np.sqrt(self.mu*self.epsilon)


    def GetK(self):
        _X, _Y, _Ex = Fourier2D(self.X, self.Y, self.Ex)
        _X, _Y, _Ey = Fourier2D(self.X, self.Y, self.Ey)
        Kx = 2.0*np.pi*_X
        Ky = 2.0*np.pi*_Y
        self._X = _X
        self._Y = _Y
        self.Kx = Kx
        self.Ky = Ky
        self.Ex_f = _Ex
        self.Ey_f = _Ey
        
        Kz = np.sqrt(self.k0**2 - Kx**2 - Ky**2 + 0j)
        self.Kz = Kz

        _Ez = -(Kx*_Ex+Ky*_Ey)/Kz
        self.Ez_f = _Ez
        
        ## jw H = jk x E
        self.Hx_f = (Ky*self.Ez_f - Kz*self.Ey_f)/(self.w*self.mu)
        self.Hy_f = (Kz*self.Ex_f - Kx*self.Ez_f)/(self.w*self.mu)
        self.Hz_f = (Kx*self.Ey_f - Ky*self.Ex_f)/(self.w*self.mu)
        
        
            
    def GetFields(self, z):
        """ finds the fields in the yz plane for [z0, z1] range
        """
        x = self.x
        y = self.y
        X = self.X
        Y = self.Y
        _X = self._X
        _Y = self._Y
        Kz = self.Kz
        
        X_3D, Y_3D, Z_3D = np.meshgrid(x, y, z, indexing='ij')
        
        Nx, Ny, Nz = len(x), len(y), len(z)
        Ex_xyz_f = np.zeros((Nx, Ny, Nz), dtype=complex)
        Ey_xyz_f = np.zeros((Nx, Ny, Nz), dtype=complex)
        Ez_xyz_f = np.zeros((Nx, Ny, Nz), dtype=complex)
        Hx_xyz_f = np.zeros((Nx, Ny, Nz), dtype=complex)
        Hy_xyz_f = np.zeros((Nx, Ny, Nz), dtype=complex)
        Hz_xyz_f = np.zeros((Nx, Ny, Nz), dtype=complex)
        for i in range(Nz):
            Ex_xyz_f[:,:,i] = self.Ex_f*np.exp(1j*Kz*z[i])
            Ey_xyz_f[:,:,i] = self.Ey_f*np.exp(1j*Kz*z[i])
            Ez_xyz_f[:,:,i] = self.Ez_f*np.exp(1j*Kz*z[i])
            Hx_xyz_f[:,:,i] = self.Hx_f*np.exp(1j*Kz*z[i])
            Hy_xyz_f[:,:,i] = self.Hy_f*np.exp(1j*Kz*z[i])
            Hz_xyz_f[:,:,i] = self.Hz_f*np.exp(1j*Kz*z[i])

        Ex_xyz = np.zeros((Nx, Ny, Nz), dtype=complex)
        Ey_xyz = np.zeros((Nx, Ny, Nz), dtype=complex)
        Ez_xyz = np.zeros((Nx, Ny, Nz), dtype=complex)
        Hx_xyz = np.zeros((Nx, Ny, Nz), dtype=complex)
        Hy_xyz = np.zeros((Nx, Ny, Nz), dtype=complex)
        Hz_xyz = np.zeros((Nx, Ny, Nz), dtype=complex)
        for i in range(Nz):
            Ex_xyz[:, :, i] = InvFourier2D(X, Y, _X, _Y, Ex_xyz_f[:,:,i])
            Ey_xyz[:, :, i] = InvFourier2D(X, Y, _X, _Y, Ey_xyz_f[:,:,i])
            Ez_xyz[:, :, i] = InvFourier2D(X, Y, _X, _Y, Ez_xyz_f[:,:,i])
            Hx_xyz[:, :, i] = InvFourier2D(X, Y, _X, _Y, Hx_xyz_f[:,:,i])
            Hy_xyz[:, :, i] = InvFourier2D(X, Y, _X, _Y, Hy_xyz_f[:,:,i])
            Hz_xyz[:, :, i] = InvFourier2D(X, Y, _X, _Y, Hz_xyz_f[:,:,i])
        
        
        return [X_3D, Y_3D, Z_3D], [Ex_xyz, Ey_xyz, Ez_xyz], [Hx_xyz, Hy_xyz, Hz_xyz]




