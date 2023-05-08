import numpy as np
from math import pi, sqrt
import cmath
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from os import _exit
import time 
start_time = time.time()

N_x = 120
N_y = 120
N_z = 120
L_x = 40E-6
L_y = 40E-6
L_z = 40E-6
delta_t = 0.05E-15
permittivity = 8.85E-12
permeability = 1.25E-6
pml = np.zeros((N_x, N_y, N_z))

x = np.zeros((N_x))
y = np.zeros((N_y))
z = np.zeros((N_z))
index_x = np.zeros((N_x))
index_y = np.zeros((N_y))
index_z = np.zeros((N_z))

for i in range(N_x):
    x[i] = -L_x/2 + L_x*i/N_x
    index_x[i] = -N_x//2 + i
for i in range(N_y):
    y[i] = -L_y/2 + L_y*i/N_y
    index_y[i] = -N_y//2 + i
for i in range(N_z):
    z[i] = -L_z/2 + L_z*i/N_z
    index_z[i] = -N_z//2 + i
    
xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
index_xx, index_yy, index_zz =  np.meshgrid(index_x, index_y, index_z, indexing='ij')

def plot(figure,title):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect('equal', adjustable='box')
    ax.set_title(title)
    contourf = ax.contourf(np.absolute(figure), 400, vmin=figure.min(), vmax=figure.max(), levels = 50)
    cbar = fig.colorbar(contourf)
    plt.show()

def plot3d(figure,title):
    mpl.rcParams['figure.dpi'] = 160
    mpl.rcParams['figure.figsize'] = [8, 8]  #[x,y] 
    mpl.rcParams['font.size'] = 13
    fig = plt.figure()
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    ax.set_title(title)
    X = np.arange(0, np.shape(figure)[0], 1)
    Y = np.arange(0, np.shape(figure)[1], 1)
    X, Y = np.meshgrid(X, Y)
    Z = figure
    surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=True)
    ax.set_zlim(-0.5, 1)
    plt.show()

def PML(r):
    pml_thickness = 4E-6
    Q = 1 / delta_t
    if r < -L_x/2 + pml_thickness:
        return Q * ((-L_x/2 + pml_thickness - r)/pml_thickness)**3
    elif r >= L_x/2 - pml_thickness:
        return Q * ((-L_x/2 + pml_thickness + r)/pml_thickness)**3
    else:
        return 0

for i in range(N_x):
    for j in range(N_y):
        for k in range(N_z):
            pml[i,j,k] = PML(x[i]) + PML(y[j]) + PML(z[k])
# DD=(1 - pml * delta_t / 2)/(1 + pml * delta_t / 2)
# plot(np.real(DD[:,:,N_z//2]),"vdsg")
# os.stop


Hx = np.zeros((N_x,N_y,N_z), dtype = complex)
Hy = np.zeros((N_x,N_y,N_z), dtype = complex)
Hz = np.zeros((N_x,N_y,N_z), dtype = complex)
Ex = np.zeros((N_x,N_y,N_z), dtype = complex)
Ey = np.zeros((N_x,N_y,N_z), dtype = complex)
Ez = np.zeros((N_x,N_y,N_z), dtype = complex)


kx = 2*pi/L_x * index_xx[:,:,:]
ky = 2*pi/L_y * index_yy[:,:,:]
kz = 2*pi/L_z * index_zz[:,:,:]

def Ez_f(x,y,z):
    return np.exp(-((np.square(x)+np.square(y))/(2E-6)**2))
    # return np.where(np.sqrt(np.square(x)+np.square(y)+np.square(z))<5E-6,1,0)


Ez = Ez_f(xx,yy,zz)
# plot(np.absolute(Ez[:,:,N//2]),"Ez")

def fold(M):
    folded_M = np.zeros((N_x,N_y,N_z), dtype = complex)
    folded_M[0:N_x//2, 0:N_y//2, 0:N_z//2] = M[N_x//2:N_x, N_y//2:N_y, N_z//2:N_z]
    folded_M[0:N_x//2, 0:N_y//2, N_z//2:N_z] = M[N_x//2:N_x, N_y//2:N_y, 0:N_z//2]
    folded_M[0:N_x//2, N_y//2:N_y, 0:N_z//2] = M[N_x//2:N_x, 0:N_y//2, N_z//2:N_z]
    folded_M[0:N_x//2, N_y//2:N_y, N_z//2:N_z] = M[N_x//2:N_x, 0:N_y//2, 0:N_z//2]
    folded_M[N_x//2:N_x, 0:N_y//2, 0:N_z//2] = M[0:N_x//2, N_y//2:N_y, N_z//2:N_z]
    folded_M[N_x//2:N_x, 0:N_y//2, N_z//2:N_z] = M[0:N_x//2, N_y//2:N_y, 0:N_z//2]
    folded_M[N_x//2:N_x, N_y//2:N_y, 0:N_z//2] = M[0:N_x//2, 0:N_y//2, N_z//2:N_z]
    folded_M[N_x//2:N_x, N_y//2:N_y, N_z//2:N_z] = M[0:N_x//2, 0:N_y//2, 0:N_z//2]
    return folded_M

for t in range(2101):
    fft_y_Ex = np.fft.fft(Ex, axis = 1)
    fft_z_Ex = np.fft.fft(Ex, axis = 2)
    fft_x_Ey = np.fft.fft(Ey, axis = 0)
    fft_z_Ey = np.fft.fft(Ey, axis = 2)
    fft_x_Ez = np.fft.fft(Ez, axis = 0)
    fft_y_Ez = np.fft.fft(Ez, axis = 1)
    
    fft_y_Ex = fold(fft_y_Ex)
    fft_z_Ex = fold(fft_z_Ex)
    fft_x_Ey = fold(fft_x_Ey)
    fft_z_Ey = fold(fft_z_Ey)
    fft_x_Ez = fold(fft_x_Ez)
    fft_y_Ez = fold(fft_y_Ez)
    
    
    ky_fft_y_Ex = 1j * ky * fft_y_Ex
    kz_fft_z_Ex = 1j * kz * fft_z_Ex
    kx_fft_x_Ey = 1j * kx * fft_x_Ey
    kz_fft_z_Ey = 1j * kz * fft_z_Ey
    kx_fft_x_Ez = 1j * kx * fft_x_Ez
    ky_fft_y_Ez = 1j * ky * fft_y_Ez
    
    ky_fft_y_Ex = fold(ky_fft_y_Ex)
    kz_fft_z_Ex = fold(kz_fft_z_Ex)
    kx_fft_x_Ey = fold(kx_fft_x_Ey)
    kz_fft_z_Ey = fold(kz_fft_z_Ey)
    kx_fft_x_Ez = fold(kx_fft_x_Ez)
    ky_fft_y_Ez = fold(ky_fft_y_Ez)
    
    ifft_ky_fft_y_Ex = np.fft.ifft(ky_fft_y_Ex, axis = 1)
    ifft_kz_fft_z_Ex = np.fft.ifft(kz_fft_z_Ex, axis = 2)
    ifft_kx_fft_x_Ey = np.fft.ifft(kx_fft_x_Ey, axis = 0)
    ifft_kz_fft_z_Ey = np.fft.ifft(kz_fft_z_Ey, axis = 2)
    ifft_kx_fft_x_Ez = np.fft.ifft(kx_fft_x_Ez, axis = 0)
    ifft_ky_fft_y_Ez = np.fft.ifft(ky_fft_y_Ez, axis = 1)
    
    Hx = (1 - pml * delta_t / 2)/(1 + pml * delta_t / 2) * Hx - (delta_t/permeability) * (1 / (1 + pml * delta_t / 2)) * (ifft_ky_fft_y_Ez - ifft_kz_fft_z_Ey)
    Hy = (1 - pml * delta_t / 2)/(1 + pml * delta_t / 2) * Hy - (delta_t/permeability) * (1 / (1 + pml * delta_t / 2)) * (ifft_kz_fft_z_Ex - ifft_kx_fft_x_Ez)
    Hz = (1 - pml * delta_t / 2)/(1 + pml * delta_t / 2) * Hz - (delta_t/permeability) * (1 / (1 + pml * delta_t / 2)) * (ifft_kx_fft_x_Ey - ifft_ky_fft_y_Ex)
    
    fft_y_Hx = np.fft.fft(Hx, axis = 1)
    fft_z_Hx = np.fft.fft(Hx, axis = 2)
    fft_x_Hy = np.fft.fft(Hy, axis = 0)
    fft_z_Hy = np.fft.fft(Hy, axis = 2)
    fft_x_Hz = np.fft.fft(Hz, axis = 0)
    fft_y_Hz = np.fft.fft(Hz, axis = 1)
    
    fft_y_Hx = fold(fft_y_Hx)
    fft_z_Hx = fold(fft_z_Hx)
    fft_x_Hy = fold(fft_x_Hy)
    fft_z_Hy = fold(fft_z_Hy)
    fft_x_Hz = fold(fft_x_Hz)
    fft_y_Hz = fold(fft_y_Hz)
    
    ky_fft_y_Hx = 1j * ky * fft_y_Hx
    kz_fft_z_Hx = 1j * kz * fft_z_Hx
    kx_fft_x_Hy = 1j * kx * fft_x_Hy
    kz_fft_z_Hy = 1j * kz * fft_z_Hy
    kx_fft_x_Hz = 1j * kx * fft_x_Hz
    ky_fft_y_Hz = 1j * ky * fft_y_Hz
    
    ky_fft_y_Hx = fold(ky_fft_y_Hx)
    kz_fft_z_Hx = fold(kz_fft_z_Hx)
    kx_fft_x_Hy = fold(kx_fft_x_Hy)
    kz_fft_z_Hy = fold(kz_fft_z_Hy)
    kx_fft_x_Hz = fold(kx_fft_x_Hz)
    ky_fft_y_Hz = fold(ky_fft_y_Hz)
    
    ifft_ky_fft_y_Hx = np.fft.ifft(ky_fft_y_Hx, axis = 1)
    ifft_kz_fft_z_Hx = np.fft.ifft(kz_fft_z_Hx, axis = 2)
    ifft_kx_fft_x_Hy = np.fft.ifft(kx_fft_x_Hy, axis = 0)
    ifft_kz_fft_z_Hy = np.fft.ifft(kz_fft_z_Hy, axis = 2)
    ifft_kx_fft_x_Hz = np.fft.ifft(kx_fft_x_Hz, axis = 0)
    ifft_ky_fft_y_Hz = np.fft.ifft(ky_fft_y_Hz, axis = 1)
    
    Ex = (1 - pml * delta_t / 2)/(1 + pml * delta_t / 2) * Ex + (delta_t/permittivity) * (1 / (1 + pml * delta_t / 2)) * (ifft_ky_fft_y_Hz - ifft_kz_fft_z_Hy)
    Ey = (1 - pml * delta_t / 2)/(1 + pml * delta_t / 2) * Ey + (delta_t/permittivity) * (1 / (1 + pml * delta_t / 2)) * (ifft_kz_fft_z_Hx - ifft_kx_fft_x_Hz)
    Ez = (1 - pml * delta_t / 2)/(1 + pml * delta_t / 2) * Ez + (delta_t/permittivity) * (1 / (1 + pml * delta_t / 2)) * (ifft_kx_fft_x_Hy - ifft_ky_fft_y_Hx)
    
    print(t)
    if t%300 == 0:
        plot3d(np.real(Ez[:,:,N_z//2]),t)
# plot3d(np.absolute(Ez),"Ez")
# plot(np.absolute(Ez[:,:,N//2]),"Ez")
print("process time: %s seconds ---" % (time.time() - start_time))
    