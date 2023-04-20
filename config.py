import numpy as np
from math import pi, sqrt
import cmath
import matplotlib.pyplot as plt

N = 6
NA = 0.6
lambda_ = 248E-9
z = 1E-2
k = 2*pi/lambda_
Lx = 2E-6
Nmax = 10 #use in Kernal function, must be less than 2**(2*N)
x = np.zeros((2**N))
fx = np.zeros((2**N))
x_double = np.zeros((2**(N+1)))
fx_double = np.zeros((2**(N+1)))

for i in range(2**N):
    x[i] = -Lx/2 + Lx*i/2**N + Lx/2**(N+1)
    fx[i] = (i-(2**N-1)/2)/Lx
for i in range(2**(N+1)):
    x_double[i] = -1*Lx + 2*Lx*i/2**(N+1) + Lx/2**(N+1)
    fx_double[i] = (i-(2**(N+1)-1)/2)/Lx #same resolution as fx but double range
    
xx, yy = np.meshgrid(x, x)
fxx, fyy = np.meshgrid(fx, fx)
xx_double, yy_double = np.meshgrid(x_double, x_double)
fxx_double, fyy_double = np.meshgrid(fx_double, fx_double)

def ft(f,x,fx):
    temp = np.zeros((len(fx), len(fx)), dtype=complex)
    for i in range(len(fx)):
        for j in range(len(fx)):
            for ii in range(len(x)):
                for jj in range(len(x)):
                    temp[i,j] += f[ii,jj]*cmath.exp(-1j*2*pi*(x[ii]*fx[i]+x[jj]*fx[j])) * (x[1]-x[0])**2
    return temp * (1/(2*pi))

def ift(f,x,fx):
    temp = np.zeros((len(x), len(x)), dtype=complex)
    for i in range(len(x)):
        for j in range(len(x)):
            for ii in range(len(fx)):
                for jj in range(len(fx)):
                    temp[i,j] += f[ii,jj]*cmath.exp(1j*2*pi*(x[ii]*fx[i]+x[jj]*fx[j])) * (fx[1]-fx[0])**2
    return temp * (2*pi)

def plot(figure,title):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect('equal', adjustable='box')
    ax.set_title(title)
    contourf = ax.contourf(np.real(figure), 400, vmin=figure.min(), vmax=figure.max(), levels = 50)
    cbar = fig.colorbar(contourf)
    plt.show()

def enclosure(x,y,x1,y1,x2,y2):
    if ((x>x1) and (x<x2)) and ((y>y1) and (y<y2)):
        return 1
    else:
        return 0
    
def t_f(x,y):
    # return np.where(np.sqrt(np.square(x) + np.square(y)) < 0.6E-6, 1, 0)
    mask = np.zeros((len(x),len(y)))
    for i in range(len(x)):
        for j in range(len(x)):
            positionX = x[0,i]
            positionY = y[j,0]
            mask[i,j] += enclosure(positionX,positionY,0.5E-6, -0.7E-6, 0.7E-6, 0.7E-6)
            mask[i,j] += enclosure(positionX,positionY,0.1E-6, -0.3E-6, 0.3E-6, 0.7E-6)
            mask[i,j] += enclosure(positionX,positionY,-0.7E-6, -0.7E-6, 0.7E-6, -0.5E-6)
            mask[i,j] += enclosure(positionX,positionY,-0.7E-6, -0.3E-6, 0.3E-6, -0.1E-6)
            mask[i,j] += enclosure(positionX,positionY,-0.3E-6, 0.1E-6, -0.1E-6, 0.7E-6)
            mask[i,j] += enclosure(positionX,positionY,-0.7E-6, 0.1E-6, -0.1E-6, 0.3E-6)
    return np.where(mask !=0, 1, 0)

def P_f(fx,fy):
    return np.where(np.sqrt(np.square(fx) + np.square(fy)) < NA/lambda_, 1, 0)

def S_f(fx,fy):
    return np.where(np.sqrt(np.square(fx) + np.square(fy)) < 0.2 * NA/lambda_, 1, 0)