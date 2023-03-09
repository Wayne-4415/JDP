import numpy as np
from math import exp, pi, sqrt, ceil, sin, cos
import cmath
import matplotlib.pyplot as plt

N = 5
NA = 0.6
lambda_ = 248E-9
z = 1E-2
k = 2*pi/lambda_
Lx = 2E-6
Nmax = 10 #use in Kernal function, must be less than N square
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

def ft(f,x,fx):
    temp = np.zeros((len(fx), len(fx)), dtype=complex)
    for i in range(len(fx)):
        for j in range(len(fx)):
            for ii in range(len(x)):
                for jj in range(len(x)):
                    temp[i,j] += f[ii,jj]*cmath.exp(-1j*2*pi*(x[ii]*fx[i]+x[jj]*fx[j])) *(Lx/len(x))
    return temp
def ift(f,x,fx):
    temp = np.zeros((len(x), len(x)), dtype=complex)
    for i in range(len(x)):
        for j in range(len(x)):
            for ii in range(len(fx)):
                for jj in range(len(fx)):
                    temp[i,j] += f[ii,jj]*cmath.exp(1j*2*pi*(x[ii]*fx[i]+x[jj]*fx[j])) *(0.5/Lx)
    return temp
def plot(figure,title):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect('equal', adjustable='box')
    ax.set_title(title)
    contourf = ax.contourf(np.real(figure), 400, vmin=figure.min(), vmax=figure.max(), levels = 50)
    cbar = fig.colorbar(contourf)
    plt.show()

def t_f(x,y):
    if sqrt(x**2 + y**2) < 0.6E-6:
        return 1
    else:
        return 0

def P_f(fx,fy): #1 Conventional illumination
    if sqrt(fx**2 + fy**2) < NA/lambda_:
        return 1
    else:
        return 0

def S_f(fx,fy): #1 Conventional illumination
    if sqrt(fx**2 + fy**2) < 0.6 * NA/lambda_:
        return 1
    else:
        return 0