import numpy as np
import matplotlib.pyplot as plt
import math 
import cmath
from math import exp, pi, sqrt, ceil, sin, cos
import config
N, NA,lambda_, z, k, Lx, Nmax, x, fx, x_double, fx_double = config.parameter()

def plot(figure,title):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect('equal', adjustable='box')
    ax.set_title(title)
    contourf = ax.contourf(np.real(figure), 400, vmin=figure.min(), vmax=figure.max(), levels = 50)
    cbar = fig.colorbar(contourf)
    plt.show()
    
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
TCC=(np.loadtxt("a.txt"))


def svd():
    phi = np.zeros((2**N,2**N), dtype=complex)
    I_total = np.zeros((2**N,2**N))
    t = np.zeros((2**N,2**N), dtype=complex)
    for i in range(2**N):
        for j in range(2**N):
            t[i][j] = t_f(x[i],x[j])
    plot(np.absolute(t),"t")
    T = ft(t,x,fx)
    for state in range(Nmax):
        for i in range(2**N):
            for j in range(2**N):
                phi[i][j] = kernal[state][i][j]
        plot(np.absolute(phi),"kernal")
        phi = phi * T
        phi = ift(phi,x,fx)
        I_total += np.real(np.square(np.absolute(phi)) * kernal_val[state])
    return I_total

def t_f(x,y):
    # if (abs(abs(x) - 0.25E-6) < 0.125E-6) & (abs(abs(y) - 0.25E-6) < 0.125E-6) :
    if (abs(abs(x) - 0.0E-6) < 0.4E-6) & (abs(abs(y) - 0.24E-6) < 0.08E-6) :
        return 1
    else:
        return 0
def kernal():
    u, s, vh = np.linalg.svd(TCC)
    kernal = np.zeros((Nmax,2**N,2**N), dtype=complex)
    kernal_val = np.zeros((Nmax), dtype=complex)
    for state in range(Nmax):
        for i in range(2**N):
            for j in range(2**N):
                kernal[state][i][j] = vh[state][i*2**N + j]
                kernal_val[state] = s[state]
    return kernal ,kernal_val

kernal ,kernal_val = kernal()
plot(np.absolute(svd()),"I_total")