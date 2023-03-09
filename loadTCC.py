import numpy as np
import matplotlib.pyplot as plt
import math 
import cmath
from math import exp, pi, sqrt, ceil, sin, cos
from config import *  #N, NA, lambda_, z, k, Lx, Nmax, x, fx, x_double, fx_double, ft(f,x,fx), ift(f,x,fx), t_f(x,y), P_f(fx,fy), S_f(fx,fy)

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
        # plot(np.absolute(phi),"kernal")
        phi = phi * T
        phi = ift(phi,x,fx)
        I_total += np.real(np.square(np.absolute(phi)) * kernal_val[state])
    return I_total

def t_f(x,y):
    if sqrt(x**2 + y**2) < 0.6E-6:
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