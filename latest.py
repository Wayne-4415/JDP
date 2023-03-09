import matplotlib.colors as col
import math 
import numpy as np
from math import exp, pi, sqrt, ceil, sin, cos
import time
from config import *  #N, NA, lambda_, z, k, Lx, Nmax, x, fx, x_double, fx_double, ft(f,x,fx), ift(f,x,fx), t_f(x,y), P_f(fx,fy), S_f(fx,fy)

#%%

t=np.zeros((2**N,2**N), dtype = complex)
P=np.zeros((2**N,2**N), dtype = complex)
S=np.zeros((2**N,2**N), dtype = complex)
for i in range(2**N):
    for j in range(2**N):
            t[i][j] = t_f(x[i],x[j])
            P[i][j] = P_f(fx[i],fx[j])
            S[i][j] = S_f(fx[i],fx[j])
plot(np.absolute(t),"t")
plot(np.absolute(P),"P")
plot(np.absolute(S),"S")

#%%
#Huygens–Fresnel
f=1E-6
A=np.zeros((2**N,2**N), dtype = complex)
B=np.zeros((2**N,2**N), dtype = complex)
C=np.zeros((2**N,2**N), dtype = complex)
D=np.zeros((2**N,2**N), dtype = complex)

for i in range(2**N):
    for j in range(2**N):
            A[i][j] = t_f(x[i],x[j])
# plot(np.absolute(A),"A")

for i in range(2**N):
    for j in range(2**N):
        for ii in range(2**N):
            for jj in range(2**N):  
                r = sqrt((x[i]-x[ii])**2 +  (x[j]-x[jj])**2 + f**2)
                B[i][j] += A[ii][jj]*cmath.exp(1j*k*r)/r
# plot(np.absolute(B),"B")

for i in range(2**N):
    for j in range(2**N):
            B[i][j] = B[i][j] * cmath.exp(-1j*k*((x[i])**2 + (x[j])**2)/(2*f))

for i in range(2**N):
    for j in range(2**N):
        for ii in range(2**N):
            for jj in range(2**N):  
                r = sqrt((x[i]-x[ii])**2 + (x[j]-x[jj])**2 + (2*f)**2)
                C[i][j] += B[ii][jj]*cmath.exp(1j*k*r)/r
# plot(np.absolute(C),"C")

for i in range(2**N):
    for j in range(2**N):
            C[i][j] = C[i][j] * cmath.exp(-1j*k*((x[i])**2 + (x[j])**2)/(2*f))

for i in range(2**N):
    for j in range(2**N):
        for ii in range(2**N):
            for jj in range(2**N):  
                r = sqrt((x[i]-x[ii])**2 + (x[j]-x[jj])**2 + f**2)
                D[i][j] += C[ii][jj]*cmath.exp(1j*k*r)/r
plot(np.absolute(D),"D")

#%%
#4f
t = np.zeros((2**N,2**N), dtype=complex)
P = np.zeros((2**N,2**N), dtype=complex)
I_4f = np.zeros((2**N,2**N), dtype=complex)

for i in range(2**N):
    for j in range(2**N):
        t[i][j] = t_f(x[i],x[j])
        P[i][j] = P_f(fx[i],fx[j])

T = ft(t,x,fx)
T = P * T
I_4f = ift(T,x,fx)
plot(np.absolute(I_4f),"I_4f")

#%%
#6f
S = np.zeros((2**N,2**N), dtype=complex)
t = np.zeros((2**N,2**N), dtype=complex)
P = np.zeros((2**N,2**N), dtype=complex)
I_6f = np.zeros((2**N,2**N), dtype=complex)


def P_f(fx,fy): #6F socs case's pupil
    if sqrt(fx**2 + fy**2) < NA/lambda_:
        return 1
    elif (sqrt(fx**2 + fy**2) >= NA/lambda_) and (sqrt(fx**2 + fy**2) <1.6*NA/lambda_) :
        return (-1/(0.6*NA/lambda_)*sqrt(fx**2 + fy**2) + 8/3)
    else:
        return 0

for i in range(2**N):
    for j in range(2**N):
        t[i][j] = t_f(x[i],x[j])
        S[i][j] = S_f(fx[i],fx[j])   
        P[i][j] = P_f(fx[i],fx[j])   
        
T = ft(t,x,fx)
T = P * T
I_6f = ift(T,x,fx) 
plot(np.absolute(I_6f),"I_6f")

