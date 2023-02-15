#%%
import matplotlib.colors as col
import math 
import cmath
import numpy as np
from math import exp, pi, sqrt, ceil, sin, cos
import matplotlib as mpl
import matplotlib.pyplot as plt
import inspect
import time

start_time = time.time()
#全域變數EEEEEE
N = 5
NA = 0.6
lambda_ = 248E-9
#sigma = 0.1
z = 1*lambda_
k = 2*pi/lambda_
Lx = 1.5E-6
Nmax = N**2 #use in Kernal function, must be less than N square
print("cell 1 finished")

#%%
x = np.zeros((2**N))
fx = np.zeros((2**N))
x_double = np.zeros((2**(N+1))) #for double the source area  //same as double xii
fx_double = np.zeros((2**(N+1)))
x_quater = np.zeros((2**(N+2)))
fx_quater = np.zeros((2**(N+2)))

for i in range(2**N):
    x[i] = -1*Lx/2 + Lx*i/2**N + Lx/2**(N+1)
    fx[i] = (i-(2**N-1)/2)/Lx

for i in range(2**(N+1)):
    x_double[i] = -1*Lx + 2*Lx*i/2**(N+1) + Lx/2**(N+1)
    fx_double[i] = (i-(2**(N+1)-1)/2)/2/Lx

for i in range(2**(N+2)):
    x_quater[i] = -2*Lx + 4*Lx*i/2**(N+2) + Lx/2**(N+1)
    fx_quater[i] = (i-(2**(N+2)-1)/2)/2/2/Lx
print("cell 2 finished")

#%%
def ft(f,x,fx):
    temp = np.zeros((len(fx), len(fx)), dtype=complex)
    for i in range(len(fx)):
        for j in range(len(fx)):
            for ii in range(len(x)):
                for jj in range(len(x)):
                    temp[i,j] += f[ii,jj]*cmath.exp(-1j*2*pi*(x[ii]*fx[i]+x[jj]*fx[j]))
    return temp

def ift(f,x,fx):
    temp = np.zeros((len(x), len(x)), dtype=complex)
    for i in range(len(x)):
        for j in range(len(x)):
            for ii in range(len(fx)):
                for jj in range(len(fx)):
                    temp[i,j] += f[ii,jj]*cmath.exp(1j*2*pi*(x[ii]*fx[i]+x[jj]*fx[j]))
    return temp
print("cell 3 finished")

#%%
def plot(figure,title):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect('equal', adjustable='box')
    ax.set_title(title)
    contourf = ax.contourf(np.real(figure), 400, vmin=figure.min(), vmax=figure.max(), levels = 50)
    cbar = fig.colorbar(contourf)
    plt.show()
print("cell 4 finished")

#%%
#pupil and source function
#pupil function P(fx,fy) in f space
def P_f(fx,fy):
    if sqrt(fx**2 + fy**2) < NA/lambda_:
        #temp1 = 1j*2*pi*z*lambda_*(fx**2 + fy**2)/2
        #result = cm.exp(temp1)
        #return result
        return 1
    else:
        return 0

#source in f space

def S_f(fx,fy): #1 Conventional illumination
    if sqrt(fx**2 + fy**2) < 0.6*NA/lambda_:
        return 1
    else:
        return 0

# def S_f(fx,fy): #2 Quadrupole illumination
#     nor = Lx**2*NA/(lambda_**2*z*2**N)
#     if sqrt((abs(fx) - 0.5*nor)**2 + (abs(fy) - 0.5*nor)**2) < 0.1*nor:
#         return 1
#     else:
#         return 0

# def S_f(fx,fy): #3 Annular illumination
#     nor = Lx**2*NA/(lambda_**2*z*2**N)
#     if abs(sqrt(fx**2 + fy**2) - 0.65*nor) < 0.05*nor:
#         return 1
#     else:
#         return 0

# def S_f(fx,fy): #4 Subresolution serifs at the contact corners
#     nor = xx / lambda_ / z
#     if ((abs(abs(fx) - 12.5E-9 * nor) < 125E-9 * nor) & (abs(abs(fy) - 12.5E-9 * nor) < 125E-9 * nor)) or ((abs(fx) < 125E-9 * nor) and (abs(fy) < 125E-9 * nor)):
#         return 1
#     else:
#         return 0

# def S_f(fx,fy): #5 Alternate PSM
#     if :
#         return 1
#     else:
#         return 0

# def S_f(fx,fy): #6 Rim PSM with a shifter width of 25 nm
#     nor = xx / lambda_ / z
#     if (abs(fx) < 125E-9 * nor) and (abs(fy) < 125E-9 * nor):
#         return 1
#     elif (abs(fx) < 150E-9 * nor) and (abs(fy) < 150E-9 * nor):
#         return -1
#     else:
#         return 0

# def S_f(fx,fy): #7 Outrigger PSM with a shifter distance of C.5 nm and a shifter width of 25 nm
#     nor = xx / lambda_ / z
#     if ((abs(fx) - 175E-9 * nor < 12.5E-9 * nor) and (abs(fy) < 150E-9 * nor)) or ((abs(fy) - 175E-9 * nor < 12.5E-9 * nor) and (abs(fx) < 150E-9 * nor)):
#         return -1
#     elif (abs(fx) < 150E-9 * nor) and (abs(fy) < 150E-9 * nor):
#         return 1
#     else:
#         return 0

# def S_f(fx,fy): #8 Attenuated PSM with a transmittance of 10%
#     nor = xx / lambda_ / z
#     if (abs(fx) < 125E-9 * nor) and (abs(fy) < 125E-9 * nor):
#         return -0.1
#     elif (abs(fx) < 250E-9 * nor) and (abs(fy) < 250E-9 * nor):
#         return 1
#     else:
#         return 0
print("cell 5 finished")

#%%
#mask t(x,y) in x space 
def t_f(x,y):
    if (abs(abs(x) - 0.25E-6) < 0.125E-6) & (abs(abs(y) - 0.25E-6) < 0.125E-6) :
        return 1
    else:
        return 0
    
t = np.zeros((2**N,2**N), dtype=complex)
for i in range(2**N):
    for j in range(2**N):
        t[i][j] = t_f(x[i],x[j])
T = ft(t,x,fx)
plot(np.square(np.absolute(t)),"t")
print("cell 6 finished")

#%%
#calculate fft of pupul, source
def mesh():
    P = np.zeros((2**(N+1),2**(N+1)), dtype=complex)
    S = np.zeros((2**(N+2),2**(N+2)), dtype=complex)

    for i in range(2**(N+1)):
        for j in range(2**(N+1)):
            P[i][j] = P_f(fx_double[i], fx_double[j])
    
    for i in range(2**(N+2)):
        for j in range(2**(N+2)):
            S[i][j] = S_f(fx_quater[i], fx_quater[j])
    
    return P,S
P, S = mesh()
h = ft(P, x_double, fx_double)
J0 = ft(S, x_quater, fx_quater)
# plot(np.square(np.absolute(J0)),"J0")
# plot(np.square(np.absolute(h)),"h")
print("cell 7 finished")

#%%
#calculate fft of pupul, source
#calculate kernal
def calculate():
    W = np.zeros((2**(N+1),2**(N+1),2**(N+1),2**(N+1)), dtype=complex)
    M = np.zeros((2**(N), 2**(N)), dtype = complex)  
    W_squashed = np.zeros((2**(2*N),2**(2*N)), dtype=complex)
    for i in range(2**(N+2)):
        for j in range(2**(N+2)):
            try:
                J0[-i-1][-j-1] = (J0[-i-1][-j-1] + J0[-i-1][-j-2] + J0[-i-2][-j-1] + J0[-i-2][-j-2]) / 4
            except:
                J0[-i-1][-j-1] = J0[-i-1][-j-1]
    
    for x1 in range(2**(N+1)):
        for y1 in range(2**(N+1)):
            for x2 in range(2**(N+1)):
                for y2 in range(2**(N+1)):
                    W[x1][y1][x2][y2] = J0[x2 - x1 + 2**(N+1)][y2 - y1 + 2**(N+1)] * h[x1][y1] * np.conjugate(h[x2][y2])

    for x1 in range(2**(N)):
        for y1 in range(2**(N)):
            for x2 in range(2**(N)):
                for y2 in range(2**(N)):
                    M[x2][y2] = J0[x2 - x1 + 2**(N+1)][y2 - y1 + 2**(N+1)] * h[x1+2**(N-1)][y1+2**(N-1)] * np.conjugate(h[x2+2**(N-1)][y2+2**(N-1)])
                    W_squashed[x1*2**(N) + x2][y1*2**(N) + y2] = M[x2][y2]
        
    return W, W_squashed
W, W_squashed = calculate()
plot(np.square(np.absolute(W_squashed)),"W")    
print("cell 8 finished")

#%%
def intensity():
    I_total = np.zeros((2**N,2**N), dtype=complex)

    for x1 in range(2**(N+1)):
        for y1 in range(2**(N+1)):
            for x2 in range(2**(N+1)):
                for y2 in range(2**(N+1)):
                    try:
                        W[-x1-1][-y1-1][-x2-1][-y2-1]=\
                            (W[-x1-1][-y1-1][-x2-1][-y2-1]+\
                            W[-x1-1][-y1-1][-x2-1][-y2-2]+\
                            W[-x1-1][-y1-1][-x2-2][-y2-1]+\
                            W[-x1-1][-y1-1][-x2-2][-y2-2]+\
                            W[-x1-1][-y1-2][-x2-1][-y2-1]+\
                            W[-x1-1][-y1-2][-x2-1][-y2-2]+\
                            W[-x1-1][-y1-2][-x2-2][-y2-1]+\
                            W[-x1-1][-y1-2][-x2-2][-y2-2]+\
                            W[-x1-2][-y1-1][-x2-1][-y2-1]+\
                            W[-x1-2][-y1-1][-x2-1][-y2-2]+\
                            W[-x1-2][-y1-1][-x2-2][-y2-1]+\
                            W[-x1-2][-y1-1][-x2-2][-y2-2]+\
                            W[-x1-2][-y1-2][-x2-1][-y2-1]+\
                            W[-x1-2][-y1-2][-x2-1][-y2-2]+\
                            W[-x1-2][-y1-2][-x2-2][-y2-1]+\
                            W[-x1-2][-y1-2][-x2-2][-y2-2]) / 16
                            
                    except:
                        W[x1][y1][x2][y2] = W[x1][y1][x2][y2]
                    

    for i in range(2**N):
        for j in range(2**N):
            for x1 in range(2**N):
                for y1 in range(2**N):
                    for x2 in range(2**N):
                        for y2 in range(2**N):
                            I_total[i][j] += W[i-x1+2**N][j-y1+2**N][i-x2+2**N][j-y2+2**N] * t[x1][y1] * np.conjugate(t[x2][y2]) / 2**(4*N)
    return np.square(np.absolute(I_total))
print("cell 9 finished")

#%%
def svd():
    phi = np.zeros((2**N,2**N), dtype=complex)
    I_total = np.zeros((2**N,2**N))
    
    u, s, vh = np.linalg.svd(W_squashed)
    
    for state in range(Nmax):
        for i in range(2**N):
            for j in range(2**N):
                phi[i][j] = vh[state][i*2**N + j]
        phi = ft(phi,x,fx)
        phi = phi * T
        phi = ift(phi,x,fx)
        I_total += np.real(np.square(np.absolute(phi)) * s[state])
    return I_total

#%%
def comfirm4f():
    t = np.zeros((2**N,2**N), dtype=complex)
    P = np.zeros((2**N,2**N), dtype=complex)
    I = np.zeros((2**N,2**N), dtype=complex)

    for i in range(2**N):
        for j in range(2**N):
            t[i][j] = t_f(x[i],x[j])
            P[i][j] = P_f(fx[i],fx[j])
    T = ft(t,x,fx)
    
    # plot(np.square(np.absolute(T)),"T")
    # plot(np.square(np.absolute(P)),"P")
    
    for i in range(2**N):
        for j in range(2**N):
            T[i][j] = P[i][j] * T[i][j]
    
    # plot(np.square(np.absolute(T)),"P*T")
    I = ift(T,x,fx)
    return np.square(np.absolute(I))
print("cell 10 finished")

#%%
def comfirm6f():
    S = np.zeros((2**N,2**N), dtype=complex)
    t = np.zeros((2**N,2**N), dtype=complex)
    P = np.zeros((2**N,2**N), dtype=complex)
    I = np.zeros((2**N,2**N), dtype=complex)

    for i in range(2**N):
        for j in range(2**N):
            t[i][j] = t_f(x[i],x[j])
            S[i][j] = S_f(fx[i],fx[j])
            P[i][j] = P_f(fx[i],fx[j])

    # plot(np.square(np.absolute(t)),"t")
    # plot(np.square(np.absolute(S)),"S")
    # plot(np.square(np.absolute(P)),"P")

    J0 = ift(S,x,fx)
    # plot(np.square(np.absolute(J0)),"J0")

    for i in range(2**N):
        for j in range(2**N):
            J0[i][j] = J0[i][j] * t[i][j]
    # plot(np.square(np.absolute(J0)),"J0*t")
    T = ft(J0,x,fx)
    # plot(np.square(np.absolute(T)),"T")
    for i in range(2**N):
        for j in range(2**N):
            T[i][j] = P[i][j] * T[i][j]
    # plot(np.square(np.absolute(T)),"T*P")
    I = ift(T,x,fx)
    return np.square(np.absolute(I))
print("cell 11 finished")

#%%9
def main():
    Intensity_4f = comfirm4f()
    Intensity_6f = comfirm6f()
    Intensity_W = intensity()
    Intensity_SOCS = svd()
    
    plot(Intensity_4f,"I_4f")
    plot(Intensity_6f,"I_6f")
    plot(Intensity_W,"I_W")
    plot(Intensity_SOCS,"I_SOCS")

    print("--- %s seconds ---" % (time.time() - start_time))
    # print('\a')
main()