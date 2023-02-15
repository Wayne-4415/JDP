#%%
import matplotlib.colors as col
import math 
import numpy as np
from math import exp, pi, sqrt, ceil, sin, cos
import matplotlib as mpl
import matplotlib.pyplot as plt

#全域變數
N = 5
NA = 0.6
lambda_ = 248E-9
#sigma = 0.1
z = 1*lambda_
k = 2*pi/lambda_
Lx = 1E-6
Nmax = N**2 #use in Kernal function, must be less than N square
xx = 0.5 #adjustment of f-space scalxxe

#%%
#matrix fold function
def reshape(M):
    M_len = math.log(len(M),2)
    M_len = int(M_len)
    temp = np.zeros((2**M_len,2**M_len),dtype=complex)
    for i in range(2**(M_len - 1)):
        for j in range(2**(M_len - 1)):
            temp[i][j] = M[2**(M_len - 1) - i - 1][2**(M_len - 1) - j - 1]
            temp[i + 2**(M_len - 1)][j] = M[2**(M_len) - i - 1][2**(M_len - 1) - j - 1]
            temp[i][j + 2**(M_len - 1)] = M[2**(M_len - 1) - i - 1][2**(M_len) - j - 1]
            temp[i + 2**(M_len - 1)][j + 2**(M_len - 1)] = M[2**(M_len) - i - 1][2**(M_len) - j - 1]
    for i in range(2**(M_len)):
        for j in range(2**(M_len)):
            temp[i][j] = temp[i][j]*(-1)**(i + j + 1)
    return temp

'''
#%%
def plot(figure):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect('equal', adjustable='box')
    plt.contourf(np.real(figure), levels = 50)
    plt.show()
'''

#%%
def plot(figure):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect('equal', adjustable='box')
    contourf = ax.contourf(np.real(figure), 400, vmin=figure.min(), vmax=figure.max(), levels = 50)
    cbar = fig.colorbar(contourf)
    plt.show()

#%%
#pupil and source function
#pupil function P(fx,fy) in f space
def P_f(fx,fy):
    if sqrt(fx**2 + fy**2) < Lx**2*NA/(lambda_**2*z*2**N):
        #temp1 = 1j*2*pi*z*lambda_*(fx**2 + fy**2)/2
        #result = cm.exp(temp1)
        #return result
        return 1
    else:
        return 0

#source in f space

def S_f(fx,fy): #1 Conventional illumination
    if sqrt(fx**2 + fy**2) < 0.7*Lx**2*NA/(lambda_**2*z*2**N):
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


#%%
#calculate fft of pupul, source
#calculate kernal
def calculate():
    x = np.zeros((2**N))
    fx = np.zeros((2**N))
    x_double = np.zeros((2**(N+1))) #for double the source area  //same as double xii
    fx_double = np.zeros((2**(N+1)))

    P = np.zeros((2**N,2**N), dtype=complex)
    S = np.zeros((2**(N+1),2**(N+1)), dtype=complex)

    for i in range(2**N):
        x[i] = -1*Lx/2 + Lx*i/2**N + Lx/2**(N+1)
        fx[i] = x[i]*xx/lambda_/z

    for i in range(2**N):
        for j in range(2**N):
            P[i][j] = P_f(fx[i],fx[j]) 

    for i in range(2**(N+1)):
        x_double[i] = -1*Lx + 2*Lx*i/2**(N+1) + Lx/2**(N+1)
        fx_double[i] = x_double[i]*xx/lambda_/z

        for i in range(2**(N+1)):
            for j in range(2**(N+1)):
                S[i][j] = S_f(fx_double[i],fx_double[j])
    h = reshape(np.fft.fft2(P, norm = "ortho"))
    J0 = reshape(np.fft.fft2(S, norm = "ortho"))

#h J0 is the fft result from P,S

#calculate kernal
    M = np.zeros((2**N,2**N), dtype=complex)
    W = np.zeros((2**(2*N),2**(2*N)), dtype=complex)
    phi = np.zeros((2**N,2**N), dtype=complex)
    temp = np.zeros((2**N,2**N), dtype=complex)
    kernal_state = np.zeros((Nmax,2**N,2**N), dtype=complex)
    kernal_Val = np.zeros((Nmax), dtype=complex)
    
    for x1 in range(2**N):
        for y1 in range(2**N):
            for x2 in range(2**N):
                for y2 in range(2**N):
                    M[x2][y2] = J0[x2 - x1 + 2**N -1][y2 - y1 + 2**N -1] * h[x1][y1] * np.conjugate(h[x2][y2])
                    W[x1*2**N + x2][y1*2**N + y2] = M[x2][y2]

    u, s, vh = np.linalg.svd(W) #s = Val u = Vec
    for i in range(Nmax):
        print(s[i])

    for state in range(Nmax):
        for i in range(2**N):
            for j in range(2**N):
                phi[i][j] = vh[state][i*2**N + j]
        temp = reshape(np.fft.fft2(phi, norm = "ortho")) 
        for i in range(2**N):
            for j in range(2**N):
                kernal_state[state][i][j] = temp[i][j]        
        kernal_Val[state] = s[state]
    return kernal_state, kernal_Val

kernal_state, kernal_Val = calculate()

#%%
def intensity():
    x = np.zeros((2**N))
    t = np.zeros((2**N,2**N), dtype=complex)
    phi_t = np.zeros((2**N,2**N), dtype=complex)
    I_total = np.zeros((2**N,2**N))
    
    for i in range(2**N):
        x[i] = -1*Lx/2 + Lx*i/2**N + Lx/2**(N+1)

    for i in range(2**N):
        for j in range(2**N):
            t[i][j] = t_f(x[i],x[j])
    T = reshape(np.fft.ifft2(t, norm = "ortho"))

    for state in range(Nmax):
        for i in range(2**(N)):
            for j in range(2**(N)):
                phi_t[i][j] = kernal_state[state][i][j] * T[i][j]

        phi_t = np.fft.ifft2(reshape(phi_t), norm = "ortho")
        I_total += np.real(np.square(np.absolute(phi_t)) * kernal_Val[state])
    return I_total

#%%
#mask t(x,y) in x space 
def t_f(x,y):
    if (abs(abs(x) - 0.25E-6) < 0.125E-6) & (abs(abs(y) - 0.25E-6) < 0.125E-6) :
    #if (sqrt(x**2+y**2) < 2.0 * lambda_):    
        return 1
    else:
        return 0

#%%
def comfirm6f():
    x = np.zeros((2**N))
    fx = np.zeros((2**N))
    S = np.zeros((2**N,2**N), dtype=complex)
    t = np.zeros((2**N,2**N), dtype=complex)
    P = np.zeros((2**N,2**N), dtype=complex)
    I = np.zeros((2**N,2**N), dtype=complex)
    
    for i in range(2**N):
        x[i] = -1*Lx/2 + Lx*i/2**N + Lx/2**(N+1)
        fx[i] = x[i]*xx/lambda_/z

    for i in range(2**N):
        for j in range(2**N):
            t[i][j] = t_f(x[i],x[j])
            S[i][j] = S_f(fx[i],fx[j])
            P[i][j] = P_f(fx[i],fx[j])
    J0 = np.fft.fft2(reshape(S), norm = "ortho")
    for i in range(2**N):
        for j in range(2**N):
            J0[i][j] = J0[i][j] * t[i][j]
    
    T = reshape(np.fft.ifft2(J0, norm = "ortho"))
    for i in range(2**N):
        for j in range(2**N):
            T[i][j] = P[i][j] * T[i][j]
    I = np.fft.fft2(reshape(T), norm = "ortho")
    plot(np.square(np.absolute(I)))

#%%
def comfirm4f():
    x = np.zeros((2**N))
    fx = np.zeros((2**N))
    t = np.zeros((2**N,2**N), dtype=complex)
    P = np.zeros((2**N,2**N), dtype=complex)
    I = np.zeros((2**N,2**N), dtype=complex)
    
    for i in range(2**N):
        x[i] = -1*Lx/2 + Lx*i/2**N + Lx/2**(N+1)
        fx[i] = x[i]*xx/lambda_/z

    for i in range(2**N):
        for j in range(2**N):
            t[i][j] = t_f(x[i],x[j])
            P[i][j] = P_f(fx[i],fx[j])
    T = reshape(np.fft.ifft2(t, norm = "ortho"))
    for i in range(2**N):
        for j in range(2**N):
            T[i][j] = P[i][j] * T[i][j]
    I = np.fft.fft2(reshape(T), norm = "ortho")
    plot(np.square(np.absolute(I)))

#%%9
def main():
    Intensity = intensity()
    plot(Intensity)

    comfirm4f()
    comfirm6f()
main()


