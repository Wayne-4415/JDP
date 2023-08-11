import numpy as np
from math import pi, sqrt, exp
import cmath
import matplotlib.pyplot as plt

N=6
Lx = 1E-6
Lfx = 2**4/Lx
x = np.zeros((2**N))
fx = np.zeros((2**N))

for i in range(2**N):
    x[i] = -Lx/2 + Lx*i/(2**N)
    fx[i] = -Lfx/2 + Lfx*i/(2**N)

def plot(figure,title):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # ax.set_aspect('equal', adjustable='box')
    ax.set_title(title)
    contourf = ax.contourf(np.real(figure), 400, vmin=figure.min(), vmax=figure.max(), levels = 50, cmap = "gray")
    cbar = fig.colorbar(contourf)
    plt.show()

def fold(M):
    M_fold = np.zeros((2**N,2**N), dtype =complex)
    M_fold[0:2**(N-1),0:2**(N-1)] = M[2**(N-1):2**N,2**(N-1):2**N]
    M_fold[0:2**(N-1),2**(N-1):2**N] = M[2**(N-1):2**N,0:2**(N-1)]
    M_fold[2**(N-1):2**N,0:2**(N-1)] = M[0:2**(N-1),2**(N-1):2**N]
    M_fold[2**(N-1):2**N,2**(N-1):2**N] = M[0:2**(N-1),0:2**(N-1)]
    return M_fold

xx, yy = np.meshgrid(x, x)
fxx, fyy = np.meshgrid(fx, fx)
def ft(f,x,fx):
    temp = np.zeros((len(fx), len(fx)), dtype=complex)
    for i in range(len(fx)):
        for j in range(len(fx)):
            temp[i,j] = np.sum(f*np.exp(-1j*2*pi*(xx*fx[i]+yy*fx[j])))
    return temp

def ift(f,x,fx):
    temp = np.zeros((len(x), len(x)), dtype=complex)
    for i in range(len(x)):
        for j in range(len(x)):
            temp[j,i] = np.sum(f*np.exp(1j*2*pi*(fxx*x[i]+fyy*x[j])))
    return temp

def P_shift(x_shift,y_shift):
    P_shift = np.zeros((2**N,2**N))
    P_shift[max(x_shift,0):min(2**N,2**N+x_shift),max(y_shift,0):min(2**N,2**N+y_shift)] = P[max(0,-x_shift):min(2**N-x_shift,2**N),max(0,-y_shift):min(2**N-y_shift,2**N)]
    return P_shift

"""
define source, pupil adn mask
"""
S = np.zeros((2**N,2**N))
for i in range(2**N):
    for j in range(2**N):
        if (abs(fx[i])-0.5*(0.25*Lfx))**2+(abs(fx[j])-0.6*(0.25*Lfx))**2 < 0.017*(0.25*Lfx)**2:
            S[i,j] = 1
        if (abs(fx[i])-0.6*(0.25*Lfx))**2+(abs(fx[j])-0.5*(0.25*Lfx))**2 < 0.017*(0.25*Lfx)**2:
            S[i,j] = 1
        
P = np.zeros((2**N,2**N))
for i in range(2**N):
    for j in range(2**N):
        if fx[i]**2+fx[j]**2 < (0.25*Lfx)**2:
            P[i,j] = 1
            
            
a = np.zeros((2**N,2**N))
for i in range(2**N):
    for j in range(2**N):
        if abs(abs(x[i])-100E-9) < 50E-9 and abs(x[j]) < 50E-9:
            a[i,j] = 1

eig_num = int(np.sum(S))
print("havs " + str(eig_num) + " eigenvalues")

P_operator = np.zeros((eig_num,2**(2*N)))
s=0
for i in range(2**N):
    for j in range(2**N):
        if S[i,j] == 1:
            P_operator[s:s+1,:] = np.transpose(P_shift(-i+2**(N-1),-j+2**(N-1))).reshape(1,2**(2*N))
            s+=1

U, lamb, Vh = np.linalg.svd(P_operator, full_matrices=False)
lamb = lamb/max(lamb)
print("first 12 square of singular values are \n" + str(np.square(lamb[0:12])))

phi = np.zeros((eig_num,2**N,2**N))
for i in range(eig_num):
    phi[i,:,:] = Vh[i,:].reshape(2**N,2**N)
    
ift_a = np.zeros((2**N,2**N), dtype = complex)
ift_a = ift(a, x, fx)

eig_I = np.zeros((eig_num,2**N,2**N), dtype = complex)
for i in range(9):# for i in range(eig_num):
    eig_I[i,:,:] = ft(ift_a * phi[i,:,:], x, fx)

I = np.zeros((2**N,2**N))
for i in range(9):# for i in range(eig_num):
    I += lamb[i]**2 * np.square(np.absolute(eig_I[i,:,:]))


plot(P_operator,"P_operator")
plot(S,"S")
plot(P,"P")
plot(a,"a")
plot(np.real(ift_a),"ift_a")
plot(phi[0,:,:],"phi0")
plot(phi[1,:,:],"phi1")
plot(phi[2,:,:],"phi2")
plot(np.real(eig_I[0,:,:]),"eig_I0")
plot(np.imag(eig_I[1,:,:]),"eig_I1")
plot(np.imag(eig_I[2,:,:]),"eig_I2")
plot(I,"I")