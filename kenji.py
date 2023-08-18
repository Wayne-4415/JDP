import numpy as np
from math import pi, sqrt, exp
import cmath
import matplotlib.pyplot as plt
import time 

####################################
'''
define physic parameters and commonly used functions.

'''
N=8
Lx = 1E-6
Lfx = 4
x = np.zeros((2**N))
fx = np.zeros((2**N))

for i in range(2**N):
    x[i] = -Lx/2 + Lx*i/(2**N)
    fx[i] = -Lfx/2 + Lfx*i/(2**N)

def plot(figure,title):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect('equal', adjustable='box')
    ax.set_title(title)
    contourf = ax.contourf(np.real(figure), 400, vmin=figure.min(), vmax=figure.max(), levels = 50, cmap = "gray")
    # fig.colorbar(contourf)
    plt.show()

def paddle(M,n):
    temp = np.zeros((n, n), dtype=complex)
    temp[n//2-len(M)//2:n//2+len(M)//2, n//2-len(M)//2:n//2+len(M)//2] = M
    return temp

def capture(M,n):
    temp = np.zeros((n, n), dtype=complex)
    temp = M[len(M)//2-n//2:len(M)//2+n//2, len(M)//2-n//2:len(M)//2+n//2]
    return temp

def P_shift(x_shift,y_shift):
    P_shift = np.zeros((2**N,2**N))
    P_shift[max(x_shift,0):min(2**N,2**N+x_shift),max(y_shift,0):min(2**N,2**N+y_shift)] = P[max(0,-x_shift):min(2**N-x_shift,2**N),max(0,-y_shift):min(2**N-y_shift,2**N)]
    return P_shift

####################################
'''
define source S, pupil P and mask a.

'''
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

####################################
'''
define parameters in main program.

th : the threshold of eigenvalues that calculated
P_operator : P_operator
ift_a : inverse fourier transform of a
eig : eigenfunction, produced by P_operator using SVD
eig_I : eigenfunction of aerial image
I : aerial image

'''
eig_num = int(np.sum(S))
print("havs " + str(eig_num) + " eigenvalues\n")

th = 0.1
P_operator = np.zeros((eig_num,2**(2*N)))
ift_a = np.zeros((2**N,2**N), dtype = complex)
eig = np.zeros((eig_num,2**N,2**N))
eig_I = np.zeros((eig_num,2**N,2**N), dtype = complex)
I = np.zeros((2**N,2**N))

####################################
'''
main program.

'''
start_time = time.time()
s=0
for i in range(2**N):
    for j in range(2**N):
        if S[i,j] == 1:
            P_operator[s:s+1,:] = np.transpose(P_shift(-i+2**(N-1),-j+2**(N-1))).reshape(1,2**(2*N))
            # P_operator[s:s+1,:] = P_shift(-i+2**(N-1),-j+2**(N-1)).reshape(1,2**(2*N))
            s+=1
print("calculate p_orerator time: %s seconds ---\n" % (time.time() - start_time))
start_time = time.time()

U, lamb, Vh = np.linalg.svd(P_operator, full_matrices=False)
lamb = lamb/max(lamb)
print("first 12 square of singular values are \n" + str(np.square(lamb[0:12])) + "\n")
print("threshold = " + str(th) +", calculate " + str(sum(a > th for a in lamb)) + " eigenvalues\n")
print("calculate SVD time: %s seconds ---" % (time.time() - start_time))
start_time = time.time()

ift_a = capture(np.fft.ifftshift(np.fft.ifft2(paddle(a, 2**(2*N-4)))), 2**N)
print("calculate ft of mask time: %s seconds ---" % (time.time() - start_time))
start_time = time.time()

for i in range(sum(a > th for a in lamb)):
    eig[i,:,:] = Vh[i,:].reshape(2**N,2**N)
    eig_I[i,:,:] = capture(np.fft.fft2(paddle(ift_a * eig[i,:,:], 2**(2*N-4))), 2**N)
    I += lamb[i]**2 * np.square(np.absolute(eig_I[i,:,:]))

####################################
'''
plot.

'''
# plot(P_operator,"P_operator")
plot(S,"S")
# plot(P,"P")
plot(a,"a")
plot(np.absolute(ift_a),"ift_a")
plot(eig[0,:,:],"phi0")
plot(eig[1,:,:],"phi1")
plot(eig[2,:,:],"phi2")
plot(np.absolute(eig_I[0,:,:]),"eig_I0")
plot(np.absolute(eig_I[1,:,:]),"eig_I1")
plot(np.absolute(eig_I[2,:,:]),"eig_I2")
plot(I,"I")

####################################
print("calculate aerial image time: %s seconds ---" % (time.time() - start_time))