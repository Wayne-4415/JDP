import math 
import cmath
import numpy as np
from math import exp, pi, sqrt
import time
import multiprocessing as mp
import matplotlib.pyplot as plt
import os

import config
N, NA,lambda_, z, k, Lx, Nmax, x, fx, x_double, fx_double = config.parameter()

def P_f(fx,fy):
    if sqrt((fx)**2 + (fy)**2) < NA/lambda_:
        return 1
    else:
        return 0

#source in f space

def S_f(fx,fy): #1 Conventional illumination
    if sqrt(fx**2 + fy**2) < 0.6*NA/lambda_:
        return 1
    else:
        return 0

def plot(figure,title):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect('equal', adjustable='box')
    ax.set_title(title)
    contourf = ax.contourf(np.real(figure), 400, vmin=figure.min(), vmax=figure.max(), levels = 50)
    cbar = fig.colorbar(contourf)
    plt.show()


P = np.zeros((2**(N+1),2**(N+1)))
S = np.zeros((2**(N),2**(N)))
for i in range(2**(N+1)):
    for j in range(2**(N+1)):
        P[i][j] = P_f(fx_double[i]+1/Lx/2,fx_double[j]+1/Lx/2)
for i in range(2**(N)):
    for j in range(2**(N)):
        S[i][j] = S_f(fx[i],fx[j])



def TCC(var):
    f1 = var[0]
    f2 = var[1]
    g1 = var[2]
    g2 = var[3]
    result = 0
    for i in range(2**N):
        for j in range(2**N):
            result += P[f1 -i + 2**N -1][g1 - j + 2**N -1]*P[f2 -i + 2**N -1][g2 - j + 2**N -1]*S[i][j]
    return result

#################################################
if __name__ == '__main__':
    try:
        os.remove(r'a.txt')
    except:
        pass
    start_time = time.time()
    variation = []
    for i in range(2**(4*N)):
        variation.append([(i//2**(2*N))-2**N*((i//2**(2*N))//2**N), (i-2**(2*N)*(i//2**(2*N)))-2**N*((i-2**(2*N)*(i//2**(2*N)))//2**N), (i//2**(2*N))//2**N, (i-2**(2*N)*(i//2**(2*N)))//2**N])
    
    pool = mp.Pool(mp.cpu_count())
    pool_outputs = pool.map(TCC, variation)

    TCC = np.zeros((2**(2*N),2**(2*N)))
    for i in range(2**(2*N)):
        for j in range(2**(2*N)):
            TCC[i][j] = pool_outputs[2**(2*N)*i + j]

    plot(TCC,"TCC")
    
    print("--- %s seconds ---" % (time.time() - start_time))
    
    np.savetxt("a.txt",TCC)
    


#################################################