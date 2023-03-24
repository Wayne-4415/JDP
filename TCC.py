from math import exp, pi
import cmath
import time
import multiprocessing as mp
import os
from config import *

P = np.zeros((2**(N+1), 2**(N+1)))
S = np.zeros((2**(N), 2**(N)))

P = P_f(fxx_double + 1/Lx/2, fyy_double + 1/Lx/2)
S = S_f(fxx, fyy)

def TCC(var):
    f1 = var[0]
    f2 = var[1]
    g1 = var[2]
    g2 = var[3]
    
    return np.sum(P[f1: f1+2**N, g1: g1+2**N] * P[f2: f2+2**N, g2: g2+2**N] * np.transpose(S))
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
    
    print("--- %s seconds ---" % (time.time() - start_time))
    #%%
    TCC4D = np.zeros((2**N,2**N,2**N,2**N))

    for i in range(2**N):
        for j in range(2**N):
            for k in range(2**N):
                for l in range(2**N):
                    TCC4D[j][l][i][k] = pool_outputs[l+2**N*k+2**(2*N)*j+2**(3*N)*i]
                    
    #%%
    t=np.zeros((2**N,2**N), dtype = complex)
    t = t_f(xx, yy)
    T = ft(t, x, fx)
    
    I = np.zeros((2**N,2**N), dtype = complex)
    
    for i in range(2**N):
        for j in range(2**N):
            for fx_1 in range(2**N):
                for fy_1 in range(2**N):
                    for fx_2 in range(2**N):
                        for fy_2 in range(2**N):
                            I[i,j] += TCC4D[fx_1,fy_1,fx_2,fy_2]*T[fx_1,fy_1]*T[fx_2,fy_2]*cmath.exp(1j*2*pi*((fx_1-fx_2)*x[i]+(fy_1-fy_2)*x[j]))

    
    plot(np.square(np.absolute(I)), "I")

#################################################