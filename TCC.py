import numpy as np
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

    TCC = np.zeros((2**(2*N), 2**(2*N)))
    for i in range(2**(2*N)):
        for j in range(2**(2*N)):
            TCC[i][j] = pool_outputs[2**(2*N)*i + j]

    plot(TCC,"TCC")
    
    print("--- %s seconds ---" % (time.time() - start_time))
    
    np.savetxt("a.txt",TCC)
    


#################################################