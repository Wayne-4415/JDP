import time
import multiprocessing as mp
import os
from config import *

P = np.zeros((2**(N+1), 2**(N+1)))
S = np.zeros((2**(N), 2**(N)))

P = P_f(fxx_double, fyy_double)
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

    TCC = np.zeros((2**(2*N), 2**(2*N)))
    TCC=np.reshape(pool_outputs / np.sum(S), (2**(2*N), 2**(2*N)))
    
    plot(TCC,"TCC")
    
    print("--- %s seconds ---" % (time.time() - start_time))
    
    np.savetxt("a.txt",TCC)
    
    


#################################################