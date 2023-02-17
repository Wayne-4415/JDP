#%%
import matplotlib.colors as col
import math 
import cmath
import numpy as np
from numpy import linalg as LA
from math import exp, pi, sqrt, ceil, sin, cos
import matplotlib as mpl
import matplotlib.pyplot as plt
import inspect
import time
from progressbar import *
import multiprocessing as mp
import csv


#%%
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

#calculate kernal
def kernal():

    def T(f1,f2,g1,g2):  
        result = 0
        for i in range(2**N):
            for j in range(2**N):
                # pupil1 = P_f(f1 - i/Lx + 2**(N-1)/Lx , g1 - j/Lx + 2**(N-1)/Lx )
                # P[f1 -i + 2**(N-1)][g1 - j + 2**(N-1)]
                # pupil2 = P_f(f2 - i/Lx + 2**(N-1)/Lx , g2 - j/Lx + 2**(N-1)/Lx )
                # P[f2 -i + 2**(N-1)][g2 - j + 2**(N-1)]
                # result += S_f(i/Lx - 2**(N-1)/Lx, j/Lx - 2**(N-1)/Lx)*pupil1*np.conjugate(pupil2)
                # S[i - 2**(N-1)][j - 2**(N-1)]
                result += P[f1 -i + 2**N -1][g1 - j + 2**N -1]*P[f2 -i + 2**N -1][g2 - j + 2**N -1]*S[i][j]
        return result

    # def TT():
    #     TT = np.zeros((2**(2*N),2**(2*N)))
    #     for i in progress(range(2**(2*N))):
    #         for j in range(2**(2*N)):
    #             # TT[i][j] = T(fx[i-2**N*(ceil(i//2**N))], fx[j-2**N*(ceil(j//2**N))], fx[ceil(i//2**N)], fx[ceil(j//2**N)])
    #             TT[i][j] = T(i-2**N*(ceil(i//2**N)), j-2**N*(ceil(j//2**N)), ceil(i//2**N), ceil(j//2**N))
    #     return TT
    def TT():
        TT = np.zeros((2**(2*N),2**(2*N)))

        # for i in progress(range(2**(4*N))):
        #     TT[i//2**(2*N)][i-2**(2*N)*(i//2**(2*N))] = T((i//2**(2*N))-2**N*((i//2**(2*N))//2**N), (i-2**(2*N)*(i//2**(2*N)))-2**N*((i-2**(2*N)*(i//2**(2*N)))//2**N), (i//2**(2*N))//2**N, (i-2**(2*N)*(i//2**(2*N)))//2**N)
        # return TT
    
def TCC(*var):
    np.savetxt('a.txt', "0", fmt='%f',delimiter=',')
        
#%%
if __name__ == '__main__':
    N = 2
    NA = 0.6
    lambda_ = 248E-9
    z = 1*lambda_
    k = 2*pi/lambda_
    Lx = 2E-6
    Nmax = 10 #use in Kernal function, must be less than N square
    
    x = np.zeros((2**N))
    fx = np.zeros((2**N))
    x_double = np.zeros((2**(N+1))) #for double the source area  //same as double xii
    fx_double = np.zeros((2**(N+1)))
    for i in range(2**N):
        x[i] = -Lx/2 + Lx*i/2**N + Lx/2**(N+1)
        fx[i] = (i-(2**N-1)/2)/Lx
    for i in range(2**(N+1)):
        x_double[i] = -1*Lx + 2*Lx*i/2**(N+1) + Lx/2**(N+1)
        fx_double[i] = (i-(2**(N+1)-1)/2)/Lx
        
    P = np.zeros((2**(N+1),2**(N+1)))
    S = np.zeros((2**(N),2**(N)))
    for i in range(2**(N+1)):
        for j in range(2**(N+1)):
            P[i][j] = P_f(fx_double[i]+1/Lx/2,fx_double[j]+1/Lx/2)
    for i in range(2**(N)):
        for j in range(2**(N)):
            S[i][j] = S_f(fx[i]+1/Lx/2,fx[j]+1/Lx/2)
                
    variation = []
    sep_variation = []
    for i in range(2**(4*N)):
        variation.append([(i//2**(2*N))-2**N*((i//2**(2*N))//2**N), (i-2**(2*N)*(i//2**(2*N)))-2**N*((i-2**(2*N)*(i//2**(2*N)))//2**N), (i//2**(2*N))//2**N, (i-2**(2*N)*(i//2**(2*N)))//2**N])
    
    for i in range(mp.cpu_count()):
        temp=[]
        for j in range(2**(4*N)//mp.cpu_count()):
            temp.append(variation[j + i * (2**(4*N)//mp.cpu_count())])
        sep_variation.append(temp)
            
    for i in range(2**(4*N) - mp.cpu_count() * (2**(4*N)//mp.cpu_count())):
        sep_variation[mp.cpu_count()-1].append(variation[-(2**(4*N) - mp.cpu_count() * (2**(4*N)//mp.cpu_count())) +i])
    
    
    process_list = []
    for i in range(mp.cpu_count()):
        process_list.append(mp.Process(target = TCC, args = sep_variation[i]))
        process_list[i].start()

    for i in range(mp.cpu_count()):
        process_list[i].join()
        

