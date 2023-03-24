import numpy as np
from config import *

#%%
#plot
t=np.zeros((2**N,2**N), dtype = complex)
P=np.zeros((2**N,2**N), dtype = complex)
S=np.zeros((2**N,2**N), dtype = complex)
t = t_f(xx, yy)
P = P_f(fxx, fyy)
S = S_f(fxx, fyy)
plot(np.absolute(t),"t")
plot(np.absolute(P),"P")
plot(np.absolute(S),"S")

#%%
#4f
I_4f = np.zeros((2**N,2**N), dtype=complex)

T = ft(t,x,fx)
T = P * T
I_4f = ift(T,x,fx)
plot(np.square(np.absolute(I_4f)),"I_4f")

#%%
#6f
I_6f = np.zeros((2**N,2**N), dtype=complex)
P_double = np.zeros((2**(N+1),2**(N+1)), dtype = complex)
TP = np.zeros((2**N,2**N), dtype=complex)
sigma=0.6

P_double = P_f(fxx_double,fyy_double)

T = ft(t,x,fx)

for i in range(2**N):
    for j in range(2**N):
        if S[i][j] != 0:
            P[:, :] = P_double[i:2**N+i, j:2**N+j]
            TP = P * T
            I_6f += ift(TP,x,fx) 
        else:
            continue
            
plot(np.square(np.absolute(I_6f)),"I_6f")

