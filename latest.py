from config import *

#%%
#plot
t=np.zeros((2**N,2**N), dtype = complex)
P=np.zeros((2**N,2**N), dtype = complex)
S=np.zeros((2**N,2**N), dtype = complex)
P_double = np.zeros((2**(N+1), 2**(N+1)), dtype = complex)

t = t_f(xx, yy)
P = P_f(fxx, fyy)
S = S_f(fxx, fyy)
P_double = P_f(fxx_double, fyy_double)

plot(np.absolute(t), "t")
plot(np.absolute(P), "P")
plot(np.absolute(S), "S")

#%%
T = ft(t, x, fx)

#%%
#4f
I_4f = np.zeros((2**N, 2**N), dtype=complex)

I_4f = ift(P * T, x, fx)
plot(np.square(np.absolute(I_4f)), "I_4f")

#%%
#6f
I_6f = np.zeros((2**N, 2**N), dtype=complex)
P_shift = np.zeros((2**N, 2**N), dtype=complex)
normalize = 0

for i in range(2**N):
    for j in range(2**N):
        if S[i][j] != 0:
            P_shift[:, :] = P_double[i:2**N+i, j:2**N+j]
            I_6f += ift(P_shift * T, x, fx)
            normalize += 1
        else:
            continue
            
plot(np.square(np.absolute(I_6f/normalize)), "I_6f")

