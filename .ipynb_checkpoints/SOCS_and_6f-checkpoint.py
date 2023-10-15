from config import *
import time

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


# plot(np.absolute(t), "t")
# plot(np.absolute(P), "P")
# plot(np.absolute(S), "S")

#%%
T = ft(t, x, fx)
# plot(np.square(np.absolute(T)), "T")

#%%
#6f
start_time = time.time()
I_6f_moveS = np.zeros((2**N, 2**N))
I_6f_moveP = np.zeros((2**N, 2**N))

normalize = 0

for i in range(2**N):
    for j in range(2**N):
        if S[i, j] != 0:
            S_part = np.zeros((2**N, 2**N), dtype=complex)
            S_part[i, j] = 1 
            # I_6f_moveS += np.square(np.absolute(ift(P * ft(t * ift(S_part, x, fx), x, fx), x, fx)))
            
            I_6f_moveP += np.square(np.absolute(ift(P_double[i:2**N+i, j:2**N+j] * T, x, fx)))

            normalize += 1
        else:
            continue

# I_6f_moveS = I_6f_moveS/normalize
I_6f_moveP = I_6f_moveP/normalize

# plot(I_6f_moveS, "I_6f_moveS")
plot(I_6f_moveP, "I_6f_moveP")
print("---6f process time: %s seconds ---" % (time.time() - start_time))
#%%
#SOCS
start_time = time.time()
TCC=(np.loadtxt("a.txt"))

I_SOCS = np.zeros((2**N, 2**N))
state = 0

u, s, vh = np.linalg.svd(TCC)

while s[state] > 0.01 * s[0]:
    I_SOCS += np.square(np.absolute(ift(np.reshape(vh[state,:], (2**N,2**N)) * T, x, fx))) * s[state]
    state += 1

plot(I_SOCS,"I_SOCS")
print("---SOCS process time: %s seconds ---" % (time.time() - start_time))
#%%

err = np.zeros((2**N, 2**N))
err = np.square(np.absolute(I_SOCS - I_6f_moveP))

plot(err,"err_SOCS-6f")

