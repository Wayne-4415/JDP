import numpy as np
import matplotlib.pyplot as plt
import config
N, NA,lambda_,z,k,Lx,Nmax = config.parameter()

def plot(figure,title):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect('equal', adjustable='box')
    ax.set_title(title)
    contourf = ax.contourf(np.real(figure), 400, vmin=figure.min(), vmax=figure.max(), levels = 50)
    cbar = fig.colorbar(contourf)
    plt.show()

A=np.ndarray.tolist(np.loadtxt("a.txt"))

sorted(A, key = lambda s: s[2])
A.sort(key = lambda s: s[2])

for i in range(2**(N)):
    temp = A[2**(3*N)*i:2**(3*N)*(i+1)]
    temp.sort(key = lambda s: s[0])
    A[2**(3*N)*i:2**(3*N)*(i+1)] = temp
for i in range(2**(2*N)):
    temp = A[2**(2*N)*i:2**(2*N)*(i+1)]
    temp.sort(key = lambda s: s[3])
    A[2**(2*N)*i:2**(2*N)*(i+1)] = temp
for i in range(2**(3*N)):
    temp = A[2**(1*N)*i:2**(1*N)*(i+1)]
    temp.sort(key = lambda s: s[1])
    A[2**(1*N)*i:2**(1*N)*(i+1)] = temp

TCC = np.zeros((2**(2*N),2**(2*N)))
for i in range(2**(2*N)):
    for j in range(2**(2*N)):
        try:
            TCC[i][j] = A[i*2**(2*N)+j][4]
        except:
            TCC[i][j] = 1000
plot(TCC,"TCC")

