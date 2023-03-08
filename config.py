import numpy as np
from math import exp, pi, sqrt, ceil, sin, cos
def parameter():
    N = 4
    NA = 0.6
    lambda_ = 248E-9
    z = 1E-2
    k = 2*pi/lambda_
    Lx = 1E-5
    Nmax = 10 #use in Kernal function, must be less than N square
    x = np.zeros((2**N))
    fx = np.zeros((2**N))
    x_double = np.zeros((2**(N+1)))
    fx_double = np.zeros((2**(N+1)))
    for i in range(2**N):
        x[i] = -Lx/2 + Lx*i/2**N + Lx/2**(N+1)
        fx[i] = (i-(2**N-1)/2)/Lx
    for i in range(2**(N+1)):
        x_double[i] = -1*Lx + 2*Lx*i/2**(N+1) + Lx/2**(N+1)
        fx_double[i] = (i-(2**(N+1)-1)/2)/Lx #same resolution as fx but double range

    return N, NA,lambda_, z, k, Lx, Nmax, x, fx, x_double, fx_double

