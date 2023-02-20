from math import exp, pi, sqrt, ceil, sin, cos
def parameter():
    N = 4
    NA = 0.6
    lambda_ = 248E-9
    z = 1*lambda_
    k = 2*pi/lambda_
    Lx = 2E-6
    Nmax = 10 #use in Kernal function, must be less than N square
    return N, NA,lambda_,z,k,Lx,Nmax