from config import *

TCC=(np.loadtxt("a.txt"))

def kernal():
    u, s, vh = np.linalg.svd(TCC)
    kernal = np.zeros((Nmax, 2**N, 2**N), dtype=complex)
    kernal_val = np.zeros((Nmax), dtype=complex)
    for state in range(Nmax):
        kernal[state,:,:] = np.reshape(vh[state,:], (2**N,2**N))
        kernal_val[state] = s[state]
    return kernal ,kernal_val

def socs():
    phi = np.zeros((2**N, 2**N), dtype=complex)
    I_total = np.zeros((2**N, 2**N))
    t = np.zeros((2**N, 2**N), dtype=complex)
    
    t = t_f(xx, yy)
    T = ft(t, x, fx)
    
    for state in range(Nmax):
        phi[:,:] = kernal[state, :, :]
        # plot(np.absolute(phi),"kernal")
        phi = ift(phi * T, x, fx)
        I_total += np.real(np.square(np.absolute(phi)) * kernal_val[state])
    return I_total

kernal ,kernal_val = kernal()
I_SOCS = socs()
plot(np.absolute(I_SOCS),"I_SOCS")