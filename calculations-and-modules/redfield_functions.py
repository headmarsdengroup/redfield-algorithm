import numpy as np
import scipy

sigz = np.matrix([[1, 0],[0, -1]]) 
sigy = np.matrix([[0, -1j], [1j, 0]])
sigx = np.matrix([[0, 1], [1, 0]])
iden = np.matrix([[1, 0], [0, 1]])

e = -1.602176634E-19 #Coulombs
hbar = 1.054571817E-34 #J*s
kb = 1.380649E-23 #J/K
c_ = 2.998E10 #cm/s
mass = 9.1093837E-31 #kg
mu = (e*hbar)/(2*mass) #J/T
g_zz = 1.988
alph = 10**(-7)

def Hamiltonian(B):
    return (mu*B*g_zz)*(1/(hbar))*sigz*2*np.pi

def coth(theta):
    x = -np.sinh(2*theta)/(1 - np.cosh(2*theta))
    return x

def ohmic_spectrum(w, T): 
    if w == 0.0: 
        return alph*4 
    if w < 0:
        return np.exp((-hbar*(-w))/(kb*T))*alph*(-w)*coth((hbar*(-w))/(2*kb*T))
    elif w > 0:
        return alph*w*coth((hbar*w)/(2*kb*T))
    
def q_jk(j, k, C, H, s):
    dim = len(H)
    ev, V = np.linalg.eig(H)
    s_eg = V.conjugate().T@s[k]@V
    Q = np.zeros((dim, dim), dtype=complex)
    if j == k:
        for n in range(dim):
            for m in range(dim):
                Q[n, m] = (V[:, n].conjugate().T@s_eg@V[:, m])*1/2*C(ev[m] - ev[n])
    else:
        for n in range(dim):
            for m in range(dim):
                Q[n, m] = 0
    return Q

def q_jk_hat(j, k, C, H, s):
    dim = len(H)
    ev, V = np.linalg.eig(H)
    s_eg = V.conjugate().T@s[k]@V
    Q = np.zeros((dim, dim), dtype=complex)
    if j == k:
        for n in range(dim):
            for m in range(dim):
                Q[n, m] = (V[:, n].conjugate().T@s_eg@V[:, m])*1/2*C(ev[n] - ev[m])
    else:
        for n in range(dim):
            for m in range(dim):
                Q[n, m] = 0
    return Q

def redfield(H, s, C):
    ev, V = np.linalg.eig(H)
    dim = len(H)
    R1 = 1j*(np.kron(H.T, np.identity(dim)) - np.kron(np.identity(dim), H))
    R2 = np.zeros((dim**2, dim**2), dtype=complex)
    R3 = np.zeros((dim**2, dim**2), dtype=complex)
    for j in range(dim-1):
        for k in range(dim-1):
            q_term = V@q_jk(j, k, C, H, s)@V.conjugate().T
            q_hat_term = V@q_jk_hat(j, k, C, H, s).T@V.T
            s_term = s[j]
            R2 += np.kron(-np.identity(dim), s_term@q_term) + np.kron(s_term.T, q_term)
            R3 += np.kron(s_term.T@q_hat_term, np.identity(dim)) - np.kron(q_hat_term, s_term)
    return R1 + R2 - R3