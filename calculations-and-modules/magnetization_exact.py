import warnings 
warnings.filterwarnings("ignore")

from redfield_functions import *
from utils import *
from circuits import *

tmin = 0
tmax = 10**(-5)
n_t = 20
t_list = tlist(tmin, tmax, n_t)

# NOTE: change magnetic field B and/or temperature on your own choice
B = 1
temperature = 25

H = Hamiltonian(B)
ev, V = np.linalg.eig(H)

s = np.zeros((1, 2, 2))
s[0] = sigx

spectrum = lambda w: ohmic_spectrum(w, temperature)

R = redfield(H, s, spectrum)

magnetization_exact = []
for i in range(len(t_list)):
    ori = scipy.linalg.expm(R*t_list[i])@np.array([[0, 0, 0, 1]]).reshape(-1, 1)
    ori = ori.reshape(H.shape)
    ori /= np.trace(ori)
    mag = np.trace(ori@sigz)
    magnetization_exact.append(mag)
    
np.save('magnetization_exact'+str(B)+'T'+str(temperature)+'K.npy', magnetization_exact)