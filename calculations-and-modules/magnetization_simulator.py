import warnings 
warnings.filterwarnings("ignore")

from redfield_functions import *
from utils import *
from circuits import *

from qiskit import QuantumCircuit
from qiskit.compiler import transpile
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel
from qiskit.providers.fake_provider import GenericBackendV2
from qiskit_ionq import IonQProvider

n_shots = 2**10


# NOTE: comment out depending on if you want to use noise-less or noisy simulator

simulator = AerSimulator(seed_simulator=42)

# For noisy Aer simulator comment out the following three lines, default is noise-less Aer

# coupling_map = [(0, 1), (1, 2), (2, 3), (3, 4)]
# device = GenericBackendV2(num_qubits=5, coupling_map=coupling_map, seed=54)
# noise_model = NoiseModel.from_backend(device)
# simulator = AerSimulator(seed_simulator=42, noise_model=noise_model)

# For IONQ uncomment/comment out the following three lines depending on the preference of noisy of noise-less simulator

# provider = IonQProvider('XXX') # enter your token to access IONQ's simulators
# simulator.set_options(noise_model="aria-1") 
# simulator = provider.get_backend("ionq_simulator")
# simulator.set_options(sampler_seed=12345) 

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
eigR, T = np.linalg.eig(R)

R_diag = np.diag(eigR)

sigmap = upper_block_dilation(R_diag, t_list)

T_inv = np.linalg.inv(T)
back_gs = T_inv@np.array([[0, 0, 0, 1]]).reshape(-1, 1)
back_gs_norm = np.linalg.norm(back_gs)
back_gs_normalized = back_gs/back_gs_norm

soln_for_state_init= scipy.optimize.minimize(opt, np.zeros(6), args=(back_gs_normalized, ), method='L-BFGS-B')

magnetization_circuits = magnetization_circuit_func(sigmap, T, soln_for_state_init.x[3])

counts_list = []
for i in range(len(t_list)):
    compiled_circuit = transpile(magnetization_circuits[i], simulator, basis_gates=['rz', 'ry', 'rx', 'h', 'cx'])
    result = simulator.run(compiled_circuit, shots=n_shots).result()
    counts = result.get_counts()
    counts_list.append(counts)

magnetization_simulator = magnetization_expval(counts_list)

np.save('magnetization_simulator'+str(B)+'T'+str(temperature)+'K.npy', magnetization_simulator)