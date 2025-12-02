# Mapping Bloch-Redfield dynamics into a unitary gate-based quantum algorithm

- A quantum algorithm to simulate the Redfield dynamics of a single molecule magnets.

## Dependencies

-  **Python:**
  - NumPy == 1.26.4
  - qiskit == 1.2.2

## Repository structure

├── calculations-and-modules
│   ├── __pycache__
│   ├── circuits.py
│   ├── dynamics_exact.py
│   ├── dynamics_simulator.py
│   ├── magnetization_exact.py
│   ├── magnetization_simulator.py
│   ├── redfield_functions.py
│   └── utils.py
├── raw-data-and-figures
│   ├── figure-3
│   ├── figure-4
│   ├── figure-5
│   └── figure-6
└── README.md

- **`calculations-and-modules/`** 
  - *.py* Python functions to perform various calculations
  - *`dynamics_X.py`* scripts that run dynamics calculations
  - *`magnetization_X.py`* scripts that run magnetization calculations
- **`raw-data-and-figures/`** 
  - Folder per figure containing
    - *.npy* raw data for manuscript figure
    - *.ipynb* Jupyter notebooks for visualization

---

## Notes on the calculation scripts

### Note - 1 

- Calculation scripts (calculations-and-modules/XXX.py) **require** IONQ token to use IONQ's Aria-1 noisy and noise-less simulator.
- To obtain a free token visit: https://ionq.com/quantum-cloud
- To insert your token:
  1. Open the file.
  2. Locate the section marked with:
     ```python
     provider = IonQProvider('XXX') # enter your token to access IONQ's simulators
     ```
  3. Modify the relevant line.


### Note - 2
- Calculation scripts contain optional lines to use different **external magnetic field** and **temperature values**.
- By default, the temperature and the external magnetic field are fixed to 25K and 1T respectively.
- To change the values:
  1. Open the file.
  2. Locate the section marked with:
     ```python
     # NOTE: change magnetic field B and/or temperature on your own choice
     # B = 1
     # temperature = 25
     ```
  3. Modify the relevant line.


### Note - 3

- Calculation scripts contain another optional lines to switch between **noise-less** and **noisy** simulators of **Aer** and **IONQ**.
- By default, the **Aer noise-less simulator** is used.
- To run with the noise-less simulator:
  1. Open the file.
  2. Locate the section marked with:
     ```python
      # NOTE: comment out depending on if you want to use noise-less or noisy simulator

      # simulator = AerSimulator(seed_simulator=42)
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
     ```
  3. Comment/uncomment the relevant line.

  
## Running calculations or doing visualization

- python3 calculations-and-modules/XXX.py 
- jupyter notebook raw-data-and-figures/figure#/XXX.ipynb
