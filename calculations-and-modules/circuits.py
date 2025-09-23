from qiskit import QuantumCircuit
import numpy as np

def propagation_circuit_func(sigmap, state_init_angle):
    propagation_circuits = []
    global_angle_for_ry = state_init_angle
    
    for i in range(len(sigmap)):
        angles = [np.angle(sigmap[i][j])*(-2) for j in range(4)]
        qc = QuantumCircuit(3)
        phi1 = (angles[0] + angles[1])/2
        phi2 = (angles[0] - angles[1])/2

        psi1 = (angles[2] + angles[3])/2
        psi2 = (angles[2] - angles[3])/2
        qc.h(0)
        qc.ry(global_angle_for_ry, 1)
    
        qc.rz((phi1 + psi1)/2, 0)
        qc.cx(2, 0)

        qc.rz((phi1 - psi1)/2, 0)
        qc.cx(1, 0)
        qc.cx(2, 0)

        qc.rz((phi2 + psi2)/2, 0)
        qc.cx(2, 0)

        qc.rz((phi2 - psi2)/2, 0)
        qc.cx(1, 0)
        qc.cx(2, 0)
        qc.h(0)
        propagation_circuits.append(qc)
    return propagation_circuits

from qiskit.circuit.library import UnitaryGate    
def magnetization_circuit_func(sigmap, T, state_init_angle):
    U_svd1 = np.linalg.svd(T)[0]
    D_svd1 = np.linalg.svd(T)[1]
    V_svd1 = np.linalg.svd(T)[2]

    dilp_svd1 = []
    diln_svd1 = []
    for i in range(len(D_svd1)):
        d = np.diag(D_svd1).astype(complex)[i, i]/max(D_svd1)
        p = d+1j*np.sqrt((1-np.abs(d)**2+0j)/np.abs(d)**2)*d
        n = d-1j*np.sqrt((1-np.abs(d)**2+0j)/np.abs(d)**2)*d
        dilp_svd1.append(p)
        diln_svd1.append(n)

    angles_svd1 = [np.angle(dilp_svd1[i])*(-2) for i in range(4)]
    phi1_svd1 = (angles_svd1[0] + angles_svd1[1])/2
    phi2_svd1 = (angles_svd1[0] - angles_svd1[1])/2

    psi1_svd1 = (angles_svd1[2] + angles_svd1[3])/2
    psi2_svd1 = (angles_svd1[2] - angles_svd1[3])/2
    
    global_angle_for_ry = state_init_angle
    
    circuits = []
    for i in range(len(sigmap)):
        qc = QuantumCircuit(4)
        angles = [np.angle(sigmap[i][j])*(-2) for j in range(4)]
        phi1 = (angles[0] + angles[1])/2
        phi2 = (angles[0] - angles[1])/2

        psi1 = (angles[2] + angles[3])/2
        psi2 = (angles[2] - angles[3])/2
        qc.h(0)
        qc.ry(global_angle_for_ry, 1)
    
        qc.rz((phi1 + psi1)/2, 0)
        qc.cx(2, 0)

        qc.rz((phi1 - psi1)/2, 0)
        qc.cx(1, 0)
        qc.cx(2, 0)

        qc.rz((phi2 + psi2)/2, 0)
        qc.cx(2, 0)

        qc.rz((phi2 - psi2)/2, 0)
        qc.cx(1, 0)
        qc.cx(2, 0)
        qc.h(0)

        gate_V1 = UnitaryGate(V_svd1)
        qc.append(gate_V1, [1, 2])

        qc.h(3)

        qc.rz((phi1_svd1 + psi1_svd1)/2, 3)
        qc.cx(1, 3)

        qc.rz((phi2_svd1 + psi2_svd1)/2, 3)
        qc.cx(2, 3)
        qc.cx(1, 3)

        qc.rz((phi1_svd1 - psi1_svd1)/2, 3)
        qc.cx(1, 3)

        qc.rz((phi2_svd1 - psi2_svd1)/2, 3)
        qc.cx(2, 3)
        qc.cx(1, 3)  

        gate_U1 = UnitaryGate(U_svd1)
        qc.append(gate_U1, [1, 2])
        qc.h(3)
        qc.measure_all()
        circuits.append(qc)
        
    return circuits
    
    
def tomography_circuit_func(propagation_circuits):
    zz_circuits = []
    zy_circuits = []
    zx_circuits = []
    yz_circuits = []
    yy_circuits = []
    yx_circuits = []
    xz_circuits = []
    xy_circuits = []
    xx_circuits = []
    zi_circuits = []
    yi_circuits = []
    xi_circuits = []
    iz_circuits = []
    iy_circuits = []
    ix_circuits = []
    
    for i in range(len(propagation_circuits)):
        qczz = propagation_circuits[i].copy()
        qczz.measure_all()
        zz_circuits.append(qczz)
    
        qczy = propagation_circuits[i].copy()
        qczy.sdg(1)
        qczy.h(1)
        qczy.measure_all()
        zy_circuits.append(qczy)
    
        qczx = propagation_circuits[i].copy()
        qczx.h(1)
        qczx.measure_all()
        zx_circuits.append(qczx)
    
        qcyz = propagation_circuits[i].copy()
        qcyz.sdg(2)
        qcyz.h(2)
        qcyz.measure_all()
        yz_circuits.append(qcyz)
    
        qcyy = propagation_circuits[i].copy()
        qcyy.sdg(1)
        qcyy.h(1)
        qcyy.sdg(2)
        qcyy.h(2)
        qcyy.measure_all()
        yy_circuits.append(qcyy)
    
        qcyx = propagation_circuits[i].copy()
        qcyx.h(1)
        qcyx.sdg(2)
        qcyx.h(2)
        qcyx.measure_all()
        yx_circuits.append(qcyx)
    
        qcxz = propagation_circuits[i].copy()
        qcxz.h(2)
        qcxz.measure_all()
        xz_circuits.append(qcxz)
    
        qcxy = propagation_circuits[i].copy()
        qcxy.sdg(1)
        qcxy.h(1)
        qcxy.h(2)
        qcxy.measure_all()
        xy_circuits.append(qcxy)
    
        qcxx = propagation_circuits[i].copy()
        qcxx.h(1)
        qcxx.h(2)
        qcxx.measure_all()
        xx_circuits.append(qcxx)
    
        qczi = propagation_circuits[i].copy()
        qczi.measure_all()
        zi_circuits.append(qczi)
    
        qcyi = propagation_circuits[i].copy()
        qcyi.sdg(2)
        qcyi.h(2)
        qcyi.measure_all()
        yi_circuits.append(qcyi)
    
        qcxi = propagation_circuits[i].copy()
        qcxi.h(2)
        qcxi.measure_all()
        xi_circuits.append(qcxi)
    
        qciz = propagation_circuits[i].copy()
        qciz.measure_all()
        iz_circuits.append(qciz)
    
        qciy = propagation_circuits[i].copy()
        qciy.sdg(1)
        qciy.h(1)
        qciy.measure_all()
        iy_circuits.append(qciy)
    
        qcix = propagation_circuits[i].copy()
        qcix.h(1)
        qcix.measure_all()
        ix_circuits.append(qcix)
        
    return zz_circuits, zy_circuits, zx_circuits, yz_circuits, yy_circuits, yx_circuits, xz_circuits, xy_circuits,  xx_circuits, zi_circuits, xi_circuits, yi_circuits, iz_circuits, iy_circuits, ix_circuits