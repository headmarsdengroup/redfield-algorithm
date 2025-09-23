import numpy as np
import scipy
import cvxpy as cp

sigz = np.matrix([[1, 0],[0, -1]]) 
sigy = np.matrix([[0, -1j], [1j, 0]])
sigx = np.matrix([[0, 1], [1, 0]])
iden = np.matrix([[1, 0], [0, 1]])

def fidel(rho, sigma):
    x = scipy.linalg.sqrtm(rho)
    y = scipy.linalg.sqrtm(x@sigma@x)
    return y.trace()**2

n_qubit = 2
n_shots = 2**10

def tomo_full(counts_list):
    for i in range(len(counts_list)):
        for j in ['000', '010', '100', '110']:
            if j not in list(counts_list[i].keys()):
                counts_list[i][j] = 0
    zz = (counts_list[0]['110']+counts_list[0]['000']-counts_list[0]['010']-counts_list[0]['100'])
    zz /= (counts_list[0]['110']+counts_list[0]['000']+counts_list[0]['010']+counts_list[0]['100'])
    mat_zz = zz*np.kron(sigz, sigz)
    
    zy = (counts_list[1]['110']+counts_list[1]['000']-counts_list[1]['010']-counts_list[1]['100'])
    zy /= (counts_list[1]['110']+counts_list[1]['000']+counts_list[1]['010']+counts_list[1]['100'])
    mat_zy = zy*np.kron(sigz, sigy)
    
    zx = (counts_list[2]['110']+counts_list[2]['000']-counts_list[2]['010']-counts_list[2]['100'])
    zx /= (counts_list[2]['110']+counts_list[2]['000']+counts_list[2]['010']+counts_list[2]['100'])
    mat_zx = zx*np.kron(sigz, sigx)
    
    yz = (counts_list[3]['110']+counts_list[3]['000']-counts_list[3]['010']-counts_list[3]['100'])
    yz /= (counts_list[3]['110']+counts_list[3]['000']+counts_list[3]['010']+counts_list[3]['100'])
    mat_yz = yz*np.kron(sigy, sigz)
    
    yy = (counts_list[4]['110']+counts_list[4]['000']-counts_list[4]['010']-counts_list[4]['100'])
    yy /= (counts_list[4]['110']+counts_list[4]['000']+counts_list[4]['010']+counts_list[4]['100'])
    mat_yy = yy*np.kron(sigy, sigy)
    
    yx = (counts_list[5]['110']+counts_list[5]['000']-counts_list[5]['010']-counts_list[5]['100'])
    yx /= (counts_list[5]['110']+counts_list[5]['000']+counts_list[5]['010']+counts_list[5]['100'])
    mat_yx = yx*np.kron(sigy, sigx)
    
    xz = (counts_list[6]['110']+counts_list[6]['000']-counts_list[6]['010']-counts_list[6]['100'])
    xz /= (counts_list[6]['110']+counts_list[6]['000']+counts_list[6]['010']+counts_list[6]['100'])
    mat_xz = xz*np.kron(sigx, sigz)
    
    xy = (counts_list[7]['110']+counts_list[7]['000']-counts_list[7]['010']-counts_list[7]['100'])
    xy /= (counts_list[7]['110']+counts_list[7]['000']+counts_list[7]['010']+counts_list[7]['100'])
    mat_xy = xy*np.kron(sigx, sigy)
    
    xx = (counts_list[8]['110']+counts_list[8]['000']-counts_list[8]['010']-counts_list[8]['100'])
    xx /= (counts_list[8]['110']+counts_list[8]['000']+counts_list[8]['010']+counts_list[8]['100'])
    mat_xx = xx*np.kron(sigx, sigx)
    
    zi = (counts_list[9]['010']+counts_list[9]['000']-counts_list[9]['110']-counts_list[9]['100'])
    zi /= (counts_list[9]['010']+counts_list[9]['000']+counts_list[9]['110']+counts_list[9]['100'])
    mat_zi = zi*np.kron(sigz, iden)
    
    yi = (counts_list[10]['010']+counts_list[10]['000']-counts_list[10]['110']-counts_list[10]['100'])
    yi /= (counts_list[10]['010']+counts_list[10]['000']+counts_list[10]['110']+counts_list[10]['100'])
    mat_yi = yi*np.kron(sigy, iden)
    
    xi = (counts_list[11]['010']+counts_list[11]['000']-counts_list[11]['110']-counts_list[11]['100'])
    xi /= (counts_list[11]['010']+counts_list[11]['000']+counts_list[11]['110']+counts_list[11]['100'])
    mat_xi = xi*np.kron(sigx, iden)
    
    iz = (counts_list[12]['100']+counts_list[12]['000']-counts_list[12]['010']-counts_list[12]['110'])
    iz /= (counts_list[12]['100']+counts_list[12]['000']+counts_list[12]['010']+counts_list[12]['110'])
    mat_iz = iz*np.kron(iden, sigz)
    
    iy = (counts_list[13]['100']+counts_list[13]['000']-counts_list[13]['010']-counts_list[13]['110'])
    iy /= (counts_list[13]['100']+counts_list[13]['000']+counts_list[13]['010']+counts_list[13]['110'])
    mat_iy = iy*np.kron(iden, sigy)
    
    ix = (counts_list[14]['100']+counts_list[14]['000']-counts_list[14]['010']-counts_list[14]['110'])
    ix /= (counts_list[14]['100']+counts_list[14]['000']+counts_list[14]['010']+counts_list[14]['110'])
    mat_ix = ix*np.kron(iden, sigx)
    
    mat = mat_zz + mat_zy + mat_zx + mat_yz + mat_yy + mat_yx + mat_xz + mat_xy +mat_xx + mat_zi + mat_yi + mat_xi + mat_iz + mat_iy + mat_ix
    mat += np.kron(iden, iden)
    mat /= n_qubit**2
    mat /= np.trace(mat)  
    return mat

e1 = np.array([[0, 1]])
e0 = np.array([[1, 0]])
def mag_expval_from_bitstring(counts_list):
    for j in ['0000', '0001', '0010', '0011', '0100', '0101', '0110', '0111', 
                  '1000', '1001', '1010', '1011', '1100', '1101', '1110', '1111']:
        if j not in list(counts_list.keys()):
            counts_list[j] = 0
                
    zzzz = (counts_list['0000']/n_shots)*np.kron(e0, np.kron(e0, np.kron(e0, e0)))
    zzzo = (counts_list['0001']/n_shots)*np.kron(e0, np.kron(e0, np.kron(e0, e1)))
    zzoz = (counts_list['0010']/n_shots)*np.kron(e0, np.kron(e0, np.kron(e1, e0)))
    zzoo = (counts_list['0011']/n_shots)*np.kron(e0, np.kron(e0, np.kron(e1, e1)))
    
    zozz = (counts_list['0100']/n_shots)*np.kron(e0, np.kron(e1, np.kron(e0, e0)))
    zozo = (counts_list['0101']/n_shots)*np.kron(e0, np.kron(e1, np.kron(e0, e1)))
    zooz = (counts_list['0110']/n_shots)*np.kron(e0, np.kron(e1, np.kron(e1, e0)))
    zooo = (counts_list['0111']/n_shots)*np.kron(e0, np.kron(e1, np.kron(e1, e1)))
    
    ozzz = (counts_list['1000']/n_shots)*np.kron(e1, np.kron(e0, np.kron(e0, e0)))
    ozzo = (counts_list['1001']/n_shots)*np.kron(e1, np.kron(e0, np.kron(e0, e1)))
    ozoz = (counts_list['1010']/n_shots)*np.kron(e1, np.kron(e0, np.kron(e1, e0)))
    ozoo = (counts_list['1011']/n_shots)*np.kron(e1, np.kron(e0, np.kron(e1, e1)))
    
    oozz = (counts_list['1100']/n_shots)*np.kron(e1, np.kron(e1, np.kron(e0, e0)))
    oozo = (counts_list['1101']/n_shots)*np.kron(e1, np.kron(e1, np.kron(e0, e1)))
    oooz = (counts_list['1110']/n_shots)*np.kron(e1, np.kron(e1, np.kron(e1, e0)))
    oooo = (counts_list['1111']/n_shots)*np.kron(e1, np.kron(e1, np.kron(e1, e1)))
    
    sv = zzzz + zzzo + zzoz + zzoo + zozz + zozo + zooz + zooo + ozzz + ozzo + ozoz + ozoo + oozz + oozo + oooz + oooo
    return sv

def postprocessor(counts_list, T):
    post_list = []
    for i in range(len(counts_list)):
        tomo = tomo_full(counts_list[i])
        
        rho00 = np.sqrt(tomo[0, 0])
        rho11 = np.sqrt(tomo[3, 3])
        rho01 = (tomo[0, 1]/rho00).conjugate()
        rho10 = (tomo[0, 2]/rho00).conjugate()
        mat = np.array([[rho00, rho01], [rho10, rho11]])
        
        back = mat.reshape(-1, 1)
        back = T@(back.reshape(-1, 1))
        back /= np.linalg.norm(back)
        
        back = back.reshape(2, 2)
        back /= np.trace(back)
        
        X = cp.Variable((back.shape), hermitian=True)
        constraints = [cp.constraints.psd.PSD(X)]
        constraints +=[cp.trace(X)==1]
        prob = cp.Problem(cp.Minimize(cp.norm(X-back, 'fro')), constraints)
        prob.solve(solver=cp.SCS)
        post_list.append(X.value)
        
    return post_list

def upper_block_dilation(R_diag, t_list):
    expR_diag_list = []
    dilp_list = []
    diln_list = []
    
    for i in range(len(t_list)):
        expR_diag_list.append(scipy.linalg.expm(R_diag*t_list[i]))
        
    for i in range(len(expR_diag_list)):
        expR_diag_list[i]/= max(np.abs(np.diag(expR_diag_list[i])))
        
    for i in range(len(expR_diag_list)):
        dilp = []
        diln = []
        for j in range(len(expR_diag_list[i])):
            d = expR_diag_list[i][j, j]
            p = d+1j*np.sqrt((1-np.abs(d)**2)/np.abs(d)**2)*d
            n = d-1j*np.sqrt((1-np.abs(d)**2)/np.abs(d)**2)*d
            dilp.append(p)
            diln.append(n)
        dilp_list.append(dilp)
        diln_list.append(diln)
        
    return dilp_list

def magnetization_expval(counts_list):
    mag_list = []
    for i in range(len(counts_list)):
        sv = mag_expval_from_bitstring(counts_list[i])
        a = np.array([sv[0][j] for j in range(8) if j%2==0])
        b = np.sqrt(a)
        dm = b.reshape(2, 2)
        dm /= np.trace(dm)
        mag = np.trace(dm@sigz) 
        mag_list.append(mag)
    return mag_list

def U3(theta, phi, lam):
    return np.array([[np.cos(theta/2), -np.exp(1j*lam)*np.sin(theta/2)], [np.exp(1j*phi)*np.sin(theta/2), np.exp(1j*(phi + lam))*np.cos(theta/2)]])

def opt(params, back_gs_normalized):
    u1 = U3(params[0], params[1], params[2])
    u2 = U3(params[3], params[4], params[5])
    target = np.kron(u1, u2)@np.array([[1, 0, 0, 0]]).reshape(-1, 1) 
    return 1-np.linalg.norm(target.conjugate().T@back_gs_normalized)

def tlist(tmin, tmax, n_t):
    pow = 1.2
    a = tmax/(pow**(n_t-1)-1)
    b = -a
    t_list = np.zeros((n_t))
    for i in range(n_t):
        t_list[i]= a*pow**(i)+b
    return t_list