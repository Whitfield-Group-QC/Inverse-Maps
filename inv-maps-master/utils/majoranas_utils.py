import scipy as sp
import numpy as np

from sympy import SparseMatrix, mod_inverse
from sympy.matrices.sparsetools import _doktocsr
from scipy.sparse import hstack, vstack, identity


def efficient_mod_inverse(G, N, modulus=2, dtype=int):
    to_rref = sp.sparse.hstack([G, sp.sparse.identity(N)])
    A = SparseMatrix(N, 2*N, to_rref.todok())
    A_rref = A.rref(iszerofunc=lambda x: x % modulus==0)
    inv = A_rref[0][:,N:]
    data = _doktocsr(inv)
    vals = [int(_mod_custom(x, modulus)) for x in data[0]]
    return sp.sparse.csr_array((vals, data[1], data[2]), shape=data[3], dtype=dtype)


def _mod_custom(x, modulus):
    numer, denom = x.as_numer_denom()
    return numer*mod_inverse(denom,modulus) % modulus


def BSM_matrix_from_G_matrix(G, N):
    U = sp.sparse.lil_array((N,N))
    for i in range(N):
        for j in range(i+1,N):
            U[i,j] = 1

    T = U + identity(N)

    G_invT = efficient_mod_inverse(G, N).transpose()

    # The binary simplectic form in tableu order - rows are majoranas
    binary_symp_majs = vstack([hstack([G, G], dtype=int), hstack([G_invT @ U, G_invT @ T], dtype=int)], dtype=int, format='csr').transpose()
    binary_symp_majs.data %= 2

    return binary_symp_majs


def avg_pauli_weight_gmatrix(G, N):
    
    maj_arr = BSM_arrs(G, N)
    x_arr = maj_arr[0:N]
    z_arr = maj_arr[N:]

    weight_arr = np.logical_or(x_arr, z_arr)

    total_parity_x = np.sum(x_arr, axis=0) % 2
    total_parity_z = np.sum(z_arr, axis=0) % 2

    weight_parity = np.sum(np.logical_or(total_parity_x, total_parity_z))
    

    return (np.sum(weight_arr, axis=(0,1)) + weight_parity) / (2*N + 1)


def BSM_arrs(G, N):
    binary_symp_majs = BSM_matrix_from_G_matrix(G, N)

    x_arr = np.array(binary_symp_majs[:,:N].todense())
    z_arr = np.array(binary_symp_majs[:,N:].todense())
    maj_arr = []
    for i in range(2*N):
        maj_arr.append(np.append(x_arr[i], z_arr[i]))
    return maj_arr


def get_majoranas_from_G_matrix(G, N):

    maj_arr = BSM_arrs(G, N)

    majs = []
    ops = {(0,0): 'I', (1,0): 'X', (0,1): 'Z', (1,1): 'Y'}

    for m in range(2*N):
        maj = ''
        for i in range(N):
            maj += ops[(maj_arr[m][i], maj_arr[m][i+N])]
        majs.append(maj)
    return majs

def majs_to_Pauli(maj_arr, N, Majs):
    
    '''
    Takes input of the majoranas of an encoding, outputs a single Pauli
    
    example input Majs: c_0 c_4' [1,0,0,0,0,0,0,1]: corresponds to c_1 c_1'
    '''
    
    ops = {(0,0): 'I', (1,0): 'X', (0,1): 'Z', (1,1): 'Y'}
    if len(Majs) != 2*N:
        return 'Error: Incorrect vector size'
    else:
        arr = np.array([0]*2*N)
        use_maj = np.array(maj_arr)
        for p in range(len(Majs)):
            if Majs[p] == 1:
                arr = (use_maj[p] + arr) % 2
        maj = ''
        for i in range(N):
            maj += ops[(arr[i], arr[i+N])]
        return maj
            
    
    