from .majoranas_utils import BSM_matrix_from_G_matrix, efficient_mod_inverse
import numpy as np


def J_from_G_matrix(G, N):
    return BSM_matrix_from_G_matrix(G, N).transpose()


def J_inv_from_G_matrix(G, N):
    return efficient_mod_inverse(J_from_G_matrix(G, N), 2*N)


def paulis_maj_string_from_J_inv(J_inv, N):
    
    x_paulis, z_paulis = [], []
    
    for p in range(N):
        x_p, z_p = '', ''
        for i in range(N):
            if J_inv[i, p] == 1:
                x_p += f'c_{i} '
            if J_inv[i+N, p] == 1:
                x_p += f'c\'_{i} '

            if J_inv[i, p+N] == 1:
                z_p += f'c_{i} '
            if J_inv[i+N, p+N] == 1:
                z_p += f'c\'_{i} '

        x_paulis.append(x_p)
        z_paulis.append(z_p)
    
    return x_paulis, z_paulis, J_inv


def paulis_maj_string_from_G_matrix(G, N):
    J_inv = J_inv_from_G_matrix(G, N)
    return paulis_maj_string_from_J_inv(J_inv, N)

def p_2_mat(J_inv, N, pstr):
    
    ops_r = {'I': (0,0), 'X': (1,0), 'Z': (0,1), 'Y': (1,1)}
    p_list = np.array([0]*2*N)
    for let in range(len(pstr)):
        p_list[let] = ops_r[pstr[let]][0]
        p_list[let+N] = ops_r[pstr[let]][1]
    J_inv_t = J_inv.T
    maj_l = np.array([0]*2*N)
    for p in range(2*N):
        if p_list[p] == 1:
            maj_l = (maj_l + J_inv_t[[p], :]) % 2
    return maj_l