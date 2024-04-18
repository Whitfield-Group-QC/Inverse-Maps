"""
For Gus & Tom, an example really
"""

from g_matrix import pruned_sierpinski_G_matrix, optimized_sierpinski_G_matrix
from utils.majoranas_utils import get_majoranas_from_G_matrix
from utils.dirac_utils import paulis_maj_string_from_G_matrix
from utils.majoranas_utils import majs_to_Pauli
from utils.majoranas_utils import BSM_arrs
from utils.dirac_utils import p_2_mat


if __name__ == "__main__":
    N = 5

    G_pruned = pruned_sierpinski_G_matrix(N)
    G_opt = optimized_sierpinski_G_matrix(N)

    majs_pruned = get_majoranas_from_G_matrix(G_pruned, N)
    majs_opt = get_majoranas_from_G_matrix(G_opt, N)

    x_pauli_majs_pruned, z_pauli_majs_pruned, J_inv_p = paulis_maj_string_from_G_matrix(G_pruned, N)
    x_pauli_majs_opt, z_pauli_majs_opt, J_inv_o = paulis_maj_string_from_G_matrix(G_opt, N)
    
    maj_mat = BSM_arrs(G_opt, N)
    p = majs_to_Pauli(maj_mat, N, [1,0,0,1,0,0,0,0,0,0])
    
    maj = p_2_mat(J_inv_o, N, 'IYIII')
    
    #print(x_pauli_majs_pruned)