"""
For Gus & Tom, an example really
"""

from g_matrix import pruned_sierpinski_G_matrix, optimized_sierpinski_G_matrix
from utils.majoranas_utils import get_majoranas_from_G_matrix
from utils.dirac_utils import paulis_maj_string_from_G_matrix
from utils.signutils import sign_check
                

if __name__ == "__main__":
    N = 4

    G_pruned = pruned_sierpinski_G_matrix(N)
    G_opt = optimized_sierpinski_G_matrix(N)

    majs_pruned = get_majoranas_from_G_matrix(G_pruned, N)
    majs_opt = get_majoranas_from_G_matrix(G_opt, N)

    x_pauli_majs_pruned, z_pauli_majs_pruned, J_inv = paulis_maj_string_from_G_matrix(G_pruned, N)
    x_pauli_majs_opt, z_pauli_majs_opt, J_inv2 = paulis_maj_string_from_G_matrix(G_opt, N)

    xps, zps = sign_check(J_inv.todense(), majs_pruned, x_pauli_majs_pruned, z_pauli_majs_pruned, N)
    print(xps, majs_pruned)

                
    