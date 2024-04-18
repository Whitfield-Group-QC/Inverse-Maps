"""
For Gus & Tom, an example really
"""

from g_matrix import pruned_sierpinski_G_matrix, optimized_sierpinski_G_matrix
from utils.majoranas_utils import get_majoranas_from_G_matrix
from utils.dirac_utils import paulis_maj_string_from_G_matrix


if __name__ == "__main__":
    N = 17

    G_pruned = pruned_sierpinski_G_matrix(N)
    G_opt = optimized_sierpinski_G_matrix(N)

    majs_pruned = get_majoranas_from_G_matrix(G_pruned, N)
    majs_opt = get_majoranas_from_G_matrix(G_opt, N)

    x_pauli_majs_pruned, z_pauli_majs_pruned = paulis_maj_string_from_G_matrix(G_pruned, N)
    x_pauli_majs_opt, z_pauli_majs_opt = paulis_maj_string_from_G_matrix(G_opt, N)

    print(x_pauli_majs_pruned)