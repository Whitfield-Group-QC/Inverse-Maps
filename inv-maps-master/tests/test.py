from scipy.sparse import identity

from ..utils.majoranas_utils import get_majoranas_from_G_matrix
from ..utils.dirac_utils import J_inv_from_G_matrix, paulis_maj_string_from_J_inv


def test_JW_majs_from_paulis(N):

    _majs_exact = []
    for i in range(N):
        _majs_exact.append('Z'*i + 'X' + 'I'*(N-i-1))
        _majs_exact.append('Z'*i + 'Y' + 'I'*(N-i-1))
    majs_exact = set(_majs_exact)

    G = identity(N)
    majs = set(get_majoranas_from_G_matrix(G, N))
    
    assert majs == majs_exact
    return(majs)


def test_JW_paulis_from_majs(N):
    
    x_paulis_exact, z_paulis_exact = [], []

    for p in range(N):
        x_p_exact, z_p_exact = '', f'c_{p} c\'_{p} '
        for i in range(p):
            x_p_exact += f'c_{i} c\'_{i} '
        x_p_exact += f'c_{p} '

        x_paulis_exact.append(x_p_exact)
        z_paulis_exact.append(z_p_exact)

    J_inv = J_inv_from_G_matrix(identity(N), N)
    x_strings, z_strings = paulis_maj_string_from_J_inv(J_inv, N)

    assert x_paulis_exact == x_strings
    assert z_paulis_exact == z_strings

    return(x_strings, z_strings)


if __name__ == "__main__":
    N = 7

    print(test_JW_majs_from_paulis(N))
    print('')
    print(test_JW_paulis_from_majs(N))