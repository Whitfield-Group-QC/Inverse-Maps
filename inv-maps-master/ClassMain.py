from g_matrix import pruned_sierpinski_G_matrix, optimized_sierpinski_G_matrix
from utils.EncodingClasses import *

N = 4
G_pruned = pruned_sierpinski_G_matrix(N)
Enc = Encoding(N, G_pruned)
#Spin_H = Spin_Hamiltonian(N, [['XXII', 'IXYX'], [2, 1.78]])
#my_encoding = Encoding(N, G_pruned)
#FullPackage = Encoding_with_Hamiltonians(N, G_pruned, Spin_H = Spin_H)

#H = Fermionic_Hamiltonian(N, FullPackage.Fermionic_H)
#FullPackage2 = Encoding_with_Hamiltonians(N, G_pruned, Fermionic_H = H)

###There is an error with signs I have to track down... but multiplication seems to work