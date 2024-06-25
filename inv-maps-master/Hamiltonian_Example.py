from g_matrix import pruned_sierpinski_G_matrix, optimized_sierpinski_G_matrix
from utils.EncodingClasses import *

'''
This is code written by Andrew Projansky, Jason Neicase, and Brent Harrison

This is an example file for more extended use of this repository. Specifically, 
this example is to showcase class functions which keep track of encodings, and 
Hamiltonians of interest. 

There are two classes of interest: Encoding, and Encoding_with_Hamiltonian. 
Encoding as a class just bookkeeps everything in the basic use example file. 
With system size and G matrix as input, it creates a object with stores N, G, 
J_inv, majoranas, and X and Z paulis. 

Encoding_with_Hamiltonian takes in inputs N, G, as well as information about 
either a Spin or Fermionic Hamiltonian. It then makes the corresponding other 
Hamiltonian using the mapping. 

For list of class attributes, see EncodingClasses.py
'''

N = 4
G_pruned = pruned_sierpinski_G_matrix(N)
Spin_H = Spin_Hamiltonian(N, [['XXII', 'IXYX'], [2, 1.78]])
Hamiltonian_1 = Encoding_with_Hamiltonians(N, G_pruned, Spin_H = Spin_H)

'''
Prints the fermionic hamiltonian corresponding to the toy 
Spin Hamiltonian defined above
'''

print(Hamiltonian_1.Fermionic_H)


'''
Consistency check to show that you can define either the spin or fermionic 
Hamiltonian as starting points
'''
H = Fermionic_Hamiltonian(N, FullPackage.Fermionic_H)
Hamiltonian_2 = Encoding_with_Hamiltonians(N, G_pruned, Fermionic_H = H)
#sign difference on second coefficient from splitting up X and Z
print(Hamiltonian_2.Spin_H)