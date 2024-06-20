from g_matrix import pruned_sierpinski_G_matrix, optimized_sierpinski_G_matrix, random_enc_G_matrix
from utils.majoranas_utils import get_majoranas_from_G_matrix
from utils.dirac_utils import paulis_maj_string_from_G_matrix
from utils.signutils import sign_check
import numpy as np
                


class Encoding:
    
    def __init__(self, N, G):
        
        self.N = N
        self.G = G
        self.majoranas = get_majoranas_from_G_matrix(self.G, self.N)
        xp, zp, self.J_inv = paulis_maj_string_from_G_matrix(self.G, self.N)
        self.xps, self.zps, self.phases = sign_check(self.J_inv.todense(),
                                                     self.majs, xp, zp, self.N)
        
        
class Encoding_with_Hamiltonians:
    
    def __init__(self, N, G, Fermionic_H = None, Spin_H = None):
        
        super().__init__(N, G)
        if Fermionic_H == None:
            self.Spin_H = Spin_H
            self.Fermionic_H = self.Spin_to_Fermion(self)
        if Spin_H == None:
            self.Fermionic_H = Fermionic_H
            self.Spin_H = self.Fermion_to_Spin(self)
            
    def Spin_to_Fermion():
        pass
        '''
        SpinHL = self.Spin_H.Paulis
        SpinCL = self.Spin_H.Coeffs
        for val in range(len(SpinHL)):
            multiply out
        '''
        
    def Fermion_to_Spin():
        pass
        
class Spin_Hamiltonian:
    
    def __init__(self, N, Hamiltonian):
        
        self.N = N
        if type(Hamiltonian[0][0] == str):
            self.Hamiltonian = self.Make_Sympletic(Hamiltonian)
            self.Written_H = Hamiltonian
        if type(Hamiltonian[0][0] == list):
            self.Hamiltonian = Hamiltonian
            self.Written_H = self.Make_written
        
    def Make_Sympletic(self, Hamiltonian):
        
        Symp_L = []
        phase_L = Hamiltonian[1]
        for st_pauli in range(len(Hamiltonian[0])):
            p = Hamiltonian[0][st_pauli]
            Symp_l = np.zeros(2*self.N)
            for let in range(self.N):
                if p[let] == 'X':
                    Symp_l[let] = 1
                elif p[let] == 'Z':
                    Symp_l[let+self.N] = 1
                elif p[let] == 'Y':
                    Symp_l[let] = 1
                    Symp_l[let+self.N] = 1
                    phase_L[st_pauli] = phase_L[st_pauli] * 1j
            Symp_L.append(Symp_l)
        return [Symp_L, phase_L]
    
    def Make_Written(self, Hamiltonian):
        
        String_List = []
        phase_L = Hamiltonian[1]
        for symp_pauli in range(len(Hamiltonian[0])):
            p = Hamiltonian[0][symp_pauli]
            stp = ''
            for let in range(self.N):
                if p[let] == 1:
                    if p[let+self.N] == 1:
                        stp = stp + 'Y'
                        phase_L[symp_pauli] = phase_L[symp_pauli] / 1j
                    else:
                        stp = stp + 'X'
                elif p[let+self.N] == 1:
                    stp = stp + 'Z'
                else:
                    stp = stp + 'I'
            String_List.append(stp)
        return [String_List, phase_L]
    
class Fermionic_Hamiltonian:
    
    def __init__(self, N, Hamiltonian):
        self.N = N
        self.vec_Hamiltonian = Hamiltonian
        
    def print_strs(self):
        
        print_l = []
        for term in self.vec_Hamiltonian[0]:
            print_str = ''
            for maj in range(self.N):
                if term[maj] == 1:
                    print_str = print_str + "c'_" + str(maj) +  ' '
            print_l.append(print_str)
        print([print_l, self.vec_Hamiltonian[1]])