from g_matrix import pruned_sierpinski_G_matrix, optimized_sierpinski_G_matrix, random_enc_G_matrix
from utils.majoranas_utils import get_majoranas_from_G_matrix
from utils.dirac_utils import paulis_maj_string_from_G_matrix
from utils.signutils import sign_check
import numpy as np
                


class Encoding:
    
    def __init__(self, N, G):
        
        self.N = N
        self.G = G
        self.majs = get_majoranas_from_G_matrix(self.G, self.N)
        xp, zp, J_inv = paulis_maj_string_from_G_matrix(self.G, self.N)
        self.xps, self.zps, self.phases = sign_check(J_inv.todense(),
                                                     self.majs, xp, zp, self.N)
        self.J_inv = J_inv.todense()
        
        
class Encoding_with_Hamiltonians(Encoding):
    
    def __init__(self, N, G, Fermionic_H = None, Spin_H = None):
        
        super().__init__(N, G)
        if Fermionic_H == None:
            self.Spin_H = Spin_H.Hamiltonian
            self.Fermionic_H = self.Spin_to_Fermion().Hamiltonian
        elif Spin_H == None:
            self.Fermionic_H = Fermionic_H.Hamiltonian
            self.Spin_H = self.Fermion_to_Spin().Hamiltonian
            
    def Spin_to_Fermion(self):
    
        spinHL = self.Spin_H[0]
        spinCL = self.Spin_H[1]
        Fermion_L = []
        Fermion_C = []
        for i in range(len(spinHL)):
            term = spinHL[i]
            base_l = list(np.zeros(2*self.N))
            phase=1
            for k in range(self.N):
                if term[k] == 1:
                    if term[k+self.N] == 1:
                        ind_l, phase_c = self.multiply(self.J_inv.T[k], self.J_inv.T[k+self.N], self.phases[k]*self.phases[k+self.N])
                    else:
                        ind_l = self.J_inv.T[k]
                        phase_c = self.phases[k]
                    base_l, phase = self.multiply(base_l, ind_l, phase_c*phase)
                elif term[k+self.N] == 1:
                    ind_l = self.J_inv.T[k+self.N]
                    phase_c = self.phases[k+self.N]
                    base_l, phase = self.multiply(base_l, ind_l, phase_c*phase)
            Fermion_L.append(base_l)
            Fermion_C.append(spinCL[i]*phase)
        return Fermionic_Hamiltonian(self.N, [Fermion_L, Fermion_C])
        
    def multiply(self, l1, l2, phase):
        
        parity = 1
        for ind in range(2*self.N):
            if l2[ind] == 1:
                for ind_2 in range(2*self.N-1,ind,-1):
                    if l1[ind_2] == 1:
                        parity = parity * -1
                l1[ind] = (l1[ind]+l2[ind])%2
        return l1, parity*phase
    
    def Pauli_mult(self, p1, p2):
        
        m_dict = {('I', 'I'): ['I', 1], ('X', 'I'): ['X', 1], ('I', 'X'): ['X', 1], 
                  ('Y', 'I'): ['Y', 1], ('I', 'Y'): ['Y', 1], ('Z', 'I'): ['Z', 1], 
                  ('I', 'Z'): ['Z', 1], ('X', 'X'): ['I', 1], ('Y', 'X'): ['Z', -1j], 
                  ('X', 'Y'): ['Z', 1j], ('Z', 'Y'): ['X', -1j], ('Y', 'Z'): ['X', 1j], 
                  ('X', 'Z'): ['Y', 1j], ('Z', 'X'): ['Y', -1j], ('Y', 'Y'): ['I', 1],
                  ('Z', 'Z'): ['I', 1]}
        
        new_word = ''
        new_phase = 1
        for lett in range(len(p1)):
            key = (p1[lett], p2[lett])
            new_word += m_dict[key][0]
            new_phase = new_phase * np.conj(m_dict[key][1])
        return new_word, new_phase
        
        
    def Fermion_to_Spin(self):
        
        FermiHL = self.Fermionic_H[0]
        FermiCL = self.Fermionic_H[1]
        Spin_L = []
        Spin_C = []
        for i in range(len(FermiHL)):
            term= FermiHL[i]
            base_str = 'I'*self.N
            phase = 1
            for coeff in range(len(term)):
                if term[coeff] == 1:
                    base_str, p_phase = self.Pauli_mult(base_str, self.majoranas[coeff])
                    phase=p_phase*phase
            Symp_list = list(np.zeros(2*self.N))
            for lett in range(self.N):
                if base_str[lett] == 'X':
                    Symp_list[lett] = 1
                elif base_str[lett] == 'Y':
                    Symp_list[lett] = 1
                    Symp_list[lett+self.N] = 1
                elif base_str[lett] == 'Z':
                    Symp_list[lett+self.N] = 1
            Spin_L.append(Symp_list)
            Spin_C.append(FermiCL[i]*phase)
        return Spin_Hamiltonian(self.N, [Spin_L, Spin_C])
        
class Spin_Hamiltonian:
    
    def __init__(self, N, Hamiltonian):
        
        self.N = N
        if type(Hamiltonian[0][0] == type('II')):
            self.Hamiltonian = self.Make_Sympletic(Hamiltonian)
            self.Written_H = Hamiltonian
        elif type(Hamiltonian[0][0] == list):
            self.Hamiltonian = Hamiltonian
            self.Written_H = self.Make_Written
        
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
            Symp_L.append(list(Symp_l))
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
        self.Hamiltonian = Hamiltonian
        
    def print_strs(self):
        
        print_l = []
        for term in self.Hamiltonian[0]:
            print_str = ''
            for maj in range(self.N):
                if term[maj] == 1:
                    print_str = print_str + "c'_" + str(maj) +  ' '
            print_l.append(print_str)
        print([print_l, self.Hamiltonian[1]])
    