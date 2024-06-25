import numpy as np

def multiply(majc, majp):
    
    m_dict = {('I', 'I'): ['I', 1], ('X', 'I'): ['X', 1], ('I', 'X'): ['X', 1], 
              ('Y', 'I'): ['Y', 1], ('I', 'Y'): ['Y', 1], ('Z', 'I'): ['Z', 1], 
              ('I', 'Z'): ['Z', 1], ('X', 'X'): ['I', 1], ('Y', 'X'): ['Z', -1j], 
              ('X', 'Y'): ['Z', 1j], ('Z', 'Y'): ['X', -1j], ('Y', 'Z'): ['X', 1j], 
              ('X', 'Z'): ['Y', 1j], ('Z', 'X'): ['Y', -1j], ('Y', 'Y'): ['I', 1],
              ('Z', 'Z'): ['I', 1]}
    
    majc_str = majc[0]
    new_word = ''
    new_phase = majc[1]
    for lett in range(len(majc_str)):
        key = (majc_str[lett], majp[lett])
        #print(majc[0][lett], m_dict[key][0])
        new_word += m_dict[key][0]
        #print(majc[1])
        new_phase = new_phase * np.conj(m_dict[key][1])
        #print(new_phase)
    majc = (new_word, new_phase)
    return majc
        
    
def sign_check(J_inv, majs, xps, zps, N):
    
    phase_v = []
    majf = majs[:N]
    maje = majs[N:]
    majs = [val for pair in zip(majf, maje) for val in pair]
    #x-check_first 
    for row in range(len(J_inv.T)//2):
        majc = ('I'*N, 1)
        new_row = J_inv.T[row]
        new_rowup = new_row[:N]
        new_rowp = new_row[N:]
        new_rowf = [val for pair in zip(new_rowup, new_rowp) for val in pair]
        for maj in range(len(new_rowf)):
            if new_rowf[maj] == 1:
                majc = multiply(majc, majs[maj])
        xps[row] = str(majc[1]) + ' ' + xps[row]
        phase_v.append(majc[1])
    
    for row in range(len(J_inv.T)//2):
        majc = ('I'*N, 1)
        new_row = J_inv.T[row+N]
        new_rowup = new_row[:N]
        new_rowp = new_row[N:]
        new_rowf = [val for pair in zip(new_rowup, new_rowp) for val in pair]
        for maj in range(len(new_rowf)):
            if new_rowf[maj] == 1:
                majc = multiply(majc, majs[maj])
        zps[row] = str(majc[1]) + ' '  + zps[row]
        phase_v.append(majc[1])
        
    return xps, zps, phase_v

