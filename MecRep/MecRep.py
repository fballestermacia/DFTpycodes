import numpy as np
import spglib as spg
import ase.io as io
import ase as ase
import ase.spacegroup as asespg
from ase.visualize import view
from spgrep.irreps import decompose_representation


def chech_permu(permu):
    dim=len(permu[0])
    for i in range(dim):
        sum=0
        for j in range(dim):
            sum=sum+permu[i,j]
        if sum!=1:
            return False
    return True

   

if __name__ == '__main__':
    scf_file = r'AgP2\Phonons_V2\AgP2_scf.in'
    structure = io.read(scf_file)
    
    basisvec = structure.get_cell()
    positions = structure.get_scaled_positions()
    atmnum = structure.get_atomic_numbers()
    recbasis = structure.cell.reciprocal()
    
    #print(recbasis)
    
    '''
    group = asespg.get_spacegroup(structure)
    print(group.symmetry_normalised_sites(structure.get_scaled_positions()))
    print(group.no)
    print(group.get_op())
    print(group.rotations)
    print(group.translations)
    
    print(group.equivalent_sites([0,0.1,0.5]))'''
    
    sym_data = spg.get_symmetry_dataset([basisvec,positions,atmnum])
    rotations = sym_data.get("rotations")
    translations = sym_data.get("translations")
    print(sym_data.get("wyckoffs"))
    
    
    i = 1
    sym_prec=4
    kpoint = np.array([0,0,0])
    num_atoms = len(atmnum)
    permu = np.zeros([num_atoms,num_atoms], dtype='complex')
    cells = np.zeros([num_atoms,3])
    
    '''phase = np.exp(1j*np.matmul(np.transpose(np.matmul(rotations[i],kpoint)), translations[i]))
    phase2 = np.exp(1j*np.matmul(np.transpose(kpoint), translations[i]))
    print(phase)
    print(phase2)'''
    
    for j in range(num_atoms):
        new_pos=np.matmul(rotations[i],positions[j]) + translations[i]
        #print(positions[j],new_pos)
        
        #Store phases
        
        #Move everything back to Unit Cell
        if np.any(new_pos<0) or np.any(new_pos>=1):
            for k in range(3):
                while new_pos[k]<0:
                    new_pos[k]=new_pos[k]+1
                    cells[j,k]=cells[j,k]-1
                    
                while new_pos[k]>=1:
                    new_pos[k]=new_pos[k]-1
                    cells[j,k]=cells[j,k]+1
          
        #print(positions[j],new_pos)
       
        #Find permutation
        for k in range(num_atoms):
            if np.array_equal(np.round(new_pos,decimals=sym_prec),np.round(positions[k],decimals=sym_prec)):
                permu[k,j]=1*np.round(np.exp(2j*np.pi*np.matmul(np.transpose(kpoint), cells[k])),3)
                
    if chech_permu(permu):
        print(np.matmul(np.matmul(np.transpose(basisvec),rotations[i]),np.linalg.inv(np.transpose(basisvec))))
        #Check if maybe the basisvec should be transposed
        
        print(cells)
        print(permu)


    

