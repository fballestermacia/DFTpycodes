import numpy as np
import spglib as spg
import ase.io as io
import ase as ase
import ase.spacegroup as asespg
from spgrep import get_spacegroup_irreps
from spgrep.representation import get_character, project_to_irrep,check_spacegroup_representation
from spgrep.irreps import decompose_representation
from utilsQE import readModesatKpoin


def get_displacements_representation(
    lattice,
    positions,
    little_rotations,
    little_translations,
    qpoint,
):
    r"""Compute representation matrix for fourier-transformed displacements.

    .. math::
       \\Gamma_{\\kappa'\\mu'; \\kappa\\mu}^{\\mathbf{q}}(g) := \\exp \\left( -i \\mathbf{R}_{g} \\mathbf{q} \\cdot \\mathbf{h}_{g}(\\kappa) \\right) [\\mathbf{R}_{g}]_{\\mu'\\mu} \\delta_{ g\\kappa, \\kappa' }
    """
    little_order = len(little_rotations)
    num_atoms = len(positions)

    # Operation-`i` moves atom-`kappa` to `permutations[i, kappa]`
    permutations = np.zeros((little_order, num_atoms), dtype=int)
    for i, (Ri, vi) in enumerate(zip(little_rotations, little_translations)):
        for kappa, position in enumerate(positions):
            new_pos = np.remainder(Ri @ position + vi, 1)
            for kappa2, position2 in enumerate(positions):
                if np.allclose(position2, new_pos):
                    permutations[i, kappa] = kappa2
                    break

    shifts = np.zeros((little_order, num_atoms, 3))
    for i, (Ri, vi) in enumerate(zip(little_rotations, little_translations)):
        perm_i = permutations[i]
        shifts[i] = positions @ Ri.T + vi[None, :] - positions[perm_i]

    perm_rep = np.zeros((little_order, num_atoms, num_atoms), dtype=np.complex_)
    for i, Ri in enumerate(little_rotations):
        for kappa in range(num_atoms):
            kappa2 = permutations[i, kappa]
            perm_rep[i, kappa2, kappa] = np.exp(
                -2j * np.pi * np.dot(Ri.T @ qpoint, shifts[i, kappa])
            )

    # Rotation matrix in cartesian (order, 3, 3)
    A = np.transpose(lattice)  # column-wise lattice vectors
    Ainv = np.linalg.inv(A)
    rotation_rep = np.array([A @ r @ Ainv for r in little_rotations], dtype=np.complex_)

    rep = np.einsum("ipq,iab->ipaqb", perm_rep, rotation_rep, optimize="greedy")
    return rep.reshape(-1, num_atoms * 3, num_atoms * 3)


def permutation_matrix(pos,r,t,q,sym_prec=4):
    kpoint = q
    num_atoms = len(pos)
    permu = np.zeros([num_atoms,num_atoms], dtype='complex')
    cells = np.zeros([num_atoms,3])
    
    '''phase = np.exp(1j*np.matmul(np.transpose(np.matmul(rotations[i],kpoint)), translations[i]))
    phase2 = np.exp(1j*np.matmul(np.transpose(kpoint), translations[i]))
    print(phase)
    print(phase2)'''
    
    for j in range(num_atoms):
        new_pos=np.matmul(r,pos[j]) + t
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
            if np.array_equal(np.round(new_pos,decimals=sym_prec),np.round(pos[k],decimals=sym_prec)):
                permu[k,j]=1*np.round(np.exp(2j*np.pi*np.matmul(np.transpose(kpoint), cells[k])),5)
    return permu 


def transformmode(pos,r,t,q,mode): 
    permu = permutation_matrix(pos,r,t,q)
    newmode = permu @ mode
    
    for i in range(len(newmode)):
        newmode[i] = r @ newmode[i] 
    
    return newmode
    
scf_file = r'AgP2\Phonons_V2\AgP2_scf.in'
structure = io.read(scf_file)

basisvec = structure.get_cell()
atmpos= structure.get_scaled_positions()
atmnum = structure.get_atomic_numbers()
recbasis = structure.cell.reciprocal()



kpoint = [-1/2, 1/2, 0]  
irreps, rotations, translations, mapping_little_group = get_spacegroup_irreps(
    basisvec, atmpos, atmnum, kpoint
)

modes = readModesatKpoin(kpoint,r'AgP2\Phonons_V2\matdyn.modes', scffile=r'AgP2\Phonons_V2\AgP2.scf.pwo')


#print(modes[0])
#print(transformmode(atmpos,rotations[1],translations[1],kpoint,modes[0]))

#print(len(rotations), len(translations), mapping_little_group,len(irreps))
#print(irreps)

#print(get_character(irreps[1]))
s =check_spacegroup_representation(rotations[mapping_little_group], translations[mapping_little_group], kpoint,irreps[0])


ds = get_displacements_representation(basisvec, atmpos, rotations[mapping_little_group], translations[mapping_little_group], kpoint)
print(np.shape(ds))
print(np.trace(ds[0]))

#print(len(ds))
#print(np.shape(project_to_irrep(ds, irreps[0])))
#print(get_character(ds))
#print(decompose_representation(rotations, max_num_random_generations=4))
