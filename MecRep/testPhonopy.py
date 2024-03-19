import sys
from phonopy.interface.qe import read_pwscf, PH_Q2R
import numpy as np
import spglib as spg
import ase.io as io
import ase as ase
import ase.spacegroup as asespg
import atomman as am
import phonopy as ph
from phonopy.structure.cells import get_reduced_bases, determinant

from jarvis.io.phonopy.fcmat2hr import get_phonon_hr
from jarvis.io.phonopy.outputs import get_phonon_tb, read_fc
from jarvis.core.atoms import ase_to_atoms

from phonopy.units import PwscfToTHz, VaspToTHz
from phonopy.interface.qe import write_supercells_with_displacements

factor = PwscfToTHz
print(factor**2*0.284552630079000/30.973761)

cell, _ = read_pwscf(r'data/AgP2/Phonons/AgP2.scf.pwi')
#print(cell)
q2r = PH_Q2R(r'data/AgP2/Phonons/AgP2.fc')
q2r.run(cell, is_full_fc=True)#, parse_fc=True)
q2r.write_force_constants(fc_format='text')




#system = am.load('phonopy_Atoms', cell)
#aseatm = system.dump('ase_Atoms')
#jarvisatm =ase_to_atoms(aseatm)

phonon = ph.load(force_constants_filename='FORCE_CONSTANTS',unitcell_filename='data/AgP2/Phonons/AgP2.scf.pwi',store_dense_svecs=False, calculator='qe') #ph.load('phonopy_params.yaml')




supercell = phonon.get_supercell()
#print(supercell)
primitive = phonon.get_primitive()


num_atom = int(cell.get_number_of_atoms())
num_satom = int(supercell.get_number_of_atoms())

dmat = phonon._dynamical_matrix
smallest_vectors, multi = phonon._primitive.get_smallest_vectors()#dmat._pcell.get_smallest_vectors()
mass = dmat._pcell.get_masses()
#print(phonon._primitive.get_smallest_vectors())
reduced_bases = get_reduced_bases(supercell.get_cell(), 1e-5)
positions = np.dot(supercell.get_positions(), np.linalg.inv(reduced_bases))

relative_scale = np.dot(reduced_bases, np.linalg.inv(primitive.get_cell()))
super_pos = np.zeros((num_satom, 3), dtype=np.float64)
for i in range(num_satom):
    super_pos[i] = np.dot(positions[i], relative_scale)
p2s_map = dmat._p2s_map = primitive.get_primitive_to_supercell_map()
s2p_map = dmat._s2p_map = primitive.get_supercell_to_primitive_map()
num_satom = supercell.get_number_of_atoms()
num_patom = primitive.get_number_of_atoms()
#print(smallest_vectors[1])
get_phonon_hr(
        phonon.force_constants*factor**2,#read_fc("FORCE_CONSTANTS")*factor**2,
        smallest_vectors,
        mass,
        multi,
        super_pos,
        p2s_map,
        s2p_map,
        num_satom,
        num_patom,
        'phonopyTB_hr.dat'
    )

print("donete")
#get_phonon_tb(atoms=jarvisatm, fc=read_fc("FORCE_CONSTANTS"))



