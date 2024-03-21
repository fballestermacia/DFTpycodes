import sys
from phonopy.interface.qe import read_pwscf, PH_Q2R
import numpy as np
import spglib as spg
import ase.io as io
import ase as ase
import ase.spacegroup as asespg
import atomman as am
import phonopy as ph
from phonopy.structure.cells import get_reduced_bases, determinant, get_smallest_vectors

from jarvis.io.phonopy.fcmat2hr import get_phonon_hr
from jarvis.io.phonopy.outputs import get_phonon_tb, read_fc
from jarvis.core.atoms import ase_to_atoms

from phonopy.units import PwscfToTHz, VaspToTHz
from phonopy.interface.qe import write_supercells_with_displacements

import cellconstructor as CC
import cellconstructor.Phonons
from cellconstructor.Units import MASS_RY_TO_UMA


def s2p(spos, cell, p2s):
    symprec = 1e-5
    frac_pos = spos

    p2s_positions = frac_pos[p2s]
    s2p_map = []
    for s_pos in frac_pos:
        # Compute distances from s_pos to all positions in _p2s_map.
        frac_diffs = p2s_positions - s_pos
        frac_diffs -= np.rint(frac_diffs)
        cart_diffs = np.dot(frac_diffs, cell)
        distances = np.sqrt((cart_diffs**2).sum(axis=1))
        indices = np.where(distances < symprec)[0]
        assert len(indices) == 1
        s2p_map.append(p2s[indices[0]])

    s2p_map = np.array(s2p_map, dtype="intc")
    

    return s2p_map


factor = PwscfToTHz
#print(factor**2*0.284552630079000/30.973761)

cellphpy, _ = read_pwscf(r'data/232/AgP2.scf.pwi')
#print(cell)
q2r = PH_Q2R(r'data/232/AgP2.fc')

q2r.run(cellphpy, is_full_fc=True)#, parse_fc=True)
#print(q2r.primitive)
q2r.write_force_constants(fc_format='text')




#system = am.load('phonopy_Atoms', cell)
#aseatm = system.dump('ase_Atoms')
#jarvisatm =ase_to_atoms(aseatm)

phonon = ph.load(force_constants_filename='FORCE_CONSTANTS',unitcell_filename='data/232/AgP2.scf.pwi',store_dense_svecs=False, calculator='qe') #ph.load('phonopy_params.yaml')

dyn = CC.Phonons.Phonons()
dyn.LoadFromQE(fildyn_prefix="data/232/dynmats/AgP2.dyn", nqirr=8)

super_dyn = dyn.GenerateSupercellDyn(dyn.GetSupercell())


#supercell = phonon.supercell
#print(supercell.get_cell(),'\n', super_dyn.structure.unit_cell)
#primitive = phonon.get_primitive()
#print(primitive)


num_atom = int(len(dyn.structure.coords))
num_satom = int(len(super_dyn.structure.coords))
#print(num_satom)


smallest_vectors, multi = get_smallest_vectors(super_dyn.structure.unit_cell,super_dyn.structure.coords,dyn.structure.coords)#phonon._primitive.get_smallest_vectors()#dmat._pcell.get_smallest_vectors()
mass = dyn.structure.get_masses_array()*MASS_RY_TO_UMA
#print(phonon._primitive.get_smallest_vectors())
print(mass)
print()
reduced_bases = dyn.structure.unit_cell#get_reduced_bases(super_dyn.structure.unit_cell, 1e-5)
positions = CC.Methods.covariant_coordinates(reduced_bases, dyn.structure.coords) #np.dot(dyn.structure.coords, np.linalg.inv(reduced_bases))

super_pos = CC.Methods.covariant_coordinates(reduced_bases, super_dyn.structure.coords) #np.dot(super_dyn.structure.coords, np.linalg.inv(reduced_bases))

p2s_map = np.arange(num_atom)#dmat._p2s_map #= primitive.get_primitive_to_supercell_map()
s2p_map = s2p(super_pos,dyn.structure.unit_cell,p2s_map)#dmat._s2p_map #= primitive.get_supercell_to_primitive_map()

print(dyn.structure.unit_cell,positions)

print(p2s_map)
print(s2p_map)
num_patom = num_atom
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
        'DFTpycodes/WTProj/AgP2_phonopyTB_hr.dat'
    )

print("donete")
#get_phonon_tb(atoms=jarvisatm, fc=read_fc("FORCE_CONSTANTS"))



