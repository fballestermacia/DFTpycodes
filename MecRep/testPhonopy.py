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

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



def smallest_vectors(spos, pos, sbasis, symprec=1e-5):

    snatoms = len(spos)
    natoms = len(pos)
    
    smallestvectors = np.zeros((snatoms, natoms, 27, 3))
    multi = np.zeros((snatoms,natoms))
    
    '''
    This next part is forked from phonopy, i do not really understand it but it works

    ' Si non confectus, no reficiat ' 

    '''
    # For each vector, we will need to consider all nearby images in the
    # reduced bases. The lattice points at which supercell images are
    # searched are composed by linear combinations of three vectors in (0,
    # a, b, c, -a-b-c, -a, -b, -c, a+b+c). There are finally 65 lattice
    # points. There is no proof that this is enough.
    lattice_1D = (-1, 0, 1)
    lattice_4D = np.array(
        [
            [i, j, k, ll]
            for i in lattice_1D
            for j in lattice_1D
            for k in lattice_1D
            for ll in lattice_1D
        ]
    )
    bases = [[1, 0, 0], [0, 1, 0], [0, 0, 1], [-1, -1, -1]]
    lattice_points = np.dot(lattice_4D, bases)
    lattice_points = np.array(
        np.unique(lattice_points, axis=0),
    )

    npoints = len(lattice_points)
    vecs = np.empty((npoints,3))

    for i in range(snatoms):
        for j in range(natoms):
            length = np.empty(npoints)
            for k in range(npoints):
                vecs[k] = spos[i]-pos[j]+lattice_points[k]
                length[k] = np.linalg.norm(vecs[k])
            
            minimum = None
            for k in range(npoints):
                if minimum is None:
                    minimum = length[k]
                elif length[k] < minimum:
                    minimum = length[k]
            
            count = 0
            for k in range(npoints):
                if np.abs(length[k] - minimum) < symprec:
                    transformedvec = np.matmul(sbasis,vecs[k])
                    smallestvectors[i,j,count,:] = transformedvec
                    count += 1
            if count > 27:
                print("As Eminem said: 'something's wrong'")
                return None, None
            else: 
                multi[i,j] = count
    return smallestvectors, multi






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

'''primcell_filename = 'data/232/AgP2.scf.pwi'
fc_filename = 'data/232/AgP2.fc'
cell, _ = read_pwscf(primcell_filename)
q2r = PH_Q2R(fc_filename)
q2r.run(cell)
q2r.write_force_constants(fc_format='txt')'''


dyn = CC.Phonons.Phonons()
dyn.LoadFromQE(fildyn_prefix="data/232/dynmats/AgP2.dyn", nqirr=8)

dyn.InterpolateMesh([7,7,7])

#dyn.Symmetrize(verbose=True)

super_dyn = dyn.GenerateSupercellDyn(dyn.GetSupercell())
#super_dyn.Symmetrize(verbose=True, asr='crystal')

dyn.save_phonopy(units_ev_ang2=False, write_unitcell=True, write_poscar=False)#, supercell_size=dyn.GetSupercell())


phonon = ph.load(force_constants_filename='FORCE_CONSTANTS',unitcell_filename='unitcell.in',store_dense_svecs=False, calculator='qe') #ph.load('phonopy_params.yaml')




primitive_basis, primitive_pos, primitive_types = phonon._primitive.cell, phonon._primitive.get_positions(), phonon._primitive.get_atomic_numbers()#spg.find_primitive((dyn.structure.unit_cell,dyn.structure.coords, dyn.structure.get_atomic_types()))

#print(primitive_basis, dyn.structure.unit_cell)


num_atom = int(len(dyn.structure.coords))
num_satom = int(len(super_dyn.structure.coords))
#print(num_satom)


svecs, multi = smallest_vectors(
    CC.Methods.covariant_coordinates(super_dyn.structure.unit_cell, super_dyn.structure.coords),
    CC.Methods.covariant_coordinates(super_dyn.structure.unit_cell, dyn.structure.coords),
    super_dyn.structure.unit_cell
    )
    

#smallest_vectorsphpy, multi = get_smallest_vectors(super_dyn.structure.unit_cell,
#                                               CC.Methods.covariant_coordinates(super_dyn.structure.unit_cell, super_dyn.structure.coords),
#                                               CC.Methods.covariant_coordinates(super_dyn.structure.unit_cell, dyn.structure.coords))#phonon._primitive.get_smallest_vectors()#dmat._pcell.get_smallest_vectors()

mass = super_dyn.structure.get_masses_array()*MASS_RY_TO_UMA

print(multi)

reduced_bases = dyn.structure.unit_cell#get_reduced_bases(super_dyn.structure.unit_cell, 1e-5)
positions = CC.Methods.covariant_coordinates(reduced_bases, dyn.structure.coords) #np.dot(dyn.structure.coords, np.linalg.inv(reduced_bases))

super_pos = CC.Methods.covariant_coordinates(reduced_bases, super_dyn.structure.coords) #np.dot(super_dyn.structure.coords, np.linalg.inv(reduced_bases))

p2s_map = np.arange(len(primitive_types))#dmat._p2s_map #= primitive.get_primitive_to_supercell_map()
s2p_map = s2p(super_pos,primitive_basis,p2s_map)#dmat._s2p_map #= primitive.get_supercell_to_primitive_map()



num_patom = num_atom

get_phonon_hr(
        phonon.force_constants*factor**2,
        svecs,
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




