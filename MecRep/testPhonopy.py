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

from datetime import datetime



def smallest_vectors(spos, pos, sbasis, symprec=1e-8):

    snatoms = len(spos)
    natoms = len(pos)
    
    smallestvectors = np.zeros((snatoms, snatoms, 27, 3))
    multi = np.zeros((snatoms,snatoms))
    
    '''
    This next part is forked from phonopy, i do not really understand it but it works

    ' Si non confectus, no reficiat ' 

    '''
    
    #THIS PART IS FROM PHONOPY, IT IS MORE EFFICIENT BUT I DOUBT IF IT'S CORRECT
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
    #print(lattice_points)
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
                    transformedvec = CC.Methods.covariant_coordinates(np.linalg.inv(sbasis), vecs[k])#np.matmul(sbasis,vecs[k])
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


def write_phonon_hr(fcmat, svecs, mass, multi, super_pos, p2s_map, s2p_map, num_satom, num_patom,sbasis, outfile):
    '''
    Write phonon tight-binding Hamiltonian in the format of wannier90_hr.dat
    Copied and modified from Jarvis-Tools, mainly because of issues with the change of basis
    '''
    # maximum dimension for hr matrix
    ndim = 51
    sdim = 20
    nrpt_max = 51 ** 3
    # hr matrix
    norbs = num_patom * 3
    hr_mat = np.zeros((ndim, ndim, ndim, norbs, norbs), dtype="complex128")
    # hr_mat = np.zeros((ndim, ndim, ndim, norbs, norbs), dtype=np.complex128)
    hr_mat0 = np.zeros((nrpt_max, norbs, norbs), dtype="complex128")
    # hr_mat0 = np.zeros((nrpt_max, norbs, norbs), dtype=np.complex128)
    # WS points
    rpts = np.zeros((nrpt_max, 3), dtype="int32")
    # rpts = np.zeros((nrpt_max, 3), dtype=np.int32)
    # degeneracy
    dege = np.zeros((nrpt_max), dtype="int32")
    # dege = np.zeros((nrpt_max), dtype=np.int32)

    for iatom in range(num_patom):  # atoms in primitive cell
        for jatom in range(num_patom):  # atoms in primitive cell
            mass_sqrt = np.sqrt(mass[iatom] * mass[jatom])
            for katom in range(num_satom):  # atoms in supercell
                if s2p_map[katom] != p2s_map[jatom]:
                    continue
                for rr in range(np.int32(multi[katom, iatom])):
                    # find which rpt
                    rvec = (
                        svecs[katom, iatom, rr]
                        + super_pos[p2s_map[iatom]]
                        - super_pos[s2p_map[katom]]
                    )
                    print(svecs[katom, iatom, rr],super_pos[p2s_map[iatom]],super_pos[s2p_map[katom]])
                    print(iatom, jatom, katom, rvec,CC.Methods.covariant_coordinates(np.linalg.inv(sbasis), rvec))
                    rvec = CC.Methods.covariant_coordinates(sbasis, rvec)
                    for ii in range(3):
                        if abs(rvec[ii]) < 1.0e-5:
                            rvec[ii] = 0.0
                    rx = np.int32(np.round(rvec[0]))#np.int32(np.round(np.dot(rvec,sbasis[0])/np.linalg.norm(sbasis[0])))
                    ry = np.int32(np.round(rvec[1]))#np.int32(np.round(np.dot(rvec,sbasis[1])/np.linalg.norm(sbasis[1])))#np.int32(np.round(rvec[1]))
                    rz = np.int32(np.round(rvec[2]))#np.int32(np.round(np.dot(rvec,sbasis[2])/np.linalg.norm(sbasis[2])))#np.int32(np.round(rvec[2]))
                    print(rx,ry,rz)
                    idx = iatom * 3
                    idy = jatom * 3
                    rx = rx + sdim
                    ry = ry + sdim
                    rz = rz + sdim
                    hr_mat[rx, ry, rz, idx : idx + 3, idy : idy + 3] = fcmat[
                        p2s_map[iatom], katom
                    ] / (multi[katom, iatom] * mass_sqrt)
                    #print(multi[katom, iatom])

    # collect all hr at rpts with none-zero data
    irpt = 0  # count for useful rpts
    for rx in range(-sdim, sdim + 1):
        ixo = rx + sdim
        for ry in range(-sdim, sdim + 1):
            iyo = ry + sdim
            for rz in range(-sdim, sdim + 1):
                izo = rz + sdim
                #print(abs(hr_mat[ixo, iyo, izo, :, :]).sum())
                if (
                    abs(hr_mat[ixo, iyo, izo, :, :]).sum() < 1.0e-17
                ):  # ommit too small
                    continue
                dege[irpt] += 1
                rpts[irpt, 0] = rx
                rpts[irpt, 1] = ry
                rpts[irpt, 2] = rz
                hr_mat0[irpt, :, :] = hr_mat[ixo, iyo, izo, :, :]
                irpt = irpt + 1
    nrpt = irpt
    dege_rpts = dege[0:nrpt]
    norbs = num_patom * 3
    with open(outfile, "w") as f:
        line = (
            " Writen on "
            + str(datetime.now())
            + "\n"
            + "          "
            + str(norbs)
            + "\n"
            + "        "
            + str(nrpt)
            + "\n"
        )
        f.write(line)
        nl = np.int32(np.ceil(nrpt / 15.0))
        for n in range(nl):
            line = "    " + "    ".join(
                [str(np.int32(i)) for i in dege_rpts[n * 15 : (n + 1) * 15]]
            )
            f.write(line)
            f.write("\n")
        for irpt in range(nrpt):
            rx = rpts[irpt, 0]
            ry = rpts[irpt, 1]
            rz = rpts[irpt, 2]
            for jatomorb in range(norbs):
                for iatomorb in range(norbs):
                    rp = hr_mat0[irpt, iatomorb, jatomorb].real
                    ip = hr_mat0[irpt, iatomorb, jatomorb].imag
                    tmp = "{:8d}{:8d}{:8d}{:8d}{:8d}{:20.10f}{:20.10f}\n"
                    line = tmp.format(
                        rx, ry, rz, iatomorb + 1, jatomorb + 1, rp, ip
                    )
                    f.write(line)
    
    



factor = PwscfToTHz

'''primcell_filename = 'data/232/AgP2.scf.pwi'
fc_filename = 'data/232/AgP2.fc'
cell, _ = read_pwscf(primcell_filename)
q2r = PH_Q2R(fc_filename)
q2r.run(cell)
q2r.write_force_constants(fc_format='txt')'''


dyn = CC.Phonons.Phonons()
dyn.LoadFromQE(fildyn_prefix="data/232/dynmats/AgP2.dyn", nqirr=8)


'''print(dyn.structure.atoms)
print(dyn.structure.coords)
print(dyn.structure.unit_cell)'''

#dyn.InterpolateMesh([10,10,10])

#dyn.Symmetrize(verbose=True) #FC is already symmetric, although it migh use another asr

super_dyn = dyn.GenerateSupercellDyn(dyn.GetSupercell())
#super_dyn.Symmetrize(verbose=True, asr='crystal')

#print(np.shape(dyn.GetRealSpaceFC(supercell_array=dyn.GetSupercell())))


dyn.save_phonopy(units_ev_ang2=False, write_unitcell=True, write_poscar=False)#, supercell_size=dyn.GetSupercell())


phonon = ph.load(force_constants_filename='FORCE_CONSTANTS',unitcell_filename='unitcell.in',store_dense_svecs=False, calculator='qe') #ph.load('phonopy_params.yaml')




primitive_basis, primitive_pos, primitive_types = phonon._primitive.cell, phonon._primitive.get_positions(), phonon._primitive.get_atomic_numbers()#spg.find_primitive((dyn.structure.unit_cell,dyn.structure.coords, dyn.structure.get_atomic_types()))

#print(primitive_basis, dyn.structure.unit_cell)


num_atom = int(len(dyn.structure.coords))
num_satom = int(len(super_dyn.structure.coords))

cell = dyn.structure.unit_cell#np.array([[1.000000000,0.000000000,0.000000000],[0.000000000,0.817188282,0.000000000],[-0.506471883,0.000000000,1.156604712]])#dyn.structure.unit_cell
scell = super_dyn.structure.unit_cell

print('hello there!')

svecs, multi = smallest_vectors(
    CC.Methods.covariant_coordinates(scell, super_dyn.structure.coords),
    CC.Methods.covariant_coordinates(scell, super_dyn.structure.coords),
    scell #scell
    )

print('general Kenobi')

#smallest_vectorsphpy, multi = get_smallest_vectors(super_dyn.structure.unit_cell,
#                                               CC.Methods.covariant_coordinates(super_dyn.structure.unit_cell, super_dyn.structure.coords),
#                                               CC.Methods.covariant_coordinates(super_dyn.structure.unit_cell, dyn.structure.coords))#phonon._primitive.get_smallest_vectors()#dmat._pcell.get_smallest_vectors()

mass = super_dyn.structure.get_masses_array()*MASS_RY_TO_UMA




positions = CC.Methods.covariant_coordinates(cell, dyn.structure.coords) #np.dot(dyn.structure.coords, np.linalg.inv(reduced_bases))

super_pos = CC.Methods.covariant_coordinates(cell, super_dyn.structure.coords) #np.dot(super_dyn.structure.coords, np.linalg.inv(reduced_bases))

p2s_map = np.arange(len(primitive_types))#dmat._p2s_map #= primitive.get_primitive_to_supercell_map()
s2p_map = s2p(super_pos,primitive_basis,p2s_map)#dmat._s2p_map #= primitive.get_supercell_to_primitive_map()





num_patom = num_atom

#print(np.transpose(np.reshape(dyn.GetRealSpaceFC(supercell_array=dyn.GetSupercell()),(num_satom,3,num_satom,3)),axes=(0,2,1,3))[70,130])

write_phonon_hr(
        np.transpose(np.reshape(dyn.GetRealSpaceFC(supercell_array=dyn.GetSupercell()),(num_satom,3,num_satom,3)),axes=(0,2,1,3))*factor**2,#phonon.force_constants*factor**2,
        svecs,
        mass,
        multi,
        super_dyn.structure.coords,
        p2s_map,
        s2p_map,
        num_satom,
        num_patom,
        cell,
        'DFTpycodes/WTProj/AgP2_phonopyTB_hr.dat'
    )

print("donete")




