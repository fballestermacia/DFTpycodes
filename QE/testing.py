import numpy as np
import spglib as spg
import ase.io as io
import ase as ase
import ase.spacegroup as asespg
from spgrep import get_spacegroup_irreps
#from ase.visualize import view
import phonopy as ph
import utilsQE
import atomman as am
from phonopy.interface.qe import read_pwscf, PH_Q2R


def similarity_transformation(rot, mat):
    """ R x M x R^-1 """
    return np.dot(rot, np.dot(mat, np.linalg.inv(rot)))


        

def findlg(rk,t,q,symprec=1e-5):
    '''tag=[]
    rotatq=[]
    tatq = []

    #print(rk)
    #print(len(rk))
    for i,r in enumerate(rk):
        diff = np.dot(q,r) - q
        if (np.abs(diff) < symprec).all():
            rotatq.append(r)
            for j in range(3):
                if np.abs(t[i][j] - 1) < symprec:
                    t[i][j] = 0.0
            tatq.append(t[i])
            tag.append(i)'''
    irreps, rotationslgroup, translationslgroup, mapping_little_group = get_spacegroup_irreps(
                basisvec, atmpos, atmnum, symk[ik]
            )
    tag = mapping_little_group
    rotatq = rotationslgroup
    tatq = translationslgroup
    return tag, np.array(rotatq), np.array(tatq)

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

def ground_matrix(q,rot,t,lat,pos):
    matrices = []
    tags, rotatq, tatq = findlg(rot,t,q)
    lat,rec = utilsQE.readcellvec(r'AgP2\Phonons_V2\AgP2.scf.pwo')
    for i in range(len(rotatq)):
        r_cart = similarity_transformation(rotatq[i],lat)
        
        perm_mat = permutation_matrix(pos,rotatq[i],tatq[i],q)
        
        matrices.append(np.kron(perm_mat, r_cart))
    return np.array(matrices)
        

def band_indices(bandsatqpoint, tol=1e-1):
    indices = []
    currentden = []
    empty =True
    for i in range(len(bandsatqpoint)):
        
        if np.abs(bandsatqpoint[i]-bandsatqpoint[i-1])<tol:
            currentden.append(i+1)
            empty = False
        
        else:
            if empty:
                currentden.append(i+1)
                empty = False
            else:
                indices.append(currentden)
                currentden = []
                currentden.append(i+1)
                empty =False
        
    indices.append(currentden)
    return indices

def get_irreps(dynmat, bandindices,qpoint,rots,trans,atmpos,basisvec):
        eigvecs = []
        phases = np.kron(
            [np.exp(1j * np.pi * np.dot(qpoint, pos))
             for pos in atmpos], [1, 1, 1])
        
        eigvecs2 = [x for x in np.linalg.eig(dynmat)[1]]
        
        for vec in np.transpose(eigvecs2): 
            eigvecs.append(vec * phases)
        
        irrep = []
        for indices in bandindices:
            irrep_Rs = []
    
            for mat in ground_matrix(qpoint,rots,trans,basisvec,atmpos):
                l = len(indices)

                if l == 1:
                    vec = eigvecs[indices[0]-1]
                    irrep_Rs.append([[np.vdot(vec, np.dot(mat, vec))]])
                    continue

                irrep_R = np.zeros((l, l), dtype=complex)
                for i, b_i in enumerate(indices):
                    
                    vec_i = eigvecs[b_i]
                    for j, b_j in enumerate(indices):
                        vec_j = eigvecs[b_j]
                        irrep_R[i, j] = np.vdot(vec_i, np.dot(mat, vec_j))
                
                irrep_Rs.append(irrep_R)

            irrep.append(irrep_Rs)

        return irrep

def get_irreps2(dynmat, bandindices,qpoint,rots,trans,atmpos,basisvec):
        eigvecs = []
        phases = np.kron(
            [np.exp(2j * np.pi * np.dot(qpoint, pos))
             for pos in atmpos], [1, 1, 1])
    
        
        eigvecs2 = [x for x in np.linalg.eig(dynmat)[1]]
        
        for vec in np.transpose(eigvecs2): 
            eigvecs.append(vec * phases)
        
        irrep = []
        for indices in bandindices:
            irrep_Rs = []
    
            for mat in ground_matrix(qpoint,rots,trans,basisvec,atmpos):
                l = len(indices)

                if l == 1:
                    vec = eigvecs[indices[0]-1]
                    irrep_Rs.append([[np.vdot(vec, np.dot(mat, vec))]])
                    continue

                irrep_R = np.zeros((l, l), dtype=complex)
                for i, b_i in enumerate(indices):
                    
                    vec_i = eigvecs[b_i]
                    for j, b_j in enumerate(indices):
                        vec_j = eigvecs[b_j]
                        irrep_R[i, j] = np.vdot(vec_i, np.dot(mat, vec_j))
                
                irrep_Rs.append(irrep_R)

            irrep.append(irrep_Rs)

        return irrep

if __name__ == '__main__':
    
    factor = 0.123983#ph.units.PwscfToTHz*4.136 # meV

    
    # Symmetry info and primitive cell
    scf_file = r'AgP2\Phonons_V2\AgP2_scf.in'
    structure = io.read(scf_file)
    
    basisvec = structure.get_cell()
    atmpos= structure.get_scaled_positions()
    atmnum = structure.get_atomic_numbers()
    recbasis = structure.cell.reciprocal()
    
    sym_data = spg.get_symmetry_dataset((basisvec,atmpos,atmnum))
    rotations = sym_data.get("rotations")
    translations = sym_data.get("translations")
    SG = sym_data.get("number")
    nsym = len(rotations)
    
    # K-points and frequencies
    qlabels, positions, qpoints = utilsQE.readHighSymPointsPhonon(r"AgP2\Phonons_V2\matdyn_AgP2.in", retKpoints=True)
    notqpoints, bands = utilsQE.readPhononbandFreq(r"AgP2\Phonons_V2\AgP2.freq.gp")
    print(qlabels)
    #print(bands[:,positions])
    
    asg = asespg.Spacegroup(SG)
    
    gammaindex = qpoints.index([0,0,0])
    
    dummyqo = []
    dummyqo.append([0.,0.,0.])
    dummypos = []
    dummypos.append(positions[gammaindex])
    for i in range(len(qpoints)):
        if i == gammaindex:
            continue
        dummyqo.append(qpoints[i])
        dummypos.append(positions[i])
    
    dummyqo2, indices = np.unique(qpoints, axis=0, return_index=True)
    newqpoints = []
    newpos = []
    
    for i in np.sort(indices):
        newqpoints.append(dummyqo[i])
        newpos.append(dummypos[i])
        
    
    #print(newpos)
    newbands = np.sort(np.transpose(bands[:,newpos]), axis=0)
    
    
    bandindex = []
    
    for i in range(len(newbands)):
        bandindex.append(band_indices(newbands[i]))
    
    
    newbands *= factor 
    
    for i in range(len(newbands)):
        for j in range(len(newbands[i])):
            if abs(newbands[i][j]) < 5e-2:
                newbands[i][j] = 0.
    
    syst = am.load('ase_Atoms', structure)
    phpyAT = am.dump('phonopy_Atoms', syst)
    phonon = ph.Phonopy(phpyAT, np.eye(3),asg.scaled_primitive_cell)
    phpysym = phonon.get_symmetry()
    #cell = phonon.get_unitcell()
    cell, _ = read_pwscf(r'AgP2\Phonons_V2\AgP2.scf.pwi')
    phonon = ph.load(force_constants_filename=r"DFTpycodes-main\QE\FORCE_CONSTANTS", unitcell=cell, symmetrize_fc=True)
    
    phonon.run_qpoints(newqpoints)
    
    '''phonon.run_band_structure(qpoints)
    phonon.get_band_structure()'''
    print(np.shape(phonon.get_dynamical_matrix_at_q(newqpoints[0])))
    #print(np.linalg.eig(phonon.get_dynamical_matrix_at_q(newqpoints[0])))
    #print('this',get_irreps(dmat, bandindex[0], newqpoints[0],rotations,translations,atmpos,basisvec))
    
    with open(r"DFTpycodes-main\QE\testTRACESC.txt",'w') as f:
        f.write(str(len(bands))+'\n')
        f.write('0'+'\n')
        
        f.write(str(len(rotations))+'\n')
        for i in range(len(rotations)):
            for j in range(3):
                f.write(" {:3d} {:3d} {:3d}".format(rotations[i,j,0],rotations[i,j,1],rotations[i,j,2]))
            f.write(" {:10.5f} {:10.5f} {:10.5f}".format(translations[i,0],translations[i,1],translations[i,2]))
            f.write(" 0 0 0 0 0 0 0 0" + '\n') #TODO: WRITE THE SPINOR ROTATIONS ALTHOUGH THEY ARE IRRELEVANT IN PHONONS

        f.write(str(len(newqpoints))+'\n')
        for i in range(len(newqpoints)):
            for j in range(3):
                f.write("{:14.7}".format(newqpoints[i][j]))
            f.write('\n')
            
            
        #######COPIED FROM phonopy2trace######### 
        symk=newqpoints
        nsymk = len(newqpoints)
        for ik in range(nsymk):
            dmat = phonon.get_dynamical_matrix_at_q(newqpoints[ik])
            irrs = get_irreps(dmat, bandindex[ik], newqpoints[ik],rotations,translations,atmpos,basisvec)
            
            kpt=symk[ik]
         #   little group
            '''rk=phonon.set_irreps(symk[0], 1e-4) #TODO FIND A WAY TO FIND THE LITTLE GROUP AT A GIVEN KPOINT
            cct = phonon.get_irreps()
            ct.get_ground_matrices()
            setirep=phonon.set_irreps(kpt, 1e-5)
            ct = phonon.get_irreps()
            print(ct._get_characters())'''
            
            
            lgtag,_,_=findlg(rotations,translations,kpt)
            f.write(str(len(lgtag))+'\n')
            for opr in lgtag:
                f.write(str(opr+1)+'  ')
            f.write('\n')
        #   phonon.show_irreps()
            #band_indices = ct.get_band_indices()
            #print(band_indices)
            #ek1 = ct.get_freq()
            #ek = []
            #for eig in ek1:
            #    ek.append(eig*factor)
        #   characters = np.rint(ct.get_characters())#.real
        
            for l, bi in enumerate(bandindex[ik]):
                f.write('{:5s}{:3s}{:12.6f}'.format(str(bi[0]),str(len(bi)),newbands[ik][bi[0]-1]))
                for iopr in range(len(lgtag)):
                    #print(irrs[l][iopr])
                    f.write('{:10.5f}{:10.5f}'.format(np.round(np.real(np.trace(irrs[l][iopr])),7) ,np.round(np.imag(np.trace(irrs[l][iopr])),7)))
                f.write('\n')
            