import numpy as np
import spglib as spg
import ase.io as io
import ase as ase
from spgrep import get_spacegroup_irreps
import utilsQE
#from cellconstructor.Methods import covariant_coordinates
import su2rot
#np.set_printoptions(threshold=np.inf)

def band_indices(bandsatqpoint, tol=1e-3):
    '''indices = []
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
        
    indices.append(currentden)'''

    indices = []
    done = []
    for i in range(len(bandsatqpoint)):
        if i in done:
            continue
        else:
            f_set = [i]
            done.append(i)
        for j in range(i + 1, len(bandsatqpoint)):
            if (np.abs(bandsatqpoint[f_set] - bandsatqpoint[j]) < tol).any():
                f_set.append(j)
                done.append(j)
        indices.append(f_set[:])

    return indices


def permutation_matrix(pos,r,t,q,sym_prec=4):
    
    num_atoms = len(pos)
    permu = np.zeros([num_atoms,num_atoms], dtype='complex')
    
    for j in range(num_atoms):
        new_pos=np.matmul(r,pos[j]) + t

        for j2 in range(num_atoms):
            diff = new_pos-pos[j2]
            if (abs(diff - np.rint(diff)) < 10**(-1*sym_prec)).all():
                phase = np.dot(q, new_pos-pos[j2])

                permu[j2,j] = np.exp(-2j*np.pi*phase)

    return permu 


def similarity_transformation(U, mat):
    """ U x M x U^-1 """
    return np.matmul(U, np.matmul(mat, np.linalg.inv(U)))


def symmetryeigval(modefunc,r,t, latticevec, permatq, bandindices):
    U = latticevec
    O = similarity_transformation(U,r)
    symoperator = np.kron(permatq,O)
    
    l = len(bandindices)
    if l == 1:
        eigval =np.matmul(np.transpose(np.conjugate(modefunc[bandindices].flatten())),np.matmul(symoperator,modefunc[bandindices].flatten()))
    
    else:
        mat = np.zeros((l,l), dtype=np.complex128)
        for i, b_i in enumerate(bandindices):
            vec_i = modefunc[b_i].flatten()
            for j, b_j in enumerate(bandindices):
                vec_j = modefunc[b_j].flatten()
                mat[i, j] = np.vdot(vec_i, np.dot(symoperator, vec_j))
        
        eigval = np.trace(mat)
    return eigval


def findlg(rk,t,q,symprec=1e-5):
    tag=[]
    rotatq=[]
    tatq = []

    #print(rk)
    #print(len(rk))
    for i,r in enumerate(rk):
        diff = np.dot(r,q) - q
        if (np.abs(diff) < symprec).all():
            rotatq.append(r)
            for j in range(3):
                if np.abs(t[i][j] - 1) < symprec:
                    t[i][j] = 0.0
            tatq.append(t[i])
            tag.append(i)
    return tag, np.array(rotatq), np.array(tatq)


if __name__ == '__main__':
    
    factor = 0.123983
    
    scf_file = r'data/AgP2/Phonons/noLOTO/AgP2.scf.pwi'
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
    SU2=np.zeros((nsym,2,2),dtype=np.complex128)
    for i in range(nsym):
        sopr_c=similarity_transformation(np.transpose(basisvec),rotations[i,:,:])
        SU2[i,:,:]=su2rot.get_su2rotation(sopr_c)

    
    
    qlabels, positions, qpoints = utilsQE.readHighSymPointsPhonon(r"data/AgP2/Phonons/444/matdyn.in", retKpoints=True)
    notqpoints, bands = utilsQE.readPhononbandFreq(r"data/AgP2/Phonons/444/AgP2.newHSP.freq.gp")
    

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
    
    
    
    newbands = np.transpose(bands[:,newpos])

    
    
    bandindex = []
    
    for i in range(len(newbands)):
        #print(newbands)
        bandindex.append(band_indices(newbands[i],tol=5e-3))
        #stop = stop
    
    newbands *= factor
    
    for i in range(len(newbands)):
        for j in range(len(newbands[i])):
            if abs(newbands[i][j]) < 5e-2:
                newbands[i][j] = 0.
    
    
    
    
    atmposcov = atmpos#covariant_coordinates(basisvec[:],atmpos)
    #print(atmposcov)

    with open(r"data/AgP2/Phonons/444/traces.txt",'w') as f:
        f.write(str(len(bands))+'\n')
        f.write('0'+'\n')
        
        f.write(str(len(rotations))+'\n')
        for i in range(len(rotations)):
            for j in range(3):
                f.write(" {:3d} {:3d} {:3d}".format(rotations[i,j,0],rotations[i,j,1],rotations[i,j,2]))
            f.write(" {:10.5f} {:10.5f} {:10.5f}".format(translations[i,0],translations[i,1],translations[i,2]))
            for ii in range(2):
                for jj in range(2):
                    f.write('{:14.7f}{:14.7f}'.format(SU2[i,ii,jj].real,SU2[i,ii,jj].imag))
            f.write('\n')
        f.write(str(len(newqpoints))+'\n')
        for i in range(len(newqpoints)):
            for j in range(3):
                f.write("{:14.7}".format(newqpoints[i][j]))
            f.write('\n')
        
        symk=newqpoints
        nsymk = len(newqpoints)
        for ik in range(nsymk):
            irreps, rotationslgroup, translationslgroup, mapping_little_group = get_spacegroup_irreps(
                basisvec, atmpos, atmnum, symk[ik]
            )
            
            #mapping_little_group, rotationslgroup, translationslgroup = findlg(rotations, translations, symk[ik])


            modes = utilsQE.readModesatKpoin(symk[ik],r'data/AgP2/Phonons/444/matdyn.modes', scffile=r'data/AgP2/Phonons/444/AgP2.scf.pwo')
            
            '''for i, m in enumerate(modes):
                for j, m2 in enumerate(modes):
                    print(i, j, np.dot(np.array(np.conjugate(m2)).flatten(),np.array(m).flatten()))'''

            lgtag = mapping_little_group
            f.write(str(len(lgtag))+'\n')
            for opr in lgtag:
                f.write(str(opr+1)+'  ')
            f.write('\n')
            
            for l, bi in enumerate(bandindex[ik]):
                f.write('{:5s}{:3s}{:12.6f}'.format(str(bi[0]+1),str(len(bi)),newbands[ik][bi[0]]))
                #print(bi)
                for iopr in range(len(lgtag)):

                    symeigval = symmetryeigval(np.array(modes),rotationslgroup[iopr], translationslgroup[iopr], basisvec[:], 
                                                    permutation_matrix(atmpos,rotationslgroup[iopr], translationslgroup[iopr], symk[ik]), bi)

                    f.write('{:10.5f}{:10.5f}'.format(np.real(symeigval) ,np.imag(symeigval)))
                f.write('\n')
        