import numpy as np
import spglib as spg
import ase.io as io
import ase as ase
from spgrep import get_spacegroup_irreps,get_spacegroup_irreps_from_primitive_symmetry
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


def permutation_matrix(pos,r,t,q,sym_prec=7):
    
    num_atoms = len(pos)
    permu = np.zeros([num_atoms,num_atoms], dtype='complex')
    
    for j in range(num_atoms):
        new_pos=np.matmul(r,pos[j]) + t

        for j2 in range(num_atoms):
            diff = new_pos-pos[j2]
            if (abs(diff - np.rint(diff)) < 10**(-1*sym_prec)).all():
                phase = round(np.dot(q, new_pos-pos[j2]),3)
                
                permu[j2,j] = np.exp(-2j*np.pi*phase)

    return permu 


def similarity_transformation(U, mat):
    """ U x M x U^-1 """
    return np.matmul(U, np.matmul(mat, np.linalg.inv(U)))


def symmetryeigval(modefunc,r,t, latticevec, permatq, bandindices):
    U = np.transpose(latticevec)
    
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
    
    tracesfile = r'data/PbTe/444/PbTetraces.txt'

    scf_file = r'data/PbTe/444/PbTe.scf.pwi'
    scfoutfile = r'data/PbTe/444/PbTe.scf.pwo'

    matdyninfile = r"data/PbTe/444/tracesHSP/matdyn.in"
    frecgpfile = r"data/PbTe/444/tracesHSP/PbTe.freq.gp"
    modesfile = r'data/PbTe/444/tracesHSP/matdyn.modes'

    factor = 0.123983/1000 #cm^-1 to eV
    
    
    structure = io.read(scf_file)

    basisvec = structure.get_cell()
    atmpos= structure.get_scaled_positions()
    atmnum = structure.get_atomic_numbers()
    recbasis = structure.cell.reciprocal()
    
    

    refinedbasis, refinedatmpos,refinedatmnum=spg.standardize_cell(cell=(basisvec,atmpos,atmnum))
    sym_data = spg.get_symmetry_dataset((basisvec, atmpos, atmnum))

    #print(spg.get_symmetry_dataset((refinedbasis,refinedatmpos,refinedatmnum)).get('number'))

    fullrotations = sym_data.get("rotations")
    fulltranslations = sym_data.get("translations")
    
    tmat = sym_data.get('transformation_matrix')
    
    maskrot = np.arange(len(fullrotations))#np.unique(fullrotations,return_index=True, axis=0)
    maskrot = np.sort(maskrot)
    #print(rotations, maskrot)
    rotations = fullrotations[maskrot]
    translations = fulltranslations[maskrot]

    vectorsym = spg.get_symmetry_from_database(sym_data.get("hall_number"))
    vectorrot = vectorsym.get('rotations')
    vectortrans = vectorsym.get('translations')

    SG = sym_data.get("number")
    #print('Space group ',SG)
    nsym = len(rotations)
    SU2=np.zeros((nsym,2,2),dtype=np.complex128)
    for i in range(nsym):
        sopr_c=similarity_transformation(np.transpose(basisvec),rotations[i,:,:])
        SU2[i,:,:]=su2rot.get_su2rotation(sopr_c)

    
    
    qlabels, positions, qpoints = utilsQE.readHighSymPointsPhonon(matdyninfile, retKpoints=True)
    notqpoints, bands = utilsQE.readPhononbandFreq(frecgpfile)
    
    

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

    newbands *= factor
    
    bandindex = []
    
    for i in range(len(newbands)):
        #print(newbands)
        bandindex.append(band_indices(newbands[i],tol=2.5e-15))
        #stop = stop
    
    
    
    for i in range(len(newbands)):
        for j in range(len(newbands[i])):
            if abs(newbands[i][j]) < 5e-5:
                newbands[i][j] = 0.
    
    
    
    
    atmposcov = atmpos

    

    with open(tracesfile,'w') as f:
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
        
        nsymk = len(symk)
        
        for ik in range(nsymk):
            
            
            irreps, rotationslgroup, translationslgroup, dummymapping_little_group = get_spacegroup_irreps(
                refinedbasis[:,:], refinedatmpos, refinedatmnum, symk[ik]
            )
            
            dummyirreps, mapping_little_group = get_spacegroup_irreps_from_primitive_symmetry(rotations,translations,symk[ik])
            
            
            #print(mapping_little_group)
            #
            #print(findlg(rotations,translations,symk[ik]))
            #print(similarity_transformation(basisvec,rotationslgroup[-1]),rotations[-3])
            #print(dummymapping_little_group)
            #print(len(rotationslgroup))
            
            #mapping_little_group, rotationslgroup, translationslgroup = findlg(rotations, translations, symk[ik])
            
            

            lmaskrot = mapping_little_group
            lmaskrot = np.sort(lmaskrot)
            rotationslgroup = rotations[lmaskrot]
            translationslgroup = translations[lmaskrot]
            
            
            
            


            modes = utilsQE.readModesatKpoin(symk[ik],modesfile, scffile=scfoutfile)
            
            

            lgtag = mapping_little_group
            f.write(str(len(lgtag))+'\n')
            for opr in lgtag:
                f.write(str(opr+1)+'  ')
            f.write('\n')
            
            for l, bi in enumerate(bandindex[ik]):
                f.write('{:5s}{:3s}{:12.6f}'.format(str(bi[0]+1),str(len(bi)),newbands[ik][bi[0]]))
                #print(bi)
                for iopr in range(len(lgtag)):

                    symeigval = symmetryeigval(np.array(modes),rotationslgroup[iopr], translationslgroup[iopr], basisvec, 
                                                    permutation_matrix(atmpos,rotationslgroup[iopr], translationslgroup[iopr], symk[ik]), bi)

                    f.write('{:10.5f}{:10.5f}'.format(np.round(np.real(symeigval),10) ,np.round(np.imag(symeigval),10)))
                f.write('\n')
        