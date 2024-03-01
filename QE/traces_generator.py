import numpy as np
import spglib as spg
import ase.io as io
import ase as ase
import ase.spacegroup as asespg
from spgrep import get_spacegroup_irreps
from spgrep.representation import get_character
import utilsQE


def band_indices(bandsatqpoint, tol=5e-3):
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


def symmetryeigval(modefunc,rs,ts):
    whatist=ts



if __name__ == '__main__':
    
    factor = 0.123983
    
    scf_file = r'data/AgP2/Phonons/AgP2.scf.pwi'
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
    
    qlabels, positions, qpoints = utilsQE.readHighSymPointsPhonon(r"data/AgP2/Phonons/matdyn.in", retKpoints=True)
    notqpoints, bands = utilsQE.readPhononbandFreq(r"data/AgP2/Phonons/AgP2.freq.gp")
    
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
    
    newbands = np.sort(np.transpose(bands[:,newpos]), axis=0)
    
    
    bandindex = []
    
    for i in range(len(newbands)):
        bandindex.append(band_indices(newbands[i]))
    
    
    newbands *= factor
    for i in range(len(newbands)):
        for j in range(len(newbands[i])):
            if abs(newbands[i][j]) < 5e-2:
                newbands[i][j] = 0.
    
    
    with open(r"DFTpycodes/QE/testTRACESC.txt",'w') as f:
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
        
        symk=newqpoints
        nsymk = len(newqpoints)
        for ik in range(nsymk):
            irreps, rotationslgroup, translationslgroup, mapping_little_group = get_spacegroup_irreps(
                basisvec, atmpos, atmnum, symk[ik]
            )
            modes = utilsQE.readModesatKpoin(symk[ik],r'data/AgP2/Phonons/matdyn.modes', scffile=r'data/AgP2/Phonons/AgP2.scf.pwo')
            
            lgtag = mapping_little_group
            f.write(str(len(lgtag))+'\n')
            for opr in lgtag:
                f.write(str(opr+1)+'  ')
            f.write('\n')
            #print(irreps)
            for l, bi in enumerate(bandindex[ik]):
                f.write('{:5s}{:3s}{:12.6f}'.format(str(bi[0]),str(len(bi)),newbands[ik][bi[0]-1]))
                #for iopr in range(len(lgtag)):
                    #print(len(irreps), iopr)
                    #print(np.trace(irreps[iopr])[0])
                    #f.write('{:10.5f}{:10.5f}'.format(np.real(np.trace(irreps[iopr])[0]) ,np.imag(np.trace(irreps[iopr])[0])))
                f.write('\n')
        