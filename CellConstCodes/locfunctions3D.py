import cellconstructor as CC
import cellconstructor.Phonons
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from datetime import datetime, timedelta

def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()



def getModesforalq(dynmat):
    
    # Returns array of polarization vectors for each qvec in the
    # shape (N_Qpoints,3*N_atoms, N_modes)
    qvecs = dynmat.q_tot
    
    modes = []
    
    for i in range(len(qvecs)):
        dummy, pols = dynmat.DyagDinQ(iq=i)
        modes.append(np.array(pols))

    return np.array(modes)
        

def findposinunitcell(modq, l, superatmpos, atmpos, reduced_cell):
    for lx in range(l[0]):
        for ly in range(l[1]):
            for lz in range(l[2]):
                translate = np.sum(reduced_cell*[lx,ly,lz], axis=1)
                if (np.abs(np.sum(np.abs(np.subtract(atmpos,superatmpos-translate)),axis=1))<1e-5).any():
                    ldummy = [lx,ly,lz]
                    mod = modq[np.abs(np.sum(np.abs(np.subtract(atmpos,superatmpos-translate)),axis=1))<1e-5]
                    return mod, ldummy
    return None, None    

def LocalizedTransformastion(modesperq, qvecs, superatmpos, atmpos, reduced_cell, m,l):
    
    newpolvecs = np.empty((len(qvecs),len(superatmpos),3), dtype=np.complex128) 
    quasiVs = np.empty((len(qvecs),len(superatmpos),3), dtype=np.complex128) 
    # INDICES ARE: q-point, atom, {x,y,z}
    coloretes = np.empty(len(superatmpos))
    for i,q in enumerate(qvecs): #qpoint
        #print('Current q:', i, '/', len(qvecs))
        printProgressBar(i, len(qvecs))
        modq = np.reshape(modesperq[i],(len(atmpos),3))
        for b in range(len(superatmpos)): #atom
            

            mod, ldummy = findposinunitcell(modq, l, superatmpos[b], atmpos, reduced_cell)
            #coloretes[b] = ldummy[2]
            #print(ldummy,np.subtract(ldummy,m))

            newpolvecs[i,b] = mod*np.exp(1j*np.dot(q,superatmpos[b]))
            #quasiVs[i,b] = mod*np.exp(1j*np.dot(q,np.subtract(superatmpos[b],m)))
            quasiVs[i,b] = newpolvecs[i,b]*np.exp(-1j*np.dot(q,m))

            #newpolvecs[i,b] = mod*np.exp(1j*np.dot(q,np.sum(reduced_cell*ldummy,axis=0)))
            #quasiVs[i,b] = mod*np.exp(1j*np.dot(q,np.sum(reduced_cell*np.subtract(ldummy,m),axis=0)))
    
    for i in range(len(qvecs)):
        newpolvecs[i] /= np.linalg.norm(newpolvecs[i])
    
    printProgressBar(1,1)
    locmodes = np.sum(quasiVs,axis=0)
    
    locmodes /= np.linalg.norm(locmodes)

    return newpolvecs, locmodes #, coloretes
        

def berrypposop(qvec, locmode, positions):
    center = np.empty(3)
    for i in range(3):
        exponential = np.array([np.exp(1j*np.dot(qvec,pos[i])) for pos in positions])
        #print(exponential)
        #print(np.kron(exponential))
        matelement = np.dot(np.conjugate(locmode.flatten()),np.multiply(exponential.flatten(),locmode.flatten()))
        center[i] = 1/np.linalg.norm(qvec)*np.imag(np.log(matelement))
    return center    


def ModAlongDirections(mode, positions, dirs=np.eye(3)):
    dirx, diry, dirz = dirs
    xvals = np.empty(len(positions))
    yvals = np.empty(len(positions))
    zvals = np.empty(len(positions))
    
    xmods = np.empty(len(positions))
    ymods = np.empty(len(positions))
    zmods = np.empty(len(positions))

    totmods = np.empty(len(positions))
    totvals = np.empty(len(positions))

    for i in range(len(positions)):
        xvals[i] = np.dot(positions[i], dirx)/np.linalg.norm(dirx)
        yvals[i] = np.dot(positions[i], diry)/np.linalg.norm(diry)
        zvals[i] = np.dot(positions[i], dirz)/np.linalg.norm(dirz)

        xmods[i] = np.linalg.norm(mode[i])
        ymods[i] = np.linalg.norm(mode[i])
        zmods[i] = np.linalg.norm(mode[i])

        totmods[i] = np.linalg.norm(mode[i])
        totvals[i] = np.linalg.norm(positions[i])

    return np.array([xvals,yvals,zvals]), np.array([xmods,ymods,zmods]), totvals, totmods






if __name__ == '__main__':
    

    inittime = datetime.now()
    

    dyn = CC.Phonons.Phonons()
    dyn.LoadFromQE(fildyn_prefix="data/AgP2/Phonons/dynmats/AgP2.dyn", nqirr=30)
    
    #dyn = dyn.InterpolateMesh([6,6,6])


    super_dyn = dyn.GenerateSupercellDyn(dyn.GetSupercell())


    dyn2 = super_dyn#dyn.InterpolateMesh([2,2,2])

    #let's do only one m to save memory
    l = dyn.GetSupercell()
    
    ##########################################
    #IMPORTANT NOTE: i think CellConstructor gives the positions in reversed order with respect to the y coordinate. 
    #Thus, the center must have a -1 in the y direction so as to center the function in the ''first'' unit cell 
    ##########################################
    mx, my, mz = np.sum(np.multiply([0,-1,0],[dyn.structure.get_reciprocal_vectors(),dyn.structure.get_reciprocal_vectors(),dyn.structure.get_reciprocal_vectors()]),axis=0)
    
    #print(dyn.structure.unit_cell)
    #print(dyn2.structure.unit_cell)

    
    modes = getModesforalq(dyn)

    qindex = 5
    modeindex = 5
    modes = modes[:,:,modeindex]
    
    
    '''xcoords = np.zeros(np.shape(dyn2.structure.coords))
    for i in range(dyn2.structure.N_atoms):
        xcoords[i,:] = CC.Methods.covariant_coordinates(dyn2.structure.unit_cell, dyn2.structure.coords[i,:])
    
    NOSCxcoords = np.zeros(np.shape(dyn.structure.coords))
    for i in range(dyn.structure.N_atoms):
        NOSCxcoords[i,:] = CC.Methods.covariant_coordinates(dyn.structure.unit_cell, dyn.structure.coords[i,:])'''

    
    
    qvecs = np.zeros(np.shape(dyn.q_tot))
    for i in range(len(dyn.q_tot)):
        qvecs[i,:] = CC.Methods.covariant_coordinates(dyn.structure.get_reciprocal_vectors(), dyn.q_tot[i])
    #print(dyn.structure.get_reciprocal_vectors())
    #qvecs = dyn.q_tot
    modes, locmodes = LocalizedTransformastion(modes,qvecs,dyn2.structure.coords,dyn.structure.coords, np.transpose(dyn.structure.unit_cell), [mx,my,mz],l)
    #print(np.shape(modes), np.shape(modes[0]))
    #print(np.shape(locmodes), np.shape(locmodes[0]))

    finishtime = datetime.now()
    print(finishtime-inittime)

    

    q = 1e-12
    center = berrypposop([q,q,q],locmodes,dyn2.structure.coords) #TODO: Check this function
    

    rvals, rmods, totvals, totmods = ModAlongDirections(locmodes, dyn2.structure.coords, dirs=np.eye(3))
    
    
    xcenter = center#CC.Methods.covariant_coordinates(dyn2.structure.unit_cell, center)
    
    #xcenter = [-1,17.5,2.5]

    #xcenter[0] = np.min(dyn2.structure.coords[:,0])+xcenter[0]
    #xcenter[2] = np.min(dyn2.structure.coords[:,2])+xcenter[2]

    print(xcenter)

    atypes = [x-1 for x in dyn2.structure.get_atomic_types()]
    
    col = ['g','k']
    colors = [col[i] for i in atypes]  #[int(c) for c in coloretes] 

    natoms = len(atypes)
    
    latpos = dyn2.structure.coords 

    normfactor = 40

    polvec = modes[qindex]
    polvec *= normfactor
    '''polvec = np.reshape(polvec,(natoms,3))'''
    
    '''normalize = np.max(np.linalg.norm(polvec,axis=0))
    polvec /= normalize'''

    locpolvec = locmodes
    locpolvec *= normfactor
    '''locpolvec = np.reshape(locpolvec,(natoms,3))'''
    

    '''locnormalize = np.max(np.linalg.norm(locpolvec,axis=0))
    locpolvec /= locnormalize'''



    ################################
    #Plotting
    ################################


    fig1 = plt.figure()
    ax = fig1.add_subplot(projection='3d')
    
    ax.scatter(latpos[:,0],latpos[:,1],latpos[:,2],marker='o',linewidths=1, c=colors)
    ax.quiver(latpos[:,0],latpos[:,1],latpos[:,2], polvec[:,0],polvec[:,1],polvec[:,2], length=0.5)
    ax.set_title("Polarization vector, real part")
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(projection='3d')
    ax2.scatter(latpos[:,0],latpos[:,1],latpos[:,2],marker='o',linewidths=1, c=colors)
    ax2.quiver(latpos[:,0],latpos[:,1],latpos[:,2], np.imag(polvec[:,0]),np.imag(polvec[:,1]),np.imag(polvec[:,2]), length=0.5, color='r')
    ax2.set_title("Polarization vector, imaginary part")
    ax2.set_xlabel('X')
    ax2.set_ylabel('Y')
    ax2.set_zlabel('Z')




    fig3 = plt.figure()
    ax3 = fig3.add_subplot(projection='3d')

    ax3.scatter(latpos[:,0],latpos[:,1],latpos[:,2],marker='o',linewidths=1, c=colors)
    ax3.quiver(latpos[:,0],latpos[:,1],latpos[:,2], locpolvec[:,0],locpolvec[:,1],locpolvec[:,2], length=0.5)
    ax3.scatter(xcenter[0],xcenter[1],xcenter[2],s=20,c='r')
    ax3.set_title("Localized polarization vector, real part")
    ax3.set_xlabel('X')
    ax3.set_ylabel('Y')
    ax3.set_zlabel('Z')
    #ax3.invert_yaxis()

    fig4 = plt.figure()
    ax4 = fig4.add_subplot(projection='3d')
    ax4.scatter(latpos[:,0],latpos[:,1],latpos[:,2],marker='o',linewidths=1, c=colors)
    ax4.quiver(latpos[:,0],latpos[:,1],latpos[:,2], np.imag(locpolvec[:,0]),np.imag(locpolvec[:,1]),np.imag(locpolvec[:,2]), length=0.5, color='r')
    ax4.set_title("Localized polarization vector, imaginary part")
    ax4.set_xlabel('X')
    ax4.set_ylabel('Y')
    ax4.set_zlabel('Z')


    fig5 = plt.figure()
    ax5 = fig5.add_subplot()
    for i in range(3):
        ax5.scatter(rvals[i], rmods[i], c=['C0','C1','C2'][i])
    
    ax5.vlines(xcenter,0,np.max(rmods), colors=['C0','C1','C2'], linestyles='--')

    fig6 = plt.figure()
    ax6 = fig6.add_subplot()
    ax6.scatter(totvals, totmods, c='k')
    ax6.vlines(np.linalg.norm(xcenter),0,np.max(totmods), colors='r', linestyles='--')

    
    '''fig7 = plt.figure()
    ax7 = fig7.add_subplot(projection='3d')
    ax7.scatter(latpos[0,0],latpos[0,1],latpos[0,2], c='r')
    ax7.scatter(latpos[-1,0],latpos[-1,1],latpos[-1,2], c='b')'''

    plt.show()

