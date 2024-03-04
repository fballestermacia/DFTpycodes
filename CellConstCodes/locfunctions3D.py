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
    
    modes = np.empty(len(qvecs), dtype='object')
    
    for i in range(len(qvecs)):
        dummy, pols = dynmat.DyagDinQ(iq=i)
        modes[i]= np.array(pols)

    return np.array(modes)
        




def LocalizedTransformastion(modesperq, qvecs, superatmpos, atmpos, reduced_cell, m,l):
    
    newpolvecs = np.empty((len(qvecs),3*len(atmpos),3*len(superatmpos)), dtype=np.complex128) 
    quasiVs = np.empty((len(qvecs),3*len(atmpos),3*len(superatmpos)), dtype=np.complex128) 
    # INDICES ARE: q-point, mode, kappa_{x,y,z}

    #TODO: MAYBE CHANGE IT SO WE ONLY CALCULATE 1 MODE TO SAVE MEMORY
    #print(l)
    for i,q in enumerate(qvecs): #qpoint
        print('Current q:', i, '/', len(qvecs))
        for a in range(3*len(atmpos)): #mode
            for b in range(len(atmpos)): #atom
                for lz in range(l[2]):
                    for ly in range(l[1]):
                        for lx in range(l[0]):
                            #print(i,a,b,np.dot(q,[reduced_cell[i]*(np.dot(atmpos[b],reduced_cell[i])//1-m[i]) for i in range(3)]))
                            ldummy = [lx,ly,lz]
                            indexX=3*len(atmpos)*lx 
                            indexY=3*len(atmpos)*l[0]*ly
                            indexZ=3*len(atmpos)*l[1]*l[0]*lz
                            indextot = indexX+indexY+indexZ
                            #print(np.subtract(ldummy,m))
                            
                            #print(a+indextot, '/', 3*len(superatmpos))
                            newpolvecs[i,3*b:3*b+3,a+indextot] = modesperq[i][3*b:3*b+3,a]*np.exp(1j*np.dot(q,np.sum(reduced_cell*ldummy,axis=0)))
                            quasiVs[i,3*b:3*b+3,a+indextot] = modesperq[i][3*b:3*b+3,a]*np.exp(1j*np.dot(q,np.sum(reduced_cell*np.subtract(ldummy,m),axis=0)))
        #print('here',reduced_cell,ldummy,'\n 2here',reduced_cell*ldummy,'\n 3here',np.sum(reduced_cell*ldummy,axis=0),'\n 4here',np.dot(q,np.sum(reduced_cell*ldummy,axis=0)))

    '''for i in range(len(modesperq)):
        phasesatq =np.empty(len(atmpos), dtype=np.complex128)
        for j in range(len(atmpos)):
            #print(np.dot(qvecs[i],atmpos[j]))
            phasesatq[j] = np.exp(2j*np.pi*np.dot(qvecs[i],atmpos[j]))
        
        
        phases3D = np.kron(phasesatq,np.ones(3))


        

        modesperq[i] = np.transpose(np.multiply(np.transpose(modesperq[i]), phases3D)) #THESE ARE THE u VECTORS IN THE NOTES ???
        
        for k in range(len(atmpos)):
            quasiVs[i,:,:,k] = modesperq[i]*np.conjugate(phases3D[3*k])
    '''


    
    ##############
    #TODO: NORMALIZE THE VECTORS 
    ##############

    

    locmodes = np.sum(quasiVs,axis=0)
    #print(locmodes)
    return newpolvecs, locmodes
        

def berrypposop(q, locmode, positions):
    locmode = np.reshape(locmode,(len(positions),3))
    exponential = [np.exp(1j*np.dot(q,pos)) for pos in positions]
    
    center = np.empty(3)
    for i in range(3):
        matelement = np.dot(np.conjugate(locmode[:,i]),np.multiply(exponential,locmode[:,i]))
        center[i] = 1/np.linalg.norm(q)*np.imag(np.log(matelement))
    return center    




if __name__ == '__main__':
    

    inittime = datetime.now()
    

    dyn = CC.Phonons.Phonons()
    dyn.LoadFromQE(fildyn_prefix="data/AgP2/Phonons/dynmats/AgP2.dyn", nqirr=30)
    dyn = dyn.InterpolateMesh([1,1,10])
    #print(dyn.GetSupercell())
    super_dyn = dyn.GenerateSupercellDyn(dyn.GetSupercell())
    #print(super_dyn.q_tot)
    #print(dyn.q_tot[63])
    #print(dyn.structure.coords)

    #freqs, pols = super_dyn.DiagonalizeSupercell(verbose=True)
    #freqs2, pols2 = dyn.DyagDinQ(iq=60)
    ####
    # pols -> polarization vectors, shape (3*N_atoms, N_modes)
    ####
    
    #print(np.shape(pols2[:,0]))
    #print(pols2[:,0])

    dyn2 = super_dyn#dyn.InterpolateMesh([2,2,2])
    #dyn2.save_qe('testmatdyn')
    #print(dyn.structure.unit_cell)
    
    #print(np.shape(dyn2.structure.coords))
    
    
    #let's do only one m to save memory
    l = dyn.GetSupercell()
    mx, my, mz = 0,0,0 #np.arange(4),np.arange(4),np.arange(4)
    
    #print(dyn.structure.unit_cell)
    #print(dyn2.structure.unit_cell)

    
    modes = getModesforalq(dyn)


    xcoords = np.zeros(np.shape(dyn2.structure.coords))
    for i in range(dyn2.structure.N_atoms):
        xcoords[i,:] = CC.Methods.covariant_coordinates(dyn2.structure.unit_cell, dyn2.structure.coords[i,:])
    
    NOSCxcoords = np.zeros(np.shape(dyn.structure.coords))
    for i in range(dyn.structure.N_atoms):
        NOSCxcoords[i,:] = CC.Methods.covariant_coordinates(dyn.structure.unit_cell, dyn.structure.coords[i,:])

    
    qvecs = np.zeros(np.shape(dyn.q_tot))
    for i in range(len(dyn.q_tot)):
        qvecs[i,:] = CC.Methods.covariant_coordinates(dyn.structure.get_reciprocal_vectors(), dyn.q_tot[i])
    print(dyn.structure.get_reciprocal_vectors())

    modes, locmodes = LocalizedTransformastion(modes,qvecs,xcoords,NOSCxcoords,dyn.structure.unit_cell, [mx,my,mz],l)
    #print(np.shape(modes), np.shape(modes[0]))
    #print(np.shape(locmodes), np.shape(locmodes[0]))

    finishtime = datetime.now()
    print(finishtime-inittime)

    modeindex = 0
    print(qvecs[-1])

    center = berrypposop(qvecs[-1],locmodes[modeindex],xcoords) #TODO: Check this function
    print(center)

    xcenter = center#CC.Methods.covariant_coordinates(dyn2.structure.unit_cell, center)

    atypes = [x-1 for x in dyn2.structure.get_atomic_types()]
    
    col = ['g','k']
    colors = [col[i] for i in atypes]

    natoms = len(atypes)
    
    latpos = xcoords

    polvec = modes[-1][modeindex]
    polvec = np.reshape(polvec,(natoms,3))
    
    normalize = np.max(np.linalg.norm(polvec,axis=0))
    polvec /= normalize

    locpolvec = locmodes[modeindex]
    locpolvec = np.reshape(locpolvec,(natoms,3))
    

    locnormalize = np.max(np.linalg.norm(locpolvec,axis=0))
    locpolvec /= locnormalize



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

    fig4 = plt.figure()
    ax4 = fig4.add_subplot(projection='3d')
    ax4.scatter(latpos[:,0],latpos[:,1],latpos[:,2],marker='o',linewidths=1, c=colors)
    ax4.quiver(latpos[:,0],latpos[:,1],latpos[:,2], np.imag(locpolvec[:,0]),np.imag(locpolvec[:,1]),np.imag(locpolvec[:,2]), length=0.5, color='r')
    ax4.set_title("Localized polarization vector, imaginary part")
    ax4.set_xlabel('X')
    ax4.set_ylabel('Y')
    ax4.set_zlabel('Z')

    
    plt.show()

