import cellconstructor as CC
import cellconstructor.Phonons
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def getModesforalq(alpha,beta,M1,M2,ks):
    
    # Returns array of polarization vectors for each qvec in the
    # shape (N_Qpoints,N_atoms, N_modes)

    
    modes = np.empty(len(ks), dtype='object')
    
    for i in range(len(ks)):
        dummy, pols = np.linalg.eig(dynmat(alpha,beta,M1,M2,ks[i]))
        modes[i]= np.array(pols)
        
    return np.array(modes)
        

def LocalizedTransformastion(modesperq,ks,ms,ls):
    
    quasiVs = np.empty((len(ks),2,len(ms),len(ls),2),dtype=np.complex128) 
    # INDICES ARE: q-point, j, m, l, kappa,

    for i,k in enumerate(ks): #qpoint
        for a in range(2): # mode
            for j,mcell in enumerate(ms): #m vector
                for l in range(len(ls)): #l vector
                    for b in range(2): #atom
                        quasiVs[i,a,j,l,b] = modesperq[i,a,l,b]*np.exp(-1j*k*mcell)

    
    ##############
    #TODO: NORMALIZE THE VECTORS 
    ##############

    

    locmodes = np.sum(quasiVs,axis=0)
    #locmodes = locmodes*np.exp(-1j*np.angle(locmodes))
    #print(locmodes)
    return locmodes


def DelocalizeTransformation(modesperRs, ks, ms, ls):
    quasiUs = np.empty((len(ks),2,len(ms),len(ls),2),dtype=np.complex128)
    quasiVs = np.empty((len(ks),2,len(ms),len(ls),2),dtype=np.complex128) 
    # INDICES ARE: q-point, j, m, l, kappa,

    for i,k in enumerate(ks): #qpoint
        for a in range(2): # mode
            for j,mcell in enumerate(ms): #m vector
                for l in range(len(ls)): #l vector
                    for b in range(2): #atom
                        quasiUs[i,a,j,l,b] = modesperRs[a,j,l,b]*np.exp(1j*k*mcell)

    
    ##############
    #TODO: NORMALIZE THE VECTORS 
    ##############

    

    delocmodes = np.sum(quasiUs,axis=2)
    #locmodes = locmodes*np.exp(-1j*np.angle(locmodes))
    #print(locmodes)
    return delocmodes


def dynmat(alpha,beta,M1,M2,k):
    Dyn = np.array([
        [(alpha+beta)/M1, (alpha+beta*np.exp(-1j*k))/(M1*M2)**0.5],
        [(alpha+beta*np.exp(1j*k))/(M1*M2)**0.5, (alpha+beta)/M2]
    ])
    return Dyn


def berrypposop(q, lmode, positions):

    #################################
    #NOTE: q needs to be sufficiently small, otherwise terms of order O(q\sigma) are significant and the operator does not return the center of the localized function
    # the bigger the system, the smaller the q vector must be
    #################################

    exponential = np.exp(1j*q*positions)

    matelement = np.dot(np.conjugate(lmode),np.multiply(exponential,lmode))
    center = 1/q*np.imag(np.log(matelement))
    return center
    



if __name__ == '__main__':
    alpha = 1
    beta = 1
    M1 = 2
    M2 = 1
    N = 50
    ls = np.arange(N)-1/2*N
    ks = np.pi/N*ls
    


    ms = np.copy(ls)-0.25

    xcoords = np.array([-0.25,0.25])
    
    dummymodes = getModesforalq(alpha,beta,M1,M2,ks)
    
    modes = np.empty((len(ks), 2, len(ls), 2), dtype=np.complex128)

    for i,k in enumerate(ks):
        for a in range(2): #mode
            for j,cell in enumerate(ls):
                for b in range(2): #atom
                    modes[i,a,j,b] = dummymodes[i][b,a]*np.exp(1j*k*(cell)) 

    if M1 == M2:
        dummycopy = modes.copy()
        phase0a = np.exp(1j*np.angle(dummycopy[:,0,N//2,0] + dummycopy [:,1,N//2,0]))
        phase0b = np.exp(1j*np.angle(dummycopy[:,0,N//2,1] + dummycopy [:,1,N//2,1]))
        phase1a = np.exp(1j*np.angle(dummycopy[:,0,N//2,0] - dummycopy [:,1,N//2,0]))
        phase1b = np.exp(1j*np.angle(dummycopy[:,0,N//2,1] - dummycopy [:,1,N//2,1]))
        modes[:,0,:,0] *= phase0a 
        modes[:,0,:,1] *= phase0b 
        modes[:,1,:,0] *= phase1a 
        modes[:,1,:,1] *= phase1b 
        

    locmodes = LocalizedTransformastion(modes,ks,ms,ls)
    nolocmodes = DelocalizeTransformation(locmodes,ks,ms,ls)
    
    natoms = 2

    qpointindex = N//2 #QPOINT MUST BE CLOSE TO GAMMA otherwise O(q\sigma) becomes significant
    
    modeindex = 1
    mcellindex = N//2-10
    
    
    latpos = np.kron(np.ones(len(ls)),xcoords)+np.kron(ls,np.array([1,1]))
    
    

    polvec = modes[qpointindex,modeindex]
    polvec = np.reshape(polvec,(N,natoms))
    
    norm = [np.linalg.norm(vec) for vec in polvec.flatten()]
    normalize = np.max(norm)
    polvec /= normalize
    norm = [np.linalg.norm(vec) for vec in polvec.flatten()]
    polvec = polvec.flatten()
    
    

    locpolvec = locmodes[modeindex,mcellindex,:,:]
    locpolvec = np.reshape(locpolvec,(N,natoms))
    
    
    locnorm = [np.linalg.norm(vec) for vec in locpolvec.flatten()] #i think this is wrong but whatever
    locnormalize = np.max(locnorm)
    locpolvec /= locnormalize
    locnorm = [np.linalg.norm(vec) for vec in locpolvec.flatten()]
    locpolvec = locpolvec.flatten()

    

    pos =  berrypposop(1e-12, locmodes[modeindex,mcellindex].flatten(),latpos[:])
    print(pos)


    
    polvec2 = nolocmodes[qpointindex,modeindex]
    polvec2 = np.reshape(polvec2,(N,natoms))

    norm2 = [np.linalg.norm(vec) for vec in polvec2.flatten()]
    normalize = np.max(norm2)
    polvec2 /= normalize
    norm2 = [np.linalg.norm(vec) for vec in polvec2.flatten()]
    polvec2 = polvec2.flatten()
    ################################
    #Plotting
    ################################

    fig1 = plt.figure()
    ax = fig1.add_subplot()
    
    
    ax.scatter(latpos,np.zeros(len(latpos)),marker='o',linewidths=1)
    ax.quiver(latpos,np.zeros(len(latpos)), np.zeros(len(polvec)),np.real(polvec))
    ax.plot(latpos, norm)
    ax.set_title("Polarization vector, real part")
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    

    fig2 = plt.figure()
    ax2 = fig2.add_subplot()
    ax2.scatter(latpos,np.zeros(len(latpos)),marker='o',linewidths=1)
    ax2.quiver(latpos,np.zeros(len(latpos)), np.zeros(len(polvec)),np.imag(polvec), color='r')
    ax2.set_title("Polarization vector, imaginary part")
    ax2.set_xlabel('X')
    ax2.set_ylabel('Y')
    




    fig3 = plt.figure()
    ax3 = fig3.add_subplot()

    ax3.scatter(latpos,np.zeros(len(latpos)),marker='o',linewidths=1)
    ax3.quiver(latpos,np.zeros(len(latpos)), np.zeros(len(locpolvec)),np.real(locpolvec))
    ax3.plot(latpos, locnorm)
    ax3.vlines(pos, 0,1, 'r')
    ax3.set_title("Localized polarization vector, real part")
    ax3.set_xlabel('X')
    ax3.set_ylabel('Y')
    

    fig4 = plt.figure()
    ax4 = fig4.add_subplot()
    ax4.scatter(latpos,np.zeros(len(latpos)),marker='o',linewidths=1)
    ax4.quiver(latpos,np.zeros(len(latpos)), np.zeros(len(locpolvec)), np.imag(locpolvec),  color='r')
    ax4.set_title("Localized polarization vector, imaginary part")
    ax4.set_xlabel('X')
    ax4.set_ylabel('Y')
    
    fig5 = plt.figure()
    ax5 = fig5.add_subplot()
    ax5.scatter(latpos,np.zeros(len(latpos)),c='k',marker='o',linewidths=1)
    ax5.plot(latpos, locnorm)
    ax5.vlines(pos, 0,1, 'r')
    ax5.set_title("Localized polarization vector, magnitude")
    ax5.set_xlabel('X')
    ax3.set_ylabel('Y')

    fig6 = plt.figure()
    ax6 = fig6.add_subplot()
    
    
    ax6.scatter(latpos,np.zeros(len(latpos)),marker='o',linewidths=1)
    ax6.quiver(latpos,np.zeros(len(latpos)), np.zeros(len(polvec2)),np.real(polvec2))
    ax6.plot(latpos, norm2)
    ax6.set_title("Polarization vector, real part")
    ax6.set_xlabel('X')
    ax6.set_ylabel('Y')
    
    plt.show()
    