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


def dynmat(alpha,beta,M1,M2,k):
    Dyn = np.array([
        [(alpha+beta)/M1, (alpha+beta*np.exp(-1j*k))/(M1*M2)**0.5],
        [(alpha+beta*np.exp(1j*k))/(M1*M2)**0.5, (alpha+beta)/M2]
    ])
    return Dyn

if __name__ == '__main__':
    alpha = 1
    beta = 1
    M1 = 1
    M2 = 1
    N = 4
    ls = np.arange(N)-1/2*N
    ks = np.pi/N*ls
    
    ms = [-0.1]#np.copy(ls)-0.25

    xcoords = np.array([-0.25,0.25])
    
    dummymodes = getModesforalq(alpha,beta,M1,M2,ks)
    
    modes = np.empty((len(ks), 2, len(ls), 2), dtype=np.complex128)

    for i,k in enumerate(ks):
        for a in range(2): #mode
            for j,cell in enumerate(ls):
                for b in range(2): #atom
                    modes[i,a,j,b] = dummymodes[i][b,a]*np.exp(1j*k*(cell)) 


    locmodes = LocalizedTransformastion(modes,ks,ms,ls)

    
    natoms = 2

    qpointindex = N//2
    modeindex = 0
    mcellindex = 0#N//2
    
    
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
    
    locnorm = [np.linalg.norm(vec) for vec in locpolvec.flatten()]
    locnormalize = np.max(locnorm)
    locpolvec /= locnormalize
    locnorm = [np.linalg.norm(vec) for vec in locpolvec.flatten()]
    locpolvec = locpolvec.flatten()


    ################################
    #Plotting
    ################################

    fig1 = plt.figure()
    ax = fig1.add_subplot()
    
    
    ax.scatter(latpos,np.zeros(len(latpos)),marker='o',linewidths=1)
    ax.quiver(latpos,np.zeros(len(latpos)), np.real(polvec),np.zeros(len(polvec)))
    ax.plot(latpos, norm)
    ax.set_title("Polarization vector, real part")
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    

    fig2 = plt.figure()
    ax2 = fig2.add_subplot()
    ax2.scatter(latpos,np.zeros(len(latpos)),marker='o',linewidths=1)
    ax2.quiver(latpos,np.zeros(len(latpos)), np.imag(polvec),np.zeros(len(polvec)), color='r')
    ax2.set_title("Polarization vector, imaginary part")
    ax2.set_xlabel('X')
    ax2.set_ylabel('Y')
    




    fig3 = plt.figure()
    ax3 = fig3.add_subplot()

    ax3.scatter(latpos,np.zeros(len(latpos)),marker='o',linewidths=1)
    ax3.quiver(latpos,np.zeros(len(latpos)), np.real(locpolvec),np.zeros(len(locpolvec)))
    ax3.plot(latpos, locnorm)
    ax3.set_title("Localized polarization vector, real part")
    ax3.set_xlabel('X')
    ax3.set_ylabel('Y')
    

    fig4 = plt.figure()
    ax4 = fig4.add_subplot()
    ax4.scatter(latpos,np.zeros(len(latpos)),marker='o',linewidths=1)
    ax4.quiver(latpos,np.zeros(len(latpos)), np.imag(locpolvec),np.zeros(len(locpolvec)),  color='r')
    ax4.set_title("Localized polarization vector, imaginary part")
    ax4.set_xlabel('X')
    ax4.set_ylabel('Y')
    

    
    plt.show()
    