import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from readKgrid import fromOUTCARtoplot

if __name__ == '__main__':
    #CONSTANTS
    eV2Hartree = 1/27.2113845#0.0367493

    qz = np.array([-0.1495968000,0,0.1495966000])
    qy = np.array([0,1,0])
    qx = np.array([0.1495968000,0,0.1495966000])

    qz /= np.linalg.norm(qz)
    qx /= np.linalg.norm(qx)

    bands, efermi, dummyqxvals, dummyqyvals, dummyqzvals, occupation = fromOUTCARtoplot(outcarfile='Ag2Te\\ArchivosMartin\\GM_MBJ\\OUTCAR_MBJ',kpointsfile="Ag2Te\\ArchivosMartin\\GM_MBJ\\KPOINTS_MBJ_BS", qx=qx, qy=qy, qz=qz, weightfilter = 0)
    
    dummy = np.transpose(bands[0])
    dummyoc = np.transpose(occupation[0])

    index = np.sum(dummyoc, axis=0)

    band = ((dummy)[int(index[0])])
    bandXY = (band-np.min(band))*1000#*eV2Hartree

    ####################################
    qxvals = dummyqxvals[dummyqyvals == 0]
    qyvals = dummyqyvals[dummyqxvals == 0]
    
    bandX = bandXY[dummyqyvals == 0]
    bandY = bandXY[dummyqzvals == 0]
    ####################################

    bands, efermi, dummyqxvals, dummyqyvals, dummyqzvals, occupation = fromOUTCARtoplot(outcarfile='Ag2Te\\PGHR4\\OUTCAR',kpointsfile="Ag2Te\\PGHR4\\KPOINTS", qx=qx, qy=qy, qz=qz, weightfilter = 0)
    
    dummy = np.transpose(bands[0])
    dummyoc = np.transpose(occupation[0])

    index = np.sum(dummyoc, axis=0)

    band = ((dummy)[int(index[0])])
    bandYZ = (band-np.min(band))*1000#*eV2Hartree

    qzvals = dummyqzvals[dummyqyvals == 0]
    
    bandZ = bandYZ[dummyqyvals == 0]


    plt.figure()
    plt.plot(qxvals,bandX, 'b', label = '$\\Gamma \\rightarrow q_x$')
    plt.plot(qyvals,bandY, 'r', label = '$\\Gamma \\rightarrow q_y$')
    plt.plot(qzvals,bandZ, 'g', label = '$\\Gamma \\rightarrow q_z$')

    plt.ylabel('$E$ (meV)')
    plt.xlabel('Position along direction in the first BZ')

    plt.legend()
    plt.show()