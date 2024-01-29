import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from readKgrid import fromOUTCARtoplot
from get_effectiveMass_FiniteDiff import secondDer_FiniteDiff
from numpy.polynomial.polynomial import polyfit, polyder, polyval

if __name__ == '__main__':
    #CONSTANTS
    eV2Hartree = 1/27.2113845#0.0367493

    qz = np.array([-4.9231649190557025-8.0899999999999999, 0. ,7.4862572210537408])
    qy = np.array([0,4.4800000000000004,0])
    qx = np.cross(qy, qz)

    qz *= 1/np.linalg.norm(qz)
    qx *= 1/np.linalg.norm(qx)
    qy *= 1/np.linalg.norm(qy)


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

    
    degree = 18

    polx = polyfit(qxvals,bandX,degree)
    poly = polyfit(qyvals,bandY,degree)
    polz = polyfit(qzvals,bandZ,degree)

    ddpolx = polyder(polx,2)
    ddpoly = polyder(poly,2)
    ddpolz = polyder(polz,2)


    fig1 = plt.figure()

    ax1 = plt.subplot(311)
    ax1.plot(qxvals, bandX, 'b')
    ax1.plot(qxvals, polyval(qxvals, polx),'--k')

    ax2 = plt.subplot(312)
    ax2.plot(qyvals, bandY, 'r')
    ax2.plot(qyvals, polyval(qyvals, poly),'--k')

    ax3 = plt.subplot(313)
    ax3.plot(qzvals, bandZ, 'g')
    ax3.plot(qzvals, polyval(qzvals, polz),'--k')

    ax2.set_ylabel('$E$ (meV)')
    ax3.set_xlabel('Position along direction in the first BZ')

    ax1.text(0.9,0.9,'$q_x$',size=20,ha='left', va='top', transform=ax1.transAxes)
    ax2.text(0.9,0.9,'$q_y$',size=20,ha='left', va='top', transform=ax2.transAxes)
    ax3.text(0.9,0.9,'$q_z$',size=20,ha='left', va='top', transform=ax3.transAxes)



    fig2 = plt.figure()

    ax1 = plt.subplot(311)
    ax1.plot(qxvals, polyval(qxvals, ddpolx),'b')

    ax2 = plt.subplot(312)
    ax2.plot(qyvals, polyval(qyvals, ddpoly),'r')

    ax3 = plt.subplot(313)
    ax3.plot(qzvals, polyval(qzvals, ddpolz),'g')

    ax1.hlines([0], np.min(qxvals), np.max(qxvals))
    ax2.hlines([0], np.min(qyvals), np.max(qyvals))
    ax3.hlines([0], np.min(qzvals), np.max(qzvals))

    ax2.set_ylabel('Second derivative, polynomial fit')
    ax3.set_xlabel('Position along direction in the first BZ')


    ax1.text(0.9,0.9,'$q_x$',size=20,ha='left', va='top', transform=ax1.transAxes)
    ax2.text(0.9,0.9,'$q_y$',size=20,ha='left', va='top', transform=ax2.transAxes)
    ax3.text(0.9,0.9,'$q_z$',size=20,ha='left', va='top', transform=ax3.transAxes)



    fig3 = plt.figure()

    ax1 = plt.subplot(311)
    ax1.plot(qxvals, polyval(qxvals, ddpolx)**-1/eV2Hartree*1000,'b')

    ax2 = plt.subplot(312)
    ax2.plot(qyvals, polyval(qyvals, ddpoly)**-1/eV2Hartree*1000,'r')

    ax3 = plt.subplot(313)
    ax3.plot(qzvals, polyval(qzvals, ddpolz)**-1/eV2Hartree*1000,'g')

    ax1.set_ylim(0,1)
    ax2.set_ylim(0,1)
    ax3.set_ylim(0,1)

    ax2.set_ylabel('Effective mass, polynomial fit')
    ax3.set_xlabel('Position along direction in the first BZ')



    ddx = secondDer_FiniteDiff(bandX,qxvals)
    ddy = secondDer_FiniteDiff(bandY,qyvals)
    ddz = secondDer_FiniteDiff(bandZ,qzvals)


    ax1.text(0.9,0.9,'$q_x$',size=20,ha='left', va='top', transform=ax1.transAxes)
    ax2.text(0.9,0.9,'$q_y$',size=20,ha='left', va='top', transform=ax2.transAxes)
    ax3.text(0.9,0.9,'$q_z$',size=20,ha='left', va='top', transform=ax3.transAxes)



    fig4 = plt.figure()

    ax1 = plt.subplot(311)
    ax1.plot(qxvals, ddx,'b')

    ax2 = plt.subplot(312)
    ax2.plot(qyvals, ddy,'r')

    ax3 = plt.subplot(313)
    ax3.plot(qzvals, ddz,'g')

    ax1.hlines([0], np.min(qxvals), np.max(qxvals))
    ax2.hlines([0], np.min(qyvals), np.max(qyvals))
    ax3.hlines([0], np.min(qzvals), np.max(qzvals))

    ax2.set_ylabel('Second derivative, finite differences')
    ax3.set_xlabel('Position along direction in the first BZ')

    ax1.text(0.9,0.9,'$q_x$',size=20,ha='left', va='top', transform=ax1.transAxes)
    ax2.text(0.9,0.9,'$q_y$',size=20,ha='left', va='top', transform=ax2.transAxes)
    ax3.text(0.9,0.9,'$q_z$',size=20,ha='left', va='top', transform=ax3.transAxes)



    fig5 = plt.figure()

    ax1 = plt.subplot(311)
    ax1.plot(qxvals, ddx**-1/eV2Hartree*1000,'b')

    ax2 = plt.subplot(312)
    ax2.plot(qyvals, ddy**-1/eV2Hartree*1000,'r')

    ax3 = plt.subplot(313)
    ax3.plot(qzvals, ddz**-1/eV2Hartree*1000,'g')
    ax3.set_ylim(0,8)

    ax2.set_ylabel('Effective mass, finite differences')
    ax3.set_xlabel('Position along direction in the first BZ')

    ax1.text(0.9,0.9,'$q_x$',size=20,ha='left', va='top', transform=ax1.transAxes)
    ax2.text(0.9,0.9,'$q_y$',size=20,ha='left', va='top', transform=ax2.transAxes)
    ax3.text(0.9,0.9,'$q_z$',size=20,ha='left', va='top', transform=ax3.transAxes)





    
    plt.show()