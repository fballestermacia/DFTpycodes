import numpy as np
import matplotlib.pyplot as plt
import utilsVASP



if __name__ == '__main__':
    qz = np.array([-13.01,0,7.48])
    #qz = np.array([-0.1495968, 0., 0.1495966])
    qz = qz/np.linalg.norm(qz)

    qy = np.array([0,1.,0])

    qx = np.cross(qy, qz)#np.cross(qy,qz/np.linalg.norm(qz))

    qz *= 1/np.linalg.norm(qz)
    qx *= 1/np.linalg.norm(qx)
    qy *= 1/np.linalg.norm(qy)

    npoints = 10

    kpath = utilsVASP.fromgridtoKPOINTS(qx,qz, [-0.2,0.2], [-0.2,0.2], 51,71)#np.array(Kpath([-qz, np.array([0.,0.,0.]), qz],npoints))
    
    
    with open('Ag2Te\MBJ_thirdtry\IBZKPT') as f:
        ibzkpt = f.readlines()
    
    nkibz = int(ibzkpt[1])
    
    with open('Ag2Te\\PGHR_XZ\\KPOINTS','w') as f:
        f.write('Kpoints generated automatically by F.B.\n')
        f.write('{}\n'.format(len(kpath)+nkibz))
        f.write('Reciprocal lattice\n')
        for i in range(3,len(ibzkpt)):
            f.write(ibzkpt[i])
        for i in range(len(kpath)):
                #f.write('{:.5f} {:.5f} {:.5f}\n'.format(kpath[i,0],kpath[i,1],kpath[i,2]))
                f.write('{:.14f} {:.14f} {:.14f} {}\n'.format(kpath[i,0],kpath[i,1],kpath[i,2],'0.00'))
    
    '''grid = fromgridtoKPOINTS(qx,qy,[-1,1],[-1,1],10,30)

    fig = plt.figure()
    fig.set_size_inches(15, 8)
    ax = fig.add_subplot(projection='3d')
    ax.scatter(grid[:,0], grid[:,1], grid[:,2])
    plt.show()'''

    '''bands, efermi, qxvals, qyvals, occupation = utilsVASP.fromOUTCARtoplot(outcarfile='Ag2Te\PerpGrid\OUTCAR',kpointsfile="Ag2Te\PerpGrid\KPOINTS", qx=qz, qy=qy, weightfilter = 0)

    dummy = np.transpose(bands[0])
    dummyoc = np.transpose(occupation[0])

    index = np.sum(dummyoc, axis=0)



    band = ((dummy)[int(index[0])])*1000#-efermi)*1000

    #band -= np.min(band)
    
    #band = np.where(band>400, 400, band)
    
    nlev = 30

    fig = plt.figure()
    fig.set_size_inches(15, 8)
    ax = plt.subplot(111)
    drawing = ax.tricontourf(qxvals,qyvals,band, levels = nlev)
    colorbar = fig.colorbar(drawing, ax = ax)
    ax.tricontour(qxvals,qyvals,band, colors='k', levels = nlev)
    ax.scatter(qxvals,qyvals,s=3,c='k')
    ax.scatter(0,0,s=10,c='r')

    fermilevel = ax.tricontour(qxvals, qyvals, band, [0], colors='r', zorder=10)'''


    plt.show()

