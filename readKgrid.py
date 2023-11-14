import numpy as np
import matplotlib.pyplot as plt

def fromKPOINTStogrid(kpointsfile='KPOINTS', qx = np.array([1,0,0]), qy = np.array([0,1,0]), weightfilter = None):
    
    kp = open(kpointsfile).readlines()
    
    RECLAT = False
    if kp[2][0].upper() == 'L':
        startline = 4
    elif kp[2][0].upper() == 'R':
        startline = 3
        RECLAT = True
    else: startline = 3
    
    qxvals = []
    qyvals = []
    
    if weightfilter is not None:
        readornot = []
    else:
        readornot = None
    
    for i in range(startline, len(kp)):
        if weightfilter is not None:
            if float(kp[i].strip().split()[3]) == weightfilter:
                readornot.append(True)
            else: 
                readornot.append(False)
                continue
                
            
        kvalx, kvaly, kvalz = [float(x) for x in kp[i].strip().split()[:3]]
        
        if not RECLAT:
            coefx = np.dot(np.array([kvalx,kvaly,kvalz]),qx)/np.linalg.norm(qx)**2
            coefy = np.dot(np.array([kvalx,kvaly,kvalz]),qy)/np.linalg.norm(qy)**2
        else: 
            coefx = kvalx
            coefy = kvaly
        
        qxvals.append(coefx)
        qyvals.append(coefy)
    
    return np.array(qxvals), np.array(qyvals), readornot


def fromgridtoKPOINTS(qx, qy, qxlims, qylims, nxpoints, nypoints, origin = np.array([0,0,0])):
    NX = np.linspace(qxlims[0], qxlims[1], nxpoints)
    NY = np.linspace(qylims[0], qylims[1], nypoints)
    
    QX, QY = np.meshgrid(NX,NY)
    
    KX = np.outer(QX.flatten(),qx) + origin
    KY = np.outer(QY.flatten(),qy) + origin

    KGRID = KX + KY
    
    return KGRID


def fromOUTCARtoplot(outcarfile = 'OUTCAR', kpointsfile = 'KPOINTS', qx = np.array([1,0,0]), qy = np.array([0,1,0]), weightfilter = None):
    
    qxvals, qyvals, readornot = fromKPOINTStogrid(kpointsfile=kpointsfile, qx = qx, qy = qy, weightfilter=weightfilter)
    
    # Separate outcar file in lines and remove trailing and heading whitespaces
    outcar = [line for line in open(outcarfile) if line.strip()] 
    
    for i, line in enumerate(outcar):
        # Check for needed constant values and parameters
        if 'NKPTS =' in line:
            nkpts = int(line.split()[3])
            nband = int(line.split()[-1])
            
        if 'ISPIN  =' in line:
            ispin = int(line.split()[2])
        
        if "k-points in reciprocal lattice and weights" in line:
            Lvkpts = i + 1

        if 'reciprocal lattice vectors' in line:
            ibasis = i + 1

        if 'E-fermi' in line:
            efermi = float(line.split()[2])
            LineEfermi = i + 1
            

        if 'NELECT' in line:
            nelect = float(line.split()[2])
        
    # Check how many lines contain information about the bands
    # For ispin = 2, there are two extra lines "spin component..."
    N = (nband + 2) * nkpts * ispin + (ispin - 1) * 2
    bands = []
    occupation = []
    
    kpointcounter = -1
    
    for line in outcar[LineEfermi:LineEfermi + N]:
        # Skip non-relevant lines
        if 'spin component' in line or 'band No.' in line:
            continue
        if 'k-point' in line:
            kpointcounter += 1
            continue
        if weightfilter is not None:
            if readornot[kpointcounter]:
                bands.append(float(line.split()[1]))
                occupation.append(float(line.split()[2]))
        else: 
            bands.append(float(line.split()[1]))
            occupation.append(float(line.split()[2]))
    
    if weightfilter is not None:
        nkpts = np.sum(readornot)
    
    bands = np.array(bands, dtype=float).reshape((ispin, nkpts, nband)) 
    occupation = np.array(occupation, dtype=float).reshape((ispin, nkpts, nband)) 
       
    return bands, efermi, qxvals, qyvals, occupation


def Kpath(path,n):
    kpath = []
    for i in range(len(path)-1):
        kp1 = np.array(path[i])
        kp2 = np.array(path[i+1])
        kp  = np.array(kp2-kp1)
        for j in range(n+1):
            kpath.append(kp1+kp*j/n)
    return kpath



if __name__ == '__main__':
    qz = np.array([-13.01,0,7.48])
    qz = qz/np.linalg.norm(qz)
    qy = np.array([0,0.5,0])
    qx = np.array([0.5,0,0])#np.cross(qy,qz/np.linalg.norm(qz))

    npoints = 10

    kpath = fromgridtoKPOINTS(qz,qy, [-0.40,0.40], [0,0], 1000,1)#np.array(Kpath([-qz, np.array([0.,0.,0.]), qz],npoints))
    
    with open('Ag2Te\MBJ_thirdtry\IBZKPT') as f:
        ibzkpt = f.readlines()
    
    nkibz = int(ibzkpt[1])
    
    with open('Ag2Te\perpLineHighRez\KPOINTS','w') as f:
        f.write('Kpoints generated automatically by F.B.\n')
        f.write('{}\n'.format(len(kpath)+nkibz))
        f.write('Reciprocal lattice\n')
        for i in range(3,len(ibzkpt)):
            f.write(ibzkpt[i])
        for i in range(len(kpath)):
                #f.write('{:.5f} {:.5f} {:.5f}\n'.format(kpath[i,0],kpath[i,1],kpath[i,2]))
                f.write('{:.5f} {:.5f} {:.5f} {}\n'.format(kpath[i,0],kpath[i,1],kpath[i,2],'0.00'))
    
    '''grid = fromgridtoKPOINTS(qx,qy,[-1,1],[-1,1],10,30)

    fig = plt.figure()
    fig.set_size_inches(15, 8)
    ax = fig.add_subplot(projection='3d')
    ax.scatter(grid[:,0], grid[:,1], grid[:,2])
    plt.show()'''

    '''bands, efermi, qxvals, qyvals, occupation = fromOUTCARtoplot(outcarfile='Ag2Te\PerpGrid\OUTCAR',kpointsfile="Ag2Te\PerpGrid\KPOINTS", qx=qz, qy=qy, weightfilter = 0)

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

