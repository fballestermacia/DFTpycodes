import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
from scipy.interpolate import splev, splrep
from numpy.polynomial.polynomial import polyder, polyval, polyfit

def fromKPOINTStoline(kpointsfile='KPOINTS', vec=np.array([0.5,0.,0.]), weightfilter = None):
    
    kp = open(kpointsfile).readlines()
    
    RECLAT = False
    if kp[2][0].upper() == 'L':
        startline = 4
    else: startline = 3
    '''elif kp[2][0].upper() == 'R':
        startline = 3
        RECLAT = True'''
    
    
    coefs = []
    
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
        
        coefx = np.dot(np.array([kvalx,kvaly,kvalz]),vec)/np.linalg.norm(vec)**2
        
        coefs.append(coefx)
    
    return np.array(coefs), readornot


def fromOUTCARtoplot1D(outcarfile = 'OUTCAR', kpointsfile = 'KPOINTS', vec=np.array([0.5,0.,0.]), weightfilter = None):
    
    coefs, readornot = fromKPOINTStoline(kpointsfile=kpointsfile, vec=vec, weightfilter=weightfilter)
    
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
       
    return bands, efermi, coefs, occupation


def gausiana(x,center,sigma,a0):
    return a0*np.exp(-(x-center)**2/sigma**2/2)


def get_mzeffs(coefs, band, polyorder):
    pol = polyfit(coefs,band, polyorder)
    
    derz = polyder(pol,2)
    
    derzvals = polyval(coefs,derz)

    
    factor = 1#1000*(6.582e-16)**2/0.510e6*3e8**2/1e-10**2*4*np.pi**2
    
    
    mz = np.power(derzvals,-1)*factor
    #mz = np.where(mz<0, 0, mz)
    
    enval = efermi#np.min(band)
    sigma = 0.0005
    a0=1/(sigma*np.sqrt(2*np.pi))
    
    envals = np.linspace(np.min(band), np.min(band)+250/1000*0.0367493, 10000)
    
    #bandgaussian = gausiana(band,enval,sigma,a0)
    #mzeff = np.sum(mz*bandgaussian)/simps(bandgaussian)
    
    #print("m_effZ = ", mzeff)
    
    mzeffs = []
    iderzvals = []
    
    for enval2 in envals:
        bandgaussian = gausiana(band,enval2,sigma,a0)
        iderzvals.append(np.sum(derzvals*bandgaussian)/simps(bandgaussian))
        mzeffs.append(np.power(np.sum(derzvals*bandgaussian)/simps(bandgaussian),-1))
    
    
    return mzeffs, envals, iderzvals, pol, enval, mz

if __name__ == '__main__':
    
    qz = np.array([-13.01,0,7.48])/np.linalg.norm(np.array([-13.01,0,7.48]))
    
    bands, efermi, coefs, occupation = fromOUTCARtoplot1D(outcarfile = 'Ag2Te\\perpLineHighRezLong\\OUTCAR', kpointsfile = 'Ag2Te\\perpLineHighRezLong\\KPOINTS', vec=qz, weightfilter = 0.)
    #print(efermi)
    #coefs = np.abs(coefs)
    #print(occupation)
    dummy = np.transpose(bands[0])
    dummyoc = np.transpose(occupation[0])

    index = np.sum(dummyoc, axis=0)
    
    
    zeroenergy = np.min(((dummy)[int(index[0])]))
    
    
    efermidoped10 = 124.49207603931427/1000#236.33401095867157/1000
    efermidoped = 26.821041107177734/1000#50.916579365730286/1000
    
    
    efermidoped10a = 124.49207603931427/1000#236.33401095867157/1000
    efermidopeda = 26.821041107177734/1000#50.916579365730286/1000
    
    efermidoped10b = 236.33401095867157/1000
    efermidopedb = 50.916579365730286/1000
    
    band = ((dummy)[int(index[0])])*0.0367493 #hartree
    
    polyorders = [20]#(np.arange(3,10))*2
    
    scaling = 1/1.88973 #Bohr
    coefs *= scaling
    
    
    '''Enlim = 0.4
    trues = np.where(band<Enlim)
    coefs = coefs[trues[0]]
    band = band[trues[0]]'''
    
    
    '''pol = polyfit(coefs,band, polyorder)
    
    derz = polyder(pol,2)
    
    derzvals = polyval(coefs,derz)
    
    factor = 1#1000*(6.582e-16)**2/0.510e6*3e8**2/1e-10**2*4*np.pi**2
    
    
    mz = np.power(derzvals,-1)*factor
    
    
    enval = efermi*0.0367493#np.min(band)
    sigma = 0.001
    a0=1
    
    envals = np.linspace(np.min(band), np.min(band)+0.01, 100)
    
    bandgaussian = gausiana(band,enval,sigma,a0)
    mzeff = np.sum(mz*bandgaussian)/simps(bandgaussian)
    
    print("m_effZ = ", mzeff)
    
    mzeffs = []
    
    for enval2 in envals:
        bandgaussian = gausiana(band,enval2,sigma,a0)
        mzeffs.append(np.sum(mz*bandgaussian)/simps(bandgaussian))'''
    
    fig1 = plt.figure()
    fig1.set_size_inches(10, 6)
    ax1 = fig1.add_subplot(221)
    ax2 = fig1.add_subplot(222)
    ax3 = fig1.add_subplot(223)
    ax4 = fig1.add_subplot(224)
    
    fig2 = plt.figure()
    fig2.set_size_inches(10, 6)
    
    ax5 = fig2.add_subplot(111)
       
    for polyorder in polyorders:
        mzeffs, envals, derzvals, pol, enval, mz = get_mzeffs(coefs, band, polyorder)
        
        zero = zeroenergy#np.min(band/0.0367493)
        envals = envals/0.0367493-zero
        
        ax1.plot(coefs, band)
        ax1.plot(coefs, polyval(coefs,pol), label='C = {}'.format(polyorder))
        
        
        ax2.plot(envals,derzvals, label='C = {}'.format(polyorder))
        
        #ax3.plot(coefs,mz, label='C = {}'.format(polyorder))
        ax3.plot(envals,np.power(derzvals,-1), label='C = {}'.format(polyorder))

        ax4.plot(envals, mzeffs, label='C = {}'.format(polyorder))
        #print(mzeffs[0])
        ax4.vlines([efermidoped10,  efermidoped], np.min(mzeffs), np.max(mzeffs), colors=['r','k'])
        
        ax5.plot(envals*1000, mzeffs, label='C = {}'.format(polyorder))
        ax5.vlines([efermidoped10a*1000,  efermidopeda*1000], 0.15,0.35, colors=['r','k'])
        ax5.vlines([efermidoped10b*1000,  efermidopedb*1000], 0.15,0.35, colors=['r','k'], linestyles='--')
    
    #ax1.set_ylim(0.275,0.365)
    #ax2.set_ylim(-5,7.5)
    #ax3.set_ylim(-5,5)
    #ax4.set_ylim(0,0.35)
    #ax5.set_ylim(0,0.35)
    
    
    plt.legend()
    plt.show()
    
    
    
    
    