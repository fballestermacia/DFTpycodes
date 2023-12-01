import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
from scipy.interpolate import splev, splrep
from numpy.polynomial.polynomial import polyder, polyval, polyfit
from get_effectiveMass1D import fromOUTCARtoplot1D, gausiana
from get_effectiveMass_FiniteDiff import secondDer_FiniteDiff

def get_mzeffs(coefs, band):
    '''pol = polyfit(coefs,band, polyorder)
    
    derz = polyder(pol,2)
    
    derzvals = polyval(coefs,derz)'''

    derzvals = secondDer_FiniteDiff(band, coefs)

    
    factor = 1#1000*(6.582e-16)**2/0.510e6*3e8**2/1e-10**2*4*np.pi**2
    
    
    mz = np.power(derzvals,-1)*factor
    #mz = np.where(mz<0, 0, mz)
    
    enval = efermi#np.min(band)
    sigma = 0.0005
    a0=1/(sigma*np.sqrt(2*np.pi))
    
    envals = np.linspace(np.min(band), np.min(band)+250/1000*0.0367493, 1000)
    
    #bandgaussian = gausiana(band,enval,sigma,a0)
    #mzeff = np.sum(mz*bandgaussian)/simps(bandgaussian)
    
    #print("m_effZ = ", mzeff)
    
    mzeffs = []
    iderzvals = []
    
    for enval2 in envals:
        bandgaussian = gausiana(band,enval2,sigma,a0)
        iderzvals.append(np.sum(derzvals*bandgaussian)/simps(bandgaussian))
        mzeffs.append(np.power(np.sum(derzvals*bandgaussian)/simps(bandgaussian),-1))
    
    return mzeffs, envals, iderzvals, enval, mz

if __name__ == '__main__':
    
    qz = np.array([-13.01,0,7.48])/np.linalg.norm(np.array([-13.01,0,7.48]))
    
    bands, efermi, coefs, occupation = fromOUTCARtoplot1D(outcarfile = 'Ag2Te\\perpLineHighRezShort\\OUTCAR', kpointsfile = 'Ag2Te\\perpLineHighRezShort\\KPOINTS', vec=qz, weightfilter = 0.)
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
    
    scaling = 1/1.88973 #Bohr
    coefs *= scaling
    
    
    '''Enlim = 0.5*0.0367493+zeroenergy*0.0367493
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

   

        
    mzeffs, envals, derzvals, enval, mz = get_mzeffs(coefs, band)
    
    zero = zeroenergy#np.min(band/0.0367493)
    envals = envals/0.0367493-zero
    
    ax1.plot(coefs, band/0.0367493-zero)
    
    ax2.plot(envals,derzvals)
    
    #ax3.plot(coefs,mz, label='C = {}'.format(polyorder))
    ax3.plot(envals,np.power(derzvals,-1))

    ax4.plot(envals, mzeffs)
    #print(mzeffs[0])
    ax4.vlines([efermidoped10,  efermidoped], np.min(mzeffs), np.max(mzeffs), colors=['r','k'])
    
    ax5.plot(envals*1000, mzeffs)
    ax5.vlines([efermidoped10a*1000,  efermidopeda*1000], 0.15,0.35, colors=['r','k'])
    ax5.vlines([efermidoped10b*1000,  efermidopedb*1000], 0.15,0.35, colors=['r','k'], linestyles='--')

    #ax1.set_ylim(0.275,0.365)
    #ax2.set_ylim(-5,7.5)
    #ax3.set_ylim(-5,5)
    #ax4.set_ylim(0,0.35)
    #ax5.set_ylim(0,0.35)
    
    
    plt.legend()
    plt.show()