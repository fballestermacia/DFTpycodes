import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

def fromDOSCARtoarray(doscarfile = 'DOSCAR'):
    dfile = open(doscarfile).readlines()
    
    initTemp = float(dfile[2].strip())
    
    infoline = [float(x) for x in dfile[5].strip().split()]
    
    #NEDOS is not always present, if not, it will be equal to Efermi.
    Emax, Emin, Erange, NEDOS, Efermi = infoline[0], infoline[1], infoline[2], infoline[3], infoline[-2] 
    
    Energy = []
    DOS = []
    iDOS = []
    
    
    for line in dfile[6:]:
        dummy = [float(x) for x in line.strip().split()]
        dummyEnergy, dummyDOS, dummyiDOS = dummy
        Energy.append(dummyEnergy)
        DOS.append(dummyDOS)
        iDOS.append(dummyiDOS)
    
    return np.array(Energy), np.array(DOS), np.array(iDOS), Emax, Emin, Erange, NEDOS, Efermi
 
 
if __name__ == '__main__':       
    Energy, DOS, iDOS, Emax, Emin, Erange, NEDOS, Efermi = fromDOSCARtoarray( doscarfile='Ag2Te\MBJ_thirdtry\DOSCAR')


    #carrier concentration n = 1.65·10^18 cm-3

    #Efermi *= 4.64158883361  #factor to the fermi energy for a 10 x carrier density
    
    #dummydos = CubicSpline(Energy, DOS)    

    fig = plt.figure()
    fig.set_size_inches(15, 8)
    ax = plt.subplot(111)

    ens = np.linspace(np.min(Energy), np.max(Energy), 1000)

    ax.plot(DOS, Energy-Efermi,'--r', label='DOS')
    ax.plot(iDOS, Energy-Efermi, 'k', label='iDOS')
    #ax.plot(dummydos.__call__(ens), ens-Efermi)
    ax.axhline(Efermi-Efermi, np.min(DOS), np.max(DOS), c='b', linestyle='--')
    ax.legend()
    plt.show() 
