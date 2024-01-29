import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
import utilsVASP
 
 
if __name__ == '__main__':       
    Energy, DOS, iDOS, Emax, Emin, Erange, NEDOS, Efermi = utilsVASP.fromDOSCARtoarray( doscarfile='Ag2Te\MBJ_thirdtry\DOSCAR')


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
