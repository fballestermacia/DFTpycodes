import numpy as np
import matplotlib.pyplot as plt
from readDOSCAR import fromDOSCARtoarray
import scipy.constants as cte
from scipy.integrate import simps, quad
from scipy.interpolate import CubicSpline, PPoly

def g(E,m):
    """
    E = energy in meV
    m = mass in kg
    Gives the DOS per cm**3 and per meV"""
    
    jul2eV = 6.241509e18
    
    E=E/(jul2eV*1000) #change to jul
    DOS=(m/(np.pi**2 * cte.hbar**3))*np.sqrt(2*m*E)
    DOS=DOS/(100**3) #change to cm3
    DOS=DOS/(jul2eV*1000) #change to "per meV"
    return DOS 

def integrated_dos(E,m):
    res, err = quad(g,0,E,args=(m,))
    return res


def fermitoDoping(efermi, meff=1):
    
    nel = integrated_dos(efermi,cte.m_e*meff)
    
    return nel

def Dopingtofermi(nel, meff=1, iguess = None, Elims = None):    
    maxiter = 1e4
    precission = 1e-4
    
    iter = 0
    

    if iguess is not None:
        Efpost = iguess
    else: Efpost = 1
    
    if Elims is not None:
        liminf, limsup = Elims
    else: liminf, limsup = 0, 1e3
    
    Efant = liminf
    
    nup=0
    while iter < maxiter and np.abs(Efant-Efpost) > precission:
        iter += 1
        nup = fermitoDoping(Efpost, meff=meff)
        if nup < nel:
            liminf = Efpost
            dummy = Efpost
            Efpost = (Efpost+limsup)/2
            Efant = dummy
        elif nup > nel:
            
            limsup = Efpost
            dummy = Efpost
            Efpost = (Efpost+liminf)/2
            Efant = dummy
        elif nup == nel:
            Efant = Efpost 
    
    return Efpost

doscar='Ag2Te\MBJ_thirdtry\DOSCAR'

Energy, DOS, iDOS, Emax, Emin, Erange, NEDOS, Efermi = fromDOSCARtoarray( doscarfile=doscar)

#print(Efermi)


dfile = open(doscar).readlines()

vol = float(dfile[1].strip().split()[0])
print(vol)  
initTemp = float(dfile[2].strip())

infoline = [float(x) for x in dfile[5].strip().split()]

#Erange is not always present, if not, it will be equal to NEDOS.
Emax, Emin, Erange, NEDOS, Efermidoscar = infoline[0], infoline[1], infoline[2], infoline[-3], infoline[-2] 

'''  8.0899999999999999    0.0000000000000000    0.0000000000000000
     0.0000000000000000    4.4800000000000004    0.0000000000000000
    -4.9231649190557025    0.0000000000000000    7.4862572210537408'''
    
a1 = np.array([8.0899999999999999,0,0])
a2 = np.array([0,4.4800000000000004,0])
a3 = np.array([-4.9231649190557025, 0.0000000000000000, 7.4862572210537408])

vol = np.dot(a1,np.cross(a2,a3))

n = (1.65e18)*10 #/cm^3 #/(1e8)**3

alat = 13.5188898664
#vol *= 1/alat 

print('vol=',vol)



doping_vol=  n*vol/((1e8)**3)

print('doping*vol = ', doping_vol)

Es = np.linspace(0, 450, 5000)

#print(integrated_dos(5.1, cte.m_e))


meff = 0.12#0.1

#print('112/vol = ', 112/vol)
print('efermi = ',Dopingtofermi(n,meff=meff, iguess=100))
print('dop_fermi = ',fermitoDoping(Dopingtofermi(n,meff=meff, iguess=100),meff = meff))
#print('sumDOS/vol = ', np.sum((DOS[Energy <= Efermi+0.1])[-1])/vol)

ig = []
for E in Es:
    ig.append(integrated_dos(E,cte.m_e*meff))

plt.figure()
plt.plot(Es,ig)
plt.show()