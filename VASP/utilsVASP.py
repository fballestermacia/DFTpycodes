import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import re
from matplotlib.ticker import AutoMinorLocator

def fromKPOINTStogrid(kpointsfile='KPOINTS', qx = np.array([1,0,0]), qy = np.array([0,1,0]), qz = np.array([0,0,1]), weightfilter = None):
    
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
    qzvals = []

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
        
        if RECLAT:
            coefx = np.dot(np.array([kvalx,kvaly,kvalz]),qx)/np.linalg.norm(qx)**2
            coefy = np.dot(np.array([kvalx,kvaly,kvalz]),qy)/np.linalg.norm(qy)**2
            coefz = np.dot(np.array([kvalx,kvaly,kvalz]),qz)/np.linalg.norm(qz)**2
        else: 
            coefx = kvalx
            coefy = kvaly
            coefz = kvalz
        
        qxvals.append(coefx)
        qyvals.append(coefy)
        qzvals.append(coefz)
    
    return np.array(qxvals), np.array(qyvals), np.array(qzvals), readornot


def fromgridtoKPOINTS(qx, qy, qxlims, qylims, nxpoints, nypoints, origin = np.array([0,0,0])):
    NX = np.linspace(qxlims[0], qxlims[1], nxpoints)
    NY = np.linspace(qylims[0], qylims[1], nypoints)
    
    QX, QY = np.meshgrid(NX,NY)
    
    KX = np.outer(QX.flatten(),qx) + origin
    KY = np.outer(QY.flatten(),qy) + origin

    KGRID = KX + KY
    
    return KGRID


def fromOUTCARtoplot(outcarfile = 'OUTCAR', kpointsfile = 'KPOINTS', qx = np.array([1,0,0]), qy = np.array([0,1,0]), qz = np.array([0,0,1]), weightfilter = None):
    
    qxvals, qyvals, qzvals, readornot = fromKPOINTStogrid(kpointsfile=kpointsfile, qx = qx, qy = qy, qz=qz, weightfilter=weightfilter)
    
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
       
    return bands, efermi, qxvals, qyvals, qzvals, occupation


def Kpath(path,n):
    kpath = []
    for i in range(len(path)-1):
        kp1 = np.array(path[i])
        kp2 = np.array(path[i+1])
        kp  = np.array(kp2-kp1)
        for j in range(n+1):
            kpath.append(kp1+kp*j/n)
    return kpath


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


def readOutcar(outcarfile='OUTCAR', kpointsfile='KPOINTS'):
    """Obtain the band energies from OUTCAR file generated by VASP.

    Args:
        outcarfile (str, optional): Path to OUTCAR file. Defaults to 'OUTCAR'.
        kpointsfile (str, optional): Path to KPOINTS file. Defaults to 'KPOINTS'.
    """
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
      
    # Get basis.       
    # When the supercell is too large, spaces are missing between real space
    # lattice constants. A bug found out by Wei Xie (weixie4@gmail.com) and
    # code provided by Irián ¿Sánchez?
    B = np.array([line.split()[-3:] for line in outcar[ibasis:ibasis+3]],
                 dtype=float)
    
    # k-points vectors and weights
    dummy = np.array([line.split() for line in outcar[Lvkpts:Lvkpts+nkpts]],
                   dtype=float)
    vkpts = dummy[:, :3]
    wkpts = dummy[:, -1]
    
    # Check how many lines contain information about the bands
    # For ispin = 2, there are two extra lines "spin component..."
    N = (nband + 2) * nkpts * ispin + (ispin - 1) * 2
    bands = []
    filledfactor = []
    
    for line in outcar[LineEfermi:LineEfermi + N]:
        # Skip non-relevant lines
        if 'spin component' in line or 'band No.' in line:
            continue
        if 'k-point' in line:
            continue
        bands.append(float(line.split()[1]))
        filledfactor.append(float(line.split()[2]))
        
    bands = np.array(bands, dtype=float).reshape((ispin, nkpts, nband))
    filledfactor = np.array(filledfactor, dtype=float).reshape((ispin, nkpts, nband))
    
    kp = open(kpointsfile).readlines()
    
    if kp[2][0].upper() == 'L':
        Nk_in_seg = int(kp[1].split()[0]) # Number of k-points per line
        Nseg = nkpts // Nk_in_seg # Number of lines
        vkpt_diff = np.zeros_like(vkpts, dtype=float)

        for ii in range(Nseg):
            start = ii * Nk_in_seg
            end = (ii + 1) * Nk_in_seg
            vkpt_diff[start:end, :] = vkpts[start:end, :] - vkpts[start, :]

        kpt_path = np.linalg.norm(np.dot(vkpt_diff, B), axis=1)
        
        for ii in range(1, Nseg):
            start = ii * Nk_in_seg
            end = (ii + 1) * Nk_in_seg
            kpt_path[start:end] += kpt_path[start-1]

        
        kpt_bounds = np.concatenate((kpt_path[0::Nk_in_seg], [kpt_path[-1], ]))

    return kpt_path, bands, efermi, kpt_bounds, wkpts, nelect, filledfactor

def bandplot(kpath, bands, efermi, kpt_bounds, nelect, kpointsfile =  "KPOINTS",fig_size = (15,8),filledvsempty=False,filledfactor = None,e_window = None,line_w = 0.5,fig_name = "bands.pdf",fig_title=None):
    """Plot Band Structure using matplotlib

    Args:
        kpath: Array of kpoints.
        bands: Array of energy values for each kpoint.
        efermi: Fermi energy.
        kpt_bounds: Values of high-symmetry points.
        nelect: Number ofelectrons
        kpointsfile: Path to KPOINTS file. Defaults to 'KPOINTS'. Defaults to "KPOINTS".
        fig_size: Size of matplotlib figure. Defaults to (15,8).
        filledvsempty: Wether to plot filled bands blue and empty bands red. If True, a filledfactor variable must be provided. Defaults to False.
        filledfactor: Occupation coeficient for each band and each kpoint. Defaults to None.
        e_window: Energy window for the plot (yaxis limits). Defaults to None.
        line_w: Line width for the plot. Defaults to 0.5.
        fig_name: Name for the figure file. Defaults to "bands.pdf".
        fig_title: Title of the matplotlib figure. Defaults to None.
    """
    
    width, height = fig_size
    
    if e_window != None:
        ymin, ymax = e_window[0], e_window[1]
    else:
        ymin, ymax = -4.0,4.0
        
    fig = plt.figure()
    fig.set_size_inches(width, height)
    ax = plt.subplot(111)
    
    nspin, nkpts, nbands = bands.shape
    
    clrs = ['r', 'b']
    
    for Ispin in range(nspin):
        for Iband in range(nbands):
            color = 'k'
            if filledvsempty and filledfactor is not None:
                if filledfactor[Ispin, 0,Iband].any() == 1.:
                    color = 'b'
                if filledfactor[Ispin, 0,Iband].any() == 0.:
                    color = 'r'

            if float(Iband) >= nelect:
                line, = ax.plot(kpath, bands[Ispin, :, Iband], lw=line_w, zorder=0,
                                alpha=0.8, color=color)
            else:
                line, = ax.plot(kpath, bands[Ispin, :, Iband], lw=line_w, zorder=0,
                                alpha=0.8, color=color)
    
    # add extra horizontal/vertical lines
    
    for bd in kpt_bounds:
        ax.axvline(x=bd, ls='-', color='k', lw=0.5, alpha=0.5)
        
    ax.set_ylabel('$E - E_f$ [eV]',  # fontsize='small',
                  labelpad=5)
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(kpath.min(), kpath.max())
    
    ax.set_xticks(kpt_bounds)
    
    # Read and use KPOINTS file for tick labels.

    with open(kpointsfile,'r') as KPointsFile:
        TmpFlag = 0;
        TmpLabels = [];
        for TmpLine in KPointsFile:
            TmpLine = TmpLine.strip()
            if TmpFlag == 1:                    
                TmpLine = re.sub(r'^.*\!\s?', '', TmpLine)
                TmpLine.strip()
                if TmpLine != "":
                    if TmpLine == "G" or TmpLine == "Gamma" or TmpLine == "GM":
                        TmpLabels.append(r'$\mathrm{{\mathsf{\Gamma}}}$')
                    else:
                        TmpLabels.append(r'$\mathrm{\mathsf{'+TmpLine+'}}$')
            if (TmpLine == "reciprocal") | (TmpLine == "rec"):
                TmpFlag = 1
        TmpLabels2 = [TmpLabels[0]]
        TmpIndex = 1
        while TmpIndex < (len(TmpLabels) - 1):
            if TmpLabels[TmpIndex + 1] == TmpLabels[TmpIndex]:
                TmpLabels2.append(TmpLabels[TmpIndex])
            else:
                TmpLabels2.append(TmpLabels[TmpIndex]+'|'+TmpLabels[TmpIndex + 1])
            TmpIndex += 2
        
        TmpLabels2.append(TmpLabels[len(TmpLabels) - 1])
        ax.set_xticklabels(TmpLabels2)
        
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.axhline(y=0, xmin=0, xmax=1, linestyle='dotted', color='black', linewidth=0.5)
        ax.set_title(fig_title)
        plt.tight_layout(pad=0.20)
        plt.savefig(fig_name)
        plt.show()
        return kpath, bands
    

def plot_bands_from_VASP(outcarfile = "OUTCAR", kpointsfile = "KPOINTS",fig_size = (10,7.5),ewindow =  None,filledvsempty=False,fig_name = "bands.pdf",fig_title=None,efermi=None):
    """ Plot the band structure from KPOINTS and OUTCAR files.

    Args:
        outcarfile (str, optional): Path to OUTCAR file. Defaults to 'OUTCAR'.
        kpointsfile (str, optional): Path to KPOINTS file. Defaults to 'KPOINTS'.
        fig_size: Size of matplotlib figure. Defaults to (10,7.5).
        e_window: Energy window for the plot (yaxis limits). Defaults to None.
        filledvsempty: Wether to plot filled bands blue and empty bands red. Defaults to False.
        fig_name: Name for the figure file. Defaults to "bands.pdf".
        fig_title: Title of the matplotlib figure. Defaults to None.
        efermi: Fermi energy if previously known.. Defaults to None.
    """
    
    # Use non-interactive backend in case there is no display
    mpl.use('agg')
    mpl.rcParams['axes.unicode_minus'] = False
    
    kpath, bands, efermi_, kpt_bounds, wkpts, nelect, filledfactor = readOutcar(outcarfile,kpointsfile)
    
    if efermi == None:
        efermi = efermi_
    
    
    if filledvsempty:
        bandplot(kpath, bands-efermi, efermi,kpt_bounds, nelect,kpointsfile=kpointsfile, e_window=ewindow, fig_size=fig_size, fig_name=fig_name,fig_title=fig_title, filledvsempty=filledvsempty, filledfactor=filledfactor)
    else:
        bandplot(kpath, bands-efermi, efermi,kpt_bounds, nelect,kpointsfile=kpointsfile, e_window=ewindow, fig_size=fig_size, fig_name=fig_name,fig_title=fig_title)


