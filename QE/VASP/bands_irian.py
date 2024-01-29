from pymatgen.io.vasp.outputs import Vasprun
import pymatgen.core.structure as pst
import pymatgen.symmetry.analyzer as psa
from pymatgen.electronic_structure.plotter import BSPlotter
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib as mpl
from matplotlib.ticker import AutoMinorLocator
from scipy.interpolate import make_interp_spline
import itertools as it
import re
import bz2

rc('text', usetex=True)
rc('font', size=14)
rc('legend', fontsize=13)
rc('text.latex')#, preamble=r'\usepackage{cmbright}')

def get_bandInfo1(inFile='OUTCAR',kpointfile = "KPOINTS"):

    """
    extract band energies from OUTCAR
    """

    outcar = [line for line in open(inFile) if line.strip()]

    for ii, line in enumerate(outcar):
        if 'NKPTS =' in line:
            nkpts = int(line.split()[3])
            nband = int(line.split()[-1])

        if 'ISPIN  =' in line:
            ispin = int(line.split()[2])

        if "k-points in reciprocal lattice and weights" in line:
            Lvkpts = ii + 1

        if 'reciprocal lattice vectors' in line:
            ibasis = ii + 1

        if 'E-fermi' in line:
            efermi = float(line.split()[2])
            LineEfermi = ii + 1
            # break

        if 'NELECT' in line:
            nelect = float(line.split()[2])
            # break

    # basis vector of reciprocal lattice
    # B = np.array([line.split()[3:] for line in outcar[ibasis:ibasis+3]],

    # When the supercell is too large, spaces are missing between real space
    # lattice constants. A bug found out by Wei Xie (weixie4@gmail.com).
    B = np.array([line.split()[-3:] for line in outcar[ibasis:ibasis+3]],
                 dtype=float)
    # k-points vectors and weights
    tmp = np.array([line.split() for line in outcar[Lvkpts:Lvkpts+nkpts]],
                   dtype=float)
    vkpts = tmp[:, :3]
    wkpts = tmp[:, -1]

    # for ispin = 2, there are two extra lines "spin component..."
    N = (nband + 2) * nkpts * ispin + (ispin - 1) * 2
    bands = []
    # vkpts = []
    for line in outcar[LineEfermi+1:LineEfermi + N+1]:
        if 'spin component' in line or 'band No.' in line:
            continue
        if 'k-point' in line:
            # vkpts += [line.split()[3:]]
            continue
        bands.append(float(line.split()[1]))

    bands = np.array(bands, dtype=float).reshape((ispin, nkpts, nband))


    kp = open(kpointfile).readlines()

    if kp[2][0].upper() == 'L':
        Nk_in_seg = int(kp[1].split()[0])
        Nseg = nkpts // Nk_in_seg
        vkpt_diff = np.zeros_like(vkpts, dtype=float)

        for ii in range(Nseg):
            start = ii * Nk_in_seg
            end = (ii + 1) * Nk_in_seg
            vkpt_diff[start:end, :] = vkpts[start:end, :] - vkpts[start, :]

        kpt_path = np.linalg.norm(np.dot(vkpt_diff, B), axis=1)
        # kpt_path = np.sqrt(np.sum(np.dot(vkpt_diff, B)**2, axis=1))
        for ii in range(1, Nseg):
            start = ii * Nk_in_seg
            end = (ii + 1) * Nk_in_seg
            kpt_path[start:end] += kpt_path[start-1]

        # kpt_path /= kpt_path[-1]
        kpt_bounds = np.concatenate((kpt_path[0::Nk_in_seg], [kpt_path[-1], ]))

    return kpt_path, bands, efermi, kpt_bounds, wkpts, nelect

def get_bandInfo2(inFile='OUTCAR',kpointfile = "KPOINTS"):

    """
    extract band energies from OUTCAR
    """

    outcar = [line for line in open(inFile) if line.strip()]

    for ii, line in enumerate(outcar):
        if 'NKPTS =' in line:
            nkpts = int(line.split()[3])
            nband = int(line.split()[-1])

        if 'ISPIN  =' in line:
            ispin = int(line.split()[2])

        if "k-points in reciprocal lattice and weights" in line:
            Lvkpts = ii + 1

        if 'reciprocal lattice vectors' in line:
            ibasis = ii + 1

        if 'E-fermi' in line:
            efermi = float(line.split()[2])
            LineEfermi = ii + 1
            # break

        if 'NELECT' in line:
            nelect = float(line.split()[2])
            # break

    # basis vector of reciprocal lattice
    # B = np.array([line.split()[3:] for line in outcar[ibasis:ibasis+3]],

    # When the supercell is too large, spaces are missing between real space
    # lattice constants. A bug found out by Wei Xie (weixie4@gmail.com).
    B = np.array([line.split()[-3:] for line in outcar[ibasis:ibasis+3]],
                 dtype=float)
    # k-points vectors and weights
    tmp = np.array([line.split() for line in outcar[Lvkpts:Lvkpts+nkpts]],
                   dtype=float)
    vkpts = tmp[:, :3]
    wkpts = tmp[:, -1]

    # for ispin = 2, there are two extra lines "spin component..."
    N = (nband + 2) * nkpts * ispin + (ispin - 1) * 2
    bands = []
    # vkpts = []
    for line in outcar[LineEfermi:LineEfermi + N]:
        if 'spin component' in line or 'band No.' in line:
            continue
        if 'k-point' in line:
            # vkpts += [line.split()[3:]]
            continue
        bands.append(float(line.split()[1]))

    bands = np.array(bands, dtype=float).reshape((ispin, nkpts, nband))


    kp = open(kpointfile).readlines()

    if kp[2][0].upper() == 'L':
        Nk_in_seg = int(kp[1].split()[0])
        Nseg = nkpts // Nk_in_seg
        vkpt_diff = np.zeros_like(vkpts, dtype=float)

        for ii in range(Nseg):
            start = ii * Nk_in_seg
            end = (ii + 1) * Nk_in_seg
            vkpt_diff[start:end, :] = vkpts[start:end, :] - vkpts[start, :]

        kpt_path = np.linalg.norm(np.dot(vkpt_diff, B), axis=1)
        # kpt_path = np.sqrt(np.sum(np.dot(vkpt_diff, B)**2, axis=1))
        for ii in range(1, Nseg):
            start = ii * Nk_in_seg
            end = (ii + 1) * Nk_in_seg
            kpt_path[start:end] += kpt_path[start-1]

        # kpt_path /= kpt_path[-1]
        kpt_bounds = np.concatenate((kpt_path[0::Nk_in_seg], [kpt_path[-1], ]))
        
    return kpt_path, bands, efermi, kpt_bounds, wkpts, nelect


def bandplot(kpath, bands, efermi, kpt_bounds, nelect, kpointfile =  "KPOINTS",fig_size = (15,8),e_window = None,line_w = 0.5,fig_name = "bands.pdf",fig_title=None):
    '''
    Use matplotlib to plot band structure
    '''

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
            
            # if Iband == 0 else line.get_color()
            lc = None if Iband == 0 else line.get_color()
            
            #if nspin == 1:
            #    new_nelect = nelect/2
            if float(Iband) >= nelect:
                line, = ax.plot(kpath, bands[Ispin, :, Iband], color='k',lw=line_w, zorder=0,
                                alpha=0.8)
            else:
                line, = ax.plot(kpath, bands[Ispin, :, Iband], lw=line_w, zorder=0,
                                alpha=0.8,
                                color='k')

    for bd in kpt_bounds:
        ax.axvline(x=bd, ls='-', color='k', lw=0.5, alpha=0.5)

    # add extra horizontal/vertical lines

    ax.set_ylabel('$E - E_f$ [eV]',  # fontsize='small',
                  labelpad=5)
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(kpath.min(), kpath.max())

    ax.set_xticks(kpt_bounds)
    
    # Read and use kpoint file

    with open(kpointfile,'r') as KPointsFile:
        TmpFlag = 0;
        TmpLabels = [];
        for TmpLine in KPointsFile:
            TmpLine = TmpLine.strip()
            if TmpFlag == 1:                    
                TmpLine = re.sub(r'^.*\!\s?', '', TmpLine)
                TmpLine.strip()
                if TmpLine != "":
                    if TmpLine == "G":
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
    plt.show()
    plt.savefig(fig_name)
    return kpath, bands

def plot_vaspbands(outcar = "OUTCAR", kpoints = "KPOINTS",fig_size = (10,7.5),ewindow =  None,fig_name = "bands.pdf",fig_title=None,efermi=None):

    """ Plot vasp bands using the aforedefined functions """

    # Use non-interactive backend in case there is no display
    mpl.use('agg')
    mpl.rcParams['axes.unicode_minus'] = False
    try:    
        kpath, bands, efermi_, kpt_bounds, wkpts, nelect = get_bandInfo1(outcar,kpoints)
    except IndexError:
        kpath, bands, efermi_, kpt_bounds, wkpts, nelect = get_bandInfo2(outcar,kpoints)
    if efermi == None:
        efermi = efermi_
    bandplot(kpath,bands-efermi,efermi,kpt_bounds,nelect,kpointfile=kpoints,e_window = ewindow,fig_size=fig_size,fig_name=fig_name,fig_title=fig_title)
    
plot_vaspbands(outcar="Ag2Te\MBJ\OUTCAR", kpoints="Ag2Te\MBJ\KPOINTS",fig_name = "Ag2Te\\bands_MBJ.pdf", ewindow=[-1,1])
