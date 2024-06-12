import matplotlib.pyplot as plt
import numpy as np
import re
import glob
import platform
import utilsQE as utqe
import constantsQE as ctqe

def readfolder(folder, extension=None, shiftEnergy=True):
    """
    Folder must contain subfolders labeled by the kgrid number of points,
    each subfolder must contain subfolders labeled by the cutoff energy.

    If extension is None, the program detects the extension used for the outputs
    but it can be specified to maybe filter only some files and not others.

    If shiftEnergy is True, the code shifts the origin of energies to the minimum obtained.

    Returns a numpy array where each element is an array of data structured as:
    [GRID, CUTOFF, TOTAL ENERGY PER ATOM (in meV), FERMI ENERGY (in eV), TOTAL RUNTIME (in hours), TOTAL RAM (in GB), TOTAL FORCES PER ATOM (in meV)]
    """

    #####################################
    #CHECK OS FOR DIRECTORY MANIPULATION#
    #####################################

    if platform.system() == 'Windows':
        dirp = '\\'
    else: dirp ='/'

    kgrids = [i.split(dirp)[-1] for i in glob.glob(folder+dirp+'*')]
    kgrids = utqe.orderedKgrid(kgrids)

    data = []

    ######################################
    #NOTE: this asumes all subdirs have the same files!!
    ######################################
    subdirs = glob.glob(folder+dirp+kgrids[0]+dirp+'*')
    files = glob.glob(subdirs[0]+dirp+'*')
    
    if extension is None:
        for file in files:
            if re.search('o', file.split('.')[-1]):
                extension = file.split('.')[-1]
                break
    
    for i in range(len(kgrids)):
        currentgrid = kgrids[i]
        files = glob.glob(folder+dirp+currentgrid+dirp+'*'+dirp+'*'+extension)
        natoms = utqe.readNAtoms(files[0])
        for j, file in enumerate(files):
            cutoff = files[j].split(dirp)[-2]
            forces = utqe.readForces(file)
            energy = utqe.readTotalEnergy(file)
            efermi = utqe.readEfermi(file)
            time = utqe.readTime(file)
            ram = utqe.readRAM(file)
            data.append([currentgrid, float(cutoff), float(energy)*ctqe.Ry2meV/natoms,float(efermi),float(time),float(ram),float(forces)*ctqe.Ry2meV/natoms])
    
    data = np.array(data)

    if shiftEnergy:
        min = np.min([float(x) for x in data[:,2]])
        data[:,2] = [float(x) for x in data[:,2]] - min
    return np.array(data)


def sortforCutoff(arrays):
    """
    Returns data in ascending order of the cutoff energy.
    """

    for i,arr in enumerate(arrays):
        arrays[i] = arr[arr[:,0].argsort()]
    
    return arrays


def dataforeachKgrid(data):
    """
    Returns the data split in two arrays, one for the kgrids and one with the rest of the data ordered in the same way as the kgrids.
    """

    kgrids, indexes = np.unique(data[:,0], return_index=True)
    kgrids = kgrids[indexes.argsort()]
    
    data4kgrid = []
    for kgrid in kgrids:
        dummy = []
        for i in range(len(data[:,0])):
            if data[i,0] == kgrid:
                dummy.append([float(x) for x in data[i,1:]])
        data4kgrid.append(dummy)
    return np.array(kgrids, dtype=str), sortforCutoff(np.array(data4kgrid))






if __name__ == '__main__':
#[GRID, CUTOFF, TOTAL ENERGY PER ATOM (in meV), FERMI ENERGY (in eV), TOTAL RUNTIME (in hours), TOTAL RAM (in Gb), TOTAL FORCES PER ATOM (in meV)]
    data = readfolder(r'data/TaCo2Te2/RoomTemp/Convergence')
    kgrids, data4kgrid = dataforeachKgrid(data)
    

    fig1 = plt.figure(1)
    fig2 = plt.figure(2)
    fig3 = plt.figure(3)
    fig4 = plt.figure(4)
    fig5 = plt.figure(5)

    plt.figure(1)
    ax1 = plt.axes()
    
    plt.figure(2)
    ax2 = plt.axes()

    plt.figure(3)
    ax3 = plt.axes()

    plt.figure(4)
    ax4 = plt.axes()

    plt.figure(5)
    ax5 = plt.axes()



    for i, kg in enumerate(kgrids):
        ax1.plot(data4kgrid[i,:,0],data4kgrid[i,:,1],'.-', label=kg)
        ax2.plot(data4kgrid[i,:,0],data4kgrid[i,:,-1],'.-', label=kg)
        ax3.plot(data4kgrid[i,:,0],data4kgrid[i,:,2],'.-', label=kg)
        ax4.plot(data4kgrid[i,:,0],data4kgrid[i,:,3],'.-', label=kg)
        ax5.plot(data4kgrid[i,:,0],data4kgrid[i,:,4],'.-', label=kg)

    ax1.legend()
    ax2.legend()
    ax3.legend()
    ax4.legend()
    ax5.legend()

    ax1.set_xlabel("Cutoff Energy (Ry)")
    ax2.set_xlabel("Cutoff Energy (Ry)")
    ax3.set_xlabel("Cutoff Energy (Ry)")
    ax4.set_xlabel("Cutoff Energy (Ry)")
    ax5.set_xlabel("Cutoff Energy (Ry)")

    ax1.set_ylabel("Total Energy per Atom (meV)")
    ax2.set_ylabel("Total Force per Atom (meV)")
    ax3.set_ylabel("Fermi Energy (eV)")
    ax4.set_ylabel("Runtime (hours)")
    ax5.set_ylabel("Ram usage (Gb)")


    plt.show()
