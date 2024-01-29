import matplotlib.pyplot as plt
import numpy as np

def read_filband(file='filband'):
    filelines = [line for line in open(file) if line.strip()]
    
    nbnd = int(filelines[0].split()[2].split(',')[0])
    nks = int(filelines[0].split()[4])
    
    kpointline = True
    for i, line in enumerate(filelines[1:]):
        if kpointline:
            try:
                kpoints=np.vstack([kpoints,np.array(line.split(), dtype='float')])
            except NameError:
                kpoints=np.array(np.array(line.split(), dtype='float'))
            
            bandval = []
            kpointline=False
        else:
            for dummy in line.split():
                bandval.append(float(dummy))
            if len(bandval) == nbnd:
                try:
                    bands=np.vstack([bands,np.array(bandval)])
                except NameError:
                    bands=np.array(bandval)
                kpointline=True
    
    return kpoints, bands

def read_bandsin(bandsinfile='bands.in'):
    filelines = [line for line in open(bandsinfile) if line.strip()]
    
    inKpoints=False
    inCellPar = False
    highsympoints = []
    pointsinline=[]
    
    for i, line in enumerate(filelines):
        if not inCellPar:
            if line.split()[0].upper() == 'CELL_PARAMETERS':
                inCellPar = True
                counter = 0
        else:
            try:
                B=np.vstack([B,np.array(line.split()[-3:],dtype=float)])
            except NameError:
                B=np.array(np.array(line.split()[-3:],dtype=float))
            counter += 1
            if counter == 3: 
                B = np.array(B)
                inCellPar=False
            
            
        
        if not inKpoints:
            if line.split()[0].upper() == 'K_POINTS':
                inKpoints = True
        else:
            if len(line.split()) == 1:
                nlines = int(line)
            else:
                try:
                    kpoints=np.vstack([kpoints,np.array(line.split()[0:3], dtype='float')])
                except NameError:
                    kpoints=np.array(np.array(line.split()[0:3], dtype='float'))
                pointsinline.append(int(line.split()[3]))
                hsp = line.split('!')[1].split('\n')[0]
                if hsp.upper() == 'GAMMA': hsp = '\Gamma'
                highsympoints.append('${}$'.format(hsp))
                
                if len(highsympoints) == nlines: inKpoints=False
    
    return kpoints, nlines, highsympoints, B, pointsinline
                    



    
if __name__ == '__main__':
    
    fig_size = (15,8)
    line_w = 1
    ymin = 4
    ymax = 18
    
    kpoints, bands = read_filband('QE_projects\\AgP2\\bands_norelaxed\\filbandAgP2')
    kpoints2, nlines, hsp, B, pointsinline = read_bandsin('QE_projects\\AgP2\\bands_norelaxed\\AgP2_bands.in')
    nlines = nlines-1
    width, height = fig_size
    
    ticks = [0]
    
    for i in range(nlines):
        ticks.append(int(np.sum(pointsinline[:i+1])))
    
    Nk_in_seg = pointsinline # Number of k-points per line
    Nseg = nlines # Number of lines
    vkpt_diff = np.zeros_like(kpoints, dtype=float)
    
    start = 0
    end = Nk_in_seg[0]
    for ii in range(Nseg):
        if ii != 0:
            start += Nk_in_seg[ii-1]
            end += Nk_in_seg[ii]
        vkpt_diff[start:end+1, :] = kpoints[start:end+1, :] - kpoints[start, :]

    kpt_path = np.linalg.norm(np.dot(vkpt_diff, B), axis=1)
    
    start = 0
    end = Nk_in_seg[0]
    for ii in range(1, Nseg):
        start += Nk_in_seg[ii-1]
        end += Nk_in_seg[ii]
        kpt_path[start:end] += kpt_path[start-1]

    #if the start and end of the kpath is the same, a jump occurs in the graph, this line solves it
    kpt_path[-1] = 2*kpt_path[-2] - kpt_path[-3] 
    
    kpt_bounds = kpt_path[ticks]
    

    
    
    fig = plt.figure()
    fig.set_size_inches(width, height)
    ax = plt.subplot(111)
    
    color = 'b'
    
    for band in np.transpose(bands):
        line, = ax.plot(kpt_path,band, lw = line_w, zorder=0, color=color)
    
    for t in kpt_bounds:
        ax.axvline(x=t, ls='-', color='k', lw=0.5, alpha=0.5)
    
    ax.set_xticks(kpt_bounds)
    ax.set_xticklabels(hsp)
    ax.set_ylim(ymin, ymax)
    
    plt.show()