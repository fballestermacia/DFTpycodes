import matplotlib.pyplot as plt
import numpy as np
import utilsQE
                
    
if __name__ == '__main__':
    
    fig_size = (15,8)
    line_w = 1
    ymin = 4
    ymax = 18
    
    kpoints, bands = utilsQE.read_filband('QE_projects\\AgP2\\bands_norelaxed\\filbandAgP2')
    kpoints2, nlines, hsp, B, pointsinline = utilsQE.read_bandsin('QE_projects\\AgP2\\bands_norelaxed\\AgP2_bands.in')
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