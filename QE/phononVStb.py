import numpy as np
import matplotlib.pyplot as plt
import utilsQE
import cellconstructor as CC, cellconstructor.Phonons
import cellconstructor.ForceTensor

def readbulkek(bulkekfile = 'bulkek.dat', linestart = 2):
    filelines = [line for line in open(bulkekfile) if line.strip()]
    isqpointfull = False
    bands = []
    qpoints = []
    band = []
    prevq = -1
    for fline in filelines[linestart:]:
        #print(fline)
        if prevq > float(fline.split()[0]):
            isqpointfull = True
            bands.append(np.array(band))
            band = [float(fline.split()[1])]
            prevq = float(fline.split()[0])
            continue
        if not isqpointfull:
            qpoints.append(float(fline.split()[0]))
        prevq = float(fline.split()[0])
        band.append(float(fline.split()[1]))
    bands.append(np.array(band))
    return np.array(qpoints), np.array(bands)
        
def readpathbulkeks(pathlabels,bulkekfile = 'bulkek.dat'):
    segments = []
    bandspersegment = []
    for i in range(len(pathlabels)-1):
        a = pathlabels[i].split('$')[1].split('_')
        b = pathlabels[i+1].split('$')[1].split('_')
        if len(a)>1:
            a = a[0] + a[1]
        else:
            a = a[0]
        if len(b)>1:
            b = b[0] + b[1]
        else:
            b = b[0]
        if a == '\\Gamma' or a == 'Gamma' or a.lower() == 'gamma':
            a = 'G'
        if b == '\\Gamma' or b == 'Gamma' or b.lower() == 'gamma':
            b = 'G'
        seg, bands = readbulkek(bulkekfile=bulkekfile + '-' + a +'-' +b, linestart=3)
        segments.append(seg)
        bandspersegment.append(bands)
    return np.array(segments), np.array(bandspersegment)

factor = 0.123983#*0.24180
thztomev = 4.15665538536

# Let us define the PATH in the brilluin zone and the total number of points
PATH = "XWLGXBKUG"
N_POINTS = 1000
# Here we define the position of the special points
SPECIAL_POINTS = {"G": [0,0,0],"W": [.5, .25, .75],"X": [.5, 0, .5],"L": [.5, .5, .5],"B": [.5, .5, 1],"K": [.375, .375, .75],"U": [.125, .125, .25]}

SSCHA_DYN = r'data/Al2Te4Zn/444/Al2Te4Zn1.dyn'
NQIRR2 = 14

sscha_dyn = CC.Phonons.Phonons(SSCHA_DYN, NQIRR2)
qpath, data = CC.Methods.get_bandpath(sscha_dyn.structure.unit_cell,PATH,SPECIAL_POINTS,N_POINTS)
xaxis, xticks, xlabels = data # Info to plot correclty the x axis
sscha_dispersion = CC.ForceTensor.get_phonons_in_qpath(sscha_dyn, qpath)*factor
nmodes = sscha_dyn.structure.N_atoms * 3



qpoints, bands = utilsQE.readPhononbandFreq(r"data/Al2Te4Zn/444/Al2Te4Zn1.freq.gp")

qlabels, positions = utilsQE.readHighSymPointsPhonon(r"data/Al2Te4Zn/444/matdyn.in")



bands *= factor


tbqs, tbbandsperline = readpathbulkeks(qlabels, r'data/Al2Te4Zn/444/bulkek.dat')

tbbandsperline *= thztomev



#resize the qpoints of the TB to fit those of QE
for i in range(len(positions)-1):
    tbqs[i] = tbqs[i]*(qpoints[positions[i+1]]-qpoints[positions[i]])/tbqs[i,-1] + qpoints[positions[i]]

newxaxis = np.array([])

for i in range(len(positions)-1):
    if i == 0:
        mask = ((xaxis[:]>=xticks[i]) & (xaxis[:]<=xticks[i+1]))
        xmasked = xaxis[mask]
        prevmask = xmasked
        newxaxis = np.append(newxaxis,np.array(((xmasked-xmasked[0])/(xmasked[-1]-xmasked[0]))*(qpoints[positions[i+1]]-qpoints[positions[i]]) + qpoints[positions[i]]).flatten())
        
    else: 
        mask = ((xaxis[:]>xticks[i]) & (xaxis[:]<=xticks[i+1]))
        xmasked = np.append(prevmask[-1],xaxis[mask])
        prevmask = xmasked
        newpoints = (xmasked-xmasked[0])/(xmasked[-1]-xmasked[0])*(qpoints[positions[i+1]]-qpoints[positions[i]]) + qpoints[positions[i]]
        newxaxis = np.append(newxaxis,np.array(newpoints).flatten()[1:])
        
    newxaxis = np.array(newxaxis).flatten()
    #xaxis2[i] = xaxis2[i]*(qpoints[positions[i+1]]-qpoints[positions[i]])/xaxis2[i,-1] + qpoints[positions[i]]

xaxis = np.array(newxaxis).flatten()


plt.figure()
topocolors = 'k'*36#'g'*10+'g'*4+'k'*12+'g'*6+'k'*4
ax1 = plt.axes()

for i,band in enumerate(bands):
    ax1.plot(qpoints, band, linewidth=1, alpha=1, color=topocolors[i], label='QE')


for i, tbands in enumerate(tbbandsperline):
    for j, tband in enumerate(tbands):
        ax1.plot(tbqs[i], tband,linewidth = 1, alpha = 1, color = 'r', label='WannierTools')

for pos in positions:
    ax1.axvline(x=qpoints[pos], linewidth=0.5, color='k')

for i in range(nmodes):
    ax1.plot(xaxis, sscha_dispersion[:,i], color = 'g', label = 'cellconstructor', linestyle='-.')

plt.xticks(ticks=qpoints[positions[:]], labels=qlabels)

ax1.axhline(y=0, linewidth=0.5, color='b', linestyle='--')

handles, labels = plt.gca().get_legend_handles_labels()
labels, ids = np.unique(labels, return_index=True)
handles = [handles[i] for i in ids]

ax1.legend(handles, labels, loc = 'best')

plt.ylabel("Frequency (meV)" )#(cm$^{-1}$)")
plt.xlim(qpoints[0], qpoints[-1])

plt.title('phonon bands')

#plt.ylim(0, )
plt.show()