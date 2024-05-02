import numpy as np
import matplotlib.pyplot as plt
import utilsQE

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
        if a == '\\Gamma':
            a = 'G'
        if b == '\\Gamma':
            b = 'G'
        seg, bands = readbulkek(bulkekfile=bulkekfile + '-' + a +'-' +b, linestart=3)
        segments.append(seg)
        bandspersegment.append(bands)
    return np.array(segments), np.array(bandspersegment)


factor = 0.123983#*0.24180
thztomev = 4.15665538536
qpoints, bands = utilsQE.readPhononbandFreq(r"data/AgP2/Phonons/444/AgP2.newHSP.freq.gp")

qlabels, positions = utilsQE.readHighSymPointsPhonon(r"data/AgP2/Phonons/444/matdyn.newHSP.in")


bands *= factor


tbqs, tbbandsperline = readpathbulkeks(qlabels, 'DFTpycodes/WTProj/bulkek.dat')

tbbandsperline *= thztomev



#resize the qpoints of the TB to fit those of QE
for i in range(len(positions)-1):
    tbqs[i] = tbqs[i]*(qpoints[positions[i+1]]-qpoints[positions[i]])/tbqs[i,-1] + qpoints[positions[i]]




plt.figure()
topocolors = 'k'*36#'g'*10+'g'*4+'k'*12+'g'*6+'k'*4


for i,band in enumerate(bands):
    plt.plot(qpoints, band, linewidth=1, alpha=1, color=topocolors[i], label='QE')


for i, tbands in enumerate(tbbandsperline):
    for j, tband in enumerate(tbands):
        plt.plot(tbqs[i], tband,linewidth = 1, alpha = 1, color = 'r', label='WannierTools')

for pos in positions:
    plt.axvline(x=qpoints[pos], linewidth=0.5, color='k')

plt.xticks(ticks=qpoints[positions[:]], labels=qlabels)

plt.axhline(y=0, linewidth=0.5, color='b', linestyle='--')

handles, labels = plt.gca().get_legend_handles_labels()
labels, ids = np.unique(labels, return_index=True)
handles = [handles[i] for i in ids]

plt.legend(handles, labels, loc = 'best')

plt.ylabel("Frequency (meV)" )#(cm$^{-1}$)")
plt.xlim(qpoints[0], qpoints[-1])

plt.title('phonon bands')

#plt.ylim(0, )
plt.show()