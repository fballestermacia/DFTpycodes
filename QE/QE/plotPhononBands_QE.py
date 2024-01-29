import numpy as np
import matplotlib.pyplot as plt


def readPhononbandFreq(freqgpfile): #freqgpfile=SYSTEM.freq.gp

    data = np.loadtxt(freqgpfile)

    nbands = data.shape[1] - 1
    qpoints = data[:, 0]
    bands = np.transpose(data)[1:,]

    return qpoints, bands


def readHighSymPointsPhonon(matdynfile):  #matdynfile=matdyn.in
    filelines = [line for line in open(matdynfile) if line.strip()]
    qbandform=False
    for i,fline in enumerate(filelines):
        
        if fline.split()[0]=='readtau' and fline.split()[-1]=='.true.':
            continue #TODO: INCLUDE THE CASE WHERE ATOMIC POSITIONS ARE READ

        if fline.split()[0]=='q_in_band_form' and fline.split()[-1]=='.true.':
            qbandform=True

        if fline.split()[0]=='/':
            if qbandform:
                start = i+1

            #TODO: make the case where qbandform=False
    
    npoints = int(filelines[start])
    labels = []
    positions = np.empty(npoints, dtype='int')
    counter = 0
    for i,kline in enumerate(filelines[start+1:]):
        dummylabel = str(kline.split('!')[-1].split()[0])
        if dummylabel[0].upper() == 'G':
            labels.append('$\\Gamma$')  
        else: labels.append('${}$'.format(dummylabel))
        positions[i] = counter
        counter += int(kline.split()[-2])
    return labels, positions






qpoints, bands = readPhononbandFreq("AgP2\phonons\AgP2.freq.gp")

qlabels, positions = readHighSymPointsPhonon("AgP2\phonons\matdyn_AgP2.in")

plt.figure()

for band in bands:
    plt.plot(qpoints, band, linewidth=1, alpha=1, color='k')

for pos in positions:
    plt.axvline(x=qpoints[pos], linewidth=0.5, color='k')
print(qlabels)
plt.xticks(ticks=qpoints[positions[:]], labels=qlabels)

plt.axhline(y=0, linewidth=0.5, color='b', linestyle='--')

plt.ylabel("Frequency (cm$^{-1}$)")
plt.xlim(qpoints[0], qpoints[-1])
#plt.ylim(0, )
plt.show()