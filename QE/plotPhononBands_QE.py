import numpy as np
import matplotlib.pyplot as plt
import utilsQE


factor = 0.123983#*0.24180
qpoints, bands = utilsQE.readPhononbandFreq(r"data/232/AgP2.freq.gp")

qlabels, positions = utilsQE.readHighSymPointsPhonon(r"data/232/matdyn.in")

bands *= factor

plt.figure()
topocolors = 'k'*36#'g'*10+'g'*4+'k'*12+'g'*6+'k'*4


for i,band in enumerate(bands):
    plt.plot(qpoints, band, linewidth=1, alpha=1, color=topocolors[i])

for pos in positions:
    plt.axvline(x=qpoints[pos], linewidth=0.5, color='k')

plt.xticks(ticks=qpoints[positions[:]], labels=qlabels)

plt.axhline(y=0, linewidth=0.5, color='b', linestyle='--')

plt.ylabel("Frequency (meV)" )#(cm$^{-1}$)")
plt.xlim(qpoints[0], qpoints[-1])
#plt.ylim(0, )
plt.show()