import numpy as np
import matplotlib.pyplot as plt
import utilsQE



qpoints, bands = utilsQE.readPhononbandFreq(r"AgP2\phonons\AgP2.freq.gp")

qlabels, positions = utilsQE.readHighSymPointsPhonon(r"AgP2\phonons\matdyn_AgP2.in")

plt.figure()

for band in bands:
    plt.plot(qpoints, band, linewidth=1, alpha=1, color='k')

for pos in positions:
    plt.axvline(x=qpoints[pos], linewidth=0.5, color='k')

plt.xticks(ticks=qpoints[positions[:]], labels=qlabels)

plt.axhline(y=0, linewidth=0.5, color='b', linestyle='--')

plt.ylabel("Frequency (cm$^{-1}$)")
plt.xlim(qpoints[0], qpoints[-1])
#plt.ylim(0, )
plt.show()