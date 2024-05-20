# Import the CellConstructor library to plot the dispersion
import cellconstructor as CC, cellconstructor.Phonons
import cellconstructor.ForceTensor
# Import the numerical libraries and those for plotting
import numpy as np
import matplotlib.pyplot as plt
import sys, os
# Let us define the PATH in the brilluin zone and the total number of points
PATH = "EAGBDCZGYCE"
N_POINTS = 1000
# Here we define the position of the special points
SPECIAL_POINTS = {"G": [0,0,0],"E": [.5, .5, .5],"A": [.5, 0, .5],"B": [0, 0, .5],"D": [0, .5, .5],"C": [.5, .5, 0],"Z": [0, .5, 0],"Y": [.5, 0, 0]}
# The two dynamical matrix to be compared
HARM_DYN = 'data/AgP2/Phonons/444/dynmatnoloto/AgP2.dyn'
SSCHA_DYN = 'data/AgP2/Phonons/444/dynmats/AgP2.dyn'
# The number of irreducible q points
# i.e., the number of files in which the phonons are stored
NQIRR1 = 30
NQIRR2 = 30

# --------------------- THE SCRIPT FOLLOWS ---------------------
# Load the harmonic and sscha phonons
harmonic_dyn = CC.Phonons.Phonons(HARM_DYN, NQIRR1)
sscha_dyn = CC.Phonons.Phonons(SSCHA_DYN, NQIRR2)
# Get the band path
qpath, data = CC.Methods.get_bandpath(harmonic_dyn.structure.unit_cell,PATH,SPECIAL_POINTS,N_POINTS)
xaxis, xticks, xlabels = data # Info to plot correclty the x axis
# Get the phonon dispersion along the path
harmonic_dispersion = CC.ForceTensor.get_phonons_in_qpath(harmonic_dyn, qpath)
sscha_dispersion = CC.ForceTensor.get_phonons_in_qpath(sscha_dyn, qpath)
nmodes = harmonic_dyn.structure.N_atoms * 3
# Plot the two dispersions
plt.figure(dpi = 150)
ax = plt.axes()
for i in range(nmodes):
    lbl=None
    lblsscha = None
    if i == 0:
        lbl = 'no LOTO'
        lblsscha = 'LOTO'

    ax.plot(xaxis, harmonic_dispersion[:,i], color = 'k', ls = 'dashed',label = lbl)
    ax.plot(xaxis, sscha_dispersion[:,i], color = 'r', label = lblsscha)
# Plot vertical lines for each high symmetry points
for x in xticks:
    ax.axvline(x, 0, 1, color = "k", lw = 0.4)
ax.axhline(0, 0, 1, color = 'k', ls = ':', lw = 0.4)
# Set the x labels to the high symmetry points
ax.set_xticks(xticks)
ax.set_xticklabels(xlabels)
ax.set_xlabel("Q path")
ax.set_ylabel("Phonons [cm-1]")
plt.tight_layout()
plt.show()