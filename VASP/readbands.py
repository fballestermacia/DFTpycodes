import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib as mpl
import re
from matplotlib.ticker import AutoMinorLocator
import utilsVASP


rc('text', usetex=True)
rc('font', size=15)
rc('legend', fontsize=13)
rc('text.latex')
        


utilsVASP.plot_bands_from_VASP(outcarfile="Ag2Te\MBJ_thirdtry\OUTCAR", kpointsfile="Ag2Te\MBJ_thirdtry\KPOINTS",fig_name = "Ag2Te\\bands_mbjthirdtry_newcode.pdf", ewindow=[-1,1], filledvsempty=True)