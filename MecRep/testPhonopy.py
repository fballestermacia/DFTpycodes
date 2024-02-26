import sys
from phonopy.interface.qe import read_pwscf, PH_Q2R
import numpy as np
import spglib as spg
import ase.io as io
import ase as ase
import ase.spacegroup as asespg

import phonopy as ph

cell, _ = read_pwscf(r'AgP2\Phonons_V2\AgP2.scf.pwi')
print(cell)
q2r = PH_Q2R(r'AgP2\Phonons_V2\AgP2.fc')
q2r.run(cell, is_full_fc=True)
q2r.write_force_constants(fc_format='text')


