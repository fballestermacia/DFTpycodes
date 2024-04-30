import ase.io as io
from ase.spacegroup import get_spacegroup
import sys


scfinfile = str(sys.argv[1])
structure = io.read(scfinfile)
sp = get_spacegroup(structure,symprec=1e-7)
print(sp)
