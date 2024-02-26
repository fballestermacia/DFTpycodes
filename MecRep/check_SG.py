import ase.io as io
from ase.spacegroup import get_spacegroup

structure = io.read(r'AgP2\Phonons_V2\AgP2_scf.in')
sp = get_spacegroup(structure,symprec=1e-7)
print(sp)
