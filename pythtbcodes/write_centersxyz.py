import cellconstructor as CC
import cellconstructor.Phonons
import numpy as np
import datetime

dyn = CC.Phonons.Phonons()
dyn.LoadFromQE(fildyn_prefix="data/AgP2/Phonons/444/dynmats/AgP2.dyn", nqirr=30)

positions = dyn.structure.coords

with open("DFTpycodes/pythtbcodes/AgP2_phononTB_centres.xyz",'w') as f:
    f.write("    "+str(len(positions))+"\n")
    f.write("Wannier centres, written by Francesc Ballester on {}\n".format(datetime.datetime.now()))
    for pos in positions:
        for i in range(3):
            f.write('X         {}      {}      {}\n'.format(pos[0],pos[1],pos[2]))
    for i, label in enumerate(dyn.structure.atoms):
        f.write('{}         {}      {}      {}\n'.format(label,positions[i,0],positions[i,1],positions[i,2]))