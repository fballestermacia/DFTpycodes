#!/usr/bin/env python

# Copyright under GNU General Public License 2010, 2012, 2016
# by Sinisa Coh and David Vanderbilt (see gpl-pythtb.txt)

from pythtb import * # import TB model class
import matplotlib.pyplot as plt
import numpy as np

# read output from Wannier90 that should be in folder named "example_a"
#   see instructions above for how to obtain the example output from 
#   Wannier90 for testing purposes
agp2=w90(r"DFTpycodes/WTProj",r"AgP2_phononTB")



my_model=agp2.model(min_hopping_norm=0.01)

'''# solve model on a path and plot it
path=[[0.5,0.5,0.5],[0.5,0.0, 0.5],[0.0,0.0,0.0], [0.0,0.0,0.5], [0.0, 0.5, 0.5], [0.5,0.5,0.0], [0.0, 0.5, 0.0],[0.0,0.0,0.0], [0.5, 0.0, 0.0],[0.5,0.5,0.0], [0.5, 0.5, 0.5]]
# labels of the nodes
k_label=(r'$E$', r'$A$', r'$\Gamma$',r'$B$', r'$D$',r'$C$', r'$Z$', r'$\Gamma$', r'$Y$',r'$C$', r'$E$')
# call function k_path to construct the actual path
(k_vec,k_dist,k_node)=my_model.k_path(path,51)
#
evals, evecs=my_model.solve_all(k_vec, eig_vectors=True)
fig, ax = plt.subplots()
for i in range(evals.shape[0]):
    ax.plot(k_dist,np.sign(evals[i])*abs(evals[i])**0.5,"k-")
for n in range(len(k_node)):
    ax.axvline(x=k_node[n],linewidth=0.5, color='k')
ax.set_xlabel("Path in k-space")
ax.set_ylabel("Band energy (eV)")
ax.set_xlim(k_dist[0],k_dist[-1])
ax.set_xticks(k_node)
ax.set_xticklabels(k_label)
fig.tight_layout()

mode = 9
(fig2,ax2)=my_model.visualize(0,2,eig_dr=evecs[mode,0,:],draw_hoppings=False)
ax2.set_title("OABR state(?)")
ax2.set_xlabel("x coordinate")
ax2.set_ylabel("y coordinate")
fig2.tight_layout()
'''
# cutout finite model first along direction z with PBC along x and y
finite_model=my_model.cut_piece(1,2,glue_edgs=False)

# solve model on a path and plot it
path=[[0,0],[0.5,0.0],[0.5,0.5], [0.0,0.5], [0.0, 0.0], [0.5,0.5]]
# labels of the nodes
k_label=(r'$\Gamma$', r'$X$', r'$M$',r'$Y$', r'$\Gamma$',r'$M$')
# call function k_path to construct the actual path
(k_vec,k_dist,k_node)=finite_model.k_path(path,51)

evals, evecs=finite_model.solve_all(k_vec, eig_vectors=True)
fig, ax = plt.subplots()
for i in range(evals.shape[0]):
    ax.plot(k_dist,np.sign(evals[i])*abs(evals[i])**0.5,"k-")
for n in range(len(k_node)):
    ax.axvline(x=k_node[n],linewidth=0.5, color='k')
ax.set_xlabel("Path in k-space")
ax.set_ylabel("Band energy (eV)")
ax.set_xlim(k_dist[0],k_dist[-1])
ax.set_xticks(k_node)
ax.set_xticklabels(k_label)
fig.tight_layout()

# solve finite models
(evalsf,evecsf)=finite_model.solve_one(k_point=[0,0],eig_vectors=True)
edgemode = 25
otheredgemode = edgemode+1

'''for i,e in enumerate(evalsf):
    print(i,e)'''


(fig3,ax3)=finite_model.visualize(1,2,eig_dr=evecsf[edgemode,:],draw_hoppings=False)
ax3.set_title("Surface State (?)")
ax3.set_xlabel("x coordinate")
ax3.set_ylabel("y coordinate")
fig3.tight_layout()

(fig4,ax4)=finite_model.visualize(1,2,eig_dr=evecsf[otheredgemode,:],draw_hoppings=False)
ax4.set_title("Surface State (?)")
ax4.set_xlabel("x coordinate")
ax4.set_ylabel("y coordinate")
fig4.tight_layout()

plt.show()
#fig.savefig("silicon_quick.pdf")