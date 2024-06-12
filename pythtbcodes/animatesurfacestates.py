import matplotlib.pyplot as plt
import numpy as np


#REMEMBER PROJECTIONS ARE ORDERES AS PZ. PX, PY

positions = np.loadtxt('DFTpycodes/pythtbcodes/positions441.txt')

mode = np.loadtxt('DFTpycodes/pythtbcodes/surfacemode441.txt',dtype='complex')

print(np.shape(mode), np.shape(np.unique(positions,axis=0)))



#_, mask = np.unique(positions,axis=0, return_index=True)
newpos = positions[::3]#np.sort(mask)]
natms = len(newpos)
polvec = np.real(mode)

polvec=np.reshape(polvec,(natms,3))
a = np.arange(natms)/1000
b = np.zeros(natms)
fig1 = plt.figure()
ax = fig1.add_subplot(projection='3d')

ax.scatter(newpos[:,0],newpos[:,1],newpos[:,2],marker='o',linewidths=1, c='r')

#ax.quiver(newpos[:,0],newpos[:,1],newpos[:,2], polvec[:],b,b, length=0.5)

#ax.quiver(newpos[:,0],newpos[:,1],newpos[:,2], polvec[natms:2*natms],polvec[2*natms:],polvec[0:natms], length=0.5)
ax.quiver(newpos[:,0],newpos[:,1],newpos[:,2], polvec[:,1],polvec[:,2],polvec[:,0], length=0.1)

ax.set_title("Surface mode, real part")
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
#ax.set_zlim(0,25)
plt.show()