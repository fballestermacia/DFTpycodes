import numpy as np
import matplotlib.pyplot as plt

def readNodes(prefix):
    data = np.transpose(np.loadtxt(prefix+'.dat'))
    kcart = np.array([data[0],data[1],data[2]])
    kdirect = np.array([data[-3],data[-2],data[-1]])

    return kcart, kdirect


kc, kd = readNodes(r'data/Al2Te4Zn/444/Nodes')

colors = ['m','m','r','r','b','b','g','g']
charge = ['-1','-1','+1','+1','-1','-1','+1','+1']

offset = [0.05,0.05,0.05]

fig1 = plt.figure(figsize=(8,8))
ax = fig1.add_subplot(projection='3d')




ax.scatter(kc[0],kc[1],kc[2],s=200, c=colors[:])

for i in range(len(kc[0])):
    ax.text(kc[0,i]+offset[0],kc[1,i]+offset[1],kc[2,i]+offset[2],charge[i],color=colors[i],fontsize=20)

ax.set_xticks([-0.5,0,0.5],[-0.5,0,0.5],fontsize=20)
ax.set_xlabel("k$_x$",fontsize=20 )

ax.set_yticks([-0.5,0,0.5],[-0.5,0,0.5],fontsize=20)
ax.set_ylabel("k$_y$",fontsize=20 )

ax.set_zticks([-0.5,0,0.5],[-0.5,0,0.5],fontsize=20)
ax.set_zlabel("k$_z$",fontsize=20 )

plt.title('Al$_2$Te$_4$Zn Weyl Nodes', fontsize=20)

plt.show()