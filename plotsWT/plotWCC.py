import numpy as np
import matplotlib.pyplot as plt

def readwcc(prefix):
    data = np.transpose(np.loadtxt(prefix+'.dat'))
    
    ks = data[0]
    largestgap = data[1]
    sumwcc = data[2]
    
    wcc = data[3:]

    return ks, largestgap, sumwcc, wcc


def readWeyl3D(prefix,nnodes):
    data = np.transpose(np.loadtxt(prefix+'.dat'))
    ks = data[0]
    weyls = data[1:nnodes]

    return ks, weyls


#k, lg, swcc, w = readwcc(r'DFTpycodes/WTProj/wcc')

k, w = readWeyl3D(r'data/Al2Te4Zn/444/wanniercenter3D_Weyl',8)


plt.figure()


colors = ['m','m','r','r','b','b','g','g']
ax1 = plt.axes()

for i in range(len(w)):
    ax1.scatter(np.array(k), np.array(w[i]), c = colors[i], label='WannierTools',s=20)


plt.yticks([0,0.5,1],fontsize=20)
plt.ylabel("WCC",fontsize=20 )#(cm$^{-1}$)")
plt.xlim(np.min(k), np.max(k))
plt.xticks(ticks=[np.min(k),(np.max(k)-np.min(k))/2, np.max(k)], labels=np.round([np.min(k),(np.max(k)-np.min(k))/2, np.max(k)],2),fontsize=20)

plt.title('Al$_2$Te$_4$Zn Weyl Chirality Calculation', fontsize=20)

#plt.ylim(35, 50)

plt.show()