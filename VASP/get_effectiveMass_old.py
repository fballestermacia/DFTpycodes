import numpy as np
import matplotlib.pyplot as plt
from utilsVASP import fromOUTCARtoplot
from scipy.integrate import simps
from numpy.polynomial.polynomial import polyder, polyval, polyfit


def directionalPolyFit(band,xs, ys, polyorder):
    
    
    polyxs = []#polyfit(xs,band[:,0], polyorder)
    polyys = []
    for i in range(len(ys)):
        polyxs.append(polyfit(xs,band[:,i], polyorder))
    for j in range(len(xs)):
        polyys.append(polyfit(ys,band[j], polyorder)) 
    
    
    
    return polyxs, polyys


def gausiana(x,center,sigma,a0):
    return a0*np.exp(-(np.power(x-center,2))/(2*sigma**2))

if __name__ == '__main__':
    qz = np.array([-13.01,0,7.48])
    qz = qz/np.linalg.norm(qz)
    qy = np.array([0,1,0])
    qx = np.cross(qy, qz/np.linalg.norm(qz))#np.array([0.5,0,0])
    
    bands, efermi, qxvals, qyvals, occupation = fromOUTCARtoplot(outcarfile='Ag2Te\\ArchivosMartin\\GM_MBJ\\OUTCAR_MBJ',kpointsfile="Ag2Te\\ArchivosMartin\\GM_MBJ\\KPOINTS_MBJ_BS", qx=qx, qy=qy, weightfilter = 0)
    #bands, efermi, qxvals, qyvals, occupation = fromOUTCARtoplot(outcarfile='Ag2Te\ArchivosMartin\GM_MBJ\OUTCAR_MBJ',kpointsfile="Ag2Te\ArchivosMartin\GM_MBJ\KPOINTS_MBJ_BS", qx=qx, qy=qy, weightfilter = 0)

    dummy = np.transpose(bands[0])
    dummyoc = np.transpose(occupation[0])

    index = np.sum(dummyoc, axis=0)



    band = ((dummy)[int(index[0])])*0.0367493#((dummy)[int(index[0]+1)]-efermi)/0.0367493 #in hartree#*1000 #in meV
    
    #band -= np.min(band)
    polyorder = 8
    
    scaling = 1/1.88973 #Bohr
    qxvals *= scaling
    qyvals *= scaling
    
    xs, ylen = np.unique(qxvals, return_counts=True) 
    ys, xlen = np.unique(qyvals, return_counts=True) #!!!!Notice the change
    
    
    band2 = band.reshape(xlen[0], ylen[0])
    
    
    polyxs, polyys = directionalPolyFit(band2,xs, ys, polyorder)
    
    derx = polyder(polyxs,2, axis=1)
    dery = polyder(polyys,2, axis=1)
    
    
    derxvals = polyval(xs, np.transpose(derx), tensor=True)
    deryvals = np.transpose(polyval(ys, np.transpose(dery), tensor=True))
    
    '''trues = np.where(band<0.285)
    print(qxvals[trues])'''

    '''derxvals = []
    deryvals = []
    
    for dx in derx:
        derxvals.append(np.array(polyval(xs,dx)))
        
    for dy in dery:
        deryvals.append(np.array(polyval(ys,dy)))
    
    derxvals = np.array(derxvals)
    deryvals = np.array(deryvals) '''
    
    '''derxvals = np.where(derxvals<0, 0, derxvals)
    deryvals = np.where(deryvals<0, 0, deryvals)'''
    
    factor = 1#1000*(6.582e-16)**2/0.510e6*3e8**2/1e-10**2*4*np.pi**2
    
    
    mx = np.power(derxvals,-1)*factor
    my = np.power(deryvals,-1)*factor
    

    mx = np.where(mx<0, 0, mx)
    #mx = np.where(mx>0.5, 0.5, mx)
    
    #my = np.where(my<0, 0, my)
    #my = np.where(my>0.4, 0.4, my)
    
    
    efermi = 26.821041107177734/1000
    fermitohartree = 0.0367493
    enval = (efermi*0.0367493+np.min(band))/fermitohartree#0.283
    sigmax = 0.01
    sigmay = 0.01
    a0x=1/(sigmax*np.sqrt(2*np.pi))
    a0y=1/(sigmay*np.sqrt(2*np.pi))

    envals = np.linspace(np.min(band/fermitohartree), np.min(band/fermitohartree)+0.2, 2000)#np.min(band)+0.01, 1000)
    
    bandgaussianx = gausiana(band/fermitohartree,enval,sigmax,a0x)
    bandgaussiany = gausiana(band/fermitohartree,enval,sigmay,a0y)
    
    
    mxeff = np.power(np.sum(derxvals.flatten()*bandgaussianx.flatten())/simps(bandgaussianx),-1)#np.sum(mx*bandgaussianx.reshape([len(mx[:,0]),len(mx[0])]))/simps(bandgaussianx)
    myeff = np.power(np.sum(deryvals.flatten()*bandgaussiany.flatten())/simps(bandgaussiany),-1)#np.sum(my*bandgaussiany.reshape([len(my[:,0]),len(my[0])]))/simps(bandgaussiany)
    
    ratio = np.true_divide(mxeff,myeff, out=np.zeros_like(mxeff), where=myeff!=0)
    
    print("m_effX = ", mxeff)
    print("m_effY = ", myeff)
    print("ratio = ", ratio)  
    
    
    mxeffs = []
    myeffs = []
    
    for enval2 in envals:
        bandgaussian = gausiana(band/fermitohartree,enval2,sigmax,a0x)
        mxeffs.append(np.power(np.sum(np.transpose(derxvals)*bandgaussian.reshape([61,31]))/simps(bandgaussian),-1))#np.sum(mx*bandgaussian.reshape([len(mx[:,0]),len(mx[0])]))/simps(bandgaussian))
        
        bandgaussian = gausiana(band/fermitohartree,enval2,sigmay,a0y)
        myeffs.append(np.power(np.sum(deryvals*bandgaussian.reshape([31,61]))/simps(bandgaussian),-1))#np.sum(my*bandgaussian.reshape([len(my[:,0]),len(my[0])]))/simps(bandgaussian))
    

    
    plt.figure()
    plt.plot(envals-np.min(band/fermitohartree), mxeffs)
    plt.plot(envals-np.min(band/fermitohartree),myeffs)
    plt.vlines(enval-np.min(band/fermitohartree), 0, 0.25, 'r')

    plt.figure()
    plt.plot(envals-np.min(band/fermitohartree),np.array(mxeffs)/np.array(myeffs))
    
    '''dxticks = np.array([0,5e3,1e4,1.5e4,2e4,2.5e4,3e4])
    dyticks = np.array([0,2e4,4e4,6e4,8e4])
    
    mxticks = np.array([0,0.1,0.2,0.3,0.4,0.5])
    myticks = np.array([0,0.5,0.1,0.15,0.2,0.25,0.3,0.35,0.4])
    
    '''
    enval=enval*fermitohartree
    X,Y = np.meshgrid(xs,ys)
    
    nlev = 20
    
    fig = plt.figure()
    fig.set_size_inches(15, 8)
    ax = plt.subplot(111)
    drawing = ax.tricontourf(qxvals,qyvals,band, levels = nlev)
    colorbar = fig.colorbar(drawing, ax = ax)
    ax.tricontour(qxvals,qyvals,band, colors='k', levels = nlev)
    ax.scatter(qxvals,qyvals,s=3,c='k')
    ax.scatter(0,0,s=10,c='r')

    fermilevel = ax.tricontour(qxvals, qyvals, band, [enval], colors='r', zorder=10)
    
    
    fig = plt.figure()
    fig.set_size_inches(15, 8)
    ax = plt.subplot(111)
    plt.title("d2E/dKx2")
    drawing = ax.contourf(X,Y,derxvals, levels = np.linspace(0, np.max(derxvals), 10*nlev), extend = 'both')
    colorbar = fig.colorbar(drawing, ax = ax)#, ticks=dxticks)
    ax.contour(X,Y,derxvals, colors='k', levels = nlev)
    ax.scatter(qxvals,qyvals,s=3,c='k')
    ax.scatter(0,0,s=10,c='r')
    
    fermilevel = ax.tricontour(qxvals, qyvals, band, [enval], colors='r', zorder=100)
    
    fig = plt.figure()
    fig.set_size_inches(15, 8)
    ax = plt.subplot(111)
    plt.title("d2E/dKy2")
    drawing = ax.contourf(X,Y,deryvals, levels = np.linspace(0, np.max(deryvals), 10*nlev), extend = 'both')
    colorbar = plt.colorbar(drawing, ax = ax)#, ticks=dyticks)
    ax.contour(X,Y,deryvals, colors='k', levels = nlev)
    ax.scatter(qxvals,qyvals,s=3,c='k')
    ax.scatter(0,0,s=10,c='r')
    
    fermilevel = ax.tricontour(qxvals, qyvals, band, [enval], colors='r', zorder=10)
    
    fig = plt.figure()
    fig.set_size_inches(15, 8)
    ax = plt.subplot(111)
    plt.title("mx")
    drawing = ax.contourf(X,Y,mx, levels = np.linspace(0, np.max(mx), 5*nlev))
    colorbar = fig.colorbar(drawing, ax = ax)#,ticks=mxticks)
    ax.contour(X,Y,mx, colors='k', levels = nlev)
    ax.scatter(qxvals,qyvals,s=3,c='k')
    ax.scatter(0,0,s=10,c='r')
    fermilevel = ax.tricontour(qxvals, qyvals, band, [enval], colors='r', zorder=10)
    
    fig = plt.figure()
    fig.set_size_inches(15, 8)
    ax = plt.subplot(111)
    plt.title("my")
    drawing = ax.contourf(X,Y,my, levels = np.linspace(0, np.max(my), 5*nlev), extend = 'both')
    colorbar = fig.colorbar(drawing, ax = ax)#,ticks=mxticks)
    ax.contour(X,Y,my, colors='k', levels = nlev)
    ax.scatter(qxvals,qyvals,s=3,c='k')
    ax.scatter(0,0,s=10,c='r')
    
    fermilevel = ax.tricontour(qxvals, qyvals, band, [enval], colors='r', zorder=10)

    
    fig = plt.figure()
    fig.set_size_inches(15, 8)
    ax = plt.subplot(111)
    plt.title("ratio")
    drawing = ax.contourf(X,Y,mx/my, levels = np.linspace(0, np.max(mx/my), 5*nlev), extend = 'both')
    colorbar = fig.colorbar(drawing, ax = ax)
    ax.contour(X,Y,mx/my, colors='k', levels = nlev)
    ax.scatter(qxvals,qyvals,s=3,c='k')
    ax.scatter(0,0,s=10,c='r')
    
    fermilevel = ax.tricontour(qxvals, qyvals, band, [enval], colors='r', zorder=10)
    
    plt.show()