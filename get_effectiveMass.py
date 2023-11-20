import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from readKgrid import fromOUTCARtoplot
from scipy.integrate import simps
from numpy.polynomial.polynomial import polyder, polyval, polyfit
from scipy.interpolate import griddata

def directionalPolyFit(band,xs, ys, polyorder):
    
    
    polyxs = []
    polyys = []
    for i in range(len(ys)):
        polyxs.append(polyfit(xs,band[:,i], polyorder))
    for j in range(len(xs)):
        polyys.append(polyfit(ys,band[j], polyorder)) 
    
    
    
    return polyxs, polyys


def gausiana(x,center,sigma):
    a0 = 1/(sigma*np.sqrt(2*np.pi))
    return a0*np.exp(-(np.power(x-center,2))/(2*sigma**2))

if __name__ == '__main__':
    #CONSTANTS
    eV2Hartree = 1/27.2113845#0.0367493
    polyorder = 8
    sigx = 0.0002
    sigy = sigx
    #########


    qz = np.array([-13.01,0,7.48])
    qz = qz/np.linalg.norm(qz)

    qy = np.array([0,1,0])

    qx = np.cross(qy, qz/np.linalg.norm(qz))

    bands, efermi, qxvals, qyvals, occupation = fromOUTCARtoplot(outcarfile='Ag2Te\\ArchivosMartin\\GM_MBJ\\OUTCAR_MBJ',kpointsfile="Ag2Te\\ArchivosMartin\\GM_MBJ\\KPOINTS_MBJ_BS", qx=qx, qy=qy, weightfilter = 0)
    
    dummy = np.transpose(bands[0])
    dummyoc = np.transpose(occupation[0])

    index = np.sum(dummyoc, axis=0)

    band = ((dummy)[int(index[0])])
    band = (band-np.min(band))*eV2Hartree

    
    

    xs, ylen = np.unique(qxvals, return_counts=True) 
    ys, xlen = np.unique(qyvals, return_counts=True)
    

    shape = [ylen[0], xlen[0]]

    band2 = band.reshape(shape[::-1])
    polyxs, polyys = directionalPolyFit(band2,xs, ys, polyorder)

    derx = polyder(polyxs,2, axis=1)
    dery = polyder(polyys,2, axis=1)

    derxvals = polyval(xs, np.transpose(derx), tensor=True)
    deryvals = polyval(ys, np.transpose(dery), tensor=True)

    #print(derxvals)

    efermi = 26.821041107177734/1000 #26 meVs
    
    envals = np.linspace(np.min(band2), np.min(band2)+0.15*eV2Hartree, 200)

    
    xx = np.linspace(np.min(xs), np.max(xs), 61)
    yy = np.linspace(np.min(ys), np.max(ys), 31)

    X,Y = np.meshgrid(xx,yy)

    gridddx = griddata(np.transpose([qxvals,qyvals]),np.transpose(derxvals).flatten(),(X,Y), method='cubic')
    gridddy = griddata(np.transpose([qxvals,qyvals]),deryvals.flatten(),(X,Y), method='cubic')
    gridE = griddata(np.transpose([qxvals,qyvals]),band2.flatten(),(X,Y), method='cubic')

    for en in envals:
        #Cicle through all data in the grgrid    x_idev2=0
        x_idev2=0
        x_weight=0
        y_idev2=0
        y_weight=0
        for i in range(np.shape(gridE)[0]):
            for j in range(np.shape(gridE)[1]):
                #sum velocities
                x_idev2=x_idev2+gausiana(gridE[i,j], en, sigx)*gridddx[i,j]
                y_idev2=y_idev2+gausiana(gridE[i,j], en, sigy)*gridddy[i,j]
                #count summed terms
                x_weight=x_weight+gausiana(gridE[i,j], en, sigx)
                y_weight=y_weight+gausiana(gridE[i,j], en, sigy)
        try:
            idevs=np.vstack([idevs,np.array([en,x_idev2/x_weight,y_idev2/y_weight])])
        except NameError:
            idevs=np.array([en,x_idev2/x_weight,y_idev2/y_weight])

    Energies, ddxvals, ddyvals = np.transpose(idevs)

    x_idev2=0
    x_weight=0
    y_idev2=0
    y_weight=0
    for i in range(np.shape(gridE)[0]):
        for j in range(np.shape(gridE)[1]):
            #sum velocities
            x_idev2=x_idev2+gausiana(gridE[i,j], efermi*eV2Hartree, sigx)*gridddx[i,j]
            y_idev2=y_idev2+gausiana(gridE[i,j], efermi*eV2Hartree, sigy)*gridddy[i,j]
            #count summed terms
            x_weight=x_weight+gausiana(gridE[i,j], efermi*eV2Hartree, sigx)
            y_weight=y_weight+gausiana(gridE[i,j], efermi*eV2Hartree, sigy)
    try:
        idevsf=np.vstack([idevsf,np.array([en,x_idev2/x_weight,y_idev2/y_weight])])
    except NameError:
        idevsf=np.array([en,x_idev2/x_weight,y_idev2/y_weight])
    Energiesf, ddxvalsf, ddyvalsf = np.transpose(idevsf)
    
    print(ddxvalsf**-1,ddyvalsf**-1)

    ens = Energies/eV2Hartree
    mxeffs = ddxvals**-1
    myeffs = ddyvals**-1

    plt.figure()
    plt.plot(ens, mxeffs)
    plt.plot(ens, myeffs)
    plt.vlines(efermi, 0, 0.35, 'r')

    plt.figure()
    plt.plot(envals,np.array(mxeffs)/np.array(myeffs))
    
    band2 = band2/eV2Hartree*1000

    enval=efermi*eV2Hartree
    X,Y = qxvals,qyvals
    cmap = mpl.cm.viridis
    norm = mpl.colors.Normalize(vmin=np.min(band2), vmax= np.max(band2))


    nlev = 10
    
    '''fig = plt.figure()
    fig.set_size_inches(15, 8)
    ax = plt.subplot(111)
    drawing = ax.imshow(band2,interpolation='bilinear',cmap=cmap,norm=norm,
                        extent=(np.min(qxvals),np.max(qxvals),np.min(qyvals),np.max(qyvals)))#ax.tricontourf(qxvals,qyvals,band, levels = nlev)
    colorbar = fig.colorbar(drawing, ax = ax)
    ax.tricontour(qxvals,qyvals,band, colors='k', levels = nlev)
    ax.scatter(qxvals,qyvals,s=3,c='k')
    ax.scatter(0,0,s=10,c='r')

    fermilevel = ax.tricontour(qxvals, qyvals, band, [enval], colors='r', zorder=10)'''

    plt.show()