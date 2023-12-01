import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from readKgrid import fromOUTCARtoplot
from scipy.integrate import simps
from numpy.polynomial.polynomial import polyder, polyval, polyfit
from scipy.interpolate import griddata
import scipy.constants as cte
from get_effectiveMass import gausiana

def secondDer_FiniteDiff(fvals, xvals, boundary = 'closed'): #boundary = zero, closed

    dervdummy = np.empty(len(fvals))

    num = fvals[2:] - 2*fvals[1:-1] + fvals[:-2]
    den = (xvals[2:] - xvals[1:-1])*(xvals[1:-1]-xvals[:-2])

    dervdummy[1:-1] = num/den

    if boundary == 'zero':
        dervdummy[0] = 0
        dervdummy[-1] = 0
    elif boundary == 'closed':
        dervdummy[0] = dervdummy[1]
        dervdummy[-1] = dervdummy[-2]
    
    return dervdummy



if __name__ == '__main__':
    #CONSTANTS
    eV2Hartree = 1/27.2113845#0.0367493
    sigx = 0.0001
    sigy = sigx
    #########

    

    qz = np.array([-13.01,0,7.48])
    qz = qz/np.linalg.norm(qz)
    

    qy = np.array([0,1,0])

    qx = np.cross(qy, qz)

    '''qz = np.array([-0.1495968000,0,0.1495966000])
    qy = np.array([0,1,0])
    qx = np.array([0.1495966000,0,0.1495968000])

    qz /= np.linalg.norm(qz)
    qx /= np.linalg.norm(qx)'''

    bands, efermi, qxvals, qyvals, qzvals, occupation = fromOUTCARtoplot(outcarfile='Ag2Te\\PGHR_DENSER\\OUTCAR',kpointsfile="Ag2Te\\PGHR_DENSER\\KPOINTS",
                                                                          qx=qz, qy=qy, qz=qx, weightfilter = 0)
    
    dummy = np.transpose(bands[0])
    dummyoc = np.transpose(occupation[0])

    index = np.sum(dummyoc, axis=0)

    band = ((dummy)[int(index[0])])
    band = (band-np.min(band))*eV2Hartree

    
    #qyvals = qzvals

    
    xs, ylen = np.unique(qxvals, return_counts=True) 
    ys, xlen = np.unique(qyvals, return_counts=True)
    

    shape = [ylen[0], xlen[0]]

    band2 = band.reshape(shape[::-1])
    
    for xband in np.transpose(band2):
        dummyder = secondDer_FiniteDiff(xband, xs)
        try:
            derxvals=np.vstack([derxvals,np.array(dummyder)])
        except NameError:
            derxvals=np.array(dummyder)

    
    for yband in band2:
        dummyder = secondDer_FiniteDiff(yband, ys)
        try:
            deryvals=np.vstack([deryvals,np.array(dummyder)])
        except NameError:
            deryvals=np.array(dummyder)


    #print(derxvals)

    efermi_no10times = 42/1000 # meVs, fermi energi with previous doping
    efermi = 196/1000 # meVs, fermi energi with 10 times the doping
    envals = np.linspace(np.min(band2), np.min(band2)+0.25*eV2Hartree, 200)

    
    xx = np.linspace(np.min(xs), np.max(xs), len(xs))
    yy = np.linspace(np.min(ys), np.max(ys), len(ys))

    X,Y = np.meshgrid(xx,yy)

    gridddx = griddata(np.transpose([qxvals,qyvals]),np.transpose(derxvals).flatten(),(X,Y), method='linear')
    gridddy = griddata(np.transpose([qxvals,qyvals]),deryvals.flatten(),(X,Y), method='linear')
    gridE = griddata(np.transpose([qxvals,qyvals]),band2.flatten(),(X,Y), method='linear')

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
    
    
    me=9.10e-31
    alat=15.25 
    Ry2eV = 13.6056923
    Ry2jul = 2.179872e-18
    au2metre = 5.29177e-11
    factor = 1#(cte.hbar**2/me)*( (((2*np.pi)**2)*1000/Ry2eV)/(Ry2jul*(alat*au2metre)**2) )
    ens = Energies/eV2Hartree*1000
    mxeffs = ddxvals**-1*factor
    myeffs = ddyvals**-1*factor

    plt.figure()
    plt.plot(ens, mxeffs, label='$m_{z, eff}$')
    plt.plot(ens, myeffs, label='$m_{y, eff}$')
    plt.vlines(efermi*1000, 0, np.max(mxeffs), 'r', label = '$E_f(10\\cdot n)$')
    plt.vlines(efermi_no10times*1000, 0, np.max(mxeffs), 'k', '--', label = '$E_f(n)$')
    plt.xlabel('$E$ (meV)')
    plt.ylabel('$m_{eff}$ $(m_e)$')
    plt.legend()

    plt.figure()
    plt.plot(ens,np.array(mxeffs)/np.array(myeffs))
    plt.vlines(efermi*1000, 0, np.max(np.array(mxeffs)/np.array(myeffs)), 'r')
    plt.vlines(efermi_no10times*1000, 0, np.max(np.array(mxeffs)/np.array(myeffs)), 'k', '--')
    plt.xlabel('$E$ (meV)')
    plt.ylabel('$m_z/m_y$')
    plt.vlines(efermi*1000, 0, np.max(mxeffs), 'r', label = '$E_f(10\\cdot n)$')
    plt.vlines(efermi_no10times*1000, 0, np.max(mxeffs), 'k', '--', label = '$E_f(n)$')
    plt.legend()


    band2 = band2/eV2Hartree*1000

    enval=efermi*1000
    enval_no10times = efermi_no10times*1000
    cmap = mpl.cm.viridis
    


    nlev = 15


    ##########################################
    norm = mpl.colors.Normalize(vmin=np.min(band2), vmax= np.max(band2))
    fig = plt.figure()
    fig.set_size_inches(12, 6)
    ax = plt.subplot(111)
    drawing = ax.imshow(band2,interpolation='bilinear',cmap=cmap,norm=norm,
                        extent=(np.min(qxvals),np.max(qxvals),np.min(qyvals),np.max(qyvals)))#ax.tricontourf(qxvals,qyvals,band, levels = nlev)
    colorbar = fig.colorbar(drawing, ax = ax)
    ax.contour(X,Y,np.transpose(band2), colors='k', levels = nlev)
    ax.scatter(qxvals,qyvals,s=3,c='k', label = 'Grid')
    ax.scatter(0,0,s=10,c='r',marker='x', label = '$\\Gamma$')
    plt.legend()

    fermilevel = ax.contour(X,Y,np.transpose(band2), [enval_no10times], colors='b', zorder=10)
    fermilevel = ax.contour(X,Y,np.transpose(band2), [enval], colors='r', zorder=10)
    
    ax.set_xlabel('$q_z$')
    ax.set_ylabel('$q_y$')
    ax.set_title('$E(q_z,q_y)$')

    ##########################################
    nlev = 30

    norm = mpl.colors.Normalize(vmin=np.min(gridddx), vmax= np.max(gridddx))
    fig = plt.figure()
    fig.set_size_inches(12, 6)
    ax = plt.subplot(111)
    drawing2 = ax.imshow(gridddx,interpolation='bilinear',cmap=cmap,norm=norm,
                        extent=(np.min(qxvals),np.max(qxvals),np.min(qyvals),np.max(qyvals)))#ax.tricontourf(qxvals,qyvals,band, levels = nlev)
    colorbar = fig.colorbar(drawing2, ax = ax)
    ax.contour(X,Y,gridddx, colors='k', levels = nlev)
    ax.scatter(qxvals,qyvals,s=3,c='k', label = 'Grid')
    ax.scatter(0,0,s=10,c='r',marker='x', label = '$\\Gamma$')
    plt.legend()

    fermilevel = ax.contour(X,Y,np.transpose(band2), [enval_no10times], colors='b', zorder=10)
    fermilevel = ax.contour(X,Y,np.transpose(band2), [enval], colors='r', zorder=10)
    ax.set_xlabel('$q_z$')
    ax.set_ylabel('$q_y$')
    ax.set_title('$\\frac{\\partial^2 E}{\\partial q_z ^2}$')
    ##########################################

    norm = mpl.colors.Normalize(vmin=np.min(gridddy), vmax= np.max(gridddy))
    fig = plt.figure()
    fig.set_size_inches(12, 6)
    ax = plt.subplot(111)
    drawing2 = ax.imshow(gridddy,interpolation='bilinear',cmap=cmap,norm=norm,
                        extent=(np.min(qxvals),np.max(qxvals),np.min(qyvals),np.max(qyvals)))#ax.tricontourf(qxvals,qyvals,band, levels = nlev)
    colorbar = fig.colorbar(drawing2, ax = ax)
    ax.contour(X,Y,gridddy, colors='k', levels = nlev)
    ax.scatter(qxvals,qyvals,s=3,c='k', label = 'Grid')
    ax.scatter(0,0,s=10,c='r',marker='x', label = '$\\Gamma$')
    plt.legend()

    fermilevel = ax.contour(X,Y,np.transpose(band2), [enval_no10times], colors='b', zorder=10)
    fermilevel = ax.contour(X,Y,np.transpose(band2), [enval], colors='r', zorder=10)
    ax.set_xlabel('$q_z$')
    ax.set_ylabel('$q_y$')
    ax.set_title('$\\frac{\\partial^2 E}{\\partial q_y ^2}$')
    ##########################################

    nlev = 10
    norm = mpl.colors.Normalize(vmin=np.min(0), vmax= np.max(10))
    levelsm = np.linspace(0,1, nlev)

    fig = plt.figure()
    fig.set_size_inches(12, 6)
    ax = plt.subplot(111)
    drawing2 = ax.imshow(np.power(gridddx,-1),interpolation='bilinear',cmap=cmap,norm=norm,
                        extent=(np.min(qxvals),np.max(qxvals),np.min(qyvals),np.max(qyvals)))#ax.tricontourf(qxvals,qyvals,band, levels = nlev)
    colorbar = fig.colorbar(drawing2, ax = ax)
    ax.contour(X,Y,np.power(gridddx,-1), colors='k', levels = levelsm)
    ax.scatter(qxvals,qyvals,s=3,c='k', label = 'Grid')
    ax.scatter(0,0,s=10,c='r',marker='x', label = '$\\Gamma$')
    plt.legend()

    fermilevel = ax.contour(X,Y,np.transpose(band2), [enval_no10times], colors='b', zorder=10)
    fermilevel = ax.contour(X,Y,np.transpose(band2), [enval], colors='r', zorder=10)
    ax.set_xlabel('$q_z$')
    ax.set_ylabel('$q_y$')
    ax.set_title('$m_z$')
    ##########################################

    norm = mpl.colors.Normalize(vmin=np.min(0), vmax= np.max(5))
    levelsm = np.linspace(0,1, nlev)

    fig = plt.figure()
    fig.set_size_inches(12, 6)
    ax = plt.subplot(111)
    drawing2 = ax.imshow(np.power(gridddy,-1),interpolation='bilinear',cmap=cmap,norm=norm,
                        extent=(np.min(qxvals),np.max(qxvals),np.min(qyvals),np.max(qyvals)))#ax.tricontourf(qxvals,qyvals,band, levels = nlev)
    colorbar = fig.colorbar(drawing2, ax = ax)
    ax.contour(X,Y,np.power(gridddy,-1), colors='k', levels = levelsm)
    ax.scatter(qxvals,qyvals,s=3,c='k', label = 'Grid')
    ax.scatter(0,0,s=10,c='r',marker='x', label = '$\\Gamma$')
    plt.legend()

    fermilevel = ax.contour(X,Y,np.transpose(band2), [enval_no10times], colors='b', zorder=10)
    fermilevel = ax.contour(X,Y,np.transpose(band2), [enval], colors='r', zorder=10)
    ax.set_xlabel('$q_z$')
    ax.set_ylabel('$q_y$')
    ax.set_title('$m_y$')
    plt.show()