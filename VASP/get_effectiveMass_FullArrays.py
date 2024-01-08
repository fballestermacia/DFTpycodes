import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from readKgrid import fromOUTCARtoplot
from scipy.integrate import simps
from numpy.polynomial.polynomial import polyder, polyval, polyfit
from scipy.interpolate import griddata
import scipy.constants as cte

def directionalPolyFitFULLARRAY(band,xs, ys, polyorder):
    
    
    polyxs = polyfit(xs,band, polyorder)
    polyys = polyfit(ys,np.transpose(band), polyorder)
    
    return polyxs, polyys


def gausiana(x,center,sigma):
    a0 = 1/(sigma*np.sqrt(2*np.pi))
    return a0*np.exp(-(np.power(x-center,2))/(2*sigma**2))

if __name__ == '__main__':
    #CONSTANTS
    eV2Hartree = 1/27.2113845#0.0367493
    polyorder = 10
    sigx = 0.0003
    sigy = 1*sigx
    alat=15.25
    #########
    

    qz = np.array([-4.9231649190557025-8.0899999999999999, 0. ,7.4862572210537408])
    #qz = np.array([-0.1495968, 0., 0.1495966])
    qz = qz#/np.linalg.norm(qz)

    qy = np.array([0,4.4800000000000004,0])

    qx = np.cross(qy, qz)
    

    '''qz = np.array([-0.1495968000,0,0.1495966000])
    qz = qz
    qy = np.array([0,1.,0])
    qx = np.array([0.1495966000,0,0.1495968000])'''

    qz *= 1/np.linalg.norm(qz)
    qx *= 1/np.linalg.norm(qx)
    qy *= 1/np.linalg.norm(qy)

    bands, efermi, qxvals, qyvals, qzvals, occupation = fromOUTCARtoplot(outcarfile='Ag2Te\\PGHR_XZ\\OUTCAR',kpointsfile="Ag2Te\\PGHR_XZ\\KPOINTS", qx=qx, qy=qz, qz=qy, weightfilter = 0)

    dummy = np.transpose(bands[0])
    dummyoc = np.transpose(occupation[0])

    index = np.sum(dummyoc, axis=0)

    band = ((dummy)[int(index[0])])
    band = (band-np.min(band))*eV2Hartree
    

    shape = [51, 71]

    NX = np.linspace(-0.20, 0.20, 61)
    NY = np.linspace(-0.20, 0.20, 61)
    
    QX, QY = np.meshgrid(NX,NY)
    
    KX = np.outer(QX.flatten(),qx)
    KY = np.outer(QY.flatten(),qy)

    xs = NX
    ys = NY

    band2 = band.reshape(shape[::-1])
    polyxs, polyys = directionalPolyFitFULLARRAY(band2,xs, ys, polyorder)

    derx = polyder(polyxs,2, axis=0)
    dery = polyder(polyys,2, axis=0)

    derxvals = polyval(xs, derx, tensor=True)
    deryvals = polyval(ys, dery, tensor=True)
    
    efermi_no10times = 42/1000 # meVs, fermi energi with previous doping
    efermi = 196/1000 # meVs, fermi energi with 10 times the doping
    envals = np.linspace(np.min(band2), np.min(band2)+0.25*eV2Hartree, 200)

    
    xx = np.linspace(np.min(xs), np.max(xs), len(xs))
    yy = np.linspace(np.min(ys), np.max(ys), len(ys))

    X,Y = np.meshgrid(xx,yy)

    '''gridddx = griddata(np.transpose([qxvals,qyvals]),np.transpose(derxvals).flatten(),(X,Y), method='linear')
    gridddy = griddata(np.transpose([qxvals,qyvals]),deryvals.flatten(),(X,Y), method='linear')
    gridE = griddata(np.transpose([qxvals,qyvals]),band2.flatten(),(X,Y), method='linear')'''

    gridddx = np.transpose(derxvals)
    gridddy = deryvals
    gridE = band2


    for en in envals:
        #Cicle through all data in the grid
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
    plt.plot(ens, mxeffs, label='$m_{x, eff}$')
    plt.plot(ens, myeffs, label='$m_{z, eff}$')
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
    plt.ylabel('$m_x/m_z$')
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
    fig1 = plt.figure()
    fig1.set_size_inches(12, 6)
    ax1 = plt.subplot(111)
    drawing = ax1.imshow(band2,interpolation='bilinear',cmap=cmap,norm=norm,
                        extent=(np.min(qxvals),np.max(qxvals),np.min(qyvals),np.max(qyvals)))#ax.tricontourf(qxvals,qyvals,band, levels = nlev)
    colorbar = fig1.colorbar(drawing, ax = ax1)
    ax1.contour(X,Y,np.transpose(band2), colors='k', levels = nlev)
    ax1.scatter(qxvals,qyvals,s=3,c='k', label = 'Grid')
    ax1.scatter(0,0,s=10,c='r',marker='x', label = '$\\Gamma$')
    plt.legend()

    fermilevel = ax1.contour(X,Y,np.transpose(band2), [enval_no10times], colors='b', zorder=10)
    fermilevel = ax1.contour(X,Y,np.transpose(band2), [enval], colors='r', zorder=10)
    
    ax1.set_xlabel('$q_x$')
    ax1.set_ylabel('$q_z$')
    ax1.set_title('$E(q_x,q_z)$')

    ##########################################
    nlev = 30

    norm = mpl.colors.Normalize(vmin=np.min(gridddx), vmax= np.max(gridddx))
    fig2 = plt.figure()
    fig2.set_size_inches(12, 6)
    ax2 = plt.subplot(111)
    drawing2 = ax2.imshow(gridddx,interpolation='bilinear',cmap=cmap,norm=norm,
                        extent=(np.min(qxvals),np.max(qxvals),np.min(qyvals),np.max(qyvals)))#ax.tricontourf(qxvals,qyvals,band, levels = nlev)
    colorbar = fig2.colorbar(drawing2, ax = ax2)
    ax2.contour(X,Y,gridddx, colors='k', levels = nlev)
    ax2.scatter(qxvals,qyvals,s=3,c='k', label = 'Grid')
    ax2.scatter(0,0,s=10,c='r',marker='x', label = '$\\Gamma$')
    plt.legend()

    fermilevel = ax2.contour(X,Y,np.transpose(band2), [enval_no10times], colors='b', zorder=10)
    fermilevel = ax2.contour(X,Y,np.transpose(band2), [enval], colors='r', zorder=10)
    ax2.set_xlabel('$q_x$')
    ax2.set_ylabel('$q_z$')
    ax2.set_title('$\\frac{\\partial^2 E}{\\partial q_x ^2}$')
    
    ##########################################
    #   BE EXTREMELY CAREFUL THERE'S A BUG THAT CHANGES THE NEXT PLOT
    #   WITH THE PREVIOUS ONE AND I HAVE NO IDEA WHY
    ##########################################

    norm = mpl.colors.Normalize(vmin=np.min(gridddy), vmax= np.max(gridddy))
    fig3 = plt.figure()
    fig3.set_size_inches(12, 6)
    ax3 = plt.subplot(111)
    drawing2 = ax3.imshow(gridddy,interpolation='bilinear',cmap=cmap,norm=norm,
                        extent=(np.min(qxvals),np.max(qxvals),np.min(qyvals),np.max(qyvals)))#ax.tricontourf(qxvals,qyvals,band, levels = nlev)
    colorbar = fig3.colorbar(drawing2, ax = ax3)
    ax3.contour(X,Y,gridddy, colors='k', levels = nlev)
    ax3.scatter(qxvals,qyvals,s=3,c='k', label = 'Grid')
    ax3.scatter(0,0,s=10,c='r',marker='x', label = '$\\Gamma$')
    plt.legend()

    fermilevel = ax3.contour(X,Y,np.transpose(band2), [enval_no10times], colors='b', zorder=10)
    fermilevel = ax3.contour(X,Y,np.transpose(band2), [enval], colors='r', zorder=10)
    ax3.set_xlabel('$q_x$')
    ax3.set_ylabel('$q_z$')
    ax3.set_title('$\\frac{\\partial^2 E}{\\partial q_z ^2}$')
    
    ##########################################
    

    nlev = 10
    norm = mpl.colors.Normalize(vmin=np.min(0), vmax= np.max(2))
    levelsm = np.linspace(0,2, nlev)

    fig4 = plt.figure()
    fig4.set_size_inches(12, 6)
    ax4 = plt.subplot(111)
    drawing2 = ax4.imshow(np.power(gridddx,-1),interpolation='bilinear',cmap=cmap,norm=norm,
                        extent=(np.min(qxvals),np.max(qxvals),np.min(qyvals),np.max(qyvals)))#ax.tricontourf(qxvals,qyvals,band, levels = nlev)
    colorbar = fig4.colorbar(drawing2, ax = ax4)
    ax4.contour(X,Y,np.power(gridddx,-1), colors='k', levels = levelsm)
    ax4.scatter(qxvals,qyvals,s=3,c='k', label = 'Grid')
    ax4.scatter(0,0,s=10,c='r',marker='x', label = '$\\Gamma$')
    plt.legend()

    fermilevel = ax4.contour(X,Y,np.transpose(band2), [enval_no10times], colors='b', zorder=10)
    fermilevel = ax4.contour(X,Y,np.transpose(band2), [enval], colors='r', zorder=10)
    ax4.set_xlabel('$q_x$')
    ax4.set_ylabel('$q_z$')
    ax4.set_title('$m_x$')

    ##########################################
    #   BE EXTREMELY CAREFUL THERE'S A BUG THAT CHANGES THE NEXT PLOT
    #   WITH THE PREVIOUS ONE AND I HAVE NO IDEA WHY
    ##########################################

    norm = mpl.colors.Normalize(vmin=np.min(0), vmax= np.max(1))
    levelsm = np.linspace(0,1, nlev)

    fig5 = plt.figure()
    fig5.set_size_inches(12, 6)
    ax5 = plt.subplot(111)
    drawing2 = ax5.imshow(np.power(gridddy,-1),interpolation='bilinear',cmap=cmap,norm=norm,
                        extent=(np.min(qxvals),np.max(qxvals),np.min(qyvals),np.max(qyvals)))#ax.tricontourf(qxvals,qyvals,band, levels = nlev)
    colorbar = fig5.colorbar(drawing2, ax = ax5)
    ax5.contour(X,Y,np.power(gridddy,-1), colors='k', levels = levelsm)
    ax5.scatter(qxvals,qyvals,s=3,c='k', label = 'Grid')
    ax5.scatter(0,0,s=10,c='r',marker='x', label = '$\\Gamma$')
    plt.legend()

    fermilevel = ax5.contour(X,Y,np.transpose(band2), [enval_no10times], colors='b', zorder=10)
    fermilevel = ax5.contour(X,Y,np.transpose(band2), [enval], colors='r', zorder=10)
    ax5.set_xlabel('$q_x$')
    ax5.set_ylabel('$q_z$')
    ax5.set_title('$m_z$')
    plt.show()