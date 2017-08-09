import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline as IUSpline

def test():

    xdata = np.loadtxt('ELVIS_Data/Etracks/AllTrees/Kek&Kauket/X.txt')
    ydata = np.loadtxt('ELVIS_Data/Etracks/AllTrees/Kek&Kauket/Y.txt')
    zdata = np.loadtxt('ELVIS_Data/Etracks/AllTrees/Kek&Kauket/Z.txt')
    adata = np.loadtxt('ELVIS_Data/Etracks/AllTrees/Kek&Kauket/scale.txt')
    virdata = np.loadtxt('ELVIS_Data/Etracks/AllTrees/Kek&Kauket/Rvir.txt')

    x = xdata[25]
    y = ydata[25]
    z = zdata[25]
    a = adata[25]
    
    subcut = a > 0.
    x_clean = x[subcut]
    y_clean = y[subcut]
    z_clean = z[subcut]
    a_clean = a[subcut]
    x_cleansort = x_clean[::-1]
    y_cleansort = y_clean[::-1]
    z_cleansort = z_clean[::-1]
    a_cleansort = a_clean[::-1]
    
    interval = 750./len(a)
    afinal = np.linspace(min(a_clean),1.,len(a_clean)*interval)
    
    xinterp = IUSpline(a_cleansort, x_cleansort, k=3)
    yinterp = IUSpline(a_cleansort, y_cleansort, k=3)
    zinterp = IUSpline(a_cleansort, z_cleansort, k=3)
        
    xpos = xinterp(afinal)
    ypos = yinterp(afinal)
    zpos = zinterp(afinal)
    
    Rvir1 = virdata[0]
    vir1cut = Rvir1 > 0
    Rvir1_clean = Rvir1[vir1cut]
    Rvir1_cleansort = Rvir1_clean[::-1]
    a1 = adata[0]
    a1_clean = a1[vir1cut]
    a1_cleansort = a1_clean[::-1]
    a1final = np.linspace(min(a1_clean),1.,len(a1_clean)*interval)
            
    Rvir2 = virdata[1]
    vir2cut = Rvir2 > 0
    Rvir2_clean = Rvir2[vir2cut]
    Rvir2_cleansort = Rvir2_clean[::-1]
    a2 = adata[1]
    a2_clean = a2[vir2cut]
    a2_cleansort = a2_clean[::-1]
    a2final = np.linspace(min(a2_clean),1.,len(a2_clean)*interval)
            
    Rvir1_interp = IUSpline(a1_cleansort, Rvir1_cleansort, k=3)
    Rvir2_interp = IUSpline(a2_cleansort, Rvir2_cleansort, k=3)
            
    virial1 = Rvir1_interp(a1final)
    virial2 = Rvir2_interp(a2final)
            
    print 'Host 1 x position'
    xHost1 = xdata[0]
    xHost1_clean = xHost1[subcut]
    xHost1_cleansort = xHost1_clean[::-1]
    xHost1interp = IUSpline(a_cleansort, xHost1_cleansort, k=3)
    xHost1pos = xHost1interp(afinal)
            
    print 'Host 2 x position'
    xHost2 = xdata[1]
    xHost2_clean = xHost2[subcut]
    xHost2_cleansort = xHost2_clean[::-1]
    xHost2interp = IUSpline(a_cleansort, xHost2_cleansort, k=3)
    xHost2pos = xHost2interp(afinal)
            
    print 'Host 1 y position'
    yHost1 = ydata[0]
    yHost1_clean = yHost1[subcut]
    yHost1_cleansort = yHost1_clean[::-1]
    yHost1interp = IUSpline(a_cleansort, yHost1_cleansort, k=3)
    yHost1pos = yHost1interp(afinal)
            
    print 'Host 2 y position'
    yHost2 = ydata[1]
    yHost2_clean = yHost2[subcut]
    yHost2_cleansort = yHost2_clean[::-1]
    yHost2interp = IUSpline(a_cleansort, yHost2_cleansort, k=3)
    yHost2pos = yHost2interp(afinal)
            
    print 'Host 1 z position'
    zHost1 = zdata[0]
    zHost1_clean = zHost1[subcut]
    zHost1_cleansort = zHost1_clean[::-1]
    zHost1interp = IUSpline(a_cleansort, zHost1_cleansort, k=3)
    zHost1pos = zHost1interp(afinal)
            
    print 'Host 2 z position'
    zHost2 = zdata[1]
    zHost2_clean = zHost2[subcut]
    zHost2_cleansort = zHost2_clean[::-1]
    zHost2interp = IUSpline(a_cleansort, zHost2_cleansort, k=3)
    zHost2pos = zHost2interp(afinal)
            
            
    dx1 = xHost1pos - xpos
    dy1 = yHost1pos - ypos
    dz1 = zHost1pos - zpos
    dx2 = xHost2pos - xpos
    dy2 = yHost2pos - ypos
    dz2 = zHost2pos - zpos
    
    distx1 = xHost1_clean - x_clean
    disty1 = yHost1_clean - y_clean
    distz1 = zHost1_clean - z_clean
    distx2 = xHost2_clean - x_clean
    disty2 = yHost2_clean - y_clean
    distz2 = zHost2_clean - z_clean
            
    comoving1 = np.sqrt(np.power(dx1,2) + np.power(dy1,2) + np.power(dz1,2))
    comoving2 = np.sqrt(np.power(dx2,2) + np.power(dy2,2) + np.power(dz2,2))
    
    distance1 = np.sqrt(np.power(distx1,2) + np.power(disty1,2) + np.power(distz1,2))
    distance2 = np.sqrt(np.power(distx2,2) + np.power(disty2,2) + np.power(distz2,2))

    plt.figure(1)
    plt.axis([0,1,30,40])

    plt.plot(a, x, 'bo')
    plt.plot(afinal, xpos, 'b-')
    plt.plot(a, y, 'go')
    plt.plot(afinal, ypos, 'g-')
    plt.plot(a, z, 'ro')
    plt.plot(afinal, zpos, 'r-')
    
    plt.plot(a_clean, xHost1_clean, 'co')
    plt.plot(afinal, xHost1pos, 'c-')
    plt.plot(a_clean, yHost1_clean, 'mo')
    plt.plot(afinal, yHost1pos, 'm-')
    plt.plot(a_clean, zHost1_clean, 'yo')
    plt.plot(afinal, zHost1pos, 'y-')

    
    plt.figure(2)
    plt.axis([0,1,0,5])
    
    plt.plot(afinal, comoving1, 'k-')
    plt.plot(a_clean, distance1, 'ko')
    

    plt.show()
