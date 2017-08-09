#plot the Wetzel data for Cooper proposal and summer research 2014

import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

def plotA():

    inputdata = Table.read('WetzelData_formatted.dat', format = 'ascii')

    mstar = np.array(inputdata['\xef\xbb\xbfmstar'])
    mhalo1213lo = np.array(inputdata['M_halo_12_13_lo'])
    mhalo1213hi = np.array(inputdata['M_halo_12_13_hi'])
    mhalo1314lo = np.array(inputdata['M_halo_13_14_lo'])
    mhalo1314hi = np.array(inputdata['M_halo_13_14_hi'])
    mhalo1415lo = np.array(inputdata['M_halo_14_15_lo'])
    mhalo1415hi = np.array(inputdata['M_halo_14_15_hi'])

    cut = mhalo1213lo != -1
    clean_mhalo1213lo = mhalo1213lo[cut]
    clean_mhalo1213hi = mhalo1213hi[cut]
    clean_mstar = mstar[cut]

    mass = np.array([7.0, 9.0])
    masserror = np.array([0.5, 1.0])
    tq = np.array([1.5, 8.0])
    tqerror = np.array([0.5, 1.0])
    linex = np.array([9.0, 11.2])
    liney = np.array([8.0, 1.5])
    
    axwidth = 3
    axlength = 10
    fontsize=28

    plt.rc('axes',linewidth=axwidth)  #make sure to do this before making any figures

    fig = plt.figure(figsize=(13,8))
    plt.subplot2grid((1,1),(0,0))
    plt.subplots_adjust(left=0.12, bottom=0.12, right=0.97, top=0.97, wspace=0.1, hspace=0.1)
    ax = plt.gca()
    minorLocatorx = AutoMinorLocator(2)
    minorLocatory = AutoMinorLocator(2)
    ax.xaxis.set_minor_locator(minorLocatorx)
    ax.yaxis.set_minor_locator(minorLocatory)
    #<do plotting>

    xerror = np.empty([2,1])
    xerror[0,0]=0.27
    xerror[1,0]=0.5
    yerror = np.empty([2,1])
    yerror[0,0]=1.4
    yerror[1,0]=1.0
    #plt.figure(figsize = (10,5), dpi = 100)
    p1 = plt.errorbar(6.5, 1.4, xerr = None, yerr=None, fmt='ko', markersize = 15)
    plt.errorbar(6.5, 1.4, xerr = xerror, yerr=yerror, fmt='ko', markersize = 15, elinewidth = 2.5, capsize = 10.0)
    plt.plot(6.5,1.4, 'ko', markersize = 20)
    p2 = plt.errorbar(9.0, 8.0, xerr = None, yerr=None, fmt='ks', markersize = 15)
    plt.plot(9.0,8.0, 'ks', markersize = 25)
    
    #p2 = plt.plot(9.0,8.0, color='k', marker = 'd', markersize = 20)
    #plt.plot(mass, tq, 'k--')
    #plt.plot(linex, liney, 'k--')
    plt.fill_between(mstar, mhalo1415hi, mhalo1415lo, color="c", edgecolor = "c")
    plt.fill_between(mstar, mhalo1314hi, mhalo1314lo, color="Gold", edgecolor = "Gold")
    plt.fill_between(clean_mstar, clean_mhalo1213hi, clean_mhalo1213lo, color="m", edgecolor = "m")

    pp1213 = plt.Rectangle((0, 0), 1, 1, fc="m", ec = "m")
    pp1314 = plt.Rectangle((0, 0), 1, 1, fc="Gold", ec = "Gold")
    pp1415 = plt.Rectangle((0, 0), 1, 1, fc="c", ec = "c")
    pp0 = plt.Rectangle((0, 0), 1, 1, fc="w", ec="w")
    

    plt.legend([pp0,pp1213,pp1314,pp1415, p2, pp0, p1, pp0], [r"$\rm Wetzel \, et \, al. \, (2013)$", r"$\rm log(M_{vir}/M_{\odot}) = [12,13]$", r"$\rm log(M_{vir}/M_{\odot}) = [13,14]$", r"$\rm log(M_{vir}/M_{\odot}) = [14,15]$", r"$\rm Wheeler \, et \, al. \, (2014)$", r"$\rm log(M_{vir}/M_{\odot}) = 13.5$", r"$\rm Fillingham \, et \, al. \, (2015)$", r"$\rm log(M_{vir}/M_{\odot}) = 12.5$"], loc = (0.05,0.58), frameon=False, numpoints=1, prop={'size':20})
    
    #plt.legend(loc = (0.05,0.58), frameon=False, numpoints=1, prop={'size':13})
    
    plt.xlabel(r'$\rm Satellite \, Stellar \, Mass \, (M_{\odot})$', fontsize = 28)
    plt.ylabel(r'$\rm Quenching \ Timescale \ (Gyr)$', fontsize = 28)
    plt.axis([5.68,11.5,-0.5,10])
    
    xtickloc = [6,7,8,9,10,11]
    xtickstr = [r'$\rm 10^{6}$',r'$\rm 10^{7}$',r'$\rm 10^{8}$',r'$\rm 10^{9}$',r'$\rm 10^{10}$',r'$\rm 10^{11}$']
    #ticklabels = np.array([r'$\rm 10^{7}$',r'$\rm 10^{8}$',r'$\rm 10^{9}$',r'$\rm 10^{10}$',r'$\rm 10^{11}$'])
    #plt.xticks([7,8,9,10,11], ticklabels)
    ytickloc = [0,2,4,6,8,10]
    ytickstr = ['$'+str(kk)+'$' for kk in ytickloc]

    ax.set_xticks(xtickloc)
    ax.set_xticklabels(xtickstr, position = (0,-0.01))
    ax.set_yticks(ytickloc)
    ax.set_yticklabels(ytickstr)

    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for line in ax.get_xticklines() + ax.get_yticklines():
        line.set_markersize(18)
        line.set_markeredgewidth(3)
    for tick in ax.xaxis.get_minor_ticks():
        tick.label1.set_fontsize(fontsize/2)
    for tick in ax.yaxis.get_minor_ticks():
        tick.label1.set_fontsize(fontsize/2)

    ax.tick_params(which='major',width=axwidth,length=axlength+5)
    ax.tick_params(which='minor',width=axwidth,length=axlength)
    
    #plt.savefig('Proposal_plot.pdf')
    plt.show()

#############################################################################
#Include 2 data points
#############################################################################

def plotB():
    
    inputdata = Table.read('WetzelData_formatted.dat', format = 'ascii')
    
    mstar = np.array(inputdata['\xef\xbb\xbfmstar'])
    mhalo1213lo = np.array(inputdata['M_halo_12_13_lo'])
    mhalo1213hi = np.array(inputdata['M_halo_12_13_hi'])
    mhalo1314lo = np.array(inputdata['M_halo_13_14_lo'])
    mhalo1314hi = np.array(inputdata['M_halo_13_14_hi'])
    mhalo1415lo = np.array(inputdata['M_halo_14_15_lo'])
    mhalo1415hi = np.array(inputdata['M_halo_14_15_hi'])
    
    cut = mhalo1213lo != -1
    clean_mhalo1213lo = mhalo1213lo[cut]
    clean_mhalo1213hi = mhalo1213hi[cut]
    clean_mstar = mstar[cut]
    
    mass = np.array([7.0, 9.0])
    masserror = np.array([0.5, 1.0])
    tq = np.array([1.5, 8.0])
    tqerror = np.array([0.5, 1.0])
    linex = np.array([9.0, 11.2])
    liney = np.array([8.0, 1.5])
    
    axwidth = 3
    axlength = 10
    fontsize=28
    
    plt.rc('axes',linewidth=axwidth)  #make sure to do this before making any figures
    
    fig = plt.figure(figsize=(8,8))
    plt.subplot2grid((1,1),(0,0))
    plt.subplots_adjust(left=0.12, bottom=0.13, right=0.97, top=0.97, wspace=0.1, hspace=0.1)
    ax = plt.gca()
    minorLocatorx   = AutoMinorLocator(2)
    minorLocatory   = AutoMinorLocator(2)
    ax.xaxis.set_minor_locator(minorLocatorx)
    ax.yaxis.set_minor_locator(minorLocatory)
    #<do plotting>
    
    xerror1 = np.empty([2,1])
    xerror1[0,0]=0.17
    xerror1[1,0]=0.25
    xerror2 = np.empty([2,1])
    xerror2[0,0]=0.19
    xerror2[1,0]=0.23
    yerror1 = np.empty([2,1])
    yerror1[0,0]=1.2
    yerror1[1,0]=0.73
    yerror2 = np.empty([2,1])
    yerror2[0,0]=1.4
    yerror2[1,0]=1.8
    yerrorc = np.empty([2,1])
    yerrorc[0,0]=1.76
    yerrorc[1,0]=0.87
    #plt.figure(figsize = (10,5), dpi = 100)
    
    #### 100% low, 80% high
    #p1a = plt.errorbar(6.3, 0.1, xerr = None, yerr=None, fmt='ko', markersize = 15)
    #plt.errorbar(6.3, 0.1, xerr = xerror1, yerr=None, fmt='ko', markersize = 15, elinewidth = 2.5, capsize = 10.0)
    #plt.plot(6.3,0.1, 'ko', markersize = 20)
    #p1b = plt.errorbar(7.4, 3.2, xerr = None, yerr=None, fmt='ko', markersize = 15)
    #plt.errorbar(7.4, 3.2, xerr = xerror2, yerr=None, fmt='ko', markersize = 15, elinewidth = 2.5, capsize = 10.0)
    #plt.plot(7.4, 3.2, 'ko', markersize = 20)
    
    ### 90% quenched in both
    p1a = plt.errorbar(6.3, 1.29, xerr = None, yerr=None, fmt='ko', markersize = 15)
    plt.errorbar(6.3, 1.29, xerr = xerror1, yerr=yerror1, fmt='ko', markersize = 15, elinewidth = 2.5, capsize = 10.0)
    plt.plot(6.3,1.29, 'ko', markersize = 20)
    p1b = plt.errorbar(7.4, 1.44, xerr = None, yerr=None, fmt='ko', markersize = 15)
    plt.errorbar(7.4, 1.44, xerr = xerror2, yerr=yerror2, fmt='ko', markersize = 15, elinewidth = 2.5, capsize = 10.0)
    plt.plot(7.4, 1.44, 'ko', markersize = 20)
    
    #### Corals point
    p2 = plt.errorbar(9.0, 7.76, xerr = None, yerr=None, fmt='ks', markersize = 15)
    plt.errorbar(9.0, 7.76, xerr = 0.3, yerr=yerrorc, fmt='ks', markersize = 15, elinewidth = 2.5, capsize = 10.0)
    plt.plot(9.0, 7.76, 'ks', markersize = 25)
    
    #p2 = plt.plot(9.0,8.0, color='k', marker = 'd', markersize = 20)
    #plt.plot(mass, tq, 'k--')
    #plt.plot(linex, liney, 'k--')
    plt.fill_between(mstar, mhalo1415hi, mhalo1415lo, color="c", edgecolor = "c")
    plt.fill_between(mstar, mhalo1314hi, mhalo1314lo, color="Gold", edgecolor = "Gold")
    plt.fill_between(clean_mstar, clean_mhalo1213hi, clean_mhalo1213lo, color="m", edgecolor = "m")
    
    pp1213 = plt.Rectangle((0, 0), 1, 1, fc="m", ec = "m")
    pp1314 = plt.Rectangle((0, 0), 1, 1, fc="Gold", ec = "Gold")
    pp1415 = plt.Rectangle((0, 0), 1, 1, fc="c", ec = "c")
    pp0 = plt.Rectangle((0, 0), 1, 1, fc="w", ec="w")
    
    
    #plt.legend([pp0,pp1213,pp1314,pp1415, p2, pp0, p1a, pp0], [r"$\rm Wetzel \, et \, al. \, (2013)$", r"$\rm log(M_{vir}/M_{\odot}) = [12,13]$", r"$\rm log(M_{vir}/M_{\odot}) = [13,14]$", r"$\rm log(M_{vir}/M_{\odot}) = [14,15]$", r"$\rm Wheeler \, et \, al. \, (2014)$", r"$\rm log(M_{vir}/M_{\odot}) = 13.5$", r"$\rm Fillingham \, et \, al. \, (2015)$", r"$\rm log(M_{vir}/M_{\odot}) = 12.5$"], loc = (0.05,0.4), frameon=False, numpoints=1, prop={'size':21})
    
    
    
    #first_legend = plt.legend([p2, pp0, p1a, pp0], [r"$\rm Wheeler \, et \, al. \, (2014)$", r"$\rm log(M_{vir}/M_{\odot}) = 13.5$", r"$\rm Fillingham \, et \, al. \, (2015)$", r"$\rm log(M_{vir}/M_{\odot}) = 12.5$"], loc = (0.02,0.68), frameon=False, numpoints=1, prop={'size':20}, handletextpad = 0.1)
    
    #plt.gca().add_artist(first_legend)
    
    #plt.legend([pp0,pp1213,pp1314,pp1415], [r"$\rm Wetzel \, et \, al. \, (2013)$", r"$\rm log(M_{vir}/M_{\odot}) = [12,13]$", r"$\rm log(M_{vir}/M_{\odot}) = [13,14]$", r"$\rm log(M_{vir}/M_{\odot}) = [14,15]$"], loc = (0.64,0.69), frameon=False, numpoints=1, prop={'size':20})
    plt.legend([pp0,pp1213,pp1314,pp1415], [r"$\rm Wetzel \, et \, al. \, (2013)$", r"$\rm log(M_{vir}/M_{\odot}) = [12,13]$", r"$\rm log(M_{vir}/M_{\odot}) = [13,14]$", r"$\rm log(M_{vir}/M_{\odot}) = [14,15]$"], loc = (0.48,0.71), frameon=False, numpoints=1, prop={'size':18})
    
    
    
    #plt.plot(6.05, 11, 'ks', markersize = 20)
    #lt.plot(6.05, 9, 'ko', markersize = 20)
    
    plt.xlabel(r'$\rm Stellar \, Mass \, (M_{\odot})$', fontsize = 28)
    plt.ylabel(r'$\rm Quenching \ Timescale \ (Gyr)$', fontsize = 28)
    #plt.axis([5.68,12,-0.5,12])
    plt.axis([7.8,12,-0.5,12])
    
    #xtickloc = [6,7,8,9,10,11]
    #xtickstr = [r'$\rm 10^{6}$',r'$\rm 10^{7}$',r'$\rm 10^{8}$',r'$\rm 10^{9}$',r'$\rm 10^{10}$',r'$\rm 10^{11}$']
    xtickloc = [8,9,10,11]
    xtickstr = [r'$\rm 10^{8}$',r'$\rm 10^{9}$',r'$\rm 10^{10}$',r'$\rm 10^{11}$']
    #ticklabels = np.array([r'$\rm 10^{7}$',r'$\rm 10^{8}$',r'$\rm 10^{9}$',r'$\rm 10^{10}$',r'$\rm 10^{11}$'])
    #plt.xticks([7,8,9,10,11], ticklabels)
    ytickloc = [0,2,4,6,8,10]
    ytickstr = ['$'+str(kk)+'$' for kk in ytickloc]
    
    ax.set_xticks(xtickloc)
    ax.set_xticklabels(xtickstr, position = (0,-0.01))
    ax.set_yticks(ytickloc)
    ax.set_yticklabels(ytickstr)
    
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for line in ax.get_xticklines() + ax.get_yticklines():
        line.set_markersize(18)
        line.set_markeredgewidth(3)
    for tick in ax.xaxis.get_minor_ticks():
        tick.label1.set_fontsize(fontsize/2)
    for tick in ax.yaxis.get_minor_ticks():
        tick.label1.set_fontsize(fontsize/2)
    
    ax.tick_params(which='major',width=axwidth,length=axlength+5)
    ax.tick_params(which='minor',width=axwidth,length=axlength)

    plt.savefig('plot_for_Tim_advancement.png')
    plt.savefig('plot_for_Tim_advancement.pdf')
    plt.show()


