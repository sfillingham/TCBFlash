import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import astropy.cosmology as cosmo
import astropy.units as u


def radial1():

    
    data_input_ver1 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.6virial_interp.dat', format = 'ascii')
    data_input_ver1a = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.6virial_interp_ver1a.dat', format = 'ascii')
    #data_input_ver2 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.6virial_interp_ver2.dat', format = 'ascii')
    #data_input_ver3 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.6virial_interp_ver3.dat', format = 'ascii')
    #data_input_ver4 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.6virial_interp_ver4.dat', format = 'ascii')
    data_input_BAM = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.6virial_interp_BehrooziAM.dat', format = 'ascii')
    
    tq_ver1 = np.array(data_input_ver1['QuenchTime'])*u.Gyr
    tq_ver1a = np.array(data_input_ver1a['QuenchTime'])*u.Gyr
    #tq_ver2 = np.array(data_input_ver2['QuenchTime'])*u.Gyr
    tq_BAM = np.array(data_input_BAM['QuenchTime'])*u.Gyr
    
    z_ver1 = np.array(data_input_ver1['zInfall'])
    z_ver1a = np.array(data_input_ver1a['zInfall'])
    #z_ver2 = np.array(data_input_ver2['zInfall'])
    z_BAM = np.array(data_input_BAM['zInfall'])
    
    t_ver1 = cosmo.lookback_time(z_ver1)
    t_ver1a = cosmo.lookback_time(z_ver1a)
    #t_ver2 = cosmo.lookback_time(z_ver2)
    t_BAM = cosmo.lookback_time(z_BAM)
    
    tau_ver1 = t_ver1 - tq_ver1
    tau_ver1a = t_ver1a - tq_ver1a
    #tau_ver2 = t_ver2 - tq_ver2
    tau_BAM = t_BAM - tq_BAM

    
    dist = np.array([50, 100, 150, 200, 250, 300])
    distfrac = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    quench = np.array([0.44, 0.72, 0.85, 0.94, 0.98, 0.99])
    quench_interp = np.array([0.469, 0.7579, 0.8704, 0.9535, 0.98, 0.99])
    quench_interpfrac_ver1 = np.array([0.6, 0.722, 0.795, 0.867, 0.910, 0.936, 0.965, 0.979, 1.0])
    quench_interpfrac_ver1a = np.array([0.707, 0.894, 0.954, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
    #quench_interpfrac_ver2 = np.array([0.558, 0.693, 0.774, 0.853, 0.901, 0.929, 0.962, 0.977, 1.0])
    quench_interpfrac_BAM = np.array([0.618, 0.754, 0.833, 0.886, 0.924, 0.949, 0.973, 0.986, 1.0])

#####################################
#Plot work
#####################################
    axwidth = 3
    axlength = 10
    fontsize=24

    plt.rc('axes',linewidth=axwidth)  #make sure to do this before making any figures

    #<do plotting>

    plt.figure(2,figsize = (26,8))
    plt.tight_layout()
    plt.subplot(121)
    plt.axis([0, 1.0, 0, 1.1])
    ax = plt.gca()
    p1, = plt.plot(distfrac, quench_interpfrac_ver1, 'k-', linewidth = 4.5)
    #p1a, = plt.plot(distfrac, quench_interpfrac_ver1a, color = 'y', linestyle = '-', linewidth = 3.0)
    #p2, = plt.plot(distfrac, quench_interpfrac_ver2, 'm-', linewidth = 3.0)
    pBAM, = plt.plot(distfrac, quench_interpfrac_BAM, 'c-', linewidth = 3.0)
    
    pp0 = plt.Rectangle((0, 0), 1, 1, fc="k", ec = "k")
    pp1 = plt.Rectangle((0, 0), 1, 1, fc="Gold", ec = "Gold")
    pp2 = plt.Rectangle((0, 0), 1, 1, fc="m", ec = "m")
    pp3 = plt.Rectangle((0, 0), 1, 1, fc="c", ec="c")
    
    plt.legend([pp0,pp2,pp3], [r'$\rm GK14\ AM$', r'$\rm GK14\ AM\ +\ scatter$',  r'$\rm Behroozi\ AM$'], loc = 4, frameon=False, numpoints=1, prop={'size':16})
    
    
    ax.set_xlabel('Distance from Host (Rvir)', fontsize = 24)
    ax.set_ylabel('Quenched Fraction', fontsize = 24)
    
    
    xtickloc = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    xtickstr = ['$'+str(kk)+'$' for kk in xtickloc]
    ytickloc = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    ytickstr = ['$'+str(kk)+'$' for kk in ytickloc]
    
    ax.set_xticks(xtickloc)
    ax.set_xticklabels(xtickstr)
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
    
    plt.subplot(122)
    plt.axis([-1,6,0,1.1])
    ax = plt.gca()
    ax.set_xlabel('Quenching Time Scale (Gyr)', fontsize = 24)
    ax.set_ylabel('Fraction of Dwarfs', fontsize = 24)
    
    n, bins, patches = plt.hist(tau_ver1, bins = range(-1,6,1), histtype = 'step', color = 'black', linewidth = 4.5, normed = True)
    #n, bins, patches = plt.hist(tau_ver1a, bins = range(-1,6,1), histtype = 'step', color = 'Gold', linewidth = 3.0, normed = True)
    #n, bins, patches = plt.hist(tau_ver2, bins = range(-1,6,1), histtype = 'step', color = 'magenta', linewidth = 1.5, normed = True)
    n, bins, patches = plt.hist(tau_BAM, bins = range(-1,6,1), histtype = 'step', color = 'cyan', linewidth = 3.0, normed = True)
    
    
    #plt.text(1.0, 1.0, r'$\rm Garrison-Kimmel\ AM$', color = 'k')
    #plt.text(1.0, 0.90, r'$\rm GK\ AM\ (<\ 0.5\ R_{vir})$', color = 'Gold')
    #plt.text(1.0, 0.95, r'$\rm Behroozi\ AM$', color = 'c')

    pp0 = plt.Rectangle((0, 0), 1, 1, fc="k", ec = "k")
    pp1 = plt.Rectangle((0, 0), 1, 1, fc="Gold", ec = "Gold")
    pp2 = plt.Rectangle((0, 0), 1, 1, fc="m", ec = "m")
    pp3 = plt.Rectangle((0, 0), 1, 1, fc="c", ec="c")
    
    plt.legend([pp0,pp2,pp3], [r'$\rm GK14\ AM$', r'$\rm GK14\ AM\ +\ scatter$',  r'$\rm Behroozi\ AM$'], loc = 1, frameon=False, numpoints=1, prop={'size':16})
    # r'$\rm GK14\ AM\ (<\ 0.5\ R_{vir})$',

    #plt.legend([p1,p2,p3,p4,p5], [r'$\rm Vmax/Vpeak\ >\ 0.3$', r'$\rm Vmax/Vpeak\ >\ 0.4$', r'$\rm Vmax/Vpeak\ >\ 0.5$', r'$\rm Vmax/Vpeak\ >\ 0.6$', r'$\rm Vmax/Vpeak\ >\ 0.7$'], loc = 4, frameon=False, numpoints=1, prop={'size':11})
    
    xtickloc = [0,1,2,3,4,5,6]
    xtickstr = ['$'+str(kk)+'$' for kk in xtickloc]
    ytickloc = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    ytickstr = ['$'+str(kk)+'$' for kk in ytickloc]

    ax.set_xticks(xtickloc)
    ax.set_xticklabels(xtickstr)
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
    
    
    plt.savefig('RadialModel/elvis_paperplot_radialquench_final1.pdf')
    plt.show()




def radial2():

    data_input_ver1_2 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.2virial_interp.dat', format = 'ascii')
    data_input_ver1_3 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.3virial_interp.dat', format = 'ascii')
    data_input_ver1_4 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.4virial_interp.dat', format = 'ascii')
    data_input_ver1_5 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.5virial_interp.dat', format = 'ascii')
    data_input_ver1_6 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.6virial_interp.dat', format = 'ascii')
    data_input_ver1_7 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.7virial_interp.dat', format = 'ascii')
    data_input_ver1_8 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.8virial_interp.dat', format = 'ascii')
    data_input_ver1_9 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.9virial_interp.dat', format = 'ascii')
    data_input_ver1_10 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_1.0virial_interp.dat', format = 'ascii')
    
    #data_input_ver1a_2 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.2virial_interp_ver1a.dat', format = 'ascii')
    #data_input_ver1a_3 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.3virial_interp_ver1a.dat', format = 'ascii')
    #data_input_ver1a_4 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.4virial_interp_ver1a.dat', format = 'ascii')
    #data_input_ver1a_5 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.5virial_interp_ver1a.dat', format = 'ascii')
    #data_input_ver1a_6 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.6virial_interp_ver1a.dat', format = 'ascii')
    #data_input_ver1a_7 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.7virial_interp_ver1a.dat', format = 'ascii')
    #data_input_ver1a_8 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.8virial_interp_ver1a.dat', format = 'ascii')
    #data_input_ver1a_9 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.9virial_interp_ver1a.dat', format = 'ascii')
    #data_input_ver1a_10 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_1.0virial_interp_ver1a.dat', format = 'ascii')
    
    data_input_ver2_2 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.2virial_interp_ver2rand.dat', format = 'ascii')
    data_input_ver2_3 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.3virial_interp_ver2rand.dat', format = 'ascii')
    data_input_ver2_4 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.4virial_interp_ver2rand.dat', format = 'ascii')
    data_input_ver2_5 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.5virial_interp_ver2rand.dat', format = 'ascii')
    data_input_ver2_6 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.6virial_interp_ver2rand.dat', format = 'ascii')
    data_input_ver2_7 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.7virial_interp_ver2rand.dat', format = 'ascii')
    data_input_ver2_8 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.8virial_interp_ver2rand.dat', format = 'ascii')
    data_input_ver2_9 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.9virial_interp_ver2rand.dat', format = 'ascii')
    data_input_ver2_10 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_1.0virial_interp_ver2rand.dat', format = 'ascii')
    
    data_input_BAM_2 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.2virial_interp_BehrooziAM.dat', format = 'ascii')
    data_input_BAM_3 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.3virial_interp_BehrooziAM.dat', format = 'ascii')
    data_input_BAM_4 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.4virial_interp_BehrooziAM.dat', format = 'ascii')
    data_input_BAM_5 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.5virial_interp_BehrooziAM.dat', format = 'ascii')
    data_input_BAM_6 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.6virial_interp_BehrooziAM.dat', format = 'ascii')
    data_input_BAM_7 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.7virial_interp_BehrooziAM.dat', format = 'ascii')
    data_input_BAM_8 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.8virial_interp_BehrooziAM.dat', format = 'ascii')
    data_input_BAM_9 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.9virial_interp_BehrooziAM.dat', format = 'ascii')
    data_input_BAM_10 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_1.0virial_interp_BehrooziAM.dat', format = 'ascii')
    
    #data_input_ver1a = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.6virial_interp_ver2.dat', format = 'ascii')
    #data_input_ver2 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.6virial_interp_ver2rand.dat', format = 'ascii')
    #data_input_ver3 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.6virial_interp_ver3.dat', format = 'ascii')
    #data_input_ver4 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.6virial_interp_ver4.dat', format = 'ascii')
    #data_input_BAM = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.6virial_interp_BehrooziAM.dat', format = 'ascii')
    
    #quenching time data############################################
    tq_ver1_2 = np.array(data_input_ver1_2['QuenchTime'])*u.Gyr
    tq_ver1_3 = np.array(data_input_ver1_3['QuenchTime'])*u.Gyr
    tq_ver1_4 = np.array(data_input_ver1_4['QuenchTime'])*u.Gyr
    tq_ver1_5 = np.array(data_input_ver1_5['QuenchTime'])*u.Gyr
    tq_ver1_6 = np.array(data_input_ver1_6['QuenchTime'])*u.Gyr
    tq_ver1_7 = np.array(data_input_ver1_7['QuenchTime'])*u.Gyr
    tq_ver1_8 = np.array(data_input_ver1_8['QuenchTime'])*u.Gyr
    tq_ver1_9 = np.array(data_input_ver1_9['QuenchTime'])*u.Gyr
    tq_ver1_10 = np.array(data_input_ver1_10['QuenchTime'])*u.Gyr
    
    #tq_ver1a_2 = np.array(data_input_ver1a_2['QuenchTime'])*u.Gyr
    #tq_ver1a_3 = np.array(data_input_ver1a_3['QuenchTime'])*u.Gyr
    #tq_ver1a_4 = np.array(data_input_ver1a_4['QuenchTime'])*u.Gyr
    #tq_ver1a_5 = np.array(data_input_ver1a_5['QuenchTime'])*u.Gyr
    #tq_ver1a_6 = np.array(data_input_ver1a_6['QuenchTime'])*u.Gyr
    #tq_ver1a_7 = np.array(data_input_ver1a_7['QuenchTime'])*u.Gyr
    #tq_ver1a_8 = np.array(data_input_ver1a_8['QuenchTime'])*u.Gyr
    #tq_ver1a_9 = np.array(data_input_ver1a_9['QuenchTime'])*u.Gyr
    #tq_ver1a_10 = np.array(data_input_ver1a_10['QuenchTime'])*u.Gyr
    
    tq_ver2_2 = np.array(data_input_ver2_2['QuenchTime'])*u.Gyr
    tq_ver2_3 = np.array(data_input_ver2_3['QuenchTime'])*u.Gyr
    tq_ver2_4 = np.array(data_input_ver2_4['QuenchTime'])*u.Gyr
    tq_ver2_5 = np.array(data_input_ver2_5['QuenchTime'])*u.Gyr
    tq_ver2_6 = np.array(data_input_ver2_6['QuenchTime'])*u.Gyr
    tq_ver2_7 = np.array(data_input_ver2_7['QuenchTime'])*u.Gyr
    tq_ver2_8 = np.array(data_input_ver2_8['QuenchTime'])*u.Gyr
    tq_ver2_9 = np.array(data_input_ver2_9['QuenchTime'])*u.Gyr
    tq_ver2_10 = np.array(data_input_ver2_10['QuenchTime'])*u.Gyr
    
    tq_BAM_2 = np.array(data_input_BAM_2['QuenchTime'])*u.Gyr
    tq_BAM_3 = np.array(data_input_BAM_3['QuenchTime'])*u.Gyr
    tq_BAM_4 = np.array(data_input_BAM_4['QuenchTime'])*u.Gyr
    tq_BAM_5 = np.array(data_input_BAM_5['QuenchTime'])*u.Gyr
    tq_BAM_6 = np.array(data_input_BAM_6['QuenchTime'])*u.Gyr
    tq_BAM_7 = np.array(data_input_BAM_7['QuenchTime'])*u.Gyr
    tq_BAM_8 = np.array(data_input_BAM_8['QuenchTime'])*u.Gyr
    tq_BAM_9 = np.array(data_input_BAM_9['QuenchTime'])*u.Gyr
    tq_BAM_10 = np.array(data_input_BAM_10['QuenchTime'])*u.Gyr
    
    #redshift data################################################
    z_ver1_2 = np.array(data_input_ver1_2['zInfall'])
    z_ver1_3 = np.array(data_input_ver1_3['zInfall'])
    z_ver1_4 = np.array(data_input_ver1_4['zInfall'])
    z_ver1_5 = np.array(data_input_ver1_5['zInfall'])
    z_ver1_6 = np.array(data_input_ver1_6['zInfall'])
    z_ver1_7 = np.array(data_input_ver1_7['zInfall'])
    z_ver1_8 = np.array(data_input_ver1_8['zInfall'])
    z_ver1_9 = np.array(data_input_ver1_9['zInfall'])
    z_ver1_10 = np.array(data_input_ver1_10['zInfall'])
    
    #z_ver1a_2 = np.array(data_input_ver1a_2['zInfall'])
    #z_ver1a_3 = np.array(data_input_ver1a_3['zInfall'])
    #z_ver1a_4 = np.array(data_input_ver1a_4['zInfall'])
    #z_ver1a_5 = np.array(data_input_ver1a_5['zInfall'])
    #z_ver1a_6 = np.array(data_input_ver1a_6['zInfall'])
    #z_ver1a_7 = np.array(data_input_ver1a_7['zInfall'])
    #z_ver1a_8 = np.array(data_input_ver1a_8['zInfall'])
    #z_ver1a_9 = np.array(data_input_ver1a_9['zInfall'])
    #z_ver1a_10 = np.array(data_input_ver1a_10['zInfall'])
    
    z_ver2_2 = np.array(data_input_ver2_2['zInfall'])
    z_ver2_3 = np.array(data_input_ver2_3['zInfall'])
    z_ver2_4 = np.array(data_input_ver2_4['zInfall'])
    z_ver2_5 = np.array(data_input_ver2_5['zInfall'])
    z_ver2_6 = np.array(data_input_ver2_6['zInfall'])
    z_ver2_7 = np.array(data_input_ver2_7['zInfall'])
    z_ver2_8 = np.array(data_input_ver2_8['zInfall'])
    z_ver2_9 = np.array(data_input_ver2_9['zInfall'])
    z_ver2_10 = np.array(data_input_ver2_10['zInfall'])
    
    z_BAM_2 = np.array(data_input_BAM_2['zInfall'])
    z_BAM_3 = np.array(data_input_BAM_3['zInfall'])
    z_BAM_4 = np.array(data_input_BAM_4['zInfall'])
    z_BAM_5 = np.array(data_input_BAM_5['zInfall'])
    z_BAM_6 = np.array(data_input_BAM_6['zInfall'])
    z_BAM_7 = np.array(data_input_BAM_7['zInfall'])
    z_BAM_8 = np.array(data_input_BAM_8['zInfall'])
    z_BAM_9 = np.array(data_input_BAM_9['zInfall'])
    z_BAM_10 = np.array(data_input_BAM_10['zInfall'])
    
    #Determine the infall time from the redshift##############################
    t_ver1_2 = cosmo.lookback_time(z_ver1_2)
    t_ver1_3 = cosmo.lookback_time(z_ver1_3)
    t_ver1_4 = cosmo.lookback_time(z_ver1_4)
    t_ver1_5 = cosmo.lookback_time(z_ver1_5)
    t_ver1_6 = cosmo.lookback_time(z_ver1_6)
    t_ver1_7 = cosmo.lookback_time(z_ver1_7)
    t_ver1_8 = cosmo.lookback_time(z_ver1_8)
    t_ver1_9 = cosmo.lookback_time(z_ver1_9)
    t_ver1_10 = cosmo.lookback_time(z_ver1_10)
    
    #t_ver1a_2 = cosmo.lookback_time(z_ver1a_2)
    #t_ver1a_3 = cosmo.lookback_time(z_ver1a_3)
    #t_ver1a_4 = cosmo.lookback_time(z_ver1a_4)
    #t_ver1a_5 = cosmo.lookback_time(z_ver1a_5)
    #t_ver1a_6 = cosmo.lookback_time(z_ver1a_6)
    #t_ver1a_7 = cosmo.lookback_time(z_ver1a_7)
    #t_ver1a_8 = cosmo.lookback_time(z_ver1a_8)
    #t_ver1a_9 = cosmo.lookback_time(z_ver1a_9)
    #t_ver1a_10 = cosmo.lookback_time(z_ver1a_10)
    
    t_ver2_2 = cosmo.lookback_time(z_ver2_2)
    t_ver2_3 = cosmo.lookback_time(z_ver2_3)
    t_ver2_4 = cosmo.lookback_time(z_ver2_4)
    t_ver2_5 = cosmo.lookback_time(z_ver2_5)
    t_ver2_6 = cosmo.lookback_time(z_ver2_6)
    t_ver2_7 = cosmo.lookback_time(z_ver2_7)
    t_ver2_8 = cosmo.lookback_time(z_ver2_8)
    t_ver2_9 = cosmo.lookback_time(z_ver2_9)
    t_ver2_10 = cosmo.lookback_time(z_ver2_10)
    
    t_BAM_2 = cosmo.lookback_time(z_BAM_2)
    t_BAM_3 = cosmo.lookback_time(z_BAM_3)
    t_BAM_4 = cosmo.lookback_time(z_BAM_4)
    t_BAM_5 = cosmo.lookback_time(z_BAM_5)
    t_BAM_6 = cosmo.lookback_time(z_BAM_6)
    t_BAM_7 = cosmo.lookback_time(z_BAM_7)
    t_BAM_8 = cosmo.lookback_time(z_BAM_8)
    t_BAM_9 = cosmo.lookback_time(z_BAM_9)
    t_BAM_10 = cosmo.lookback_time(z_BAM_10)
    
    #determine the quenching timescale#################################
    tau_ver1_2 = t_ver1_2 - tq_ver1_2
    tau_ver1_3 = t_ver1_3 - tq_ver1_3
    tau_ver1_4 = t_ver1_4 - tq_ver1_4
    tau_ver1_5 = t_ver1_5 - tq_ver1_5
    tau_ver1_6 = t_ver1_6 - tq_ver1_6
    tau_ver1_7 = t_ver1_7 - tq_ver1_7
    tau_ver1_8 = t_ver1_8 - tq_ver1_8
    tau_ver1_9 = t_ver1_9 - tq_ver1_9
    tau_ver1_10 = t_ver1_10 - tq_ver1_10
    
    #tau_ver1a_2 = t_ver1a_2 - tq_ver1a_2
    #tau_ver1a_3 = t_ver1a_3 - tq_ver1a_3
    #tau_ver1a_4 = t_ver1a_4 - tq_ver1a_4
    #tau_ver1a_5 = t_ver1a_5 - tq_ver1a_5
    #tau_ver1a_6 = t_ver1a_6 - tq_ver1a_6
    #tau_ver1a_7 = t_ver1a_7 - tq_ver1a_7
    #tau_ver1a_8 = t_ver1a_8 - tq_ver1a_8
    #tau_ver1a_9 = t_ver1a_9 - tq_ver1a_9
    #tau_ver1a_10 = t_ver1a_10 - tq_ver1a_10
    
    tau_ver2_2 = t_ver2_2 - tq_ver2_2
    tau_ver2_3 = t_ver2_3 - tq_ver2_3
    tau_ver2_4 = t_ver2_4 - tq_ver2_4
    tau_ver2_5 = t_ver2_5 - tq_ver2_5
    tau_ver2_6 = t_ver2_6 - tq_ver2_6
    tau_ver2_7 = t_ver2_7 - tq_ver2_7
    tau_ver2_8 = t_ver2_8 - tq_ver2_8
    tau_ver2_9 = t_ver2_9 - tq_ver2_9
    tau_ver2_10 = t_ver2_10 - tq_ver2_10
    
    tau_BAM_2 = t_BAM_2 - tq_BAM_2
    tau_BAM_3 = t_BAM_3 - tq_BAM_3
    tau_BAM_4 = t_BAM_4 - tq_BAM_4
    tau_BAM_5 = t_BAM_5 - tq_BAM_5
    tau_BAM_6 = t_BAM_6 - tq_BAM_6
    tau_BAM_7 = t_BAM_7 - tq_BAM_7
    tau_BAM_8 = t_BAM_8 - tq_BAM_8
    tau_BAM_9 = t_BAM_9 - tq_BAM_9
    tau_BAM_10 = t_BAM_10 - tq_BAM_10
    

    
    mean_ver1 = np.array([np.mean(tau_ver1_2)/u.Gyr,np.mean(tau_ver1_3)/u.Gyr,np.mean(tau_ver1_4)/u.Gyr,np.mean(tau_ver1_5)/u.Gyr,np.mean(tau_ver1_6)/u.Gyr,np.mean(tau_ver1_7)/u.Gyr,np.mean(tau_ver1_8)/u.Gyr,np.mean(tau_ver1_9)/u.Gyr,np.mean(tau_ver1_10)/u.Gyr])
    #mean_ver1a = np.array([np.mean(tau_ver1a_2),np.mean(tau_ver1a_3),np.mean(tau_ver1a_4),np.mean(tau_ver1a_5),np.mean(tau_ver1a_6),np.mean(tau_ver1a_7),np.mean(tau_ver1a_8),np.mean(tau_ver1a_9),np.mean(tau_ver1a_10)])
    
    mean_ver2 = np.array([np.mean(tau_ver2_2)/u.Gyr,np.mean(tau_ver2_3)/u.Gyr,np.mean(tau_ver2_4)/u.Gyr,np.mean(tau_ver2_5)/u.Gyr,np.mean(tau_ver2_6)/u.Gyr,np.mean(tau_ver2_7)/u.Gyr,np.mean(tau_ver2_8)/u.Gyr,np.mean(tau_ver2_9)/u.Gyr,np.mean(tau_ver2_10)/u.Gyr])
    
    mean_BAM = np.array([np.mean(tau_BAM_2)/u.Gyr,np.mean(tau_BAM_3)/u.Gyr,np.mean(tau_BAM_4)/u.Gyr,np.mean(tau_BAM_5)/u.Gyr,np.mean(tau_BAM_6)/u.Gyr,np.mean(tau_BAM_7)/u.Gyr,np.mean(tau_BAM_8)/u.Gyr,np.mean(tau_BAM_9)/u.Gyr,np.mean(tau_BAM_10)/u.Gyr])
    
    median_ver1 = np.array([np.median(tau_ver1_2)/u.Gyr,np.median(tau_ver1_3)/u.Gyr,np.median(tau_ver1_4)/u.Gyr,np.median(tau_ver1_5)/u.Gyr,np.median(tau_ver1_6)/u.Gyr,np.median(tau_ver1_7)/u.Gyr,np.median(tau_ver1_8)/u.Gyr,np.median(tau_ver1_9)/u.Gyr,np.median(tau_ver1_10)/u.Gyr])
    
    #median_ver1a = np.array([np.median(tau_ver1a_2),np.median(tau_ver1a_3),np.median(tau_ver1a_4),np.median(tau_ver1a_5),np.median(tau_ver1a_6),np.median(tau_ver1a_7),np.median(tau_ver1a_8),np.median(tau_ver1a_9),np.median(tau_ver1a_10)])
    
    median_ver2 = np.array([np.median(tau_ver2_2)/u.Gyr,np.median(tau_ver2_3)/u.Gyr,np.median(tau_ver2_4)/u.Gyr,np.median(tau_ver2_5)/u.Gyr,np.median(tau_ver2_6)/u.Gyr,np.median(tau_ver2_7)/u.Gyr,np.median(tau_ver2_8)/u.Gyr,np.median(tau_ver2_9)/u.Gyr,np.median(tau_ver2_10)/u.Gyr])
    
    median_BAM = np.array([np.median(tau_BAM_2)/u.Gyr,np.median(tau_BAM_3)/u.Gyr,np.median(tau_BAM_4)/u.Gyr,np.median(tau_BAM_5)/u.Gyr,np.median(tau_BAM_6)/u.Gyr,np.median(tau_BAM_7)/u.Gyr,np.median(tau_BAM_8)/u.Gyr,np.median(tau_BAM_9)/u.Gyr,np.median(tau_BAM_10)/u.Gyr])
    
    error = np.array([np.std(tau_ver1_2)/u.Gyr,np.std(tau_ver1_3)/u.Gyr,np.std(tau_ver1_4)/u.Gyr,np.std(tau_ver1_5)/u.Gyr,np.std(tau_ver1_6)/u.Gyr,np.std(tau_ver1_7)/u.Gyr,np.std(tau_ver1_8)/u.Gyr,np.std(tau_ver1_9)/u.Gyr,np.std(tau_ver1_10)/u.Gyr])
    
    upper_error = mean_ver1 + error
    lower_error = mean_ver1 - error
    print upper_error
    print lower_error
    
    
    dist = np.array([50, 100, 150, 200, 250, 300])
    distfrac = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    quench = np.array([0.44, 0.72, 0.85, 0.94, 0.98, 0.99])
    quench_interp = np.array([0.469, 0.7579, 0.8704, 0.9535, 0.98, 0.99])
    quench_interpfrac_ver1 = np.array([0.6, 0.722, 0.795, 0.867, 0.910, 0.936, 0.965, 0.979, 1.0])
    #quench_interpfrac_ver1a = np.array([0.707, 0.894, 0.954, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
    quench_interpfrac_ver2 = np.array([0.627, 0.751, 0.836, 0.896, 0.927, 0.946, 0.974, 0.991, 1.0])
    quench_interpfrac_BAM = np.array([0.618, 0.754, 0.833, 0.886, 0.924, 0.949, 0.973, 0.986, 1.0])
    
    #####################################
    #Plot work
    #####################################
    axwidth = 3
    axlength = 10
    fontsize=24
    
    plt.rc('axes',linewidth=axwidth)  #make sure to do this before making any figures
    
    #<do plotting>
    
    plt.figure(figsize = (26,8))
    #plt.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9, wspace=0.2, hspace=0.2)
    #plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)    #plt.Subplotparams(top = 0.95)
    #plt.subplot(121)
    plt.subplot2grid((1,2),(0,0), colspan = 1)
    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.97, top=0.97, wspace=0.6, hspace=0.1)
    plt.axis([0, 1.0, 0, 1.1])
    ax = plt.gca()
    p1, = plt.plot(distfrac, quench_interpfrac_ver1, color='k', linestyle = '-', linewidth = 4.5, label = r'$\rm GK14\ AM$')
    #p1a, = plt.plot(distfrac, quench_interpfrac_ver1a, color = 'y', linestyle = '-', linewidth = 3.0)
    p2, = plt.plot(distfrac, quench_interpfrac_ver2, color='m', linestyle = '--', linewidth = 3.0, label = r'$\rm GK14\ AM\ +\ scatter$')
    pBAM, = plt.plot(distfrac, quench_interpfrac_BAM, color='c', linestyle = '--', linewidth = 3.0, label = r'$\rm Behroozi\ AM$')
    
    pp0 = plt.Rectangle((0, 0), 1, 1, fc="k", ec = "k")
    #pp1 = plt.Rectangle((0, 0), 1, 1, fc="Gold", ec = "Gold")
    pp2 = plt.Rectangle((0, 0), 1, 1, fc="m", ec = "m")
    pp3 = plt.Rectangle((0, 0), 1, 1, fc="c", ec="c")
    
    plt.legend(loc = 4, frameon=False, numpoints=1, prop={'size':16})
    
    #plt.legend([pp0,pp2,pp3], [r'$\rm GK14\ AM$', r'$\rm GK14\ AM\ +\ scatter$',  r'$\rm Behroozi\ AM$'], loc = 4, frameon=False, numpoints=1, prop={'size':16})
    
    
    ax.set_xlabel(r'$\rm Distance\ from\ Host\ (Rvir)$', fontsize = 24)
    ax.set_ylabel(r'$\rm Quenched\ Fraction$', fontsize = 24)
    
    xtickloc = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    xtickstr = ['$'+str(kk)+'$' for kk in xtickloc]
    ytickloc = [0.2, 0.4, 0.6, 0.8, 1.0]
    ytickstr = ['$'+str(kk)+'$' for kk in ytickloc]
    
    ax.set_xticks(xtickloc)
    ax.set_xticklabels(xtickstr)
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
    

    #plt.subplot(122)
    plt.subplot2grid((1,2),(0,1), colspan = 1)
    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.97, top=0.97, wspace=0.6, hspace=0.1)
    plt.axis([0,1,0,5])
    ax = plt.gca()

    plt.fill_between(distfrac, upper_error, lower_error, facecolor = 'grey', alpha = 0.4)
    p1_mean, = plt.plot(distfrac, mean_ver1, color='k', linestyle='-', linewidth = 4.5, label = r'$\rm GK14\ AM$')
    p1_med, = plt.plot(distfrac, median_ver1, color='k', linestyle=':', linewidth = 4.5)
    #p1a_mean, = plt.plot(distfrac, mean_ver1a, color = 'Gold', linestyle = '-', linewidth = 3.0)
    #p1a_median, = plt.plot(distfrac, median_ver1a, color = 'Gold', linestyle = '--', linewidth = 3.0)
    p2_mean, = plt.plot(distfrac, mean_ver2, color='m', linestyle = '--', linewidth = 3.0, label = r'$\rm GK14\ AM\ +\ scatter$')
    p2_median, = plt.plot(distfrac, median_ver2, color='m', linestyle = ':', linewidth = 3.0)
    pBAM_mean, = plt.plot(distfrac, mean_BAM, color='c', linestyle = '--', linewidth = 3.0, label = r'$\rm Behroozi\ AM$')
    pBAM_median, = plt.plot(distfrac, median_BAM, color='c', linestyle = ':', linewidth = 3.0)

    #plt.text(1.0, 1.0, r'$\rm Garrison-Kimmel\ AM$', color = 'k')
    #plt.text(1.0, 0.90, r'$\rm GK\ AM\ (<\ 0.5\ R_{vir})$', color = 'Gold')
    #plt.text(1.0, 0.95, r'$\rm Behroozi\ AM$', color = 'c')

    pp0 = plt.Rectangle((0, 0), 1, 1, fc="k", ec = "k")
    #pp1 = plt.Rectangle((0, 0), 1, 1, fc="Gold", ec = "Gold")
    pp2 = plt.Rectangle((0, 0), 1, 1, fc="m", ec = "m")
    pp3 = plt.Rectangle((0, 0), 1, 1, fc="c", ec="c")
    
    plt.legend(loc = 1, frameon=False, numpoints=1, prop={'size':16})
    #plt.legend([pp0,pp2,pp3], [r'$\rm GK14\ AM$', r'$\rm GK14\ AM\ +\ scatter$',  r'$\rm Behroozi\ AM$'], loc = 1, frameon=False, numpoints=1, prop={'size':16})
    
    #plt.legend([p1,p2,p3,p4,p5], [r'$\rm Vmax/Vpeak\ >\ 0.3$', r'$\rm Vmax/Vpeak\ >\ 0.4$', r'$\rm Vmax/Vpeak\ >\ 0.5$', r'$\rm Vmax/Vpeak\ >\ 0.6$', r'$\rm Vmax/Vpeak\ >\ 0.7$'], loc = 4, frameon=False, numpoints=1, prop={'size':11})
    
    ax.set_xlabel(r'$\rm Quenching\ Radius\ (R_{quench} / R_{vir})$', fontsize = 24)
    ax.set_ylabel(r'$\rm Quenching\ Time\ Scale\ (Gyr)$', fontsize = 24)
    
    xtickloc = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    xtickstr = ['$'+str(kk)+'$' for kk in xtickloc]
    ytickloc = [1,2,3,4,5]
    ytickstr = ['$'+str(kk)+'$' for kk in ytickloc]
    
    ax.set_xticks(xtickloc)
    ax.set_xticklabels(xtickstr)
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
    
    #plt.savefig('RadialModel/elvis_paperplot_radialquench_final2.pdf')
    
    
    plt.show()
    
