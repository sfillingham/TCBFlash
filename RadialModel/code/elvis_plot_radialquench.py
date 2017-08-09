import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import astropy.cosmology as cosmo
import astropy.units as u


def radial():

    data50_input = Table.read('RadialModel/OutputData/elvis_alldwarfs_quenchtime_r50_1.0virial_interp.dat', format = 'ascii')
    data100_input = Table.read('RadialModel/OutputData/elvis_alldwarfs_quenchtime_r100_1.0virial_interp.dat', format = 'ascii')
    data150_input = Table.read('RadialModel/OutputData/elvis_alldwarfs_quenchtime_r150_1.0virial_interp.dat', format = 'ascii')
    data200_input = Table.read('RadialModel/OutputData/elvis_alldwarfs_quenchtime_r200_1.0virial_interp.dat', format = 'ascii')
    data250_input = Table.read('RadialModel/OutputData/elvis_alldwarfs_quenchtime_r250_1.0virial_interp.dat', format = 'ascii')
    data300_input = Table.read('RadialModel/OutputData/elvis_alldwarfs_quenchtime_r300_1.0virial_interp.dat', format = 'ascii')

    #shea50_input = Table.read('RadialModel/OutputData/elvis_alldwarfs_quenchtime_r50.dat', format = 'ascii')
    #shea100_input = Table.read('RadialModel/OutputData/elvis_alldwarfs_quenchtime_r100.dat', format = 'ascii')
    #shea150_input = Table.read('RadialModel/OutputData/elvis_alldwarfs_quenchtime_r150.dat', format = 'ascii')
    #shea200_input = Table.read('RadialModel/OutputData/elvis_alldwarfs_quenchtime_r200.dat', format = 'ascii')
    #shea250_input = Table.read('RadialModel/OutputData/elvis_alldwarfs_quenchtime_r250.dat', format = 'ascii')
    #shea300_input = Table.read('RadialModel/OutputData/elvis_alldwarfs_quenchtime_r300.dat', format = 'ascii')

    data50 = np.array(data50_input['QuenchTime'])
    data100 = np.array(data100_input['QuenchTime'])
    data150 = np.array(data150_input['QuenchTime'])
    data200 = np.array(data200_input['QuenchTime'])
    data250 = np.array(data250_input['QuenchTime'])
    data300 = np.array(data300_input['QuenchTime'])
    
    data_input_ver1 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.6virial_interp.dat', format = 'ascii')
    #data_input_ver2 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.6virial_interp_ver2rand.dat', format = 'ascii')
    #data_input_ver3 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.6virial_interp_ver3.dat', format = 'ascii')
    #data_input_ver4 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.6virial_interp_ver4.dat', format = 'ascii')
    data_input_BAM = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.6virial_interp_BehrooziAM.dat', format = 'ascii')
    
    tq_ver1 = np.array(data_input_ver1['QuenchTime'])*u.Gyr
    tq_BAM = np.array(data_input_BAM['QuenchTime'])*u.Gyr
    
    z_ver1 = np.array(data_input_ver1['zInfall'])
    z_BAM = np.array(data_input_BAM['zInfall'])
    
    t_ver1 = cosmo.lookback_time(z_ver1)
    t_BAM = cosmo.lookback_time(z_BAM)

    #data2_input = Table.read('RadialModel/OutputData/elvis_alldwarfs_quenchtime_r1000_0.2virial_interp.dat', format = 'ascii')
    #data3_input = Table.read('RadialModel/OutputData/elvis_alldwarfs_quenchtime_r1000_0.3virial_interp.dat', format = 'ascii')
    #data4_input = Table.read('RadialModel/OutputData/elvis_alldwarfs_quenchtime_r1000_0.4virial_interp.dat', format = 'ascii')
    #data5_input = Table.read('RadialModel/OutputData/elvis_alldwarfs_quenchtime_r1000_0.5virial_interp.dat', format = 'ascii')
    #data6_input = Table.read('RadialModel/OutputData/elvis_alldwarfs_quenchtime_r1000_0.6virial_interp.dat', format = 'ascii')
    #data7_input = Table.read('RadialModel/OutputData/elvis_alldwarfs_quenchtime_r1000_0.7virial_interp.dat', format = 'ascii')
    #data8_input = Table.read('RadialModel/OutputData/elvis_alldwarfs_quenchtime_r1000_0.8virial_interp.dat', format = 'ascii')
    #data9_input = Table.read('RadialModel/OutputData/elvis_alldwarfs_quenchtime_r1000_0.9virial_interp.dat', format = 'ascii')
    #data10_input = Table.read('RadialModel/OutputData/elvis_alldwarfs_quenchtime_r1000_1.0virial_interp.dat', format = 'ascii')

    #data2 = np.array(data2_input['QuenchTime'])
    #data3 = np.array(data3_input['QuenchTime'])
    #data4 = np.array(data4_input['QuenchTime'])
    #data5 = np.array(data5_input['QuenchTime'])
    #data6 = np.array(data6_input['QuenchTime'])
    #data7 = np.array(data7_input['QuenchTime'])
    #data8 = np.array(data8_input['QuenchTime'])
    #data9 = np.array(data9_input['QuenchTime'])
    #data10 = np.array(data10_input['QuenchTime'])
    
    dist = np.array([50, 100, 150, 200, 250, 300])
    distfrac = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    quench = np.array([0.44, 0.72, 0.85, 0.94, 0.98, 0.99])
    quench_interp = np.array([0.469, 0.7579, 0.8704, 0.9535, 0.98, 0.99])
    quench_interpfrac_ver1 = np.array([0.6, 0.722, 0.795, 0.867, 0.910, 0.936, 0.965, 0.979, 1.0])
    #quench_interpfrac_04 = np.array([0.584, 0.711, 0.789, 0.861, 0.905, 0.932, 0.963, 0.978, 1.0])
    #quench_interpfrac_05 = np.array([0.558, 0.693, 0.774, 0.853, 0.901, 0.929, 0.962, 0.977, 1.0])
    #quench_interpfrac_06 = np.array([0.525, 0.670, 0.758, 0.843, 0.896, 0.926, 0.961, 0.975, 1.0])
    #quench_interpfrac_07 = np.array([0.448, 0.603, 0.710, 0.808, 0.872, 0.909, 0.953, 0.969, 1.0])
    quench_interpfrac_BAM = np.array([0.618, 0.754, 0.833, 0.886, 0.924, 0.949, 0.973, 0.986, 1.0])

    #plt.figure(1)
    #plt.axis([0, 350, 0, 1.1])
    #plt.plot(dist, quench, 'bo-')
    #plt.plot(dist, quench_interp, 'ro-')
    #plt.xlabel('Distance from Host (kpc)')
    #plt.ylabel('Quenched Fraction')
    #plt.text(200.0, 0.5, r'Round 1', color = 'b')
    #plt.text(200.0, 0.45, r'Round 2: Cubic Interpolation', color = 'r')
    
    
    #plt.savefig('RadialModel/elvis_plot_radialquench.pdf')

    plt.figure(2)
    plt.subplot(211)
    plt.axis([0, 1.0, 0, 1.1])
    p1, = plt.plot(distfrac, quench_interpfrac_ver1, 'k-')
    #p2, = plt.plot(distfrac, quench_interpfrac_04, 'g-')
    #p3, = plt.plot(distfrac, quench_interpfrac_05, 'r-')
    p4, = plt.plot(distfrac, quench_interpfrac_BAM, 'c-')
    
    
    plt.xlabel('Distance from Host (Rvir)')
    plt.ylabel('Quenched Fraction')

    #plt.legend([p1,p2,p3,p4,p5], [r'$\rm Vmax/Vpeak\ >\ 0.3$', r'$\rm Vmax/Vpeak\ >\ 0.4$', r'$\rm Vmax/Vpeak\ >\ 0.5$', r'$\rm Vmax/Vpeak\ >\ 0.6$', r'$\rm Vmax/Vpeak\ >\ 0.7$'], loc = 4, frameon=False, numpoints=1, prop={'size':11})
    
    
    #plt.savefig('RadialModel/elvis_plot_radialquench_virialratio.pdf')


#plt.figure(3)
#   plt.axis([0,12,0,1.1])
#   plt.xlabel('Quenching Time (Gyr)')
#  plt.ylabel('Fraction of Dwarfs')

# n, bins, patches = plt.hist(data50, bins = range(0,13,1), histtype = 'step', color = 'green', linewidth = 1.5, cumulative = True, normed = True)
#   n, bins, patches = plt.hist(data100, bins = range(0,13,1), histtype = 'step', color = 'blue', linewidth = 1.5, cumulative = True, normed = True)
#   n, bins, patches = plt.hist(data150, bins = range(0,13,1), histtype = 'step', color = 'red', linewidth = 1.5, cumulative = True, normed = True)
#   n, bins, patches = plt.hist(data200, bins = range(0,13,1), histtype = 'step', color = 'DarkOrange', linewidth = 1.5, cumulative = True, normed = True)
#   n, bins, patches = plt.hist(data250, bins = range(0,13,1), histtype = 'step', color = 'cyan', linewidth = 1.5, cumulative = True, normed = True)
#   n, bins, patches = plt.hist(data300, bins = range(0,13,1), histtype = 'step', color = 'magenta', linewidth = 1.5, cumulative = True, normed = True)

#   tq_eye = np.array([10.0, 1.5, 3.5, 1.0, 7.0, 5.0, 12.0, 3.0, 2.5])
#  n, bins, patches = plt.hist(tq_eye, bins = range(0,13,1), histtype = 'step', color = 'black', linewidth = 2, cumulative = True, normed = True)

#  plt.text(8.0, 0.45, r'Weisz', color = 'k')
#  plt.text(8.0, 0.41, r'50 kpc', color = 'g')
#  plt.text(8.0, 0.37, r'100 kpc', color = 'b')
#  plt.text(8.0, 0.33, r'150 kpc', color = 'r')
#  plt.text(8.0, 0.29, r'200 kpc', color = 'DarkOrange')
#  plt.text(8.0, 0.25, r'250 kpc', color = 'c')
#  plt.text(8.0, 0.21, r'300 kpc', color = 'm')

    #plt.savefig('RadialModel/elvis_cumhistogram_probDF_interp_all.pdf')

#  plt.figure(4)
#  plt.axis([0,12,0,1.1])
#  plt.xlabel('Quenching Time (Gyr)')
#  plt.ylabel('Fraction of Dwarfs')

    #n, bins, patches = plt.hist(data2, bins = range(0,13,1), histtype = 'step', color = 'pink', linewidth = 1.5, cumulative = True, normed = True)
    #n, bins, patches = plt.hist(data3, bins = range(0,13,1), histtype = 'step', color = 'DarkBlue', linewidth = 1.5, cumulative = True, normed = True)
    #n, bins, patches = plt.hist(data4, bins = range(0,13,1), histtype = 'step', color = 'yellow', linewidth = 1.5, cumulative = True, normed = True)
    #n, bins, patches = plt.hist(data5, bins = range(0,13,1), histtype = 'step', color = 'green', linewidth = 1.5, cumulative = True, normed = True)
    #n, bins, patches = plt.hist(data6, bins = range(0,13,1), histtype = 'step', color = 'blue', linewidth = 1.5, cumulative = True, normed = True)
    #n, bins, patches = plt.hist(data7, bins = range(0,13,1), histtype = 'step', color = 'red', linewidth = 1.5, cumulative = True, normed = True)
    #n, bins, patches = plt.hist(data8, bins = range(0,13,1), histtype = 'step', color = 'DarkOrange', linewidth = 1.5, cumulative = True, normed = True)
    #n, bins, patches = plt.hist(data9, bins = range(0,13,1), histtype = 'step', color = 'cyan', linewidth = 1.5, cumulative = True, normed = True)
    #n, bins, patches = plt.hist(data10, bins = range(0,13,1), histtype = 'step', color = 'magenta', linewidth = 1.5, cumulative = True, normed = True)

#  tq_eye = np.array([10.0, 1.5, 3.5, 1.0, 7.0, 5.0, 12.0, 3.0, 2.5])
#   n, bins, patches = plt.hist(tq_eye, bins = range(0,13,1), histtype = 'step', color = 'black', linewidth = 2, cumulative = True, normed = True)

#  plt.text(8.0, 0.45, r'Weisz', color = 'k')
#  plt.text(8.0, 0.41, r'0.5 Rvir', color = 'g')
#  plt.text(8.0, 0.37, r'0.6 Rvir', color = 'b')
#   plt.text(8.0, 0.33, r'0.7 Rvir', color = 'r')
#  plt.text(8.0, 0.29, r'0.8 Rvir', color = 'DarkOrange')
    #plt.text(8.0, 0.25, r'0.9 Rvir', color = 'c')
    #plt.text(8.0, 0.21, r'1.0 Rvir', color = 'm')

    #plt.savefig('RadialModel/elvis_cumhistogram_probDF_interp_virialratio.pdf')
    plt.show()


def tau_pdf():

    #data50_input = Table.read('RadialModel/OutputData/elvis_alldwarfs_quenchtime_r50_1.0virial_interp.dat', format = 'ascii')
    #data100_input = Table.read('RadialModel/OutputData/elvis_alldwarfs_quenchtime_r100_1.0virial_interp.dat', format = 'ascii')
    #data150_input = Table.read('RadialModel/OutputData/elvis_alldwarfs_quenchtime_r150_1.0virial_interp.dat', format = 'ascii')
    #data200_input = Table.read('RadialModel/OutputData/elvis_alldwarfs_quenchtime_r200_1.0virial_interp.dat', format = 'ascii')
    #data250_input = Table.read('RadialModel/OutputData/elvis_alldwarfs_quenchtime_r250_1.0virial_interp.dat', format = 'ascii')
    #data300_input = Table.read('RadialModel/OutputData/elvis_alldwarfs_quenchtime_r300_1.0virial_interp.dat', format = 'ascii')

    data3_input = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.6virial_interp.dat', format = 'ascii')
    data4_input = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.4_0.6virial_interp.dat', format = 'ascii')
    data5_input = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.5_0.6virial_interp.dat', format = 'ascii')
    data6_input = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.6_0.6virial_interp.dat', format = 'ascii')
    data7_input = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.7_0.6virial_interp.dat', format = 'ascii')

    data_input_ver1 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.6virial_interp.dat', format = 'ascii')
    #data_input_ver2 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.6virial_interp_ver2rand.dat', format = 'ascii')
    #data_input_ver3 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.6virial_interp_ver3.dat', format = 'ascii')
    #data_input_ver4 = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.6virial_interp_ver4.dat', format = 'ascii')
    data_input_BAM = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.6virial_interp_BehrooziAM.dat', format = 'ascii')
    

    #tq50 = np.array(data50_input['QuenchTime'])*u.Gyr
    #tq100 = np.array(data100_input['QuenchTime'])*u.Gyr
    #tq150 = np.array(data150_input['QuenchTime'])*u.Gyr
    #tq200 = np.array(data200_input['QuenchTime'])*u.Gyr
    #tq250 = np.array(data250_input['QuenchTime'])*u.Gyr
    #tq300 = np.array(data300_input['QuenchTime'])*u.Gyr

    tq3 = np.array(data3_input['QuenchTime'])*u.Gyr
    tq4 = np.array(data4_input['QuenchTime'])*u.Gyr
    tq5 = np.array(data5_input['QuenchTime'])*u.Gyr
    tq6 = np.array(data6_input['QuenchTime'])*u.Gyr
    tq7 = np.array(data7_input['QuenchTime'])*u.Gyr
    
    tq_ver1 = np.array(data_input_ver1['QuenchTime'])*u.Gyr
    tq_BAM = np.array(data_input_BAM['QuenchTime'])*u.Gyr

    #z50 = np.array(data50_input['zInfall'])
    #z100 = np.array(data100_input['zInfall'])
    #z150 = np.array(data150_input['zInfall'])
    #z200 = np.array(data200_input['zInfall'])
    #z250 = np.array(data250_input['zInfall'])
    #z300 = np.array(data300_input['zInfall'])

    z3 = np.array(data3_input['zInfall'])
    z4 = np.array(data4_input['zInfall'])
    z5 = np.array(data5_input['zInfall'])
    z6 = np.array(data6_input['zInfall'])
    z7 = np.array(data7_input['zInfall'])
    
    z_ver1 = np.array(data_input_ver1['zInfall'])
    z_BAM = np.array(data_input_BAM['zInfall'])

    #t50 = cosmo.lookback_time(z50)
    #t100 = cosmo.lookback_time(z100)
    #t150 = cosmo.lookback_time(z150)
    #t200 = cosmo.lookback_time(z200)
    #t250 = cosmo.lookback_time(z250)
    #t300 = cosmo.lookback_time(z300)

    t3 = cosmo.lookback_time(z3)
    t4 = cosmo.lookback_time(z4)
    t5 = cosmo.lookback_time(z5)
    t6 = cosmo.lookback_time(z6)
    t7 = cosmo.lookback_time(z7)
    
    t_ver1 = cosmo.lookback_time(z_ver1)
    t_BAM = cosmo.lookback_time(z_BAM)

    #tau50 = t50/u.Gyr - tq50/u.Gyr
    #tau100 = t100/u.Gyr - tq100/u.Gyr
    #tau150 = t150/u.Gyr - tq150/u.Gyr
    #tau200 = t200/u.Gyr - tq200/u.Gyr
    #tau250 = t250/u.Gyr - tq250/u.Gyr
    #tau300 = t300/u.Gyr - tq300/u.Gyr

    tau3 = t3 - tq3
    tau4 = t4 - tq4
    tau5 = t5 - tq5
    tau6 = t6 - tq6
    tau7 = t7 - tq7
    
    tau_ver1 = t_ver1 - tq_ver1
    tau_BAM = t_BAM - tq_BAM

    #tau50_median = np.median(tau50)
    #tau100_median = np.median(tau100)
    #tau150_median = np.median(tau150)
    #tau200_median = np.median(tau200)
    #tau250_median = np.median(tau250)
    #tau300_median = np.median(tau300)

    #tau3_median = np.median(tau3)
    #tau4_median = np.median(tau4)
    #tau5_median = np.median(tau5)
    #tau6_median = np.median(tau6)
    #tau7_median = np.median(tau7)

    #ratio model

    plt.figure(10)
    plt.axis([-1,6,0,1.1])
    plt.xlabel('Quenching Time Scale (Gyr)')
    plt.ylabel('Fraction of Dwarfs')
    
    n, bins, patches = plt.hist(tau3, bins = range(-1,13,1), histtype = 'step', color = 'green', linewidth = 1.5, cumulative = True, normed = True)
    n, bins, patches = plt.hist(tau4, bins = range(-1,13,1), histtype = 'step', color = 'red', linewidth = 1.5, cumulative = True, normed = True)
    n, bins, patches = plt.hist(tau5, bins = range(-1,13,1), histtype = 'step', color = 'blue', linewidth = 1.5, cumulative = True, normed = True)
    n, bins, patches = plt.hist(tau6, bins = range(-1,13,1), histtype = 'step', color = 'black', linewidth = 1.5, cumulative = True, normed = True)
    n, bins, patches = plt.hist(tau7, bins = range(-1,13,1), histtype = 'step', color = 'cyan', linewidth = 1.5, cumulative = True, normed = True)
    

    plt.text(3.0, 0.49, r'$\rm V_{max}/V_{peak}\ >\ 0.3$', color = 'g')
    plt.text(3.0, 0.45, r'$\rm V_{max}/V_{peak}\ >\ 0.4$', color = 'b')
    plt.text(3.0, 0.41, r'$\rm V_{max}/V_{peak}\ >\ 0.5$', color = 'r')
    plt.text(3.0, 0.37, r'$\rm V_{max}/V_{peak}\ >\ 0.6$', color = 'k')
    plt.text(3.0, 0.33, r'$\rm V_{max}/V_{peak}\ >\ 0.7$', color = 'c')

    
    #plt.savefig('RadialModel/elvis_cumhistogram_tauquench_r0.6_radialmodel.pdf')
    plt.show()

    plt.figure(9)
    plt.axis([-1,6,0,1.1])
    plt.xlabel('Quenching Time Scale (Gyr)')
    plt.ylabel('Fraction of Dwarfs')

    n, bins, patches = plt.hist(tau3, bins = range(-1,6,1), histtype = 'step', color = 'green', linewidth = 1.5, normed = True)
    n, bins, patches = plt.hist(tau4, bins = range(-1,6,1), histtype = 'step', color = 'red', linewidth = 1.5, normed = True)
    n, bins, patches = plt.hist(tau5, bins = range(-1,6,1), histtype = 'step', color = 'blue', linewidth = 1.5, normed = True)
    n, bins, patches = plt.hist(tau6, bins = range(-1,6,1), histtype = 'step', color = 'black', linewidth = 1.5, normed = True)
    n, bins, patches = plt.hist(tau7, bins = range(-1,6,1), histtype = 'step', color = 'cyan', linewidth = 1.5, normed = True)

    plt.text(3.0, 0.49, r'$\rm V_{max}/V_{peak}\ >\ 0.3$', color = 'g')
    plt.text(3.0, 0.45, r'$\rm V_{max}/V_{peak}\ >\ 0.4$', color = 'b')
    plt.text(3.0, 0.41, r'$\rm V_{max}/V_{peak}\ >\ 0.5$', color = 'r')
    plt.text(3.0, 0.37, r'$\rm V_{max}/V_{peak}\ >\ 0.6$', color = 'k')
    plt.text(3.0, 0.33, r'$\rm V_{max}/V_{peak}\ >\ 0.7$', color = 'c')

    #plt.savefig('RadialModel/elvis_histogram_tauquench_r0.6_radialmodel.pdf')
    plt.show()

    plt.figure(11)
    plt.axis([0,12,0,70])
    plt.xlabel('Infall Time and Quenching Time (Gyr)')
    plt.ylabel('Number of Dwarfs')

    n, bins, patches = plt.hist(t3, bins = range(0,13,1), histtype = 'step', color = 'green', linewidth = 1.5)
    n, bins, patches = plt.hist(tq3, bins = range(0,13,1), histtype = 'step', color = 'green', linewidth = 1.5, linestyle = 'dashed')

    plt.text(1.0, 65, r'$\rm V_{max}/V_{peak}\ >\ 0.3$', color = 'g')
    plt.text(1.0, 60, r'solid = infall time', color = 'k')
    plt.text(1.0, 55, r'dashed = quenching lookback time', color = 'k')
    
    #plt.savefig('RadialModel/elvis_taudistribution_r0.6_v0.3_radialmodel.pdf')
    plt.show()

    plt.figure(99)
    plt.axis([-1,6,0,1.1])
    plt.xlabel('Quenching Time Scale (Gyr)')
    plt.ylabel('Fraction of Dwarfs')
    
    n, bins, patches = plt.hist(tau_ver1, bins = range(-1,6,1), histtype = 'step', color = 'green', linewidth = 1.5, normed = True)
    n, bins, patches = plt.hist(tau_BAM, bins = range(-1,6,1), histtype = 'step', color = 'red', linewidth = 1.5, normed = True)

    
    plt.text(3.0, 0.49, r'$\rm GK14AM$', color = 'g')
    plt.text(3.0, 0.41, r'$\rm BAM$', color = 'r')
    
    #plt.savefig('RadialModel/elvis_histogram_tauquench_r0.6_radialmodel.pdf')
    plt.show()


#plt.figure(12)
#plt.axis([0,12,0,70])
# plt.xlabel('Infall Time and Quenching Time (Gyr)')
# plt.ylabel('Number of Dwarfs')

#n, bins, patches = plt.hist(t4, bins = range(0,13,1), histtype = 'step', color = 'red', linewidth = 1.5)
#n, bins, patches = plt.hist(tq4, bins = range(0,13,1), histtype = 'step', color = 'red', linewidth = 1.5, linestyle = 'dashed')

#plt.text(1.0, 65, r'$\rm V_{max}/V_{peak}\ >\ 0.4$', color = 'r')
#plt.text(1.0, 60, r'solid = infall time', color = 'k')
#plt.text(1.0, 55, r'dashed = quenching lookback time', color = 'k')

#plt.figure(13)
#plt.axis([0,12,0,70])
#plt.xlabel('Infall Time and Quenching Time (Gyr)')
# plt.ylabel('Number of Dwarfs')

#n, bins, patches = plt.hist(t5, bins = range(0,13,1), histtype = 'step', color = 'blue', linewidth = 1.5)
#n, bins, patches = plt.hist(tq5, bins = range(0,13,1), histtype = 'step', color = 'blue', linewidth = 1.5, linestyle = 'dashed')

#plt.text(1.0, 65, r'$\rm V_{max}/V_{peak}\ >\ 0.5$', color = 'blue')
#plt.text(1.0, 60, r'solid = infall time', color = 'k')
#plt.text(1.0, 55, r'dashed = quenching lookback time', color = 'k')

    #plt.show()

#plt.figure(14)
#plt.axis([0,12,0,70])
#  plt.xlabel('Infall Time and Quenching Time (Gyr)')
#   plt.ylabel('Number of Dwarfs')

#  n, bins, patches = plt.hist(t6, bins = range(0,13,1), histtype = 'step', color = 'black', linewidth = 1.5)
#  n, bins, patches = plt.hist(tq6, bins = range(0,13,1), histtype = 'step', color = 'black', linewidth = 1.5, linestyle = 'dashed')

#   plt.text(1.0, 65, r'$\rm V_{max}/V_{peak}\ >\ 0.6$', color = 'k')
#  plt.text(1.0, 60, r'solid = infall time', color = 'k')
#   plt.text(1.0, 55, r'dashed = quenching lookback time', color = 'k')

    #plt.savefig('RadialModel/elvis_histogram_quenchtime_r0.6_radialmodel.pdf')

    #plt.show()

#   plt.figure(15)
#  plt.axis([0,12,0,70])
# plt.xlabel('Infall Time and Quenching Time (Gyr)')
# plt.ylabel('Number of Dwarfs')

#   n, bins, patches = plt.hist(t7, bins = range(0,13,1), histtype = 'step', color = 'cyan', linewidth = 1.5)
#   n, bins, patches = plt.hist(tq7, bins = range(0,13,1), histtype = 'step', color = 'cyan', linewidth = 1.5, linestyle = 'dashed')

#  plt.text(1.0, 65, r'$\rm V_{max}/V_{peak}\ >\ 0.7$', color = 'c')
#  plt.text(1.0, 60, r'solid = infall time', color = 'k')
# plt.text(1.0, 55, r'dashed = quenching lookback time', color = 'k')

#  plt.show()


def fqtq():

    data2_input = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.2virial_interp.dat', format = 'ascii')
    data3_input = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.3virial_interp.dat', format = 'ascii')
    data4_input = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.4virial_interp.dat', format = 'ascii')
    data5_input = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.5virial_interp.dat', format = 'ascii')
    data6_input = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.6virial_interp.dat', format = 'ascii')
    data7_input = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.7virial_interp.dat', format = 'ascii')
    data8_input = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.8virial_interp.dat', format = 'ascii')
    data9_input = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_0.9virial_interp.dat', format = 'ascii')
    data10_input = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r1000_v0.3_1.0virial_interp.dat', format = 'ascii')

    data2 = np.array(data2_input['QuenchTime'])
    data3 = np.array(data3_input['QuenchTime'])
    data4 = np.array(data4_input['QuenchTime'])
    data5 = np.array(data5_input['QuenchTime'])
    data6 = np.array(data6_input['QuenchTime'])
    data7 = np.array(data7_input['QuenchTime'])
    data8 = np.array(data8_input['QuenchTime'])
    data9 = np.array(data9_input['QuenchTime'])
    data10 = np.array(data10_input['QuenchTime'])

    data50_input = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r50_v0.3_1.0virial_interp.dat', format = 'ascii')
    data100_input = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r100_v0.3_1.0virial_interp.dat', format = 'ascii')
    data150_input = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r150_v0.3_1.0virial_interp.dat', format = 'ascii')
    data200_input = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r200_v0.3_1.0virial_interp.dat', format = 'ascii')
    data250_input = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r250_v0.3_1.0virial_interp.dat', format = 'ascii')
    data300_input = Table.read('RadialModel/elvis_alldwarfs_quenchtime_r300_v0.3_1.0virial_interp.dat', format = 'ascii')

    tq50 = np.array(data50_input['QuenchTime'])*u.Gyr
    tq100 = np.array(data100_input['QuenchTime'])*u.Gyr
    tq150 = np.array(data150_input['QuenchTime'])*u.Gyr
    tq200 = np.array(data200_input['QuenchTime'])*u.Gyr
    tq250 = np.array(data250_input['QuenchTime'])*u.Gyr
    tq300 = np.array(data300_input['QuenchTime'])*u.Gyr
    
    dist = np.array([50, 100, 150, 200, 250, 300])
    distfrac = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    quench_interp_v03 = np.array([0.485, 0.722, 0.839, 0.926, 0.956, 0.982])
    quench_interpfrac_v03 = np.array([0.584, 0.711, 0.789, 0.861, 0.905, 0.932, 0.963, 0.977, 1.0])

    tq2 = np.array(data2_input['QuenchTime'])*u.Gyr
    tq3 = np.array(data3_input['QuenchTime'])*u.Gyr
    tq4 = np.array(data4_input['QuenchTime'])*u.Gyr
    tq5 = np.array(data5_input['QuenchTime'])*u.Gyr
    tq6 = np.array(data6_input['QuenchTime'])*u.Gyr
    tq7 = np.array(data7_input['QuenchTime'])*u.Gyr
    tq8 = np.array(data8_input['QuenchTime'])*u.Gyr
    tq9 = np.array(data9_input['QuenchTime'])*u.Gyr
    tq10 = np.array(data10_input['QuenchTime'])*u.Gyr

    z50 = np.array(data50_input['zInfall'])
    z100 = np.array(data100_input['zInfall'])
    z150 = np.array(data150_input['zInfall'])
    z200 = np.array(data200_input['zInfall'])
    z250 = np.array(data250_input['zInfall'])
    z300 = np.array(data300_input['zInfall'])

    z2 = np.array(data2_input['zInfall'])
    z3 = np.array(data3_input['zInfall'])
    z4 = np.array(data4_input['zInfall'])
    z5 = np.array(data5_input['zInfall'])
    z6 = np.array(data6_input['zInfall'])
    z7 = np.array(data7_input['zInfall'])
    z8 = np.array(data8_input['zInfall'])
    z9 = np.array(data9_input['zInfall'])
    z10 = np.array(data10_input['zInfall'])

    t50 = cosmo.lookback_time(z50)
    t100 = cosmo.lookback_time(z100)
    t150 = cosmo.lookback_time(z150)
    t200 = cosmo.lookback_time(z200)
    t250 = cosmo.lookback_time(z250)
    t300 = cosmo.lookback_time(z300)

    t2 = cosmo.lookback_time(z2)
    t3 = cosmo.lookback_time(z3)
    t4 = cosmo.lookback_time(z4)
    t5 = cosmo.lookback_time(z5)
    t6 = cosmo.lookback_time(z6)
    t7 = cosmo.lookback_time(z7)
    t8 = cosmo.lookback_time(z8)
    t9 = cosmo.lookback_time(z9)
    t10 = cosmo.lookback_time(z10)

    tau50 = t50 - tq50
    tau100 = t100 - tq100
    tau150 = t150 - tq150
    tau200 = t200 - tq200
    tau250 = t250 - tq250
    tau300 = t300 - tq300

    tau2 = t2 - tq2
    tau3 = t3 - tq3
    tau4 = t4 - tq4
    tau5 = t5 - tq5
    tau6 = t6 - tq6
    tau7 = t7 - tq7
    tau8 = t8 - tq8
    tau9 = t9 - tq9
    tau10 = t10 - tq10

    tau50_median = np.median(tau50)/u.Gyr
    tau100_median = np.median(tau100)/u.Gyr
    tau150_median = np.median(tau150)/u.Gyr
    tau200_median = np.median(tau200)/u.Gyr
    tau250_median = np.median(tau250)/u.Gyr
    tau300_median = np.median(tau300)/u.Gyr

    tau2_median = np.median(tau2)/u.Gyr
    tau3_median = np.median(tau3)/u.Gyr
    tau4_median = np.median(tau4)/u.Gyr
    tau5_median = np.median(tau5)/u.Gyr
    tau6_median = np.median(tau6)/u.Gyr
    tau7_median = np.median(tau7)/u.Gyr
    tau8_median = np.median(tau8)/u.Gyr
    tau9_median = np.median(tau9)/u.Gyr
    tau10_median = np.median(tau10)/u.Gyr

    taulist_frac = np.array([tau2_median, tau3_median, tau4_median, tau5_median, tau6_median, tau7_median, tau8_median, tau9_median, tau10_median])
    taulist_dist = np.array([tau50_median, tau100_median, tau150_median, tau200_median, tau250_median, tau300_median])

    #data from time model
    timedata = Table.read('TimeModel/Output/elvis_fqtq_r1.0_v0.3.dat', format = 'ascii')
    timefq = np.array(timedata['fq'])
    timetq = np.array(timedata['tq'])

    plt.figure(1)
    plt.axis([0.2,1,0,12])
    plt.xlabel('quenched fraction')
    plt.ylabel('quenching timescale (Gyr)')
    plt.title(r'$\rm GK14AM:\ V_{max}/V_{peak}\ >\ 0.3$')
    plt.text(0.5,10.5, r'$\rm Time \, Model$', color = 'b')
    plt.text(0.5,10.1, r'$\rm Radial \, Model: \, virial \, frac$', color = 'g')
    plt.text(0.5,9.7, r'$\rm Radial \, Model: \, absolute \, distance$', color = 'r')
    
    plt.plot(timefq, timetq, 'b-', quench_interpfrac_v03, taulist_frac, 'g-', quench_interp_v03, taulist_dist, 'r-')

    plt.savefig('elvis_fqtqplot_allmodels.pdf')
    plt.show()

    
