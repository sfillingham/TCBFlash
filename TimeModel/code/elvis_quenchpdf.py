#This script determines the lookback time when a subhalo has been quenched, then plots a traditional or cumulative histogram of the lookback times.
#All subhalos need to be read in and which host halo they belong to (if they have fallen in at all).

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import astropy.cosmology as cosmo
import astropy.units as u

def traditional_pdf(subhalo_input, quench_timescale, v_ratio):

    #read in files and select relevant columns 
    subhalo_data = Table.read(subhalo_input, format = 'ascii')
    
    h1dist = np.array(subhalo_data['h1dist(Mpc)'])
    h2dist = np.array(subhalo_data['h2dist(Mpc)'])
    h1firstz = np.array(subhalo_data['h1firstz_interacted'])
    h2firstz = np.array(subhalo_data['h2firstz_interacted'])
    vpeak = np.array(subhalo_data['Vpeak(km/s)'])
    vmax = np.array(subhalo_data['Vmax(km/s)'])
    hostname = np.array(subhalo_data['HostName'])
    R1vir = np.array(subhalo_data['R1vir(kpc)'])
    R2vir = np.array(subhalo_data['R2vir(kpc)'])

    #define empty velocity and distance columns
    master_distfrac = np.array([])
    master_velfrac = np.array([])
    master_zinfall = np.array([])

    #create fractional velocity, fractional distance, and infall time columns
    for i in range(len(h1dist)):

        dist1 = h1dist[i]*1000 #multiplying by 1000 changes units from Mpc to kpc
        dist2 = h2dist[i]*1000 #multiplying by 1000 changes units from Mpc to kpc
        
        cond_a = (h1firstz[i] > h2firstz[i]) & (dist1 < dist2) & (dist1 < R1vir[i])
        cond_b = (h1firstz[i] > h2firstz[i]) & (dist2 < dist1) & (dist2 < R2vir[i])
        cond_c = (h2firstz[i] > h1firstz[i]) & (dist2 < dist1) & (dist2 < R2vir[i])
        cond_d = (h2firstz[i] > h1firstz[i]) & (dist1 < dist2) & (dist1 < R1vir[i])
            
        if (cond_a or cond_d):
            distance = dist1
            z_infall = h1firstz[i]
            virial = R1vir[i]

        elif (cond_b or cond_c):
            distance = dist2
            z_infall = h2firstz[i]
            virial = R2vir[i]

        else:
            continue

        dist_frac = distance/virial
        vel_frac = vmax[i]/vpeak[i]

        master_distfrac = np.append(master_distfrac, dist_frac)
        master_velfrac = np.append(master_velfrac, vel_frac)
        master_zinfall = np.append(master_zinfall, z_infall)

        #print(master_zinfall)

    #select velocity ratio and distance cuts
    vel_cut = master_velfrac > v_ratio
    dist_cut = master_distfrac < 1.0

    finalcut = vel_cut & dist_cut & (master_zinfall != -1)

    t_infall = cosmo.lookback_time(master_zinfall[finalcut])
    #print t_infall

    quench_cut = t_infall > quench_timescale*u.Gyr

    time_infall = t_infall[quench_cut]

    quenchtime = np.arange(len(time_infall))
    quenchtime[:] = quench_timescale
    tau_q = quenchtime*u.Gyr

    t_quench = time_infall - tau_q

    # the histogram of the data
    #n, bins, patches = plt.hist(t_quench, 10, facecolor='blue', alpha=0.75)
    n, bins, patches = plt.hist(t_quench, bins = range(0,11,1))

    plt.xlabel('Quenching Time')
    plt.ylabel('Number of Dwarfs')
    histtitle = r'Histogram of Quenching Time: Vmax/Vpeak > 0.4'
    plt.title(histtitle)
    plt.axis([0, 10, 0, 40])
    #plt.savefig('CheckBug/elvis_histogram_probDF.pdf')

    plt.show()

        
def cumulative_pdf(subhalo_input, quench_timescale, v_ratio):

    #read in files and select relevant columns 
    subhalo_data = Table.read(subhalo_input, format = 'ascii')
    
    h1dist = np.array(subhalo_data['h1dist(Mpc)'])
    h2dist = np.array(subhalo_data['h2dist(Mpc)'])
    h1firstz = np.array(subhalo_data['h1firstz_interacted'])
    h2firstz = np.array(subhalo_data['h2firstz_interacted'])
    vpeak = np.array(subhalo_data['Vpeak(km/s)'])
    vmax = np.array(subhalo_data['Vmax(km/s)'])
    hostname = np.array(subhalo_data['HostName'])
    R1vir = np.array(subhalo_data['R1vir(kpc)'])
    R2vir = np.array(subhalo_data['R2vir(kpc)'])

    #define empty velocity and distance columns
    master_distfrac = np.array([])
    master_velfrac = np.array([])
    master_zinfall = np.array([])

    #create fractional velocity, fractional distance, and infall time columns
    for i in range(len(h1dist)):

        dist1 = h1dist[i]*1000 #multiplying by 1000 changes units from Mpc to kpc
        dist2 = h2dist[i]*1000 #multiplying by 1000 changes units from Mpc to kpc
        
        cond_a = (h1firstz[i] > h2firstz[i]) & (dist1 < dist2) & (dist1 < R1vir[i])
        cond_b = (h1firstz[i] > h2firstz[i]) & (dist2 < dist1) & (dist2 < R2vir[i])
        cond_c = (h2firstz[i] > h1firstz[i]) & (dist2 < dist1) & (dist2 < R2vir[i])
        cond_d = (h2firstz[i] > h1firstz[i]) & (dist1 < dist2) & (dist1 < R1vir[i])
            
        if (cond_a or cond_d):
            distance = dist1
            z_infall = h1firstz[i]
            virial = R1vir[i]

        elif (cond_b or cond_c):
            distance = dist2
            z_infall = h2firstz[i]
            virial = R2vir[i]

        else:
            continue

        dist_frac = distance/virial
        vel_frac = vmax[i]/vpeak[i]

        master_distfrac = np.append(master_distfrac, dist_frac)
        master_velfrac = np.append(master_velfrac, vel_frac)
        master_zinfall = np.append(master_zinfall, z_infall)

        #print(master_zinfall)

    #select velocity ratio and distance cuts
    vel_cut = master_velfrac > v_ratio
    dist_cut = master_distfrac < 1.0

    finalcut = vel_cut & dist_cut & (master_zinfall != -1)

    t_infall = cosmo.lookback_time(master_zinfall[finalcut])
    #print t_infall

    quench_cut = t_infall > quench_timescale*u.Gyr

    time_infall = t_infall[quench_cut]

    quenchtime = np.arange(len(time_infall))
    quenchtime[:] = quench_timescale
    tau_q = quenchtime*u.Gyr

    t_quench = time_infall - tau_q
    
    print t_quench
    print len(t_quench)
    
    #set the values of the bin edges for the histogram
    bins = np.array([14.0 ,  12.5,  12.0 ,  11.5,  11.0 ,  10.5,  10.0 ,   9.5,   9.0 ,
                     8.5,   8.0 ,   7.5,   7.0 ,   6.5,   6.0 ,   5.5,   5.0 ,   4.5,
                     4.0 ,   3.5,   3.0 ,   2.5,   2.0 ,   1.5,   1.0 ,   0.5,   -1.0 ])

    #Introduce the Weisz data
    tq_eye = np.array([10.0, 1.5, 3.5, 1.0, 7.0, 5.0, 12.0, 3.0, 2.5])
    n, bins, patches = plt.hist(tq_eye, bins = range(0,13,1), histtype = 'step', color = 'green', linewidth = 4.0, cumulative = -1, normed = True)

    tq_10 = np.array([1.41, 1.26, 0.71, 0.71, 0.71, 0.79, 0.71, 3.16, 2.51])
    n, bins, patches = plt.hist(tq_10, bins = range(0,13,1), histtype = 'step', color = 'blue', linewidth = 4.0,cumulative = -1, normed = True)

    tq_95 = np.array([3.55, 2.51, 12.5, 1.67, 1.45, 2.54, 2.63, 3.28, 10.1])
    n, bins, patches = plt.hist(tq_95, bins = range(0,13,1), histtype = 'step', color = 'red', linewidth = 4.0,cumulative = -1, normed = True)

    tq_90 = np.array([7.12, 5.62, 12.6, 2.24, 1.67, 2.66, 3.35, 3.40, 10.6])
    n, bins, patches = plt.hist(tq_90, bins = range(0,13,1), histtype = 'step', color = 'DarkOrange', linewidth = 4.0,cumulative = -1, normed = True)

    tq_80 = np.array([7.33, 5,91, 12.6, 3.10, 2.04, 4.37, 5.90, 3.67, 11.3])
    n, bins, patches = plt.hist(tq_80, bins = range(0,13,1), histtype = 'step', color = 'magenta', linewidth = 4.0,cumulative = -1, normed = True)

    # the cumulative histogram of the data
    n, bins, patches = plt.hist(t_quench, bins = range(0,13,1), normed = True, color = 'black', linewidth = 4.0,cumulative = -1, histtype = 'step')
    
    final_array_eye = np.empty([1,len(bins)-1])
    final_array_10 = np.empty([1,len(bins)-1])
    final_array_95 = np.empty([1,len(bins)-1])
    final_array_90 = np.empty([1,len(bins)-1])
    final_array_80 = np.empty([1,len(bins)-1])
    final_array_elvis = np.empty([1,len(bins)-1])

    #initialize the number of subhalos in the first histogram bin
    cumtotal1 = 0
    cumtotal2 = 0
    cumtotal3 = 0
    cumtotal4 = 0
    cumtotal5 = 0
    cumtotal6 = 0

        
    for j in range(len(bins)-1):
            
        upper = bins[j]
        lower = bins[j+1]
        
        cut1 = (tq_eye < upper) & (tq_eye > lower)
        cut2 = (tq_10 < upper) & (tq_10 > lower)
        cut3 = (tq_95 < upper) & (tq_95 > lower)
        cut4 = (tq_90 < upper) & (tq_90 > lower)
        cut5 = (tq_80 < upper) & (tq_80 > lower)
        cut6 = (t_quench < upper) & (t_quench > lower)
            
        number1 = len(t_infall1[cut1])
        number2 = len(t_infall2[cut2])
        number3 = len(t_infall3[cut3])
        number4 = len(t_infall4[cut4])
        number5 = len(t_infall5[cut5])
        number6 = len(t_infall6[cut6])
            
        cumtotal1 = cumtotal1 + number1
        cumtotal2 = cumtotal2 + number2
        cumtotal3 = cumtotal3 + number3
        cumtotal4 = cumtotal4 + number4
        cumtotal5 = cumtotal5 + number5
        cumtotal6 = cumtotal6 + number6
            
        cumfrac1 = cumtotal1/float(len(t_infall1))
        cumfrac2 = cumtotal2/float(len(t_infall2))
        cumfrac3 = cumtotal3/float(len(t_infall3))
        cumfrac4 = cumtotal4/float(len(t_infall4))
        cumfrac5 = cumtotal5/float(len(t_infall5))
        cumfrac6 = cumtotal6/float(len(t_infall6))
            
        final_array_eye[0][j] = cumfrac1
        final_array_10[0][j] = cumfrac2
        final_array_95[0][j] = cumfrac3
        final_array_90[0][j] = cumfrac4
        final_array_80[0][j] = cumfrac5
        final_array_elvis[0][j] = cumfrac6

    plotbins = np.array([12.5,  12.0 ,  11.5,  11.0 ,  10.5,  10.0 ,   9.5,   9.0 ,
                     8.5,   8.0 ,   7.5,   7.0 ,   6.5,   6.0 ,   5.5,   5.0 ,   4.5,
                     4.0 ,   3.5,   3.0 ,   2.5,   2.0 ,   1.5,   1.0 ,   0.5,   0.0 ])
                     
                     
    p1 = plt.plot(plotbins, final_array_eye[0], color = 'g', linestyle = '--', linewidth = 4.0, label = r'$\rm Behroozi13$')
    p2 = plt.plot(plotbins, final_array_10[0], color = 'b', linestyle = '--', linewidth = 4.0, label = r'$\rm Behroozi13$')
    p3 = plt.plot(plotbins, final_array_95[0], color = 'r', linestyle = '--', linewidth = 4.0, label = r'$\rm Behroozi13$')
    p4 = plt.plot(plotbins, final_array_90[0], color = 'DarkOrange', linestyle = '--', linewidth = 4.0, label = r'$\rm Behroozi13$')
    p5 = plt.plot(plotbins, final_array_80[0], color = 'm', linestyle = '--', linewidth = 4.0, label = r'$\rm Behroozi13$')
    p6 = plt.plot(plotbins, final_array_elvis[0], color = 'k', linestyle = '--', linewidth = 4.0, label = r'$\rm Behroozi13$')
    
    plt.xlabel('Quenching Time (Gyr)')
    plt.ylabel('Fraction of Dwarfs')
    histtitle = r'Cumulative Histogram of Quenching Time: Vmax/Vpeak > 0.3'
    #plt.title(histtitle)
    plt.axis([0, 12, 0, 1.1])
    plt.text(8.0, 0.95, r'ELVIS data',color = 'k')
    #plt.text(8.0, 0.45, r'Weisz et al. 2014a data', color = 'k')
    plt.text(8.0, 0.91, r'by eye', color = 'g')
    plt.text(8.0, 0.87, r'stellar_frac = 1.0', color = 'b')
    plt.text(8.0, 0.83, r'stellar_frac = 0.95', color = 'r')
    plt.text(8.0, 0.79, r'stellar_frac = 0.90', color = 'DarkOrange')
    plt.text(8.0, 0.75, r'stellar_frac = 0.80', color = 'm')
    
    plt.savefig('paper2/elvis_cumhistogram_probDF_ver2.pdf')

    plt.show()

    
    
