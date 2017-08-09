#creates a histogram of the elvis data
#number of dwarf galaxies in a clean file as a function of infall time

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy import cosmology as cosmo
from astropy import units as u

def histogram():

    subhalo_data = Table.read('ELVIS_Data_Clean/elvis_alldwarfs_clean.dat', format = 'ascii')

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
    master_distfrac1 = np.array([])
    master_velfrac1 = np.array([])
    master_zinfall1 = np.array([])

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

        master_distfrac1 = np.append(master_distfrac1, dist_frac)
        master_velfrac1 = np.append(master_velfrac1, vel_frac)
        master_zinfall1 = np.append(master_zinfall1, z_infall)

    subhalo_data_b = Table.read('ELVIS_Data_Clean/elvis_alldwarfs_clean_ver2rand.dat', format = 'ascii')

    h1dist_b = np.array(subhalo_data_b['h1dist(Mpc)'])
    h2dist_b = np.array(subhalo_data_b['h2dist(Mpc)'])
    h1firstz_b = np.array(subhalo_data_b['h1firstz_interacted'])
    h2firstz_b = np.array(subhalo_data_b['h2firstz_interacted'])
    vpeak_b = np.array(subhalo_data_b['Vpeak(km/s)'])
    vmax_b = np.array(subhalo_data_b['Vmax(km/s)'])
    hostname_b = np.array(subhalo_data_b['HostName'])
    R1vir_b = np.array(subhalo_data_b['R1vir(kpc)'])
    R2vir_b = np.array(subhalo_data_b['R2vir(kpc)'])

    #define empty velocity and distance columns
    master_distfrac2 = np.array([])
    master_velfrac2 = np.array([])
    master_zinfall2 = np.array([])

    #create fractional velocity, fractional distance, and infall time columns
    for i in range(len(h1dist_b)):

        dist1 = h1dist_b[i]*1000 #multiplying by 1000 changes units from Mpc to kpc
        dist2 = h2dist_b[i]*1000 #multiplying by 1000 changes units from Mpc to kpc
        
        cond_a = (h1firstz_b[i] > h2firstz_b[i]) & (dist1 < dist2) & (dist1 < R1vir_b[i])
        cond_b = (h1firstz_b[i] > h2firstz_b[i]) & (dist2 < dist1) & (dist2 < R2vir_b[i])
        cond_c = (h2firstz_b[i] > h1firstz_b[i]) & (dist2 < dist1) & (dist2 < R2vir_b[i])
        cond_d = (h2firstz_b[i] > h1firstz_b[i]) & (dist1 < dist2) & (dist1 < R1vir_b[i])
            
        if (cond_a or cond_d):
            distance = dist1
            z_infall = h1firstz_b[i]
            virial = R1vir_b[i]

        elif (cond_b or cond_c):
            distance = dist2
            z_infall = h2firstz_b[i]
            virial = R2vir_b[i]

        else:
            continue

        dist_frac_b = distance/virial
        vel_frac_b = vmax[i]/vpeak[i]

        master_distfrac2 = np.append(master_distfrac2, dist_frac_b)
        master_velfrac2 = np.append(master_velfrac2, vel_frac_b)
        master_zinfall2 = np.append(master_zinfall2, z_infall)
        
    #conditionlist = np.array([0.3,0.4,0.5,0.7])
    #conditionstring = np.array(['0.3','0.4','0.5','0.7'])

    #for i in range(len(conditionlist)):
        
    #eliminate the dwarfs that haven't fallen in, ie z = -1, and meet the other selection requirements
    #condition = conditionlist[i]
    #string = conditionstring[i]
    condition = 0.3
    string = '0.3'
    distcondition = 1.0
    distcondition2 = 0.5

    vel_cut1 = master_velfrac1 > condition
    vel_cut2 = master_velfrac2 > condition
    dist_cut1 = master_distfrac1 < distcondition
    dist_cut2 = master_distfrac2 < distcondition
    dist_cut1a = master_distfrac1 < distcondition2

    finalcut1 = vel_cut1 & dist_cut1 & (master_zinfall1 != -1)
    finalcut2 = vel_cut2 & dist_cut2 & (master_zinfall2 != -1)
    finalcut1a = vel_cut1 & dist_cut1a & (master_zinfall1 != -1)

    t_infall1 = cosmo.lookback_time(master_zinfall1[finalcut1])
    t_infall1a = cosmo.lookback_time(master_zinfall1[finalcut1a])
    t_infall2 = cosmo.lookback_time(master_zinfall2[finalcut2])

    plt.figure(figsize = (10,5), dpi = 100)
    # the histogram of the data
    #n, bins, patches = plt.hist(t_infall, 12, facecolor='blue', alpha=0.75)
    n, bins, patches = plt.hist(t_infall1, histtype = 'step', color = 'blue', bins = range(0,13,1), cumulative = -1, normed = True)
    n, bins, patches = plt.hist(t_infall1a, histtype = 'step', color = 'red', bins = range(0,13,1), cumulative = -1, normed = True)
    n, bins, patches = plt.hist(t_infall2, histtype = 'step', color = 'green', bins = range(0,13,1), cumulative = -1, normed = True)

    l = plt.axhline(y = 0.9, linewidth=1, linestyle = '--', color='m')

    
    plt.xlabel('Lookback Time (Gyr)')
    #plt.ylabel('Fraction of Dwarfs')
    histtitle = r'Cumulative Histogram of Time Since Infall'
    plt.title(histtitle)
    plt.axis([0, 12, 0, 1.1])
    plt.text(0.5,1.5, r'$\rm GK14AM,\ r = R_{vir}$', color = 'b')
    plt.text(0.5,1.45, r'$\rm Random,\ r = R_{vir}$', color = 'g')
    plt.text(0.5,1.4, r'$\rm GK14AM,\ r = 0.5R_{vir}$', color = 'r')
    plt.savefig('elvis_histogram_v'+string+'_paperfigure1.pdf')

    plt.show()
