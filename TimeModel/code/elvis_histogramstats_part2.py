#This script gets the 1-sigma region for the cumulative histogram of infall times for the ELVIS suite.
#It combines many features of the histogram.py and getstats.py scripts.
#The inputlist is '~/elvis/TimeModel/elvis_infalltime_input.dat'

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy import cosmology as cosmo
from astropy import units as u
from matplotlib.ticker import AutoMinorLocator

def stats(input_list):

    inputlist = Table.read(input_list, format = 'ascii')
    galname = np.array(inputlist['Galname'])

    #set the values of the bin edges for the histogram
    bins = np.array([14.0 ,  12.5,  12.0 ,  11.5,  11.0 ,  10.5,  10.0 ,   9.5,   9.0 ,
         8.5,   8.0 ,   7.5,   7.0 ,   6.5,   6.0 ,   5.5,   5.0 ,   4.5,
         4.0 ,   3.5,   3.0 ,   2.5,   2.0 ,   1.5,   1.0 ,   0.5,   -1.0 ])
    
    ###########################
    #Garrison-Kimmel Abundance Matching routine
    ###########################
    
    subhalo_data = Table.read('ELVIS_Data_Clean/elvis_alldwarfs_clean.dat', format = 'ascii')
    red_data = Table.read('TimeModel/elvis_allsubhalos_infallz_interp.dat', format = 'ascii')

    h1dist = np.array(subhalo_data['h1dist(Mpc)'])
    h2dist = np.array(subhalo_data['h2dist(Mpc)'])
    #h1firstz = np.array(subhalo_data['h1firstz_interacted'])
    #h2firstz = np.array(subhalo_data['h2firstz_interacted'])
    h1firstz = np.array(red_data['h1z_infall'])
    h2firstz = np.array(red_data['h2z_infall'])
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

    ####################################################################
    #Garrison-Kimmel Abundance Matching routine for the INDIVIDUAL hosts
    #####################################################################



    subhalo_data = Table.read('ELVIS_Data_Clean/elvis_'+host+'_clean.dat', format = 'ascii')
    red_data = Table.read('TimeModel/elvis_allsubhalos_infallz_interp.dat', format = 'ascii')
    
    h1dist = np.array(subhalo_data['h1dist(Mpc)'])
    h2dist = np.array(subhalo_data['h2dist(Mpc)'])
    #h1firstz = np.array(subhalo_data['h1firstz_interacted'])
    #h2firstz = np.array(subhalo_data['h2firstz_interacted'])
    h1firstz = np.array(red_data['h1z_infall'])
    h2firstz = np.array(red_data['h2z_infall'])
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

    ############################
    #Extending the Halo Mass range and randomly selecting 9 subhalos per host
    ############################
    subhalo_data_b = Table.read('ELVIS_Data_Clean/randomtest/elvis_alldwarfs_clean_ver2rand.dat', format = 'ascii')
    red_data_b = Table.read('TimeModel/elvis_allsubhalos_infallz_interp_ver2rand.dat', format = 'ascii')

    h1dist_b = np.array(subhalo_data_b['h1dist(Mpc)'])
    h2dist_b = np.array(subhalo_data_b['h2dist(Mpc)'])
    #h1firstz_b = np.array(subhalo_data_b['h1firstz_interacted'])
    #h2firstz_b = np.array(subhalo_data_b['h2firstz_interacted'])
    h1firstz_b = np.array(red_data_b['h1z_infall'])
    h2firstz_b = np.array(red_data_b['h2z_infall'])
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
        vel_frac_b = vmax_b[i]/vpeak_b[i]

        master_distfrac2 = np.append(master_distfrac2, dist_frac_b)
        master_velfrac2 = np.append(master_velfrac2, vel_frac_b)
        master_zinfall2 = np.append(master_zinfall2, z_infall)

        
    ############################
    #Behroozi Abundance Matching routine
    ############################
    subhalo_data_bz = Table.read('ELVIS_Data_Clean/elvis_alldwarfs_clean_BehrooziAM.dat', format = 'ascii')

    h1dist_bz = np.array(subhalo_data_bz['h1dist(Mpc)'])
    h2dist_bz = np.array(subhalo_data_bz['h2dist(Mpc)'])
    h1firstz_bz = np.array(subhalo_data_bz['h1firstz_interacted'])
    h2firstz_bz = np.array(subhalo_data_bz['h2firstz_interacted'])
    vpeak_bz = np.array(subhalo_data_bz['Vpeak(km/s)'])
    vmax_bz = np.array(subhalo_data_bz['Vmax(km/s)'])
    hostname_bz = np.array(subhalo_data_bz['HostName'])
    R1vir_bz = np.array(subhalo_data_bz['R1vir(kpc)'])
    R2vir_bz = np.array(subhalo_data_bz['R2vir(kpc)'])

    #define empty velocity and distance columns
    master_distfrac3 = np.array([])
    master_velfrac3 = np.array([])
    master_zinfall3 = np.array([])

    #create fractional velocity, fractional distance, and infall time columns
    for i in range(len(h1dist_bz)):

        dist1 = h1dist_bz[i]*1000 #multiplying by 1000 changes units from Mpc to kpc
        dist2 = h2dist_bz[i]*1000 #multiplying by 1000 changes units from Mpc to kpc
        
        cond_a = (h1firstz_bz[i] > h2firstz_bz[i]) & (dist1 < dist2) & (dist1 < R1vir_bz[i])
        cond_b = (h1firstz_bz[i] > h2firstz_bz[i]) & (dist2 < dist1) & (dist2 < R2vir_bz[i])
        cond_c = (h2firstz_bz[i] > h1firstz_bz[i]) & (dist2 < dist1) & (dist2 < R2vir_bz[i])
        cond_d = (h2firstz_bz[i] > h1firstz_bz[i]) & (dist1 < dist2) & (dist1 < R1vir_bz[i])
            
        if (cond_a or cond_d):
            distance = dist1
            z_infall = h1firstz_bz[i]
            virial = R1vir_bz[i]

        elif (cond_b or cond_c):
            distance = dist2
            z_infall = h2firstz_bz[i]
            virial = R2vir_bz[i]

        else:
            continue

        dist_frac_bz = distance/virial
        vel_frac_bz = vmax_bz[i]/vpeak_bz[i]

        master_distfrac3 = np.append(master_distfrac3, dist_frac_bz)
        master_velfrac3 = np.append(master_velfrac3, vel_frac_bz)
        master_zinfall3 = np.append(master_zinfall3, z_infall)


    ############################
    #Higher Mass Range, check on Coral's point
    ############################
    subhalo_data_4 = Table.read('ELVIS_Data_Clean/elvis_alldwarfs_clean_ver4.dat', format = 'ascii')
    red_data_4 = Table.read('TimeModel/elvis_allsubhalos_infallz_interp_ver4.dat', format = 'ascii')

    h1dist_4 = np.array(subhalo_data_4['h1dist(Mpc)'])
    h2dist_4 = np.array(subhalo_data_4['h2dist(Mpc)'])
    #h1firstz_4 = np.array(subhalo_data_4['h1firstz_interacted'])
    #h2firstz_4 = np.array(subhalo_data_4['h2firstz_interacted'])
    h1firstz_4 = np.array(red_data_4['h1z_infall'])
    h2firstz_4 = np.array(red_data_4['h2z_infall'])
    vpeak_4 = np.array(subhalo_data_4['Vpeak(km/s)'])
    vmax_4 = np.array(subhalo_data_4['Vmax(km/s)'])
    hostname_4 = np.array(subhalo_data_4['HostName'])
    R1vir_4 = np.array(subhalo_data_4['R1vir(kpc)'])
    R2vir_4 = np.array(subhalo_data_4['R2vir(kpc)'])

    #define empty velocity and distance columns
    master_distfrac4 = np.array([])
    master_velfrac4 = np.array([])
    master_zinfall4 = np.array([])

    #create fractional velocity, fractional distance, and infall time columns
    for i in range(len(h1dist_4)):

        dist1 = h1dist_4[i]*1000 #multiplying by 1000 changes units from Mpc to kpc
        dist2 = h2dist_4[i]*1000 #multiplying by 1000 changes units from Mpc to kpc
        
        cond_a = (h1firstz_4[i] > h2firstz_4[i]) & (dist1 < dist2) & (dist1 < R1vir_4[i])
        cond_b = (h1firstz_4[i] > h2firstz_4[i]) & (dist2 < dist1) & (dist2 < R2vir_4[i])
        cond_c = (h2firstz_4[i] > h1firstz_4[i]) & (dist2 < dist1) & (dist2 < R2vir_4[i])
        cond_d = (h2firstz_4[i] > h1firstz_4[i]) & (dist1 < dist2) & (dist1 < R1vir_4[i])
            
        if (cond_a or cond_d):
            distance = dist1
            z_infall = h1firstz_4[i]
            virial = R1vir_4[i]

        elif (cond_b or cond_c):
            distance = dist2
            z_infall = h2firstz_4[i]
            virial = R2vir_4[i]

        else:
            continue

        dist_frac_4 = distance/virial
        vel_frac_4 = vmax_4[i]/vpeak_4[i]

        master_distfrac4 = np.append(master_distfrac4, dist_frac_4)
        master_velfrac4 = np.append(master_velfrac4, vel_frac_4)
        master_zinfall4 = np.append(master_zinfall4, z_infall)

            
    #eliminate the dwarfs that haven't fallen in, ie z = -1, and meet the other selection requirements
    condition = 0.3
    condition2 = 0.7
    string = '0.3'
    distcondition = 1.0
    distcondition2 = 0.5

    vel_cut1 = master_velfrac1 > condition
    vel_cut2 = master_velfrac2 > condition
    vel_cut3 = master_velfrac3 > condition
    vel_cut1b = master_velfrac1 > condition2
    vel_cut4 = master_velfrac4 > condition
    dist_cut1 = master_distfrac1 < distcondition
    dist_cut1a = master_distfrac1 < distcondition2
    dist_cut2 = master_distfrac2 < distcondition
    dist_cut3 = master_distfrac3 < distcondition
    dist_cut4 = master_distfrac4 < distcondition

    finalcut1 = vel_cut1 & dist_cut1 & (master_zinfall1 != -1)
    finalcut1a = vel_cut1 & dist_cut1a & (master_zinfall1 != -1)
    finalcut1b = vel_cut1b & dist_cut1 & (master_zinfall1 != -1)
    finalcut2 = vel_cut2 & dist_cut2 & (master_zinfall2 != -1)
    finalcut3 = vel_cut3 & dist_cut3 & (master_zinfall3 != -1)
    finalcut4 = vel_cut4 & dist_cut4 & (master_zinfall4 != -1)

    t_infall1 = cosmo.lookback_time(master_zinfall1[finalcut1])/u.Gyr
    t_infall1a = cosmo.lookback_time(master_zinfall1[finalcut1a])/u.Gyr
    t_infall1b = cosmo.lookback_time(master_zinfall1[finalcut1b])/u.Gyr
    t_infall2 = cosmo.lookback_time(master_zinfall2[finalcut2])/u.Gyr
    t_infall3 = cosmo.lookback_time(master_zinfall3[finalcut3])/u.Gyr
    t_infall4 = cosmo.lookback_time(master_zinfall4[finalcut4])/u.Gyr


    #########################
    #The script now moves onto the statistics of individual cumulative distributions.
    #########################
    #Initialize final array which will contain the cumulative quenched fraction for each host.
    final_array = np.empty([len(galname),len(bins)-1])
    final_array1 = np.empty([1,len(bins)-1])
    final_array1a = np.empty([1,len(bins)-1])
    final_array1b = np.empty([1,len(bins)-1])
    final_array2 = np.empty([1,len(bins)-1])
    final_array3 = np.empty([1,len(bins)-1])
    final_array4 = np.empty([1,len(bins)-1])

    for i in range(len(galname)):

        name = galname[i]

        inputfile = Table.read('TimeModel/Output/'+name, format = 'ascii')
        time = np.array(inputfile['InfallTime'])

        #initialize the number of subhalos in the first histogram bin
        cumtotal = 0
        cumtotal1 = 0
        cumtotal1a = 0
        cumtotal1b = 0
        cumtotal2 = 0
        cumtotal3 = 0
        cumtotal4 = 0

        for j in range(len(bins)-1):

            upper = bins[j]
            lower = bins[j+1]

            cut = (time < upper) & (time > lower)
            cut1 = (t_infall1 < upper) & (t_infall1 > lower)
            cut1a = (t_infall1a < upper) & (t_infall1a > lower)
            cut1b = (t_infall1b < upper) & (t_infall1b > lower)
            cut2 = (t_infall2 < upper) & (t_infall2 > lower)
            cut3 = (t_infall3 < upper) & (t_infall3 > lower)
            cut4 = (t_infall4 < upper) & (t_infall4 > lower)
            

            number = len(time[cut])
            cumtotal = cumtotal + number
            cumfrac = cumtotal/float(len(time))
            final_array[i][j] = cumfrac

            number1 = len(t_infall1[cut1])
            number1a = len(t_infall1a[cut1a])
            number1b = len(t_infall1b[cut1b])
            number2 = len(t_infall2[cut2])
            number3 = len(t_infall3[cut3])
            number4 = len(t_infall4[cut4])

            cumtotal1 = cumtotal1 + number1
            cumtotal1a = cumtotal1a + number1a
            cumtotal1b = cumtotal1b + number1b
            cumtotal2 = cumtotal2 + number2
            cumtotal3 = cumtotal3 + number3
            cumtotal4 = cumtotal4 + number4

            cumfrac1 = cumtotal1/float(len(t_infall1))
            cumfrac1a = cumtotal1a/float(len(t_infall1a))
            cumfrac1b = cumtotal1b/float(len(t_infall1b))
            cumfrac2 = cumtotal2/float(len(t_infall2))
            cumfrac3 = cumtotal3/float(len(t_infall3))
            cumfrac4 = cumtotal4/float(len(t_infall4))

            final_array1[0][j] = cumfrac1
            final_array1a[0][j] = cumfrac1a
            final_array1b[0][j] = cumfrac1b
            final_array2[0][j] = cumfrac2
            final_array3[0][j] = cumfrac3
            final_array4[0][j] = cumfrac4

    ####The histogram of the data
    #n, bins, patches = plt.hist(t_infall1, histtype = 'step', color = 'blue', bins = [0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,10.5,11.0,11.5,12.0,12.5,13.0], cumulative = -1, normed = True)
    #n, bins, patches = plt.hist(t_infall1a, histtype = 'step', color = 'red', bins = [0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,10.5,11.0,11.5,12.0,12.5,13.0], cumulative = -1, normed = True)
    #n, bins, patches = plt.hist(t_infall2, histtype = 'step', color = 'black', bins = [0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,10.5,11.0,11.5,12.0,12.5,13.0], cumulative = -1, normed = True)

    #determine the range and plot the histogram

    standard = np.std(final_array, axis=0)
    standard2 = np.std(final_array, axis=0)
    median = np.median(final_array, axis = 0)

    plotbins = np.array([12.5,  12.0 ,  11.5,  11.0 ,  10.5,  10.0 ,   9.5,   9.0 ,
         8.5,   8.0 ,   7.5,   7.0 ,   6.5,   6.0 ,   5.5,   5.0 ,   4.5,
         4.0 ,   3.5,   3.0 ,   2.5,   2.0 ,   1.5,   1.0 ,   0.5,   0.0 ])

    upper_error = final_array1[0]+standard
    lower_error = final_array1[0]-standard

    sanity_cut = upper_error > 1.0
    sanity_index = np.arange(len(upper_error))
    input_index = sanity_index[sanity_cut]

    for i in range(len(input_index)):
        index_num = input_index[i]
        upper_error[index_num] = 1.0

    print upper_error

    axwidth = 3
    axlength = 10
    fontsize=24

    plt.rc('axes',linewidth=axwidth)  #make sure to do this before making any figures

    fig = plt.figure(figsize=(13,8))
    plt.subplot2grid((1,1),(0,0))
    plt.subplots_adjust(left=0.15, bottom=0.15, right=0.9, top=0.95, wspace=0.6, hspace=0.1)
    ax = plt.gca()
    minorLocatorx   = AutoMinorLocator(2)
    minorLocatory   = AutoMinorLocator(2)
    ax.xaxis.set_minor_locator(minorLocatorx)
    ax.yaxis.set_minor_locator(minorLocatory)
    #<do plotting>
    
    #plt.figure(figsize = (10,6), dpi = 100)
    plt.fill_between(plotbins, upper_error, lower_error, facecolor = 'grey', alpha = 0.4)
    p0 = plt.plot(plotbins, final_array1[0], color = 'k', linestyle = '-', linewidth = 5.5, label = r'$\rm GK14\ AM$')
    p1 = plt.plot(plotbins, final_array1a[0], color = 'Gold', linestyle ='-', linewidth = 4.0, label = r'$\rm GK14\ AM\ (<\ 0.5\ R_{vir})$')
    #plt.plot(plotbins, final_array1b[0], 'b-', linewidth = 1.5)
    p2 = plt.plot(plotbins, final_array2[0], color = 'm', linestyle = '-', linewidth = 4.0, label = r'$\rm GK14\ AM\ +\ scatter$')
    p3 = plt.plot(plotbins, final_array3[0], color = 'c', linestyle = '-', linewidth = 4.0, label = r'$\rm Behroozi\ AM$')
    #plt.plot(plotbins, final_array4[0], 'r-', linewidth = 1.5)

    plt.axhline(y = 0.90, xmin=0.0, xmax=0.109, linewidth=2.5, linestyle = '--', color='k')
    plt.axvline(x = 1.361, ymin=0.0, ymax=0.818, linewidth=2.5, linestyle = '--', color='k')

    pp0 = plt.Rectangle((0, 0), 1, 1, fc="k", ec = "k")
    pp1 = plt.Rectangle((0, 0), 1, 1, fc="Gold", ec = "Gold")
    pp2 = plt.Rectangle((0, 0), 1, 1, fc="m", ec = "m")
    pp3 = plt.Rectangle((0, 0), 1, 1, fc="c", ec="c")

    plt.legend(loc = 1, frameon=False, numpoints=1, prop={'size':16})
    
    #plt.legend([pp0,pp1,pp2,pp3], [r'$\rm GK14\ AM$', r'$\rm GK14\ AM\ (<\ 0.5\ R_{vir})$', r'$\rm GK14\ AM\ +\ scatter$',  r'$\rm Behroozi\ AM$'], loc = 1, frameon=False, numpoints=1, prop={'size':16})

    plt.xlabel('Lookback Time (Gyr)', fontsize = 24)
    plt.ylabel(r'$f\ (\rm t_{infall}>t_{lookback})$', fontsize = 24)
    plt.axis([0, 12.5, 0, 1.1])
    #plt.text(8.0,1.03, r'$\rm Garrison-Kimmel\ AM$', color = 'k')
    #plt.text(8.0,0.99, r'$\rm Behroozi\ AM$', color = 'c')
    #plt.text(8.0,0.95, r'$\rm Modified\ GK\ AM$', color = 'm')
    #plt.text(8.0,0.91, r'$\rm GK\ AM\ (<\ 0.5\ R_{vir})$', color = 'Gold')
    #plt.text(8.0,0.87, r'$\rm GK\ AM\ (V_{max}/V_{peak} >\ 0.7)$', color = 'b')
    #plt.text(8.0,0.83, r'$\rm High\ Mass$', color = 'r')

    plt.text(1.7,0.16, r'$\tau_{\rm quench}\ \sim\ 1.4\ \rm Gyr$', color = 'k', fontsize = 22)

    xtickloc = [0,2,4,6,8,10,12]
    xtickstr = ['$'+str(kk)+'$' for kk in xtickloc]
    ytickloc = [0.0,0.2,0.4,0.6,0.8,1.0]
    ytickstr = ['$'+str(kk)+'$' for kk in ytickloc]

    plt.xticks(xtickloc,xtickstr)
    plt.yticks(ytickloc,ytickstr)

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


    plt.savefig('distribution_of_infall_times_final.pdf')

    plt.show()
    

