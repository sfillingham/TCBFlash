#This file will get the quenching timescale using the ELVIS simulation.
#Input files need to be generated for the subhalos, quenched fraction, and fractional distance from host
#Subhalo file needs to be cleaned using mass cuts based on a specific abundance matching routine.

import numpy as np
import astropy.cosmology as cosmo
from astropy.table import Table


def quenchtime_all(subhalo_input, red_input, fq_input, r_input, v_ratio):

    #read in files and select relevant columns 
    subhalo_data = Table.read(subhalo_input, format = 'ascii')
    red_data = Table.read(red_input, format = 'ascii')
    fq_data = Table.read(fq_input, format = 'ascii')
    r_data = Table.read(r_input, format = 'ascii')

    h1dist = np.array(subhalo_data['h1dist(Mpc)'])
    h2dist = np.array(subhalo_data['h2dist(Mpc)'])
    h1firstz = np.array(red_data['h1z_infall'])
    h2firstz = np.array(red_data['h2z_infall'])
    #h1firstz = np.array(subhalo_data['h1firstz_interacted'])
    #h2firstz = np.array(subhalo_data['h2firstz_interacted'])
    vpeak = np.array(subhalo_data['Vpeak(km/s)'])
    vmax = np.array(subhalo_data['Vmax(km/s)'])
    hostname = np.array(subhalo_data['HostName'])
    R1vir = np.array(subhalo_data['R1vir(kpc)'])
    R2vir = np.array(subhalo_data['R2vir(kpc)'])
    fq = np.array(fq_data['fq'])
    r_ratio = np.array(r_data['radius'])

    r_string = np.array(map(str, r_ratio))
    v_string = str(v_ratio)

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
    
    #select velocity ratio cuts
    vel_cut = master_velfrac > v_ratio
           
    #loop through fractional distance
    for j in range(len(r_ratio)):
        
        #select distance cuts
        dist_ratio = r_ratio[j]
        dist_cut = master_distfrac < dist_ratio

        #define quenching time array
        master_tq = np.array([])

        #loop through quenched fraction
        for k in range(len(fq)):

            quenched_frac = fq[k]

            #eliminate dwarfs that do not satisfy necessary conditions
            cut = vel_cut & dist_cut & (master_zinfall != -1)

            #determine the sorted list of infall times 
            t = cosmo.lookback_time(master_zinfall[cut])
            t_infall = np.sort(t)

            #given quenched fraction, determine the quenching time
            t_slice = np.round(len(t_infall)*quenched_frac)
            search_slice = len(t_infall) - t_slice

            #t_quench = (t_infall[search_slice] + t_infall[search_slice - 1]) / 2
            t_quench = t_infall[search_slice]
            
            master_tq = np.append(master_tq, t_quench)

        output = Table([fq, master_tq], names=('fq', 'tq'), meta={'fq vs tq': 'ELVIS'})

        outputfile = 'TimeModel/mass68/elvis_fqtq_r'+r_string[j]+'_v'+v_string+'_m68_alt.dat'
        output.write(outputfile, format = 'ascii')
            


def quenchtime_individual(input_list, red_input, fq_input, r_input, v_ratio):

    #loop through each input file in order to determine the 
    #quenching timescale for each host.

    inputlist = Table.read(input_list, format = 'ascii')
    fq_data = Table.read(fq_input, format = 'ascii')
    r_data = Table.read(r_input, format = 'ascii')
    red_data = Table.read(red_input, format = 'ascii')
    
    galname = np.array(inputlist['Galname'])
    fq = np.array(fq_data['fq'])
    r_ratio = np.array(r_data['radius'])
    redID = np.array(red_data['ID'])
    red1 = np.array(red_data['h1z_infall'])
    red2 = np.array(red_data['h2z_infall'])
    
    r_string = np.array(map(str, r_ratio))
    v_string = str(v_ratio)

    for s in range(len(galname)):
        
        print galname[s]
        
        #read in files and select relevant columns 
        subhalo_data = Table.read('ELVIS_Data_Clean/elvis_'+galname[s]+'_clean_m58.dat', format = 'ascii')

        subID = np.array(subhalo_data['ID'])
        h1dist = np.array(subhalo_data['h1dist(Mpc)'])
        h2dist = np.array(subhalo_data['h2dist(Mpc)'])
        #h1firstz = np.array(subhalo_data['h1firstz_interacted'])
        #h2firstz = np.array(subhalo_data['h2firstz_interacted'])
        vpeak = np.array(subhalo_data['Vpeak(km/s)'])
        vmax = np.array(subhalo_data['Vmax(km/s)'])
        hostname = np.array(subhalo_data['HostName'])
        R1vir = np.array(subhalo_data['R1vir(kpc)'])
        R2vir = np.array(subhalo_data['R2vir(kpc)'])
        

        #define empty velocity and distance columns
        master_distfrac1 = np.array([])
        master_distfrac2 = np.array([])
        master_velfrac1 = np.array([])
        master_velfrac2 = np.array([])
        master_zinfall1 = np.array([])
        master_zinfall2 = np.array([])

        #determine whether isolated or double
        if R2vir[0] > 0:
        
            #create fractional velocity, fractional distance, and infall timescale columns for each host in the system
            for i in range(len(h1dist)):
                #print 100+i

                redcut = subID[i] == redID
                h1firstz = red1[redcut]
                h2firstz = red2[redcut]
                
                dist1 = h1dist[i]*1000 #multiplying by 1000 changes units from Mpc to kpc
                dist2 = h2dist[i]*1000 #multiplying by 1000 changes units from Mpc to kpc

                cond_a = (h1firstz[0] > h2firstz[0]) & (dist1 < dist2) & (dist1 < R1vir[i]) & (dist1 != -1000)
                cond_b = (h1firstz[0] > h2firstz[0]) & (dist2 < dist1) & (dist2 < R2vir[i]) & (dist2 != -1000)
                cond_c = (h2firstz[0] > h1firstz[0]) & (dist2 < dist1) & (dist2 < R2vir[i]) & (dist2 != -1000)
                cond_d = (h2firstz[0] > h1firstz[0]) & (dist1 < dist2) & (dist1 < R1vir[i]) & (dist1 != -1000)
            
                if (cond_a or cond_d):
                    distance = dist1
                    z_infall = h1firstz[0]
                    virial = R1vir[i]

                    #append data for host1
                    dist_frac = distance/virial
                    vel_frac = vmax[i]/vpeak[i]

                    master_distfrac1 = np.append(master_distfrac1, dist_frac)
                    master_velfrac1 = np.append(master_velfrac1, vel_frac)
                    master_zinfall1 = np.append(master_zinfall1, z_infall)

                elif (cond_b or cond_c):
                    distance = dist2
                    z_infall = h2firstz[0]
                    virial = R2vir[i]

                    #append data for host2
                    dist_frac = distance/virial
                    vel_frac = vmax[i]/vpeak[i]

                    master_distfrac2 = np.append(master_distfrac2, dist_frac)
                    master_velfrac2 = np.append(master_velfrac2, vel_frac)
                    master_zinfall2 = np.append(master_zinfall2, z_infall)

                else:
                    continue

            vel_cut1 = master_velfrac1 > v_ratio
            vel_cut2 = master_velfrac2 > v_ratio
           
            #loop through fractional distance
            for j in range(len(r_ratio)):
        
                #select distance cuts
                dist_ratio = r_ratio[j]
                #print dist_ratio
                
                dist_cut1 = (master_distfrac1 < dist_ratio) & (master_distfrac1 > 0)
                dist_cut2 = (master_distfrac2 < dist_ratio) & (master_distfrac2 > 0)

                #define quenching time array
                master_tq1 = np.array([])
                master_tq2 = np.array([])
                master_time1 = np.array([])
                master_time2 = np.array([])
                masterz1 = np.array([])
                masterz2 = np.array([])
            
                #loop through quenched fraction
                for k in range(len(fq)):

                    quenched_frac = fq[k]

                    print quenched_frac

                    #eliminate dwarfs that do not satisfy necessary conditions
                    cut1 = vel_cut1 & dist_cut1 & (master_zinfall1 != -1)
                    cut2 = vel_cut2 & dist_cut2 & (master_zinfall2 != -1)

                    #determine the sorted list of infall times 
                    t1 = cosmo.lookback_time(master_zinfall1[cut1])
                    t2 = cosmo.lookback_time(master_zinfall2[cut2])

                    #print len(t1)
                    #print len(t2)
                    
                    t_infall1 = np.sort(t1)
                    t_infall2 = np.sort(t2)
                    
                    #given quenched fraction, determine the quenching time
                    t_slice1 = np.round(len(t_infall1)*quenched_frac)
                    search_slice1 = len(t_infall1) - t_slice1
                    t_slice2 = np.round(len(t_infall2)*quenched_frac)
                    search_slice2 = len(t_infall2) - t_slice2

                    #t_quench1 = (t_infall1[search_slice1] + t_infall1[search_slice1 - 1]) / 2
                    t_quench1 = t_infall1[search_slice1]
                    master_tq1 = np.append(master_tq1, t_quench1)
                    #t_quench2 = (t_infall2[search_slice2] + t_infall2[search_slice2 - 1]) / 2
                    t_quench2 = t_infall2[search_slice2]
                    master_tq2 = np.append(master_tq2, t_quench2)

                master_time1 = np.append(master_time1, t1)
                master_time2 = np.append(master_time2, t2)
                masterz1 = np.append(masterz1, master_zinfall1[cut1])
                masterz2 = np.append(masterz2, master_zinfall2[cut2])

                #print master_time1
                #print masterz1
                #print master_time2
                #print masterz2

                output1 = Table([fq, master_tq1], names=('fq', 'tq'), meta={'fq vs tq': 'ELVIS'})
                output2 = Table([fq, master_tq2], names=('fq', 'tq'), meta={'fq vs tq': 'ELVIS'})

                timeoutput1 = Table([masterz1, master_time1], names=('Infallz', 'InfallTime'), meta={'ELVIS' : 'Infall Time'})
                timeoutput2 = Table([masterz2, master_time2], names=('Infallz', 'InfallTime'), meta={'ELVIS' : 'Infall Time'})

                timeoutputfile1 = 'TimeModel/mass58/elvis_'+galname[s]+'_host1_infalltime_r'+r_string[j]+'_v'+v_string+'_m58.dat'
                timeoutput1.write(timeoutputfile1, format = 'ascii')
                outputfile1 = 'TimeModel/mass58/elvis_'+galname[s]+'_host1_fqtq_r'+r_string[j]+'_v'+v_string+'_m58.dat'
                output1.write(outputfile1, format = 'ascii')
                
                timeoutputfile2 = 'TimeModel/mass58/elvis_'+galname[s]+'_host2_infalltime_r'+r_string[j]+'_v'+v_string+'_m58.dat'
                timeoutput2.write(timeoutputfile2, format = 'ascii')
                outputfile2 = 'TimeModel/mass58/elvis_'+galname[s]+'_host2_fqtq_r'+r_string[j]+'_v'+v_string+'_m58.dat'
                output2.write(outputfile2, format = 'ascii')

        else : # R2vir[i] < 0

            #create fractional velocity, fractional distance, and infall timescale columns for each host in the system
            for i in range(len(h1dist)):

                redcut = subID[i] == redID
                h1firstz = red1[redcut]
                h2firstz = red2[redcut]

                dist1 = h1dist[i]*1000 #multiplying by 1000 changes units from Mpc to kpc

                distance = dist1
                z_infall = h1firstz[0]
                virial = R1vir[i]

                dist_frac = distance/virial
                vel_frac = vmax[i]/vpeak[i]

                master_distfrac1 = np.append(master_distfrac1, dist_frac)
                master_velfrac1 = np.append(master_velfrac1, vel_frac)
                master_zinfall1 = np.append(master_zinfall1, z_infall)

            vel_cut = master_velfrac1 > v_ratio
           
            #loop through fractional distance
            for j in range(len(r_ratio)):
        
                #select distance cuts
                dist_ratio = r_ratio[j]
                dist_cut = (master_distfrac1 < dist_ratio) & (master_distfrac1 > 0)

                #define quenching time array
                master_tq = np.array([])
                master_time = np.array([])
                masterz = np.array([])

                #loop through quenched fraction
                for k in range(len(fq)):

                    quenched_frac = fq[k]

                    #eliminate dwarfs that do not satisfy necessary conditions
                    cut = vel_cut & dist_cut & (master_zinfall1 != -1)
                
                    #determine the sorted list of infall times 
                    t = cosmo.lookback_time(master_zinfall1[cut])
                    t_infall = np.sort(t)
                    #print len(t_infall)
                    #given quenched fraction, determine the quenching time
                    t_slice = np.round(len(t_infall)*quenched_frac)
                    search_slice = len(t_infall) - t_slice
                    #print t_infall
                    #t_quench = (t_infall[search_slice] + t_infall[search_slice - 1]) / 2
                    t_quench = t_infall[search_slice]
                    #f_quench = (len(t_infall)-1)/np.float(len(t_infall))
            
                    master_tq = np.append(master_tq, t_quench)

                master_time = np.append(master_time, t)
                masterz = np.append(masterz, master_zinfall1[cut])

                output = Table([fq, master_tq], names=('fq', 'tq'), meta={'fq vs tq': 'ELVIS'})

                outputfile = 'TimeModel/mass58/elvis_'+galname[s]+'_fqtq_r'+r_string[j]+'_v'+v_string+'_m58.dat'
                output.write(outputfile, format = 'ascii')

                timeoutput1 = Table([masterz, master_time], names=('Infallz', 'InfallTime'), meta={'ELVIS' : 'Infall Time'})

                timeoutputfile1 = 'TimeModel/mass58/elvis_'+galname[s]+'_infalltime_r'+r_string[j]+'_v'+v_string+'_m58.dat'
                timeoutput1.write(timeoutputfile1, format = 'ascii')
                
