# This will determine the quenching time, then return the ID, quenching time, zfirst_interacted, and which host it belongs to, for the subhalos in the ELVIS suite.  The input need to be generated using elvis_getdistance.py, along with the 'elvis_alldwarfs_clean.dat' file.
# Input distance needs to be in kpc
#
#
#

import numpy as np
from astropy.table import Table
import astropy.cosmology as cosmo

def time(dist, ratio, vel, AM):

    host1data = np.loadtxt('RadialModel/elvis_alldwarfs_v'+vel+'_distance_time_interp_host1_'+AM+'.dat')
    host2data = np.loadtxt('RadialModel/elvis_alldwarfs_v'+vel+'_distance_time_interp_host2_'+AM+'.dat')
    scale = np.loadtxt('RadialModel/elvis_alldwarfs_v'+vel+'_scalefactor_interp_'+AM+'.dat')
    vir1data = np.loadtxt('RadialModel/elvis_alldwarfs_v'+vel+'_virial1_interp_'+AM+'.dat')
    vir2data = np.loadtxt('RadialModel/elvis_alldwarfs_v'+vel+'_virial2_interp_'+AM+'.dat')
    subhalolist = Table.read('RadialModel/InputData/elvis_alldwarfs_Etracks_distanceinput_r1.0_v'+vel+'_'+AM+'.dat', format = 'ascii')
    totalsublist = Table.read('ELVIS_Data_Clean/elvis_alldwarfs_clean_'+AM+'.dat', format = 'ascii')
    
    subhaloID = np.array(subhalolist['ID'])
    hostlist = np.array(subhalolist['ELVISname'])
    subhalotest = np.array(totalsublist['ID'])
    #zin1 = np.array(totalsublist['h1firstz_interacted'])
    #zin2 = np.array(totalsublist['h2firstz_interacted'])
    
    starform = np.array([])
    quench = np.array([])
    quenchID = np.array([])
    zinfall = np.array([])
    hostnum = np.array([])

    for i in range(len(host1data)):
        

        #read in host virial radius data 
        hostname = hostlist[i]
        Rvir1 = vir1data[i]
        Rvir2 = vir2data[i]

        #select distance data for subhalo, multiplying by 1000 converts units from Mpc to kpc
        dist1_data = host1data[i]*1000 
        dist2_data = host2data[i]*1000
        a_data = scale[i]

        #eliminate data from before subhalo formed, ie a = 0
        acut = a_data > 0
        a = a_data[acut]
        dist1 = dist1_data[acut]
        dist2 = dist2_data[acut]
        vir1 = Rvir1[acut]
        vir2 = Rvir2[acut]
        

        #determine when the subhalo crosses the virial radius, this is the infall redshift
        redcut1 = (dist1 < vir1) & (vir1 > 0)
        redcut2 = (dist2 < vir2) & (vir2 > 0)

        z = (1-a)/a
        
        zinfall1_cut = z[redcut1]
        zinfall2_cut = z[redcut2]

        if ((len(zinfall1_cut) !=0) & (len(zinfall2_cut) !=0)):
            zinfall1 = np.amax(zinfall1_cut)
            zinfall2 = np.amax(zinfall2_cut)

        elif ((len(zinfall1_cut) ==0) & (len(zinfall2_cut) !=0)):
            zinfall1 = -1
            zinfall2 = np.amax(zinfall2_cut)

        elif ((len(zinfall1_cut) !=0) & (len(zinfall2_cut) ==0)):
            zinfall1 = np.amax(zinfall1_cut)
            zinfall2 = -1

        else:
            zinfall1 = -1
            zinfall2 = -1
    
    
        #determine if the subhalo is within the input distance
        test1 = (dist1 < dist) & (dist1 > 0)
        test2 = (dist2 < dist) & (dist2 > 0)
        
        if (len(dist1[test1]) == 0) & (len(dist2[test2]) == 0):
            star_num = subhaloID[i]
            starform = np.append(starform, star_num)

        elif (len(dist1[test1]) == 0) & (len(dist2[test2]) != 0):
            
            vir2_test = vir2[test2]*ratio
            dist2_test = dist2[test2]
            z2_test = z[test2]

            cut = dist2_test < vir2_test
            d2 = dist2_test[cut]

            if (len(d2) != 0):

                z2 = np.amax(z2_test[cut])

                if (zinfall2 > z2):

                    #choose the redshift when the subhalo crosses the distance cut
                    qt2 = cosmo.lookback_time(z2)
                    
                    quench = np.append(quench, qt2)
                    quenchID = np.append(quenchID,subhaloID[i])
                    zinfall = np.append(zinfall, zinfall2)
                    hostnum = np.append(hostnum, 'host2')

                else:

                    #choose the redshift when the subhalo crosses the virial radius
                    qt2 = cosmo.lookback_time(zinfall2)

                    quench = np.append(quench, qt2)
                    quenchID = np.append(quenchID,subhaloID[i])
                    zinfall = np.append(zinfall, zinfall2)
                    hostnum = np.append(hostnum, 'host2')
                    
            else:
                
                star_num = subhaloID[i]
                starform = np.append(starform, star_num)
            
        elif (len(dist1[test1]) != 0) & (len(dist2[test2]) == 0):  
            
            vir1_test = vir1[test1]*ratio
            
            dist1_test = dist1[test1]
            z1_test = z[test1]

            cut = dist1_test < vir1_test
            d1 = dist1_test[cut]

            if (len(d1) != 0):

                z1 = np.amax(z1_test[cut])

                if (zinfall1 > z1):

                    #choose the redshift when the subhalo crosses the distance cut
                    qt1 = cosmo.lookback_time(z1)
                    
                    quench = np.append(quench, qt1)
                    quenchID = np.append(quenchID,subhaloID[i])
                    zinfall = np.append(zinfall, zinfall1)
                    hostnum = np.append(hostnum, 'host1')

                else:

                    #choose the redshift when the subhalo crosses the virial radius
                    qt1 = cosmo.lookback_time(zinfall1)

                    quench = np.append(quench, qt1)
                    quenchID = np.append(quenchID,subhaloID[i])
                    zinfall = np.append(zinfall, zinfall1)
                    hostnum = np.append(hostnum, 'host1')
                    
            else:
                
                star_num = subhaloID[i]
                starform = np.append(starform, star_num)
            
        else:
            
            vir1_test = vir1[test1]*ratio
            dist1_test = dist1[test1]
            z1_test = z[test1]
            vir2_test = vir2[test2]*ratio
            dist2_test = dist2[test2]
            z2_test = z[test2]
            
            cut1 = dist1_test < vir1_test
            d1 = dist1_test[cut1]
            cut2 = dist2_test < vir2_test
            d2 = dist2_test[cut2]
            

            if ((len(d1) != 0) & (len(d2) != 0)):

                z1 = np.amax(z1_test[cut1])
                z2 = np.amax(z2_test[cut2])
                redlist = np.array([z1, z2])
                redmax = np.amax(redlist)

                if ((z1 == redmax) & (z1 <= zinfall1)):

                    qt1 = cosmo.lookback_time(z1)
                    
                    quench = np.append(quench, qt1)
                    quenchID = np.append(quenchID,subhaloID[i])
                    zinfall = np.append(zinfall, zinfall1)
                    hostnum = np.append(hostnum, 'host1')

                elif ((z1 == redmax) & (z1 > zinfall1)):
                 
                    qt1 = cosmo.lookback_time(zinfall1)

                    quench = np.append(quench, qt1)
                    quenchID = np.append(quenchID,subhaloID[i])
                    zinfall = np.append(zinfall, zinfall1)
                    hostnum = np.append(hostnum, 'host1')

                elif ((z2 == redmax) & (z2 <= zinfall2)):

                    qt2 = cosmo.lookback_time(z2)
                    
                    quench = np.append(quench, qt2)
                    quenchID = np.append(quenchID,subhaloID[i])
                    zinfall = np.append(zinfall, zinfall2)
                    hostnum = np.append(hostnum, 'host2')

                elif ((z2 == redmax) & (z2 > zinfall2)):

                    qt2 = cosmo.lookback_time(zinfall2)

                    quench = np.append(quench, qt2)
                    quenchID = np.append(quenchID,subhaloID[i])
                    zinfall = np.append(zinfall, zinfall2)
                    hostnum = np.append(hostnum, 'host2')

                else: 

                    star_num = subhaloID[i]
                    starform = np.append(starform, star_num)

            elif ((len(d1) == 0) & (len(d2) != 0)):

                z2 = np.amax(z2_test[cut2])

                if (zinfall2 > z2):

                    #choose the redshift when the subhalo crosses the distance cut
                    qt2 = cosmo.lookback_time(z2)
                    
                    quench = np.append(quench, qt2)
                    quenchID = np.append(quenchID,subhaloID[i])
                    zinfall = np.append(zinfall, zinfall2)
                    hostnum = np.append(hostnum, 'host2')

                else:

                    #choose the redshift when the subhalo crosses the virial radius
                    qt2 = cosmo.lookback_time(zinfall2)

                    quench = np.append(quench, qt2)
                    quenchID = np.append(quenchID,subhaloID[i])
                    zinfall = np.append(zinfall, zinfall2)
                    hostnum = np.append(hostnum, 'host2')
                

            elif ((len(d1) != 0) & (len(d2) == 0)):

                z1 = np.max(z1_test[cut1])

                if (zinfall1 > z1):

                    #choose the redshift when the subhalo crosses the distance cut
                    qt1 = cosmo.lookback_time(z1)
                    
                    quench = np.append(quench, qt1)
                    quenchID = np.append(quenchID,subhaloID[i])
                    zinfall = np.append(zinfall, zinfall1)
                    hostnum = np.append(hostnum, 'host1')

                else:

                    
                    #choose the redshift when the subhalo crosses the virial radius
                    qt1 = cosmo.lookback_time(zinfall1)

                    quench = np.append(quench, qt1)
                    quenchID = np.append(quenchID,subhaloID[i])
                    zinfall = np.append(zinfall, zinfall1)
                    hostnum = np.append(hostnum, 'host1')
                
            else:
                
                star_num = subhaloID[i]
                starform = np.append(starform, star_num)
                    
    total = float(len(host1data))

    quenchfrac = len(quench) / total
    
    print len(quench)
    print len(starform)
    print total
    print quenchfrac
    
    dist_string = str(dist)
    ratio_string = str(ratio)
    outputfile = 'RadialModel/elvis_alldwarfs_quenchtime_r'+dist_string+'_v'+vel+'_'+ratio_string+'virial_interp_'+AM+'.dat'

    mastertable = Table([quenchID, quench, zinfall, hostnum], names=('ID', 'QuenchTime', 'zInfall', 'Host'), meta={'Quenchtime': 'ELVIS'})    
    mastertable.write(outputfile, format = 'ascii')
