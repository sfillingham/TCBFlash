#This will select subhalos from the ELVIS suite based on the input criteria
#Input Criteria:
#      mass min, mass max, r_ratio, v_ratio
#  
#Such that: mass min < log(Mstar) < mass max; distance(z=0) < (r_ratio)*(Rvir); (Vmax / Vpeak) > v_ratio
#
#The input files are taken from the email sent by SGK to MC, UCI Summer 2014
#
#Be sure to input the correct abundance matching prescription in the work function
#Notes:
#  ver1 = GK14AM and inside 1.0 Rvir
#  ver1a = GK14AM and inside 0.5 Rvir
#  ver2 = RandomAM and inside 1.0 Rvir
#  ver3 = Extension down to lower mass, Mstar ~ 1e5 - 1e6
#  ver4 = GK14AM, High Mass, and inside 1.0 Rvir
#  BehrooziAM = Behroozi abundance matching in the same mass range as version 1
#  I am redo-ing Behroozi in order to ensure the analysis is correct.
#

import numpy as np
from astropy.table import Table

def work(Mstar_min, Mstar_max, r_ratio, v_ratio, AM):

    #check the input to ensure its valid
    if Mstar_min >= Mstar_max:
        print 'Fail! Try again...Min first, Max second'

    else:
    
        #read input files and select necessary columns to begin loop below
        inputlist = Table.read('elvis_inputlist.dat', format = 'ascii')
        
        #abundance matching prescription
        ab_match = Table.read('GK14AM.txt', format = 'ascii')

        inputdata = np.array(inputlist['Galname'])
        R1vir = np.array(inputlist['R1vir(kpc)'])
        R2vir = np.array(inputlist['R2vir(kpc)'])
        elvisname = np.array(inputlist['Newname'])
        Mpeak = np.array(ab_match['Mpeak(Msun)'])
        Mstar = np.array(ab_match['Mstar(Msun)'])

        #determine the range of halo masses based on interpolation of abundance matching file
        Mpeak_min = np.interp(Mstar_min, Mstar, Mpeak)
        Mpeak_max = np.interp(Mstar_max, Mstar, Mpeak)

        #define the final arrays that will contain entire selected dwarf population
        masterID = np.array([])
        masterVmax = np.array([])
        masterVpeak = np.array([])
        masterMpeak = np.array([])
        masterh1dist = np.array([])
        masterh2dist = np.array([])
        mastertag = np.array([])
        masterELVIS = np.array([])
        masterR1 = np.array([])
        masterR2 = np.array([])

        r_string = str(r_ratio)
        v_string = str(v_ratio)

        #loop through all the files in the ELVIS suite 

        for t in range(len(inputdata)):

            data = Table.read('ELVIS_Data/Subhalos/'+inputdata[t]+'.txt', format = 'ascii')

            #Select relevant columns from the files
            ID = np.array(data['(0)ID'])
            Vmax = np.array(data['(2)Vmax(km/s)'])
            Vpeak = np.array(data['(4)Vpeak(km/s)'])
            Mpeak = np.array(data['(5)Mpeak(Msun)'])
            h1dist = np.array(data['(6)h1dist(Mpc)'])
            h2dist = np.array(data['(9)h2dist(Mpc)'])

            #define empty velocity and distance columns
            master_distfrac = np.array([])
            master_velfrac = np.array([])
            
        
            #create fractional velocity and fractional distance columns in order to make cuts
            for i in range(len(h1dist)):
        
                dist1 = h1dist[i]*1000 #multiplying by 1000 changes units from Mpc to kpc
                dist2 = h2dist[i]*1000 #multiplying by 1000 changes units from Mpc to kpc
        
                cond_a = (dist1 < dist2) & (dist1 < R1vir[t]) & (dist1 != -1000)
                cond_b = (dist2 < dist1) & (dist2 < R2vir[t]) & (dist2 != -1000)
                cond_c = (dist1 < dist2) & (dist1 < R1vir[t]) & (dist1 == -1000)
                cond_d = (dist2 < dist1) & (dist2 < R2vir[t]) & (dist2 == -1000)
            
                if (cond_a or cond_d):
                    distance = dist1
                    virial = R1vir[t]
                    dist_frac = distance/virial

                elif (cond_b or cond_c):
                    distance = dist2
                    virial = R2vir[t]
                    dist_frac = distance/virial

                else:
                    dist_frac = -1

                vel_frac = Vmax[i]/Vpeak[i]

                master_distfrac = np.append(master_distfrac, dist_frac)
                master_velfrac = np.append(master_velfrac, vel_frac)

            #select velocity ratio cut
            vel_cut = master_velfrac > v_ratio
            print len(vel_cut)

            #select distance ratio cut
            dist_cut = (master_distfrac < r_ratio) & (master_distfrac != -1)
            print len(dist_cut)

            #create the mass cuts based on the interpolation above
            masscuts = (Mpeak > Mpeak_min) & (Mpeak < Mpeak_max)
            print len(masscuts)

            #Total cuts
            cuts = vel_cut & dist_cut & masscuts

            #select galaxies based on the above cuts
            ID_sel = ID[cuts]
            Vmax_sel = Vmax[cuts]
            Vpeak_sel = Vpeak[cuts]
            Mpeak_sel = Mpeak[cuts]
            h1dist_sel = h1dist[cuts]
            h2dist_sel = h2dist[cuts]

            #tag dwarfs with host name and respective virial radii
            tag = np.chararray(len(ID_sel))
            tag = np.chararray(tag.shape, itemsize = 20)
            tag[:] = inputdata[t]
            etag = np.chararray(len(ID_sel))
            etag = np.chararray(tag.shape, itemsize = 20)
            etag[:] = elvisname[t]
            radius1 = np.arange(len(ID_sel))
            radius2 = np.arange(len(ID_sel))
            radius1[:] = R1vir[t]
            radius2[:] = R2vir[t]

            #append each column into master columns 
            masterID = np.append(masterID, ID_sel)
            masterVmax = np.append(masterVmax, Vmax_sel)
            masterVpeak = np.append(masterVpeak, Vpeak_sel)
            masterMpeak = np.append(masterMpeak, Mpeak_sel)
            masterh1dist = np.append(masterh1dist, h1dist_sel)
            masterh2dist = np.append(masterh2dist, h2dist_sel)
            mastertag = np.append(mastertag, tag)
            masterR1 = np.append(masterR1, radius1)
            masterR2 = np.append(masterR2, radius2)
            masterELVIS = np.append(masterELVIS, etag)
            
            #put individual columns into table
            output = Table([ID_sel, Vmax_sel, Vpeak_sel, Mpeak_sel, h1dist_sel, h2dist_sel, tag, radius1, radius2, etag], names=('ID', 'Vmax(km/s)', 'Vpeak(km/s)', 'Mpeak(Msun)', 'h1dist(Mpc)', 'h2dist(Mpc)', 'HostName', 'R1vir(kpc)', 'R2vir(kpc)', 'ELVISname'), meta={'Mass Cuts': 'ELVIS'})

            outputfile = 'RadialModel/InputData/elvis_'+inputdata[t]+'_Etracks_distanceinput_r'+r_string+'_v'+v_string+'_'+AM+'.dat'
            output.write(outputfile, format = 'ascii')

        print len(masterID)
        #put master columns into table
        mastertable = Table([masterID, masterVmax, masterVpeak, masterMpeak, masterh1dist, masterh2dist, mastertag, masterR1, masterR2, masterELVIS], names=('ID', 'Vmax(km/s)', 'Vpeak(km/s)', 'Mpeak(Msun)', 'h1dist(Mpc)', 'h2dist(Mpc)', 'HostName', 'R1vir(kpc)', 'R2vir(kpc)', 'ELVISname'), meta={'Mass Cuts': 'ELVIS'})    

        masteroutputfile = 'RadialModel/InputData/elvis_alldwarfs_Etracks_distanceinput_r'+r_string+'_v'+v_string+'_'+AM+'.dat'
        mastertable.write(masteroutputfile, format = 'ascii')
