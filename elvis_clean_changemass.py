#This code will take a file from the Elvis_Data directory and select the columns 
#needed for the Summer 2014 project

#Column Names in the ELVIS files:
#(0)ID	(1)interacted	(2)Vmax(km/s)	(3)Mvir(Msun)	(4)Vpeak(km/s)	(5)Mpeak(Msun)	(6)h1dist(Mpc)	(7)h1vradphysical(km/s)	(8)h1vtan(km/s)	(9)h2dist(Mpc)	(10)h2vradphysical(km/s)	(11)h2vtan(km/s)	(12)h1mindist(comovingMpc)	(13)h2mindist(comovingMpc)	(14)h1firstz_interacted	(15)h2firstz_interacted	(16)treenum	(17)normh1dist	(18)normh1Vrad	(19)normh1Vtan	(20)normh2dist	(21)normh2Vrad	(22)normh2Vtan	(23)Angleh1toh2(rads)	(24)Angleh2toh1(rads)	(25)h1lastz_interacted	(26)h2lastz_interacted	(27)x(Mpc)	(28)y(Mpc)	(29)z(Mpc)	(30)peculiarVx(km/s)	(31)peculiarVy(km/s)	(32)peculiarVz(km/s)

import numpy as np
from astropy.table import Table

#use inputlist file and all data files from ELVIS_Data directory

def work(Mhalo_min, Mhalo_max):

    #check the input to ensure its valid
    if Mhalo_min >= Mhalo_max:
        print 'Fail! Try again...'

    else:
    
        #read input files and select necessary columns to begin loop below
        inputlist = Table.read('elvis_inputlist.dat', format = 'ascii')
        ab_match = Table.read('GK14AM.txt', format = 'ascii')

        inputdata = np.array(inputlist['Galname'])
        R1vir = np.array(inputlist['R1vir(kpc)'])
        R2vir = np.array(inputlist['R2vir(kpc)'])
        Mpeak = np.array(ab_match['Mpeak(Msun)'])
        Mstar = np.array(ab_match['Mstar(Msun)'])

        #determine the range of halo masses based on interpolation of abundance matching file or explicit selection (see below)
        # Mpeak_min = np.interp(Mstar_min, Mstar, Mpeak)
        # Mpeak_max = np.interp(Mstar_max, Mstar, Mpeak)
        Mpeak_min = Mhalo_min
        Mpeak_max = Mhalo_max
        

        #define the final arrays that will contain entire selected dwarf population
        masterID = np.array([])
        masterVmax = np.array([])
        masterVpeak = np.array([])
        masterMpeak = np.array([])
        masterh1dist = np.array([])
        masterh2dist = np.array([])
        masterh1firstz = np.array([])
        masterh2firstz = np.array([])
        masterx = np.array([])
        mastery = np.array([])
        masterz = np.array([])
        mastertag = np.array([])
        masterR1 = np.array([])
        masterR2 = np.array([])
    
        #loop through all the files in the ELVIS suite 

        for i in range(len(inputdata)):

            data = Table.read('ELVIS_Data/Subhalos/'+inputdata[i]+'.txt', format = 'ascii')

            #Select relevant columns from the files
            ID = np.array(data['(0)ID'])
            Vmax = np.array(data['(2)Vmax(km/s)'])
            Vpeak = np.array(data['(4)Vpeak(km/s)'])
            Mpeak = np.array(data['(5)Mpeak(Msun)'])
            h1dist = np.array(data['(6)h1dist(Mpc)'])
            h2dist = np.array(data['(9)h2dist(Mpc)'])
            h1firstz_interact = np.array(data['(14)h1firstz_interacted'])
            h2firstz_interact = np.array(data['(15)h2firstz_interacted'])
            xpos = np.array(data['(27)x(Mpc)'])
            ypos = np.array(data['(28)y(Mpc)'])
            zpos = np.array(data['(29)z(Mpc)'])

            #create the mass cuts based on the interpolation above
            masscuts = (Mpeak > Mpeak_min) & (Mpeak < Mpeak_max)

            #select galaxies based on the above cuts
            ID_sel = ID[masscuts]
            Vmax_sel = Vmax[masscuts]
            Vpeak_sel = Vpeak[masscuts]
            Mpeak_sel = Mpeak[masscuts]
            h1dist_sel = h1dist[masscuts]
            h2dist_sel = h2dist[masscuts]
            h1firstz_sel = h1firstz_interact[masscuts]
            h2firstz_sel = h2firstz_interact[masscuts]
            xpos_sel = xpos[masscuts]
            ypos_sel = ypos[masscuts]
            zpos_sel = zpos[masscuts]

            #tag dwarfs with host name and respective virial radii
            tag = np.chararray(len(ID_sel))
            tag = np.chararray(tag.shape, itemsize = 20)
            tag[:] = inputdata[i]
            radius1 = np.arange(len(ID_sel))
            radius2 = np.arange(len(ID_sel))
            radius1[:] = R1vir[i]
            radius2[:] = R2vir[i]

            #append each column into master columns 
            masterID = np.append(masterID, ID_sel)
            masterVmax = np.append(masterVmax, Vmax_sel)
            masterVpeak = np.append(masterVpeak, Vpeak_sel)
            masterMpeak = np.append(masterMpeak, Mpeak_sel)
            masterh1dist = np.append(masterh1dist, h1dist_sel)
            masterh2dist = np.append(masterh2dist, h2dist_sel)
            masterh1firstz = np.append(masterh1firstz, h1firstz_sel)
            masterh2firstz = np.append(masterh2firstz, h2firstz_sel)
            masterx = np.append(masterx, xpos_sel)
            mastery = np.append(mastery, ypos_sel)
            masterz = np.append(masterz, zpos_sel)
            mastertag = np.append(mastertag, tag)
            masterR1 = np.append(masterR1, radius1)
            masterR2 = np.append(masterR2, radius2)
            
            #put individual columns into table
            output = Table([ID_sel, Vmax_sel, Vpeak_sel, Mpeak_sel, h1dist_sel, h2dist_sel, h1firstz_sel, h2firstz_sel, xpos_sel, ypos_sel, zpos_sel, tag, radius1, radius2], names=('ID', 'Vmax(km/s)', 'Vpeak(km/s)', 'Mpeak(Msun)', 'h1dist(Mpc)', 'h2dist(Mpc)', 'h1firstz_interacted', 'h2firstz_interacted', 'x(Mpc)', 'y(Mpc)', 'z(Mpc)', 'HostName', 'R1vir(kpc)', 'R2vir(kpc)'), meta={'Mass Cuts': 'ELVIS'})

            outputfile = 'ELVIS_Data_Clean/elvis_'+inputdata[i]+'_clean_ver2.dat'
            output.write(outputfile, format = 'ascii')

        #put master columns into table
        mastertable = Table([masterID, masterVmax, masterVpeak, masterMpeak, masterh1dist, masterh2dist, masterh1firstz, masterh2firstz, masterx, mastery, masterz, mastertag, masterR1, masterR2], names=('ID', 'Vmax(km/s)', 'Vpeak(km/s)', 'Mpeak(Msun)', 'h1dist(Mpc)', 'h2dist(Mpc)', 'h1firstz_interacted', 'h2firstz_interacted', 'x(Mpc)', 'y(Mpc)', 'z(Mpc)', 'HostName', 'R1vir(kpc)', 'R2vir(kpc)'), meta={'Mass Cuts': 'ELVIS'})    

        masteroutputfile = 'ELVIS_Data_Clean/elvis_alldwarfs_clean_ver2.dat'
        mastertable.write(masteroutputfile, format = 'ascii')
