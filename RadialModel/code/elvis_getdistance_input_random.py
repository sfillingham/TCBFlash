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

def work(v_ratio, AM):
    
    #read input files and select necessary columns to begin loop below
    inputlist = Table.read('elvis_inputlist.dat', format = 'ascii')
    inputfile = Table.read('ELVIS_Data_Clean/randomtest/elvis_alldwarfs_clean_ver2rand.dat', format = 'ascii')

    inputdata = np.array(inputlist['Galname'])
    R1vir = np.array(inputlist['R1vir(kpc)'])
    R2vir = np.array(inputlist['R2vir(kpc)'])
    elvisname = np.array(inputlist['Newname'])
    
    subID = np.array(inputfile['ID'])
    subVmax = np.array(inputfile['Vmax(km/s)'])
    subVpeak = np.array(inputfile['Vpeak(km/s)'])
    subMpeak = np.array(inputfile['Mpeak(Msun)'])
    subh1dist = np.array(inputfile['h1dist(Mpc)'])
    subh2dist = np.array(inputfile['h2dist(Mpc)'])
    subtag = np.array(inputfile['HostName'])
    subR1 = np.array(inputfile['R1vir(kpc)'])
    subR2 = np.array(inputfile['R2vir(kpc)'])
    
    
    v_string = str(v_ratio)
    
    #Eliminate dwarfs which do not meet the velocity ratio cuts
    vratio = subVmax / subVpeak
    
    cut = vratio > v_ratio

    #define the final arrays that will contain entire selected dwarf population
    masterID = subID[cut]
    masterVmax = subVmax[cut]
    masterVpeak = subVpeak[cut]
    masterMpeak = subMpeak[cut]
    masterh1dist = subh1dist[cut]
    masterh2dist = subh2dist[cut]
    mastertag = subtag[cut]
    masterR1 = subR1[cut]
    masterR2 = subR2[cut]
    masterELVIS = np.array([])


    #loop through the subhalos attaching the elvis host halo name

    for t in range(len(mastertag)):
        name = mastertag[t]
        namecut = (inputdata == name)
        etag = elvisname[namecut]
        #etag = tag[0]
        #print tag
        print etag
        masterELVIS = np.append(masterELVIS, etag)
        


    print masterELVIS
    #put master columns into table
    mastertable = Table([masterID, masterVmax, masterVpeak, masterMpeak, masterh1dist, masterh2dist, mastertag, masterR1, masterR2, masterELVIS], names=('ID', 'Vmax(km/s)', 'Vpeak(km/s)', 'Mpeak(Msun)', 'h1dist(Mpc)', 'h2dist(Mpc)', 'HostName', 'R1vir(kpc)', 'R2vir(kpc)', 'ELVISname'), meta={'Mass Cuts': 'ELVIS'})

    masteroutputfile = 'RadialModel/InputData/elvis_alldwarfs_Etracks_distanceinput_r1.0_v'+v_string+'_'+AM+'.dat'
    mastertable.write(masteroutputfile, format = 'ascii')






