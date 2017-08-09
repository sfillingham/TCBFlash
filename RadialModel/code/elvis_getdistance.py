# This will determine the distance from the host as a function of time for a list of subhalos in the ELVIS suite
#
#
# Notes:
#  distanceinput requires an ascii file that contains one column of ID numbers that correspond to subhalos in ELVIS at z = 0, the header should be labeled 'ID'
#  I typically use the '*alldwarfs*' file in the Input folder for whatever cuts I am looking to run.
#
#

import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as IUSpline
from astropy.table import Table

def work(distanceinput, velfrac, AM):

    #read input files
    inputfile = Table.read(distanceinput, format = 'ascii')

    #assign relevant columns from input files to numpy arrays
    inputlist = np.array(inputfile['ID'])
    hostlist = np.array(inputfile['ELVISname'])
    R2vir_test = np.array(inputfile['R2vir(kpc)'])

    #define output arrays and empty pathname list
    h1out = np.empty([len(inputlist),750])
    h2out = np.empty([len(inputlist),750])
    v1out = np.empty([len(inputlist),750])
    v2out = np.empty([len(inputlist),750])
    scale_out = np.empty([len(inputlist),750])
    con_out = np.empty([len(inputlist),750])
    mvir = np.empty([len(inputlist),750])
    Rvir1_out = np.empty([len(inputlist),750])
    Rvir2_out = np.empty([len(inputlist),750])
    loaded = []
    
    #loop through the subhalos
    for i in range(len(inputlist)):

        #select new subhalo ID
        ID = int(inputlist[i])
        
        print inputlist[i]
        print i

        #select host name 
        hostname = hostlist[i]

        IDpath = 'ELVIS_Data/Etracks/AllTrees/'+hostname+'/ID.txt'
        xpath = 'ELVIS_Data/Etracks/AllTrees/'+hostname+'/X.txt'
        ypath = 'ELVIS_Data/Etracks/AllTrees/'+hostname+'/Y.txt'
        zpath = 'ELVIS_Data/Etracks/AllTrees/'+hostname+'/Z.txt'
        vxpath = 'ELVIS_Data/Etracks/AllTrees/'+hostname+'/Vx.txt'
        vypath = 'ELVIS_Data/Etracks/AllTrees/'+hostname+'/Vy.txt'
        vzpath = 'ELVIS_Data/Etracks/AllTrees/'+hostname+'/Vz.txt'
        apath = 'ELVIS_Data/Etracks/AllTrees/'+hostname+'/scale.txt'
        virpath = 'ELVIS_Data/Etracks/AllTrees/'+hostname+'/Rvir.txt'
        rspath = 'ELVIS_Data/Etracks/AllTrees/'+hostname+'/Rs.txt'

        if (IDpath == loaded):

            #use the same files, which are already open
            pass 

        else:

            #open the new files 
            ID_input = np.loadtxt(IDpath)
            x_data = np.loadtxt(xpath)
            y_data = np.loadtxt(ypath)
            z_data = np.loadtxt(zpath)
            a_input = np.loadtxt(apath)
            vir_input = np.loadtxt(virpath)

            loaded = IDpath

        #select subhalo cut based on input list
        ID_data = ID_input[:,0]
        subhalo_cut = (ID == ID_data)

        index = np.arange(len(ID_data))
        index_input = index[subhalo_cut]

        #select Rvir2 in order to determine if the subhalo belongs to a single or double host
        virial2_test = R2vir_test[i]
        
        x = x_data[index_input[0]]
        y = y_data[index_input[0]]
        z = z_data[index_input[0]]
        a = a_input[index_input[0]]
        
        #clean and sort the position data for interpolation
        
        subcut = a > 0.
        x_clean = x[subcut]
        y_clean = y[subcut]
        z_clean = z[subcut]
        a_clean = a[subcut]
        x_cleansort = x_clean[::-1]
        y_cleansort = y_clean[::-1]
        z_cleansort = z_clean[::-1]
        a_cleansort = a_clean[::-1]
        
        #determine how many steps are needed for the interpolation
        interval = 750./len(a)
        
        
        #expand data by interpolating
        afinal = np.linspace(min(a_clean),1.,len(a_clean)*interval)

        xinterp = IUSpline(a_cleansort, x_cleansort, k=3)
        yinterp = IUSpline(a_cleansort, y_cleansort, k=3)
        zinterp = IUSpline(a_cleansort, z_cleansort, k=3)

        xpos = xinterp(afinal)
        ypos = yinterp(afinal)
        zpos = zinterp(afinal)
        
        #determine whether the simulation has 1 or 2 hosts and determine distances accordingly
        if virial2_test != -1: #two hosts
            
            print '2 hosts'

            Rvir1 = vir_input[0]
            vir1cut = Rvir1 > 0
            Rvir1_clean = Rvir1[vir1cut]
            Rvir1_cleansort = Rvir1_clean[::-1]
            a1 = a_input[0]
            a1_clean = a1[vir1cut]
            a1_cleansort = a1_clean[::-1]
            a1final = np.linspace(min(a1_clean),1.,len(a1_clean)*interval)
            
            Rvir2 = vir_input[1]
            vir2cut = Rvir2 > 0
            Rvir2_clean = Rvir2[vir2cut]
            Rvir2_cleansort = Rvir2_clean[::-1]
            a2 = a_input[1]
            a2_clean = a2[vir2cut]
            a2_cleansort = a2_clean[::-1]
            a2final = np.linspace(min(a2_clean),1.,len(a2_clean)*interval)
            
            Rvir1_interp = IUSpline(a1_cleansort, Rvir1_cleansort, k=3)
            Rvir2_interp = IUSpline(a2_cleansort, Rvir2_cleansort, k=3)
            
            virial1 = Rvir1_interp(a1final)
            virial2 = Rvir2_interp(a2final)
            
            print 'Host 1 x position'
            xHost1 = x_data[0]
            xHost1_clean = xHost1[subcut]
            xHost1_cleansort = xHost1_clean[::-1]
            xHost1interp = IUSpline(a_cleansort, xHost1_cleansort, k=3)
            xHost1pos = xHost1interp(afinal)

            print 'Host 2 x position'
            xHost2 = x_data[1]
            xHost2_clean = xHost2[subcut]
            xHost2_cleansort = xHost2_clean[::-1]
            xHost2interp = IUSpline(a_cleansort, xHost2_cleansort, k=3)
            xHost2pos = xHost2interp(afinal)
            
            print 'Host 1 y position'
            yHost1 = y_data[0]
            yHost1_clean = yHost1[subcut]
            yHost1_cleansort = yHost1_clean[::-1]
            yHost1interp = IUSpline(a_cleansort, yHost1_cleansort, k=3)
            yHost1pos = yHost1interp(afinal)
            
            print 'Host 2 y position'
            yHost2 = y_data[1]
            yHost2_clean = yHost2[subcut]
            yHost2_cleansort = yHost2_clean[::-1]
            yHost2interp = IUSpline(a_cleansort, yHost2_cleansort, k=3)
            yHost2pos = yHost2interp(afinal)
            
            print 'Host 1 z position'
            zHost1 = z_data[0]
            zHost1_clean = zHost1[subcut]
            zHost1_cleansort = zHost1_clean[::-1]
            zHost1interp = IUSpline(a_cleansort, zHost1_cleansort, k=3)
            zHost1pos = zHost1interp(afinal)
            
            print 'Host 2 z position'
            zHost2 = z_data[1]
            zHost2_clean = zHost2[subcut]
            zHost2_cleansort = zHost2_clean[::-1]
            zHost2interp = IUSpline(a_cleansort, zHost2_cleansort, k=3)
            zHost2pos = zHost2interp(afinal)
            
            
            dx1 = xHost1pos - xpos
            dy1 = yHost1pos - ypos
            dz1 = zHost1pos - zpos
            dx2 = xHost2pos - xpos
            dy2 = yHost2pos - ypos
            dz2 = zHost2pos - zpos
            
            comoving1 = np.sqrt(np.power(dx1,2) + np.power(dy1,2) + np.power(dz1,2))
            comoving2 = np.sqrt(np.power(dx2,2) + np.power(dy2,2) + np.power(dz2,2))
        
        
        else: #one host
            
            print '1 host'
            
            Rvir1 = vir_input[0]
            vir1cut = Rvir1 > 0
            Rvir1_clean = Rvir1[vir1cut]
            #print vir1cut
            #print Rvir1_clean
            Rvir1_cleansort = Rvir1_clean[::-1]
            a1 = a_input[0]
            a1_clean = a1[vir1cut]
            a1_cleansort = a1_clean[::-1]
            a1final = np.linspace(min(a1_clean),1.,len(a1_clean)*interval)
            
            Rvir1_interp = IUSpline(a1_cleansort, Rvir1_cleansort, k=3)

            virial1 = Rvir1_interp(a1final)
            virial2 = np.arange(len(virial1))
            virial2[:] = -1
            
            print 'Host 1 x position'
            xHost1 = x_data[0]
            xHost1_clean = xHost1[subcut]
            xHost1_cleansort = xHost1_clean[::-1]
            xHost1interp = IUSpline(a_cleansort, xHost1_cleansort, k=3)
            xHost1pos = xHost1interp(afinal)
            
            print 'Host 1 y position'
            yHost1 = y_data[0]
            yHost1_clean = yHost1[subcut]
            yHost1_cleansort = yHost1_clean[::-1]
            yHost1interp = IUSpline(a_cleansort, yHost1_cleansort, k=3)
            yHost1pos = yHost1interp(afinal)
            
            print 'Host 1 z position'
            zHost1 = z_data[0]
            zHost1_clean = zHost1[subcut]
            zHost1_cleansort = zHost1_clean[::-1]
            zHost1interp = IUSpline(a_cleansort, zHost1_cleansort, k=3)
            zHost1pos = zHost1interp(afinal)
            
            dx1 = xHost1pos - xpos
            dy1 = yHost1pos - ypos
            dz1 = zHost1pos - zpos
            
            comoving1 = np.sqrt(np.power(dx1,2) + np.power(dy1,2) + np.power(dz1,2))
            comoving2 = np.arange(len(xpos))
            comoving2[:] = -1

        
        #turn comoving into physical distances (Mpc)
        physical1 = comoving1*afinal
        physical2 = comoving2*afinal
        print len(physical1)
        print len(physical2)
        print len(virial1)
        print len(virial2)
        
        #determine the additional length necessary to get the distance array length to 750. Fill addarray with 0
        addon = 750 - len(comoving1)
        print addon
        addarray = np.arange(addon)
        addarray[:] = 0.0
        
        #append addarray onto respective distance arrays to get the length up to 750
        physical1 = np.append(addarray, physical1)
        physical2 = np.append(addarray, physical2)
        afinal = np.append(addarray, afinal)
    
        #determine the additional length necessary to get the virial radius array length to 750. Fill virialaddarray with 0
        virialaddon1 = 750 - len(virial1)
        virialaddon2 = 750 - len(virial2)
        virialaddarray1 = np.arange(virialaddon1)
        virialaddarray1[:] = 0.0
        virialaddarray2 = np.arange(virialaddon2)
        virialaddarray2[:] = 0.0
        
        #append virialaddarray onto respective virial arrays to get the length up to 750
        virial1 = np.append(virialaddarray1, virial1)
        virial2 = np.append(virialaddarray2, virial2)
        
        #put the subhalo arrays into the final output arrays
        print 'output'
        h1out[i] = physical1
        h2out[i] = physical2
        scale_out[i] = afinal
        Rvir1_out[i] = virial1
        Rvir2_out[i] = virial2
        

    #final output file pathnames
    outputfile1 = 'RadialModel/elvis_alldwarfs_v'+velfrac+'_distance_time_interp_host1_'+AM+'.dat'
    outputfile2 = 'RadialModel/elvis_alldwarfs_v'+velfrac+'_distance_time_interp_host2_'+AM+'.dat'
    outputfile_scale = 'RadialModel/elvis_alldwarfs_v'+velfrac+'_scalefactor_interp_'+AM+'.dat'
    outputfile_vir1 = 'RadialModel/elvis_alldwarfs_v'+velfrac+'_virial1_interp_'+AM+'.dat'
    outputfile_vir2 = 'RadialModel/elvis_alldwarfs_v'+velfrac+'_virial2_interp_'+AM+'.dat'

    #save the output arrays to the pathnames listed above
    np.savetxt(outputfile1, h1out)
    np.savetxt(outputfile2, h2out)
    np.savetxt(outputfile_scale, scale_out)
    np.savetxt(outputfile_vir1, Rvir1_out)
    np.savetxt(outputfile_vir2, Rvir2_out) 
