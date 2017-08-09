#This script takes the elvis_clean_changemass.py output and selects 15 random subhalos for each host

import numpy as np
from astropy.table import Table

#number = number of subhalos to randomly select
def random(number):

    data1 = Table.read('elvis_inputlist.dat', format = 'ascii')
    inputdata = np.array(data1['Galname'])
    
    #define the final arrays that will contain entire selected dwarf population
    masterID = np.array([])
    masterVmax = np.array([])
    masterVpeak = np.array([])
    masterMpeak = np.array([])
    masterh1dist = np.array([])
    masterh2dist = np.array([])
    masterh1firstz = np.array([])
    masterh2firstz = np.array([])
    mastertag = np.array([])
    masterR1 = np.array([])
    masterR2 = np.array([])

    for i in range(len(inputdata)):

        data = Table.read('ELVIS_Data_clean/elvis_'+inputdata[i]+'_clean_ver2.dat', format = 'ascii')

        #read in the relevant columns of data
        ID = np.array(data['ID'])
        Vmax = np.array(data['Vmax(km/s)'])
        Vpeak = np.array(data['Vpeak(km/s)'])
        Mpeak = np.array(data['Mpeak(Msun)'])
        h1dist = np.array(data['h1dist(Mpc)'])
        h2dist = np.array(data['h2dist(Mpc)'])
        h1firstz_interact = np.array(data['h1firstz_interacted'])
        h2firstz_interact = np.array(data['h2firstz_interacted'])
        R1vir = np.array(data['R1vir(kpc)'])
        R2vir = np.array(data['R2vir(kpc)'])
        hostname = np.array(data['HostName'])

        #define a shitload of arrays that need to erase after each loop.
        ID_sel = np.array([])
        Vmax_sel = np.array([])
        Vpeak_sel = np.array([])
        Mpeak_sel = np.array([])
        h1dist_sel = np.array([])
        h2dist_sel = np.array([])
        h1firstz_sel = np.array([])
        h2firstz_sel = np.array([])
        radius1 = np.array([])
        radius2 = np.array([])
        tag = np.array([])
        
        ida = np.array([])
        vma = np.array([])
        vpa = np.array([])
        mpa = np.array([])
        h1dista = np.array([])
        h2dista = np.array([])
        h1za = np.array([])
        h2za = np.array([])
        R1a = np.array([])
        R2a = np.array([])
        hosta = np.array([])

        idb = np.array([])
        vmb = np.array([])
        vpb = np.array([])
        mpb = np.array([])
        h1distb = np.array([])
        h2distb = np.array([])
        h1zb = np.array([])
        h2zb = np.array([])
        R1b = np.array([])
        R2b = np.array([])
        hostb = np.array([]) 

        for j in range(len(ID)):

            dist1 = h1dist[j]*1000 #Mpc to kpc
            dist2 = h2dist[j]*1000 #Mpc to kpc
            R1 = R1vir[j]
            R2 = R2vir[j]

            cond1 = (dist1 < dist2)&(dist1 < R1)&(dist1 != -1000)
            cond2 = (dist1 > dist2)&(dist1 < R1)&(dist2 == -1000)
            cond3 = (dist2 < dist1)&(dist2 < R2)&(dist2 != -1000)
            cond4 = (dist2 > dist1)&(dist2 < R2)&(dist1 == -1000)

            if (cond1 or cond2):

                ida = np.append(ida, ID[j])
                vma = np.append(vma, Vmax[j])
                vpa = np.append(vpa, Vpeak[j])
                mpa = np.append(mpa, Mpeak[j])
                h1dista = np.append(h1dista, h1dist[j])
                h2dista = np.append(h2dista, h2dist[j])
                h1za = np.append(h1za, h1firstz_interact[j])
                h2za = np.append(h2za, h2firstz_interact[j])
                R1a = np.append(R1a, R1)
                R2a = np.append(R2a, R2)
                hosta = np.append(hosta, hostname[j])

            elif (cond3 or cond4):

                idb = np.append(idb, ID[j])
                vmb = np.append(vmb, Vmax[j])
                vpb = np.append(vpb, Vpeak[j])
                mpb = np.append(mpb, Mpeak[j])
                h1distb = np.append(h1distb, h1dist[j])
                h2distb = np.append(h2distb, h2dist[j])
                h1zb = np.append(h1zb, h1firstz_interact[j])
                h2zb = np.append(h2zb, h2firstz_interact[j])
                R1b = np.append(R1b, R1)
                R2b = np.append(R2b, R2)
                hostb = np.append(hostb, hostname[j])

            else:

                continue

        print i
        if (len(idb) != 0):
            
            cuta = np.random.choice(len(ida), number, replace = False) 
            cutb = np.random.choice(len(idb), number, replace = False)

            #select galaxies based on the above cuts
            ID_sel = np.append(ida[cuta],idb[cutb])
            Vmax_sel = np.append(vma[cuta],vmb[cutb])
            Vpeak_sel = np.append(vpa[cuta],vpb[cutb])
            Mpeak_sel = np.append(mpa[cuta],mpb[cutb])
            h1dist_sel = np.append(h1dista[cuta],h1distb[cutb])
            h2dist_sel = np.append(h2dista[cuta],h2distb[cutb])
            h1firstz_sel = np.append(h1za[cuta],h1zb[cutb])
            h2firstz_sel = np.append(h2za[cuta],h2zb[cutb])
            radius1 = np.append(R1a[cuta],R1b[cutb])
            radius2 = np.append(R2a[cuta],R2b[cutb])
            tag = np.append(hosta[cuta],hostb[cutb])

        else:

            cuta = np.random.choice(len(ida), number, replace = False)

            #select galaxies based on the above cuts
            ID_sel = ida[cuta]
            Vmax_sel = vma[cuta]
            Vpeak_sel = vpa[cuta]
            Mpeak_sel = mpa[cuta]
            h1dist_sel = h1dista[cuta]
            h2dist_sel = h2dista[cuta]
            h1firstz_sel = h1za[cuta]
            h2firstz_sel = h2za[cuta]
            radius1 = R1a[cuta]
            radius2 = R2a[cuta]
            tag = hosta[cuta]

        
        #append each column into master columns 
        masterID = np.append(masterID, ID_sel)
        masterVmax = np.append(masterVmax, Vmax_sel)
        masterVpeak = np.append(masterVpeak, Vpeak_sel)
        masterMpeak = np.append(masterMpeak, Mpeak_sel)
        masterh1dist = np.append(masterh1dist, h1dist_sel)
        masterh2dist = np.append(masterh2dist, h2dist_sel)
        masterh1firstz = np.append(masterh1firstz, h1firstz_sel)
        masterh2firstz = np.append(masterh2firstz, h2firstz_sel)
        mastertag = np.append(mastertag, tag)
        masterR1 = np.append(masterR1, radius1)
        masterR2 = np.append(masterR2, radius2)
            
        #put individual columns into table
        output = Table([ID_sel, Vmax_sel, Vpeak_sel, Mpeak_sel, h1dist_sel, h2dist_sel, h1firstz_sel, h2firstz_sel, tag, radius1, radius2], names=('ID', 'Vmax(km/s)', 'Vpeak(km/s)', 'Mpeak(Msun)', 'h1dist(Mpc)', 'h2dist(Mpc)', 'h1firstz_interacted', 'h2firstz_interacted', 'HostName', 'R1vir(kpc)', 'R2vir(kpc)'), meta={'Mass Cuts': 'ELVIS'})

        outputfile = 'ELVIS_Data_Clean/elvis_'+inputdata[i]+'_clean_ver2rand.dat'
        output.write(outputfile, format = 'ascii')

    #put master columns into table
    mastertable = Table([masterID, masterVmax, masterVpeak, masterMpeak, masterh1dist, masterh2dist, masterh1firstz, masterh2firstz, mastertag, masterR1, masterR2], names=('ID', 'Vmax(km/s)', 'Vpeak(km/s)', 'Mpeak(Msun)', 'h1dist(Mpc)', 'h2dist(Mpc)', 'h1firstz_interacted', 'h2firstz_interacted', 'HostName', 'R1vir(kpc)', 'R2vir(kpc)'), meta={'Mass Cuts': 'ELVIS'})    

    masteroutputfile = 'ELVIS_Data_Clean/elvis_alldwarfs_clean_ver2rand.dat'
    mastertable.write(masteroutputfile, format = 'ascii')
