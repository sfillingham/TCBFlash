#This takes the output of the get.quenchtime_individual() and plots the median points along
#with the standard deviation as the error bars.

import numpy as np
from astropy.table import Table
from astropy.stats import sigma_clip
import matplotlib.pyplot as plt

def stats(input_indiv):

    inputlist = Table.read(input_indiv, format = 'ascii')
    filename = np.array(inputlist['fqinput'])

    #loop through all files, creating master columns to do stats on
    tq2 = np.array([])
    tq3 = np.array([])
    tq4 = np.array([])
    tq5 = np.array([])
    tq6 = np.array([])
    tq7 = np.array([])
    tq8 = np.array([])
    tq9 = np.array([])

    for i in range(len(filename)):

        data = Table.read('TimeModel/Output/'+filename[i]+'.dat', format = 'ascii')

        fq = np.array(data['fq'])
        tq = np.array(data['tq'])

        tq2 = np.append(tq2,tq[0])
        tq3 = np.append(tq3,tq[1])
        tq4 = np.append(tq4,tq[2])
        tq5 = np.append(tq5,tq[3])
        tq6 = np.append(tq6,tq[4])
        tq7 = np.append(tq7,tq[5])
        tq8 = np.append(tq8,tq[6])
        tq9 = np.append(tq9,tq[7])

    filtered_tq2 = sigma_clip(tq2, 3, None)
    filtered_tq3 = sigma_clip(tq3, 3, None)
    filtered_tq4 = sigma_clip(tq4, 3, None)
    filtered_tq5 = sigma_clip(tq5, 3, None)
    filtered_tq6 = sigma_clip(tq6, 3, None)
    filtered_tq7 = sigma_clip(tq7, 3, None)
    filtered_tq8 = sigma_clip(tq8, 3, None)
    filtered_tq9 = sigma_clip(tq9, 3, None)
    
    med2 = np.median(filtered_tq2)
    med3 = np.median(filtered_tq3)
    med4 = np.median(filtered_tq4)
    med5 = np.median(filtered_tq5)
    med6 = np.median(filtered_tq6)
    med7 = np.median(filtered_tq7)
    med8 = np.median(filtered_tq8)
    med9 = np.median(filtered_tq9)

    sig2 = np.std(filtered_tq2)
    sig3 = np.std(filtered_tq3)
    sig4 = np.std(filtered_tq4)
    sig5 = np.std(filtered_tq5)
    sig6 = np.std(filtered_tq6)
    sig7 = np.std(filtered_tq7)
    sig8 = np.std(filtered_tq8)
    sig9 = np.std(filtered_tq9)

    median_tq = np.array([med2, med3, med4, med5, med6, med7, med8, med9])
    sig_tq = np.array([sig2, sig3, sig4, sig5, sig6, sig7, sig8, sig9])

    print(median_tq)
    print(sig_tq)

    plt.figure(1)

    plt.axis([0.1,1,0,12])
    plt.xlabel('quenched fraction')
    plt.ylabel('quenching timescale (Gyr)')
    plt.title('GK14AM')
    plt.text(0.7,10.5, r'$Vmax/Vpeak\ >\ 0.3$',color = 'k')
    plt.text(0.7,10.1, r'$r\ <\ 1.0\ R_{vir}$',color = 'k')
    
    plt.errorbar(fq, median_tq, yerr=sig_tq, fmt='bo')
    plt.savefig('TimeModel/MyPlots/elvis_fqtqplot_v0.3_r1.0_stats_3sigclip.pdf')
    plt.show()
