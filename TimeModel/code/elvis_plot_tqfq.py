#this takes the output files from elvis_getquenchtine.py and plots tq = tq[fq]

import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt

def plot_velratio():

    #file parameters to loop through
    radii = np.array(['0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1.0'])
    vel = np.array(['0.3','0.4', '0.5','0.6', '0.7'])

    for i in range(len(vel)):

        data1 = Table.read('TimeModel/Output/elvis_fqtq_r'+radii[0]+'_v'+vel[i]+'.dat', format = 'ascii')
        data2 = Table.read('TimeModel/Output/elvis_fqtq_r'+radii[1]+'_v'+vel[i]+'.dat', format = 'ascii')
        data3 = Table.read('TimeModel/Output/elvis_fqtq_r'+radii[2]+'_v'+vel[i]+'.dat', format = 'ascii')
        data4 = Table.read('TimeModel/Output/elvis_fqtq_r'+radii[3]+'_v'+vel[i]+'.dat', format = 'ascii')
        data5 = Table.read('TimeModel/Output/elvis_fqtq_r'+radii[4]+'_v'+vel[i]+'.dat', format = 'ascii')
        data6 = Table.read('TimeModel/Output/elvis_fqtq_r'+radii[5]+'_v'+vel[i]+'.dat', format = 'ascii')
        data7 = Table.read('TimeModel/Output/elvis_fqtq_r'+radii[6]+'_v'+vel[i]+'.dat', format = 'ascii')
        data8 = Table.read('TimeModel/Output/elvis_fqtq_r'+radii[7]+'_v'+vel[i]+'.dat', format = 'ascii')

        fq1 = np.array(data1['fq'])
        fq2 = np.array(data2['fq'])
        fq3 = np.array(data3['fq'])
        fq4 = np.array(data4['fq'])
        fq5 = np.array(data5['fq'])
        fq6 = np.array(data6['fq'])
        fq7 = np.array(data7['fq'])
        fq8 = np.array(data8['fq'])
        tq1 = np.array(data1['tq'])
        tq2 = np.array(data2['tq'])
        tq3 = np.array(data3['tq'])
        tq4 = np.array(data4['tq'])
        tq5 = np.array(data5['tq'])
        tq6 = np.array(data6['tq'])
        tq7 = np.array(data7['tq'])
        tq8 = np.array(data8['tq'])

        plt.figure(i)

        plt.axis([0.2,1,0,12])
        plt.xlabel('quenched fraction')
        plt.ylabel('quenching timescale (Gyr)')
        plt.title('GK14AM: vmax/vpeak > '+vel[i])
        plt.text(0.3,7.0, r'$r\ <\ $'+radii[0]+'$\ R_{vir}$',color = 'b')
        plt.text(0.3,6.6, r'$r\ <\ $'+radii[1]+'$\ R_{vir}$',color = 'g')
        plt.text(0.3,6.2, r'$r\ <\ $'+radii[2]+'$\ R_{vir}$',color = 'r')
        plt.text(0.3,5.8, r'$r\ <\ $'+radii[3]+'$\ R_{vir}$',color = 'k')
        plt.text(0.3,5.4, r'$r\ <\ $'+radii[4]+'$\ R_{vir}$',color = 'c')
        plt.text(0.3,5.0, r'$r\ <\ $'+radii[5]+'$\ R_{vir}$',color = 'm')
        plt.text(0.3,4.6, r'$r\ <\ $'+radii[6]+'$\ R_{vir}$',color = 'DarkOrange')
        plt.text(0.3,4.2, r'$r\ <\ $'+radii[7]+'$\ R_{vir}$',color = 'y')

        plt.plot(fq1, tq1, 'b-', fq2, tq2, 'g-', fq3, tq3, 'r-', fq4, tq4, 'k-', fq5, tq5, 'c-', fq6, tq6, 'm-', fq8, tq8, 'y-')
        plt.plot(fq7, tq7, color = 'DarkOrange', linestyle = '-')
        plt.savefig('TimeModel/MyPlots/elvis_fqtqplot_GK14AM_v'+vel[i]+'.pdf')
    plt.show()


def plot_radial():

    #file parameters to loop through
    radii = np.array(['1.0'])#'0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1.0'])
    vel = np.array(['0.3','0.4', '0.5','0.6', '0.7'])

    for i in range(len(radii)):

        data1 = Table.read('TimeModel/Output/elvis_fqtq_r'+radii[i]+'_v'+vel[0]+'.dat', format = 'ascii')
        data2 = Table.read('TimeModel/Output/elvis_fqtq_r'+radii[i]+'_v'+vel[1]+'.dat', format = 'ascii')
        data3 = Table.read('TimeModel/Output/elvis_fqtq_r'+radii[i]+'_v'+vel[2]+'.dat', format = 'ascii')
        data4 = Table.read('TimeModel/Output/elvis_fqtq_r'+radii[i]+'_v'+vel[3]+'.dat', format = 'ascii')
        data5 = Table.read('TimeModel/Output/elvis_fqtq_r'+radii[i]+'_v'+vel[4]+'.dat', format = 'ascii')

        fq1 = np.array(data1['fq'])
        fq2 = np.array(data2['fq'])
        fq3 = np.array(data3['fq'])
        fq4 = np.array(data4['fq'])
        fq5 = np.array(data5['fq'])
        tq1 = np.array(data1['tq'])
        tq2 = np.array(data2['tq'])
        tq3 = np.array(data3['tq'])
        tq4 = np.array(data4['tq'])
        tq5 = np.array(data5['tq'])

        plt.figure(i)

        plt.axis([0.2,1,0,12])
        plt.xlabel('quenched fraction')
        plt.ylabel('quenching timescale (Gyr)')
        plt.title(r'$\rm GK14AM:\ r\ <\ $'+radii[i]+r'$\rm R_{vir}$')
        plt.text(0.3,7.0, r'$\rm V_{max}/V_{peak}\ >\ $'+vel[0],color = 'b')
        plt.text(0.3,6.6, r'$\rm V_{max}/V_{peak}\ >\ $'+vel[1],color = 'g')
        plt.text(0.3,6.2, r'$\rm V_{max}/V_{peak}\ >\ $'+vel[2],color = 'r')
        plt.text(0.3,5.8, r'$\rm V_{max}/V_{peak}\ >\ $'+vel[3],color = 'k')
        plt.text(0.3,5.4, r'$\rm V_{max}/V_{peak}\ >\ $'+vel[4],color = 'c')

        plt.plot(fq1, tq1, 'b-', fq2, tq2, 'g-', fq3, tq3, 'r-', fq4, tq4, 'k-', fq5, tq5, 'c-')
        plt.savefig('TimeModel/MyPlots/elvis_fqtqplot_GK14AM_r'+radii[i]+'.pdf')
    plt.show()


def plot_abmatch():
    
    #file parameters to loop through
    #radii = np.array(['0.6', '0.8', '1.0'])
    #vel = np.array(['0.3', '0.5', '0.7'])

    #for i in range(len(vel)):

    data1 = Table.read('TimeModel/Output/elvis_fqtq_r1.0_v0.3.dat', format = 'ascii')
    data2 = Table.read('TimeModel/Output/elvis_fqtq_r1.0_v0.4.dat', format = 'ascii')
    data3 = Table.read('TimeModel/Output/elvis_fqtq_r1.0_v0.5.dat', format = 'ascii')
    data4 = Table.read('TimeModel/Output/elvis_fqtq_r1.0_v0.6.dat', format = 'ascii')
    data5 = Table.read('TimeModel/Output/elvis_fqtq_r1.0_v0.7.dat', format = 'ascii')
    data1r = Table.read('TimeModel/Output/elvis_fqtq_r1.0_v0.3_ver2rand.dat', format = 'ascii')
    data2r = Table.read('TimeModel/Output/elvis_fqtq_r1.0_v0.4_ver2rand.dat', format = 'ascii')
    data3r = Table.read('TimeModel/Output/elvis_fqtq_r1.0_v0.5_ver2rand.dat', format = 'ascii')
    data4r = Table.read('TimeModel/Output/elvis_fqtq_r1.0_v0.6_ver2rand.dat', format = 'ascii')
    data5r = Table.read('TimeModel/Output/elvis_fqtq_r1.0_v0.7_ver2rand.dat', format = 'ascii')
    
    fq1 = np.array(data1['fq'])
    fq2 = np.array(data2['fq'])
    fq3 = np.array(data3['fq'])
    fq4 = np.array(data4['fq'])
    fq5 = np.array(data5['fq'])
    fq1r = np.array(data1r['fq'])
    fq2r = np.array(data2r['fq'])
    fq3r = np.array(data3r['fq'])
    fq4r = np.array(data4r['fq'])
    fq5r = np.array(data5r['fq'])
    tq1 = np.array(data1['tq'])
    tq2 = np.array(data2['tq'])
    tq3 = np.array(data3['tq'])
    tq4 = np.array(data4['tq'])
    tq5 = np.array(data5['tq'])
    tq1r = np.array(data1r['tq'])
    tq2r = np.array(data2r['tq'])
    tq3r = np.array(data3r['tq'])
    tq4r = np.array(data4r['tq'])
    tq5r = np.array(data5r['tq'])
    
    plt.figure(1)
    
    plt.axis([0.2,1,0,12])
    plt.xlabel('quenched fraction')
    plt.ylabel('quenching timescale (Gyr)')
    plt.title(r'$\rm Various\ Abundance\ Matching:\ r\ <\ 1.0\ R_{vir}$')

    #pp1 = plt.Rectangle((0, 0), 1, 1, fc="b", ec = "b")
    #pp2 = plt.Rectangle((0, 0), 1, 1, fc="g", ec = "g")
    #pp3 = plt.Rectangle((0, 0), 1, 1, fc="r", ec = "r")
    #pp4 = plt.Rectangle((0, 0), 1, 1, fc="k", ec = "k")
    #pp5 = plt.Rectangle((0, 0), 1, 1, fc="c", ec = "c")
    pp0 = plt.Rectangle((0, 0), 1, 1, fc="w", ec="w")
    
    p1, = plt.plot(fq1, tq1, 'b-')
    p2, = plt.plot(fq2, tq2, 'g-')
    p3, = plt.plot(fq3, tq3, 'r-')
    p4, = plt.plot(fq4, tq4, 'k-')
    p5, = plt.plot(fq5, tq5, 'c-')
    p1r, = plt.plot(fq1r, tq1r, 'b--')
    p2r, = plt.plot(fq2r, tq2r, 'g--')
    p3r, = plt.plot(fq3r, tq3r, 'r--')
    p4r, = plt.plot(fq4r, tq4r, 'k--')
    p5r, = plt.plot(fq5r, tq5r, 'c--')

    plt.legend([pp0,p1,p2,p3,p4,p5,pp0,p1r,p2r,p3r,p4r,p5r], [r'$\rm GK14AM:\ 5*10^{9}\ <\ M_{halo}\ <\ 1*10^{11}$', r'$\rm Vmax/Vpeak\ >\ 0.3$', r'$\rm Vmax/Vpeak\ >\ 0.4$', r'$\rm Vmax/Vpeak\ >\ 0.5$', r'$\rm Vmax/Vpeak\ >\ 0.6$', r'$\rm Vmax/Vpeak\ >\ 0.7$',r'$\rm random\ select:\ 2*10^{9}\ <\ M_{halo}\ <\ 2*10^{11}$', r'$\rm Vmax/Vpeak\ >\ 0.3$', r'$\rm Vmax/Vpeak\ >\ 0.4$', r'$\rm Vmax/Vpeak\ >\ 0.5$', r'$\rm Vmax/Vpeak\ >\ 0.6$',r'$\rm Vmax/Vpeak\ >\ 0.7$'], loc = 3, frameon=False, numpoints=1, prop={'size':9})
    
    plt.savefig('TimeModel/MyPlots/elvis_fqtqplot_allv_r1.0.pdf')
    plt.show()
