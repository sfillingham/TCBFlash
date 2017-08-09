# Wetzel's Data

import numpy as np
import matplotlib.pyplot as plt

wz = open("WetzelData_formatted.dat", "r")
wlines = wz.readlines()
wz.close


def get_wetzel_tinf(wlines):
    """
    
    """
    N = np.size(wlines)
    mstar = np.ndarray(N-1)
    mhalo1213lo = np.ndarray(N-1)
    mhalo1213hi = np.ndarray(N-1)
    mhalo1314lo = np.ndarray(N-1)
    mhalo1314hi = np.ndarray(N-1)
    mhalo1415lo = np.ndarray(N-1)
    mhalo1415hi = np.ndarray(N-1)
    
    for t in range(N-1):
        
        mstar[t] = float(wlines[t+1].strip().split()[0])
        mhalo1213lo[t] = float(wlines[t+1].strip().split()[1])
        mhalo1213hi[t] = float(wlines[t+1].strip().split()[2])
        mhalo1314lo[t] = float(wlines[t+1].strip().split()[3])
        mhalo1314hi[t] = float(wlines[t+1].strip().split()[4])
        mhalo1415lo[t] = float(wlines[t+1].strip().split()[5])
        mhalo1415hi[t] = float(wlines[t+1].strip().split()[6])

    plt.figure(1)
    plt.axes([5,12,0,10])
    plt.plot(mstar,mhalo1213lo,'g',mstar,mhalo1213hi,'g')
    plt.show()

    return mstar, mhalo1213lo, mhalo1213hi, mhalo1314lo, mhalo1314hi, mhalo1415lo, mhalo1415hi

    
    
