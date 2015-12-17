#!/Users/chaz/anaconda/bin/python  #This is the python interpreter I use on Mac
##!/home/charles/anaconda/bin/python 

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize


rc = 100.
r_half = 265.339e3
M200 = 1.e13

r_h = np.arange(265.e3, 500.e3, 100.)

for rh in r_h:
    rho_c = M200*(rh + rc)/(2*np.pi**2*rc**2*rh**2)
    M = (4*np.pi*rc**2*rh**2*rho_c/(rh**2-rc**2))*(rh*np.arctan(r_half/rh)-rc*np.arctan(r_half/rc))
    if M/(0.5*M200) > 0.999 and M/(0.5*M200) < 1.001:
        print M/(0.5*M200), rh