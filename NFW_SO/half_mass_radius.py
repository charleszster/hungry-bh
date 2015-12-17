#!/Users/chaz/anaconda/bin/python  #This is the python interpreter I use on Mac
##!/home/charles/anaconda/bin/python 

import numpy as np
import matplotlib.pyplot as plt

xs = np.arange(0., 8.771, 0.0001)

for x in xs:
    m = np.log(1.+x) - x/(1.+x)
    if m/.99007644 < 1.001 and m/.99007644 > 0.999:
        print m, x