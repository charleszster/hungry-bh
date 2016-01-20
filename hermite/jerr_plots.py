#! /home/charles/anaconda/bin/python

import numpy as np
import matplotlib.pyplot as plt
import ast
import os

jerr_file = open('jerr001', 'r')
jerr = jerr_file.readlines()
jerr = [line.rstrip().split() for line in jerr]
jerr = [[ast.literal_eval(entry) for entry in line] for line in jerr]
jerr = np.array(jerr)

plt.figure()
t = jerr[:,0]
r = jerr[:,4]
rc = jerr[:,9]
ax1 = plt.subplot(211)
ax1.plot(t, r, 'r', label='r')
plt.xlabel('t (Myr)')
plt.ylabel('r (pc)')
ax1.legend(loc='best')

ax2 = plt.subplot(212, sharex=ax1)
ax2.plot(t, rc, 'b', label='rc')
ax2.get_yaxis().get_major_formatter().set_useOffset(False)
plt.ylabel('rc (pc)')
ax2.legend(loc='best')
plt.suptitle('initial r = %.1f; initial rc = %.1f' % (r[0], rc[0]))

savepath = os.path.join(os.path.expanduser('~'), 'Dropbox', 'Columbia', 'Ostriker', 'Cannibalism', 'hermite', 'Plots')
plt.savefig(os.path.join(savepath, 'diffusion_plot.png'))
plt.show()
