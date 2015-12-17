#!/Users/chaz/anaconda/bin/python

import os
import numpy as np
import sys
import matplotlib.pyplot as plt

directory = os.path.join(os.path.expanduser('~'), 'Dropbox', 'Columbia', 'Ostriker', 'Cannibalism', 'hermite')
plum_files = [f for f in os.listdir(directory) if 'plum0' in f]
if len(plum_files) == 0:
   sys.exit('No data files to plot!  Aborting')
data = []
for i, file in enumerate(plum_files):
    data.append([])
    full_path = os.path.join(directory, file)
    plum_data = open(full_path, 'r').read()
    data[i].append([line.split() for line in plum_data.splitlines()])
    data[i][0] = [[float(piece) for piece in line] for line in data[i][0]]
    sys.stdout.write('Done importing data from %i of %i\r' %(i+1, len(plum_files)))
    sys.stdout.flush()
data = np.array(data)

r_max = None
t_max = None
plt.figure()
for plum_file in data:
    time = [i[0] for i in plum_file[0]]
    r = [np.sqrt(i[1]**2 + i[2]**2 + i[3]**2) for i in plum_file[0]]
    if np.amax(r) > r_max:
        r_max = np.amax(r)
    plt.plot(time, r, label='mass=%.2E, x=100, vx=206' % (plum_file[0][0][7]))
    
    if np.amax(time) > t_max:
        t_max = np.amax(time)
plt.xlabel('time, Myr')
plt.ylabel('Distance from center, pc')
plt.legend(loc=4, fontsize='small')
plt.xlim(0, t_max)
plt.ylim(0, r_max)
plt.show()
