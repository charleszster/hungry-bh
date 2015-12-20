#!/home/charles/anaconda/bin/python 
##!/Users/chaz/anaconda/bin/python #This is the python interpreter I use on Mac
# -*- coding: utf-8 -*-

import time
import datetime
import ast
import numpy as np
import matplotlib.pyplot as plt
import pickle
import os
import sys

def run():
    with open('jerr001', 'r') as f:
        results = f.readlines()
    results = [line.rstrip().split() for line in results]
    results = [[ast.literal_eval(entry) for entry in line] for line in results]
    results = np.array(results)
    t = [line[0] for line in results if line[0]<=13720.]
    r = [np.sqrt(line[1]**2+line[2]**2+line[3]**2) for line in results]
    r = r[:len(t)]

    plt.figure()
    plt.plot(t, r)
    plt.xlabel('Time (Myr)')
    plt.ylabel('Radial Distance (pc)')
    plt.ylim([min(np.amin(r),0.), np.amax(r)*1.1])
    curr_time = datetime.datetime.now()
    curr_time = curr_time.strftime('%Y_%m_%d_%H_%M_%S')
    save_location = os.path.join(os.path.expanduser('~'),'Dropbox', 'Columbia', 'Ostriker', 'Cannibalism', 'hermite', 'Plots')
    filename = ''.join(['hermite_radial_plot_', curr_time, '.png'])
    full_name = os.path.join(save_location, filename)
    plt.savefig(full_name)
    plt.show()

if __name__ == '__main__':
    run()
