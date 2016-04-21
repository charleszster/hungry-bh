#! /Users/chaz/anaconda/envs/py35/bin/python
##! /home/charles/anaconda/envs/py35/bin/python

import numpy as np
import matplotlib.pyplot as plt
import itertools
import ast

def run():
    output = []
    with open('test.out', 'r') as f:
        skippedlines = itertools.islice(f, 0, None, 30)
        for line in skippedlines:
            output.append(line.split())
    output = [[ast.literal_eval(entry) for entry in line] for line in output]
    little_bh_data = []
    for line in output:
        little_bh_data.append([line[0], line[7]])
    little_bh_data = np.array(little_bh_data)
    print(little_bh_data)
    plt.figure()
    plt.plot(little_bh_data[:,0], little_bh_data[:,1])
    plt.xlabel('time (Myr)')
    plt.ylabel('R (pc)')
    plt.show()

if __name__ == '__main__':
    run()
