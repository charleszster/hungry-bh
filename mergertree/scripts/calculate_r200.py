#!/Users/chaz/anaconda/bin/python #This is the python interpreter I use on Mac
##!/home/charles/anaconda/bin/python 
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import pickle
import os
import scipy.optimize
import get_clusters
import analyze_clusters
import Constants

H0 = Constants.H0
WM = Constants.WM
WV = Constants.WV
G = Constants.G
cluster_folder = Constants.cluster_folder
galaxies_max_mass_fitted_file = '3_max_mass_galaxies_fitted.pkl'

def get_H(z):
    return H0*(WV + (1 - WV - WM)*(1 + z)**2 + WM*(1 + z)**3)**0.5

def get_r200(m200, H):
    return (m200*G/(100.*H**2))**(1./3.)

def run():
    F = open(os.path.join(cluster_folder, galaxies_max_mass_fitted_file), 'rb')
    galaxies_max_mass = pickle.load(F)
    F.close

    z = galaxies_max_mass[galaxies_max_mass.keys()[0]][0]
    time = galaxies_max_mass[galaxies_max_mass.keys()[0]][1]
    H_z = []
    for _z in z:
        H_z.append(get_H(_z))

    r200 = []
    for key, value in galaxies_max_mass.iteritems():
        for i, H in enumerate(H_z):
            r200.append(get_r200(value[2][i], H))
    print r200

if __name__ == '__main__':
    run()