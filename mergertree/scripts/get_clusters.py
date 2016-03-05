#!/Users/chaz/anaconda/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import os
import ast
import sys
import CC
import pickle
import Constants

cluster_folder = os.path.join(os.path.expanduser('~'), 'Dropbox', 'Columbia', 'Ostriker', 'Cannibalism', 'mergertree',
                              'cluster_files')

def get_cluster_file(file):
    cluster_file = open(file, 'r')
    cluster = cluster_file.readlines()
    cluster = [line.rstrip().split() for line in cluster]
    cluster = [[ast.literal_eval(entry) for entry in line] for line in cluster]
    return cluster

def process_smbh_cluster(H0, WM, WV):
    smbh_cluster_file = 'SMBH_cluster'
    smbh_cluster = get_cluster_file(os.path.join(cluster_folder, smbh_cluster_file))
    for i, line in enumerate(smbh_cluster):
        sys.stdout.write(''.join(['\rConverting z to age in SMBH Cluster for #', str(i+1), ' of ',
                                   str(len(smbh_cluster))]))
        sys.stdout.flush()
        z = float(line[0])
        input_params= [z, H0, WM, WV]
        age = CC.get_output_params(input_params)
        line.insert(1, age)
    print
    smbh_cluster = np.array(smbh_cluster, dtype=object)
    smbh_cluster = sorted(smbh_cluster, key=lambda l:l[2])
    F = open(os.path.join(cluster_folder, 'smbh_cluster.pkl'), 'wb')
    pickle.dump(smbh_cluster, F)
    F.close()
    print('Saved SMBH Cluster array to pickle')
    print
    print(smbh_cluster[0])
    print(smbh_cluster[int(len(smbh_cluster)/4)])
    print(smbh_cluster[int(len(smbh_cluster)/2)])
    print(smbh_cluster[3*int(len(smbh_cluster)/4)])
    print(smbh_cluster[len(smbh_cluster)-1])
    return smbh_cluster

def process_galaxies_cluster(H0, WM, WV):
    galaxies_cluster_file = 'galaxies_cluster'
    galaxies_cluster = get_cluster_file(os.path.join(cluster_folder, galaxies_cluster_file))
    galaxies_cluster_length = len(galaxies_cluster)
    for i, line in enumerate(galaxies_cluster):
        if len(line) >= 6:
            line[5] = line[5:]
        if len(line) >= 7:
            line[6:] = []
        sys.stdout.write(''.join(['\rConverting z to age in Galaxy Cluster for #', str(i+1), ' of ',
                                   str(galaxies_cluster_length)]))
        sys.stdout.flush()
        z = float(line[0])
        input_params= [z, H0, WM, WV]
        age = CC.get_output_params(input_params)
        line.insert(1, age)
    print
    #galaxies_cluster = np.array([np.array(line, dtype=object) for line in galaxies_cluster], dtype=object)
    galaxies_cluster = np.array(galaxies_cluster, dtype=object)
    galaxies_cluster = sorted(galaxies_cluster, key=lambda l:l[2])
#    print(galaxies_cluster)
    F = open(os.path.join(cluster_folder, 'galaxies_cluster.pkl'), 'wb')
    pickle.dump(galaxies_cluster, F)
    F.close()
    print('Saved Galaxies Cluster array to pickle')
    return galaxies_cluster

def get_pickled_file(cluster_file):
    F = open(os.path.join(cluster_folder, cluster_file), 'rb')
    cluster = pickle.load(F)
    cluster = np.array(cluster, dtype=object)
    return cluster

def fix_galaxies_cluster(galaxies_cluster):
#This function removes the mysterious gigantic mass that each galaxy temporarily attains at z=0.3
    bad_z_indices = [i for i, gal in enumerate(galaxies_cluster) if gal[0]==0.3]
    galaxies_cluster_no_bad_z = np.delete(galaxies_cluster, bad_z_indices, axis=0)
    return galaxies_cluster_no_bad_z

def save_fixed_galaxies_cluster(galaxies_cluster_no_bad_z):
    F = open(os.path.join(cluster_folder, 'galaxies_cluster_no_bad_z.pkl'), 'wb')
    pickle.dump(galaxies_cluster_no_bad_z, F)
    F.close()
