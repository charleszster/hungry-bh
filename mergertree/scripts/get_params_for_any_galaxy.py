#!/Users/chaz/anaconda/bin/python  #This is the python interpreter I use on Mac
##!/home/charles/anaconda/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import pickle
import os
import pprint
import scipy.optimize
import get_clusters
import analyze_clusters
import galaxies_star_mass
import merger_tree_plots
import Constants

plots_folder = Constants.plots_folder
cluster_folder = Constants.clusters_folder

def get_galaxy_mass_coefficients(galaxies_masses, galaxy_num):
    popt = analyze_clusters.curve_fit_any_galaxy_mass(galaxies_masses, galaxy_num)
    time = [1.e9*t for t in galaxies_masses[galaxy_num][1]]
    gm_original = galaxies_masses[galaxy_num][2]
    ydata_est = [analyze_clusters.poly_func(x, *tuple(popt)) for x in galaxies_masses[galaxy_num][1]]
    gm_fitted = [10**y for y in ydata_est]    
    plt.figure()
    plt.loglog(time, gm_original)
    plt.loglog(time, gm_fitted)
    plt.xlabel('Time (Yr)')
    plt.ylabel('Mass (Msun)')
    plt.title('Galaxy Mass for Galaxy #%s' % (galaxy_num))
    plt.show()
    return popt

def get_star_mass_coefficients(galaxies_by_id, galaxy_num):
    #time orig. in Gyr; use Myr to not have to convert it in hermite code
    time = [time_slice[1]*1.e3 for time_slice in galaxies_by_id[galaxy_num]] 
    stellar_mass = [time_slice[3] for time_slice in galaxies_by_id[galaxy_num]]
    popt = galaxies_star_mass.get_func_coeffs(time, stellar_mass)
    stellar_mass_fitted = [10.**(galaxies_star_mass.poly_func(t, *tuple(popt))) for t in time]
    galaxies_star_mass.plot_stellar_mass(time, galaxy_num, stellar_mass, stellar_mass_fitted)
    plt.show()
    return popt

def run():
    H0 = Constants.H0
    WM = Constants.WM
    WV = Constants.WV

    smbh_cluster = get_clusters.get_pickled_file('smbh_cluster.pkl')
    smbh_by_id = analyze_clusters.get_smbh_by_id(smbh_cluster)

    galaxies_cluster_no_bad_z = get_clusters.get_pickled_file('galaxies_cluster_no_bad_z.pkl')
    galaxies_by_id = analyze_clusters.get_galaxies_by_id(galaxies_cluster_no_bad_z)

    galaxies_masses, final_masses = analyze_clusters.get_galaxies_masses(galaxies_by_id)

    galaxy_num = '80'
    if galaxy_num in galaxies_by_id.keys():
        gm_popt = get_galaxy_mass_coefficients(galaxies_masses, galaxy_num)
        print 'Galaxy Mass Coefficients for Galaxy #%s:' % (galaxy_num)
        print gm_popt    

#        bh_popt, central_bh_mass_fitted = merger_tree_plots.curve_fit_central_bh_masses(galaxy_num, central_bh_mass)
#        print 'Central Black Hole Mass Coefficients for Galaxy #%s:' % (galaxy_num)
#        print sm_popt

        sm_popt = get_star_mass_coefficients(galaxies_by_id, galaxy_num)
        print 'Stellar Mass Coefficients for Galaxy #%s:' % (galaxy_num)
        print sm_popt

    else:
        sys.exit('No galaxy with such a number!')

if __name__ == '__main__':
    run()