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
import Constants

plots_folder = Constants.plots_folder

def poly_func(x, a, b, c, d, e, f, g, h):
#    return a*x**2. + b*x + c
#    return a*x**3. + b*x**2. + c*x + d
#    return a*x**4. + b*x**3. + c*x**2. + d*x + e
#    return a*x**5. + b*x**4. + c*x**3. + d*x**2. + e*x + f
#    return a*x**6. + b*x**5. + c*x**4. + d*x**3. + e*x**2. + f*x + g
    return a*x**7. + b*x**6. + c*x**5. + d*x**4. + e*x**3. + f*x**2. + g*x + h

def exp_func(x, a, b):
    return a*np.exp(b/x) #USE THIS TO FIT CENTRAL BLACK HOLE MASSES

def get_func_coeffs(time, stellar_mass):
    stellar_mass_log10 = np.log10(stellar_mass)
    sigma = np.ones(len(time))
    sigma[[0, -1]] = 0.01
    popt, pcov = scipy.optimize.curve_fit(poly_func, time, stellar_mass_log10, sigma=sigma)
    return popt

def plot_stellar_mass(time, galaxy_num, stellar_mass, stellar_mass_fitted):
#    plt.figure()
    t_plotted = [t/1000. for t in time]
    plt.semilogy(t_plotted, stellar_mass, label='Orig Data')
    plt.semilogy(t_plotted, stellar_mass_fitted, label='Fitted Data')
    plt.grid(b=True, which='major', color='g', linestyle='-')
    plt.grid(b=True, which='minor', color='r', linestyle='--')
    plt.legend(loc='best')
    plt.xlabel('Time (yr)')
    plt.ylabel('Stellar Mass (Msun)')
    plt.title('Stellar mass for Galaxy %s' %(galaxy_num))
    plt.show()
#    plt.savefig(os.path.join(plots_folder, ''.join(['Stellar_Mass_Galaxy_%s' % (galaxy_num)])))
#    plt.close()

def run():
    galaxies_cluster_no_bad_z = get_clusters.get_pickled_file('galaxies_cluster_no_bad_z.pkl')
    galaxies_by_id = analyze_clusters.get_galaxies_by_id(galaxies_cluster_no_bad_z)
    smbh_cluster = get_clusters.get_pickled_file('smbh_cluster.pkl')
    for galaxy_num in ['1']:
        #time orig. in Gyr; use Myr to not have to convert it in hermite code
        time = [time_slice[1]*1.e3 for time_slice in galaxies_by_id[galaxy_num]]
        print(time)
        stellar_mass = analyze_clusters.get_galaxy_stellar_mass(galaxies_by_id, galaxy_num)
        popt = get_func_coeffs(time, stellar_mass)
        print 'Galaxy #%s:' % (galaxy_num)
        print popt
        stellar_mass_fitted = [10.**(poly_func(t, *tuple(popt))) for t in time]
        plot_stellar_mass(time, galaxy_num, stellar_mass, stellar_mass_fitted)


if __name__ == '__main__':
    run()
