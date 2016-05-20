#!/Users/chaz/anaconda/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import itertools
import operator
import os
import scipy
import Constants

plots_folder = Constants.plots_folder

def poly_func(x, a, b, c, d, e, f, g, h):
#    return a*x**2. + b*x + c
#    return a*x**3. + b*x**2. + c*x + d
#    return a*x**4. + b*x**3. + c*x**2. + d*x + e
#    return a*x**5. + b*x**4. + c*x**3. + d*x**2. + e*x + f
#    return a*x**6. + b*x**5. + c*x**4. + d*x**3. + e*x**2. + f*x + g
    return a*x**7. + b*x**6. + c*x**5. + d*x**4. + e*x**3. + f*x**2. + g*x + h #USE THIS TO FIT GALAXY MASSES

def exp_func(x, a, b):
    return a*np.exp(b/x) #USE THIS TO FIT CENTRAL BLACK HOLE MASSES

def get_galaxies_by_id(galaxies_cluster):
    galaxies_by_id = {}
    for k, g in itertools.groupby(galaxies_cluster, key=operator.itemgetter(5)):
        group = list(g)
#        time = [t[1] for t in group]
#        galaxy_mass = [m[3]+m[4] for m in group]
        galaxies_by_id[str(k)] = group
    return galaxies_by_id

def get_smbh_by_id(smbh_cluster):
    smbh_by_id = {}
    for k, g in itertools.groupby(smbh_cluster, key=operator.itemgetter(2)):
        group = list(g)
        smbh_by_id[str(k)] = group
    return smbh_by_id

def get_cbh_accreted_plus_seed_mass(smbh_by_id, galaxy_num):
    return [[line[1], line[4]] for line in smbh_by_id[galaxy_num]]

def get_galaxies_masses(galaxies_by_id):
    galaxies_masses = {}
    final_masses = {}
    for k, g in galaxies_by_id.iteritems():
        z = [tz[0] for tz in g]
        years = [t[1] for t in g]
        galaxy_mass = [m[3]+m[4] for m in g]
        galaxies_masses[str(k)] = [z, years, galaxy_mass]
        if years[-1] >= 13.7226:
            final_masses[str(k)] = [z[-1], years[-1], galaxy_mass[-1]]
    return galaxies_masses, final_masses

def get_most_massive_galaxies(galaxies_masses, final_masses, top_masses):
    max_masses_ids = []
    for i in range(top_masses):
        galaxy_id_max_mass = max(final_masses.iteritems(), key=operator.itemgetter(1))[0]
        max_masses_ids.append(galaxy_id_max_mass)
        del final_masses[galaxy_id_max_mass]
    galaxies_max_mass = {}
    for id in max_masses_ids:
        galaxies_max_mass.update({str(id): galaxies_masses[str(id)]})
    return galaxies_max_mass

def get_galaxy_stellar_mass(galaxies_by_id, galaxy_num):
    stellar_mass = [time_slice[3] for time_slice in galaxies_by_id[galaxy_num]]
    return stellar_mass

def plot_galaxies_all_masses(galaxies_masses):
    plt.figure()
    for key, values in galaxies_masses.iteritems():
        plt.plot(values[1], np.log10(values[2]))
    plt.xlabel('Time (Gyr)')
    plt.ylabel('Galaxy Mass (Msun), log10')
    plt.savefig(os.path.join(plots_folder, 'galaxies_masses_v_time_no_bad_z.png'))
    plt.close()

def plot_max_galaxies_masses(galaxies_max_mass):
    plt.figure()
    for key, value in galaxies_max_mass.iteritems():
        time = value[1] #[t for t in galaxies_masses[max_masses_ids[i]][1]]
        mass = value[2] #[m for m in galaxies_masses[max_masses_ids[i]][2]]
        plt.plot(np.log10(time), np.log10(mass), label=''.join(['Galaxy ID ', key]))
    plt.xlabel('Time (Gyr), log10')
    plt.ylabel('Galaxy Mass (Msun), log10')
    plt.title('Galaxy Mass Evolution for ' + str(len(galaxies_max_mass)) + ' most Massive Galaxies at z=0.0; IDs = ' +
              str(galaxies_max_mass.keys()) + ' (bad z removed)', fontsize=10)
    ttl = plt.gca().title
    ttl.set_position([.5, 1.03])
    plt.legend(loc='best')
    plt.savefig(os.path.join(plots_folder, 'most_massive_v_time_no_bad_z_double_log.png'))
    plt.close()

def curve_fit_most_massive_galaxies(galaxies_max_mass):
    galaxies_max_mass_fitted = {}
    popt = {}
    for k, v in galaxies_max_mass.iteritems():
        xdata = v[1]
        ydata_log10 = np.log10(v[2])
        sigma = np.ones(len(xdata))
        sigma[[0, -1]] = 0.01
        popt[k], pcov = scipy.optimize.curve_fit(poly_func, xdata, ydata_log10, sigma=sigma)
        ydata_est = [poly_func(x, *tuple(popt[k])) for x in xdata]
        galaxies_max_mass_fitted[k] = [v[0], v[1], [10**y for y in ydata_est]]
    return popt, galaxies_max_mass_fitted

def curve_fit_any_galaxy_mass(galaxies_masses, galaxy_num):
    xdata = galaxies_masses[galaxy_num][1]
    ydata_log10 = np.log10(galaxies_masses[galaxy_num][2])
    sigma = np.ones(len(xdata))
    sigma[[0, -1]] = 0.01
    popt, pcov = scipy.optimize.curve_fit(poly_func, xdata, ydata_log10, sigma=sigma)
    return popt

def plot_max_masses_galaxies_orig_and_fitted(galaxies_max_mass, galaxies_max_mass_fitted, top_masses):
    for k, v in galaxies_max_mass.iteritems():
        print('plotting galaxy id ', k)
        xdata = v[1]
        ydata_log10 = np.log10(v[2])
        if k=='65':
            plt.figure()
            plt.plot(xdata, ydata_log10, label='Original')
            ydata_est = np.log10(galaxies_max_mass_fitted[k][2])
            plt.plot(xdata, ydata_est, label='Fitted')
            plt.legend(loc='best')
            plt.xlabel('Time (Gyr)')
            plt.ylabel('Galaxy Mass (Msun), log10')
            plt.savefig(os.path.join(plots_folder, ''.join(['curve_fit_for_galaxy_id_', k, '.png'])))
