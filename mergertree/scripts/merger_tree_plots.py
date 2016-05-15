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
cluster_folder = Constants.cluster_folder

def get_orbiting_bhs(redshifts, galaxy):
    orbiting_bhs = []
    for redshift in redshifts:
        for line in galaxy:
            if line[0] == redshift[0]:
                if len(line) >= 7:
                    orbiting_bhs.append([redshift[0], redshift[1], line[6]])
                else:
                    orbiting_bhs.append([redshift[0], redshift[1], [np.nan]])
    return orbiting_bhs

def get_masses(orbiting_bhs, smbh_cluster, key):
    galaxy_bh_mass = []
    central_bh_mass = []
    with open(os.path.join(cluster_folder, 'Central_BH_mass_galaxy_%s' % (key)), 'wb') as cbh:
        cbh.write('Time (Gyr)\tMass (Msun)\n')
#        print orbiting_bhs
        for time_slice in orbiting_bhs:
            redshift = time_slice[0]
            bhs = time_slice[2]
#            print bhs
            mass_this_redshift = 0.
            for line in smbh_cluster:
                smbh_redshift = line[0]
                smbh_bh = line[2]
                accreted_plus_seed_mass = line[4]
                if smbh_redshift == redshift and smbh_bh == int(key):
                    central_bh_accreted_mass = line[4]
                    cbh.write('%s\t%s\n' % (line[1], line[4]))
                    central_bh_mass.append([line[1], line[4]])
                    mass_this_redshift += central_bh_accreted_mass
                else:
                    for bh in bhs:
                        if smbh_redshift == redshift and bh == smbh_bh:
                            mass_this_redshift += accreted_plus_seed_mass
            galaxy_bh_mass.append([time_slice[1], mass_this_redshift])
            if redshift == 0.0:
                with open(os.path.join(cluster_folder, 
                                       'Infalling_BH_masses_galaxy_%s' % (key)), 'wb') as f:
#                    print key, redshift, bhs
                    f.write('BH ID\tInfall mass [Msun]\tInfall time [Myr]\n')
                    bh_masses = []
                    for line in smbh_cluster:
                        if line[0] == redshift:
                            for bh in bhs:
                               if bh == line[2]:
                                   f.write('%s\t%s\t%s\n' % (bh, line[4], line[7]*1000.))
                                   bh_masses.append(line[4])
                plt.figure()
                plt.hist(np.log10(bh_masses), bins=30)
                plt.xlabel('BH Mass (Msun), log10')
                plt.ylabel('Number')
                plt.title('Histogram of black hole masses for Galaxy %s at z=0.0' % (key))
                plt.savefig(os.path.join(plots_folder, 'Histogram_bh_masses_z_0_galaxy_%s' % (key)))
                plt.close()
    return galaxy_bh_mass, central_bh_mass, bh_masses

def plot_central_bh_mass(time, mass, key, label):
    plt.semilogy(time, mass, label=label)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='g', linestyle='-')
    plt.grid(b=True, which='minor', color='r', linestyle='--')
    plt.legend(loc='best')
    plt.xlabel('Time (Gyr)')
    plt.ylabel('Mass (Msun)')
    plt.title('Accreted + seed mass of Central BH for Galaxy %s' % (key))

def plot_galaxy_bh_mass(galaxy_bh_mass, key):
    plt.figure()
    plt.plot([line[0] for line in galaxy_bh_mass], [np.log10(line[1]) for line in galaxy_bh_mass])
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='g', linestyle='-')
    plt.grid(b=True, which='minor', color='r', linestyle='--')
    plt.xlabel('Time (Gyr)')
    plt.ylabel('Black Hole Mass (Msun, log10)')
    plt.title('Accreted Black Hole mass for Galaxy \'%s\'' % (key))
    plt.savefig(os.path.join(plots_folder, 'accreted_bh_mass_for_galaxy_%s' % (key)))
    plt.close()

def curve_fit_central_bh_masses(key, central_bh_mass):
    central_bh_mass_fitted = {}
    time = central_bh_mass[:, 0]
    mass_log10 = np.log10(central_bh_mass[:, 1])
    sigma = np.ones(len(time))
    sigma[[0, -1]] = 0.01
    func = analyze_clusters.exp_func
    popt, pcov = scipy.optimize.curve_fit(func, time, mass_log10, sigma=sigma)
    mass_log10_est = np.array([func(x, *tuple(popt)) for x in time])
    central_bh_mass_fitted = [time, [10**y for y in mass_log10_est]]  #[time, mass_log10_est]
    return popt, central_bh_mass_fitted

def fit_and_plot_central_bh_masses(smbh_cluster, galaxies_cluster_no_bad_z):
    f = open(os.path.join(cluster_folder, '3_max_mass_galaxies.pkl'), 'rb')
    max_mass_galaxies = pickle.load(f)
    f.close()

    max_mass_galaxies_ids = max_mass_galaxies.keys()
    galaxies = {}
    for id in max_mass_galaxies_ids:
        galaxies[str(id)] = [line for line in galaxies_cluster_no_bad_z if line[2] == int(id)]

    for key, value in galaxies.iteritems():
        print 'Plotting BH mass for galaxy %s' % (key)
        redshifts = [[line[0], line[1]] for line in galaxies[key]]

        orbiting_bhs = get_orbiting_bhs(redshifts, galaxies[key])

        galaxy_bh_mass, central_bh_mass, bh_masses = get_masses(orbiting_bhs, smbh_cluster, key)
        central_bh_mass = np.array(central_bh_mass)

        popt, central_bh_mass_fitted = curve_fit_central_bh_masses(key, central_bh_mass)
        print popt
        plt.figure()
        label = 'Orig data'
        plot_central_bh_mass(central_bh_mass[:, 0], central_bh_mass[:, 1], key, label)
        ''' THIS PORTION WAS ONLY TO TEST WHETHER THE FORTRAN-GENERATED TIME AND MASS MATCHED SAME GENERATED IN PYTHON. ONLY TESTED IT FOR GALAXY 1.
        if key=='1':
            f = open(os.path.join(os.path.expanduser('~'), 'Dropbox', 'Columbia', 'Ostriker', 'Cannibalism', 'hermite',
                              'mass.txt'), 'rb')
            time_mass = f.readlines()
            time = []
            mass = []
            for line in time_mass:
                time.append(float(line.split()[0]))
                mass.append(float(line.split()[1]))
            time = np.array(time)
            mass = np.array(mass)
            time_mass = [time, mass]
            label = 'Fortran'
            plot_central_bh_mass(time, mass, key, label)
        '''
        label = 'Fitted data'
        plot_central_bh_mass(central_bh_mass_fitted[0], central_bh_mass_fitted[1], key, label)
        plt.savefig(os.path.join(plots_folder, 'Central_bh_mass_galaxy_%s_orig_and_fitted' % (key)))
        plt.show()

def run():
    H0 = Constants.H0
    WM = Constants.WM
    WV = Constants.WV

#    smbh_cluster = get_clusters.process_smbh_cluster(H0, WM, WV)
    smbh_cluster = get_clusters.get_pickled_file('smbh_cluster.pkl')

#    galaxies_cluster = get_clusters.process_galaxies_cluster(H0, WM, WV)
    galaxies_cluster_no_bad_z = get_clusters.get_pickled_file('galaxies_cluster_no_bad_z.pkl')

#    galaxy_1 = [line for line in galaxies_cluster if line[2]==1]
#    black_holes = [time_slice[-1] for time_slice in galaxy_1 if time_slice[0]<=0.0][0]
#    print black_holes

    galaxies_by_id = analyze_clusters.get_galaxies_by_id(galaxies_cluster_no_bad_z)
    galaxies_masses, final_masses = analyze_clusters.get_galaxies_masses(galaxies_by_id)
#    print galaxies_masses['1']
    print galaxies_by_id['1']

    analyze_clusters.plot_galaxies_all_masses(galaxies_masses)
    '''
    top_masses = 3
    galaxies_max_mass = analyze_clusters.get_most_massive_galaxies(galaxies_masses, final_masses, top_masses)
#    print galaxies_max_mass

    F1 = open(os.path.join(get_clusters.cluster_folder, ''.join([str(top_masses), '_max_mass_galaxies.pkl'])), 'wb')
    pickle.dump(galaxies_max_mass, F1)
    F1.close()
#    analyze_clusters.plot_max_galaxies_masses(galaxies_max_mass)

    popt_galaxies_masses, galaxies_max_mass_fitted = analyze_clusters.curve_fit_most_massive_galaxies(galaxies_max_mass)
    print popt_galaxies_masses

    F = open(os.path.join(get_clusters.cluster_folder,
                          ''.join([str(top_masses), '_max_mass_galaxies_fitted.pkl'])), 'wb')
    pickle.dump(galaxies_max_mass_fitted, F)
    F.close()
####################################################
    f = open(os.path.join(os.path.expanduser('~'), 'Dropbox', 'Columbia', 'Ostriker', 'Cannibalism', 'hermite',
                      'galaxy_mass.txt'), 'rb')
    time_mass = f.readlines()
    time = []
    mass = []
    for line in time_mass:
        time.append(float(line.split()[0]))
        mass.append(float(line.split()[1]))
    time = np.array(time)
    garbage = []
    mass = np.array(mass)
    galaxies_max_mass_fitted['65'] = [time, garbage, mass]

    analyze_clusters.plot_max_masses_galaxies_orig_and_fitted(galaxies_max_mass, galaxies_max_mass_fitted, top_masses)

#    fit_and_plot_central_bh_masses(smbh_cluster, galaxies_cluster_no_bad_z)
    '''
if __name__ == '__main__':
    run()
