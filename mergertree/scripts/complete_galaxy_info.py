#!/Users/chaz/anaconda/bin/python  #This is the python interpreter I use on Mac
##!/home/charles/anaconda/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import pickle
import os
import ast
import get_clusters
import analyze_clusters
import merger_tree_plots
import galaxies_star_mass
import Constants

plots_folder = Constants.plots_folder
cluster_folder = Constants.cluster_folder

H0 = Constants.H0
WM = Constants.WM
WV = Constants.WV
t0 = 1.56537653+03 #Myr

def get_galaxy_mass_data(galaxies_masses, galaxy_num):
    galaxy_mass_coeffs = analyze_clusters.curve_fit_any_galaxy_mass(galaxies_masses, galaxy_num)
    print 'Galaxy Mass Coefficients (10**(a*t**7. + b*t**6. + c*t**5. + d*t**4. + e*t**3. + f*t**2. + g*t + h)), t in Gyr:'
    print galaxy_mass_coeffs, '\n\n'
    t_s = np.linspace(galaxies_masses[galaxy_num][1][0], galaxies_masses[galaxy_num][1][-1], 5000)
    masses_fitted = [analyze_clusters.poly_func(t, *tuple(galaxy_mass_coeffs)) for t in t_s]
    masses_fitted = [10**y for y in masses_fitted]
    return t_s, masses_fitted

def get_central_bh_mass_growth(galaxies_cluster_no_bad_z, smbh_cluster, galaxy_num):
    galaxy = [line for line in galaxies_cluster_no_bad_z if line[2] == int(galaxy_num)]
    redshifts = [[line[0], line[1]] for line in galaxy]
    orbiting_bhs = merger_tree_plots.get_orbiting_bhs(redshifts, galaxy)
    galaxy_bh_mass, central_bh_mass, bh_masses = merger_tree_plots.get_masses(orbiting_bhs, smbh_cluster, galaxy_num)
    central_bh_mass = np.array(central_bh_mass)
    central_bh_mass_coeffs, central_bh_mass_fitted = merger_tree_plots.curve_fit_central_bh_masses(galaxy_num,
                                                                                                   central_bh_mass)
    print 'Central Black Hole Mass Coefficients 10**((a*np.exp(b/t))), t in Gyr:'
    print central_bh_mass_coeffs, '\n\n'
    return central_bh_mass, central_bh_mass_fitted

def get_stellar_mass_growth(galaxies_by_id, galaxy_num):
    time = [time_slice[1]*1.e3 for time_slice in galaxies_by_id[galaxy_num]] 
    stellar_mass = [time_slice[3] for time_slice in galaxies_by_id[galaxy_num]]
    stellar_mass_coefficients = galaxies_star_mass.get_func_coeffs(time, stellar_mass)
    print 'Stellar Mass Growth Coefficients 10**(a*t**7. + b*t**6. + c*t**5. + d*t**4. + e*t**3. + f*t**2. + g*t + h), t in Myr:'
    print stellar_mass_coefficients
    stellar_mass_fitted = [10.**(galaxies_star_mass.poly_func(t, *tuple(stellar_mass_coefficients))) for t in time]
    return time, stellar_mass, stellar_mass_coefficients, stellar_mass_fitted

def plot_all_data(t_s, masses_fitted, galaxies_masses, galaxy_num, central_bh_mass, central_bh_mass_fitted, time,
                                                                                     stellar_mass, stellar_mass_fitted):
    fig = plt.figure()
    ax1 = fig.add_subplot(2, 2, 1)
    ax1.semilogy(np.array(galaxies_masses[galaxy_num][1]), np.array(galaxies_masses[galaxy_num][2]), label='Orig data')
    ax1.semilogy(t_s, masses_fitted, label='Fitted data')
    ax1.set_xlabel('Time (Gyr)')
    ax1.set_ylabel('Galaxy Mass (Msun)')
    ax1.grid(b=True, which='major', color='g', linestyle='-')
    ax1.grid(b=True, which='minor', color='r', linestyle='--')
    ax1.legend(loc='best')
    ax1.set_title(''.join(['Galaxy ', galaxy_num, ' Mass Evolution']))

    ax2 = fig.add_subplot(2, 2, 2)
    label = 'Orig data'
    merger_tree_plots.plot_central_bh_mass(central_bh_mass[:, 0], central_bh_mass[:, 1], galaxy_num, label)
    label = 'Fitted data'
    merger_tree_plots.plot_central_bh_mass(central_bh_mass_fitted[0], central_bh_mass_fitted[1], galaxy_num, label)

    ax3 = fig.add_subplot(2, 2, 3)
    galaxies_star_mass.plot_stellar_mass(time, galaxy_num, stellar_mass, stellar_mass_fitted)

    plt.show()

def get_z(timo):
    return 7.20196192*(timo/1000)**(-0.59331986) - 1.52145449

def get_stlr_mass(timo, smc):
    return 10.**(smc[0]*timo**7. + smc[1]*timo**6. + smc[2]*timo**5. + smc[3]*timo**4. + smc[4]*timo**3. + 
                                                                                 smc[5]*timo**2. + smc[6]*timo + smc[7])

def get_eff_rad(stlr_mass, z):
    return 2500.*(stlr_mass/1.e+11)**0.73*(1.+ z)**(-0.98)

def get_sigma_eff(stlr_mass, z):
    return 190.*(stlr_mass/1.e+11)**0.2*(1.+ z)**(0.47)

def get_t_fric(sigma_eff, MBH, eff_rad, r_drop):
    return (19./6.)*(sigma_eff/200.)*(1.e8/MBH)*(eff_rad/r_drop)**2.

def read_orbiting_bhs(galaxy_num, smc):
    with open(os.path.join(cluster_folder, 'Infalling_BH_masses_galaxy_%s' % (galaxy_num)), 'rb') as f:
        whole = f.readlines()
    whole = [part.split() for part in whole[1:]]
    whole = np.array([[ast.literal_eval(entry) for entry in line] for line in whole])
    whole.view('i8,i8,i8').sort(order=['f1'], axis=0)
    with open(os.path.join(cluster_folder, 'Infalling_BH_masses_galaxy_%s' % (galaxy_num)), 'wb') as f:
        f.write('BH ID\tInfall time [Gyr]\tInfall mass [Msun]\tt at Ctr [Gyr]\n')
        for line in whole:
            timo = t0+line[1]*1000
            z = get_z(timo)
            stlr_mass = get_stlr_mass(timo, smc)
            eff_rad = get_eff_rad(stlr_mass, z)
            sigma_eff = get_sigma_eff(stlr_mass, z)
            r_drop = eff_rad
            t_fric = get_t_fric(sigma_eff, line[2], eff_rad, r_drop)
            t_at_center = timo/1000 + t_fric
            f.write('%s\t%9.3f\t\t\t\t%11.0f\t\t%11.2f\n' % (str(int(line[0])).zfill(4), line[1], line[2],
                                                             t_at_center))

def run():
    smbh_cluster = get_clusters.get_pickled_file('smbh_cluster.pkl')
    galaxies_cluster_no_bad_z = get_clusters.get_pickled_file('galaxies_cluster_no_bad_z.pkl')

    galaxies_by_id = analyze_clusters.get_galaxies_by_id(galaxies_cluster_no_bad_z)
    galaxies_masses, final_masses = analyze_clusters.get_galaxies_masses(galaxies_by_id)
    galaxy_num = '1'
    t_s, masses_fitted = get_galaxy_mass_data(galaxies_masses, galaxy_num)

    central_bh_mass, central_bh_mass_fitted = get_central_bh_mass_growth(galaxies_cluster_no_bad_z, smbh_cluster,
                                                                                                             galaxy_num)
    time, stellar_mass, stellar_mass_coefficients, stellar_mass_fitted = get_stellar_mass_growth(galaxies_by_id,
                                                                                                             galaxy_num)
    
#    plot_all_data(t_s, masses_fitted, galaxies_masses, galaxy_num, central_bh_mass, central_bh_mass_fitted, time,
#                                                                                     stellar_mass, stellar_mass_fitted)
    read_orbiting_bhs(galaxy_num, stellar_mass_coefficients)

if __name__ == '__main__':
    run()