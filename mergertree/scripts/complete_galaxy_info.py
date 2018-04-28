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
import matplotlib as mpl

plots_folder = Constants.plots_folder
cluster_folder = Constants.cluster_folder

H0 = Constants.H0
WM = Constants.WM
WV = Constants.WV
t0 = 1.56537653e+03 #Myr
hubble_time = 13720e6 #MYr

def write_galaxies_final_masses(galaxies_masses):
    galaxies_masses_list = []
    print galaxies_masses['1']
    for galaxy_nums, galaxy_data in galaxies_masses.iteritems():
        if galaxy_data[1][-1] > 13.72:
            galaxies_masses_list.append([galaxy_nums, galaxy_data[2][-1], galaxy_data[0][0]])
    galaxies_masses_list = sorted(galaxies_masses_list, key=lambda l: l[1], reverse=True)
    with open(os.path.join(cluster_folder, 'galaxies_final_masses.txt'), 'wb') as f:
        f.write('Galaxy ID\t\tFinal Galaxy Mass (Msun)\tGalaxy Birth (z)\n')
        for line in galaxies_masses_list:
            f.write('%s\t\t\t%.0f\t\t\t\t%.2f\n' % (line[0].zfill(4), line[1], line[2]))

def get_galaxy_mass_data(galaxies_masses, galaxy_num):
    galaxy_mass_coeffs = analyze_clusters.curve_fit_any_galaxy_mass(galaxies_masses, galaxy_num)
#    print 'Galaxy Mass Coefficients (10**(a*t**7. + b*t**6. + c*t**5. + d*t**4. + e*t**3. + f*t**2. + g*t + h)), t in Gyr:'
#    print galaxy_mass_coeffs, '\n\n'
    t_s = np.linspace(galaxies_masses[galaxy_num][1][0], galaxies_masses[galaxy_num][1][-1], 5000)
    masses_fitted = [analyze_clusters.poly_func(t, *tuple(galaxy_mass_coeffs)) for t in t_s]
    masses_fitted = [10**y for y in masses_fitted]
    return t_s, masses_fitted, galaxy_mass_coeffs

def get_central_bh_mass_growth(galaxies_cluster_no_bad_z, smbh_cluster, galaxy_num):
    galaxy = [line for line in galaxies_cluster_no_bad_z if line[2] == int(galaxy_num)]
    redshifts = [[line[0], line[1]] for line in galaxy]
    orbiting_bhs = merger_tree_plots.get_orbiting_bhs(redshifts, galaxy)
    galaxy_bh_mass, central_bh_mass, bh_masses = merger_tree_plots.get_masses(orbiting_bhs, smbh_cluster, galaxy_num)
    central_bh_mass = np.array(central_bh_mass)
    central_bh_mass_coeffs, central_bh_mass_fitted = merger_tree_plots.curve_fit_central_bh_masses(galaxy_num,
                                                                                                   central_bh_mass)
#    print 'Central Black Hole Mass Coefficients 10**((a*np.exp(b/t))), t in Gyr:'
#    print central_bh_mass_coeffs, '\n\n'
    return central_bh_mass, central_bh_mass_fitted, central_bh_mass_coeffs, bh_masses

def get_stellar_mass_growth(galaxies_by_id, galaxy_num):
    time = [time_slice[1]*1.e3 for time_slice in galaxies_by_id[galaxy_num]]
    stellar_mass = [time_slice[3] for time_slice in galaxies_by_id[galaxy_num]]
    stellar_mass_coefficients = galaxies_star_mass.get_func_coeffs(time, stellar_mass)
#    print 'Stellar Mass Growth Coefficients 10**(a*t**7. + b*t**6. + c*t**5. + d*t**4. + e*t**3. + f*t**2. + g*t + h), t in Myr:'
#    print stellar_mass_coefficients
    stellar_mass_fitted = [10.**(galaxies_star_mass.poly_func(t, *tuple(stellar_mass_coefficients))) for t in time]
    return time, stellar_mass, stellar_mass_coefficients, stellar_mass_fitted

def output_coeffs_for_achain_h(galaxy_num, central_bh_mass_coeffs, galaxy_mass_coeffs, stellar_mass_coefficients):
    print 'C*********Parameters used for calculating the central black hole mass (mbch*), galaxy scale mass (mg*) and',
    print 'stellar mass (sm*) for GALAXY %s' % (galaxy_num)
    print 'C      PARAMETER(mbch1=%.8f, mbch2=%.8f) !central bh mass parameters for Galaxy %s' % \
                                                      (central_bh_mass_coeffs[0], central_bh_mass_coeffs[1], galaxy_num)
    print 'C      PARAMETER(mg1=%.8e,mg2=%.8e,' % (galaxy_mass_coeffs[0], galaxy_mass_coeffs[1])
    print 'C     &          mg3=%.8e,mg4=%.8f,' % (galaxy_mass_coeffs[2], galaxy_mass_coeffs[3])
    print 'C     &          mg5=%.8f,mg6=%.8f,'  % (galaxy_mass_coeffs[4], galaxy_mass_coeffs[5])
    print 'C     &          mg7=%.8f,mg8=%.8f)' % (galaxy_mass_coeffs[6], galaxy_mass_coeffs[7])
    print 'C      PARAMETER(sm1=%.8e, sm2=%.8e,' % (stellar_mass_coefficients[0], stellar_mass_coefficients[1])
    print 'C     &          sm3=%.8e, sm4=%.8e,' % (stellar_mass_coefficients[2], stellar_mass_coefficients[3])
    print 'C     &          sm5=%.8e, sm6=%.8e,' % (stellar_mass_coefficients[4], stellar_mass_coefficients[5])
    print 'C     &          sm7=%.8e, sm8=%.8f)' % (stellar_mass_coefficients[6], stellar_mass_coefficients[7])
    print 'C***********************************************************************************************************'

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
    with open(os.path.join(cluster_folder, 'Infalling_BH_masses_galaxy_%s.txt' % (galaxy_num)), 'rb') as f:
        whole = f.readlines()
    whole = [part.split() for part in whole[1:]]
    whole = np.array([[ast.literal_eval(entry) for entry in line] for line in whole])
    whole.view('i8,i8,i8').sort(order=['f2'], axis=0)
    ts_at_ctr = []
    with open(os.path.join(cluster_folder, 'Infalling_BH_masses_galaxy_%s.txt' % (galaxy_num)), 'wb') as f:
        f.write('BH ID\tInfall mass [Msun]\tInfall time [Myr]\tt at Ctr [Myr]\n')
        for line in whole:
            timo = t0+line[2]
            z = get_z(timo)
            stlr_mass = get_stlr_mass(timo, smc)
            eff_rad = get_eff_rad(stlr_mass, z)
            sigma_eff = get_sigma_eff(stlr_mass, z)
            r_drop = 5000.0
            t_fric = get_t_fric(sigma_eff, line[1], eff_rad, r_drop)
            t_at_center = timo + t_fric*1000.
            f.write('%s\t%11.0f\t\t\t%8.3f\t\t\t%11.0f\n' % (str(int(line[0])).zfill(4), line[1], line[2],
                                                             t_at_center))
            ts_at_ctr.append(t_at_center)
    return ts_at_ctr

def run():
    mpl.rcParams.update({'font.size': 15})
    plt.rc('text', usetex=True)
    plt.rc('font', family='Computer Modern Sans serif')
    smbh_cluster = get_clusters.get_pickled_file('smbh_cluster.pkl')
    galaxies_cluster_no_bad_z = get_clusters.get_pickled_file('galaxies_cluster_no_bad_z.pkl')

    galaxies_by_id = analyze_clusters.get_galaxies_by_id(galaxies_cluster_no_bad_z)
#    print(galaxies_by_id['217'])
    galaxies_masses, final_masses = analyze_clusters.get_galaxies_masses(galaxies_by_id)
#    write_galaxies_final_masses(galaxies_masses)

#    galaxy_num = '65'
    merged_bhs_from_simulations = {'1': [[7.63838E+06,2.35248E+06,1.23739E+09,2.14763E+08],
                                         [3146.E+06, 12560.E+06, 2974.E+06, 71209.E+06]],
                                   '65': [[8.08983E+08,7.20130E+07,2.96940E+08,7.06418E+07],
                                          [10000.E+06, 61228.E+06, 19018.E+06, 63883.E+06]],
                                   '187': [[1.21178E+07, 2.02427E+07], [72505.E+06, 50722.E+06]],
                                   '217': [[1.34011E+07, 1.11187E+06], [3248.E+06, 20166.E+06]]}
    galaxy_names = {'138': 'Ignore Name', '201': 'Ignore Name', '153': 'Ignore Name', '215': 'Ignore Name', '213': 'Ignore Name', '210':\
    'Ignore Name', '217': 'Ignore Name', '195': 'Ignore Name', '158': 'Ignore Name', '197': 'Ignore Name', '242': 'Ignore Name', '191':\
    'Ignore Name', '193': 'Ignore Name', '132': 'Ignore Name', '238': 'Ignore Name', '239': 'Ignore Name', '179': 'Ignore Name', '178':\
    'Ignore Name', '234': 'Ignore Name', '176': 'Ignore Name', '236': 'Ignore Name', '174': 'Ignore Name', '230': 'Ignore Name', '166':\
    'Ignore Name', '232': 'Ignore Name', '170': 'Ignore Name', '92': 'Ignore Name', '223': 'Ignore Name', '211': 'Ignore Name', '222':\
    'Ignore Name', '1': 'Ignore Name', '180': 'Ignore Name', '181': 'Ignore Name', '186': 'Ignore Name', '187': 'Ignore Name', '184':\
    'Ignore Name', '231': 'Ignore Name', '188': 'Ignore Name', '220': 'Ignore Name', '235': 'Ignore Name', '144': 'Ignore Name', '205':\
    'Ignore Name', '207': 'Ignore Name', '209': 'Ignore Name', '208': 'Ignore Name', '149': 'Ignore Name', '244': 'Ignore Name', '168':\
    'Ignore Name', '240': 'Ignore Name', '243': 'Ignore Name', '228': 'Ignore Name', '227': 'Ignore Name', '226': 'Ignore Name', '225':\
    'Ignore Name', '224': 'Ignore Name', '160': 'Ignore Name', '161': 'Ignore Name', '221': 'Ignore Name', '163': 'Ignore Name', '80':\
    'Ignore Name', '114': 'Ignore Name', '126': 'Ignore Name', '51': 'Ignore Name', '233': 'Ignore Name', '241': 'Ignore Name', '237':\
    'Ignore Name', '65': 'Ignore Name'}
    stellar_data = {}
    for galaxy_num, name in galaxy_names.iteritems():
        '''
        t_s, masses_fitted, galaxy_mass_coeffs = get_galaxy_mass_data(galaxies_masses, galaxy_num)
#        print(name, ':', str(galaxies_masses[galaxy_num][2][0]), str(galaxies_masses[galaxy_num][2][-1]))
        central_bh_mass, central_bh_mass_fitted, central_bh_mass_coeffs, bh_masses = \
                                                get_central_bh_mass_growth(galaxies_cluster_no_bad_z, smbh_cluster,
                                                                                                             galaxy_num)
        print(''.join([name, ' bh mass: ']), central_bh_mass[0], central_bh_mass[-1])
        '''
        time, stellar_mass, stellar_mass_coefficients, stellar_mass_fitted = get_stellar_mass_growth(galaxies_by_id,
                                                                                                             galaxy_num)
        stellar_data[galaxy_num] = stellar_mass_coefficients

        '''
        output_coeffs_for_achain_h(galaxy_num, central_bh_mass_coeffs, galaxy_mass_coeffs, stellar_mass_coefficients)
#        plot_all_data(t_s, masses_fitted, galaxies_masses, galaxy_num, central_bh_mass, central_bh_mass_fitted, time,
#                                                                                     stellar_mass, stellar_mass_fitted)
        ts_at_ctr = read_orbiting_bhs(galaxy_num, stellar_mass_coefficients)
        ts_at_ctr = np.array(ts_at_ctr)
#        print ts_at_ctr
        plt.figure()
        plt.hist(np.log10(ts_at_ctr*1.e6), bins=np.arange(9., 18.+.25, .25), alpha=0.5, color='blue',
                                                                                                   label='Orbiting BHs')
        plt.hist(np.log10(merged_bhs_from_simulations[galaxy_num][1]), bins=np.arange(9., 18.+.25, .25), alpha=0.5,
                                                                                       color='red', label='Merged BHs')
        plt.axvline(np.log10(hubble_time), color='g', linestyle='dashed', linewidth=2, label='Hubble Time')
        plt.xlabel(r'$\log_{10}\left(t_{fric}\right) [Yrs]$')
        plt.ylim([0,23])
        plt.yticks(np.arange(0, 23, 2))
        plt.ylabel('Count of black holes')
#        plt.title(''.join(['Histogram of Orbiting Black Holes Time to Reach Center for Galaxy ',
#                           name]))
        plt.legend(loc='upper left')
        plt.savefig(os.path.join(plots_folder, ''.join(['t_at_center_histogram_gal_', galaxy_num, '.png'])))
        plt.close()

        plt.figure()
        bh_masses = np.array(bh_masses)
#        print(galaxy_num, np.log10(bh_masses))
        plt.hist(np.log10(bh_masses), bins=np.arange(2., 10.+.25, .25), alpha=0.5, color='blue', label='Orbiting BHs')
        plt.hist(np.log10(merged_bhs_from_simulations[galaxy_num][0]), bins=np.arange(2., 10.+.25, .25), alpha=0.5,
                                                                                       color='red', label='Merged BHs')
#        plt.xlabel('BH Mass (M$_\odot$), log')
        plt.xlabel(r'$\log_{10}\left(M_{BH}\right) [M_\odot]$')
        plt.ylim([0,22])
        plt.yticks(np.arange(0, 23, 2))
        plt.ylabel('Count of black holes')
        plt.legend(loc='upper right')
#        plt.title(''.join(['Histogram of Masses of Orbiting Black Holes for Galaxy ',
#                           name]))
        plt.savefig(os.path.join(plots_folder, ''.join(['orbiting_bh_mass_histogram_gal_', galaxy_num, '.png'])))
        plt.close()
        '''
    print(stellar_data)

if __name__ == '__main__':
    run()
