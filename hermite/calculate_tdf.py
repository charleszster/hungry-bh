#! /home/charles/anaconda/bin/python
##! /Users/chaz/anaconda/bin/python	#Mac

import numpy as np
import matplotlib.pyplot as plt
import collections

#Stellar mass coefficients for each galaxy
g1 = np.array([5.55238047e-27, -3.16757613e-22, 7.51108077e-18, -9.59173067e-14, 7.13616762e-10,
                         -3.10946120e-06, 7.48070688e-03, 4.76681897])
g51 = np.array([9.21720403e-28, -5.85756647e-23, 1.59023638e-18, -2.40208805e-14, 2.19815322e-10, -1.23870657e-06,
                4.11467642e-03, 6.60341066])
g65 = np.array([1.40101390e-27, -7.68778493e-23, 1.72469599e-18, -2.05188123e-14, 1.41935495e-10, -6.05073236e-07,
                1.72244648e-03, 9.60556175])


def stellar_mass(t):
    return 10.**(g1[0]*t**7. + g1[1]*t**6. + g1[2]*t**5. + g1[3]*t**4. + g1[4]*t**3. + g1[5]*t**2. + g1[6]*t + g1[7])

def z(t):
    return 7.20196192*(t/1000.)**(-0.59331986) - 1.52145449

t0 = 1.5653765e+03

with open('InfallBHmasses_gal_1_biggest.txt') as f:
    lines = f.readlines()

lines.pop(0)
lines = [[el for el in line.rstrip().split('\t') if el is not ''] for line in lines]
lines_dict = []
for line in lines:
    lines_dict.append(dict({line[0]: [float(line[1]),float(line[2])]}))

for line in lines_dict:
    infall_time = line.values()[0][0]*1000.
    bh_mass = line.values()[0][1]
    t = t0+infall_time
    Re = 2500. * (stellar_mass(t)/1.e+11)**0.73 * (1.+z(t))**(-0.98)
    sigma = 190. * (stellar_mass(t)/1.e+11)**0.2 * (1.+z(t))**0.47
    t_df = (19./np.log(1.+stellar_mass(t)/bh_mass)) * (Re/5000.)**2. * (sigma/200.) * (1.e8/bh_mass) #* 0.5 #Gyr
    print '%s, %.2f, %.2f, %.3e, %.4e, %.2f, %.2f' % (line.keys()[0], z(t), t, bh_mass, Re, sigma, t_df)

'''
lines_dict = {}
for line in lines:
    lines_dict[line[0]] =  [float(line[1]),float(line[2])]
for id, value in lines_dict.iteritems():
    infall_time = value[0]*1000.
    bh_mass = value[1]
    t = t0+infall_time
    Re = 2500. * (stellar_mass(t)/1.e+11)**0.73 * (1.+z(t))**(-0.98)
    sigma = 190. * (stellar_mass(t)/1.e+11)**0.2 * (1.+z(t))**0.47
    t_df = (19./np.log(1.+stellar_mass(t)/bh_mass)) * (Re/5000.)**2 * (sigma/200.) * (1.e8/bh_mass) * 0.5 #Gyr
    print '%s, %.2f, %.2f, %.3e, %.4e, %.2f, %.2f' % (id, z(t), t, bh_mass, Re, sigma, t_df)
'''
