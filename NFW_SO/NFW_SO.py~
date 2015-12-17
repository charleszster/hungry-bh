#!/Users/chaz/anaconda/bin/python #This is the python interpreter I use on Mac
##!/home/charles/anaconda/bin/python 
# -*- coding: utf-8 -*-

import time
import numpy as np
import matplotlib.pyplot as plt
import pickle
import os
import scipy.optimize
import sys
sys.path.append('/Users/chaz/Dropbox/Columbia/Ostriker/Cannibalism/mergertree/scripts')
import sympy as sp
from sympy.solvers import solve
from sympy import Symbol
import Constants
import galaxy_mass_coefficients as gmc
from calculate_r200 import get_H, get_r200

plots_folder = Constants.plots_folder
G = Constants.G
delta_vir = Constants.delta_vir
t0 = Constants.t0 #Myr
h = Constants.h
mg = gmc.mg

def get_z(t):
    return 7.20196192*(t/1000.)**(-0.59331986) - 1.52145449

def galaxy_mass(mg, t):
    return 10.**(mg[0]*(t/1000.)**7. + mg[1]*(t/1000.)**6. + mg[2]*(t/1000.)**5. + mg[3]*(t/1000.)**4. \
            + mg[4]*(t/1000.)**3. + mg[5]*(t/1000.)**2. + mg[6]*(t/1000.) + mg[7])

def get_c(m200):
    return 10.**(1.025-0.097*np.log10(m200/((10.**12)/h)))	#Equation 7 from NFW paper

def get_RS(r200, c):
    return r200/c	#From paragraph below Equation 3 in NFW paper

def get_rho_crit(H, m200, r200):
    return 3.*(H**2.)/(8.*np.pi*G)	#From paragraph above Equation 1 in NFW paper
#    return m200*(3./4.)/(np.pi*delta_vir*r200**3)	#Equation 1 in NFW paper

def get_rho_0(rho_crit, c):
	return rho_crit*(delta_vir*c**3.)/(3.*(np.log(1.+c)-c/(1.+c))) #Equation 4 in NFW paper

def get_rho_NFW(R, rho_0, R_S):
	return rho_0/((R/R_S)*(1.+R/R_S)**2.)	#Equation 2 in NFW paper

def get_R_half_mass(c, r200=None, m200=None, rho_0=None, R_S=None, optimize=False):
    if optimize:
        eqn = lambda r_rs: 4.*np.pi*rho_0*(R_S**3.)*(np.log(1.+r_rs)-r_rs/(1.+r_rs)) - 0.5*m200
        init_guess = c
        sol = scipy.optimize.fsolve(eqn, init_guess)[0]
        R_half_mass = sol*R_S
        return R_half_mass
    else:
        return r200*(0.6082 - 0.1843*np.log10(c) - 0.1011*(np.log10(c)**2.) + 0.03918*(np.log10(c)**3.))

def get_NFW_parameters(t, H, gal):
    m200 = galaxy_mass(mg[gal], t)
    r200 = get_r200(m200, H)
    c = get_c(m200)
    R_S = get_RS(r200, c)
    rho_crit = get_rho_crit(H, m200, r200)
    rho_0 = get_rho_0(rho_crit, c)
    return (m200, r200, c, R_S, rho_crit, rho_0)

def get_R_H(m200, R_C, R_half_mass):
#    rho_c = m200*(R_H + R_C)/(2*np.pi**2*R_C**2*R_C**2)
    eqn = lambda R_H: (4.*np.pi*(R_C**2.)*(R_H**2.) * \
                      (m200*(R_H + R_C) / (2.*(np.pi**2.)*(R_C**2.)*(R_H**2.))) / ((R_H**2.)-(R_C**2.))) * \
                       (R_H*np.arctan(R_half_mass/R_H) - R_C*np.arctan(R_half_mass/R_C)) - 0.5*m200
    init_guess_R_H = R_half_mass
    sol = scipy.optimize.fsolve(eqn, init_guess_R_H)[0]
    return sol

def R_half_mass_r200_func(x, a, b, c, d):
    return a*(x**3.) + b*(x**2.) + c*x + d

def get_rho_c(m200, R_C, R_H):
    return m200*(R_H+R_C)/(2.*np.pi**2.*R_C**2.*R_H**2.)

def plot_c_vs_t(t_s, c_array):
    plt.figure()
    plt.semilogx(t_s, c_array)
    plt.xlabel('Time (Myr)')
    plt.ylabel('c')
    plt.title('Mass concentration (c) throughout simulation')
    plt.minorticks_on()
    plt.savefig(''.join([plots_folder, 'c.png']))
    plt.close()

def getYr1r1(r, r1):
    return -(4.*r1*(np.pi*r1**2. + r*r1 - 3.*np.pi*r**2. - 2.*(r1**2. - 3.*r**2.)*np.arctan(r1/r)) - \
             16.*r**3.*np.log(1. + r1**2./r**2.) + 3.*r**3.*(np.pi**2. - 4.*np.arctan(r/r1)**2.))/(24.*r**3.*r1**2.)
#    return -(4.*np.pi*(r1**3.) + 4.*r*(r1**2.) - 12.*np.pi*(r**2.)*r1 + 3.*(np.pi**2.)*(r**3.) + (24.*(r**2.)*r1 - \
#             8.*(r1**3.))*np.arctan(r1/r) - 16.*(r**3.)*np.log(1.+(r1/r)**2.) - \
#             12.*(r**3.)*np.arctan(r/r1)**2.)/(24.*(r**3.)*(r1**2.))

def getYrhrc(r, rh, rc):
    return 1./(2.*r**2.) + np.log(r**2./(r**2. + rc**2.))/(2.*rc**2.) - np.log(1. + rc**2./r**2.)/(6.*rh**2.)

def getYrcrh(r, rh, rc):
    return (np.pi*rc + r - 2.*rc*np.arctan(rc/r))/(6.*r**3.) - np.log(1. + rc**2./r**2.)/(6.*rc**2.) - \
            np.pi*rc/(2.*r*rh**2.) + rc*np.arctan(rc/r)/(r*rh**2.) - np.log(1. + rc**2./r**2.)/(2.*rh**2.)
#    return (np.pi*rc + r - 3.*np.pi*(r**2.)*rc/(rh**2.) + 2.*rc*(np.arctan(rc/r)*(3.*(r**2.)/(rh**2.) - 1.))) / \
#            (6.*(r**3.)) - (1./(6.*rc**2.) + 1./(2.*rh**2.))*np.log(1. + (rc**2.)/(r**2.))

def get_sigma_near(r, rh, rc, rho_c):
    vr2so = -4.*np.pi*G*rho_c*rc**2.*rh**2.*(r**2. + rh**2.)*(r**2. + rc**2.)*(getYrhrc(r, rh, rc) + getYr1r1(r, rc) +\
             getYr1r1(r, rh) + getYrcrh(r, rh, rc))/(rc**2. - rh**2.)**2.
    return (3*vr2so)**0.5

def get_sigma_far(r, rh, rc, m200):
    return (6.*G*m200*(r**2. + rh**2.)*(r**2. + rc**2.)*(np.pi*(rh-rc))**(-1.)*(np.pi**2./(8.*rh**4.) + \
             np.pi/(6.*rh*r**3.) + 1./(6.*rh**2.*r**2.) - np.pi/(2.*rh**3.*r) - np.arctan(r/rh)**2./(2.*rh**4.) - \
             np.arctan(rh/r)/(3.*rh*r**3.) + np.arctan(rh/r)/(rh**3.*r) - 2.*np.log(1.+rh**2./r**2.)/(3.*rh**4.) - \
             rc*np.pi*(rh**3. - 3.*rh*r**2. + 3.*r**3.*np.arctan(rh/r))/(6.*rh**5.*r**3.)))**0.5

def run():
    tend = 13.722672e3 #Myr
    gal = '1'
    r200_list = []
    m200_list = []
    c_list = []
    R_half_mass_list = []
    R_H_list = []
    rho_c_list = []
    t_s = np.arange(t0, tend, 100.)
    for t in t_s: #increment is 100 Myr
        z = get_z(t)
        H = get_H(z)
        NFW_params = get_NFW_parameters(t, H, gal)
        m200, r200, c, R_S, rho_crit, rho_0 = NFW_params
        r200_list.append(r200)
        m200_list.append(m200)
        c_list.append(c)
#        R = np.linspace(0.01, r200, num=100000)
#        rho_NFW = get_rho_NFW(R, rho_0, R_S)
        R_half_mass = get_R_half_mass(c, r200=r200)
        R_half_mass_list.append(R_half_mass)
        R_H = R_half_mass
        R_H_list.append(R_H)
        R_C = Constants.R_C
        rho_c = get_rho_c(m200, R_C, R_H)
        rho_c_list.append(rho_c)

    r200_array = np.array(r200_list)
    m200_array = np.array(m200_list)
    c_array = np.array(c_list)
    R_half_mass_array = np.array(R_half_mass_list)
    R_H_array = np.array(R_H_list)
    rho_c_array = np.array(rho_c_list)

#    plot_c_vs_t(t_s, c_array)
    '''
    plt.figure()
    plt.loglog(t_s, rho_c_array)
    plt.xlabel('Log time (Myr)')
    plt.ylabel('Log rho_c of OS profile (Msun/pc**3)')
    plt.savefig(''.join([plots_folder, 'rho_c.png']))
    plt.close()
    '''
    R_H_index = 0
#    R = np.linspace(0.1, R_H_array[R_H_index], num=10000)
    '''
    rho_rhoc = ((1+R**2/R_C**2)*(1+R**2/R_H_array[R_H_index]**2))**(-1)
    plt.figure()
    plt.loglog(R, rho_rhoc)
    plt.axvline(x=R_H_array[R_H_index], color='r')
    plt.xlabel('Log R')
    plt.ylabel('Log rho/rho_c')
    plt.title('rho/rho_c as function of R for r_h=%.f' % (R_H_array[R_H_index]))
    plt.savefig(''.join([plots_folder, 'rho_rhoc_as_func_of_r.png']))
    '''
    m200_index = R_H_index
    rho_c_index = R_H_index
##############################
#  TEMPORARY, JUST FOR NICK'S TEST CASE FROM THE PAPER
    R_H = 1.
    m200 = 1.e5
    R_C = [.002, .01, .05, .1]
#    colors = ['k', 'b', 'm', 'r']
    R = np.linspace(1.e-3, 100., num=100000)
    R_near = [R[R<=0.065*1.15], R[R<=0.12*1.15], R[R<=0.2*1.15], R[R<=0.25*1.15]]
    R_far = [R[R>=0.075*.85], R[R>=0.15*.85], R[R>=0.25*.85], R[R>=0.3*.85]]
###############################
    fig, ax = plt.subplots(2, 2)
    for i, a in enumerate(ax.flatten()):
        rho_c = get_rho_c(m200, R_C[i], R_H)
        sigma_near = get_sigma_near(R, R_H, R_C[i], rho_c)
##############################
#    sigma_near = get_sigma_near(R, R_H_array[R_H_index], R_C, rho_c_array[rho_c_index])
#        a.loglog(R_near[i], sigma_near[:len(R_near[i])], c='k', label='sigma near')
        a.loglog(R, sigma_near, c='k', label='sigma near')
##############################
#  TEMPORARY, JUST FOR NICK'S TEST CASE FROM THE PAPER
        sigma_far = get_sigma_far(R, R_H, R_C[i], m200)
##############################
#    sigma_far = get_sigma_far(R, R_H_array[R_H_index], R_C, m200_array[m200_index])
#        a.loglog(R_far[i], sigma_far[len(sigma_far)-len(R_far[i]):], c='r', label='sigma far')
        a.loglog(R, sigma_far, c='r', label='sigma far')
##############################
#  TEMPORARY, JUST FOR NICK'S TEST CASE FROM THE PAPER
        a.axvline(x=R_C[i], color='b', linestyle='--')
        a.set_ylim([8., 25.])
        a.set_xlim([.001, 10.])
##############################
        a.set_xlabel('Log r[pc]')
        a.set_ylabel('Log sigma[km/s]')
#    plt.legend(loc='best')
        a.grid(b=True, which='major', color='r', linestyle='-')
##############################
#  TEMPORARY, JUST FOR NICK'S TEST CASE FROM THE PAPER
        a.set_title('r_c=%.3f' % (R_C[i]), fontsize=10)
##############################
#    plt.title('sigma as function of R for r_h=%.f' % (R_H_array[R_H_index]))
    fig.suptitle(''.join(['sigma as f\'n of r for r_h=%.1f, m200=%.f' % (R_H, m200)]), fontsize=12, y=1.)
    fig.tight_layout()
    fig.savefig(''.join([plots_folder, 'sigma_as_func_of_r.png']))
    
    
if __name__ == '__main__':
    run()
