#! /Users/chaz/anaconda/bin/python	#Mac
##! /home/charles/anaconda/bin/python	#Linux

import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import sigma_equations as sig_eq

def test_dsigs_dr(r, rh, rc, G, Mtot, sigma_near, sigma_far):
    dsignear_dr = sig_eq.get_dsignear_dr(rh, r, rc, G, Mtot)
    dsigfar_dr =  sig_eq.get_dsigfar_dr(rh, r, rc, G, Mtot)
    f, (ax1, ax2) = plt.subplots(2, sharex=True)
    ax1.loglog(r, sigma_near, 'b', label='sigma_near')
    ax1.loglog(r, sigma_far, 'r', label='sigma_far')
    ax1.set_xlim([1.e-3, 6])
    ax1.set_xlabel('r (pc)')
    ax1.set_ylabel('sigma (km/s)')
    ax1.set_ylim([8., 25.])
    ax1.set_title('sigma near and sigma far')
    ax1.legend(loc='best')
    ax2.plot(r, dsignear_dr, 'b', label='dsignear_dr')
    ax2.plot(r[r>=1.e-2], dsigfar_dr[r>=1.e-2], 'r', label='dsigfar_dr')    
    ax2.set_title('dsigma/dr (near and far)')
    ax2.set_xlabel('r (pc)')
    ax2.set_ylabel('dsigma/dr (pc^-1*km/s)')
    ax2.legend(loc='best')
    plt.show()

def test_dsigs_drh(r, rh, rc, G, Mtot, sigma_near, sigma_far):
    dsignear_drh = sig_eq.get_dsignear_drh(rh, r, rc, G, Mtot)
    dsigfar_drh =  sig_eq.get_dsigfar_drh(rh, r, rc, G, Mtot)
    f, (ax1, ax2) = plt.subplots(2, sharex=True)
    ax1.plot(rh, sigma_near, 'b', label='sigma_near')
    ax1.plot(rh, sigma_far, 'r', label='sigma_far')
#    ax1.set_xlim([1.e-3, 6])
    ax1.set_xlabel('rh (pc)')
    ax1.set_ylabel('sigma (km/s)')
#    ax1.set_ylim([8., 25.])
    ax1.set_title('sigma near and sigma far')
    ax1.legend(loc='best')
    ax2.plot(rh, dsignear_drh, 'b', label='dsignear_drh')
    ax2.plot(rh, dsigfar_drh, 'r', label='dsigfar_drh')
    ax2.set_title('dsigma/drh (near and far)')
    ax2.set_xlabel('rh (pc)')
    ax2.set_ylabel('dsigma/drh (pc^-1*km/s)')
    ax2.legend(loc='best')
    plt.show()

def run():
    '''
    G = sp.Symbol('G')
    Mtot = sp.Symbol('Mtot')
    rh = sp.Symbol('rh')
    rc = sp.Symbol('rc')
    r = sp.Symbol('r')
    pie = sp.Symbol('pi')

    arktan = sp.atan
    lahg = sp.log
    squirt = sp.sqrt
    init_guess = 0.
    dsig_drh = sp.diff(sig_eq.get_sigma_near(rh, r, rc, G, Mtot, init_guess, arktan, lahg, squirt, pie), r)
    print dsig_drh
    '''
    arktan = np.arctan
    lahg = np.log
    squirt = np.sqrt
    pie = np.pi
    init_guess = 0.
    G = 0.0043009211
    Mtot = 1.e5
    rc = .002

    rh = 1.
    r = np.linspace(1.e-3, 10., 100000)
    sigma_near = sig_eq.get_sigma_near(rh, r, rc, G, Mtot, init_guess, arktan, lahg, squirt, pie)
    sigma_far =  sig_eq.get_sigma_far(rh, r, rc, G, Mtot, init_guess, arktan, lahg, squirt, pie)
    test_dsigs_dr(r, rh, rc, G, Mtot, sigma_near, sigma_far)

    rh = np.linspace(1., 10., 100000)
    r = 1.
    sigma_near = sig_eq.get_sigma_near(rh, r, rc, G, Mtot, init_guess, arktan, lahg, squirt, pie)
    sigma_far =  sig_eq.get_sigma_far(rh, r, rc, G, Mtot, init_guess, arktan, lahg, squirt, pie)
    test_dsigs_drh(r, rh, rc, G, Mtot, sigma_near, sigma_far)

if __name__ == '__main__':
    run()