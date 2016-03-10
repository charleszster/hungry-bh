#! /home/charles/anaconda/bin/python	
##! /Users/chaz/anaconda/bin/python	#Mac

import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

def Yr1r1(arktan, lahg, pie, r, r1):
    return -(4.*r1*(pie*r1**2. + r*r1 - 3.*pie*r**2. - 2.*(r1**2. - 3.*r**2.)*arktan(r1/r)) -\
            16.*r**3.*lahg(1. + r1**2./r**2.) + 3.*r**3.*(pie**2. - 4.*arktan(r/r1)**2.)) /\
            (24.*r**3.*r1**2.)

def Yrhrc(lahg, r, rh, rc):
    return 1./(2.*r**2.) + lahg(r**2./(r**2. + rc**2.))/(2.*rc**2.) - lahg(1. + rc**2./r**2.)/(6.*rh**2.)

def Yrcrh(arktan, lahg, pie, r, rh, rc):
    return (pie*rc + r - 2.*rc*arktan(rc/r))/(6.*r**3.) - lahg(1. + rc**2./r**2.)/(6.*rc**2.) -\
            pie*rc/(2.*r*rh**2.) + rc*arktan(rc/r)/(r*rh**2.) - lahg(1. + rc**2./r**2.)/(2.*rh**2.)

def get_sigma_near(rh, r, rc, G, Mtot, init_guess=0., arktan=None, lahg=None, squirt=None, pie=None):
    rho_c = Mtot*(rh+rc)/(2*pie**2*rc**2*rh**2)
    return squirt(-12.*pie*G*rho_c*rc**2.*rh**2.*(r**2. + rh**2.)*(r**2. + rc**2.)*(Yrhrc(lahg, r, rh, rc) +\
                   Yr1r1(arktan, lahg, pie, r, rc) + Yr1r1(arktan, lahg, pie, r, rh) + Yrcrh(arktan, lahg, pie, r, rh,\
                                                                                        rc))/(rc**2. - rh**2.)**2.) -\
                                                                                        init_guess

def get_dsignear_drh(rh, r, rc, G, Mtot, init_guess=0., arktan=None, lahg=None, squirt=None, pie=None):
    return -2.44948974278318*np.pi*np.sqrt(-G*Mtot*(r**2.0 + rc**2.0)*(r**2.0 + rh**2.0)*(rc + rh)*(rc**2.0 - rh**2.0)**(-2.0)*(-0.5*np.pi*rc*rh**(-2.0)/r + 0.0416666666666667*r**(-3.0)*rc**(-2.0)*(-3.0*r**3.0*(np.pi**2.0 - 4.0*np.arctan(r/rc)**2.0) + 16.0*r**3.0*np.log(r**(-2.0)*rc**2.0 + 1.0) - 4.0*rc*(-3.0*np.pi*r**2.0 + np.pi*rc**2.0 + r*rc - (-6.0*r**2.0 + 2.0*rc**2.0)*np.arctan(rc/r))) + 0.0416666666666667*r**(-3.0)*rh**(-2.0)*(-3.0*r**3.0*(np.pi**2.0 - 4.0*np.arctan(r/rh)**2.0) + 16.0*r**3.0*np.log(r**(-2.0)*rh**2.0 + 1.0) - 4.0*rh*(-3.0*np.pi*r**2.0 + np.pi*rh**2.0 + r*rh - (-6.0*r**2.0 + 2.0*rh**2.0)*np.arctan(rh/r))) + 0.166666666666667*r**(-3.0)*(np.pi*rc + r - 2.0*rc*np.arctan(rc/r)) + 0.5*r**(-2.0) + 0.5*rc**(-2.0)*np.log(r**2.0/(r**2.0 + rc**2.0)) - 0.166666666666667*rc**(-2.0)*np.log(r**(-2.0)*rc**2.0 + 1.0) - 0.666666666666667*rh**(-2.0)*np.log(r**(-2.0)*rc**2.0 + 1.0) + rc*rh**(-2.0)*np.arctan(rc/r)/r)/np.pi)*(rc**2.0 - rh**2.0)**2.0*(-2.0*G*Mtot*rh**1.0*(r**2.0 + rc**2.0)*(r**2.0 + rh**2.0)*(rc + rh)*(rc**2.0 - rh**2.0)**(-3.0)*(-0.5*np.pi*rc*rh**(-2.0)/r + 0.0416666666666667*r**(-3.0)*rc**(-2.0)*(-3.0*r**3.0*(np.pi**2.0 - 4.0*np.arctan(r/rc)**2.0) + 16.0*r**3.0*np.log(r**(-2.0)*rc**2.0 + 1.0) - 4.0*rc*(-3.0*np.pi*r**2.0 + np.pi*rc**2.0 + r*rc - (-6.0*r**2.0 + 2.0*rc**2.0)*np.arctan(rc/r))) + 0.0416666666666667*r**(-3.0)*rh**(-2.0)*(-3.0*r**3.0*(np.pi**2.0 - 4.0*np.arctan(r/rh)**2.0) + 16.0*r**3.0*np.log(r**(-2.0)*rh**2.0 + 1.0) - 4.0*rh*(-3.0*np.pi*r**2.0 + np.pi*rh**2.0 + r*rh - (-6.0*r**2.0 + 2.0*rh**2.0)*np.arctan(rh/r))) + 0.166666666666667*r**(-3.0)*(np.pi*rc + r - 2.0*rc*np.arctan(rc/r)) + 0.5*r**(-2.0) + 0.5*rc**(-2.0)*np.log(r**2.0/(r**2.0 + rc**2.0)) - 0.166666666666667*rc**(-2.0)*np.log(r**(-2.0)*rc**2.0 + 1.0) - 0.666666666666667*rh**(-2.0)*np.log(r**(-2.0)*rc**2.0 + 1.0) + rc*rh**(-2.0)*np.arctan(rc/r)/r)/np.pi - 1.0*G*Mtot*rh**1.0*(r**2.0 + rc**2.0)*(rc + rh)*(rc**2.0 - rh**2.0)**(-2.0)*(-0.5*np.pi*rc*rh**(-2.0)/r + 0.0416666666666667*r**(-3.0)*rc**(-2.0)*(-3.0*r**3.0*(np.pi**2.0 - 4.0*np.arctan(r/rc)**2.0) + 16.0*r**3.0*np.log(r**(-2.0)*rc**2.0 + 1.0) - 4.0*rc*(-3.0*np.pi*r**2.0 + np.pi*rc**2.0 + r*rc - (-6.0*r**2.0 + 2.0*rc**2.0)*np.arctan(rc/r))) + 0.0416666666666667*r**(-3.0)*rh**(-2.0)*(-3.0*r**3.0*(np.pi**2.0 - 4.0*np.arctan(r/rh)**2.0) + 16.0*r**3.0*np.log(r**(-2.0)*rh**2.0 + 1.0) - 4.0*rh*(-3.0*np.pi*r**2.0 + np.pi*rh**2.0 + r*rh - (-6.0*r**2.0 + 2.0*rh**2.0)*np.arctan(rh/r))) + 0.166666666666667*r**(-3.0)*(np.pi*rc + r - 2.0*rc*np.arctan(rc/r)) + 0.5*r**(-2.0) + 0.5*rc**(-2.0)*np.log(r**2.0/(r**2.0 + rc**2.0)) - 0.166666666666667*rc**(-2.0)*np.log(r**(-2.0)*rc**2.0 + 1.0) - 0.666666666666667*rh**(-2.0)*np.log(r**(-2.0)*rc**2.0 + 1.0) + rc*rh**(-2.0)*np.arctan(rc/r)/r)/np.pi - G*Mtot*(r**2.0 + rc**2.0)*(r**2.0 + rh**2.0)*(rc + rh)*(rc**2.0 - rh**2.0)**(-2.0)*(1.0*np.pi*rc*rh**(-3.0)/r - 0.0833333333333333*r**(-3.0)*rh**(-3.0)*(-3.0*r**3.0*(np.pi**2.0 - 4.0*np.arctan(r/rh)**2.0) + 16.0*r**3.0*np.log(r**(-2.0)*rh**2.0 + 1.0) - 4.0*rh*(-3.0*np.pi*r**2.0 + np.pi*rh**2.0 + r*rh - (-6.0*r**2.0 + 2.0*rh**2.0)*np.arctan(rh/r))) + 0.0416666666666667*r**(-3.0)*rh**(-2.0)*(12.0*np.pi*r**2.0 - 4.0*np.pi*rh**2.0 + 32.0*r**1.0*rh**1.0/(r**(-2.0)*rh**2.0 + 1.0) - 4.0*r*rh - 24.0*r**4.0*np.arctan(r/rh)**1.0/(rh**2*(r**2/rh**2 + 1)) - 4.0*rh*(2.0*np.pi*rh**1.0 + r - 4.0*rh**1.0*np.arctan(rh/r) - (-6.0*r**2.0 + 2.0*rh**2.0)/(r*(1 + rh**2/r**2))) + 4.0*(-6.0*r**2.0 + 2.0*rh**2.0)*np.arctan(rh/r)) + 1.33333333333333*rh**(-3.0)*np.log(r**(-2.0)*rc**2.0 + 1.0) - 2.0*rc*rh**(-3.0)*np.arctan(rc/r)/r)/(2*np.pi) - G*Mtot*(r**2.0 + rc**2.0)*(r**2.0 + rh**2.0)*(rc**2.0 - rh**2.0)**(-2.0)*(-0.5*np.pi*rc*rh**(-2.0)/r + 0.0416666666666667*r**(-3.0)*rc**(-2.0)*(-3.0*r**3.0*(np.pi**2.0 - 4.0*np.arctan(r/rc)**2.0) + 16.0*r**3.0*np.log(r**(-2.0)*rc**2.0 + 1.0) - 4.0*rc*(-3.0*np.pi*r**2.0 + np.pi*rc**2.0 + r*rc - (-6.0*r**2.0 + 2.0*rc**2.0)*np.arctan(rc/r))) + 0.0416666666666667*r**(-3.0)*rh**(-2.0)*(-3.0*r**3.0*(np.pi**2.0 - 4.0*np.arctan(r/rh)**2.0) + 16.0*r**3.0*np.log(r**(-2.0)*rh**2.0 + 1.0) - 4.0*rh*(-3.0*np.pi*r**2.0 + np.pi*rh**2.0 + r*rh - (-6.0*r**2.0 + 2.0*rh**2.0)*np.arctan(rh/r))) + 0.166666666666667*r**(-3.0)*(np.pi*rc + r - 2.0*rc*np.arctan(rc/r)) + 0.5*r**(-2.0) + 0.5*rc**(-2.0)*np.log(r**2.0/(r**2.0 + rc**2.0)) - 0.166666666666667*rc**(-2.0)*np.log(r**(-2.0)*rc**2.0 + 1.0) - 0.666666666666667*rh**(-2.0)*np.log(r**(-2.0)*rc**2.0 + 1.0) + rc*rh**(-2.0)*np.arctan(rc/r)/r)/(2*np.pi))/(G*Mtot*(r**2.0 + rc**2.0)*(r**2.0 + rh**2.0)*(rc + rh)*(-0.5*np.pi*rc*rh**(-2.0)/r + 0.0416666666666667*r**(-3.0)*rc**(-2.0)*(-3.0*r**3.0*(np.pi**2.0 - 4.0*np.arctan(r/rc)**2.0) + 16.0*r**3.0*np.log(r**(-2.0)*rc**2.0 + 1.0) - 4.0*rc*(-3.0*np.pi*r**2.0 + np.pi*rc**2.0 + r*rc - (-6.0*r**2.0 + 2.0*rc**2.0)*np.arctan(rc/r))) + 0.0416666666666667*r**(-3.0)*rh**(-2.0)*(-3.0*r**3.0*(np.pi**2.0 - 4.0*np.arctan(r/rh)**2.0) + 16.0*r**3.0*np.log(r**(-2.0)*rh**2.0 + 1.0) - 4.0*rh*(-3.0*np.pi*r**2.0 + np.pi*rh**2.0 + r*rh - (-6.0*r**2.0 + 2.0*rh**2.0)*np.arctan(rh/r))) + 0.166666666666667*r**(-3.0)*(np.pi*rc + r - 2.0*rc*np.arctan(rc/r)) + 0.5*r**(-2.0) + 0.5*rc**(-2.0)*np.log(r**2.0/(r**2.0 + rc**2.0)) - 0.166666666666667*rc**(-2.0)*np.log(r**(-2.0)*rc**2.0 + 1.0) - 0.666666666666667*rh**(-2.0)*np.log(r**(-2.0)*rc**2.0 + 1.0) + rc*rh**(-2.0)*np.arctan(rc/r)/r))

def get_dsignear_dr(rh, r, rc, G, Mtot, init_guess=0., arktan=None, lahg=None, squirt=None, pie=None):
    return -2.44948974278318*np.pi*np.sqrt(-G*Mtot*(r**2.0 + rc**2.0)*(r**2.0 + rh**2.0)*(rc + rh)*(rc**2.0 - rh**2.0)**(-2.0)*(-0.5*np.pi*rc*rh**(-2.0)/r + 0.0416666666666667*r**(-3.0)*rc**(-2.0)*(-3.0*r**3.0*(np.pi**2.0 - 4.0*np.arctan(r/rc)**2.0) + 16.0*r**3.0*np.log(r**(-2.0)*rc**2.0 + 1.0) - 4.0*rc*(-3.0*np.pi*r**2.0 + np.pi*rc**2.0 + r*rc - (-6.0*r**2.0 + 2.0*rc**2.0)*np.arctan(rc/r))) + 0.0416666666666667*r**(-3.0)*rh**(-2.0)*(-3.0*r**3.0*(np.pi**2.0 - 4.0*np.arctan(r/rh)**2.0) + 16.0*r**3.0*np.log(r**(-2.0)*rh**2.0 + 1.0) - 4.0*rh*(-3.0*np.pi*r**2.0 + np.pi*rh**2.0 + r*rh - (-6.0*r**2.0 + 2.0*rh**2.0)*np.arctan(rh/r))) + 0.166666666666667*r**(-3.0)*(np.pi*rc + r - 2.0*rc*np.arctan(rc/r)) + 0.5*r**(-2.0) + 0.5*rc**(-2.0)*np.log(r**2.0/(r**2.0 + rc**2.0)) - 0.166666666666667*rc**(-2.0)*np.log(r**(-2.0)*rc**2.0 + 1.0) - 0.666666666666667*rh**(-2.0)*np.log(r**(-2.0)*rc**2.0 + 1.0) + rc*rh**(-2.0)*np.arctan(rc/r)/r)/np.pi)*(rc**2.0 - rh**2.0)**2.0*(-1.0*G*Mtot*r**1.0*(r**2.0 + rc**2.0)*(rc + rh)*(rc**2.0 - rh**2.0)**(-2.0)*(-0.5*np.pi*rc*rh**(-2.0)/r + 0.0416666666666667*r**(-3.0)*rc**(-2.0)*(-3.0*r**3.0*(np.pi**2.0 - 4.0*np.arctan(r/rc)**2.0) + 16.0*r**3.0*np.log(r**(-2.0)*rc**2.0 + 1.0) - 4.0*rc*(-3.0*np.pi*r**2.0 + np.pi*rc**2.0 + r*rc - (-6.0*r**2.0 + 2.0*rc**2.0)*np.arctan(rc/r))) + 0.0416666666666667*r**(-3.0)*rh**(-2.0)*(-3.0*r**3.0*(np.pi**2.0 - 4.0*np.arctan(r/rh)**2.0) + 16.0*r**3.0*np.log(r**(-2.0)*rh**2.0 + 1.0) - 4.0*rh*(-3.0*np.pi*r**2.0 + np.pi*rh**2.0 + r*rh - (-6.0*r**2.0 + 2.0*rh**2.0)*np.arctan(rh/r))) + 0.166666666666667*r**(-3.0)*(np.pi*rc + r - 2.0*rc*np.arctan(rc/r)) + 0.5*r**(-2.0) + 0.5*rc**(-2.0)*np.log(r**2.0/(r**2.0 + rc**2.0)) - 0.166666666666667*rc**(-2.0)*np.log(r**(-2.0)*rc**2.0 + 1.0) - 0.666666666666667*rh**(-2.0)*np.log(r**(-2.0)*rc**2.0 + 1.0) + rc*rh**(-2.0)*np.arctan(rc/r)/r)/np.pi - 1.0*G*Mtot*r**1.0*(r**2.0 + rh**2.0)*(rc + rh)*(rc**2.0 - rh**2.0)**(-2.0)*(-0.5*np.pi*rc*rh**(-2.0)/r + 0.0416666666666667*r**(-3.0)*rc**(-2.0)*(-3.0*r**3.0*(np.pi**2.0 - 4.0*np.arctan(r/rc)**2.0) + 16.0*r**3.0*np.log(r**(-2.0)*rc**2.0 + 1.0) - 4.0*rc*(-3.0*np.pi*r**2.0 + np.pi*rc**2.0 + r*rc - (-6.0*r**2.0 + 2.0*rc**2.0)*np.arctan(rc/r))) + 0.0416666666666667*r**(-3.0)*rh**(-2.0)*(-3.0*r**3.0*(np.pi**2.0 - 4.0*np.arctan(r/rh)**2.0) + 16.0*r**3.0*np.log(r**(-2.0)*rh**2.0 + 1.0) - 4.0*rh*(-3.0*np.pi*r**2.0 + np.pi*rh**2.0 + r*rh - (-6.0*r**2.0 + 2.0*rh**2.0)*np.arctan(rh/r))) + 0.166666666666667*r**(-3.0)*(np.pi*rc + r - 2.0*rc*np.arctan(rc/r)) + 0.5*r**(-2.0) + 0.5*rc**(-2.0)*np.log(r**2.0/(r**2.0 + rc**2.0)) - 0.166666666666667*rc**(-2.0)*np.log(r**(-2.0)*rc**2.0 + 1.0) - 0.666666666666667*rh**(-2.0)*np.log(r**(-2.0)*rc**2.0 + 1.0) + rc*rh**(-2.0)*np.arctan(rc/r)/r)/np.pi - G*Mtot*(r**2.0 + rc**2.0)*(r**2.0 + rh**2.0)*(rc + rh)*(rc**2.0 - rh**2.0)**(-2.0)*(0.5*np.pi*rc*rh**(-2.0)/r**2 - 0.125*r**(-4.0)*rc**(-2.0)*(-3.0*r**3.0*(np.pi**2.0 - 4.0*np.arctan(r/rc)**2.0) + 16.0*r**3.0*np.log(r**(-2.0)*rc**2.0 + 1.0) - 4.0*rc*(-3.0*np.pi*r**2.0 + np.pi*rc**2.0 + r*rc - (-6.0*r**2.0 + 2.0*rc**2.0)*np.arctan(rc/r))) - 0.125*r**(-4.0)*rh**(-2.0)*(-3.0*r**3.0*(np.pi**2.0 - 4.0*np.arctan(r/rh)**2.0) + 16.0*r**3.0*np.log(r**(-2.0)*rh**2.0 + 1.0) - 4.0*rh*(-3.0*np.pi*r**2.0 + np.pi*rh**2.0 + r*rh - (-6.0*r**2.0 + 2.0*rh**2.0)*np.arctan(rh/r))) - 0.5*r**(-4.0)*(np.pi*rc + r - 2.0*rc*np.arctan(rc/r)) + 0.0416666666666667*r**(-3.0)*rc**(-2.0)*(-9.0*r**2.0*(np.pi**2.0 - 4.0*np.arctan(r/rc)**2.0) + 48.0*r**2.0*np.log(r**(-2.0)*rc**2.0 + 1.0) + 24.0*r**3.0*np.arctan(r/rc)**1.0/(rc*(r**2/rc**2 + 1)) - 4.0*rc*(-6.0*np.pi*r**1.0 + 12.0*r**1.0*np.arctan(rc/r) + rc + rc*(-6.0*r**2.0 + 2.0*rc**2.0)/(r**2*(1 + rc**2/r**2))) - 32.0*rc**2.0/(r**(-2.0)*rc**2.0 + 1.0)) + 1.33333333333333*r**(-3.0)*rc**2.0*rh**(-2.0)/(r**(-2.0)*rc**2.0 + 1.0) + 0.0416666666666667*r**(-3.0)*rh**(-2.0)*(-9.0*r**2.0*(np.pi**2.0 - 4.0*np.arctan(r/rh)**2.0) + 48.0*r**2.0*np.log(r**(-2.0)*rh**2.0 + 1.0) + 24.0*r**3.0*np.arctan(r/rh)**1.0/(rh*(r**2/rh**2 + 1)) - 4.0*rh*(-6.0*np.pi*r**1.0 + 12.0*r**1.0*np.arctan(rh/r) + rh + rh*(-6.0*r**2.0 + 2.0*rh**2.0)/(r**2*(1 + rh**2/r**2))) - 32.0*rh**2.0/(r**(-2.0)*rh**2.0 + 1.0)) + 0.166666666666667*r**(-3.0)*(1 + 2.0*rc**2/(r**2*(1 + rc**2/r**2))) - 1.0*r**(-3.0) + 0.333333333333333*r**(-3.0)/(r**(-2.0)*rc**2.0 + 1.0) + 0.5*r**(-2.0)*rc**(-2.0)*(r**2.0 + rc**2.0)*(2.0*r**1.0/(r**2.0 + rc**2.0) - 2.0*r**3.0/(r**2.0 + rc**2.0)**2) - rc*rh**(-2.0)*np.arctan(rc/r)/r**2 - rc**2*rh**(-2.0)/(r**3*(1 + rc**2/r**2)))/(2*np.pi))/(G*Mtot*(r**2.0 + rc**2.0)*(r**2.0 + rh**2.0)*(rc + rh)*(-0.5*np.pi*rc*rh**(-2.0)/r + 0.0416666666666667*r**(-3.0)*rc**(-2.0)*(-3.0*r**3.0*(np.pi**2.0 - 4.0*np.arctan(r/rc)**2.0) + 16.0*r**3.0*np.log(r**(-2.0)*rc**2.0 + 1.0) - 4.0*rc*(-3.0*np.pi*r**2.0 + np.pi*rc**2.0 + r*rc - (-6.0*r**2.0 + 2.0*rc**2.0)*np.arctan(rc/r))) + 0.0416666666666667*r**(-3.0)*rh**(-2.0)*(-3.0*r**3.0*(np.pi**2.0 - 4.0*np.arctan(r/rh)**2.0) + 16.0*r**3.0*np.log(r**(-2.0)*rh**2.0 + 1.0) - 4.0*rh*(-3.0*np.pi*r**2.0 + np.pi*rh**2.0 + r*rh - (-6.0*r**2.0 + 2.0*rh**2.0)*np.arctan(rh/r))) + 0.166666666666667*r**(-3.0)*(np.pi*rc + r - 2.0*rc*np.arctan(rc/r)) + 0.5*r**(-2.0) + 0.5*rc**(-2.0)*np.log(r**2.0/(r**2.0 + rc**2.0)) - 0.166666666666667*rc**(-2.0)*np.log(r**(-2.0)*rc**2.0 + 1.0) - 0.666666666666667*rh**(-2.0)*np.log(r**(-2.0)*rc**2.0 + 1.0) + rc*rh**(-2.0)*np.arctan(rc/r)/r))

def get_sigma_far(rh, r, rc, G, Mtot, init_guess=0., arktan=None, lahg=None, squirt=None, pie=None):
    return squirt(6.*G*Mtot*(r**2. + rh**2.)*(r**2. + rc**2.)*(pie*(rh-rc))**(-1.)*(pie**2./(8.*rh**4.) + \
             pie/(6.*rh*r**3.) + 1./(6.*rh**2.*r**2.) - pie/(2.*rh**3.*r) - arktan(r/rh)**2./(2.*rh**4.) - \
             arktan(rh/r)/(3.*rh*r**3.) + arktan(rh/r)/(rh**3.*r) - 2.*lahg(1.+rh**2./r**2.)/(3.*rh**4.) - \
             rc*pie*(rh**3. - 3.*rh*r**2. + 3.*r**3.*arktan(rh/r))/(6.*rh**5.*r**3.))) - init_guess

def get_dsigfar_drh(rh, r, rc, G, Mtot, init_guess=0., arktan=None, lahg=None, squirt=None, pie=None):
    return 2.44948974278318*(np.pi*(-rc + rh))**1.0*np.sqrt(G*Mtot*(np.pi*(-rc + rh))**(-1.0)*(r**2.0 + rc**2.0)*(r**2.0 + rh**2.0)*(-0.166666666666667*np.pi*r**(-3.0)*rc*rh**(-5.0)*(-3.0*r**2.0*rh + 3.0*r**3.0*np.arctan(rh/r) + rh**3.0) + 0.166666666666667*np.pi*r**(-3.0)/rh - 0.5*np.pi*rh**(-3.0)/r + 0.125*np.pi**2.0*rh**(-4.0) - 0.333333333333333*r**(-3.0)*np.arctan(rh/r)/rh + 0.166666666666667*r**(-2.0)*rh**(-2.0) - 0.666666666666667*rh**(-4.0)*np.log(r**(-2.0)*rh**2.0 + 1.0) - 0.5*rh**(-4.0)*np.arctan(r/rh)**2.0 + rh**(-3.0)*np.arctan(rh/r)/r))*(G*Mtot*rh**1.0*(np.pi*(-rc + rh))**(-1.0)*(r**2.0 + rc**2.0)*(-0.166666666666667*np.pi*r**(-3.0)*rc*rh**(-5.0)*(-3.0*r**2.0*rh + 3.0*r**3.0*np.arctan(rh/r) + rh**3.0) + 0.166666666666667*np.pi*r**(-3.0)/rh - 0.5*np.pi*rh**(-3.0)/r + 0.125*np.pi**2.0*rh**(-4.0) - 0.333333333333333*r**(-3.0)*np.arctan(rh/r)/rh + 0.166666666666667*r**(-2.0)*rh**(-2.0) - 0.666666666666667*rh**(-4.0)*np.log(r**(-2.0)*rh**2.0 + 1.0) - 0.5*rh**(-4.0)*np.arctan(r/rh)**2.0 + rh**(-3.0)*np.arctan(rh/r)/r) + G*Mtot*(np.pi*(-rc + rh))**(-1.0)*(r**2.0 + rc**2.0)*(r**2.0 + rh**2.0)*(0.833333333333333*np.pi*r**(-3.0)*rc*rh**(-6.0)*(-3.0*r**2.0*rh + 3.0*r**3.0*np.arctan(rh/r) + rh**3.0) - 0.166666666666667*np.pi*r**(-3.0)*rc*rh**(-5.0)*(-3.0*r**2.0 + 3.0*r**2.0/(1 + rh**2/r**2) + 3.0*rh**2.0) - 0.166666666666667*np.pi*r**(-3.0)/rh**2 + 1.5*np.pi*rh**(-4.0)/r - 0.5*np.pi**2.0*rh**(-5.0) - 0.333333333333333*r**(-4.0)/(rh*(1 + rh**2/r**2)) + 0.333333333333333*r**(-3.0)*np.arctan(rh/r)/rh**2 - 0.333333333333333*r**(-2.0)*rh**(-3.0) - 1.33333333333333*r**(-2.0)*rh**(-3.0)/(r**(-2.0)*rh**2.0 + 1.0) + 1.0*r*rh**(-6.0)*np.arctan(r/rh)**1.0/(r**2/rh**2 + 1) + 2.66666666666667*rh**(-5.0)*np.log(r**(-2.0)*rh**2.0 + 1.0) + 2.0*rh**(-5.0)*np.arctan(r/rh)**2.0 - 3.0*rh**(-4.0)*np.arctan(rh/r)/r + rh**(-3.0)/(r**2*(1 + rh**2/r**2)))/2 - 0.5*G*Mtot*(np.pi*(-rc + rh))**(-1.0)*(r**2.0 + rc**2.0)*(r**2.0 + rh**2.0)*(-0.166666666666667*np.pi*r**(-3.0)*rc*rh**(-5.0)*(-3.0*r**2.0*rh + 3.0*r**3.0*np.arctan(rh/r) + rh**3.0) + 0.166666666666667*np.pi*r**(-3.0)/rh - 0.5*np.pi*rh**(-3.0)/r + 0.125*np.pi**2.0*rh**(-4.0) - 0.333333333333333*r**(-3.0)*np.arctan(rh/r)/rh + 0.166666666666667*r**(-2.0)*rh**(-2.0) - 0.666666666666667*rh**(-4.0)*np.log(r**(-2.0)*rh**2.0 + 1.0) - 0.5*rh**(-4.0)*np.arctan(r/rh)**2.0 + rh**(-3.0)*np.arctan(rh/r)/r)/(-rc + rh)) / (G*Mtot*(r**2.0 + rc**2.0)*(r**2.0 + rh**2.0)*(-0.166666666666667*np.pi*r**(-3.0)*rc*rh**(-5.0)*(-3.0*r**2.0*rh + 3.0*r**3.0*np.arctan(rh/r) + rh**3.0) + 0.166666666666667*np.pi*r**(-3.0)/rh - 0.5*np.pi*rh**(-3.0)/r + 0.125*np.pi**2.0*rh**(-4.0) - 0.333333333333333*r**(-3.0)*np.arctan(rh/r)/rh + 0.166666666666667*r**(-2.0)*rh**(-2.0) - 0.666666666666667*rh**(-4.0)*np.log(r**(-2.0)*rh**2.0 + 1.0) - 0.5*rh**(-4.0)*np.arctan(r/rh)**2.0 + rh**(-3.0)*np.arctan(rh/r)/r))

def get_dsigfar_dr(rh, r, rc, G, Mtot, init_guess=0., arktan=None, lahg=None, squirt=None, pie=None):
    return 2.44948974278318*(np.pi*(-rc + rh))**1.0*np.sqrt(G*Mtot*(np.pi*(-rc + rh))**(-1.0)*(r**2.0 + rc**2.0)*(r**2.0 + rh**2.0)*(-0.166666666666667*np.pi*r**(-3.0)*rc*rh**(-5.0)*(-3.0*r**2.0*rh + 3.0*r**3.0*np.arctan(rh/r) + rh**3.0) + 0.166666666666667*np.pi*r**(-3.0)/rh - 0.5*np.pi*rh**(-3.0)/r + 0.125*np.pi**2.0*rh**(-4.0) - 0.333333333333333*r**(-3.0)*np.arctan(rh/r)/rh + 0.166666666666667*r**(-2.0)*rh**(-2.0) - 0.666666666666667*rh**(-4.0)*np.log(r**(-2.0)*rh**2.0 + 1.0) - 0.5*rh**(-4.0)*np.arctan(r/rh)**2.0 + rh**(-3.0)*np.arctan(rh/r)/r))*(G*Mtot*r**1.0*(np.pi*(-rc + rh))**(-1.0)*(r**2.0 + rc**2.0)*(-0.166666666666667*np.pi*r**(-3.0)*rc*rh**(-5.0)*(-3.0*r**2.0*rh + 3.0*r**3.0*np.arctan(rh/r) + rh**3.0) + 0.166666666666667*np.pi*r**(-3.0)/rh - 0.5*np.pi*rh**(-3.0)/r + 0.125*np.pi**2.0*rh**(-4.0) - 0.333333333333333*r**(-3.0)*np.arctan(rh/r)/rh + 0.166666666666667*r**(-2.0)*rh**(-2.0) - 0.666666666666667*rh**(-4.0)*np.log(r**(-2.0)*rh**2.0 + 1.0) - 0.5*rh**(-4.0)*np.arctan(r/rh)**2.0 + rh**(-3.0)*np.arctan(rh/r)/r) + G*Mtot*r**1.0*(np.pi*(-rc + rh))**(-1.0)*(r**2.0 + rh**2.0)*(-0.166666666666667*np.pi*r**(-3.0)*rc*rh**(-5.0)*(-3.0*r**2.0*rh + 3.0*r**3.0*np.arctan(rh/r) + rh**3.0) + 0.166666666666667*np.pi*r**(-3.0)/rh - 0.5*np.pi*rh**(-3.0)/r + 0.125*np.pi**2.0*rh**(-4.0) - 0.333333333333333*r**(-3.0)*np.arctan(rh/r)/rh + 0.166666666666667*r**(-2.0)*rh**(-2.0) - 0.666666666666667*rh**(-4.0)*np.log(r**(-2.0)*rh**2.0 + 1.0) - 0.5*rh**(-4.0)*np.arctan(r/rh)**2.0 + rh**(-3.0)*np.arctan(rh/r)/r) + G*Mtot*(np.pi*(-rc + rh))**(-1.0)*(r**2.0 + rc**2.0)*(r**2.0 + rh**2.0)*(0.5*np.pi*r**(-4.0)*rc*rh**(-5.0)*(-3.0*r**2.0*rh + 3.0*r**3.0*np.arctan(rh/r) + rh**3.0) - 0.5*np.pi*r**(-4.0)/rh - 0.166666666666667*np.pi*r**(-3.0)*rc*rh**(-5.0)*(-6.0*r**1.0*rh - 3.0*r**1.0*rh/(1 + rh**2/r**2) + 9.0*r**2.0*np.arctan(rh/r)) + 0.5*np.pi*rh**(-3.0)/r**2 + 0.333333333333333*r**(-5.0)/(1 + rh**2/r**2) + 1.0*r**(-4.0)*np.arctan(rh/r)/rh - 0.333333333333333*r**(-3.0)*rh**(-2.0) + 1.33333333333333*r**(-3.0)*rh**(-2.0)/(r**(-2.0)*rh**2.0 + 1.0) - 1.0*rh**(-5.0)*np.arctan(r/rh)**1.0/(r**2/rh**2 + 1) - rh**(-3.0)*np.arctan(rh/r)/r**2 - rh**(-2.0)/(r**3*(1 + rh**2/r**2)))/2)/(G*Mtot*(r**2.0 + rc**2.0)*(r**2.0 + rh**2.0)*(-0.166666666666667*np.pi*r**(-3.0)*rc*rh**(-5.0)*(-3.0*r**2.0*rh + 3.0*r**3.0*np.arctan(rh/r) + rh**3.0) + 0.166666666666667*np.pi*r**(-3.0)/rh - 0.5*np.pi*rh**(-3.0)/r + 0.125*np.pi**2.0*rh**(-4.0) - 0.333333333333333*r**(-3.0)*np.arctan(rh/r)/rh + 0.166666666666667*r**(-2.0)*rh**(-2.0) - 0.666666666666667*rh**(-4.0)*np.log(r**(-2.0)*rh**2.0 + 1.0) - 0.5*rh**(-4.0)*np.arctan(r/rh)**2.0 + rh**(-3.0)*np.arctan(rh/r)/r))


def run():
    rc = sp.Symbol('rc')
    r = sp.Symbol('r')
    rh = sp.Symbol('rh')
    Mtot = sp.Symbol('Mtot')
    G = sp.Symbol('G')
    eqn = 2.44948974278318*(np.pi*(-rc + rh))**1.0*sp.sqrt(G*Mtot*(np.pi*(-rc + rh))**(-1.0)*(r**2.0 + rc**2.0)*(r**2.0 + rh**2.0)*(-0.166666666666667*np.pi*r**(-3.0)*rc*rh**(-5.0)*(-3.0*r**2.0*rh + 3.0*r**3.0*sp.atan(rh/r) + rh**3.0) + 0.166666666666667*np.pi*r**(-3.0)/rh - 0.5*np.pi*rh**(-3.0)/r + 0.125*np.pi**2.0*rh**(-4.0) - 0.333333333333333*r**(-3.0)*sp.atan(rh/r)/rh + 0.166666666666667*r**(-2.0)*rh**(-2.0) - 0.666666666666667*rh**(-4.0)*sp.log(r**(-2.0)*rh**2.0 + 1.0) - 0.5*rh**(-4.0)*sp.atan(r/rh)**2.0 + rh**(-3.0)*sp.atan(rh/r)/r))*(G*Mtot*rh**1.0*(np.pi*(-rc + rh))**(-1.0)*(r**2.0 + rc**2.0)*(-0.166666666666667*np.pi*r**(-3.0)*rc*rh**(-5.0)*(-3.0*r**2.0*rh + 3.0*r**3.0*sp.atan(rh/r) + rh**3.0) + 0.166666666666667*np.pi*r**(-3.0)/rh - 0.5*np.pi*rh**(-3.0)/r + 0.125*np.pi**2.0*rh**(-4.0) - 0.333333333333333*r**(-3.0)*sp.atan(rh/r)/rh + 0.166666666666667*r**(-2.0)*rh**(-2.0) - 0.666666666666667*rh**(-4.0)*sp.log(r**(-2.0)*rh**2.0 + 1.0) - 0.5*rh**(-4.0)*sp.atan(r/rh)**2.0 + rh**(-3.0)*sp.atan(rh/r)/r) + G*Mtot*(np.pi*(-rc + rh))**(-1.0)*(r**2.0 + rc**2.0)*(r**2.0 + rh**2.0)*(0.833333333333333*np.pi*r**(-3.0)*rc*rh**(-6.0)*(-3.0*r**2.0*rh + 3.0*r**3.0*sp.atan(rh/r) + rh**3.0) - 0.166666666666667*np.pi*r**(-3.0)*rc*rh**(-5.0)*(-3.0*r**2.0 + 3.0*r**2.0/(1 + rh**2/r**2) + 3.0*rh**2.0) - 0.166666666666667*np.pi*r**(-3.0)/rh**2 + 1.5*np.pi*rh**(-4.0)/r - 0.5*np.pi**2.0*rh**(-5.0) - 0.333333333333333*r**(-4.0)/(rh*(1 + rh**2/r**2)) + 0.333333333333333*r**(-3.0)*sp.atan(rh/r)/rh**2 - 0.333333333333333*r**(-2.0)*rh**(-3.0) - 1.33333333333333*r**(-2.0)*rh**(-3.0)/(r**(-2.0)*rh**2.0 + 1.0) + 1.0*r*rh**(-6.0)*sp.atan(r/rh)**1.0/(r**2/rh**2 + 1) + 2.66666666666667*rh**(-5.0)*sp.log(r**(-2.0)*rh**2.0 + 1.0) + 2.0*rh**(-5.0)*sp.atan(r/rh)**2.0 - 3.0*rh**(-4.0)*sp.atan(rh/r)/r + rh**(-3.0)/(r**2*(1 + rh**2/r**2)))/2 - 0.5*G*Mtot*(np.pi*(-rc + rh))**(-1.0)*(r**2.0 + rc**2.0)*(r**2.0 + rh**2.0)*(-0.166666666666667*np.pi*r**(-3.0)*rc*rh**(-5.0)*(-3.0*r**2.0*rh + 3.0*r**3.0*sp.atan(rh/r) + rh**3.0) + 0.166666666666667*np.pi*r**(-3.0)/rh - 0.5*np.pi*rh**(-3.0)/r + 0.125*np.pi**2.0*rh**(-4.0) - 0.333333333333333*r**(-3.0)*sp.atan(rh/r)/rh + 0.166666666666667*r**(-2.0)*rh**(-2.0) - 0.666666666666667*rh**(-4.0)*sp.log(r**(-2.0)*rh**2.0 + 1.0) - 0.5*rh**(-4.0)*sp.atan(r/rh)**2.0 + rh**(-3.0)*sp.atan(rh/r)/r)/(-rc + rh)) / (G*Mtot*(r**2.0 + rc**2.0)*(r**2.0 + rh**2.0)*(-0.166666666666667*np.pi*r**(-3.0)*rc*rh**(-5.0)*(-3.0*r**2.0*rh + 3.0*r**3.0*sp.atan(rh/r) + rh**3.0) + 0.166666666666667*np.pi*r**(-3.0)/rh - 0.5*np.pi*rh**(-3.0)/r + 0.125*np.pi**2.0*rh**(-4.0) - 0.333333333333333*r**(-3.0)*sp.atan(rh/r)/rh + 0.166666666666667*r**(-2.0)*rh**(-2.0) - 0.666666666666667*rh**(-4.0)*sp.log(r**(-2.0)*rh**2.0 + 1.0) - 0.5*rh**(-4.0)*sp.atan(r/rh)**2.0 + rh**(-3.0)*sp.atan(rh/r)/r))
    print(sp.simplify(eqn))



if __name__ == '__main__':
    run()