#!/Users/chaz/anaconda/envs/py35/bin/python

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate

galaxy_num = 1
beta = 3./2
#Galaxy #1:
c = [5.55308668e-27, -3.16794807e-22, 7.51187513e-18, -9.59261428e-14, 7.13671421e-10, -3.10964658e-06, 7.48102061e-03,
4.76661860e+00]

def get_stellar_mass(timo):
    return 10**(c[0]*timo**7 + c[1]*timo**6 + c[2]*timo**5 + c[3]*timo**4 + c[4]*timo**3 + c[5]*timo**2 + c[6]*timo +\
            c[7])

def get_z_conv(timo):
    return 7.20196192*(timo/1000)**(-0.59331986) - 1.52145449

def get_eff_rad(timo):
    return 2500.*(get_stellar_mass(timo)/1.e+11)**0.73*(1.+ get_z_conv(timo))**(-0.98)

def integrand(x, r_c):
    return x**2/(1+(x/r_c)**2)**beta
    
def run():
    time = np.linspace(1565.3765, 13722.6728, 5000)
    rho_0 = []
    gas_mass_list = []
    r_c_list = []
    for t in time:
        R_e = get_eff_rad(t)
        r_c = 1.5 * R_e
        r_c_list.append(r_c)
        Mstar = get_stellar_mass(t)
        Mgas = 0.03*Mstar
        gas_mass_list.append(Mgas)
        dens_integral = scipy.integrate.quad(integrand, 0., 400000., args=(r_c))
        rho_0.append(Mgas/(4*np.pi*dens_integral[0]))
    plt.figure()
    plt.semilogy(time, rho_0)
    plt.xlabel('Time (Myr)')
    plt.ylabel('rho_0 (Msun/pc^3)')
    plt.title('Gas central density rho_0 as function of time for Galaxy #1')
    plt.show()
    
    plt.figure()
    plt.semilogy(time, gas_mass_list)
    plt.xlabel('Time (Myr)')
    plt.ylabel('Gas Mass (Msun)')
    plt.title('Gas mass as function of time for Galaxy #1')
    plt.show()

    plt.figure()
    plt.semilogy(time, r_c_list)
    plt.xlabel('Time (Myr)')
    plt.ylabel('r_c (pc)')
    plt.title('r_c as function of time for Galaxy #1')
    plt.show()

if __name__ == '__main__':
    run()
