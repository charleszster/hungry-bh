#!/Users/chaz/anaconda/bin/python  #This is the python interpreter I use on Mac
##!/home/charles/anaconda/bin/python 

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import Constants

H0 = Constants.H0
WM = Constants.WM
WV = Constants.WV
G = Constants.G
h = 0.7
delta_vir = 200.

M_200 = 1.e13 #Msun
_z = 0.

R_C = 100. #pc

def get_nfw_params():
    NFW_params = {}
    NFW_params['H_z'] = H0*(WV + (1 - WV - WM)*(1 + _z)**2 + WM*(1 + _z)**3)**0.5
    NFW_params['R_200'] = (M_200*G/(100.*NFW_params['H_z']**2))**(1./3.)
    NFW_params['c'] = 10**(1.025-0.097*np.log10(M_200/(1.e12/h))) #equation 7 from Zonoozi et al paper
    NFW_params['R_S'] = NFW_params['R_200']/NFW_params['c']
    NFW_params['rho_crit'] = 3*NFW_params['H_z']**2/(8*np.pi*G)
    NFW_params['rho_0'] = NFW_params['rho_crit']*delta_vir*NFW_params['c']**3/(3*np.log(1.+NFW_params['c'])- 
                                                                              NFW_params['c']/(1.+NFW_params['c'])) #equation 4 from Zonoozi et al paper 
    print 'NFW profile parameters:'
    print 'H = %.5e; R_200 = %.3f; c = %.3f; R_S = %.3f; rho_crit = %.3e; rho_0 = %.3e' % \
           (NFW_params['H_z'], NFW_params['R_200'], NFW_params['c'], NFW_params['R_S'], NFW_params['rho_crit'], 
            NFW_params['rho_0'])
    return NFW_params

def rho_c_SO(R_H, R_C):
    return M_200*(R_H+R_C)/(2*np.pi**2*R_C**2*R_H**2)

def M_NFW(r, rho_0, R_S):
    return 4*np.pi*rho_0*R_S**3*(np.log(1+r/R_S) - (r/R_S)/(1+r/R_S))

def M_SO(r, R_H, R_C):
    return 4*np.pi*R_C**2*R_H**2*rho_c_SO(R_H, R_C)*(R_H*np.arctan(r/R_H)-R_C*np.arctan(r/R_C))/(R_H**2-R_C**2)

def rho_SO_func(x, R_H):
    return rho_c_SO(R_H, R_C)/((1+x**2/R_C**2)*(1+x**2/R_H**2))

def get_so(R, rho_NFW):
    R_H_guess = 100000. #pc
    R_H = scipy.optimize.curve_fit(rho_SO_func, R[R>=R_C], rho_NFW[R>=R_C], R_H_guess)[0][0]
    rho_c = M_200*(R_H + R_C)/(2*np.pi**2*R_C**2*R_H**2) 
    rho_SO = rho_c/((1+R**2/R_C**2)*(1+R**2/R_H**2))
    print
    print 'SO profile parameters:'
    print 'R_C = %.3f; R_H optimized = %.3f; rho_c = %.3f' % (R_C, R_H, rho_c)
    return rho_SO

def run():
    NFW_params = get_nfw_params()
    R = np.linspace(0.01, NFW_params['R_200'], num=100000)
    rho_NFW = NFW_params['rho_0']/((R/NFW_params['R_S'])*(1.+R/NFW_params['R_S'])**2) #equation 2 from Zonoozi et al 
    '''
    plt.figure()
    plt.loglog(R, rho_NFW, label='rho_NFW')

    rho_SO = get_so(R, rho_NFW)
    plt.loglog(R, rho_SO, label='rho_SO')
    plt.legend(loc='best')
    plt.xlabel('Log10 Distance from Center (pc)')
    plt.ylabel('Log10 rho (Mo/pc**3)')
    plt.title('NFW and SO density profiles', fontsize=12)
    plt.savefig('NFW_and_SO_Density_Profiles.png') 
    plt.show()
    '''
    r_half_mass = 265.339e3 #pc
    R_H = 265.1e3 #pc
    half_mass = M_NFW(r_half_mass, NFW_params['rho_0'], NFW_params['R_S'])
    NFW_mass = M_NFW(R, NFW_params['rho_0'], NFW_params['R_S'])
    SO_mass = M_SO(R, R_H, R_C)

    half_mass_line = np.ones(len(R))*half_mass
    r_half_mass_line = np.ones(len(R))*r_half_mass
    plt.figure()
    plt.loglog(R, NFW_mass, label='NFW Mass')
    plt.loglog(R, SO_mass, label='SO Mass')
    plt.loglog(R, half_mass_line, 'r')
    plt.loglog(r_half_mass_line, NFW_mass, 'r')
    plt.ylim([0, 1.e14])
    plt.legend(loc='best')
    plt.xlabel('Log Distance from Center (pc)')
    plt.ylabel('Log Mass (Mo)')
    plt.title('NFW and SO Mass Profiles', fontsize=12)
    plt.savefig('NFW_and_SO_Mass_Profiles.png') 
    plt.show()
    
if __name__ == '__main__':
    run()