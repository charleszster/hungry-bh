#!/home/charles/anaconda/bin/python
# #!/Users/chaz/anaconda/bin/python  #This is the python interpreter I use on Mac

import CC
import numpy as np
import matplotlib.pyplot as plt
import os
import scipy.optimize
import Constants
import analyze_clusters

H0 = Constants.H0
WM = Constants.WM
WV = Constants.WV

def func(x, a, b, c):
#    return a*np.exp(b*x)
    return a*x**(-b)+c

def run():
#    z_s = [4, 3.1, 2.8, 2.5, 2.2, 2.0, 1.9, 1.75, 1.6, 1.5, 1.35, 1.3, 1.25, 1.2, 1.15, 1.1, 1.05, 1.0, 0.95, 0.9, 0.85,
#         0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.45, 0.4, 0.35, 0.25, 0.2, 0.15, 0.1, 0.05, 0.0]

    z_s = np.linspace(4., 0., num=1000)
    age = []
    for z in z_s:
        input_params= [z, H0, WM, WV]
        age.append(CC.get_output_params(input_params))

    sigma = np.ones(len(age))
    sigma[[0, -1]] = 0.02
    popt, pcov = scipy.optimize.curve_fit(func, age, z_s, sigma=sigma)
    print popt
    z_est = np.array([func(x, *tuple(popt)) for x in age])

    plt.figure(1)
    plt.plot(age, z_s, label='t-to-z func')
    plt.plot(age, z_est, label='curve fit')
    locs, labels = plt.xticks()
    plt.setp(labels, rotation=45)
    plt.xticks(np.arange(1.5, 14., 0.5))
    plt.yticks(np.arange(0., 4.5, 0.5))
    plt.legend()
    plt.grid(b=True, which='both', color='0.65',linestyle='-')
    plt.xlabel('Age (Gyr)')
    plt.ylabel('Z')
    
    plt.figure(2)
    diff = [(zest-zs)/zs for (zest, zs) in zip(z_est, z_s)]
    plt.plot(age, diff)
    plt.xlabel('Age (Gyr)')
    plt.ylabel('(z_est-z_s)/z_s')
    locs, labels = plt.xticks()
    plt.setp(labels, rotation=45)
    plt.xticks(np.arange(1.5, 14., 0.5))
    plt.yticks(np.arange(-0.5, 9., 0.5))
    plt.grid(b=True, which='both', color='0.65',linestyle='-')
    plt.title('Difference ratio between curve fit and z_function')
    plt.show()

if __name__ == '__main__':
    run()
