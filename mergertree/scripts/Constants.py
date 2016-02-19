#!/Users/chaz/anaconda/bin/python
# -*- coding: utf-8 -*-

import os

#ΩM = 0.28,Ωb = 0.046, ΩΛ = 0.72, σ8 = 0.82, H0 = 100 h−1 Mpc−1 = 70 km s−1 Mpc−1 , and n = 0.96
#Ω0: Total Matter Density
#ΩΛ0: Cosmological Constant

cluster_folder = os.path.join(os.path.expanduser('~'), 'Dropbox', 'Columbia', 'Ostriker', 'Cannibalism', 'mergertree',
                              'cluster_files')
plots_folder = os.path.join(os.path.expanduser('~'), 'Dropbox', 'Columbia', 'Ostriker', 'Cannibalism', 'mergertree',
                            'plots')
#H0 = 100.*0.7*1.e-6 #km*s^-1*pc^-1
H0 = 100.*0.7 #km*s^-1*Mpc^-1
WM = 0.28
WV = 0.7
G = 0.0043009211 #pc Mo**-1 (km/s)**2; https://en.wikipedia.org/wiki/Gravitational_constant
