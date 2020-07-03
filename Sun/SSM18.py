import numpy as np
import pandas as pd
import physics as ph

M = 1.99e33  # solar mass
R = 6.96e10  # solar radius
L = 3.85e33  # solar luminosity
T = pow(L / 4 / ph.pi / ph.sigma / R / R, 1 / 4)  # 5780 solar surface temperature
Z = 0.02

M18 = [0, .0099, .0385, .1038, .1620, .2100, .2580, .3100, .3900, .4700, .5500, .6900, .8300,
       .9264, .9602, .9784, .9954, 1.]

R18 = [0., .046, .076, .113, .138, .156, .173, .190, .217, .245, .275, .336, .430,
       .554, .641, .718, .849, 1.]

T18 = [15.5, 14.8, 13.8, 12.4, 11.4, 10.8, 10.2, 9.6, 8.77, 8., 7.27, 6.03,
       4.64, 3.4, 2.72, 2.12, .95, .0058]

D18 = [156.3, 133.9, 108.1, 78.9, 63.2, 53.6, 45.7, 38.5, 29.4, 22.1, 16.1,
       8.03, 2.85, 0.773, .338, .169, .050, 2.8e-7]

L18 = [0., .079, .264, .555, .718, .809, .874, .921, .964, .986, .996, 1., 1., 1., 1., 1., 1., 1.]

X18 = [.355, .417, .497, .592, .641, .668, .688, .702, .716, .724, .728, .731, .732, .732, .732, .732, .732, .732]

K18 = [1.1, 1.2, 1.3, 1.4, 1.6, 1.7, 1.8, 1.9, 2.1, 2.4, 2.8, 3.5, 4.7, 8.0, 12.2, 16.8, 10, 0.3]

# NumPy arrays
m = np.array(M18)
r = np.array(R18)
t = 1e6 * np.array(T18)
d = np.array(D18)
l = np.array(L18)
x = np.array(X18)
k = np.array(K18)
c = np.zeros(18)
c[16] = 1.

# Pressure from density and temperature
p = d * t * (2 * x + .75 * (1 - Z - x) + .5 * Z) * ph.kB / ph.m_prot


def load_data():
    return pd.DataFrame({'Mass': m, 'Radius': r, 'Temperature': t, 'Density': d, 'Luminosity': l,
                         'Hydrogen': x, 'Opacity': k, 'Convection': c, 'Pressure': p})


print('SSM18 version 1.1')
