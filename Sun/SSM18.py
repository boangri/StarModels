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
m1 = np.array(M18)
r1 = np.array(R18)
t1 = 1e6 * np.array(T18)
d1 = np.array(D18)
l1 = np.array(L18)
x1 = np.array(X18)
o1 = np.array(K18)
c1 = np.zeros(18)
c1[16] = 1.

# Pressure from density and temperature
p1 = d1 * t1 * (2 * x1 + .75 * (1 - Z - x1) + .5 * Z) * ph.kB / ph.m_prot


def load_data():
    return pd.DataFrame({'Mass': m1, 'Radius': r1, 'Temperature': t1, 'Density': d1, 'Luminosity': l1,
                         'Hydrogen': x1, 'Opacity': o1, 'Convection': c1, 'Pressure': p1})


"""
Data interpolation
"""


def interpolate(x, x0, x1, x2, y0, y1, y2):
    """
    Interpolation with 2nd degree polynomial.
    x must be between x0 and x2
    """
    u1 = x1 - x0
    u2 = x2 - x0
    w1 = y1 - y0
    w2 = y2 - y0
    a = (w1 * u2 - w2 * u1) / (u1 * u1 * u2 - u2 * u2 * u1)
    b = (w1 * u2 * u2 - w2 * u1 * u1) / (u1 * u2 * u2 - u2 * u1 * u1)
    u = x - x0
    y = y0 + (a * u + b) * u
    return y


def load_interpolated_data(K=12):
    """
    :param K: Number of intermediate (interpolated) layers
    :return: DataFrame 17*K+1 rows, 8 columns
    """
    N = 17  # number of original layers

    r = np.zeros(N * K + 1)
    m = np.zeros(N * K + 1)
    t = np.zeros(N * K + 1)
    d = np.zeros(N * K + 1)
    l = np.zeros(N * K + 1)
    x = np.zeros(N * K + 1)
    o = np.zeros(N * K + 1)
    p = np.zeros(N * K + 1)
    r[N * K] = r1[N]
    m[N * K] = m1[N]
    t[N * K] = t1[N]
    d[N * K] = d1[N]
    l[N * K] = l1[N]
    x[N * K] = x1[N]
    o[N * K] = o1[N]

    for i in range(N - 1):
        for j in range(K):
            R = r1[i] + j * (r1[i + 1] - r1[i]) / K
            r[j + K * i] = R
            m[j + K * i] = interpolate(R, r1[i], r1[i + 1], r1[i + 2], m1[i], m1[i + 1], m1[i + 2])
            t[j + K * i] = interpolate(R, r1[i], r1[i + 1], r1[i + 2], t1[i], t1[i + 1], t1[i + 2])
            d[j + K * i] = interpolate(R, r1[i], r1[i + 1], r1[i + 2], d1[i], d1[i + 1], d1[i + 2])
            l[j + K * i] = interpolate(R, r1[i], r1[i + 1], r1[i + 2], l1[i], l1[i + 1], l1[i + 2])
            x[j + K * i] = interpolate(R, r1[i], r1[i + 1], r1[i + 2], x1[i], x1[i + 1], x1[i + 2])
            o[j + K * i] = interpolate(R, r1[i], r1[i + 1], r1[i + 2], o1[i], o1[i + 1], o1[i + 2])
    # special processing for the 1st layer
    for j in range(1, K):
        d[j] = d1[0] + (d1[1] - d1[0])*j*j/K/K
        # m[j] = m1[1]*pow(r[j]/r1[1], 3)
        l[j] = l1[1]*pow(r[j]/r1[1], 3)
        m[j] = m[j-1] + (d[j] + d[j-1])*2/3*ph.pi*(pow(r[j], 3) - pow(r[j-1], 3))
    # special processing for the last layer
    for j in range(K):
        R = r1[N - 1] + (r1[N] - r1[N - 1]) * j / K
        r[j + K * (N - 1)] = R
        m[j + K * (N - 1)] = interpolate(R, r1[N - 2], r1[N - 1], r1[N], m1[N - 2], m1[N - 1], m1[N])
        t[j + K * (N - 1)] = interpolate(R, r1[N - 2], r1[N - 1], r1[N], t1[N - 2], t1[N - 1], t1[N])
        d[j + K * (N - 1)] = interpolate(R, r1[N - 2], r1[N - 1], r1[N], d1[N - 2], d1[N - 1], d1[N])
        l[j + K * (N - 1)] = interpolate(R, r1[N - 2], r1[N - 1], r1[N], l1[N - 2], l1[N - 1], l1[N])
        x[j + K * (N - 1)] = interpolate(R, r1[N - 2], r1[N - 1], r1[N], x1[N - 2], x1[N - 1], x1[N])
        o[j + K * (N - 1)] = interpolate(R, r1[N - 2], r1[N - 1], r1[N], o1[N - 2], o1[N - 1], o1[N])

    for j in range(N * K + 1):
        if m[j] > 1.:
            m[j] = 1.
        p[j] = ph.Pressure(d[j], t[j], x[j], 1 - Z - x[j], Z)

    return pd.DataFrame({'Mass': m, 'Radius': r, 'Temperature': t, 'Density': d, 'Luminosity': l,
                         'Hydrogen': x, 'Opacity': o, 'Pressure': p})


print('SSM18 version 1.7 4.07.2020')
