import math

pi = math.pi
h = 1.054E-27  # reduced Planck constant
c = 3e10  # speed of light
G = 6.67E-8  # gravitational constant
kB = 1.381E-16  # Boltzmann constant
m_prot = 1.67E-24  # proton mass
m_elec = 9.11E-28  # electron mass
alpha = 137  # e2/h/c - fine structure constant
e2 = h * c / alpha
r_elec = e2 / m_elec / c / c  # classic electron radius
sigmaT = 8 * pi / 3 * r_elec * r_elec  # Thomson electron scattering cross section
kappaT = sigmaT / m_prot
sigma = pow(pi, 2) * pow(kB, 4) / 60 * pow(h, -3) * pow(c, -2)  # sigma*T^4
gamma = 5 / 3
year = 365.25 * 24 * 3600  # seconds in a year
eVolt = 1.6e-12  # ergs
dEpp = (26.23 - 0.5) * 1e6 * eVolt  # energy output in (4p -> He) reaction minus neutrinos energy


# Equation of state
def Pressure(den, T, X, Y, Z):
    mu = 1 / (2 * X + 3 / 4 * Y + 1 / 2 * Z)
    return den / m_prot / mu * kB * T


# Zeldovich, p.70, 71:
# Energy generation rate in p-p cycle [erg/g/sm]
def Epp(den, T, X):
    T0 = (1e6, 5e6, 10e6, 15e6, 20e6, 30e6)
    e0 = (4e-9, 1.8e-3, 6.8e-2, 0.377, 1.09, 4.01)
    n = (10.6, 5.95, 4.60, 3.95, 3.64, 3.03)
    found = False
    for i in range(len(T0) - 1):
        if T < pow(T0[i] * T0[i + 1], 0.5):
            found = True
            break
        if not found:
            i = len(T0) - 1
    return den * X * X * e0[i] * pow(T / T0[i], n[i])


# Energy generation rate in CNO cycle [erg/g/sm]
def Ecno(den, T, X, XCNO):
    T0 = (6e6, 10e6, 15e6, 20e6, 30e6, 50e6, 100e6)
    e0 = (9e-10, 3.4e-4, 1.94, 4.5e2, 4.1e5, 6.2e8, 1.9e12)
    n = (27.3, 22.9, 19.9, 18.0, 15.6, 13.6, 10.2)
    found = False
    for i in range(len(T0) - 1):
        if T < pow(T0[i] * T0[i + 1], 0.5):
            found = True
            break
    if not found:
        i = len(T0) - 1
    return den * X * XCNO * e0[i] * pow(T / T0[i], n[i])


# Total energy generation
def Etot(den, T, X, Y, Z):
    return Epp(den, T, X) + Ecno(den, T, X, Z)


def MU(X, Y, Z):
    return 1 / (2 * X + 3 / 4 * Y + 1 / 2 * Z)


# Table data for opacity (taken from SSM18 model)
Tst = [2120000., 2720000., 3400000., 4640000., 6030000.,
       7270000., 8000000., 8770000., 9600000., 10200000., 10800000.,
       11400000., 12400000., 13800000., 14800000., 15500000.]
Kst = [9.94082840e+01, 3.60946746e+01,
       1.03492885e+01, 1.64912281e+00, 4.35865504e-01, 1.73913043e-01,
       1.08597285e-01, 7.14285714e-02, 4.93506494e-02, 3.93873085e-02,
       3.17164179e-02, 2.53164557e-02, 1.77439797e-02, 1.20259019e-02,
       8.96191187e-03, 7.03774792e-03]


def opacity(den, T, X, Y, Z):
    if T <= Tst[0]:
        return den * Kst[0]
    for i in range(len(Tst) - 1):
        if T <= Tst[i + 1]:
            return den * (Kst[i] + (T - Tst[i]) / (Tst[i + 1] - Tst[i]) * (Kst[i + 1] - Kst[i]))
    return den * Kst[-1]


print('Physics version 1.3 5.07.2020')
