import math

pi = math.pi
h = 1.054E-27  # reduced Planck constant
c = 3e10  # speed of light
G = 6.67E-8  # gravitational constant
kB = 1.381E-16  # Boltzmann constant
m_prot = 1.67E-24  # масса протона
m_elec = 9.11E-28  # масса электрона
alpha = 137  # e2/h/c - постоянная тонкой структуры
e2 = h * c / alpha
r_elec = e2 / m_elec / c / c  # классический радиус электрона
sigmaT = 8 * pi / 3 * r_elec * r_elec  # сечение томсоновского рассеяния на электронах
kappaT = sigmaT / m_prot
sigma = pow(pi, 2) * pow(kB, 4) / 60 * pow(h, -3) * pow(c, -2)  # sigma*T^4
gamma = 5 / 3


# Уравнение состояния
def Pressure(den, T, X, Y, Z):
    mu = 1 / (2 * X + 3 / 4 * Y + 1 / 2 * Z)
    return den / m_prot / mu * kB * T


# Зельдович, стр 70, 71:
# Скорость энерговыделения в p-p реакции [эрг/с/г] деленная на X
def E0pp(den, T, X):
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
    return den * X * e0[i] * pow(T / T0[i], n[i])


# Скорость энерговыделения в CNO цикле [эрг/с/г] деленная на X
def E0cno(den, T, XCNO):
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
    return den * XCNO * e0[i] * pow(T / T0[i], n[i])


# Полное энерговыделение - домножаем на X.
def Etot(den, T, X, Y, Z):
    return X * (E0pp(den, T, X) + E0cno(den, T, Z))


print('Physics version 1.0')
