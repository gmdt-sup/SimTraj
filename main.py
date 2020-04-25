import math
from math import exp, cos, sin, sqrt, pi, acos, asin
import numpy as np
import matplotlib.pyplot as plt

from atmosphere import *
from fusee import *
from lancement import *
from aerodynamique import *

# TODO : Prendre en compte l'altitude de décollage
rho = rho_0
g = g_0


# Fonctions de mecanique
# TODO : Prendre en compte le moment des forces aéro
# TODO : Garder phi constant tant qu'on est sur la rampe

def fx(t, x, Vx, Vy):
    if Vx != 0 or Vy != 0:
        phi = asin(Vy / sqrt(Vx ** 2 + Vy ** 2))
    else:
        phi = angle_rampe
    Ftrainee = 1/2 * rho * C_trainee * S * (Vx ** 2 + Vy ** 2)
    Fportance = 1/2 * rho * C_portance * S * (Vx ** 2 + Vy ** 2)

    f_x = - Ftrainee * cos(phi) - Fportance * sin(phi) + Fmot(t) * cos(phi)
    return f_x / m(t)


def fy(t, y, Vy, Vx):
    if Vx != 0 or Vy != 0:
        phi = asin(Vy / sqrt(Vx ** 2 + Vy ** 2))
    else:
        phi = angle_rampe

    Ftrainee = 1 / 2 * rho * C_trainee * S * (Vx ** 2 + Vy ** 2)
    Fportance = 1 / 2 * rho * C_portance * S * (Vx ** 2 + Vy ** 2)

    f_y = - Ftrainee * sin(phi) + Fportance * cos(phi) + Fmot(t) * sin(phi)
    return f_y / m(t) - g


## Integration par Runge-Kutta d'ordre 4

# Conditions Initiales

C_I = [1, 0, 0, 0, 0, 0, angle_rampe, m_seche + m_carbu]  # [ X0 Vx0 Ax0 Y0 Vy0 Ay0 Phi0 m0]
Sol = np.array([C_I])
# Initialisation de la simulation

h = 0.01  # pas de temps (en s)
t_simulation = 50  # durée de la simulation
step = int(t_simulation // h)
t = np.linspace(0, t_simulation, step)

# Integration
i = 0
while i < np.size(t) - 1:

    kx1 = fx(t[i], Sol[i, 0], Sol[i, 1], Sol[i, 4])
    kx2 = fx(t[i] + h / 2, Sol[i, 0] + h / 2 * Sol[i, 1], Sol[i, 1] + h / 2 * kx1, Sol[i, 4])
    kx3 = fx(t[i] + h / 2, Sol[i, 0] + h / 2 * Sol[i, 1] + h ** 2 / 4 * kx1, Sol[i, 1] * h / 2 * kx2, Sol[i, 4])
    kx4 = fx(t[i] + h, Sol[i, 0] + h * Sol[i, 1] + h ** 2 / 2 * kx2, Sol[i, 1] * h * kx3, Sol[i, 4])

    ky1 = fy(t[i], Sol[i, 3], Sol[i, 4], Sol[i, 1])
    ky2 = fy(t[i] + h / 2, Sol[i, 3] + h / 2 * Sol[i, 4], Sol[i, 4], Sol[i, 1])
    ky3 = fy(t[i] + h / 2, Sol[i, 3] + h / 2 * Sol[i, 4] + h ** 2 / 4 * kx1, Sol[i, 4] * h / 2 * kx2, Sol[i, 1])
    ky4 = fy(t[i] + h, Sol[i, 3] + h * Sol[i, 4] + h ** 2 / 2 * kx2, Sol[i, 4] * h * kx3, Sol[i, 1])

    X = Sol[i, 0] + h * Sol[i, 1] + h ** 2 / 6 * (kx1 + kx2 + kx3)
    Vx = Sol[i, 1] + h / 6 * (kx1 + 2 * kx2 + 2 * kx3 + kx4)
    Ax = fx(t[i], Sol[i, 0], Sol[i, 1], Sol[i, 4])

    Y = Sol[i, 3] + h * Sol[i, 4] + h ** 2 / 6 * (ky1 + ky2 + ky3)
    Vy = Sol[i, 4] + h / 6 * (ky1 + 2 * ky2 + 2 * ky3 + ky4)
    Ay = fy(t[i], Sol[i, 3], Sol[i, 4], Sol[i, 1])

    if Vx != 0 or Vy != 0:
        Phi = asin(Vy / sqrt(Vx ** 2 + Vy ** 2))
    else:
        Phi = angle_rampe

    if Y < 0:
        while i < np.size(t) - 1:
            Sol = np.append(Sol, np.reshape(np.array([0, 0, 0, 0, 0, 0, 0, m(t[i + 1])]), (1, 8)), axis=0)
            i += 1
    else:
        Sol = np.append(Sol, np.reshape(np.array([X, Vx, Ax, Y, Vy, Ay, Phi, m(t[i + 1])]), (1, 8)), axis=0)
        i += 1

# Calculs
# plt.plot(t,Sol[:,3])
plt.plot(Sol[:, 0], Sol[:, 3])
plt.show()




