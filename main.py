import math
from math import exp, cos, sin, sqrt, pi, acos, asin
import numpy as np
import matplotlib.pyplot as plt

from atmosphere import *
from fusee import *
from lancement import *
from etat_vol import *
from aerodynamique import *

# TODO : Prendre en compte l'altitude de décollage
rho = rho_0
g = g_0

etat = EtatVol.RAMPE


# Fonctions de mecanique
# TODO : Prendre en compte le moment des forces aéro
# TODO : Garder phi constant tant qu'on est sur la rampe
# TODO : La pousée est selon l'axe de la fusée, pas selon la vitesse

def fx(t, x, y, Vx, Vy):
    if Vx != 0 or Vy != 0:
        phi = asin(Vy / sqrt(Vx ** 2 + Vy ** 2))
    else:
        phi = angle_rampe

    """if x != 0 or y != 0:
            theta = asin(y / sqrt(x ** 2 + y ** 2))
        else:
            theta = phi"""

    Ftrainee = 1/2 * rho * C_trainee * S * (Vx ** 2 + Vy ** 2)
    Fportance = 1/2 * rho * C_portance * S * (Vx ** 2 + Vy ** 2)

    f_x = - Ftrainee * cos(phi) - Fportance * sin(phi) + Fmot(t) * cos(phi)

    if etat == EtatVol.SAUVEGARDE:
        Ftrainee_para = 1 / 2 * rho * C_trainee_para * S * (Vx ** 2 + Vy ** 2)
        Fportance_para = 1 / 2 * rho * C_portance_para * S * (Vx ** 2 + Vy ** 2)

        f_x += - Ftrainee_para * cos(phi) - Fportance_para * sin(phi)

    return f_x / m(t)


def fy(t, x, y, Vx, Vy):
    if Vx != 0 or Vy != 0:
        phi = asin(Vy / sqrt(Vx ** 2 + Vy ** 2))
    else:
        phi = angle_rampe

    """if x != 0 or y != 0:
        theta = asin(y / sqrt(x ** 2 + y ** 2))
    else:
        theta = phi"""

    Ftrainee = 1 / 2 * rho * C_trainee * S * (Vx ** 2 + Vy ** 2)
    Fportance = 1 / 2 * rho * C_portance * S * (Vx ** 2 + Vy ** 2)

    f_y = - Ftrainee * sin(phi) + Fportance * cos(phi) + Fmot(t) * sin(phi)

    if etat == EtatVol.SAUVEGARDE:
        Ftrainee_para = 1 / 2 * rho * C_trainee_para * S * (Vx ** 2 + Vy ** 2)
        Fportance_para = 1 / 2 * rho * C_portance_para * S * (Vx ** 2 + Vy ** 2)

        f_y += - Ftrainee_para * sin(phi) + Fportance_para * cos(phi)

    return f_y / m(t) - g


"""Integration par Runge-Kutta d'ordre 4"""

# Conditions Initiales

C_I = [1, 0, 0, 0, 0, 0, angle_rampe, m_seche + m_carbu]  # [ X0 Vx0 Ax0 Y0 Vy0 Ay0 Phi0 m0]
Sol = np.array([C_I])
# Initialisation de la simulation

h = 0.01  # pas de temps (en s)
t_simulation = 60  # durée de la simulation
step = int(t_simulation // h)
t = np.linspace(0, t_simulation, step)

# Integration
i = 0
while i < np.size(t) - 1:
    x, Vx, y, Vy = Sol[i, 0], Sol[i, 1], Sol[i, 3], Sol[i, 4]

    kx1 = fx(t[i], x, y, Vx, Vy)
    ky1 = fy(t[i], x, y, Vx, Vy)

    kx2 = fx(t[i] + h / 2, x + h / 2 * Vx, y + h / 2 * Vy, Vx + h / 2 * kx1, Vy + h / 2 * ky1)
    ky2 = fy(t[i] + h / 2, x + h / 2 * Vx, y + h / 2 * Vy, Vx + h / 2 * kx1, Vy + h / 2 * ky1)

    kx3 = fx(t[i] + h / 2, x + h / 2 * Vx + h ** 2 / 4 * kx1, y + h / 2 * Vy + h ** 2 / 4 * ky1, Vx + h / 2 * kx2, Vy + h / 2 * ky2)
    ky3 = fy(t[i] + h / 2, x + h / 2 * Vx + h ** 2 / 4 * kx1, y + h / 2 * Vy + h ** 2 / 4 * ky1, Vx + h / 2 * kx2, Vy + h / 2 * ky2)

    kx4 = fx(t[i] + h, x + h * Vx + h ** 2 / 2 * kx2, y + h * Vy + h ** 2 / 2 * ky2, Vx + h * kx3, Vy + h * ky3)
    ky4 = fy(t[i] + h, x + h * Vx + h ** 2 / 2 * kx2, y + h * Vy + h ** 2 / 2 * ky2, Vx + h * kx3, Vy + h * ky3)

    xn = x + h * Vx + h**2 / 6 * (kx1 + kx2 + kx3)
    Vxn = Vx + h / 6 * (kx1 + 2 * kx2 + 2 * kx3 + kx4)
    axn = fx(t[i], x, y, Vx, Vy)

    yn = y + h * Vy + h ** 2 / 6 * (ky1 + ky2 + ky3)
    Vyn = Vy + h / 6 * (ky1 + 2 * ky2 + 2 * ky3 + ky4)
    ayn = fx(t[i], x, y, Vx, Vy)

    if Vx != 0 or Vy != 0:
        Phi = asin(Vy / sqrt(Vx ** 2 + Vy ** 2))
    else:
        Phi = angle_rampe

    if yn < 0:
        while i < np.size(t) - 1:
            Sol = np.append(Sol, np.reshape(np.array([0, 0, 0, 0, 0, 0, 0, m(t[i + 1])]), (1, 8)), axis=0)
            i += 1
    else:
        Sol = np.append(Sol, np.reshape(np.array([xn, Vxn, axn, yn, Vyn, ayn, Phi, m(t[i + 1])]), (1, 8)), axis=0)
        i += 1

    if etat == EtatVol.RAMPE and sqrt(xn**2 + yn**2) > longueur_rampe:
        etat = EtatVol.VOL
    elif etat == EtatVol.VOL and yn < y:
        etat = EtatVol.SAUVEGARDE
    elif etat == EtatVol.SAUVEGARDE and y <= 0:
        etat = EtatVol.ATTERI

# Calculs
# plt.plot(t,Sol[:,3])
plt.plot(Sol[:, 0], Sol[:, 3])  # Affichage de la trajectoire
plt.show()
