from math import cos, sin, asin

import matplotlib.pyplot as plt
import numpy as np

from atmosphere import *
from fusee import *
from lancement import *
from etat_vol import *
from aerodynamique import *

etat = EtatVol.RAMPE


## Fonctions de mecanique
# TODO : Terme dû à l'inertie et la masse pas constante


def phi(Vx, Vy):  # Angle du vecteur vitesse dans le rerentiel terrestre
    if etat == EtatVol.RAMPE:
        return angle_rampe
    else:
        return asin(Vy / sqrt(Vx ** 2 + Vy ** 2))


def fx(t, x, y, Vx, Vy, theta, thetap):  # FIXME besoin d'un nom de fonction explicite

    f_x = - trainee(y, Vx, Vy) * cos(phi(Vx, Vy)) - portance(y, Vx, Vy) * sin(phi(Vx, Vy)) + Fmot(t) * cos(theta)

    if etat == EtatVol.SAUVEGARDE:
        Ftrainee_para = 1 / 2 * rho(y) * C_trainee_para * S * (Vx ** 2 + Vy ** 2)
        Fportance_para = 1 / 2 * rho(y) * C_portance_para * S * (Vx ** 2 + Vy ** 2)

        f_x += - Ftrainee_para * cos(phi(Vx, Vy)) - Fportance_para * sin(phi(Vx, Vy))

    return f_x / m(t)


def fy(t, x, y, Vx, Vy, theta, thetap):  # FIXME besoin d'un nom de fonction explicite

    f_y = - trainee(y, Vx, Vy) * sin(phi(Vx, Vy)) + portance(y, Vx, Vy) * cos(phi(Vx, Vy)) + Fmot(t) * sin(theta) - m(
        t) * g

    if etat == EtatVol.SAUVEGARDE:
        Ftrainee_para = 1 / 2 * rho(y) * C_trainee_para * S * (Vx ** 2 + Vy ** 2)
        Fportance_para = 1 / 2 * rho(y) * C_portance_para * S * (Vx ** 2 + Vy ** 2)

        f_y += - Ftrainee_para * sin(phi(Vx, Vy)) + Fportance_para * cos(phi(Vx, Vy))

    return f_y / m(t)


def mz(t, x, y, Vx, Vy, theta, thetap):  # FIXME besoin d'un nom de fonction explicite
    if etat == EtatVol.RAMPE:  # FIXME : Et pour EtatVol.SOL ? c'est différent ?
        return 0

    alpha = theta - phi(Vx, Vy)

    m_z = - portance(y, Vx, Vy) * cos(alpha) * marge_stat(t) - trainee(y, Vx, Vy) * sin(alpha) * marge_stat(t)
    return m_z / Iz(t)


"""Integration par Runge-Kutta d'ordre 4"""

## Conditions Initiales

C_I = [0, 0, 0, y_0, 0, 0, angle_rampe, 0, 0, m_seche + m_carbu]  # [X0 Vx0 Ax0 Y0 Vy0 Ay0 theta0 thetap0 thetapp0]
Sol = np.array([C_I])
# Initialisation de la simulation

h = 0.005  # pas de temps (en s)
t_simulation = 60  # durée de la simulation
step = int(t_simulation // h)
t = np.linspace(0, t_simulation, step)

## Integration
i = 0
while i < np.size(t) - 1:
    x, Vx, y, Vy, theta, thetap = Sol[i, 0], Sol[i, 1], Sol[i, 3], Sol[i, 4], Sol[i, 6], Sol[i, 7]

    # Calcul des coefficients pour Runge Kutta (c'est moche)
    kx1 = fx(t[i], x, y, Vx, Vy, theta, thetap)
    ky1 = fy(t[i], x, y, Vx, Vy, theta, thetap)
    km1 = mz(t[i], x, y, Vx, Vy, theta, thetap)

    kx2 = fx(t[i] + h / 2, x + h / 2 * Vx, y + h / 2 * Vy, Vx + h / 2 * kx1, Vy + h / 2 * ky1, theta + h / 2 * thetap,
             thetap + h / 2 * km1)
    ky2 = fy(t[i] + h / 2, x + h / 2 * Vx, y + h / 2 * Vy, Vx + h / 2 * kx1, Vy + h / 2 * ky1, theta + h / 2 * thetap,
             thetap + h / 2 * km1)
    km2 = mz(t[i] + h / 2, x + h / 2 * Vx, y + h / 2 * Vy, Vx + h / 2 * kx1, Vy + h / 2 * ky1, theta + h / 2 * thetap,
             thetap + h / 2 * km1)

    kx3 = fx(t[i] + h / 2, x + h / 2 * Vx + h ** 2 / 4 * kx1, y + h / 2 * Vy + h ** 2 / 4 * ky1, Vx + h / 2 * kx2,
             Vy + h / 2 * ky2, theta + h / 2 * thetap + h ** 2 / 4 * km1, thetap + h / 2 * km2)
    ky3 = fy(t[i] + h / 2, x + h / 2 * Vx + h ** 2 / 4 * kx1, y + h / 2 * Vy + h ** 2 / 4 * ky1, Vx + h / 2 * kx2,
             Vy + h / 2 * ky2, theta + h / 2 * thetap + h ** 2 / 4 * km1, thetap + h / 2 * km2)
    km3 = mz(t[i] + h / 2, x + h / 2 * Vx + h ** 2 / 4 * kx1, y + h / 2 * Vy + h ** 2 / 4 * ky1, Vx + h / 2 * kx2,
             Vy + h / 2 * ky2, theta + h / 2 * thetap + h ** 2 / 4 * km1, thetap + h / 2 * km2)

    kx4 = fx(t[i] + h, x + h * Vx + h ** 2 / 2 * kx2, y + h * Vy + h ** 2 / 2 * ky2, Vx + h * kx3, Vy + h * ky3,
             theta + h * thetap + h ** 2 / 2 * km2, thetap + h * km3)
    ky4 = fy(t[i] + h, x + h * Vx + h ** 2 / 2 * kx2, y + h * Vy + h ** 2 / 2 * ky2, Vx + h * kx3, Vy + h * ky3,
             theta + h * thetap + h ** 2 / 2 * km2, thetap + h * km3)
    km4 = mz(t[i] + h, x + h * Vx + h ** 2 / 2 * kx2, y + h * Vy + h ** 2 / 2 * ky2, Vx + h * kx3, Vy + h * ky3,
             theta + h * thetap + h ** 2 / 2 * km2, thetap + h * km3)

    # Nouvelles données
    xn = x + h * Vx + h ** 2 / 6 * (kx1 + kx2 + kx3) # FIXME noms a definir explicitement ou doc a cote
    Vxn = Vx + h / 6 * (kx1 + 2 * kx2 + 2 * kx3 + kx4)
    axn = fx(t[i], x, y, Vx, Vy, theta, thetap)

    yn = y + h * Vy + h ** 2 / 6 * (ky1 + ky2 + ky3)
    Vyn = Vy + h / 6 * (ky1 + 2 * ky2 + 2 * ky3 + ky4)
    ayn = fy(t[i], x, y, Vx, Vy, theta, thetap)

    thetan = theta + h * thetap + h ** 2 / 6 * (km1 + km2 + km3)
    thetapn = thetap + h / 6 * (km1 + 2 * km2 + 2 * km3 + km4)
    thetappn = mz(t[i], x, y, Vx, Vy, theta, thetap)

    if yn < 0:
        while i < np.size(t) - 1:
            Sol = np.append(Sol, np.reshape(np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, m(t[i + 1])]), (1, 10)), axis=0)
            i += 1
    else:
        Sol = np.append(Sol, np.reshape(np.array([xn, Vxn, axn, yn, Vyn, ayn, thetan, thetapn, thetappn, m(t[i + 1])]),
                                        (1, 10)), axis=0)
        i += 1

    # Mise à jour de l'état du vol
    if etat == EtatVol.RAMPE and sqrt(xn ** 2 + yn ** 2) > longueur_rampe:
        etat = EtatVol.VOL
    elif etat == EtatVol.VOL and yn < y:
        etat = EtatVol.SAUVEGARDE
    elif etat == EtatVol.SAUVEGARDE and y <= 0:
        etat = EtatVol.ATTERI

#np.savetxt('test.txt',Sol,delimiter=',')   #pour sauver dans un fichier texte pour analyse ensuite
#plt.plot(t, Sol[:, 6])
#plt.plot(Sol[:, 0], Sol[:, 3])  # Affichage de la trajectoire
plt.show()
