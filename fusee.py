from math import pi

"""Données mécaniques de la fusée et du moteur"""

"Données méca"
m_seche = 20
D_f = 0.134  # diametre fusee (en m)
C_trainee = 0.38  # TODO : Prendre en compte les variations avec le mach, la pression, ...
C_portance = 0.05  # TODO : J'ai mis ca totalement au pif + prendre en compte les variations ...
S = pi * (D_f / 2) ** 2  # rocket cross-sectional area (en m^2)
L = 2.5  # longueur fusee (en m)

"Données moteur"
m_carbu = 7
t_moteur = 5  # temps d'allumage moteur (en s)
F_moteur = 3500  # poussee moteur (en N)


def Fmot(t):
    if (t >= t_moteur):
        return 0
    else:
        return F_moteur


def m(t):
    if t >= t_moteur:
        return m_seche
    else:
        return m_seche + m_carbu * (1 - t / t_moteur)
