"""" Décrit l'atmosphère standard"""

# TODO : Calibrer/rectifier le modèle avec les mesures en temps réel ?

mol = 0.029  # masse molaire de l'air (en g/mol)
R = 8.314  # constante des gaz parfaits
P_0 = 101325  # pression atmosphérique standard (en Pa) au niveau de la mer
T_0 = 300  # temperature (en K) au niveau de la mer
rho_0 = (mol * P_0) / (R * T_0)  # densite de l'air (en kg/m^3) au niveau de la mer
g = 9.81  # pesanteur (en m/s^2) au niveau de la mer
gamma = 1.4


def T(y):
    if y < 11019000:
        return T_0 - y * 10**(-3) * -6.5
    else:
        return -56.5


def P(y):
    return P_0 * (1 - 0.0065 * y / 288.15)**5.255


def rho(y):
    return (mol * P(y)) / (R * T(y))
