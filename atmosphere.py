"""" Décrit l'atmosphère standard"""

mol = 0.029  # masse molaire de l'air (en g/mol)
R = 8.314  # constante des gaz parfaits
P_0 = 101325  # pression atmosphérique standard (en Pa) au niveau de la mer
T_0 = 300  # temperature (en K) au niveau de la mer
rho_0 = (mol * P_0) / (R * T_0)  # densite de l'air (en kg/m^3) au niveau de la mer
g_0 = 9.81  # pesanteur (en m/s^2) au niveau de la mer
gamma = 1.4

