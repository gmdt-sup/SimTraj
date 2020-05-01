"""Données aérodynamiques de la fusée"""
from math import exp

from atmosphere import *
from fusee import *


def C_x(mach, C_x_init):
    if mach < 0.95:
        return C_x0 / sqrt(1 - mach ** 2)
    elif 0.95 <= mach < 1.05:
        return 3.2 * C_x0
    else:
        return 0.3 * exp(-8 * (mach - 1.05)) + 1.8 * C_x_init


def trainee(y, Vx, Vy):
    mach = mach_number(sqrt(Vx ** 2 + Vy ** 2), y)
    return 1 / 2 * rho(y) * C_x(mach, C_x0) * S * (Vx ** 2 + Vy ** 2)


def portance(y, Vx, Vy):
    return 1 / 2 * rho(y) * C_l0 * S * (Vx ** 2 + Vy ** 2)
