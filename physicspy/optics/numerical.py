#!/usr/bin/env python
""" Collection of numerically evaluated optical parameters
"""

def edlen(P,T,k,f):
    """ Index of refraction of air, using the Edlen formula

    Input parameters:
    P - air pressure in Pa
    T - the temperature in C
    k - the vacuum wave number kL/(2*Pi) in um^-1
    f - partial pressure of water vapor in the air, in Pa (can be
        calculated from the relative humidity using the Goff-Gratch equation.
    """
    return 1 + ((8342.54 + 2406147 / (130 - k ** 2) + 15998 / (38.9 - k **2)) * \
        (P / 96095.43) * ((1 + 1e-8 * (0.601 - 0.00972 * T) * P) / (1 + 0.0036610 * T)) - \
        f * (0.037345 - 0.000401 * k ** 2)) * 1e-8

