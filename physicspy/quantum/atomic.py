#!/usr/bin/env python

from __future__ import division
from numpy import abs, sqrt, min, max
from scipy import factorial, zeros, prod


def threej(j1,j2,j3,m1,m2,m3):
    """ Calculate the Wigner three-j symbol of three angular momenta 
    """
    def bad_values(j1,j2,j3,m1,m2,m3):
        """ Check validity of supplied values """
        if (j1<abs(j2-j3) or j1>(j2+j3)):
            """ Braking the triangular rule """
            return 1
        if (abs(m1)>j1 or abs(m2)>j2 or abs(m3)>j3):
            """ Braking the |m| <= j rule """
            return 1
        if m1+m2+m3 !=0:
            """ Braking the sum rule """
            return 1
        return 0

    if bad_values(j1,j2,j3,m1,m2,m3):
        return 0
        
    jphase = (-1)**(j1-j2-m3)
    fac = zeros(10,long)
    fac[0] = factorial(j1+j2-j3)
    fac[1] = factorial(j1-j2+j3)
    fac[2] = factorial(-j1+j2+j3)
    fac[3] = factorial(j1+m1)
    fac[4] = factorial(j1-m1)
    fac[5] = factorial(j2+m2)
    fac[6] = factorial(j2-m2)
    fac[7] = factorial(j3+m3)
    fac[8] = factorial(j3-m3)
    fac[9] = factorial(j1+j2+j3+1)
    jprodfac = sqrt(prod(fac[0:9])/fac[9])

    kmax = min([(j1+j2-j3), (j1-m1) , (j2+m2)])
    kmin = max([0 , -(j3-j2+m1) , -(j3-j1-m2)])

    jsum=0
    for k in range(kmin,kmax+1):
        jsfac = zeros(6,long)
        jsfac[0] = factorial(k)
        jsfac[1] = factorial(j1+j2-j3-k)
        jsfac[2] = factorial(j1-m1-k)
        jsfac[3] = factorial(j2+m2-k)
        jsfac[4] = factorial(j3-j2+m1+k)
        jsfac[5] = factorial(j3-j1-m2+k)
        jsum += (-1)**k / prod(jsfac[:])
    return jphase*jprodfac*jsum


def sixj(j1,j2,j3,l1,l2,l3):
    """ Calculate the Wigner six-j symbol of six angular momenta
    """
    def bad_values(j1,j2,j3,l1,l2,l3):
        """ Check triangular rules for supplied values """
        if (j1<(abs(j2-j3)) or j1>(j2+j3)):
            return 1
        if (j1<(abs(l2-l3)) or j1>(l2+l3)): 
            return 1
        if (l1<(abs(j2-l3)) or l1>(j2+l3)):
            return 1
        if (l1<(abs(l2-j3)) or l1>(l2+j3)):
            return 1
        return 0

    def delta(a,b,c):
        """ Calculate delta """
        fac = zeros(4,long)
        fac[0] = factorial(a+b-c)
        fac[1] = factorial(a-b+c)
        fac[2] = factorial(-a+b+c)
        fac[3] = factorial(a+b+c+1)
        return sqrt(prod(fac[0:3])/fac[3]);

    if bad_values(j1,j2,j3,l1,l2,l3):
        return 0

    jphase=(-1)**(j1+j2+l1+l2);
    proddelt=delta(j1,j2,j3)*delta(l1,l2,j3)*delta(l1,j2,l3)*delta(j1,l2,l3);

    val = zeros(7,long)
    val[0] = j1+j2+l1+l2+1
    val[1] = j1+j2-j3
    val[2] = l1+l2-j3
    val[3] = j1+l2-l3
    val[4] = l1+j2-l3
    val[5] = -j1-l1+j3+l3
    val[6] = -j2-l2+j3+l3

    kmax = min(val[0:5])
    kmin = max([0, -val[5], -val[6]])

    jsum = 0
    for k in range(kmin,kmax+1):
        jsfac = zeros(8,long)
        jsfac[0] = factorial(val[0]-k);
        jsfac[1] = factorial(k);
        jsfac[2] = factorial(val[1]-k);
        jsfac[3] = factorial(val[2]-k);
        jsfac[4] = factorial(val[3]-k);
        jsfac[5] = factorial(val[4]-k);
        jsfac[6] = factorial(val[5]+k);
        jsfac[7] = factorial(val[6]+k);
        jsum += (-1)**k * jsfac[0] / prod(jsfac[1:])
    return jphase*proddelt*jsum


def landeg(gL,gS,J,S,L):
    """ Calculating the Lande factor g,
        For fine structure:      landeg(gL,gS,J,S,L)
        For hyperfine structure: landeg(gJ,gI,F,I,J)
    """
    return gL * (J * (J + 1) - S * (S + 1) + L * (L + 1)) / (2 * J * (J + 1)) + \
        gS * (J * (J + 1) + S * (S + 1) - L * (L + 1)) / (2 * J * (J + 1))
