#!/usr/bin/env python

from __future__ import division
from numpy import abs, sqrt, min, max
from scipy import factorial, zeros, prod

def threej(j1,j2,j3,m1,m2,m3):
    """ Calculate the three-j symbol of three angular momenta """

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
    fac = zeros(10,int)
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
        jsfac = zeros(6,int)
        jsfac[0] = factorial(k)
        jsfac[1] = factorial(j1+j2-j3-k)
        jsfac[2] = factorial(j1-m1-k)
        jsfac[3] = factorial(j2+m2-k)
        jsfac[4] = factorial(j3-j2+m1+k)
        jsfac[5] = factorial(j3-j1-m2+k)
        jsum += (-1)**k / prod(jsfac[:])
   
    return jphase*jprodfac*jsum
