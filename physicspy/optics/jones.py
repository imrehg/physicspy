#!/usr/bin/env python
from __future__ import division
from numpy import sqrt, cos, sin, arctan, exp, abs, pi
from scipy import array, dot, sum

def lplate(phi,g):
    """ Waveplate matrix """
    W = dot(rotatem(-phi),dot(waveplate(g),rotatem(phi)))
    return W

def waveplate(G):
    """ Waveplate with arbitrary retardance """
    W0 = array([[exp(1j*G/2), 0], \
                [0, exp(-1j*G/2)]])
    return W0

def rotatem(phi):
    """ Rotation matrix """
    R = array([[cos(phi), sin(phi)], \
               [-sin(phi), cos(phi)]])
    return R

def pol(phi):
    """ Polarizer matrix """
    P = array([[cos(phi)**2, cos(phi)*sin(phi)], \
               [sin(phi)*cos(phi), sin(phi)**2]])
    return P

def metalmirror(n,k,th):
    """ Metal mirror matrix """
    dr = mphase(n,k,th);
    W0 = array([[dr[3]*exp(-1j*dr[1]), 0],\
                [0, dr[2]*exp(-1j*dr[0])]])
    return W0

def mphase(n,k,th):
    """ Calculate phase shift and reflectance of a metal in the s and p directions"""
    u = sqrt(0.5 *((n**2 - k**2 - sin(th)**2) + sqrt( (n**2 - k**2 - sin(th)**2)**2 + 4*n**2*k**2 )))
    v = sqrt(0.5*(-(n**2 - k**2 - sin(th)**2) + sqrt( (n**2 - k**2 - sin(th)**2)**2 + 4*n**2*k**2 )))
    ds =  arctan(2*v*cos(th)/(u**2+v**2-cos(th)**2));
    dp =  arctan(2*v*cos(th)*(n**2-k**2-2*u**2)/(u**2+v**2-(n**2+k**2)**2*cos(th)**2));
    if(dp < 0):
        dp = dp+pi;
    rs = abs((cos(th) - (u+v*1j))/(cos(th) + (u+v*1j)))
    rp = abs(((n**2 + k**2)*cos(th) - (u+v*1j))/((n**2 + k**2)*cos(th) + (u+v*1j)));
    return array([ds, dp, rs, rp])

def intensity(v):
    """ Intensity from electric field vector """
    intens = sum(abs(v)**2)
    return intens

def scanrefangle(n,k,g):
    """ Plot reflectance as a function of incidence angle """
    thl = linspace(0,pi/2,1000);
    refs = []
    reff = []
    tin = array([[1],[1]])/sqrt(2)
    for y in range(0,thl.shape[0]):
        refs.append(intensity(dot(pol(pi/4),dot(metalmirror(n,k,thl[y]),dot(lplate(0,g),tin)))))
        reff.append(intensity(dot(pol(pi/4),dot(metalmirror(n,k,thl[y]),dot(lplate(pi/2,g),tin)))))
    #~ have to set return values!

def scanplateangle(n,k,g,th,alist):
    #~ alist = linspace(-pi,pi,200);
    ref = []
    tin = array([[1],[1]])/sqrt(2)
    for a in range(0,alist.shape[0]):
        ref.append(intensity(dot(pol(pi/4),dot(metalmirror(n,k,th),dot(lplate(alist[a],g),tin)))))
    return array(ref)


