#!/usr/bin/env python
""" Gaussian optics package

"""
from __future__ import division
from numpy import pi, inf, allclose, dot
from scipy import array, dot

class GaussianBeam:
    """ A Gaussian beam """
    
    def __init__(self, z0, w0, lamb):
        self.z0 = z0
        self.w0 = w0
        self.lamb = lamb

    def rrange(self):
        """ Rayleigh range """
        return pi*self.w0**2/self.lamb

    def confocal(self):
        """ Confocal parameter """
        return 2*pi*self.w0**2/self.lamb

    def rcurve(self,z):
        """ Radius of curvature """
        zz = z - self.z0
        if allclose(z,0):
            return inf
        return zz*(1 + (self.rrange()/zz)**2)

    def divergence(self):
        """ Divergence """
        return self.lamb/(pi*self.w0)

    def beamparam(self,z):
        """ Complex beam parameter """
        return z + 1j*self.rrange()

    def freespace(self,d):
        """ Free space propagation """
        q = array([[self.beamparam(0)],[1]])
        p = array([[1,d],[0,1]])
        res = dot(p,q)
        return res[0]/res[1]

    def lens(self,f):
        """ Effect of lens on complex beam parameter """
        q = array([[self.beamparam(0)],[1]])
        p = array([[1,0],[-1/f,1]])
        res = dot(p,q)
        return res[0]/res[1]
