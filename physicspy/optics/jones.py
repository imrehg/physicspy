#!/usr/bin/env python
from __future__ import division
from numpy import sqrt, cos, sin, arctan, exp, abs, pi, conj
from scipy import array, dot, sum

class JonesVector:
    """ A Jones vector class to represent polarized EM waves """
    def __init__(self,Jarray=array([1,0])):
        self.Jx = Jarray[0]
        self.Jy = Jarray[1]

    def size(self):
        """ Jones vector size """
        return sqrt(dot(self.toArray().conj(),self.toArray()).real)

    def normalize(self):
        """ Normalized Jones vector """
        result = self
        try:
            size = result.size()
            if size == 0:
                raise Exception('Zero-sized Jones vector cannot be normalized')
            result.Jx /= size
            result.Jy /= size
        except Exception as inst:
            print "Error: ",inst
        finally:
            return result

    def toArray(self):
        """ Convert into array format """
        return array([self.Jx, self.Jy])

    def rotate(self,phi):
        """ Rotated Jones vector
        
        Argument:
        phi - rotation angle in radians (clockwise is positive)
        """
        R = array([[cos(phi), sin(phi)], \
                    [-sin(phi), cos(phi)]])
        return JonesVector(dot(R, self.toArray()))

    def waveplate(self,G):
        """ Waveplate with arbitrary retardance 
        Slow axis (or "c axis") is along X 
        
        Argument:
        G - retartandance in phase units
            (e.g. one wavelength retardance is G = 2 * pi)
        """
        W0 = array([[exp(-1j*G/2), 0], \
                [0, exp(1j*G/2)]])
        return JonesVector(dot(W0, self.toArray()))

    def waveplateRot(self,phi,G):
        """ Waveplate matrix with arbitrary rotation 
        
        Arguments:
        phi - rotation angle in radians 
              (clockwise is positive)
        G - retardance in phase units
            (e.g. one wavelength retardance is G = 2 * pi)
        """
        return self.rotate(phi).waveplate(G).rotate(-phi)

    def pol(self,phi):
        """ Polarizer matrix """
        P = array([[cos(phi)**2, cos(phi)*sin(phi)], \
               [sin(phi)*cos(phi), sin(phi)**2]])
        return JonesVector(dot(P, self.toArray()))

    def mirrormetal(self,n,k,th):
        """ Reflection off a metal mirror 
        Incoming and reflected beams are assumed to be in the X plane
        """
        dr = mphase(n,k,th);
        W0 = array([[dr[3]*exp(-1j*dr[1]), 0],\
                [0, dr[2]*exp(-1j*dr[0])]])
        return JonesVector(dot(W0, self.toArray()))

    def intensity(self):
        """ Intensity from electric field vector """
        return real(self.Jx)**2 + real(self.Jy)**2


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


