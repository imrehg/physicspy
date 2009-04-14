#!/usr/bin/env python
""" Testing the atomic module """

from physicspy.quantum.atomic import *
from numpy import random, allclose
import unittest

class TestSequenceFunctions(unittest.TestCase):

    def testThreeJTriangle1(self):
        """ Triangle rule test (lower bound for J[0])  """
        J = ([1,1,3])
        M = ([1,0,-1])
        self.failUnless(threej(J[0],J[1],J[2],M[0],M[1],M[2]) == 0.0, 'three-j symbol vector triangularity check fail (lower)')

    def testThreeJTriangle2(self):
        """ Triangle rule test (upper bound for J[0]) """
        J = ([4,1,2])
        M = ([1,0,-1])
        self.failUnless(threej(J[0],J[1],J[2],M[0],M[1],M[2]) == 0.0, 'three-j symbol vector triangularity check fail (upper)')

    def testThreeJLimitM(self):
        """ Limiting value of |M| < J """
        J = ([1,1,2])
        M = ([2,0,-1])
        self.failUnless(threej(J[0],J[1],J[2],M[0],M[1],M[2]) == 0.0, 'three-j symbol m-limit check fail')

    def testThreeJM(self):
        """ Sum(M) = 0 test"""
        J = ([1,1,2])
        M = ([0,0,-1])
        self.failUnless(threej(J[0],J[1],J[2],M[0],M[1],M[2]) == 0.0, 'three-j symbol m-sum check fail')

    def testThreeJTriangle2(self):
        """ Three-J symbol value test by calculating Clebsch-Gordan coefficient """
        J = ([1,1,2])
        M = ([0,1,-1])
        self.failUnless(allclose( (2*J[0]+1)*threej(J[0],J[1],J[2],M[0],M[1],M[2])**2 , 0.3), 'wrong three-j symbol value')

if __name__ == '__main__':
    unittest.main()
