#!/usr/bin/env python
""" Testing the optics modules """
from physicspy.optics import *
from numpy import random, allclose
import unittest

class TestSequenceFunctions(unittest.TestCase):
    
    def testRandomJonesNormalization(self):
        """ Normalize random Jones vector """
        jones = JonesVector((1-random.rand())+(1-random.rand())*1j,(1-random.rand())++(1-random.rand())*1j)
        #~ print self.jones.toArray()
        jones = jones.normalize()
        self.failUnless(allclose(jones.size(),1.0), 'wrong Jones vector normalization')

    def testToArray(self):
        """ Convert random Jones vector to array """
        x = (1-random.rand())-(1-random.rand())*1j
        y = (1-random.rand())+(1-random.rand())*1j
        jones = JonesVector(x,y)
        self.failUnless(allclose(jones.toArray(),array([x,y])), 'wrong Jones vector to array conversion')

if __name__ == '__main__':
    unittest.main()
