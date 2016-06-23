'''
A class that serves as a generic catalog
'''

################

import unittest
import numpy as np

################

class MismatchedLengthException(Exception): pass

class Catalog(object):

    def __init__(self):
        super(Catalog, self).__setattr__('table', {})
        super(Catalog, self).__setattr__('header', {})
        super(Catalog, self).__setattr__('length', -1)

    def copy(self):

        newcat = Catalog()
        newcat.header.update(self.header)
        newcat.table.update(self.table)
        super(Catalog, newcat).__setattr__('length', self.length)

        return newcat


    def __getattr__(self, name):

        header = Catalog.__getattribute__(self, 'header')
        if name in header:
            return header[name]
        table = Catalog.__getattribute__(self, 'table')
        return table[name]

    def __setattr__(self, name, val):
        
        try:

            if len(val) == self.length:
                self.table[name] = val
            elif self.length == -1:
                self.table[name] = val
                super(Catalog, self).__setattr__('length', len(val))
            else:
                raise MismatchedLengthException

        except TypeError:

            self.header[name] = val


    def __len__(self):
        
        if self.length == -1:
            return 0
        return self.length

    def filter(self, mask):

        newcat = Catalog()
        for key, val in self.header.iteritems():
            newcat.__setattr__(key, val)
        for key, val in self.table.iteritems():
            newcat.__setattr__(key, val[mask])

        return newcat


#########################


class TestCatalog(unittest.TestCase):

    def testSetAttribs(self):

        zcluster = 0.5
        redshifts = np.arange(0.0, 1.0, 0.01)
        betas = np.arange(len(redshifts))

        cat = Catalog()
        cat.clusterz = 0.5
        cat.redshifts = redshifts
        cat.betas = betas

        self.assertEqual(zcluster, cat.clusterz)
        self.assertTrue((redshifts == cat.redshifts).all())
        self.assertTrue((betas ==  cat.betas).all())

    ####

    def testEnforceLength(self):

        cat = Catalog()
        cat.clusterz = 0.5
        cat.redshifts = np.arange(5)
        self.assertEqual(len(cat), 5)
        with self.assertRaises(MismatchedLengthException):
            cat.betas = np.arange(3)

    ####

    def testFilter(self):

        redshifts = np.arange(0.0, 1.0, 0.01)
        betas = np.arange(len(redshifts))

        cat = Catalog()
        cat.clusterz = 0.5
        cat.redshifts = redshifts
        cat.betas = betas
        
        newcat = cat.filter(cat.redshifts < 0.5)
        self.assertEqual(len(redshifts[redshifts < 0.5]), len(newcat))
        self.assertTrue((redshifts[redshifts < 0.5] == newcat.redshifts).all())
        self.assertTrue((betas[redshifts < 0.5] == newcat.betas).all())
        self.assertEqual(cat.clusterz, 0.5)

    ###

    def testCopy(self):

        redshifts = np.arange(0.0, 1.0, 0.01)
        betas = np.arange(len(redshifts))

        cat = Catalog()
        cat.clusterz = 0.5
        cat.redshifts = redshifts
        cat.betas = betas

        newcat = cat.copy()
        newcat.redshifts = np.zeros_like(redshifts)
        newcat.clusterz = 0.3

        self.assertTrue((cat.redshifts == redshifts).all())
        self.assertEqual(cat.clusterz, 0.5)
        
        self.assertTrue((newcat.redshifts == 0).all())
        self.assertTrue(newcat.clusterz == 0.3)
        self.assertTrue((newcat.betas == betas).all())
        



def runtests():

    unittest.main()

if __name__ == '__main__':

    runtests()



            
