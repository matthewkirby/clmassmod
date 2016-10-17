###################
# Attempt to analytically model N-body sim results for WL bias
###################



######
# Imports
######

import numpy as np
import unittest

import colossusMassCon as cmc


#######

def generateDiemerShearProfile(m200, z, radii):

    return 0.



######
# Testing
######

class TestAnalyticBiasEstimate(unittest.TestCase):

    def test_generateDiemerShearProfile(self):

        radii = np.arange(0.3, 1.5, 0.1)

        z=0.25
        m200 = 1e14

        g_profile = generateDiemerShearProfile(m200, z, radii) 

        self.assertIsInstance(g_profile, type(np.ones(5)))
        self.assertEqual(g_profile.shape, radii.shape)





def runTests():

    suite = unittest.TestLoader().loadTestsFromTestCase(TestAnalyticBiasEstimate)
    unittest.TextTestRunner(verbosity=2).run(suite)


if __name__ == '__main__':

    runTests()
