import numpy as np
#import deconvolvedlognorm as dln

#############



def createData(datasig = .15):

    truth1 = -0.05 + 0.2*np.random.standard_normal(size=100)
    truth2 = 0.1 + 0.3*np.random.standard_normal(size=45)
    truth3 = -0.02 + 0.01*np.random.standard_normal(size=8)

    true_deviates = np.hstack([truth1, truth2, truth3])

    truth=5e14

    truevals = truth + true_deviates*1e15

    masses = np.arange(-1e15, 5e15, 1e13)

    sig = datasig*1e15

    statvals = truevals + sig*np.random.standard_normal(len(truevals))

    statval_grid, mass_grid = np.meshgrid(statvals, masses, indexing='ij')
    
    probs = np.exp(-0.5*((mass_grid - statval_grid)/sig)**2)/(np.sqrt(2*np.pi)*sig)

    
    halos = [dict(id = i,
                  true_mass = truth,
                  masses = masses,
                  pdf = probs[i,:]) for i in range(len(truevals))]

    

    return halos
    

def runTest():

    halos = createData()

    parts = dln.buildGaussMixture1DModel(halos, 3)

    dln.sample(parts, 'testgaussmix', 100, singlecore=True)


if __name__ == '__main__':

    runTest()
