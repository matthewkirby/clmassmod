import numpy as np
import pylab
import scipy.ndimage

def Confidence2D(histogram, xedges, yedges, problevel, smooth=None):
    
    #first, sort the histogram from high to low
    sortOrder = np.argsort(histogram.flatten())[::-1]
    sumprob = np.cumsum(histogram.flatten()[sortOrder])
    selected = (sumprob/float(np.sum(histogram))) < problevel
    
    CLregion = np.zeros(len(sumprob))
    CLregion[sortOrder[selected]] = 1.
    CLregion = np.resize(CLregion, histogram.shape)
    
    if smooth:
        CLregion = scipy.ndimage.gaussian_filter(CLregion, smooth)
    
    #if we want to plot contours, not meshes, need to create X,Y grids
    
    xcenters = (xedges[1:] + xedges[:-1])/2.
    ycenters = (yedges[1:] + yedges[:-1])/2.
    
    Ygrid, Xgrid = np.meshgrid(ycenters, xcenters)
    
    return CLregion, Xgrid, Ygrid



def CreateSigmaSteps(histogram, xedges, yedges, problevels=[0.68,0.95,0.99], smooth = 0.75, fig = None):

    CLregion, Xgrid, Ygrid = Confidence2D(histogram, xedges,yedges, problevels[0], smooth=0.75)

    CLregion = problevels[0]*CLregion

    for i in range(1, len(problevels)):
        curRegion, X, Y = Confidence2D(histogram, xedges,yedges, problevels[i], smooth=smooth)
        CLregion += problevels[i]*curRegion

    problevels.append(1.)
    problevels = np.array(problevels)
    splits = (problevels[1:] + problevels[0:-1])/2.

    if fig is None:
        fig = pylab.figure()

    ax = pylab.gca()

    ax.contour(Xgrid, Ygrid, CLregion, splits)

    return fig

    
