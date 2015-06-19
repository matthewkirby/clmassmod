#######################
# Publication quality plots
#######################


import pylab
import numpy as np


golden_mean = (np.sqrt(5) - 1)/2.0
screen_ratio = 0.75
fig_width = 6
fig_height = fig_width*screen_ratio
fig_size = [fig_width,fig_height]


params = {'text.usetex' : True,
          'ps.usedistiller' : 'xpdf',
          'ps.distiller.res' : 6000,
          'axes.labelsize' : 16,
          'text.fontsize' : 16,
          'legend.fontsize' : 14,
          'xtick.labelsize' : 14,
          'ytick.labelsize' : 14,
          'figure.figsize' : fig_size}
pylab.rcParams.update(params)


colors = [(.9,.6,0), (.35, .7, .9), (0,.6,.5), (0.95, 0.9, 0.25)]
