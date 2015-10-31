#Utilities for use with python data analysis, specifically for analyzing Synergia and Warp outputs.
#
#Author: Nathan Cook
# 10/28/2015

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


def plot_distribution(arr, numBins, norm=False):
    '''Construct a histrogram, then plot the resulting distribution
    
    Arguments:
        arr - 1D array for binning
        numbins - # of bins for histogram
        
        (optional) norm - normalization flag (default False)
    
    Constructs a matplotlib plot and returns it. 
    Automatically displayed when using IPython backend.
    
    '''
    
    myVals, myBins = np.histogram(arr,numBins)
    bincenters = 0.5*(myBins[1:]+myBins[:-1])
    bin_width = 1.0*(myBins[1]-myBins[0])
    #normalize
    myVals_norm = myVals/(np.max(myVals)*1.0)
    
    #Set some matplotlib standards for plots 
    #mpl.rcParams['font.size': 14]
    mpl.rcParams['figure.figsize'] = 8, 6
    
    fig = plt.figure()
    ax = fig.gca()
    
    ax.set_xlabel("Array values", fontsize=14)
    ax.set_ylabel("Relative population", fontsize=14)
    
    ax.plot(bincenters,myVals_norm, c='k')
    #close the first display call
    plt.close()
    
    return fig #return the figure handle