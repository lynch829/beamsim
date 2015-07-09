import sys
import synergia
import matplotlib.pyplot as plt
from matplotlib import gridspec

def plot_bunch(bunch):
    
    '''
    Plot the coordinate space and x, y phase spaces of a synergia bunch
    '''
    
    #get particle properties from bunch
    myParticles = bunch.get_local_particles()
    xp = myParticles[:,bunch.xp]
    x = myParticles[:,bunch.x]
    yp = myParticles[:,bunch.yp]
    y = myParticles[:,bunch.y]
    
    #one way to use subplots
    #fig, (ax0, ax1, ax2)  = plt.subplots(ncols=3, figsize=(10,6))

    #another way - use gridspec
    fig = plt.figure(figsize=(9.9,3.3))
    gs = gridspec.GridSpec(1, 3, width_ratios=[1,1,1]) 
    
    ax0 = plt.subplot(gs[0])
    ax0.scatter(x, y, c='b')
    ax0.set_title('X-Y Coordinate Space')
    #ax0.set_aspect(aspect=2.0)
    
    ax1 = plt.subplot(gs[1])
    ax1.scatter(x, xp, c='r')
    ax1.set_title('X Trace Space')
    #ax1.set_aspect(aspect=2.0)
    
    ax2 = plt.subplot(gs[2])
    ax2.scatter(y, yp, c='g')
    ax2.set_title('Y Trace Space')
    #ax2.set_aspect(aspect=2.0)
    
    # Tweak spacing between subplots to prevent labels from overlapping
    #plt.subplots_adjust(hspace=2)
    
    #set figure title
    fig.canvas.set_window_title('Synergia Phase Space Distribution')
    fig.tight_layout()
    plt.show()
    
def plot_long(bunch):
    '''
    Plot the longitudinal coordinate and phase space of the bunch
    '''
    
    #get particle properties from bunch
    myParticles = bunch.get_local_particles()
    x = myParticles[:,bunch.x]
    y = myParticles[:,bunch.y]
    z = myParticles[:,bunch.z]
    zp = myParticles[:,bunch.zp]
    
    #one way to use subplots
    #fig, (ax0, ax1, ax2)  = plt.subplots(ncols=3, figsize=(10,6))
    
    #another way - use gridspec
    fig = plt.figure(figsize=(9.9,3.3))
    gs = gridspec.GridSpec(1, 3, width_ratios=[1,1,1]) 
    
    ax0 = plt.subplot(gs[0])
    ax0.scatter(x, y, c='b')
    ax0.set_title('X-Y Coordinate Space')
    #ax0.set_aspect(aspect=2.0)
    
    ax1 = plt.subplot(gs[1])
    ax1.scatter(z, zp, c='r')
    ax1.set_title('Z Trace Space')
    #ax1.set_aspect(aspect=2.0)
    
    ax2 = plt.subplot(gs[2])
    ax2.scatter(z, y, c='g')
    ax2.set_title('Z-Y coordinate space')
    #ax2.set_aspect(aspect=2.0)
    
    # Tweak spacing between subplots to prevent labels from overlapping
    #plt.subplots_adjust(hspace=2)
    
    #set figure title
    fig.canvas.set_window_title('Synergia Phase Space Distribution')
    fig.tight_layout()
    plt.show()