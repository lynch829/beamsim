import sys
import os
import tables
import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt
import options
#from mpl_toolkits.axes_grid import make_axes_locatable
#from mpl_toolkits import axes_grid

coords = {}
coords['x'] = 0
coords['xp'] = 1
coords['y'] = 2
coords['yp'] = 3
coords['cdt'] = 4
coords['dpop'] = 5
coords['id'] = 6

#opts.coords = coords


###################################### PLOTTERS ###################################

def plot_SPH(hamArray,opts, ID=0):
    
    '''
    Create plot for single particle hamiltonian for one particle over pre-specified turn #s
    
    Arguments:
    ham - numpy array of single particle hamiltonian values
    opts - options object
    
    Optional:
    ID - specify the specific particle to be plotted (default 0)
    
    '''
    
    #general plot settings
    matplotlib.rcParams['figure.autolayout'] = True
    
    #change ID only if necessary
    if opts.ID:
        ID = opts.ID
   
    h = np.arange(opts.turns+1) #plus 1 to account for initial conditions
    v = hamArray[:,ID]
    vinit = v[0] #take the iniital values as the normalization value
    
    if opts.norm:
        vScale = v/vinit
        ymin = 0
        ymax = 2
    else:
        vScale = v
        ymax = 1.5*vScale.max()
        ymin = 0.5*vScale.min()
    
    fig = plt.figure()
    plt.subplot(1,1,1)
    ax = plt.gca()
    
    #print ymin, ymax
    
    ax.scatter(h,vScale, c ='b', s=6)
    #ax.set_aspect('equal')
    
    ax.set_ylim([ymin,ymax])
    #plt.plot(h,v, 'o')
    plt.xlabel(opts.hcoord,fontsize=12)
    plt.ylabel(opts.vcoord,fontsize=12)
    title = 'Single particle hamiltonian for particle ' + str(ID)
    if not opts.lattice_name== None:
        title = title + ' for lattice ' + opts.lattice_name
    plt.title(title, y=1.05, fontsize=14)
    #plt.draw()
    #fig.tight_layout()
    plt.show()
    
    if opts.save:
        sv_title = 'SPH_0_'+ opts.lattice_name + '.pdf'
        fig.savefig(sv_title, bbox_inches='tight')


def plot_P(PArray, opts, num=10, ID=0):
    
    '''Create a Poincare plot for specified particle IDs over pre-specified turn #s
    
    
    Note: Currently, opts.plots must be a list of coordinates of length 2!
    
    '''

    #general plot settings
    matplotlib.rcParams['figure.autolayout'] = True
    
    if opts.num:
        num = opts.num
    
    #plot up to 10 particle tracks
    if PArray.shape[1] < num:
        num = PArray.shape[1]

    #plot specified # of turns
    if opts.turns:
        turns = opts.turns
    else: turns = PArray.shape[0]        
    
    plots = {'x':0, 'px':1, 'y':2, 'py':3}
    
    cropped = PArray[:turns,:num,(plots[opts.plots[0]],plots[opts.plots[1]])]
    
    reshaped = cropped.reshape(cropped.size/2,2).T
    
    #separate horizontal and vertical components
    h = reshaped[0]
    v = reshaped[1]
    
    
    if opts.scale:
        fig = plt.figure(figsize=(opts.scale*8,opts.scale*6))
    else:
        fig = plt.figure(figsize=(8,6))
        
    plt.subplot(1,1,1)
    ax = plt.gca()
    
    #print ymin, ymax
    
    ax.scatter(h,v, c ='b', s=2)
    ax.set_aspect('equal') #want equal aspect ratios for Poincare plots
    
    #ax.set_ylim([ymin,ymax])
    #plt.plot(h,v, 'o')
    plt.xlabel(opts.hcoord,fontsize=round(12*opts.scale))
    plt.ylabel(opts.vcoord,fontsize=round(12*opts.scale))
    title = opts.hcoord + '-'+ opts.vcoord+' for ' + str(turns) + ' turns'
    if not opts.lattice_name== None:
        title = title + ' for lattice ' + opts.lattice_name
    plt.title(title, y=1+0.05/opts.scale, fontsize=round(14*opts.scale))
    #plt.draw()
    #fig.tight_layout()
    plt.show()
    
    if opts.save:
        sv_title = 'Poincare'+'_' + opts.hcoord+'_' + opts.vcoord+'_'+ str(turns) + '_turns_'+  opts.lattice_name + '.pdf'
        fig.savefig(sv_title, bbox_inches='tight')     
        
 
def plot_J(JArray,opts, ID=0):
    
    '''
    Create plot for single particle invariant for one particle over pre-specified turn #s
    
    Arguments:
    JArray - numpy array of single particle hamiltonian values
    opts - options object
    
    Optional:
    ID - specify the specific particle to be plotted (default 0)
    
    '''
    
    #general plot settings
    matplotlib.rcParams['figure.autolayout'] = True
    
    #change ID only if necessary
    if opts.ID:
        ID = opts.ID
   
    h = np.arange(opts.turns+1) #plus 1 to account for initial conditions
    v = JArray[:,ID]
    vinit = v[0] #take the iniital values as the normalization value
    
    if opts.norm:
        vScale = v/vinit
        ymin = 0
        ymax = 2
    else:
        vScale = v
        if opts.variance:
            ymax = (1+opts.variance)*vScale.max()
            ymin = (1-opts.variance)*vScale.min()
        else:
            ymax = 1.05*vScale.max()
            ymin = 0.95*vScale.min()            
    
    fig = plt.figure(figsize=(8,6))
    plt.subplot(1,1,1)
    ax = plt.gca()
    
    #print ymin, ymax
    
    ax.scatter(h,vScale, c ='b', s=6)
    #ax.set_aspect('equal')
    
    ax.set_ylim([ymin,ymax])
    #plt.plot(h,v, 'o')
    plt.xlabel(opts.hcoord,fontsize=12)
    plt.ylabel(opts.vcoord,fontsize=12)
    title = 'Courant synder invariant for particle ' + str(ID)
    if not opts.lattice_name== None:
        title = title + ' for lattice ' + opts.lattice_name
    plt.title(title, y=1.05, fontsize=14)
    #plt.draw()
    #fig.tight_layout()
    plt.show()
    
    if opts.save:
        sv_title = 'J_'+str(ID)+'_'+ opts.lattice_name + '.pdf'
        fig.savefig(sv_title, bbox_inches='tight') 
 
    
#define an option to replicate the pltbunch.plot_bunch function?

def plot_tracks(tracks, opts, ID=0):
    '''
    
    Plot particle coordinates vs. 's' coordinate for a single particle.
    
    Arguments:
    tracks - array of coordinate values - organized to correspond to index of opts.plots!
    opts - options object
    
    Optional:
    ID - use to specify a particle - default None
    
    '''

    #general plot settings
    matplotlib.rcParams['figure.autolayout'] = True
    
    #change ID only if necessary
    if opts.ID:
        ID = opts.ID

    #if #turns specified, then slice on that as well
    if opts.turns:
        hcoord = np.arange(opts.turns+1) #plus 1 to account for initial conditions
        vcoords = tracks[:opts.turns+1,:,ID].transpose() 
    else:
        hcoord = np.arange(tracks.shape[0]) 
        vcoords = tracks[::,:,ID].transpose()
   
    fig = plt.figure()
    plt.subplot(1,1,1)
    ax = plt.gca()
    
    mins = []
    maxes = []
    
    for index,v in enumerate(vcoords):
        #vinit = v[0] #take the iniital values as the normalization value
        scale = v.max() #take the max value to scale v
        
        if opts.norm:
            vScale = v/scale
            mins.append(0)
            maxes.append(2)
        else:
            vScale = v
            min1 = vScale.min()
            if min1 < 0:
                mins.append(1.5*min1)
            else:
                mins.append(0.5*min1)
            maxes.append(1.5*vScale.max())
        #ax.scatter(hcoord,vScale, c ='b', s=6)
        ax.plot(hcoord,vScale, label=opts.plots[index])
        
    ax.set_ylim([min(mins),max(maxes)])
    #plt.plot(h,v, 'o')
    plt.xlabel(opts.hcoord,fontsize=12)
    
    ylabel = ', '.join(opts.plots)
    plt.ylabel(ylabel,fontsize=12)
    title = 'Particle tracks for particle ' + str(ID)
    if not opts.lattice_name== None:
        title = title + ' for lattice ' + opts.lattice_name
    plt.title(title, y=1.05, fontsize=14)
    #plt.draw()
    #fig.tight_layout()
    plt.legend(loc='best')
    plt.show()
    
    if opts.save:
        sv_title = 'Tracks_' + str(ID) + '_' + opts.lattice_name + '.pdf'
        fig.savefig(sv_title, bbox_inches='tight')
    



################################## GETTERS ####################################

def get_particles(inputfile):
    
    '''Reads an input file and returns a numpy array of particles and a dictionary of root values'''
    
    f = tables.openFile(inputfile, 'r')
    particles = f.root.particles.read()
    
    #get appropriate reference properties from file root
    npart = particles.shape[0]
    mass = f.root.mass[()]
    p_ref = f.root.pz[()]
    sn = f.root.s_n[()] #period length
    tn = f.root.tlen[()] #cumulative tracked length for this file

    f.close()
    
    header = dict()
    header['n_part'] = npart
    header['mass'] = mass
    header['p_ref'] = p_ref
    header['s_val'] = sn
    header['t_len'] = tn
    
    return header,particles
    

def get_file_list(opts):

    '''
    
    Returns a list of files of the form 'particles*.h5' from the current directory 
    
    
    Optional Arguments:
    path - link to directory storing .h5 files - default None (executes within present working directory)
    
    '''
    
    if not opts.relpath:
        #If no relative path specified, check current directory 
        files = os.listdir(os.getcwd())
    else:
        #Otherwise check specified relative path
        path = os.path.join(os.getcwd(),opts.relpath)
        files = [os.path.join(path,fn) for fn in os.listdir(path)]
        #files = os.listdir(opts.path) #only change directories if a different path is specified

    pfiles = []

    #basic filtering for particles*.h5 files
    for filename in files:
            if filename.find('particles') > -1 and filename.endswith('.h5'):
                pfiles.append(filename)
    
    return pfiles
    
def get_twiss(lattice_simulator):
    '''
    Returns an array of twiss parameters versus longitudinal coordinate 's' for a given lattice.
    
    Arguments:
    lattice_simulator - a Synergia lattice simulator
    
    Return values have array configuration: [s,betax,alphax,gammax,betay,alphay,gammay]
    '''
    
    lattice = lattice_simulator.get_lattice()
    
    twiss = []
    
    for elem in lattice.get_elements():
        temp = []
        
        lattice_functions=lattice_simulator.get_lattice_functions(elem)
        gamma_x = (1 + lattice_functions.alpha_x**2)/lattice_functions.beta_x
        gamma_y = (1 + lattice_functions.alpha_y**2)/lattice_functions.beta_y
        
        temp = [lattice_functions.arc_length, lattice_functions.beta_x, lattice_functions.alpha_x, gamma_x, lattice_functions.beta_y, lattice_functions.alpha_y, gamma_y]
        twiss.append(temp)
    
    #at the end, pre-pend final twiss values @ beginning with s-value 0
    twiss_init = twiss[-1]
    twiss_init[0]=0.0
    twiss.insert(0,twiss_init)
        
    return np.asarray(twiss)
    

def normalized_coordinates(header, particles, twiss, units=None, ID=None):
    '''Return the an array of particles (fixed s) in normalized transverse coordinates rather than trace space
    
    Input Coordinates - x, x', y, y'
    Output Coordinates - x, beta*x' + alpha*x, y, beta*y + alpha*y
    
    '''
    
    sval = header['s_val']
    svals = twiss[::,0]
    betaxvals = twiss[::,1]
    alphaxvals = twiss[::,2]
    gammaxvals = twiss[::,3]
    betayvals = twiss[::,4]
    alphayvals = twiss[::,5]
    gammayvals = twiss[::,6]
    
    #interpolate if needed
    if not sval in svals:
        betax = np.interp(sval, svals, betaxvals)
        betay = np.interp(sval, svals, betayvals)
        alphax = np.interp(sval, svals, alphaxvals)
        alphay = np.interp(sval, svals, alphayvals)
        gammax = np.interp(sval, svals, gammaxvals)
        gammay = np.interp(sval, svals, gammayvals)
    else:
        ind = list(svals).index(sval)
        betax = twiss[ind,1]
        alphax = twiss[ind,2]
        gammax = twiss[ind,3]
        betay = twiss[ind,4]
        alphay = twiss[ind,5]
        gammay = twiss[ind,6]
    
    
    x = particles[:,coords['x']] #units m
    newx = x / math.sqrt(betax) #normalized
    
    xp = particles[:,coords['xp']] #unitless
    px = (betax*xp + alphax*x)/math.sqrt(betax) #normalized
    
    y = particles[:,coords['y']] #units m
    newy = y / math.sqrt(betay) #normalized
    
    yp = particles[:,coords['yp']] #unitless
    py = (betay*yp + alphay*y) /math.sqrt(betay) #normalized    
    
    #stack arrays then transpose to get array of coordinate vectors
    particles_norm = np.vstack((newx,px,newy,py)).T
    
    return particles_norm

def get_normalized_coords(filelist, twiss, num=None, ID=None):
    
    '''
    
    Returns a numpy array of normalized coordinate vectors obtained from .h5 files in filelist
    
    Arguments:
    filelist - A list of .h5 files containing particle array information
    twiss - array of twiss parameters for the lattice
    
    
    Returns a numpy array with dimensions #turns x #particles x #transverse coordinates(4).
    norms[B][A] returns vector of coordinates for particle A at turn B.
    '''

    #opts = options.Options()
    #opts.hcoord = 'x'
    #opts.vcoord = 'xp'
    #opts.lattice_name = 'FODO'
    norms = [] #norms is a list of arrays of particle coordinates
    
    for index,fileName in enumerate(filelist):
        inputfile = fileName
        header, particles = get_particles(inputfile)
        norm_coords = normalized_coordinates(header, particles, twiss)
        #only append first num of tracked particles if not plotting all
        if num:
            norms.append(norm_coords[:num])
        else:  
            norms.append(norm_coords)
    return np.asarray(norms)


def single_particle_invariant(header, particles, twiss, units=None, ID=None):
    
    '''
    
    Returns an array of single particle invariants (Courant Synder)
    
    Arguments:
    header - a header dictionary obtained from 'get_particles()'
    particles - an array of particles obtained from .h5 files via f.root.particles.read()
    twiss - an array of twiss functions
    
    Optional:
    ID - use to specify values for a single particle - default None
    
    '''
    
    sval = header['s_val']
    svals = twiss[::,0]
    betaxvals = twiss[::,1]
    alphaxvals = twiss[::,2]
    gammaxvals = twiss[::,3]
    betayvals = twiss[::,4]
    alphayvals = twiss[::,5]
    gammayvals = twiss[::,6]
    
    #interpolate if needed
    if not sval in svals:
        betax = np.interp(sval, svals, betaxvals)
        betay = np.interp(sval, svals, betayvals)
        alphax = np.interp(sval, svals, alphaxvals)
        alphay = np.interp(sval, svals, alphayvals)
        gammax = np.interp(sval, svals, gammaxvals)
        gammay = np.interp(sval, svals, gammayvals)
    else:
        ind = list(svals).index(sval)
        betax = twiss[ind,1]
        alphax = twiss[ind,2]
        gammax = twiss[ind,3]
        betay = twiss[ind,4]
        alphay = twiss[ind,5]
        gammay = twiss[ind,6]
    
    
    x = particles[:,coords['x']] #units m
    xp = particles[:,coords['xp']] #unitless
    
    y = particles[:,coords['y']] #units m
    yp = particles[:,coords['yp']] #unitless     
    
    invariantx = gammax*x**2 + 2*alphax*x*xp + betax*xp**2
    invarianty = gammay*y**2 + 2*alphay*y*yp + betay*yp**2
    
    #inv = invariantx**2+invarianty**2
    inv2 = invariantx + invarianty
    
    return inv2
    
def get_invariants(filelist, twiss):
    
    '''
    
    Returns a numpy array of single particle invariants (Courant Synder) values obtained from .h5 files in filelist
    
    Arguments:
    filelist - A list of .h5 files containing particle array information
    twiss - array of twiss parameters for the lattice
    
    '''

    #opts = options.Options()
    #opts.hcoord = 'x'
    #opts.vcoord = 'xp'
    #opts.lattice_name = 'FODO'
    invariant = [] #hamiltonian is a list of arrays of macroparticles
    
    for index,fileName in enumerate(filelist):
        inputfile = fileName
        header, particles = get_particles(inputfile)
        invariant.append(single_particle_invariant(header, particles, twiss))
        
    return np.asarray(invariant)


def single_particle_hamiltonian(header, particles, units=None, ID=None):
    
    '''
    
    DEPRECATED: Returns an array of single particle Hamiltonian values for the particles described by 'particles' 
    
    Arguments:
    header - a header dictionary obtained from 'get_particles()'
    particles - an array of particles obtained from .h5 files via f.root.particles.read()
    
    Optional:
    ID - use to specify values for a single particle - default None
    
    '''
    
    pref = header['p_ref'] #reference momentum in GeV/c
    mass = header['mass'] #mass in GeV/c^2
    gevc = 5.34428576e-19 #mkg/s per GeV/c
    
    x = particles[:,coords['x']] #units m
    px = particles[:,coords['xp']] #unitless
    
    y = particles[:,coords['y']] #units m
    py = particles[:,coords['yp']] #unitless 
    
    if units:
        print "Units flag specified"
        px = px*pref*gevc #units kgm/s
        py = py*pref*gevc #units kgm/s
    
    if ID:
        #quadratic is the basis for calculating the hamiltonian, in absence of nonlinear couplings
        quadratic = 0.5* (px[ID]**2 + py[ID]**2) + 0.5*(x[ID]**2 + y[ID]**2)
    else:
        quadratic = 0.5* (px**2 + py**2) + 0.5*(x**2 + y**2)

    return quadratic
        
    

def get_hamiltonians(filelist):
    
    '''
    
    Returns a numpy array of single particle hamiltonian values obtained from .h5 files in filelist
    
    Arguments:
    filelist - A list of .h5 files containing particle array information
    
    '''

    #opts = options.Options()
    #opts.hcoord = 'x'
    #opts.vcoord = 'xp'
    #opts.lattice_name = 'FODO'
    hamiltonian = [] #hamiltonian is a list of arrays of macroparticles
    
    for index,fileName in enumerate(filelist):
        inputfile = fileName
        header, particles = get_particles(inputfile)
        hamiltonian.append(single_particle_hamiltonian(header, particles))
        
    return np.asarray(hamiltonian)

def get_tracks(filelist, opts):
    
    '''
    
    Returns a numpy array of single particle coordinate values obtained from .h5 files in filelist
    
    Arguments:
    filelist - A list of .h5 files containing particle array information
    
    '''

    #opts = options.Options()
    #opts.hcoord = 'x'
    #opts.vcoord = 'xp'
    #opts.lattice_name = 'FODO'
    tracks = [] #hamiltonian is a list of arrays of macroparticles
    
    for index,fileName in enumerate(filelist):
        inputfile = fileName
        header, particles = get_particles(inputfile)
        vec = []
        for coord in opts.plots:
            #make sure specified plot options are obtainable
            assert coord in coords.keys(), "Specified plot, %s is not available from: %s" %(coord, coords.keys())   
            vec.append(particles[:,coords[coord]])
        tracks.append(vec)
        
    return np.asarray(tracks)    


################################################################################################


def plot_Poincare(opts):
    
    '''Plot a poincare section in the desired normalized coordinates'''
    
    opts.hcoord = opts.plots[0]
    opts.vcoord = opts.plots[1]
    
    files = get_file_list(opts)
    twiss = get_twiss(opts.lattice_simulator)
    pArray = get_normalized_coords(files,twiss)
    plot_P(pArray, opts) 


def plot_Invariant(opts):
    
    '''
    
    Plots the single particle hamiltonian over a specified # of turns and # of particles
    
    Arguments:
    opts - an Options object specifying # of turns, particle #s, etc.
    
    '''
    
    opts.hcoord = 'turn #'
    opts.vcoord = 'J(p,q)'
    
    files = get_file_list(opts)
    twiss = get_twiss(opts.lattice_simulator)
    jArray = get_invariants(files, twiss)
    #return jArray
    plot_J(jArray,opts)


def plot_Hamiltonian(opts):
    
    '''
    
    Plots the single particle hamiltonian over a specified # of turns and # of particles
    
    Arguments:
    opts - an Options object specifying # of turns, particle #s, etc.
    
    '''
    
    opts.hcoord = 's'
    opts.vcoord = 'H(p,q)'
    
    files = get_file_list(opts)
    hamArray = get_hamiltonians(files)
    plot_SPH(hamArray,opts)
    
    
def track(opts):
    
    '''
    Plots specified coordinates for a single particle versus longitudinal coordinate `s`.
    
    Arguments:
    opts - an Options object specifying # of coordinates, turns, particle #s, etc.
    
    '''
    
    opts.hcoord = 'turn #'
    files = get_file_list(opts)
    tracks = get_tracks(files, opts)
    plot_tracks(tracks, opts)