import sys, os
import synergia
from mpi4py import MPI
import synergia_workflow

#Add options

def make_opts(name, order, outputdir, steps, steps_per_element):
    '''A quick function for defining a Synergia options object for map propagator comparisons.'''

    opts = synergia_workflow.Options(name)
    opts.add("map_order", order, "Map order", int)
    opts.add("output_dir",outputdir,"Directory for output files", str)
    
    #opts.add("map_order", 1, "Map order", int)
    opts.add("steps", steps, "Number of steps per turn", int)
    opts.add("steps_per_element",steps_per_element,"Number of steps per element", int)


    opts.add("verbosity", 1, "Verbosity of propagation", int)
    opts.add("turns", 1000, "Number of turns", int)
    opts.add("maxturns", 2000, "Maximum number of turns to run before checkpointing and quitting", int)
    opts.add("checkpointperiod", 3000, "Number of turns to run between checkpoints", int)


    opts.add("emitx", 2.5e-6, "real sigma Horizontal emittance [m rad]", float)
    opts.add("emity", 2.5e-6, "real sigma Vertical emittance [m rad]", float)
    opts.add("stdz", 10.0, "sigma read z [m]", float)
    opts.add("dpop", 0.0, "Delta-p/p spread", float)

    opts.add("macro_particles", 100, "Number of macro particles", int)
    opts.add("real_particles", 1.0e11, "Number of real particles", float)
    opts.add("tracked_particles", 100, "Number of tracked particles", int)
    opts.add("seed", 349250524, "Pseudorandom number generator seed", int)


    #space charge additions

    opts.add("gridx", 64, "grid points in x for solver", int)
    opts.add("gridy", 64, "grid points in y for solver", int)
    opts.add("gridz", 64, "grid points in z for solver", int)

    #options for controlling chef propagation vs. chef mapping!
    opts.add("use_maps", "all", "use maps for propagation either all, none, onlyrf, nonrf")
    #opts.add("allmaps", False, "Use all maps for propagation", bool)
    #opts.add("stepper", "splitoperator",
    #         "Simulation stepper, either 'independent','elements','splitoperator','soelements'", str)
    
    return opts

def make_path(dirname):
    '''Create a directory with the specified name - avoid race conditions if possible'''

    try:
        os.makedirs(dirname)
    except OSError:
        if os.path.exists(dirname):
            # We are nearly safe
            pass
        else:
            # There was an error on creation, so make sure we know about it
            raise
            
def cleanup(dirname):
    '''Cleanup files after run and move diagnostic outputs to proper directory.
    
    Arguments:
        -dirname: This is the relative path - e.g. full path = pwd + dirname
    
    '''
    curdir = os.getcwd()
    newdir = os.getcwd() + dirname
    
    for filename in os.listdir(curdir):
        if filename.endswith('.h5'):
            try:
                oldfn = '/'.join((curdir,filename))
                newfn = '/'.join((curdir,dirname,filename))
                #don't worry about shutil here, keep simple with os module
                os.rename(oldfn,newfn)
            except OSError:
                if os.path.exists(newfn):
                    #file already exists so delete and re-try
                    os.remove(newfn)
                    os.rename(oldfn,newfn)
                else:
                    #perhaps trying to move to a new disk or something that os can't handle
                    raise