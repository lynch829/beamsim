#!/usr/bin/env python
#basic imports
import sys, os, time
from math import sqrt
from mpi4py import MPI
import numpy as np
import tables
import scipy
import matplotlib.pyplot as plt

#synergia imports
import synergia
import synergia_workflow
from iota_66_1IO_nll_sc_options import opts

#my script imports
from base_diagnostics import read_bunch
from base_diagnostics import workflow
from base_diagnostics import elliptic_sp
from base_diagnostics import basic_calcs
from elliptic import EllipticBeam6D



#================== Setting up logger and MPI comunicator ============================
#try:
#if True:
# this is the communicator object that will be used for MPI operations
comm = synergia.utils.Commxx()
myrank = comm.get_rank()
mpisize = comm.get_size()
verbose = opts.verbosity>0

logger = synergia.utils.Logger(0)


if myrank == 0:
    print "my rank is 0"
else:
    print "not rank 0"

#=================  Lattice Dictionary ==================================

lattice_dir = '/home/ncook/src/beamsim/synergia/lattices/Iota6-6/' #location of repo for lattices

lattices = {}
lattices['t3_1IO'] = lattice_dir + "lattice_1IO_nll_center.madx"
lattices['t1_1IO'] = '{}lattice_1IO_proton.madx'.format(lattice_dir)

#=================  Setting up the lattice ==================================

lattice = synergia.lattice.MadX_reader().get_lattice("iota", lattices['t3_1IO'])
lattice_length = lattice.get_length()

#f_quads, d_quads = get_fd_quads(lattice)
#print "There are ", len(f_quads), " focussing quadrupoles"
#print "There are ", len(d_quads), " defocussing quadrupoles"


reference_particle = lattice.get_reference_particle()
mass_from_beam_statement = reference_particle.get_four_momentum().get_mass()


print >>logger, "reference particle mass: ", mass_from_beam_statement
print >>logger, "beam statement mass - electron mass:", mass_from_beam_statement-synergia.foundation.pconstants.me


energy = reference_particle.get_total_energy()
beta = reference_particle.get_beta()
gamma = reference_particle.get_gamma()

# calculate RF frequency.  The exact cavity frequency
# will be set automatically using the harmonic number and the
# orbit length, but it needs to be set in the element so
# that the lattice simulator knows the size of the RF bucket.
# find the harmonic number from the first RF cavity
h = None
for elem in lattice.get_elements():
    if elem.get_type() == "rfcavity":
        if elem.has_double_attribute("harmon"):
            h = elem.get_double_attribute("harmon")
            break
        else:
            raise RuntimeError, "The first RF cavity doesn't have the harmon attribute set"

# MAD language freq attribute is in MHz.
if h:   #only adjust if we've determined a harmonic for the cavity
    freq = h * beta * synergia.foundation.pconstants.c*1.0e-6/lattice_length
    for elem in lattice.get_elements():
        if elem.get_type() == "rfcavity":
            elem.set_double_attribute("freq", freq)
else:
    freq = 0

print >>logger, "lattice length: ", lattice_length
print >>logger, "energy: ", energy
print >>logger, "beta: ", beta
print >>logger, "gamma: ", gamma
print >>logger, "RF frequency: ", freq, " MHz"


# Set the same aperture radius for all elements
for elem in lattice.get_elements():
    
    #set an aperture if specified
    if opts.radius:
        elem.set_double_attribute("circular_aperture_radius", opts.radius)

    # set the propagation method.
    # extractor_type may be "chef_maps", "chef_mixed", or chef_propagate
    # Don't modify for polynomial maps
    if opts.use_maps == "taylor":
        pass
    elif opts.use_maps == "all":
        elem.set_string_attribute("extractor_type", "chef_map")
    elif opts.use_maps == "none":
        elem.set_string_attribute("extractor_type", "chef_propagate")
    elif opts.use_maps == "nonrf":
        elem.set_string_attribute("extractor_type", "chef_mixed")
    elif opts.use_maps == "onlyrf":
        if elem.get_type() == "rfcavity":
            elem.set_string_attribute("extractor_type", "chef_map")
        else:
            elem.set_string_attribute("extractor_type", "chef_propagate")
    else:
        raise RuntimeError, "bad options for use_maps: %d"%opts.use_maps

print >>logger, "use maps for: ", opts.use_maps


#====================== make space charge solver if needed ==================
requested_stepper = opts.stepper
if opts.spacecharge:
    
    solver = opts.solver

    # space charge only works with the split operator stepper, or soelements 
    if (requested_stepper != "splitoperator") and (requested_stepper != "soelements"):
        requested_stepper = "soelements"
        print >>logger, "requested stepper changed to soelements for space charge"

    #force these for test run
    gridx = 64
    gridy = 64
    gridz = 1
    grid = [gridx, gridy, gridz]

    print >>logger, "grid: ", grid

    if opts.comm_divide:
        sc_comm = synergia.utils.Commxx_divider(opts.comm_divide, False)
    else:
        sc_comm = synergia.utils.Commxx(True)

    #sc_comm = synergia.utils.Commxx(True)
    if solver == "2dopen-hockney":
        coll_operator = synergia.collective.Space_charge_2d_open_hockney(sc_comm, grid)
    elif solver == "3dopen-hockney":
        # full signature for 3d_open_hockney constructor is
        # comm, grid, long_kicks, z_periodic, period, grid_entire_period,
        # nsigma

        coll_operator = synergia.collective.Space_charge_3d_open_hockney(sc_comm, grid, opts.long_kicks, False, 0.0, False, opts.nsigma)
    elif solver == "2dbassetti-erskine":
        coll_operator = synergia.collective.Space_charge_2d_bassetti_erskine()
    else:
        raise RuntimeError, "requested space charge operator %s invalid.  Must be either 2dopen-hockney or 3dopen-hockney"%opts.solver

    print >>logger, "Using space charge solver ", solver
    print >>logger, "Grid: ", gridx, " x ", gridy, " x ", gridz

else:
    coll_operator = synergia.simulation.Dummy_collective_operator("stub")
    print >>logger, "No space charge solver used"

#==================== Finished space charge solver ======================


#==================== set up the stepper ===============================


turns = 200
order = 1
nsteps_per_element = opts.steps_per_element
n_ppc = 100 #n_ppc particles per transverse cell
n_macro = n_ppc*opts.gridx*opts.gridy 
dpop = 0.00
emittances = [9.74e-6]
tval = 0.4
cval = 0.01
name = 't3_1IO_NLL'
outputdir = "{}_turns_order{}_{}_{}PPC-14mA".format(turns,order,name,n_ppc)

nsteps = len(lattice.get_elements())*nsteps_per_element
opts.output_dir = outputdir
opts.relpath = opts.output_dir
opts.macro_particles = n_macro
opts.steps = nsteps
opts.steps_per_element = nsteps_per_element
workflow.make_path(outputdir)

print >>logger, "output directory:", opts.output_dir


requested_stepper = opts.stepper
print >>logger, "requested_stepper: ",  requested_stepper


if requested_stepper == "splitoperator":

    print >>logger, "Using split-operator stepper with ", opts.steps, " steps/turn"

    stepper = synergia.simulation.Split_operator_stepper(
        lattice, opts.map_order, coll_operator, opts.steps)

elif requested_stepper == "soelements":

    print >>logger, "Using split-operator stepper elements with ", opts.steps_per_element, " steps/element"

    stepper = synergia.simulation.Split_operator_stepper_elements(
        lattice, opts.map_order, coll_operator, opts.steps_per_element)

elif requested_stepper == "independent":

    print >>logger, "Using independent-operator stepper with ", opts.steps, " steps/turn"

    stepper = synergia.simulation.Independent_stepper(
        lattice, opts.map_order, opts.steps)

elif requested_stepper == "elements":

    print >>logger, "Using step-by-elements-operator stepper with ", opts.steps_per_element, " steps/element"

    stepper = synergia.simulation.Independent_stepper_elements(
        lattice, opts.map_order, opts.steps_per_element)

else:
    raise RuntimeError, "stepper %s invalid, must be either 'splitoperator', 'independent' or 'elements'"%requested_stepper

#================== finished setting up the stepper ========================



#=================  Setting up the lattice simulator ==================================


lattice_simulator = stepper.get_lattice_simulator()

opts.xtune = None
if opts.xtune:
    (orig_xtune, orig_ytune) = lattice_simulator.get_both_tunes()
    new_xtune = opts.xtune
    print >>logger, "Adjusting x tune to: ", new_xtune

    lattice_simulator.adjust_tunes(new_xtune, orig_ytune, f_quads, d_quads,
                                   opts.tune_tolerance)

    # after magnets have been changed I need to update the lattice simulator
    lattice_simulator.update()

#bucket_length = lattice_simulator.get_bucket_length()
#compaction_factor = lattice_simulator.get_momentum_compaction()
#slip_factor = lattice_simulator.get_slip_factor()
#print >>logger, "bucket length: ", bucket_length
#print >>logger, "compaction factor (alfa): ", compaction_factor
#print >>logger, "slip factor: ", slip_factor

opts.xml_save_lattice = None
if opts.xml_save_lattice:
    synergia.lattice.xml_save_lattice(lattice, "iota-ring.xml")


otmap = lattice_simulator.get_linear_one_turn_map()
if verbose:
    print >>logger, "one turn map from synergia2.5 infrastructure"
    print >>logger, np.array2string(otmap, max_line_width=200)

[l, v] = np.linalg.eig(otmap)

if verbose:
    print >>logger, "eigenvalues: "
    for z in l:
        print >>logger, "|z|: ", abs(z), " z: ", z, " tune: ", np.log(z).imag/(2.0*np.pi)

#[ax, bx, qx] = map2twiss(otmap[0:2,0:2])
#[ay, by, qy] = map2twiss(otmap[2:4, 2:4])
#[az, b_cdt, qz] = map2twiss(otmap[4:6,4:6])

#if verbose:
#    print >>logger, "Lattice parameters (assuming uncoupled map)"
#    print >>logger, "alpha_x: ", ax, " alpha_y: ", ay
#    print >>logger, "beta_x: ", bx, " beta_y: ", by
#    print >>logger, "q_x: ", qx, " q_y: ", qy, "q_s: ", qz
#    print >>logger, "beta_cdt: ", b_cdt, " beta_longitudinal: ", b_cdt*beta

#=================  Bunch Options Adjustments  ==================================


opts.t = tval
opts.c = cval
opts.new_tune = 0.3
opts.lnll = 1.8
opts.nseg = 20
vals = basic_calcs.get_base_nll(opts.nseg, opts.lnll, opts.new_tune, opts.t, opts.c)

#specify vals for center of the section
opts.betae = vals[3]
opts.alphae = 0 #fixed 0 alpha for center
opts.beta0 = vals[3]

opts.emits = emittances
opts.lattice = lattice
opts.dpop = dpop
opts.emitx = opts.emits[0]
opts.emity = opts.emits[0]



#================= Either generate or read in a bunch ===================================

if myrank == 0:
    #construct a bunch and make sure to add longitudinal momentum variation
    particles = EllipticBeam6D.toyEllipticalBeam6D(opts)

    for index in range(len(opts.emits)):
        bunch = particles[index]
        initialH,initialI = elliptic_sp.calc_bunch_H(bunch,opts)
        bunch_mean = np.mean(initialH)
        bunch_std = np.std(initialH)
        bunch_var = (bunch_std/bunch_mean)*100
        print "Constructed bunch with {} macroparticles, having mean H: {} and std: {}%".format(opts.macro_particles, bunch_mean,bunch_var)
        #now add longitudinal momentum variation
        #For random samples with mean = 0, sigma = sigma, use sigma*np.random.randn(...)
        #bunch[:,5] = opts.dpop*np.random.randn(1,len(bunch))
        bunch[:,5] = np.zeros(len(bunch)) #0 dpop

    print
    opts.num_total_particles = opts.macro_particles*len(opts.emits)
    opts.tracked_particles = opts.num_total_particles

    np.savetxt('myBunch.txt',bunch)         #write the bunch to a text file

    #comm = synergia.utils.Commxx(True)



#Calculate number of particles based on space charge considerations
current = 14.e-3 #mA of current 
bunch_length = 20.0 #effective bunch length

rp_perlength = current/(beta*scipy.constants.c*scipy.constants.e)
real_particles = rp_perlength*bunch_length
#real_particles = np.round(current*bunch_length/(beta*scipy.constants.c*scipy.constants.e))

#current = (14./400)*1.e-3 #mA of current
#blt = bunch_length/(beta*scipy.constants.c) #temporal length of bunch
#real_particles = np.round((current*blt)/scipy.constants.e)
opts.real_particles = real_particles
print "Approximating a current of {}A using {} particles".format(current,opts.real_particles)
#totalQ = opts.real_particles*scipy.constants.e

bucket_length = beta*lattice.get_length()/4 #RF harmonic number is 4


particles_file = 'myBunch.txt'
myBunch = read_bunch.read_bunch(particles_file, reference_particle, opts.real_particles, bucket_length, comm)

#============================ Finished acquiring a bunch =====================



#================== set up the output diagnostics by creating the bunch_simulator =============

bunch_simulator = synergia.simulation.Bunch_simulator(myBunch)

#basic diagnostics - PER STEP
basicdiag = synergia.bunch.Diagnostics_basic("basic.h5", opts.output_dir)
bunch_simulator.add_per_step(basicdiag)
print >>logger, "saving basic diagnostics each step"

#include full diagnostics
fulldiag = synergia.bunch.Diagnostics_full2("full.h5", opts.output_dir)
bunch_simulator.add_per_turn(fulldiag)
print >>logger, "saving full2 diagnostics each turn"

#tracking diagnostics - PER STEP
#diagnostics = synergia.bunch.Diagnostics_track("track.h5", opts.tracked_particles, opts.output_dir)
#Track a particle twice per turn (e.g. every ~ num_steps/2 add the diagnostic)
#bunch_simulator.add_per_step(diagnostics)
#print >>logger, "saving particle track data each step for ", opts.tracked_particles, " particles"

#particle diagnostics - PER TURN
opts.turnsPerDiag = 1
particlediag = synergia.bunch.Diagnostics_particles("particles.h5",0,0,opts.output_dir)
bunch_simulator.add_per_turn(particlediag, opts.turnsPerDiag)
print >>logger, "saving turn-by-turn particle data every {} turns".format(opts.turnsPerDiag)


#=========== finished setting up diagnostics ===========================================================

#=========== put it all together in the propagator and go!

#import time
opts.turns = 200
    
opts.checkpointperiod = 10
opts.maxturns = opts.turns+1

print "setting up propagator for rank {}".format(myrank)

propagator = synergia.simulation.Propagator(stepper)
propagator.set_checkpoint_period(opts.checkpointperiod)
propagator.set_concurrent_io(opts.concurrent_io)

print "starting simulation for rank {}".format(myrank)
if myrank == 0:
    t_start = time.time()

propagator.propagate(bunch_simulator,opts.turns, opts.maxturns,opts.verbosity)

if myrank == 0:
    t_end = time.time()
    elapsed = t_end - t_start
    print "{} turns took {} seconds, or {} seconds per turn, with a grid size of {}x{}x{}.".format(opts.turns,elapsed, elapsed*1.0/opts.turns, gridx,gridy,gridz)

if myrank == 0:
    #clean up files
    workflow.cleanup(opts.output_dir)

#print "{} turns took {} seconds, or {} seconds per turn.".format(opts.turns,elapsed, elapsed*1.0/opts.turns)


#except Exception, e:
#    sys.stderr.write(str(e) + '\n')
#    MPI.COMM_WORLD.Abort(777)