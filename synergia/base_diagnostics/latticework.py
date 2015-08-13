import sys, os
import synergia
from mpi4py import MPI
import synergia_workflow
import numpy as np

def make_chef(lattice, lattice_simulator):
    '''Make all elements of a lattice chef-propogated and update the lattice simulator.'''
    
    elems = lattice.get_elements()
    #nlelems = []
    for elem in elems:
        if elem.get_type() == 'nllens':
            #nlelems.append(elem)
            elem.set_string_attribute("extractor_type", "chef_propagate")
        else:
            elem.set_string_attribute("extractor_type", "chef_propagate")
    lattice_simulator.update()


def get_fd_quads(lattice):
    '''Return a list of focussing and defocussing quads in a lattice
    
    Arguments
        -lattice - Synergia lattice object
        
    Outputs a pair of lists -  ( [list of focussing quad elements], [list of defocussing quad elements] )
    
    '''
    
    f_quads = []
    d_quads = []
    for elem in lattice.get_elements():
        if elem.get_type() == "quadrupole":
            k1 = elem.get_double_attribute("k1")
            if k1 > 0.0:
                f_quads.append(elem)
            elif k1 < 0.0:
                d_quads.append(elem)
    return (f_quads, d_quads)
    
def get_sextupoles(lattice):
    '''Return a list of positive and negative sextupoles
    
    Arguments
        -lattice - Synergia lattice object
        
    Outputs a pair of lists -  ( [list of positive sextupoles], [list of negative sextupoles] )
    
    '''
    
    p_six = []
    n_six = []
    last = 'n_six'
    
    for elem in lattice.get_elements():
        if elem.get_type() == "sextupole":
            k2 = elem.get_double_attribute("k2")
            if k2 > 0.0:
                p_six.append(elem)
            elif k2 < 0.0:
                n_six.append(elem)
            elif k2 == 0:
                'sextupole strength is 0, so split elements between positive and negative list'
                if last == 'n_six':
                    p_six.append(elem)
                    last = 'p_six'
                elif last == 'p_six':
                    n_six.append(elem)
                    last = 'n_six'
 
    return (p_six, n_six)


#quick helper method
def print_strengths(elemslist, unique=True):
    '''Print the strengths of quads/sextupoles from a list of lattice elements
    
    Arguments:
        elemslist - the list of elements (from lattice.get_elements() usually)
        
    Optional:
        unique - (defaults to True) print only unique elements
    
    '''
    strengths = []
    for elem in elemslist:
        if elem.get_type() == "quadrupole":
            strength = elem.get_double_attribute("k1")
            if unique:
                #only print new strengths
                if strength not in strengths:
                    print elem.get_name() + ' K: ' + str(elem.get_double_attribute("k1"))
                    strengths.append(strength)
                else:
                    pass
            else:
                #print all strengths
                print elem.get_name() + ' K: ' + str(elem.get_double_attribute("k1"))
        #do the same for sextupoles
        if elem.get_type() == "sextupole":
            strength = elem.get_double_attribute("k2")
            if unique:
                #only print new strengths
                if strength not in strengths:
                    print elem.get_name() + ' K2: ' + str(elem.get_double_attribute("k2"))
                    strengths.append(strength)
                else:
                    pass
            else:
                #print all strengths
                print elem.get_name() + ' K2: ' + str(elem.get_double_attribute("k2"))             


def get_maps(lattice):
    '''Return the matrices for each step in the stepper'''
    
    for element in lattice.get_elements():
        element.set_string_attribute("extractor_type", "chef_map")
    map_order = 1
    lattice_simulator = synergia.simulation.Lattice_simulator(lattice,
                                                              map_order)
    steps_per_element = 1
    stepper = synergia.simulation.Independent_stepper_elements(lattice_simulator,
                                                                steps_per_element)
    
    maps = []
    
    stepper.force_update_operations_no_collective()
    for step in stepper.get_steps():
        for operator in step.get_operators():
            if operator.get_type() == 'independent':
                io = synergia.simulation.as_independent_operator(operator)
                slices = io.get_slices()
                #print "%s:" % slices[0].get_lattice_element().get_name()
                if len(slices) > 1:
                    raise RuntimeError, "found multiple slices in an independent operator"
                for operation in io.get_operations():
                    if operation.get_type() == 'fast_mapping':
                        fmo = synergia.simulation.as_fast_mapping_operation(operation)
                        dense_mapping = synergia.simulation.Dense_mapping(fmo.get_fast_mapping())
                        #print "order 0:"
                        #print dense_mapping.get_constant_term()
                        #print "order 1:"
                        #print dense_mapping.get_linear_term()
                        #etype = io.get_slices()[0].get_lattice_element().get_type()
                        #ename = io.get_slices()[0].get_lattice_element().get_name()
                        dmap = dense_mapping.get_linear_term()
                        dmap_copy = np.copy(dmap)
                        #dmmap = dense_mapping.get_linear_term_mad()                        
                        maps.append(np.asarray(dmap_copy))
                        #if etype == 'quadrupole':
                            #dmap = dense_mapping.get_linear_term()
                            #dmap_copy = np.copy(dmap)
                            #dmmap = dense_mapping.get_linear_term_mad()
                            #print ename + ":"
                            #print dmap
                            #if dmap.all() == dmmap.all():
                                #print "mad map is same"
                            #maps.append(np.asarray(dmap_copy))
    return maps
    
def combine_maps(maplist):
    copylist = list(maplist) #copy the list to preserve the original
    dim = copylist[0].shape[0]
    final = np.identity(dim)
    while len(copylist) > 0:
        mult = copylist.pop()
        final = np.dot(mult,final)
    return final