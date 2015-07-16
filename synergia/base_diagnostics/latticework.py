import sys, os
import synergia
from mpi4py import MPI
import synergia_workflow


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


