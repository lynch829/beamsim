import sys
import numpy as np
import matplotlib.pyplot as plt
import synergia


def get_lf_fns(lattice, lattice_simulator):
    '''
    Return the lattice functions for every element in the lattice, assuming periodicity
    '''
    #define the lattice function names
    lf_names = ("beta_x", "alpha_x", "beta_y", "alpha_y", 
            "psi_x", "psi_y","D_x", "Dprime_x", "D_y", "Dprime_y")
    #construct an empty array of dictionaries corresponding to each lattice element
    lf_info = [{}]
    #loop through lattice elements
    index=0
    for element in lattice.get_elements():
        index += 1
        #get lattice functions for a given element
        lattice_functions = lattice_simulator.get_lattice_functions(element)
        #define dictionary for this element
        lf = {}
        lf['name'] = element.get_name()
        lf['s'] = lattice_functions.arc_length
        lf['length'] = element.get_length()
        #loop through lattice functions for the element
        for lfname in lf_names:
            lf[lfname] = getattr(lattice_functions,lfname)
        #append individual element to array
        lf_info.append(lf)
        
    #construct initial element lf_info[0]
    lf_info[0]['s'] = 0.0
    lf_info[0]['name'] = lattice.get_name()
    lf_info[0]['length'] = 0.0
    lf_info[0]['psi_x'] = 0.0
    lf_info[0]['psi_y'] = 0.0
    #Take lattice functions from the final element to ensure periodicity
    for lfname in lf_names:
        lf_info[0][lfname] = lf_info[-1][lfname]
        
    return lf_info

#A basic plotting function for plotting said lattice functions
def lf_plot(lf_info, lf_fns,save=False):
    '''
        Plot lattice functions on same plot
        
        Arguments:
        lf_info - array of lattice element dictionaries, each containing lattice functions
        lf_fns - list of lattice function names (e.g. 'betax', etc.)

    '''
    
    #get s coordinate and put in numpy array
    ss = np.array([lfd['s'] for lfd in lf_info])
    
    #grab lattice functions as needed
    for fn in lf_fns:
        #create y array for each lattice function
        y = np.array([lfd[fn] for lfd in lf_info])
        #add to plot
        plt.plot(ss, y, label=fn)
    
    #add plot labels
    plt.xlabel('s', fontsize=12)
    plt.ylabel(', '.join(lf_fns), fontsize=12)
    
    #add legend and show
    plt.legend(loc='best')
    title = 'Lattice functions for Linearized Iota3-4'
    plt.title(title, y=1.05, fontsize=15)

    fig1 = plt.gcf()
    plt.show()
    
    if save:
        name = '_'.join(lf_fns) +lf_info[0]['name'] + '.pdf'
        fig1.savefig(name)
    
if __name__ == '__main__':
    
    #handle arguments
    lattice = sys.argv[0]
    #define lattice functions being plotted - s is always independent coordinate
    lf_fns = sys.argv[1:]
    
    #get lattice function array
    lf_array = get_lf_fns(lattice)
    lf_plot(lf_array, lf_fns)