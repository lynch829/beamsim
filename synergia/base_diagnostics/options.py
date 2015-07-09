import sys

#Options class definition - All other diagnostics will make use of this class or subclasses

#To Do: Make some subclasses with more specific options for different plots

class Options:
    def __init__(self):
        self.hist = False
        self.plots = 1
        self.inputfile = None
        self.outputfile = None
        self.show = True
        self.hcoord = None
        self.vcoord = None
        self.bins = 20
        self.minh = -sys.float_info.max
        self.maxh = sys.float_info.max
        self.minv = -sys.float_info.max
        self.maxv = sys.float_info.max
        self.contour = None
        self.num_contour = None
        self.lattice_name = None
        self.save = False #default don't save
        self.ID = None
        self.path = None
        self.norm = False
        self.variance = None
        self.relpath = None
        self.lattice_simulator = None
        
    def set_new(self,name,val):
        '''A class method which wraps setattr()'''
        return setattr(self,name,val)