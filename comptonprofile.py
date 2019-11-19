#!/usr/bin/env python
import math
import re
import numpy as np
from scipy.interpolate import interp1d
import scipy.integrate as integrate
import matplotlib.pyplot as plt
PROFILEFILE = '/home/wimmer/programs/coulex/comptonprofiles.dat'
idtag = {'1H': 1, '9Be': 2, '12C': 3, '179Au': 4, '208Pb':5}

class comptonprofile:
    def __init__(self,*args, **kwargs):
        self.tar = kwargs.get('tar',None)
        self.init()
        
    def init(self):
        vals  = np.loadtxt(PROFILEFILE,unpack=True, comments='#',skiprows=1)
        if self.tar not in idtag:
            print "Compton profile of %s no implemented" % self.tar
            print "check: F. Biggs, L.B. Mendelsohn, J.B. Mann, Atomic Data and Nuclear Data Tables, 16, 201 (1975)"
            return
        self.J = interp1d(vals[0],vals[idtag[self.tar]], kind='linear')
        normJ = integrate.quad(self.profile,-100,100,limit=200,full_output=1)
        self.J = interp1d(vals[0],vals[idtag[self.tar]]/normJ[0], kind='linear')
        
    def draw(self):
        q = np.linspace(-100,100,100)
        plt.plot(q,self.profile(abs(q)))
        plt.show()
 
    def profile(self,q):
        return self.J(abs(q))
