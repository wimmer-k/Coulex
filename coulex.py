#!/usr/bin/env python
import math
import cmath
import scipy.special as sp
import numpy as np
from kinematics import *

EE = 1.439964535166 #MeV*fm 
HBARC = 197.3269788 #MeV*fm 

## Winther Alder Nucl Phys A 319 (1979) 518
## equation 3.1
def relativistic(kin, bmin, mult, mat, targetexc=False):
    k = kin.exc/HBARC
    ## equation 2.24 (with 1.4 and 2.11)
    xi = k/kin.proj.blab/gamma(kin.proj.blab)*kin.brel(bmin)
    sigma = 0
    #print mult.substates()
    Pm = []
    for ll in mult.substates():
       sigma += abs(G(mult,ll,1./kin.proj.blab))**2 * g(ll,xi)
       Pm.append(abs(G(mult,ll,1./kin.proj.blab))**2 * g(ll,xi))
       #print ll,  abs(G(mult,ll,1./kin.proj.blab))**2, g(ll,xi), sigma
    if targetexc is True:
        sigma *= (kin.proj.Z*EE/HBARC)**2*k**(2*mult.L-2)*mat**2
    else:
        sigma *= (kin.targ.Z*EE/HBARC)**2*k**(2*mult.L-2)*mat**2

    print "P(m) = ", Pm/sum(Pm)
    return sigma
    
    
## Winther Alder Nucl Phys A 319 (1979) 518
## equation 3.4
def g(mu,xi):
    if mu<0:
        return g(-mu,xi)
    Kmu = sp.kn(mu,xi)
    Kmupp = sp.kn(mu+1,xi)
    return math.pi *xi**2 * (Kmupp**2 - Kmu**2 -2*mu/xi*Kmupp*Kmu)

## Winther Alder Nucl Phys A 319 (1979) 518
## Appendix B
def G(mult,mu,x):
    if mu<0:
        return G(mult,-mu,x)
    if mult.EM == "E":
        if mult.L == 1:
            if mu == 0:
                return complex(0, -4./3*cmath.sqrt(math.pi)*cmath.sqrt(x**2-1))
            if mu == 1:
                return complex(1./3*cmath.sqrt(8*math.pi)*x, 0)
        if mult.L == 2:
            if mu == 0:
                return complex(2./5*cmath.sqrt(math.pi)*x*cmath.sqrt(x**2-1), 0)
            if mu == 1:
                return complex(0, 2./5*cmath.sqrt(math.pi/6)*(2*x**2-1))
            if mu == 2:
                return complex(2./5*cmath.sqrt(math.pi/6)*x*cmath.sqrt(x**2-1), 0)
        if mult.L == 3:
            if mu == 0:
                return complex(0, 2./105*cmath.sqrt(math.pi)*(5*x**2-1)*cmath.sqrt(x**2-1))
            if mu == 1:
                return complex(-1./105*cmath.sqrt(math.pi/3)*x*(15*x**2-11), 0)
            if mu == 2:
                return complex(0, -1./21*cmath.sqrt(2./15*math.pi)*(3*x**2-1)*cmath.sqrt(x**2-1))
            if mu == 3:
                return complex(1./21*cmath.sqrt(math.pi/5)*x*cmath.sqrt(x**2-1), 0)
            
            
    elif mult.EM == "M":
        if mult.L == 1:
            if mu == 0:
                return complex(0, 0)
            if mu == 1:
                return complex(0, -1./3*cmath.sqrt(8*math.pi))
        if mult.L == 2:
            if mu == 0:
                return complex(0, 0)
            if mu == 1:
                return complex(2/5*cmath.sqrt(math.pi/6)*x, 0)
            if mu == 2:
                return complex(0, 2/5*cmath.sqrt(math.pi/6)*cmath.sqrt(x**2-1))

    return complex(0,0)


class multipole:
    def __init__(self,*args, **kwargs):
        self.EM = kwargs.get('EM',None)
        self.L = kwargs.get('L',None)
        if self.EM is None:
            self.EM = args[0]
        if self.L is None:
            self.L = args[1]
            
    def report(self):
        print str(self.EM) + str(self.L)


    # magnetic substates for L = -L , -L +1 ,... , L - 1, L
    # arange goes [a,b) -> +1
    def substates(self):
        return np.arange(-self.L,self.L+1)
    
