#!/usr/bin/env python
import math
import re
#import numpy as np
MASSFILE = '/home/wimmer/programs/coulex/ame2016.dat'
AMU = 931.4940954 # MeV
R0 = 1.25
EE = 1.439964535166

class kinematics:
    def __init__(self,*args, **kwargs):
        
        self.proj = kwargs.get('proj',None)
        self.targ = kwargs.get('targ',None)
        self.ejec = kwargs.get('ejec',None)
        self.reco = kwargs.get('reco',None)
        self.epa = kwargs.get('epa',None)
        self.epu = kwargs.get('epu',None)
        self.exc = kwargs.get('exc',None)
        if self.epa is None:
            self.epa = self.epu/self.proj.A*self.proj.mass/AMU
        if self.epu is None:
            self.epu = self.epa*self.proj.A/(self.proj.mass/AMU)

        if self.ejec is None and self.reco is None:
            print "(in)elastic scattering"
            self.ejec = self.proj
            self.reco = self.targ

        self.qval = (self.proj.mass+self.targ.mass) - (self.ejec.mass+self.reco.mass) - self.exc
            
        self.init()

    def init(self):
        self.redmass = (self.proj.mass * self.targ.mass)/(self.proj.mass + self.targ.mass)
        self.comenergy = math.sqrt(self.proj.mass**2 + self.targ.mass**2 + 2*self.targ.mass*(self.proj.mass+self.epa*self.proj.A) )

        self.proj.COM(Ecom = self.comenergy/2 + (self.proj.mass**2 - self.targ.mass**2)/(2*self.comenergy))
        self.targ.COM(Ecom = self.comenergy/2 - (self.proj.mass**2 - self.targ.mass**2)/(2*self.comenergy))
        self.proj.LAB(Tlab = self.epa*self.proj.A)
        self.targ.LAB(Tlab = 0)
            
        self.targ.bcom = - self.targ.bcom

        #velocity of the center of mass system
        self.betacm = (self.proj.Plab-self.targ.Plab)/(self.proj.Elab+self.targ.Elab)
        self.gammacm = gamma(self.betacm)
       
        self.Tcom_initial = self.comenergy - self.proj.mass - self.targ.mass

        #final state
        self.Tcom_final   = self.Tcom_initial + self.qval

        self.reco.COM(Tcom = self.Tcom_final/2*(self.Tcom_final+2*self.ejec.mass)/self.comenergy)
        self.ejec.COM(Tcom = self.Tcom_final/2*(self.Tcom_final+2*self.reco.mass)/self.comenergy)

        self.reco.bcom =  self.reco.bcom
       
        #print self.betacm
        
    def report(self):
        print "projectile: ",
        self.proj.report()
        print "target:     ",
        self.targ.report()
        print "Ebeam = %.2f AMeV = %.2f MeV/u = %.2f MeV" % (self.epa, self.epu, self.epa*self.proj.A)
        print "E_x = %.3f MeV" % self.exc
        print "beta_p(lab) = %.5f" % self.proj.blab
        
    def b_fromd(self,d):
        a = self.proj.Z*self.targ.Z*EE/self.redmass/self.proj.blab**2/gamma(self.proj.blab)
        return math.sqrt(d*d-2*a*d/gamma(self.proj.blab))

    def brel(self,b):
        a = self.proj.Z*self.targ.Z*EE/self.redmass/self.proj.blab**2/gamma(self.proj.blab)
        return b + math.pi/2*a

    def thetacom_fromb(self,b):
        a = self.proj.Z*self.targ.Z*EE/self.redmass/self.proj.blab**2/gamma(self.proj.blab)
        return 2*math.atan2(a,b)

    def thetacom(self, anglelab, part=None):
        if part is None:
            part = self.ejec
        gtan = math.tan(anglelab)**2*self.gammacm**2
        x = self.betacm/part.bcom
        if(math.tan(anglelab)>0):
            return math.acos( (-x*gtan+math.sqrt( 1+gtan*(1-x*x) ))/(1+gtan) )
        else:
            return math.acos( (-x*gtan-math.sqrt( 1+gtan*(1-x*x) ))/(1+gtan) )

    def thetalab(self, anglecom, bcm=None):
        if bcm is None:
            bcm = self.ejec.bcom
        x = self.betacm/bcm
        return math.atan2(math.sin(anglecom),self.gammacm*(math.cos(anglecom)+x));

#    def sigmacm(self, angle_lab, sigma_lab):
#        angle_com = self.thetacom(angle_lab,self.ejec.bcom)
#        ggg = self.proj.mass*self.reco.mass/self.targ.mass/self.ejec.mass * self.Tcom_initial/self.Tcom_final
#        ggg = math.sqrt(ggg)
#        wurzel=1.+ggg**2+2.*ggg*math.cos(math.pi-angle_com);
#        wurzel = math.sqrt(wurzel);
#        return sigma_lab/(1/self.gammacm*wurzel*wurzel*wurzel/(1+ggg*math.cos(math.pi-angle_com)))  

    def cm_fromlab(self, angle_lab, sigma_lab, part=None):
        if part is None:
            part = self.ejec            
        angle_com = self.thetacom(angle_lab,part)
        x = self.betacm/part.bcom
        jacobi = self.gammacm *(1+x*math.cos(angle_com))/ (math.sin(angle_com)**2+self.gammacm**2*(math.cos(angle_com)+x)**2)**(3./2)
        #print jacobi
        return angle_com, sigma_lab*jacobi
        
class nucleus:
    def __init__(self,*args, **kwargs):
        if len(args) == 1:
            self.symbol = args[0]
            self.A = None
            self.Z = None
        else:
            self.symbol = kwargs.get('iso',None)
            self.A = kwargs.get('A',None)
            self.Z = kwargs.get('Z',None)
        if self.symbol is None:
            self.symbol = str(self.A)+symbols[self.Z]
        if self.A is None or self.Z is None:
             self.A, self.Z = self.AZ_fromsymbol()

        self.mass = self.mass_fromsymbol()
        self.radius = R0*math.pow(self.A,1./3)
        #self.report()

    def report(self):     
        print "%s: A = %d, Z= %d, N= %d, M = %.4f u" % (self.symbol,self.A,self.Z,self.A-self.Z,self.mass/AMU)

    def LAB(self,**kwargs):
        self.Tlab =  kwargs.get('Tlab',None)
        self.Elab = self.Tlab + self.mass
        self.Plab = math.sqrt(self.Tlab**2+2*self.Tlab*self.mass)
        self.blab = self.Plab/self.Elab
        
    def COM(self,**kwargs):      
        self.Tcom = kwargs.get('Tcom',None)
        if self.Tcom is None:
            self.Ecom = kwargs.get('Ecom',None)
            self.Tcom = self.Ecom - self.mass
        else:
            self.Ecom = self.Tcom + self.mass
        self.Pcom = math.sqrt(self.Ecom**2-self.mass**2)
        self.bcom = self.Pcom/self.Ecom
 
    def AZ_fromsymbol(self):
        digit_pattern = re.compile(r'\D')
        alpha_pattern = re.compile(r'\d')
        A = filter(None, digit_pattern.split(self.symbol))
        sym = filter(None, alpha_pattern.split(self.symbol))
        return int(A[0]), symbols.index(sym[0])

    def mass_fromsymbol(self):
        with open(MASSFILE,'r') as f:
            for line in f:
                vals = line.split()
                if len(vals) == 6 and int(vals[2]) == self.Z and int(vals[3]) == self.A:
                    return float(vals[5])/1e6*AMU
        print "Warning! Mass for %s not found!" % self.symbol
        return 0

    
symbols = ["n", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"]


def gamma(beta):
    return 1./math.sqrt(1-beta*beta)

def beta(gamma):
    return math.sqrt(1-1/gamma/gamma)

def gammaE(ekin, mass):
    return (ekin+mass)/mass
