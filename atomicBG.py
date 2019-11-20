#!/usr/bin/env python
import sys, getopt
import math
import numpy as np
from kinematics import *
from comptonprofile import *
import matplotlib.pyplot as plt
from ROOT import TFile, TH2F

#constants
ALPHA = 0.0072973525693 # finestrcuture constant
MEC2 = 510.998950   # mass electron keV
HBARC = 197.3269788 #MeV*fm 

CP = None

def main(argv):
    global CP
    # define projectile and target, defaults
    pro = nucleus('40Ar')
    tar = nucleus('208Pb')
    #tar = nucleus('e')
    #beam energy
    E  = 300 #in MeV/u
    COMP = 'REC,PB,SEB'
    
    ## read commandline arguments
    try:
        opts, args = getopt.getopt(argv,"hv:p:t:e:c:",["proj=","targ=","ebeam=","comp="])
    except getopt.GetoptError:
        print 'atomicBG.py -p <projectile> -t <target> -e <beam energy (MeV/u)> -c <component: REC, PB, SEB, all>' 
        sys.exit(2)
#    print opts
#    print args
    for opt, arg in opts:
#        print opt, arg
        if opt == '-h':
            print 'atomicBG.py -p <projectile> -t <target> -e <beam energy (MeV/u)> -c <component: REC, PB, SEB, all>' 
            sys.exit()
        elif opt == '-v':
            tests(arg)
            sys.exit()
        elif opt in ("-p", "--proj"):
            pro = nucleus(arg)
        elif opt in ("-t", "--targ"):
            tar = nucleus(arg)
        elif opt in ("-e", "--ebeam"):
           E = float(arg)
        elif opt in ("-c", "--comp"):
             COMP = arg
    if COMP == 'all':
        COMP = 'REC,PB,SEB'

    ## setup and print kinematics
    kin = kinematics(proj=pro,targ=tar,epu=E,exc=0)
    kin.report()

    CP = comptonprofile(tar=tar.symbol)

    #energy range and stepsize
    er = [10,1010,20]
    #theta range and stepsize
    tr = [0,180,10]
    
    hfile = TFile('test.root', 'RECREATE')
    if 'REC' in COMP:
        hREC  = TH2F('hREC', 'REC', (tr[1]-tr[0])/tr[2],tr[0],tr[1], (er[1]-er[0])/er[2],er[0],er[1])
        sREC = 0
        for t in range(180):
            e,s = dSdO_REC(kin,t*math.pi/180 ,1)  # b/sr
            sREC = sREC+2*math.pi*s*math.sin(t*math.pi/180)*math.pi/180
            hREC.Fill(t,e,s)
            e,s = dSdO_REC(kin,t*math.pi/180 ,2)  # b/sr
            hREC.Fill(t,e,s)
            sREC = sREC+2*math.pi*s*math.sin(t*math.pi/180)*math.pi/180
        print "cross section REC = %.3f barn" % sREC

    if 'PB' in COMP:
        hPB  = TH2F('hPB', 'PB', (tr[1]-tr[0])/tr[2],tr[0],tr[1], (er[1]-er[0])/er[2],er[0],er[1])
        sPB = 0
        for ee in np.arange(*er):
            we = 1
            if ee==er[0] or ee == er[1]:
                we=0.5
            for t in np.arange(*tr):
                    #print ee
                wt =1
                if t==tr[0] or t ==tr[1]:
                    wt=0.5
                e,s = dSdEdO_PB(kin,ee,t*math.pi/180)  # b/sr
                sPB = sPB+2.*math.pi*s*math.sin(t*math.pi/180)*math.pi/180*tr[2]*er[2]*wt*we
                hPB.Fill(t,e,s)
        print "cross section PB = %.3f barn" % sPB

    if 'SEB' in COMP:
        hSEB  = TH2F('hSEB', 'SEB', (tr[1]-tr[0])/tr[2],tr[0],tr[1], (er[1]-er[0])/er[2],er[0],er[1])
        sSEB = 0
        print "cross section SEB = %.3f barn" % sSEB
    hfile.Write()

    return


def dSdO_REC(kin,Theta, shell):
    E_shell = (1.0-math.sqrt(1.0-math.pow(ALPHA*kin.proj.Z/shell,2)))*MEC2 # binding energy of shell K, L, ... for hydrogen-like
    E_REC = (gamma(kin.proj.blab)-1)*MEC2 + E_shell
    Theta_cm = math.acos( (kin.proj.blab - math.cos(Theta)) / (1-kin.proj.blab*math.cos(Theta)) )
    boost = 1./(gamma(kin.proj.blab)*(1-kin.proj.blab*math.cos(Theta)))    #derivation in log
    dSdO = boost**2 * BetheSalpeter(kin.proj.Z,kin.proj.blab,Theta_cm)/math.pow(shell,3)
    # correction, source unknown, seems to be small
    if shell == 1:
        dSdO = dSdO * min(max(0,1-(kin.proj.Z-Zeff(kin))/2),1)
    if shell == 2:
        dSdO = dSdO * min(max(0,1-(kin.proj.Z-Zeff(kin)-2)/8),1)
    #print kin.proj.Z , E_shell, E_REC, Theta*180/math.pi , E_REC*boost, BetheSalpeter(kin.proj.Z,kin.proj.blab,Theta_cm), dSdO*kin.targ.Z
    return E_REC*boost, dSdO*kin.targ.Z

def BetheSalpeter(Z, beta, theta_cm):
    # H. E. Bethe and E. E. Salpeter, Handb. Phys. +3, 408 (1957).
    # Bethe and Salpeter, Quantum mechanics of one and two electron atoms, Springer 1957
    # Phys. Rev. A 26 (1982) 154
    # theta is angle between trajectory of electron and emitted photon
    eta = Z*ALPHA/beta
    eta2 = eta*eta
    F = math.pow(eta*eta2/(1+eta2),2)*math.exp(-4*eta*math.atan(1/eta))/(1-math.exp(-2*math.pi*eta))
    # in fm*2
    F =  F* pow(2,5)*math.pi*EE*HBARC/math.pow(MEC2/1e3,2)
    # in b/sr
    return F*math.pow(math.sin(theta_cm),2)/math.pow(1-beta*math.cos(theta_cm),4) /100

def Zeff(kin):
    b = 0.886*math.sqrt(kin.epa*40.0)/pow(kin.proj.Z,2./3)
    a = b + 0.0378*math.sin(1.5708*b)
    return kin.proj.Z*(1.0-math.exp(-a)*(1.034-0.1777*math.exp(-0.08114*kin.proj.Z)))
    
def dSdEdO_PB(kin, e_phot, theta):    
    theta_cm = math.acos( (kin.proj.blab - math.cos(theta)) / (1-kin.proj.blab*math.cos(theta)) )
    boost = 1./(gamma(kin.proj.blab)*(1-kin.proj.blab*math.cos(theta)))    #derivation in log
    e_phot_cm = e_phot/boost
    # integrate over compton profile J(pz)
    # R. Anholt et al., Phys. Rev. A 33 (1986) 2270, eq 15
    integral = GaussLegende(dSdEdOdp_PB,[-100,100],args=(kin,e_phot_cm,theta_cm),accuracy =0.1)
    return e_phot, boost*integral*kin.targ.Z
    
def dSdEdOdp_PB(pz, kin, e_phot, theta_cm):
    # R. Anholt et al., Phys. Rev. A 33 (1986) 2270, eq 16
    # q in me^2/hbar in F. Biggs, L.B. Mendelsohn, J.B. Mann, Atomic Data and Nuclear Data Tables, 16, 201 (1975)
    e_elec = (math.sqrt((gamma(kin.proj.blab)*kin.proj.blab+pz*EE/HBARC)**2+1)-1)*MEC2
    return CP.profile(pz)*BetheHeitler(e_phot,e_elec,theta_cm,kin.proj.Z)
    
def BetheHeitler(e_phot, e_elec, theta, Z_pro):
    if e_phot < 0:
        return 0
    if e_phot > 0.99*e_elec:
        return 0
    
    # H.W. Koch and J.W. Motz, formula 2BN page 924 in Rev. Mod. Phys. 31 (1959) 920
    # also F. Sauter, Ann. Phys. 20 (1934) 404
    k = e_phot/MEC2
    #initial kinetic, total energy, momentum, and velocity of electron
    T0 = e_elec/MEC2
    E0 = T0 + 1.0
    p0 = math.sqrt(T0*(T0+2.0))
    beta0 = p0/E0
    #final kinetic, total energy, momentum, and velocity of electron
    T = T0 - k
    E = E0 - k
    p = math.sqrt(T*(T+2.0))
    beta = p/E
    Q2 = p0*p0 + k*k-2.0*p0*k*math.cos(theta)
    Q = math.sqrt(Q2)
    
    L = math.log( (E*E0-1.0+p*p0)/(E*E0-1.0-p*p0) )
    D0 = E0 - p0*math.cos(theta)
    eps = math.log( (E+p)/(E-p) )
    epsQ = math.log( (Q+p)/(Q-p) )

    k2 = k*k
    E02 = E0*E0
    E2 = E*E
    E0E = E0*E
    p02 = p0*p0
    p0p = p0*p
    D02 = D0*D0

    #d^2 sigma /dE_phot /dOmega
    # in barn/keV/sr
    
    a = 8*math.sin(theta)**2*(2*E02+1) / (p02*math.pow(D0,4))
    b = - 2.*(5.*E02+2*E0E+3.)/p02/D02
    c = - 2.*(p02-k2)/Q2/D02 + 4.*E/p02/D0
    d = 4.*E0*math.sin(theta)**2*(3.*k-p02*E) / (p02*math.pow(D0,4))
    e = 4.*E02*(E02+E2) / (p02*D02)
    f = (2.-2.*(7.*E02-3.*E0E+E2)) / (p02*D02)
    g = 2.*k*(E02+E0E-1)/(p02*D0)
    h = -4.*eps/p/D0 + epsQ/p/Q * (4./D02 - 6.*k/D0 - 2.*k*(p02-k2)/Q2/D0)
    return Z_pro*Z_pro*CoulombCorrection(beta0,beta,Z_pro) * Lambda(T0,Z_pro) * \
        p/p0 * math.pow(EE/(MEC2/1000),2)/8/math.pi*ALPHA/100/e_phot * \
        ( 8*math.sin(theta)**2*(2*E02+1) / (p02*math.pow(D0,4)) \
          - 2*(5*E02+2*E0E+3)/p02/D02 \
          - 2*(p02-k2)/Q2/D02 + 4*E/p02/D0 \
          + L/p0p * (4*E0*math.sin(theta)**2*(3*k-p02*E) / (p02*math.pow(D0,4)) \
                     + 4*E02*(E02+E2) / (p02*D02) \
                     + (2-2*(7*E02-3*E0E+E2)) / (p02*D02) \
                     + 2*k*(E02+E0E-1)/(p02*D0) ) \
          - 4*eps/p/D0 + epsQ/p/Q * (4/D02 - 6*k/D0 - 2*k*(p02-k2)/Q2/D0) )    

def CoulombCorrection(beta0, beta, Z):
    # G. Elwert, Ann Phys. 34 (1939) 178
    # equaton 52
    #print Z*ALPHA * (1/beta - 1/beta0)
    return beta0*(1-math.exp(-2*math.pi*Z*ALPHA/beta0)) / beta/(1-math.exp(-2*math.pi*Z*ALPHA/beta))

def Lambda(T0,Z):
    #correction for Born approximation, origin unknown
    t = np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0, 2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,20.0, 100.0,200.0])
    al = [1.03,1.11,1.17,1.25,1.22,1.28,1.30,1.32,1.33,1.34, 1.27,1.18,1.13,1.08,1.06,1.04,1.03,1.02,1.01,1.00, 1.00,1.00]
    au = [1.20,1.31,1.38,1.42,1.46,1.50,1.52,1.54,1.55,1.54, 1.39,1.17,1.12,1.10,1.08,1.06,1.04,1.03,1.01,0.96, 0.92,0.92]
    if T0<t[0]:
        return 0
    i = np.searchsorted(t, T0)
    #print T0 , i
    #print t[i], al[i], au[i]
    dt = T0 - t[i-1]
    fal = al[i-1] + dt*(al[i]-al[i-1])/(t[i]-t[i-1])
    fau = au[i-1] + dt*(au[i]-au[i-1])/(t[i]-t[i-1])
    return max(fal + (Z-13)*(fau-fal)/66 ,0)

def ReadComptonProfile(tar, draw=False):
    vals  = np.loadtxt("comptonprofiles.dat",unpack=True, comments='#',skiprows=1)
    idtag = {'1H': 1, '9Be': 2, '12C': 3, '179Au': 4, '208Pb':5}
    if tar.symbol not in idtag:
        print "Compton profile of %s no implemented" % tar.symbol
        print "check: F. Biggs, L.B. Mendelsohn, J.B. Mann, Atomic Data and Nuclear Data Tables, 16, 201 (1975)"
        return
    J = interp1d(vals[0],vals[idtag[tar.symbol]], kind='linear')
    if draw:
        q = np.linspace(vals[0].min(),vals[0].max(),100)
        plt.plot(q,J(q))
        plt.plot(vals[0],vals[idtag[tar.symbol]])
        #plt.yscale('log')
        plt.show()
    return J

def NormComptonProfile(tar, norm):
    vals  = np.loadtxt("comptonprofiles.dat",unpack=True, comments='#',skiprows=1)
    idtag = {'1H': 1, '9Be': 2, '12C': 3, '179Au': 4, '208Pb':5}
    if tar.symbol not in idtag:
        print "Compton profile of %s no implemented" % tar.symbol
        print "check: F. Biggs, L.B. Mendelsohn, J.B. Mann, Atomic Data and Nuclear Data Tables, 16, 201 (1975)"
        return
    return interp1d(vals[0],vals[idtag[tar.symbol]], kind='linear')


points = [-0.9894009349916499, -0.9445750230732326, -0.8656312023878318, -0.7554044083550030, -0.6178762444026438, -0.4580167776572274, -0.2816035507792589, -0.0950125098376374, 0.0950125098376374, 0.2816035507792589, 0.4580167776572274, 0.6178762444026438, 0.7554044083550030, 0.8656312023878318, 0.9445750230732326, 0.9894009349916499]
weights = [0.0271524594117541, 0.0622535239386479, 0.0951585116824928, 0.1246289712555339, 0.1495959888165767, 0.1691565193950025, 0.1826034150449236, 0.1894506104550685, 0.1894506104550685, 0.1826034150449236, 0.1691565193950025, 0.1495959888165767, 0.1246289712555339, 0.0951585116824928, 0.0622535239386479, 0.0271524594117541]

def GaussLegende(func,interval,args,accuracy = 0.1):
    #print interval, args,accuracy

    inv, a,b = interval[1] < interval[0], min(interval[0],interval[1]), max(interval[0],interval[1])
    #print inv, a, b
    nsteps = int(1.0/accuracy)+1
    dx = float(b-a)/nsteps
    integral = 0
    for x in np.linspace(a,b,nsteps, endpoint=False):
        for p,w in zip(points,weights):
            integral = integral + dx/2 *w*func(dx/2*(p+1)+x,*args)
            #print integral
    return -integral if inv else integral
if __name__ == "__main__":
   main(sys.argv[1:])

