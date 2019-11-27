#!/usr/bin/env python
# Python program to calculate the Alder-Winther cross section for electromagnetic excitations
# based on Winther and Alder, Nucl Phys A 319 (1979) 518

import sys, getopt
from kinematics import *
from coulex import *
import matplotlib.pyplot as plt

def main(argv):
    # define projectile and target, defaults
    pro = nucleus('22O')
    tar = nucleus('197Au')
    #beam energy
    E  = 55.6 #in MeV/u
    #excitation energy in Mev
    Ex = 3.17
    minb = 0
    mat = 4.5345
    targetexc = False


    ## read commandline arguments
    try:
        opts, args = getopt.getopt(argv,"hv:p:t:e:x:b:m:i",["proj=","targ=","ebeam=","exc=","bmin=","matE="])
    except getopt.GetoptError:
        print 'alderwinther.py -p <projectile> -t <target> -e <beam energy (MeV/u)> -x <excitation energy (MeV)> -b <b_min (fm)> -m <matrix element (e^2fm^4)>' 
        sys.exit(2)
#    print opts
#    print args
    for opt, arg in opts:
        if opt == '-h':
            print 'alderwinther.py -p <projectile> -t <target> -e <beam energy (MeV/u)> -x <excitation energy (MeV)> -b <b_min (fm)> -m <matrix element (e^2fm^4)>' 
            sys.exit()
        elif opt == '-v':
            tests(arg)
            sys.exit()
        elif opt == '-i':
            print "target excitation"
            targetexc = True
        elif opt in ("-p", "--proj"):
            pro = nucleus(arg)
        elif opt in ("-t", "--targ"):
            tar = nucleus(arg)
        elif opt in ("-e", "--ebeam"):
            E = float(arg)
        elif opt in ("-b", "--bmin"):
            minb = float(arg)
        elif opt in ("-m", "--matE"):
            mat = float(arg)
        elif opt in ("-x", "--exc"):
            Ex = float(arg)

    ## setup and print kinematics
    kin = kinematics(proj=pro,targ=tar,epu=E,exc=Ex)
    kin.report()

    #distance of closes approach for touchign spheres, grazing
    dgrazing = pro.radius + tar.radius
    bgrazing =  kin.b_fromd(dgrazing)
    print "grazing distance d_grazing = %.4f fm, impact parameter b_grazing = %.4f fm" % (dgrazing,bgrazing)

#    theta = np.linspace(0,5,501)
#    #diff = []
#    for t in theta:
#        tcm = kin.thetacom(t*math.pi/180)
#        b = kin.b_fromthetacom(tcm)
#        diff = relativistic(kin,b, multipole("E",2),mat,targetexc)*10
#        #test = kin.thetacom_fromb(b)
#        #print "%04f\t%04f\t%04f\t%04f" %(t,tcm,b,test)
#        print "%04f\t%04f" %(t,diff)
#    #print theta
#    return
    
    ## if a minimum impact parameter is given
    if minb > 0:
        minbrel = kin.brel(minb)
        print "minimum impact paramter = %.4f fm, relativistic modification = %.4f fm" % (minb,minbrel)
        theta_com = kin.thetacom_fromb(minbrel)
        theta_lab = kin.thetalab(theta_com)
        print "theta_com = %.4f rad, theta_lab = %.4f rad" % (theta_com,theta_lab)
        print "theta_com = %.4f deg, theta_lab = %.4f deg" % (theta_com*180/math.pi,theta_lab*180/math.pi)

        print "sigma = %.5f mb" % (relativistic(kin,minb, multipole("E",2),mat,targetexc)*10)
    ## use the touching sphere distance + 2 fm as safe distance, 
    else:
        dsafe = pro.radius + tar.radius + 2
        bsafe = kin.b_fromd(dsafe)
        minbrel = kin.brel(bsafe)
        print "safe distance d_safe (Rp + Rt + 2 fm) = %.4f, impact parameter b_safe = %.4f fm, relativistic modification = %.4f fm" % (dsafe,bsafe, minbrel)
        theta_com = kin.thetacom_fromb(minbrel)
        theta_lab = kin.thetalab(theta_com)
        print "theta_com = %.4f rad, theta_lab = %.4f rad" % (theta_com,theta_lab)
        print "theta_com = %.4f deg, theta_lab = %.4f deg" % (theta_com*180/math.pi,theta_lab*180/math.pi)
       
        print "sigma = %.5f mb" % (relativistic(kin,bsafe, multipole("E",2),mat,targetexc)*10)
  

    
def tests(arg):
    # test for figure 3 of Winther Alder Nucl Phys A 319 (1979) 518
    if arg == "fig3":
        plt.rcParams['mathtext.default']='regular'
        fig = plt.figure(figsize=(7,5))
        ax = fig.add_subplot(111)
        alpha = 0.681085
        xi = np.linspace(0, 2.5, 1000)
        plt.plot(xi, g(0, xi)/math.pi,label=r'$g_0(\xi)/\pi$')
        plt.plot(xi, g(1, xi)/math.pi/np.log(alpha**2/xi**2+1),label=r'$g_1(\xi)/(\pi ln(\alpha^2 / \xi^2+1))$')
        plt.plot(xi, g(2, xi)*xi**2/4/math.pi,label=r'$g_2(\xi)\xi^2/(4\pi)$')
        plt.plot(xi, g(3, xi)*xi**4/32/math.pi,label=r'$g_3(\xi)\xi^4/(32\pi)$')
        plt.ylim(0, 1.50)
        ax.set_xlabel(r'$\xi$',fontsize=12)
        plt.legend(fontsize=12)
        plt.title(r'$g_\mu(\xi)$ functions for the total cross section, normalized to $g_\mu(\xi=0)=1$',fontsize=12)
        plt.tight_layout()
        plt.show()

    # test G functions (appendix B) of Winther Alder Nucl Phys A 319 (1979) 518
    if arg == "G":
        mult = multipole("E",1)
        mult.report()
        print G(mult,0,0)
        print G(mult,0,0.5)
        print G(mult,0,1)
        print G(mult,0,1.5)
        print G(mult,0,2)


    return
    
if __name__ == "__main__":
   main(sys.argv[1:])

