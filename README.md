# alderwinther
Relativistic Coulomb excitation

following

A. Winther, K. Alder, Nucl. Phys. A 319 (1979) 518

usage:

chmod 744 alderwinther.py

./alderwinther.py -p <projectile> -t <target> -e <beam energy (MeV/u)> -x <excitation energy (MeV)> -b <b_min (fm)> -m <matrix element (e^2fm^4)>
  
if no b_min is given, a "safe" bmin of Rp + Rt + 2 fm is used


example:

./alderwinther.py -p 22O -e 55.6 -m 4.5345 -x 3.17


# atomicBG
Calculation of atomic background spectra with relativistic heavy ion beams

following

R. Anholt et al., Phys. Rev. A 33 (1986) 2270

based on the abkg code R. Holzmann (GSI) 1998

usage:

chmod 744 atomicBG.py

./atomicBG.py -p <projectile> -t <target> -e <beam energy (MeV/u)> -c <component: REC, PB, SEB, all>

example:

./atomicBG.py -p 110Zr -t 9Be - 205 -c all

The cross section depends on the range of the integration. The energy range is defined in lines 68 of atomicBG.py and the following and should be adjusted.
