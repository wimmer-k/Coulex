# Coulex
Relativistic Coulomb excitation

following

A. Winther, K. Alder, Nucl. Phys. A 319 (1979) 518

usage:

chmod 744 alderwinther.py

./alderwinther.py -p <projectile> -t <target> -e <beam energy (MeV/u)> -x <excitation energy (MeV)> -b <b_min (fm)> -m <matrix element (e^2fm^4)>
  
if not b_min is given, a "safe" bmin of Rp + Rt + 2 fm is used


example:

./alderwinther.py -p 22O -e 55.6 -m 4.5345 -x 3.17
