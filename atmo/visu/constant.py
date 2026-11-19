
import math

print("define physical constants in cgs")

gramme      = 1.0
centimetre  = 1.0
seconde     = 1.0
degre       = 1.0
gauss       = 1.0

# additional units
dyne    =  1.000000e+00 * gramme*centimetre/seconde**2
erg     =  1.000000e+00 * gramme*centimetre**2/seconde**2

# physical constants
c       =  2.997925e+10 * centimetre/seconde
c2      =  c*c
eV      =  1.602200e-12 * erg
G       =  6.673000e-08 * dyne*centimetre**2/gramme**2
hplanck =  6.626200e-27 * erg*seconde
kb      =  1.380620e-16 * erg/degre
uma     =  1.660531e-24 * gramme
a_R     =  7.564640e-15 * erg/centimetre**3/degre**4
Lsol    =  3.826800e+33 * erg/seconde
Msol    =  1.989100e+33 * gramme
Rsol    =  6.959900e+10 * centimetre
pc      =  3.085678e+18 * centimetre
mu0     =  4.000000e+00 * math.pi * gauss**2 * centimetre * seconde**2 / gramme
au      =  1.495980e+13 * centimetre
sig_sb  =  5.670373e-05 * erg/centimetre**2/seconde/degre**4

# useful units
cm2     =  1.000000e+00 * centimetre**2
cm3     =  1.000000e+00 * centimetre**3
ms      =  1.000000e+02 * centimetre/seconde
kms     =  1.000000e+05 * centimetre/seconde
an      =  3.155760e+07 * seconde
Mans    =  3.155760e+13 * seconde
Gans    =  3.155760e+16 * seconde
kpc     =  3.085678e+21 * centimetre
Mpc     =  3.085678e+24 * centimetre
Kelvin  =  1.000000e+00 * degre
kg      =  1.000000e+03 * gramme
metre   =  1.000000e+02 * centimetre
Pascal  =  1.000000e+00 * kg/metre   /seconde**2
joule   =  1.000000e+00 * kg*metre**2/seconde**2
H0      =  5.000000e+01 * kms/Mpc
