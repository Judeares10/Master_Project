MODULE cs2006_relaxation

IMPLICIT NONE

CONTAINS

! Calculate chemical timescales following method of 
! Tsai et al 2017
SUBROUTINE get_timescale_cs2006(neq,tscale)

USE mod_grid, ONLY: ndepth, ttf, ppf
USE mod_eos, ONLY: get_nn
USE mod_chem, ONLY: imol, molParam, nmol, Amolf
USE mod_cst, ONLY: unitL, unitP

IMPLICIT NONE

INTEGER :: nmol_relax
REAL :: neq(ndepth+1,nmol)
REAL :: tscale(ndepth+1)

INTEGER :: idepth
REAL ::  tloc, ploc, nloc, k0, kinf, k, Keq
REAL :: a, b, c, f, co
REAL :: mu(nmol), h(nmol), s(nmol), cp(nmol)
REAL :: number_density(ndepth+1)


REAL,PARAMETER :: P0 = 1.013250e+5
REAL,PARAMETER :: acarbon = 4.570e-4
REAL,PARAMETER :: aoxygen = 8.3104e-4

  ! Compute number density
  number_density = get_nn(ndepth+1,ppf,ttf)/unitL**3.
  
  DO idepth = 1, ndepth+1
  
    tloc = ttf(idepth)
    ploc = ppf(idepth)*unitP*0.1 ! in SI
    nloc = number_density(idepth)
  
    CALL get_thermodynamics(nmol,tloc,molParam,mu,h,s,cp)
    
    ! Calculate equilibrium constant
!		Keq = mu(imol('CO')) + 3.*mu(imol('H2')) - mu(imol('H2O')) - mu(imol('CH4'))
!		Keq = exp(-Keq)

!		! Calculate f
!		f = Keq/(ploc/P0)**2.
!		f = f/Amolf(idepth,imol('H2'))**3.

!		! Calculate a, b and c
!		a = f
!		b = -f*(acarbon+aoxygen)-1.
!		c = f*acarbon*aoxygen

!		! Calculate CO abundance
!		co = -b-(b**2.-(4*a*c))**0.5
!		co = co/(2.*a)

!    Amolf(idepth,imol('CO')) = co*nloc
    
		! Step 1: Calculate the rate constant of limiting reaction
		!         number 9 from C&S2006

		! Calculate k in low pressure limit, EQ13 from C&S2006
		k0 = 1.4e-6*tloc**(-1.2)
		k0 = k0*exp(-7800/tloc)

		! Calculate k in high pressure limit, EQ14 from C&S2006
		kinf = 1.5e11*tloc
		kinf = kinf*exp(-12880/tloc)

		! Calculate overall k, EQ15 from C&S2006
		k = k0*kinf/(k0*nloc+kinf)

		! Step 2: Calculate the equilibrium constant for reaction 17 of C&S 2006:
		!         CH3O <-> CO + 3/2H2
		Keq = mu(imol('CO')) + (3./2.)*mu(imol('H2')) - mu(imol('CH3O'))
		Keq = exp(-Keq)

		! Step 3: Calculate the chemical timescale
		tscale(idepth) = nloc*k*Amolf(idepth,imol('H2'))**(3./2.)*(ploc/P0)**(3./2.)
		tscale(idepth) = Keq/tscale(idepth)

  END DO

END SUBROUTINE get_timescale_cs2006

END MODULE cs2006_relaxation
