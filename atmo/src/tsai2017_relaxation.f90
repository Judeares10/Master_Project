MODULE tsai2017_relaxation

IMPLICIT NONE

CONTAINS

FUNCTION eq_co2(nmol_relax,ndim,inloc,imol_relax,neq,nrel)

USE mod_grid, ONLY: ndepth, ttf
USE mod_chem, ONLY: nmol,  imol, molParam

IMPLICIT NONE

INTEGER, INTENT(IN) :: nmol_relax, ndim
INTEGER, INTENT(IN) :: imol_relax(nmol_relax), inloc(ndepth+1,nmol_relax)
REAL, INTENT(INOUT) :: neq(ndepth+1,nmol)
REAL, INTENT(IN)    :: nrel(ndim)

INTEGER :: i, j
REAL    :: dgrt
REAL    :: eq_co2(ndepth+1)
REAL    :: mu(nmol), h(nmol), s(nmol), cp(nmol)

  ! Set CO2 equilibrium profile based on equilibrium with CO and H2O
  DO j = 1, ndepth+1
    DO i = 1, nmol_relax
      IF (imol_relax(i)==imol('CO2')) THEN
        ! Get thermodynamics quantities
        CALL get_thermodynamics(nmol,ttf(j),molParam,mu,h,s,cp)
        ! Gibbs energy of reaction
        dgrt = mu(imol('H2O')) + mu(imol('CO')) - mu(imol('CO2')) - mu(imol('H2'))
        ! CO2 number density
        eq_co2(j) = nrel(inloc(j,2))*nrel(inloc(j,3)) / &
         (neq(j,imol('H2'))*exp(-dgrt))
      END IF
    END DO
  END DO

END FUNCTION eq_co2

! Calculate chemical timescales following method of 
! Tsai et al 2017
SUBROUTINE get_timescale_tsai2017(imol_loc,neq,tscale,tfreeze_relax)

USE mod_grid, ONLY: ndepth, ttf, ppf
USE mod_eos, ONLY: get_nn
USE mod_chem, ONLY: imol, molParam, nmol, molname
USE mod_cst, ONLY: unitL

IMPLICIT NONE

INTEGER :: imol_loc
REAL :: neq(ndepth+1,nmol)
REAL :: tscale(ndepth+1)
REAL :: tfreeze_relax

INTEGER :: idepth
REAL :: tau1, tau2, tloc
REAL :: a, b, c, d
REAL :: mu(nmol), h(nmol), s(nmol), cp(nmol)
REAL :: dGRT17, dGRT18
REAL :: number_density(ndepth+1)

  ! Compute number density
  number_density = get_nn(ndepth+1,ppf,ttf)/unitL**3.

  ! Methane, carbon monoxide and water timescales
  ! Loop over levels
  DO idepth = 1,ndepth+1
   
    if (ttf(idepth) .lt. tfreeze_relax) then 
      tloc = tfreeze_relax
    else 
      tloc = ttf(idepth)
    endif 

    ! Timescales of CO, CH4 and H2O
    IF (any(imol_loc==(/imol('CO'),imol('CH4'),imol('H2O')/))) THEN
 
      d = k4()*neq(idepth,imol('CH3'))*neq(idepth,imol('O-3P'))
      c = min(k2()*neq(idepth,imol('CH2OH'))*neq(idepth,imol('H')),           &
          k3(tloc,number_density(idepth))*neq(idepth,imol('CH2OH')))
      a = k1(tloc,number_density(idepth))*neq(idepth,imol('CH3'))*neq(idepth,imol('OH')) &
          + max(c,d) + k9(tloc)*neq(idepth,imol('CH3OH'))*neq(idepth,imol('H'))
      b = min(k5(tloc)*neq(idepth,imol('C'))*neq(idepth,imol('OH')), &
        k6(tloc)*neq(idepth,imol('CH4'))*neq(idepth,imol('H')))

      ! Timescale for CO-H 
      tau2 = neq(idepth,imol('CO'))/(kh(tloc,number_density(idepth))*neq(idepth,imol('H'))**2)
     
      IF (imol_loc==imol('CH4')) tau1 = neq(idepth,imol('CH4'))/(a+b)
      IF (imol_loc==imol('H2O')) tau1 = neq(idepth,imol('H2O'))/(a+b)
      IF (imol_loc==imol('CO'))  tau1 = neq(idepth,imol('CO'))/(a+b)
      
    ! Timescale of CO
    ELSE IF (imol_loc==imol('CO2')) THEN
      
      tau1 = kco2(tloc)*neq(idepth,imol('CO'))*neq(idepth,imol('OH'))
      tau1 = neq(idepth,imol('CO2'))/tau1
      tau2 = 1./(2.*kh(tloc,number_density(idepth))*neq(idepth,imol('H')))
    
    ! Timescales of NH3 and N2
    ELSE IF (any(imol_loc==(/imol('NH3'),imol('N2')/))) THEN

      ! Thermodynamic properties 
      CALL get_thermodynamics(nmol,tloc,molParam,mu,h,s,cp)
   
      ! Calculate Gibbs energy of reaction for reactions 17 and 18
      dGRT17 = mu(imol('NH2')) + mu(imol('CO')) - mu(imol('OH')) - mu(imol('HCN'))
      dGRT18 = mu(imol('CO'))  + mu(imol('NH')) - mu(imol('HNCO'))
  
      a = k11(tloc)*neq(idepth,imol('NH2'))*neq(idepth,imol('NH2')) +          &
          k12(tloc)*neq(idepth,imol('NH'))*neq(idepth,imol('NH2')) +           &
          min(k13(tloc)*neq(idepth,imol('NO'))*neq(idepth,imol('NH2')),        &
          k14()*neq(idepth,imol('N-4S'))*neq(idepth,imol('NO'))) +             &
          k15(tloc)*neq(idepth,imol('N-4S'))*neq(idepth,imol('NH'))
    
      b = k16(tloc)*neq(idepth,imol('CH2NH2'))*neq(idepth,imol('H')) +         &
          k17(tloc,dGRT17)*neq(idepth,imol('NH2'))*neq(idepth,imol('CO'))   +         &
          k18(tloc,number_density(idepth),dGRT18)*neq(idepth,imol('NH')) &
          *neq(idepth,imol('CO')) + &
          k19(tloc)*neq(idepth,imol('HNCO'))*neq(idepth,imol('H'))   +         &
          k20(tloc)*neq(idepth,imol('N-4S'))*neq(idepth,imol('CO'))
    
      tau1 = 1./(2.*kh(tloc,number_density(idepth))*neq(idepth,imol('H')))
      
      IF (imol_loc==imol('NH3')) THEN
        a = neq(idepth,imol('NH3'))/(2.*a)
        b = neq(idepth,imol('NH3'))/b
        tau2 = (a*neq(idepth,imol('N2')) + b*neq(idepth,imol('HCN'))) / &
          (neq(idepth,imol('N2'))+neq(idepth,imol('HCN')))
      ELSE
        a = neq(idepth,imol('N2'))/(a)
        b = neq(idepth,imol('N2'))/b
        tau2 = (a*neq(idepth,imol('NH3')) + b*neq(idepth,imol('HCN'))) / &
          (neq(idepth,imol('NH3'))+neq(idepth,imol('HCN')))
      END IF
    
    ELSE
      write(*,*) 'Error: get_timescale_tsai2017'
      write(*,*) 'Molecule: ',molname(imol_loc), ' not allowed'
      stop
    END IF
    
    tscale(idepth) = max(tau1,tau2)
  
  END DO
  
END SUBROUTINE get_timescale_tsai2017

! Function to calculate reaction rate of R1
REAL FUNCTION k1(t,n)

IMPLICIT NONE

REAL,INTENT(IN) :: t, n
REAL :: k0, kinf

  k0 = k_bimol(1.93e3,-9.88,-7544.,t)
  k0 = k0 + k_bimol(5.11e-11,-6.25,-1433.,t)
  kinf = k_bimol(1.03e-10,-0.018,16.74,t)
  
  k1 = k_three_body(k0,kinf,n)

END FUNCTION k1

! Function to calculate reaction rate of R2
REAL FUNCTION k2()

IMPLICIT NONE

  k2 = 1.6e-10

END FUNCTION k2

! Function to calculate reaction rate of R3
REAL FUNCTION k3(t,n)

IMPLICIT NONE

REAL,INTENT(IN) :: t, n
REAL :: k0, kinf

  k0 = k_bimol(1.66e-10,0.,-12630.,t)
  kinf = k_bimol(3.0e9,0.,-14600.,t)
  
  k3 = k_three_body(k0,kinf,n)

END FUNCTION k3

! Function to calculate reaction rate of R4
REAL FUNCTION k4()

IMPLICIT NONE

  k4 = 1.4e-10

END FUNCTION k4

! Function to calculate reaction rate of R5
REAL FUNCTION k5(t)

IMPLICIT NONE

REAL,INTENT(IN) :: t

  k5 = k_bimol(1.05e-12,0.5,0.,t)

END FUNCTION k5

! Function to calculate reaction rate of R6
REAL FUNCTION k6(t)

IMPLICIT NONE

REAL,INTENT(IN) :: t

  k6 = k_bimol(2.2e-20,3.,-4040.,t)

END FUNCTION k6

! Function to calculate reaction rate of R7
REAL FUNCTION k7(t)

IMPLICIT NONE

REAL,INTENT(IN) :: t

  k7 = k_bimol(6.82e-20,2.685,-4643.,t)

END FUNCTION k7

! Function to calculate reaction rate of R8
REAL FUNCTION k8(t)

IMPLICIT NONE

REAL,INTENT(IN) :: t

  k8 = k_bimol(6.78e-16,1.5,-854.,t)

END FUNCTION k8

! Function to calculate reaction rate of R9
REAL FUNCTION k9(t)

IMPLICIT NONE

REAL,INTENT(IN) :: t

  k9 = k_bimol(4.91e-19,2.485,-10380.,t)

END FUNCTION k9

! Function to calculate reaction rate of R10
REAL FUNCTION k10(t)

IMPLICIT NONE

REAL,INTENT(IN) :: t

  k10 = k_bimol(8.04e-28,4.,1010.,t)

END FUNCTION k10

! Function to calculate reaction rate of R11
REAL FUNCTION k11(t)

IMPLICIT NONE

REAL,INTENT(IN) :: t

  k11 = k_bimol(2.89e-16,1.02,-5930.,t)

END FUNCTION k11

! Function to calculate reaction rate of R12
REAL FUNCTION k12(t)

IMPLICIT NONE

REAL,INTENT(IN) :: t

  k12 = k_bimol(6.98e-10,-0.27,39.,t)

END FUNCTION k12

! Function to calculate reaction rate of R13
REAL FUNCTION k13(t)

IMPLICIT NONE

REAL,INTENT(IN) :: t

  k13 = k_bimol(7.9e-9,-1.1,-98.,t)

END FUNCTION k13

! Function to calculate reaction rate of R14
REAL FUNCTION k14()

IMPLICIT NONE

  k14 = 3.7e-11

END FUNCTION k14

! Function to calculate reaction rate of R15
REAL FUNCTION k15(t)

IMPLICIT NONE

REAL,INTENT(IN) :: t

  k15 = k_bimol(1.07e-12,0.51,9.3,t)

END FUNCTION k15

! Function to calculate reaction rate of R16
REAL FUNCTION k16(t)

IMPLICIT NONE

REAL,INTENT(IN) :: t

  k16 = k_bimol(6.64e-16,1.5,450.,t)

END FUNCTION k16

! Function to calculate reaction rate of R17 from reverse rate
REAL FUNCTION k17(t,dGRT)

IMPLICIT NONE

REAL,INTENT(IN) :: t, dGRT
REAL :: krev

  ! Reverse reaction rate
  krev = k_bimol(1.07e-13,0.51,-5890.,t)
  
  k17 = krev*exp(-dGRT)
  
END FUNCTION k17

! Function to calculate reaction rate of R18 reverse
REAL FUNCTION k18(t,n,dGRT)

USE mod_chem , ONLY : P0
USE mod_cst  , ONLY : kb, unitL

IMPLICIT NONE

REAL,INTENT(IN) :: t, n, dGRT
REAL :: k0, kinf, krev

  k0 = k_bimol(3.6e4,-3.1,-51300.,t)
  kinf = k_bimol(6e13,0.,-50200.,t)  
  
  krev = k_three_body(k0,kinf,n)
  
  k18 = krev*exp(-dGRT)*P0/(kb*t)/unitL**3.

  
END FUNCTION k18

! Function to calculate reaction rate of R19
REAL FUNCTION k19(t)

IMPLICIT NONE

REAL,INTENT(IN) :: t

  k19 = k_bimol(9.3e-30,4.71,248.,t)

END FUNCTION k19

! Function to calculate reaction rate of R20
REAL FUNCTION k20(t)

IMPLICIT NONE

REAL,INTENT(IN) :: t

  k20 = k_bimol(3.84e-9,0.,-36000.,t)

END FUNCTION k20

! Function to calculate reaction rate of RH
REAL FUNCTION kh(t,n)

IMPLICIT NONE

REAL,INTENT(IN) :: t,n
REAL :: k0, kinf

  k0 = k_bimol(2.7e-31,-0.6,0.,t)
  kinf = k_bimol(3.31e-6,-1.,0.,t)
  
  kh = k_three_body(k0,kinf,n)
  
END FUNCTION kh

! Function to calculate reaction rate of RCO2
REAL FUNCTION kco2(t)

IMPLICIT NONE

REAL,INTENT(IN) :: t

  kco2 = k_bimol(1.05e-17,1.5,259.,t)
  
END FUNCTION kco2

! Function to calculate bimolecular reaction rate
REAL FUNCTION k_bimol(a,b,c,t)

IMPLICIT NONE

REAL, INTENT(IN) :: a, b, c, t

  k_bimol = a*t**(b)*exp(c/t)

END FUNCTION k_bimol

! Function to calculate three body reaction rate
REAL FUNCTION k_three_body(k0,kinf,n)

IMPLICIT NONE

REAL,INTENT(IN) :: k0, kinf, n

  k_three_body = k0*n/kinf
  k_three_body = k0*n/(1.+k_three_body)

END FUNCTION

END MODULE tsai2017_relaxation
