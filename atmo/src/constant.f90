!> \file
!! Contains modules mod_cst, subroutine mod_cst::init_cst() \n
!! initialize the physical constants and units
!<

!> Hold the code units, the physical units and the physical constants 
module mod_cst

  real, parameter :: Pi = 3.1415926535898  !< \f$  \pi \f$
  real :: kb                               !< Boltzmann constant \f$ k_{B} \f$
  real :: amu                              !< Atomic mass unit \f$ u \f$
  real :: cvel                             !< Light speed
  real :: hplanck                          !< Planck constant
  real :: sigsb                            !< Stefan-Boltzmann constant
  real :: atm                              !< Atmosphere unit for pressure
  real :: Rjup                             !< Radius of Jupiter
  real :: Rsun                             !< Radius of the sun
  real :: au                               !< astronomical unit
  real :: day                              !< day 
  real :: rgas                             !< Gas constant
  
  real :: nl                               !< Loschmidt density
  integer, parameter :: nraymax = 10       !< maximum number of molecules for Rayleigh scattering (10)
  real,dimension(nraymax) :: nr            !< Refractive indices of nraymax molecules for Rayleigh scattering
  real,dimension(nraymax) :: p_fac         !< Polarisation factor for nraymax molecules for Rayleigh scattering
  character(10),dimension(nraymax) :: rayname !< names of the nraymax molecules used for Rayleigh scattering

  real :: gram                             !< Gram unit
  real :: cm                               !< Centimetre unit
  real :: sec                              !< Second unit
  real :: kelvin                           !< Temperature unit
  real :: gcc                              !< Gram per cc unit
  real :: dyne                             !< Dyne unit
  real :: erg                              !< Erg unit
  
  real :: kg                               !< Kilogram unit
  real :: metre                            !< Metre unit
  real :: kg_per_m3                        !< Kilogram per m**3 unit
  real :: newton                           !< Newton unit
  real :: joule                            !< Joule unit
  real :: mole                             !< Mole unit

  real :: unitM                            !< Mass unit
  real :: unitL                            !< Length unit
  real :: unitT                            !< Time unit
  real :: unitK                            !< Temperature unit
  real :: unitD                            !< Density unit
  real :: unitF                            !< Force unit
  real :: unitE                            !< Energy unit
  real :: unitP                            !< Pressure unit
  real :: unitSF                           !< Spectral flux unit erg s-1 cm-2 sr-1 Hz-1
  real :: unitIF                           !< Integrated flux unit erg s-1 cm-2 sr-1
  real :: unitMo                           !< Mole unit erg K-1 mole-1
  
  real :: gas_constant_um                  !< Specific gas constant from the UM (SI)
  real :: specific_heat_um                 !< Specific heat constant from the UM (SI)
  real :: pp_ref_surf_um                   !< Reference surface pressure from the UM (SI)
  
  real :: one_over_L3                      !< inverse of unitL**3
  

contains

  !> Initialise the code units and the values of the physical constant in these units
  !<
  subroutine init_cst

    use mod_param

    implicit none

    integer :: stat

    namelist/units/unitM,unitL,unitT,unitK,gas_constant_um,specific_heat_um,pp_ref_surf_um

    ! cgs basis
    gram   = 1.0
    cm     = 1.0
    sec    = 1.0
    kelvin = 1.0
    gcc    = gram/cm**3
    dyne   = gram*cm/sec**2
    erg    = gram*cm**2/sec**2
    
    ! SI basis
    kg        = 1.0E-3*gram
    metre     = 1.0E-2*cm
    kg_per_m3 = kg/metre**3
    newton    = kg*metre/sec**2
    joule     = kg*metre**2/sec**2
    mole      = 1.0

    ! simulation basis
    unitM  = 1.0E17
    unitL  = 1.0E8
    unitT  = 1.0E2
    unitK  = 1.0
    unitMo = 1.0

    
    ! Constants used in UM simulation (given in SI)
    gas_constant_um = 3556.8
    specific_heat_um = 1.3e4
    pp_ref_surf_um = 200E5

    open(1,file=fparam)
    rewind(1) ; read(1,NML=units,iostat=stat)
     IF (stat/=0) THEN
      BACKSPACE(1)
      READ(1,fmt='(A)') line
      WRITE(*,'(A)') 'Invalid line in namelist: '//trim(line)
      STOP
     END IF
    close(1)

    if (mype == cputerm) then
       write(*,headrfmt) 'Units'
       write(*,exppmfmt) 'UnitL',unitL,'[cm] Length unit'
       write(*,exppmfmt) 'UnitM',unitM,'[g ] Mass unit'
       write(*,exppmfmt) 'UnitT',unitT,'[s ] Time unit'
       write(*,exppmfmt) 'UnitK',unitK,'[K ] Temp. unit'
       write(*,*)
       write(*,separfmt)
    endif

    unitD  = unitM/unitL**3
    unitF  = unitM*unitL/unitT**2
    unitE  = unitM*unitL**2/unitT**2
    unitP  = unitF/unitL**2
    unitSF = unitD*(unitL**4)/(unitT**3)
    unitIF = unitD*(unitL**3)/(unitT**3)

    gram   = 1.0/unitM
    cm     = 1.0/unitL
    sec    = 1.0/unitT
    kelvin = 1.0/unitK
    gcc    = 1.0/unitD
    dyne   = 1.0/unitF
    erg    = 1.0/unitE
    mole   = 1.0/unitMo
    
    kg        = 1.0E+3*gram
    metre     = 1.0E+2*cm
    kg_per_m3 = kg/metre**3
    newton    = kg*metre/sec**2
    joule     = kg*metre**2/sec**2
    
    ! physical constant

    kb      = 1.380620e-16 * erg/kelvin
    amu     = 1.660531e-24 * gram
    cvel    = 29979245800. * cm/sec
    hplanck = 6.626068e-27 * erg*sec
    sigsb   = 5.67051e-5   * erg/cm**2/kelvin**4/sec
    Atm     = 1.01325e06   * dyne/cm**2
    Rjup    = 6991100000.  * cm
    Rsun    = 69580000000. * cm
    au      = 14959787100000. * cm
    day     = 86400.       * sec
    rgas    = 8.3144598    * erg/kelvin/mole ! gas constant
    
    nl      = 2.687E19 * (1./cm**3)
    
    !! Refractive Indices and Polarisation factors for molecules used in Rayleigh scattering 
    
    rayname(1) = 'H2'  ; nr(1) = 1.0001384 ; p_fac(1) = 0.02    ! Leonard 1974 & Penndorf 1957 
    rayname(2) = 'He'  ; nr(2) = 1.000035  ; p_fac(2) = 0.02    ! Mansfield 1959 & Penndorf 1957  
    rayname(3) = 'CO'  ; nr(3) = 1.000325  ; p_fac(3) = 0.0095  ! Source Sneep & Ubachs 2005
    rayname(4) = 'N2'  ; nr(4) = 1.000284  ; p_fac(4) = 0.0205  ! Source Sneep & Ubachs 2005
    rayname(5) = 'CH4' ; nr(5) = 1.000481  ; p_fac(5) = 0.      ! Source Sneep & Ubachs 2005
    rayname(6) = 'NH3' ; nr(6) = 1.000375  ; p_fac(6) = 0.      ! Source Alan
    rayname(7) = 'H2O' ; nr(7) = 1.000254  ; p_fac(7) = 0.      ! Source Alan
    rayname(8) = 'CO2' ; nr(8) = 1.000427  ; p_fac(8) = 0.0790  ! Source Sneep & Ubachs 2005
    rayname(9) = 'H2S' ; nr(9) = 1.000634  ; p_fac(9) = 0       ! Source NPL
    rayname(10) = 'SO2' ; nr(10) = 1.000686 ; p_fac(10) = 0     ! Source NPL
    !rayname(11) = 'N2O' ; nr(11) = 1.000483  ; p_fac(11) = 0.1166  ! Source Sneep & Ubachs 2005
    !rayname(12) = 'SF6' ; nr(12) = 1.000733 ; p_fac(12) = 0.    ! Source Sneep & Ubachs 2005
    !rayname(13) = 'HCN' ; nr(13) = 1.0001384 ; p_fac(13) = 0.2 ! not found
    
    one_over_l3 = 1./unitL**3
    

  end subroutine init_cst

end module mod_cst
