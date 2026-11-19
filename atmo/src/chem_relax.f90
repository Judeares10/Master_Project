!> \file
!! Contains modules mod_chem_relax
!<

!> Holds subroutines to solve for the chemical steady-state   \n
!! using the chemical relaxation method
!<
MODULE mod_chem_relax

IMPLICIT NONE

! Number of molecules in relaxation scheme
INTEGER :: nmol_relax
! Relaxation scheme option
INTEGER :: relax_scheme
! Map of molecules in relaxation scheme in main molecules array
INTEGER,ALLOCATABLE :: imol_relax(:)
INTEGER,ALLOCATABLE :: inloc(:,:)

! Minimum time for integration
REAL :: tmin
! Convergence criteria
REAL :: dnn_conv 
REAL :: dnndt_conv
! Chemical timescales 
REAL,ALLOCATABLE :: tscale(:,:)
! Equilibrium number densities
REAL,ALLOCATABLE :: ndeq(:,:)
! Frequency at which to recalculate PT profile
INTEGER :: mod_pt
! Temperature below which the timescales are frozen at this temperature (currently just for tsai scheme)
REAL :: tfreeze_relax
! Name of the output PT profile if running self-consistently
CHARACTER(256) :: frelaxpt 

CONTAINS 

! Subroutine to initialise the chemical relaxation scheme
SUBROUTINE init_relax_chemistry

USE mod_grid, ONLY: ndepth
USE mod_chem, ONLY: imol, nmol
USE mod_param

IMPLICIT NONE

INTEGER :: stat

NAMELIST /chem_relax/ relax_scheme, tmin, dnn_conv, dnndt_conv, tfreeze_relax, mod_pt, frelaxpt

  ! Default namelist values
  relax_scheme = 1
  
  tmin       = 1.0e1
  dnn_conv   = 1.0e-2
  dnndt_conv = 1.0e-4
  tfreeze_relax = 0.0

  mod_pt = 10
  frelaxpt = 'frelaxpt.ncdf'
  
  !Read namelist
  OPEN(1,file=fparam)
  REWIND(1) ; READ(1,NML=chem_relax,iostat=stat)
   IF (stat/=0) THEN
    BACKSPACE(1)
    READ(1,fmt='(A)') line
    WRITE(*,'(A)') 'Invalid line in namelist: '//trim(line)
    STOP
   END IF
  CLOSE(1)

  ! Set the molecules included in the relaxation scheme
  ! depending on the option
  IF (relax_scheme==1) THEN
    ! Tsai et al (2017) relaxation method
    ! Six molecules: CH4, CO, H2O, CO2, N2 and NH3
    
    nmol_relax = 6
    ALLOCATE(imol_relax(nmol_relax))
    
    ! Set molecules to relax
    imol_relax(1) = imol('CH4')
    imol_relax(2) = imol('CO')
    imol_relax(3) = imol('H2O')
    imol_relax(4) = imol('CO2')
    imol_relax(5) = imol('N2')
    imol_relax(6) = imol('NH3')
  
  ELSE IF (relax_scheme==2) THEN
    ! Cooper and Showman (2006) relaxation method
    ! One molecule: CO
    
    nmol_relax = 1
    ALLOCATE(imol_relax(nmol_relax))
    
    ! Set molecules to relax
    imol_relax(1) = imol('CO')
  
  ELSE
    WRITE(*,*) 'Error: init_relax_chemistry'
    WRITE(*,*) 'relaxation scheme: ',relax_scheme, 'not supported'
    STOP
  END IF
  
  ! Allocate arrays
  ! Timescale for each molecule
  ALLOCATE(tscale(ndepth+1,nmol_relax))
  ! Equilibrium number densities of all molecules
  ALLOCATE(ndeq(ndepth+1,nmol))

  tscale = 0.
  ndeq   = 0.

END SUBROUTINE init_relax_chemistry

! ---------------------------------
! Main chemical relaxation routine
! ---------------------------------
SUBROUTINE get_relax_chemistry

IMPLICIT NONE

  ! Get the chemical timescales 
  CALL get_timescales

  ! Call the chemical relaxation solver
  CALL relaxation_solver
  
END SUBROUTINE get_relax_chemistry

! ------------------------------------------
! Subroutine to get the chemical timescales
! ------------------------------------------
SUBROUTINE get_timescales

  USE mod_chem, ONLY: nmol, Amolf, Amol, imol, get_eq_chemistry
  USE mod_grid, ONLY: ndepth,ppf,ttf
  USE mod_eos,  ONLY: get_nn
  USE mod_cst,  ONLY: unitL
  USE tsai2017_relaxation
  USE cs2006_relaxation 

  IMPLICIT NONE

  ! Local variables
  INTEGER :: i , j
  REAL :: number_density(ndepth+1)

  ! Get equilibrium abundances
  call get_eq_chemistry

  ! Calculate total number density (cgs)
  number_density = get_nn(ndepth+1,ppf,ttf)/unitL**3.

  ! Calculate equilibrium number densities
  DO j = 1, ndepth+1
    DO i = 1, nmol
      ndeq(j,i) = Amolf(j,i)*number_density(j)
    END DO
  END DO
  
  ! Calculate the chemical timescales
  IF (relax_scheme==1) THEN
    DO i = 1, nmol_relax
      CALL get_timescale_tsai2017(imol_relax(i),ndeq,tscale(:,i), tfreeze_relax)
    END DO
  ELSE IF (relax_scheme==2) THEN
    DO i = 1, nmol_relax
      CALL get_timescale_cs2006(ndeq,tscale(:,i))
    END DO
  END IF
    

END SUBROUTINE get_timescales


END MODULE mod_chem_relax

!----------------------
!Chemical relaxation solver routine
!----------------------
SUBROUTINE relaxation_solver

USE mod_grid, ONLY: ndepth, dr, drf, ttf, tt, pp, ppf, rho, rhof, kzz
USE mod_chem, ONLY: amol, amolf, imol, print_chem, molname, get_eq_chemistry, write_chemistry
USE mod_chem_relax
USE mod_eos, ONLY : get_rrgr, get_rho, get_drdrf
USE mod_cst, ONLY : unitL, unitT
USE mod_solver, ONLY: lsode_init, solve_dlsodes, rwork, solve_hydro, solve_energy, find_profile, lsode_istate
USE mod_param, ONLY: mype, cputerm, headrfmt, fout
USE tsai2017_relaxation
USE cs2006_relaxation

IMPLICIT NONE

EXTERNAL dndt_relax, jac_relax

INTEGER :: i, nit, im,idepth,node
REAL :: time, dt
REAL,ALLOCATABLE :: nd(:),nd_old(:),nd_old2(:)
REAL,PARAMETER :: acarbon = 4.570e-4
! Oxygen abundance (f_H2O + f_CO)
REAL,PARAMETER :: aoxygen = 8.3104e-4
CHARACTER(5)  :: fnit


  real :: dnn,dnnmax
  real :: dnndt,dnndtmax
  real :: time_save, time_save_old

  ! Initialise variables
  time = 0.
  nit = 0
  dt = 1e-10
  dnnmax = 1.0
  dnndtmax = 1.0
  
  ! Allocate working arrays
  node = (ndepth+1)*nmol_relax
  ALLOCATE(nd(node))
  ALLOCATE(nd_old(node))
  ALLOCATE(nd_old2(node))
  ALLOCATE(inloc(ndepth+1,nmol_relax))
  
  nd = 0. ; nd_old = 0. ; nd_old2 = 0.
  
  ! Initialise to equilibrium
  DO idepth = 1,ndepth+1
    DO im = 1, nmol_relax
      inloc(idepth,im) = (idepth-1)*nmol_relax+im
      nd(inloc(idepth,im)) = ndeq(idepth,imol_relax(im))
    END DO
  END DO
  
  nd_old  = nd
  nd_old2 = 0.
  
  ! Setup atmosphere
  ! Calculate mass densities, layer thicknesses
  rho = get_rho(ndepth,pp,tt,amol)
  rhof = get_rho(ndepth+1,ppf,ttf,amolf)
  
  ! Get radius and gravity of grid
  CALL get_rrgr
  
  !Calculate thickness of layers (cgs)
  CALL get_drdrf
  dr = dr*unitL
  drf = drf*unitL

  ! Initialise LSODE solver
  CALL lsode_init(node)

  time_save = dt*1.2
  time_save_old = dt
  
  ! Convert to cgs
  kzz   = kzz*unitL**2/unitT

  ! Loop until reached convergence criteria are met
  DO WHILE (dnnmax>dnn_conv .or. dnndtmax>dnndt_conv .or. time<tmin)

    ! Reset variations
    dnn = 0.
    dnndt = 0.
    dnnmax = 0.
    dnndtmax = 0.
    
    ! Save previous iterations
    if (time>time_save) then
      nd_old = nd
      nd_old2 = nd_old
      time_save_old = time_save
      time_save = time_save*1.2
    end if
    
    ! Update equilibrium values if required
    IF (relax_scheme==1) THEN
      ndeq(:,imol('CO2')) = eq_co2(nmol_relax,node,inloc,imol_relax,ndeq,nd)
    END IF

    ! Call dlsodes solver
    CALL solve_dlsodes(node,nd,time,time+dt,dndt_relax,jac_relax)

    if (lsode_istate<0) then 
      lsode_istate = 1
    endif

    ! Timestep in next iteration
    dt = 0.9*rwork(11)
    
    ! Advance number of iterations
    nit = nit + 1
    
    ! Calculate variations: dn/n and (dn/n)/dt
    ! and find maximum variations
    DO i = 1, node
        IF (nd(i)>0. .AND. nd_old2(i)>0.) THEN
          dnn =  abs((nd_old2(i)-nd(i)))/nd_old2(i)
          dnndt = dnn/(time-time_save_old)
          IF (dnn>dnnmax) THEN
            dnnmax = dnn
          END IF
          IF (dnndt>dnndtmax) THEN
            dnndtmax = dnndt
          END IF
        ELSE IF (nit == 1) THEN
          dnndtmax = 1.0
          dnnmax = 1.0
        END IF
    END DO
    
    ! Save solution
    nd_old = nd
    
    ! Copy equilibrium number densities
    Amolf = ndeq
    ! Copy relaxation number densities
    DO i = 1, node
      idepth = (i-1)/nmol_relax+1
      im     = i-(idepth-1)*nmol_relax
      Amolf(idepth,imol_relax(im)) = nd(i)    
    END DO

    ! Calculate mole fraction
    DO idepth = 1, ndepth+1
      Amolf(idepth,:) = Amolf(idepth,:)/sum(Amolf(idepth,:)) 
    END DO

    ! If Cooper and Showman 2006, update CH4 and H2O
    IF (relax_scheme==2) THEN
      Amolf(:,imol('CH4')) = acarbon - Amolf(:,imol('CO'))
      Amolf(:,imol('H2O')) = aoxygen - Amolf(:,imol('CO'))
    END IF

    !interpolate to get Amol at cell center, using geometric mean
    DO idepth = 1,ndepth
      Amol(idepth,:) = sqrt(Amolf(idepth,:)*Amolf(idepth+1,:))
    END DO
  
    ! Print to screen
    IF (print_chem .and. mype == cputerm) THEN
      WRITE(fnit,'(i5)') nit
      WRITE(*,headrfmt) 'Chem. Iter.: '//fnit 
      WRITE(*,999) time,dt
      WRITE(*,1000) dnnmax,dnndtmax
      WRITE(*,1001) molname(imol('CH4')),Amolf(1,imol('CH4')),molname(imol('H2O')),Amolf(1,imol('H2O')), &
      molname(imol('CO2')),Amolf(1,imol('CO2'))
      WRITE(*,*)
    END IF

999   FORMAT(5x,' time: ',es10.1,5x,' dt: ',es10.1)
1000  FORMAT(5x,' max(dn/n): ',es10.1,5x,' max((dn/n)/dt): ',es10.1)
1001  FORMAT(5x,a10,': ',es10.2,', ',a10,': ',es10.2,', ',a10,': ',es10.2)

    ! Write abundances to file
    CALL write_chemistry 

    ! If requested, reconverge the PT profile
    if (mod(nit, mod_pt)==0) then 
      if (solve_hydro .or. solve_energy) then 
        if (mype==cputerm) write(*,*) 'Reconverging the PT profile'
        call find_profile('cst')
        if (dnnmax>dnn_conv .or. dnndtmax>dnndt_conv .or. time<tmin) then 
          call get_timescales
        endif
        lsode_istate=1
      endif 
    endif 

  END DO
  
  ! Print to screen
  IF (mype == cputerm) THEN
    WRITE(*,headrfmt) 'Convergence criteria met'    
  END IF

  ! Reconverge PT profile one last time
  if (solve_hydro .or. solve_energy) then 
    if (mype==cputerm) write(*,*) 'Reconverging the PT profile'
    call find_profile('cst')
    call write_pt_profile(frelaxpt)
  endif 

END SUBROUTINE relaxation_solver

!--------------------------------------------
! Compute dn/dt for the relaxation solver
!--------------------------------------------
SUBROUTINE dndt_relax(lsode_neq,t,nd,dndt)

USE mod_chem_relax
USE mod_chem_neq, ONLY: get_phi
USE mod_grid, ONLY: ndepth, drf, mixing

IMPLICIT NONE

! Number of ODEs
INTEGER,INTENT(IN) :: lsode_neq
! Current time
REAL,INTENT(IN) :: t
! Current number densities
REAL,INTENT(IN) :: nd(lsode_neq)
! Variations
REAL,INTENT(INOUT) :: dndt(lsode_neq)

! Local variables
INTEGER :: idepth,im
REAL :: phi(ndepth,nmol_relax) 

dndt(1:lsode_neq) = 0.

  DO idepth = 1, ndepth+1
    DO im = 1, nmol_relax
      dndt(inloc(idepth,im)) = -(nd(inloc(idepth,im))-ndeq(idepth,imol_relax(im)))/tscale(idepth,im)
    END DO
  END DO
  
  ! Mixing
  IF (mixing) THEN
    ! Calculate phi
    DO im = 1, nmol_relax
      phi(:,im) = get_phi(im,imol_relax(im),lsode_neq,nmol_relax,inloc,nd)
    END DO

    ! Add diffusion term to dndt
    DO idepth = 2, ndepth-1
      dndt(inloc(idepth,1:nmol_relax)) = dndt(inloc(idepth,1:nmol_relax)) - &
        (phi(idepth,1:nmol_relax)- phi(idepth-1,1:nmol_relax))/drf(idepth)    
    END DO
  
    ! Boundaries
    dndt(inloc(1,:)) = dndt(inloc(1,:)) - (phi(1,:)-0.)/drf(1)
    dndt(inloc(ndepth+1,:)) = dndt(inloc(ndepth+1,:)) - (0.-phi(ndepth,:))/drf(ndepth+1)
  END IF


  
END SUBROUTINE dndt_relax

!----------------------
! Dummy subroutine for LSODE
!----------------------
SUBROUTINE jac_relax

IMPLICIT NONE

  ! Dummy routine

END SUBROUTINE jac_relax



