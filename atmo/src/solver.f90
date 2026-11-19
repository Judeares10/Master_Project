!> \file
!! Contains module mod_solver\n
!! Compute the jacobi matrix and the right and side of the system of equation
!<

!> hold the parameters of the solver
!<
module mod_solver

  real,dimension(:,:),allocatable :: jacobi !< jacobi matrix
  real,dimension(:  ),allocatable :: rhs    !< right hand side

  real :: tau_heateq !< optical depth below which the heat equation is solved instead of the flux equation

  real :: minstep !< minimum step for the raphson-newton
  real :: maxstep !< maximum step for the raphson-newton

  logical :: solve_energy    !< switch to solve the energy equation, if = F fix the temperature profile
  logical :: solve_hydro     !< switch to solve the hydrostatic equation, if = F fix the pressure profile
  logical :: include_dkapdpt !< include the derivative of the opacity in the jacobi matrix
  logical :: surface_spectrum !< compute the surface emission spectrum of the body
  logical :: surface_spectrum_3d !< Compute the surface emission spectrum of the body using 3D ouutput from the UM
  logical :: phase_curve !< Compute phase curves from UM output
  logical :: cloud_night_side !< Use clouds for the night side for phase curves
  logical :: transmission_spectrum !< Compute the transmission spectrum of the planet
  logical :: transmission_spectrum_3d !< Compute the transmission spectrum of the planet using 3D output from the UM
  logical :: hydrostatic   !< switch to constrain the standard tau (and change taumax) to satisfy \n
                                      !the hydrostatic equation if the solver is not used to constrain pressure
  logical :: include_uv !< include uv irradiation
  real :: accuracy !< requested accuracy for the profile
  integer :: nitmax !< Maximum number of iteration for the solver
  integer :: nbin !< number of frequency per bins for the perturbations in the jacobi matrix
  real :: psurf     !< pressure condition at the top of the atmosphere
  real :: unitPsurf !< units of the pressure condition
  logical :: print_err !< flag to print on terminal the error in the temp profile
  logical :: cst_h2o   !< flag to keep h2o cst in jacabian matrix
  real :: step
  logical :: first
  logical :: rainout_freeze !< flag to freeze the rainout model level of a given element
  integer :: rainout_freeze_nit !< number of iterations after which rainout model level is frozen
  integer,dimension(:),allocatable :: rainout_levels !< array to keep past rainout model levels
  integer :: rainout_freeze_past !< number of past iterations used to determine the rainout model level
  real :: rainout_freeze_dpt !< dpt/pt threshold below which rainout model level is frozen
  real :: pp2d_decoupling_momt !< pressure under which we decouple the pressure gradient from longitudinal momentum
  real :: pp2d_decoupling_momr !< pressure under which we decouple the velocity gradient from vertical momentum
  real :: pp2d_decoupling_err  !< pressure under which we do not check errors
  real :: pp1d_decoupling_err_hydro !< Pressure below which we do not check errors in hydrostatic equilibrium
  real :: unitPdecoupling

  real,dimension(:,:),allocatable  :: fsurf_2d     !< surface radiative flux in 2d
  integer :: np_phase_curve

  logical :: calc_cf !< Logical for computing the Contribution Function
  logical :: calc_radtime !< Logical for computing the radiative timescale

  ! LSODE parameters
  INTEGER :: lsode_order         !< order of problem for lsode
  INTEGER :: lsode_nit           !< maximum number of (internal) steps for lsode
  INTEGER :: lsode_print_warning !< maximum number of warnings printed to screen
  INTEGER :: lsode_imf           !< method flag for lsode
  INTEGER :: lsode_istate        !< Integer flag - 0 recompute Jacobian, 1 normal calc
  INTEGER :: lsode_itask         !< itask = 1
  INTEGER :: lsode_itol          !< tolerance flag
  INTEGER :: lsode_iopt          !< optional inputs
  INTEGER :: liw                 !< Length of integer work array
  INTEGER :: lrw                 !< Length of real work array
  INTEGER,ALLOCATABLE :: iwork(:)!< Integer work array
  REAL    :: lsode_atol          !< absolute tolerance
  REAL    :: lsode_rtol          !< relative tolerance
  REAL,ALLOCATABLE :: rwork(:)   !< Real work array
  LOGICAL :: lsode_dt, hminus    !< Use internally generated (lsode) timestep
  LOGICAL :: change_minstep      !< Flag to change minstep after a given number of iterations
  INTEGER :: nit_minstep         !< Number of iterations after which to change minstep
  REAL    :: minstep_factor      !< Factor by which to multiply minstep by after a given number of iterations
  REAL    :: dpt_limit           !< Limit to changes of P, T above which
                                 !< stepsize is set to minstep

  REAL    :: temp_coupling !< coupling parameter for the temperature discretization in the newton iterations
  REAL    :: ppmin_temp_coupling !< min pressure above which a temperature coupling is applied
  REAL    :: ppmax_temp_coupling !< max pressure below which a temperature coupling is applied
  REAL    :: unitP_temp_coupling

  contains

    !> Initialize the solver to find the correct pressure/temp profile
    !<
    subroutine init_solver

      use mod_grid
      use mod_param
      use mod_opacity
      use mod_cst
      use mod_radtrans
      use mod_chem
      use mod_conv

      implicit none

      integer :: stat,it


      namelist/solver/solve_energy,solve_hydro,accuracy,nitmax,tau_heateq,nbin,minstep,maxstep, &
      include_dkapdpt,surface_spectrum,surface_spectrum_3d,phase_curve,cloud_night_side,transmission_spectrum, &
      transmission_spectrum_3d,hydrostatic,psurf,unitPsurf,print_err,include_uv,cst_h2o,calc_cf, &
      np_phase_curve,first,rainout_freeze,rainout_freeze_nit,rainout_freeze_dpt,rainout_freeze_past, &
      lsode_order, lsode_nit, lsode_print_warning, lsode_imf, lsode_atol, lsode_rtol, lsode_dt, hminus, &
      change_minstep, nit_minstep, minstep_factor,pp2d_decoupling_momt,pp2d_decoupling_momr,pp2d_decoupling_err,unitPdecoupling, &
      pp1d_decoupling_err_hydro,calc_radtime, dpt_limit,&
      temp_coupling, ppmin_temp_coupling, ppmax_temp_coupling,unitP_temp_coupling

      np_phase_curve = 200
      first = .true.
      !default values
      solve_energy = .False.
      solve_hydro  = .False.
      surface_spectrum = .False.
      surface_spectrum_3d = .False.
      phase_curve = .False.
      cloud_night_side = .False.
      transmission_spectrum = .False.
      transmission_spectrum_3d = .False.
      include_uv = .false.
      accuracy = 1.0E-3
      nitmax = 1000
      tau_heateq = 1.E-40
      nbin = 1
      minstep = 0.1
      maxstep = 0.9
      include_dkapdpt = .False.
      hydrostatic = .True.
      psurf = 1.0E-6
      unitPsurf = 1.0E6 !bar
      print_err = .false.
      cst_h2o = .false.
      calc_cf = .false.
      calc_radtime = .false.
      rainout_freeze = .false.
      rainout_freeze_nit = 250
      rainout_freeze_past = 20
      rainout_freeze_dpt = 5.0E-2
      hminus = .false.
      change_minstep = .false.
      nit_minstep = 500
      minstep_factor = 10.0
      dpt_limit = 1E1
      pp2d_decoupling_momt = 1E-10
      pp2d_decoupling_momr = 1E-10
      pp2d_decoupling_err  = 1E-10
      pp1d_decoupling_err_hydro = 1E-10
      unitPdecoupling = 1.0E6  !bars

      !temperature coupling
      temp_coupling =0.
      ppmin_temp_coupling =1E-10
      ppmax_temp_coupling =1E-10
      unitP_temp_coupling = 1E6 !bar

      ! LSODE default parameters
      lsode_order         = 2
      lsode_nit           = 200
      lsode_print_warning = 0
      lsode_imf           = 222
      lsode_atol          = 1.0e-10
      lsode_rtol          = 1.0e-3
      lsode_dt            = .True.

      ! Read parameters
      open(unit=1,file=fparam)
      read(unit=1,nml=solver,iostat=stat)
       IF (stat/=0) THEN
        BACKSPACE(1)
        READ(1,fmt='(A)') line
        WRITE(*,'(A)') 'Invalid line in namelist: '//trim(line)
        STOP
       END IF
      close(1)

      ppmin_temp_coupling = ppmin_temp_coupling*unitP_temp_coupling/unitP
      ppmax_temp_coupling = ppmax_temp_coupling*unitP_temp_coupling/unitP

      pp2d_decoupling_momr = pp2d_decoupling_momr*unitPdecoupling/unitP
      pp2d_decoupling_momt = pp2d_decoupling_momt*unitPdecoupling/unitP
      pp2d_decoupling_err  = pp2d_decoupling_err *unitPdecoupling/unitP
      pp1d_decoupling_err_hydro  = pp1d_decoupling_err_hydro *unitPdecoupling/unitP

      if (ndim==1) then
         if(solve_energy .or. solve_hydro) then
            allocate(jacobi(2*(ndepth+1),2*(ndepth+1)),rhs(2*(ndepth+1)))
         endif
      else if (ndim==2) then
         allocate(jacobi(5*(ndepth+1)*(ntheta+1),5*(ndepth+1)*(ntheta+1)),rhs(5*(ndepth+1)*(ntheta+1)))
      endif

      allocate(rainout_levels(rainout_freeze_past))
      rainout_levels = 0

      !units
      Psurf = Psurf*unitPsurf/unitP

      !snu_scat (save the scattering source term to accelerate the convergence)
      if (ndim==1 .or. um_input) then
         if (nfreqkap>0) then
            allocate(snu_scat(ndepth+1,nfreqkap))
            snu_scat = get_bnu(ndepth+1,nfreqkap,ttf)
         endif

      else

         if (nfreqkap>0) then
            allocate(snu_scat_2d(ndepth+1,ntheta+1,nfreqkap))
            do it=1,ntheta+1
               snu_scat_2d(:,it,:) = get_bnu(ndepth+1,nfreqkap,ttf)
            enddo
         endif

      endif

      if (mype==cputerm) then
       write(*,headrfmt)'Solver'
       write(*,logpmfmt)'solve_hydro',solve_hydro,'Perturb the hydrostatic equation'
       write(*,logpmfmt)'solve_energy',solve_energy,'Perturb the energy equation'
       write(*,exppmfmt)'psurf',psurf*unitP/1.0E6,'[bar] surface pressure'
       write(*,exppmfmt)'accuracy',accuracy,'Requested accuracy for the Newton-Raphson method'
       write(*,intpmfmt)'nitmax',nitmax,'Maximum number of iteration'
       !write(*,intpmfmt)'nbin',nbin,'Frequency binning for the jacobi matrix'
       write(*,exppmfmt)'minstep',minstep,'Minimum step size of the Newton-Raphson'
       write(*,exppmfmt)'maxstep',maxstep,'Maximum step size of the Newton-Raphson'
       !write(*,exppmfmt)'tau_heateq',tau_heateq,'Switch to Qrad = 0 below tau_heateq'
       if (rainout_freeze) then
         write(*,*)
         write(*,logpmfmt)'rainout_freeze',rainout_freeze,'Freeze rainout of element (set in chemistry mod)'
         write(*,intpmfmt)'rainout_freeze_nit',rainout_freeze_nit,'No. of iterations after which rainout is frozen'
     	 write(*,intpmfmt)'rainout_freeze_past',rainout_freeze_past,'No. of past iterations used to set rainout freeze level'
       	 write(*,exppmfmt)'rainout_freeze_dpt',rainout_freeze_dpt,'dpt threshold below which rainout is frozen'
       endif
       write(*,*)
       write(*,separfmt)
    endif




    if(corr_k .and. include_dkapdpt) then
       if (mype==cputerm) then

          write(*,headrfmt) 'Error init_solver'
         write(*,'(A)') 'kappa derivatives with correlated-k not implemented'
         write(*,logpmfmt) 'corr_k',corr_k
         write(*,logpmfmt) 'include_dkapdpt',include_dkapdpt
         write(*,*)
         write(*,separfmt)
         stop

      endif
   endif

       if(corr_k .and. nbin>1) then
       if (mype==cputerm) then

          write(*,headrfmt) 'Error init_solver'
         write(*,'(A)') 'Cannot use freq. reduction with correlated-k'
         write(*,logpmfmt) 'corr_k',corr_k
         write(*,intpmfmt) 'nbin',nbin
         write(*,*)
         write(*,separfmt)
         stop

      endif
   endif

      if(solve_energy .and. nkap==0) then
       if (mype==cputerm) then

          write(*,headrfmt) 'Error init_solver'
          write(*,'(A)') 'Cannot solve energy equation without opacities'
          write(*,logpmfmt) 'solve_energy',solve_energy
          write(*,intpmfmt) 'nkap',nkap
          write(*,*)
          write(*,separfmt)
         stop

      endif
   endif

   ! Phase curves can only be calculated if um_input == .true.
   if (.not.um_input .and. phase_curve) then
     write(*,headrfmt) 'Error init_solver'
     write(*,'(A)') 'Phase curves can only be calculated if UM output is provided.'
     write(*,logpmfmt) 'phase_curve',phase_curve
     write(*,logpmfmt) 'um_input',um_input
     write(*,*)
     write(*,separfmt)
     stop
   endif

   ! Irradiation must be turned on if calculating phase curves
   if (.not.irrad .and. phase_curve) then
     write(*,headrfmt) 'Error init_solver'
     write(*,'(A)') 'Phase curves can only be calculated if radiation is on.'
     write(*,logpmfmt) 'phase_curve',phase_curve
     write(*,logpmfmt) 'irrad',irrad
     write(*,*)
     write(*,separfmt)
     stop
   endif

   ! calc_cp_chem must be turned on while calculating heating rates
   if (.not.calc_cp_chem .and. calc_cf) then
     write(*,'(A)') 'Heating Rates can only be calculated if calc_cp_chem is true.'
     stop
  endif
  
  ! calc_cp_chem must be turned on while calculating heating rates
  if ((.not. alpha==0.) .and. calc_radtime) then
     write(*,'(A)') 'Error solver.f90: Radiative timescale should be computed without convection, set alpha=0'
     stop
  endif

   ! calc_cp_chem must be turned on while calculating heating rates
  if ((.not.(solve_energy .and. solve_hydro)) .and. calc_radtime) then
     write(*,'(A)') 'Error solver.f90: Radiative timescales should be computed on an equilibrium profile, set solve_energy=solve_hydro=T'
     stop
  endif

  if(calc_radtime) then
     allocate(radtime(ndepth+1))
     radtime = 0.
  endif

   if (elename_freeze == 'None' .and. rainout_freeze) then
     if (mype == cputerm) then
       WRITE(*,headrfmt) 'Error init_solver'
       WRITE(*,'(A)') 'Set elename_freeze in &CHEMISTRY when using rainout_freeze'
       WRITE(*,*)
       WRITE(*,separfmt)
       stop
     endif
   endif


 end subroutine init_solver



    !> fill the jacobi matrix and the rhs for the hydrostatic equation \n
    !<
    subroutine perturb_hydro(err_hydro, err_hydro_save)

      use mod_grid
      use mod_param

      implicit none

      integer :: i
      real    :: dtau,err_hydro,kapstd_loc
     real,dimension(:),INTENT(OUT) :: err_hydro_save

      err_hydro = 0.
      err_hydro_save(1) = 0.

      !boundary condition dp0  = 0. do nothing for i == 1

      !interior, write the hydrostatic equation at cell centre i in jacobi(i+1,:)
      do i=1,ndepth

         dtau = taufstd(i+1) - taufstd(i)

         !if(i==1 .or. i==ndepth) then
            kapstd_loc = kapstd(i)
         !else
         !   kapstd_loc = (kapstd(i-1)+kapstd(i)+kapstd(i+1))/3.
         !endif

         if (kapstd_loc /=0.) then
            ! pressure perturbation
            jacobi(i+1,i  ) = -1.0/dtau + gr(i)/kapstd_loc**2* &
                 & dkapstd(i)*0.5*sqrt(ppf(i+1)/ppf(i))
            jacobi(i+1,i+1) =  1.0/dtau + gr(i)/kapstd_loc**2* &
                 & dkapstd(i)*0.5*sqrt(ppf(i)/ppf(i+1))

            !temperature perturbation
            jacobi(i+1,ndepth+1+i  ) = gr(i)/kapstd_loc**2* &
                 & dkapstd(ndepth+i)*0.5*sqrt(ttf(i+1)/ttf(i))
            jacobi(i+1,ndepth+1+i+1) = gr(i)/kapstd_loc**2* &
                 & dkapstd(ndepth+i)*0.5*sqrt(ttf(i)/ttf(i+1))


            !rhs
            rhs(i+1) = -(ppf(i+1)-ppf(i))/dtau+gr(i)/kapstd_loc

            !error hydrostatic equation
            if (ppf(i)>pp1d_decoupling_err_hydro) then
              if (i.gt.0) err_hydro = max(err_hydro,abs(((ppf(i+1)-ppf(i))/dtau-gr(i)/kapstd_loc)*kapstd_loc/gr(i)))
            endif
            err_hydro_save(i+1) = abs(((ppf(i+1)-ppf(i))/dtau-gr(i)/kapstd_loc)*kapstd_loc/gr(i))
            !err_hydro_save(i) = abs(((ppf(i+1)-ppf(i))/dtau-gr(i)/kapstd_loc)*kapstd_loc/gr(i))
          !  if (mype==0) write(*,"(i3,E10.2,E10.2)") i,taufstd(i),abs(((ppf(i+1)-ppf(i))/dtau-gr(i)/kapstd(i))*kapstd(i)/gr(i))

         endif

      enddo


    end subroutine perturb_hydro


    !> fill the jacobi matrix and the rhs for the energy equation \n
    !<
    subroutine perturb_energy(err_energy,err_energy_save,nit)

      use mod_cst
      use mod_grid
      use mod_conv
      use mod_radtrans
      use mod_opacity
      use mod_eos
      use mod_param
      use mod_chem

      implicit none

      integer :: i,i2,ii,nfreqjac,idepth,ib

      real,dimension(:,:),allocatable :: jnu,jnup,jnud,hnu,hnup,hnud,bnu,dsnu,dhnup,dhnu,dhnud,djnu,djnup,djnud,dbnu,dbnu_tot,kapjac,bjac,jjac,jpjac,jdjac,hjac,hpjac,hdjac,kapfjac,dkapfjac,epsfjac,snu,sjac,dsnu_scat
      real,dimension(:  ),allocatable :: dfrad,drhodtt,dqrad,qbrad,nujac,ah2o_save,anh3_save,afh2o_save,afnh3_save, afe_temp, affe_temp
      real,dimension(:),INTENT(OUT) :: err_energy_save
      integer,INTENT(IN) :: nit

      real,dimension(1) :: int_fsurf
      real,dimension(1,nfreqkap) :: bnusurf
      real :: fsurf,fint,err_energy,dH2,ddust,Mdust
      real,dimension(ndepth+1) :: kapfross,kapfplanck,muf,gfad,cp,gam,K_mu,K_T

      kapfross = get_kapross(ndepth+1,ttf,kapfnu)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! block computing the condition for fingering convection (brown dwarf paper)
      !!! in eos.f90 get_mu() nmolg should be replaced with nmol
      !!! to be used on a profile with no convection

      ! kapfplanck = get_kapplanck(ndepth+1,ttf,kapfnu)

      ! muf = get_mu(ndepth+1,amolf)
      ! gfad = get_gradad(ndepth+1,ppf,ttf,amolf)
      ! cp = get_cpspec(ndepth+1,ppf,ttf,amolf)/kb*muf*amu
      ! gam = cp/(cp-1.)
      ! cp = get_cpspec(ndepth+1,ppf,ttf,amolf)

      ! dH2 = 2.5E-8/unitL

      ! Mdust = Molmass(147)
      ! ddust = dmol(147)

      ! !K_T  = 1./(3.*pi)*(kb*ttf)**(3./2.)/ppf*1./(dH2**2*gam*sqrt((gam-1.)*muf*amu))
      ! K_T  = 16.*sigsb*teff**3/(3.*max(kapfross,kapfplanck)*cp*rhof**2)
      ! K_mu = 1./3.*sqrt(3.*kb*ttf/(Mdust*amu))*1./(sqrt(2.)*pi*(ddust)**2)*kb*ttf/ppf

      ! do idepth = 2,ndepth

      !    write(*,*)  idepth,ppf(idepth)*unitP/1.0E6,ttf(idepth),(((log(ttf(idepth+1))-log(ttf(idepth-1)))/(log(ppf(idepth+1))-log(ppf(idepth-1))))-gfad(idepth))/((log(muf(idepth+1))-log(muf(idepth-1)))/(log(ppf(idepth+1))-log(ppf(idepth-1)))),K_mu(idepth)*unitL**2/unitT,K_T(idepth)*unitL**2/unitT,K_T(idepth)/K_mu(idepth)
      ! enddo
      ! stop
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      err_energy = 0.

      if (mod(nfreqkap,nbin) > 0 ) then

         write(*,headrfmt) 'Error perturb_energy'
         write(*,'(A)') 'The local number of frequencies has to be a multiple of nbin'
         write(*,intpmfmt) 'nfreqkap',nfreqkap
         write(*,intpmfmt) 'nbin',nbin
         write(*,intpmfmt) 'mype',mype
         write(*,*)
         write(*,separfmt)
         stop

      endif

      nfreqjac = nfreqkap/nbin

      !allocate moments of intensity and source function
      allocate(hnu(ndepth+1,nfreqkap),jnu(ndepth+1,nfreqkap),snu(ndepth+1,nfreqkap))
      allocate(hnup(ndepth+1,nfreqkap),jnup(ndepth+1,nfreqkap))
      allocate(hnud(ndepth+1,nfreqkap),jnud(ndepth+1,nfreqkap))
      allocate(hjac(ndepth+1,nfreqjac),jjac(ndepth+1,nfreqjac),bjac(ndepth+1,nfreqjac),sjac(ndepth+1,nfreqjac))
      allocate(hpjac(ndepth+1,nfreqjac),jpjac(ndepth+1,nfreqjac))
      allocate(hdjac(ndepth+1,nfreqjac),jdjac(ndepth+1,nfreqjac))
      allocate(dhnup(ndepth+1,nfreqjac),djnup(ndepth+1,nfreqjac))
      allocate(dhnud(ndepth+1,nfreqjac),djnud(ndepth+1,nfreqjac))
      allocate(dhnu(ndepth+1,nfreqjac),djnu(ndepth+1,nfreqjac))
      allocate(dsnu(ndepth+1,nfreqjac),bnu(ndepth+1,nfreqkap))
      allocate(dbnu(ndepth+1,nfreqjac),qbrad(ndepth+1))
      allocate(dfrad(ndepth+1),drhodtt(ndepth+1),dqrad(ndepth+1))
      allocate(dbnu_tot(ndepth+1,nfreqkap),kapjac(ndepth,nfreqjac))
      allocate(dkapfjac(2*(ndepth+1),nfreqjac),kapfjac(ndepth+1,nfreqjac))
      allocate(nujac(nfreqjac+1),epsfjac(ndepth+1,nfreqjac))
      allocate(afh2o_save(ndepth+1),afnh3_save(ndepth+1),ah2o_save(ndepth),anh3_save(ndepth), afe_temp(ndepth), affe_temp(ndepth+1))

      !get convective flux at interfaces and its derivatives
      if(mype==cputerm .and. debug>0.) write(*,*) 'start get_convection'
      call get_convection
      if(mype==cputerm .and. debug>0.) write(*,*) 'end get_convection'


       do idepth=1,ndepth+1

         if (fconv(idepth)>0. .and. taufstd(idepth) < tau_heateq) then

            if (mype==cputerm) then

               write(*,headrfmt) 'Error perturb_energy'
               write(*,'(A)') 'Convection below tau_heateq, cannot solve the heat equation when convection is present'
               write(*,intpmfmt) 'idepth',idepth
               write(*,exppmfmt) 'fconv',fconv
               write(*,exppmfmt) 'tau_heateq',tau_heateq
               write(*,*)
               write(*,separfmt)
               stop

            endif

         endif


      enddo


      drhodtt = get_drhodt(ndepth+1,ppf,ttf,Amolf)

      !compute planck function derivative
      bnu      = get_bnu (ndepth+1,nfreqkap,ttf)
      dbnu_tot = get_dbnu(ndepth+1,nfreqkap,ttf)

      fstar_nu = 0.

      !compute star flux
      if (irrad) then
         call startrans(nfreqkap,kapnu,kapstd,fstar_nu,Irad)
      endif

      !source term for RT
      if (scatter) then
         snu = epsfnu*bnu+(1.-epsfnu)*abs(fstar_nu/murad)/(4.*pi)
      else
         !compute planck function
         snu      = bnu
      endif

	  if (corr_k) then
      	where(weightsfnu == 0.) snu = 0.
	  endif



      if (nbin == 1) then

         dbnu   = dbnu_tot
         kapjac = kapnu
         nujac  = nuf
         epsfjac = epsfnu
         bjac   = bnu
         dkapfjac = dkapfnu
         kapfjac = kapfnu

      else

         ii = 1
         nujac(1) = nuf(1)
         do i = 1,nfreqjac

            dbnu(:,i) = sum(dbnu_tot(:,ii:ii+nbin-1)**2,2)/sum(dbnu_tot(:,ii:ii+nbin-1),2)

            kapfjac(:,i) = sum(dbnu_tot(:,ii:ii+nbin-1)*kapfnu(:,ii:ii+nbin-1),2)/sum(dbnu_tot(:,ii:ii+nbin-1),2)

            dkapfjac(1:ndepth+1,i) = sum(dbnu_tot(:,ii:ii+nbin-1)*dkapfnu(1:ndepth+1,ii:ii+nbin-1),2)/sum(dbnu_tot(:,ii:ii+nbin-1),2)

            dkapfjac(ndepth+2:2*(ndepth+1),i) = sum(dbnu_tot(:,ii:ii+nbin-1)*dkapfnu(ndepth+2:2*(ndepth+1),ii:ii+nbin-1),2)/sum(dbnu_tot(:,ii:ii+nbin-1),2)

            bjac(:,i) = sum(dbnu_tot(:,ii:ii+nbin-1)*bnu(:,ii:ii+nbin-1),2)/sum(dbnu_tot(:,ii:ii+nbin-1),2)

            epsfjac(:,i) = sum(epsfnu(:,ii:ii+nbin-1),2)/nbin

            ii = ii + nbin
            nujac(i+1) = nuf(ii)

         enddo

         dbnu_tot = get_dbnu(ndepth,nfreqkap,tt)

         ii = 1

         do i = 1,nfreqjac

            kapjac(:,i) = sum(dbnu_tot(1:ndepth,ii:ii+nbin-1)*kapnu(:,ii:ii+nbin-1),2)/sum(dbnu_tot(1:ndepth,ii:ii+nbin-1),2)

            ii = ii + nbin

         enddo

      endif

      !compute radiative flux at interfaces
      if(mype==cputerm .and. debug>0) write(*,*) 'start radtrans'
      call radtrans(nrays,scatter,nfreqkap,kapnu,kapstd,epsfnu,snu,bnu,hnup,hnud,jnup,jnud,snu_scat)
      if(mype==cputerm .and. debug>0) write(*,*) 'end radtrans'



      hnu = hnup+hnud
      jnu = jnup+jnud

      !compute the surface flux
      bnusurf = get_bnu(1,nfreqkap,(/ Teff /))
      int_fsurf = int_freq(1,nfreqkap,nuf,bnusurf)
      fsurf = int_fsurf(1)*pi







      ! compute total irradiation and add uv
      if (irrad) then

        fstar    = int_freq(ndepth+1,nfreqkap,nuf,fstar_nu)

        if (include_uv) then

            fuv = 0.
            call uvtrans_pt
            fstar = fstar + fuv

         endif

      endif

      ! normalisation for the heating rate
      do i=1,ndepth+1
         qnu(i,:) = 4.*pi*epsfnu(i,:)*bnu(i,:)*kapfnu(i,:)*rhof(i)
      enddo

      qbrad = int_freq(ndepth+1,nfreqkap,nuf,qnu)

      !compute radiative flux and heating rate
      do i=1,ndepth+1
         qnu(i,:) = 4.*pi*epsfnu(i,:)*(jnu(i,:)-snu(i,:))*kapfnu(i,:)*rhof(i)
         fnu(i,:) = 4.*pi*hnu(i,:)
      enddo

      !integrate the flux and heating rate
      frad = int_freq(ndepth+1,nfreqkap,nuf,fnu)
      qrad = int_freq(ndepth+1,nfreqkap,nuf,qnu)

      !fill the right hand side
      err_energy = 0.

      if (irrad) then
         rhs(2*(ndepth+1)) = -(frad(ndepth+1)+fstar(ndepth+1)+fconv(ndepth+1) - fsurf)
         err_energy = max(err_energy,abs(rhs(2*(ndepth+1))/frad(ndepth+1)))

         if (print_err) then
         write(*,*) 'idepth    pres.     temp.     err.        fstar       frad        fconv       fint'
         if (mype==cputerm) write(*,2000) ndepth+1,ppf(ndepth+1)*unitP/1E6,ttf(ndepth+1),abs(rhs(2*(ndepth+1))/frad(ndepth+1)),fstar(ndepth+1),frad(ndepth+1),fconv(ndepth+1),fsurf
      endif


      else
         rhs(2*(ndepth+1)) = -(frad(ndepth+1)+fconv(ndepth+1) - fsurf)
         err_energy = max(err_energy,abs(rhs(2*(ndepth+1))/fsurf))
         err_energy_save(ndepth+1) = abs(rhs(2*(ndepth+1))/fsurf)

         if (print_err) then
         write(*,*) 'idepth    pres.     temp.     err.        fstar       frad        fconv       fint'
         if (mype==cputerm) write(*,2000) ndepth+1,ppf(ndepth+1)*unitP/1E6,ttf(ndepth+1),abs(rhs(2*(ndepth+1))/fsurf),fstar(ndepth+1),frad(ndepth+1),fconv(ndepth+1),fsurf
      endif

      endif




2000 format(3x,i3,2x,es10.1,2x,f6.0,5(2x,es10.1))

      do i = ndepth,1,-1


         if (Rp /=0.) then

            fint = fsurf*Rp**2/rrf(i)**2

         else

            fint = fsurf

         endif

         if (taufstd(i)>tau_heateq) then

            if (irrad) then

               rhs(ndepth+1+i) = -(frad(i)+fstar(i)+fconv(i)-fint)
               err_energy = max(err_energy,abs(rhs(ndepth+1+i)/frad(i)))

                if (print_err) then
               if (mype==cputerm) write(*,2000) i,ppf(i)*unitP/1E6, ttf(i),abs(rhs(ndepth+1+i)/frad(i)),fstar(i),frad(i),fconv(i),fint
            endif

            else

               rhs(ndepth+1+i) = -(frad(i)+fconv(i)-fint)
               err_energy = max(err_energy,abs(rhs(ndepth+1+i)/fint))
               err_energy_save(i) = abs(rhs(ndepth+1+i)/fint)

                if (print_err) then
               if (mype==cputerm) write(*,2000) i,ppf(i)*unitP/1E6, ttf(i),abs(rhs(ndepth+1+i)/fint),fstar(i),frad(i),fconv(i),fint
            endif

            endif


         else

            rhs(ndepth+1+i) = -qrad(i)
            err_energy = max(err_energy,abs(rhs(ndepth+1+i))/qbrad(i))
       !     if (mype==0) write(*,*) i,taufstd(i),abs(rhs(ndepth+1+i))/qbrad(i)
         endif
      enddo

      sjac = epsfjac*bjac

      !compute radiative transfer with one ray for eddington factors
      call radtrans(1,.false.,nfreqjac,kapjac,kapstd,epsfjac,bjac,bjac,hpjac,hdjac,jpjac,jdjac,bjac)

      hjac = hpjac+hdjac
      jjac = jpjac+jdjac

         !keep H2O opacity constant for the jacobian to help convergence
         !also Fe opacity when converging with rainout chemistry
      if(cst_h2o) then
         afh2o_save(:) = amolf(:,imol('H2O'))
         afnh3_save(:) = amolf(:,imol('NH3'))

         ah2o_save(:) = amol(:,imol('H2O'))
         anh3_save(:) = amol(:,imol('NH3'))

         amolf(:,imol('H2O')) = amolf(ndepth+1,imol('H2O'))
         amol(:,imol('H2O'))  = amolf(ndepth+1,imol('H2O'))

         amolf(:,imol('NH3')) = amolf(ndepth+1,imol('NH3'))
         amol(:,imol('NH3'))  = amolf(ndepth+1,imol('NH3'))

         if(debug>0 .and. mype==cputerm) write(*,*) 'start get_opacity'
         if (corr_k) then
            call get_corr_k_opacity
         else
            call get_lbl_opacity
         endif

         if (hminus) then
            call get_hm_opacity
         endif

         if (scatter) then
            call get_scat_opacity
         endif
         if(debug>0 .and. mype==cputerm) write(*,*) 'end get_opacity'

         if (nbin == 1) then

            dbnu   = dbnu_tot
            kapjac = kapnu
            nujac  = nuf
            epsfjac = epsfnu
            bjac   = bnu
            dkapfjac = dkapfnu
            kapfjac = kapfnu

         endif

         amolf(:,imol('H2O')) = afh2o_save(:)
         amolf(:,imol('NH3')) = afnh3_save(:)

         amol(:,imol('H2O'))  = ah2o_save(:)
         amol(:,imol('NH3'))  = anh3_save(:)

      endif

      do i =1,ndepth+1

         if (include_dkapdpt) then
            !reset the perturbation
            dsnu = 0. ; dhnu = 0. ; djnu = 0.

            !compute the source pressure perturbation with jnu_1ray
            where (kapfjac(i,:) /= 0.)
               dsnu(i,:) = (bjac(i,:) - jjac(i,:)) &
                    & *(dkapfjac(i,:)/kapfjac(i,:)-dkapfstd(i)/kapfstd(i))
            elsewhere

               dsnu(i,:) = 0.

            endwhere

            !compute the pressure perturbation
            call radtrans(1,.false.,nfreqjac,kapjac,kapstd,epsfjac,dsnu,dsnu,dhnup,dhnud,djnup,djnud,dsnu)

            dhnu = dhnup+dhnud
            djnu = djnup+djnud

            ! !correct with eddington factors
            ! do ii = 1,nfreqjac
            !    where (hjac(:,ii) /=0.) dhnu(:,ii) = dhnu(:,ii)/hjac(:,ii)*hnu(:,ii*nbin)
            !    where (jjac(:,ii) /=0.) djnu(:,ii) = djnu(:,ii)/jjac(:,ii)*jnu(:,ii*nbin)
            ! enddo

            do i2 =1,ndepth+1

               djnu(i2,:) = epsfjac(i2,:)*(djnu(i2,:)-dsnu(i2,:))*4.*pi*kapfjac(i2,:)*rhof(i2)
               dhnu(i2,:) = 4.*pi*dhnu(i2,:)

            enddo


            !integrate the perturbation
            dfrad =  int_freq(ndepth+1,nfreqjac,nujac,dhnu)
            dqrad =  int_freq(ndepth+1,nfreqjac,nujac,djnu)

            !fill the jacobi matrix
            do i2 = 1,ndepth

               if (taufstd(i2)>tau_heateq) then
                  jacobi(ndepth+1+i2,i) = dfrad(i2)+dfconv(i2,i)
               else
                  jacobi(ndepth+1+i2,i) = dqrad(i2)
               endif

            enddo

            jacobi(2*(ndepth+1),i) = dfrad(ndepth+1)+dfconv(ndepth+1,i)

         else

            !fill the jacobi matrix
            do i2 = 1,ndepth

               if (taufstd(i2)>tau_heateq) then
                  jacobi(ndepth+1+i2,i) = dfconv(i2,i)
               else
                  jacobi(ndepth+1+i2,i) = 0.
               endif


            enddo

            jacobi(2*(ndepth+1),i) = dfconv(ndepth+1,i)


         endif

         !reset the perturbation
         dsnu = 0. ; dhnu = 0. ; djnu = 0.

         !compute the source temperature perturbation with jnu_1ray
         where(kapfjac(i,:) /= 0.)

            dsnu(i,:) = epsfnu(i,:)*dbnu(i,:) +(bjac(i,:) - jjac(i,:)) &
                 & *(dkapfjac(ndepth+1+i,:)/kapfjac(i,:)-dkapfstd(ndepth+1+i)/kapfstd(i))
         elsewhere

            dsnu(i,:) = epsfnu(i,:)*dbnu(i,:)

         endwhere



         !compute the temperature perturbation
         call radtrans(1,.false.,nfreqjac,kapjac,kapstd,epsfjac,dsnu,dsnu,dhnup,dhnud,djnup,djnud,dsnu)

         dhnu = dhnup+dhnud
         djnu = djnup+djnud

         ! !correct with eddington factors
         ! do ii = 1,nfreqjac
         !    where (hjac(:,ii) /=0.) dhnu(:,ii) = dhnu(:,ii)/hjac(:,ii)*hnu(:,ii*nbin)
         !    where (jjac(:,ii) /=0.) djnu(:,ii) = djnu(:,ii)/jjac(:,ii)*jnu(:,ii*nbin)
         ! enddo

         do i2 =1,ndepth+1

            djnu(i2,:) = epsfjac(i2,:)*(djnu(i2,:)-dsnu(i2,:))*4.*pi*kapfjac(i2,:)*rhof(i2)
            dhnu(i2,:) = 4.*pi*dhnu(i2,:)

         enddo


         !integrate the perturbation
         dfrad =  int_freq(ndepth+1,nfreqjac,nujac,dhnu)
         dqrad =  int_freq(ndepth+1,nfreqjac,nujac,djnu)

         !ib = 24
         !if (i==1) write(*,*) sum(kapfnu(:,(ib-1)*nkmix+1:ib*nkmix)*weightsfnu(:,(ib-1)*nkmix+1:ib*nkmix),2)

         ! if(i==15) write(*,*) dfrad
         ! write(*,*) dfrad(15)


         !fill the jacobi matrix
         do i2 = 1,ndepth

            if (taufstd(i2)>tau_heateq) then
                jacobi(ndepth+1+i2,ndepth+1+i) = dfrad(i2)+dfconv(i2,ndepth+1+i)
             else
                jacobi(ndepth+1+i2,ndepth+1+i) = dqrad(i2)
             endif


         enddo

         jacobi(2*(ndepth+1),ndepth+1+i) = dfrad(ndepth+1)+dfconv(ndepth+1,ndepth+1+i)

    !     if (i<ndepth+1) then
    !     jacobi(ndepth+1+i,ndepth+1+i) = jacobi(ndepth+1+i,ndepth+1+i)+fstar(i)*(kapfross(i))/(ttf(i))*(taufstd(i+1)-taufstd(i))/kapstd(i)*10000.
    !     endif
      enddo



      !deallocate
      deallocate(hnu,jnu)
      deallocate(dhnu,djnu)
      deallocate(dsnu,bnu)
      deallocate(dbnu)
      deallocate(dfrad)

    end subroutine perturb_energy




    subroutine compute_spectrum(chem_loc,ilong)

      use mod_grid
      use mod_opacity
      use mod_param
      use mod_cst
      use mod_eos
      use mod_util
      use mod_chem
      use mod_radtrans

      implicit none


      real,dimension(:,:),allocatable :: jnu,hnu,bnu,snu,jnup,jnud,hnup,hnud
      real,dimension(:),allocatable :: rnu
      integer :: idepth,i,ib,ilong
      character(3) :: chem_loc

      if (any(chem_loc==(/'eq ','ana','man'/))) then

         if (debug>0 .and. mype==cputerm) write(*,*) 'start get_eq_chemistry'
         call get_eq_chemistry
         if (debug>0 .and. mype==cputerm) write(*,*) 'end get_eq_chemistry'

      else if(chem_loc == 'neq') then

         if (mype==cputerm) then

            write(*,headrfmt) 'Error compute_spectrum'
            write(*,'(A)') 'neq cannot be used with compute_spectrum should use cst'
            write(*,chrpmfmt) 'chem',chem
            write(*,*)
            write(*,separfmt)
            stop

         endif

      endif

      do idepth=1,ndepth+1

         do i=1,nmol

            if (Acst(i)>0.) then
               Amolf(idepth,i) = Acst(i)
            endif

            if (ttf(idepth)<Tcrit(i) .or. ppf(idepth)>Pcrit(i)) then
               Amolf(idepth,i) = 0.
            endif

         enddo

      enddo

      !interpolate the Amol at cell center, using geometric mean
      do idepth = 1,ndepth

         Amol(idepth,:) = sqrt(Amolf(idepth,:)*Amolf(idepth+1,:))

      enddo

      nn  = get_nn(ndepth,pp,tt)
      nnf = get_nn(ndepth+1,ppf,ttf)

      rho = get_rho(ndepth,pp,tt,amol)
      rhof = get_rho(ndepth+1,ppf,ttf,amolf)

      ! Compute vertical scale and gravity
      if (um_input) then
        call get_um_rrgr
      else
        call get_rrgr
      end if

      if(phase_curve .and. cloud_night_side) then

         if(cos(longitude_um(ilong))>0.) then

            cloud = .true.

         endif

      endif

      if(debug>0 .and. mype==cputerm) write(*,*) 'start get_opacity'
      if (corr_k) then
         call get_corr_k_opacity
      else
         call get_lbl_opacity
      endif

      cloud = .false.

      if (hminus) then
         call get_hm_opacity
      endif

      if (scatter) then
         call get_scat_opacity
      endif
      if(debug>0 .and. mype==cputerm) write(*,*) 'end get_opacity'

      ! assume hydrostatic equilibrium
      if (um_input) then
         do idepth=2,ndepth+1
            taufstd(idepth) = taufstd(idepth-1)+rho(idepth-1)*kapstd(idepth-1)*(rrf(idepth-1)-rrf(idepth))
         enddo

      else if (hydrostatic) then

         do idepth =2,ndepth+1

            taufstd(idepth) = taufstd(idepth-1)+(ppf(idepth)-ppf(idepth-1))*kapstd(idepth-1)/gr(idepth-1)

         enddo
      endif


      if (surface_spectrum .or. (transmission_spectrum .and. trans_spec_emission) .or. phase_curve) then

         !allocate moments of intensity and source function
         allocate(hnu(ndepth+1,nfreqkap),jnu(ndepth+1,nfreqkap),bnu(ndepth+1,nfreqkap),snu(ndepth+1,nfreqkap))
         allocate(hnup(ndepth+1,nfreqkap),hnud(ndepth+1,nfreqkap),jnup(ndepth+1,nfreqkap),jnud(ndepth+1,nfreqkap))

         fstar_nu = 0.

         if (irrad) then
            !compute star flux
            call startrans(nfreqkap,kapnu,kapstd,fstar_nu,Irad)
            fstar = int_freq(ndepth+1,nfreqkap,nuf,fstar_nu)

         endif

         bnu      = get_bnu (ndepth+1,nfreqkap,ttf)

         !source term for RT
         if (scatter) then
            if (irrad) then
               snu = epsfnu*bnu+(1.-epsfnu)*abs(fstar_nu/murad)/(4.*pi)
            else
               snu = epsfnu*bnu
            endif

         else

            !compute planck function
            snu      = bnu

         endif

         !compute radiative flux at interfaces
         if(debug>0 .and. mype==cputerm) write(*,*) 'start radtrans'
         call radtrans(nrays,scatter,nfreqkap,kapnu,kapstd,epsfnu,snu,bnu,hnup,hnud,jnup,jnud,snu_scat)
         if(debug>0 .and. mype==cputerm) write(*,*) 'end radtrans'

         hnu = hnup+hnud
         jnu = jnup+jnud

         !compute radiative flux and heating rate
         do idepth=1,ndepth+1
            qnu(idepth,:) = 4.*pi*epsfnu(idepth,:)*(jnu(idepth,:)-snu(idepth,:))*kapfnu(idepth,:)*rhof(idepth)
            fnu(idepth,:) = 4.*pi*hnu(idepth,:)
         enddo

         !integrate the flux and heating rate
         frad = int_freq(ndepth+1,nfreqkap,nuf,fnu)
         qrad = int_freq(ndepth+1,nfreqkap,nuf,qnu)



      ! do idepth = 1,ndepth

      !    if (mype==0) write(*,*) ppf(idepth)*unitP,(frad(idepth)+fstar(idepth))*unitIF,(frad(idepth+1)-frad(idepth)+fstar(idepth+1)-fstar(idepth))/(rrf(idepth+1)-rrf(idepth))*unitIF/unitL

      ! enddo
      !    stop

         ! if ((surface_spectrum .or. transmission_spectrum .or. phase_curve) .and. corr_k .and. mype==cputerm) then


         !    write(*,headrfmt) 'Error compute_spectrum'
         !    write(*,'(A)') 'compute_spectrum not implemented for correlated k'
         !    write(*,logpmfmt) 'corr_k',corr_k
         !    write(*,*)
         !    write(*,separfmt)
         !    stop

         ! endif

      endif

       if (surface_spectrum .and. .not.surface_spectrum_3d .and. &
           .not.phase_curve) then
         call compute_surf_spec
       endif

    ! Compute Contribution Function when this true

        if (calc_cf) then
           call compute_con_fun
        end if

    end subroutine compute_spectrum


    !> Wrapper routine for computing 3D spectra
    !<
    subroutine compute_spectrum_um(chem_loc)

      use mod_grid
      use mod_param
      use mod_cst
      use mod_radtrans
      implicit none

      character(3) :: chem_loc
      integer :: ilong, ilat
      integer :: ilong_1d(1), ilat_1d(1)

      ! Initialise 3D spectrum calculation using entire UM output.
      if (transmission_spectrum_3d .or. surface_spectrum_3d .or. &
          phase_curve) then

        ! When temperature field is read in, fluxes are computed
        do ilong=1,nlong
          do ilat=1,nlat
            if (debug > 0 .and. mype == cputerm) then
              write(*,*) 'Processing longitude', ilong, ' of ', nlong
              write(*,*) 'Processing latitude ', ilat, ' of ', nlat
            end if

          ! Copy 3D fields over to 1D array for use by compute_spectrum
          call copy_um_pt(ilong,ilat)
          ! Set correct upper boundary condition for irradiation
          call set_murad(ilong,ilat)

          ! Compute spectrum for current longitude and latitude
          call compute_spectrum(chem_loc,ilong)

          ! Need to save fields for calculating spectra
          call save_um_grid_vars(ilong,ilat)
          call save_um_radtrans_vars(ilong,ilat)

        end do
      end do

      ! Initialise 1D spectrum calculation using UM output at chosen longitude
      ! and latitude point.
      else

        ! Find closest longitude and latitude point on grid
        ilong_1d = minloc(abs(longitude_um - long_spec))
        ilat_1d = minloc(abs(latitude_um - lat_spec))
        if (mype == cputerm) then
          write(*,*) 'Using UM output at:'
          write(*,'(a,f5.1,a)') ' Longitude: ', longitude_um(ilong_1d)*180./pi
          write(*,'(a,f5.1,a)') ' Latitude: ', (pi/2. - latitude_um(ilat_1d))*180./pi
        endif

        ! Copy pressure and temperature from 3D to 1D arrays
        call copy_um_pt(ilong_1d(1), ilat_1d(1))

        ! Set correct upper boundary condition for irradiation
        call set_murad(ilong_1d(1), ilat_1d(1))

        ! Compute spectrum for this longitude and latitude
        call compute_spectrum(chem_loc,ilong_1d(1))

      endif

      ! Compute 3D surface spectrum from UM output
      if (surface_spectrum .and. surface_spectrum_3d .and. &
          .not.phase_curve) then
        if (debug>0 .and. mype==cputerm) write(*,*) 'start compute_surf_spec'
        call compute_surf_spec
        if (debug>0 .and. mype==cputerm) write(*,*) 'end compute_surf_spec'
      endif


    end subroutine compute_spectrum_um


    !> Calculates surface emission and writes to file - use in LbL, ATMO corrk and UM corrk
    !<
    subroutine compute_surf_spec

      use mod_grid
      use mod_opacity
      use mod_param
      use mod_cst
      use mod_eos
      use mod_util
      use mod_chem
      use mod_radtrans

      implicit none
      include 'netcdf.inc'

#if PARA==1
      include 'mpif.h'
      integer :: info,tag=0,stat
#endif

      integer :: id_file,id_var,id_nfreq,iwrite,istart,icount,i,idepth,iband,id_nfreq_sw,id_nfreq_lw, icount_sw,icount_lw
      integer :: nband_um
      real,dimension(:),allocatable    :: nu_sw,nu_lw,fnu_sw_int_um,fnu_lw_int_um
      ! Surface spectrum
      real,dimension(:),allocatable :: fnu_surf_spec,nu_surf

      ! Allocate flux and wavenumber arrays
      if (.not. corr_k) then
         !Line-by-line
         allocate(fnu_surf_spec(nfreqkap))
         allocate(nu_surf(nfreqkap))
      else
         !Correlated-k
         allocate(fnu_surf_spec(nband))
         allocate(nu_surf(nband))
      endif

      ! Calculate surface spectrum if necessary
      if (surface_spectrum_3d) then
        ! Combine ATMO fluxes
        call calc_surf_spec_3d(long_obs,lat_obs,nfreqkap,fnu(1,:))
      endif

      ! Set TOA flux and wavenumber
      if (.not. corr_k) then
         ! Line-by-line
         nu_surf = nu(:)
         fnu_surf_spec = fnu(1,:)
      else
        ! Correlated-k
        do iband = 1,nband
          nu_surf(iband) = 0.5*(nubandmin(iband)+nubandmax(iband))
          fnu_surf_spec(iband) = sum(fnu(1,(iband-1)*nkmix+1:iband*nkmix)*weightsfnu(1,(iband-1)*nkmix+1:iband*nkmix)) &
          /(nubandmax(iband)-nubandmin(iband))
        enddo
      endif

      ! Set start and end counts for netCDF
      !Line-by-line
      if (.not. corr_k) then
         istart = inukap_start
         icount = nfreqkap
      ! Correlated-k
      else
         istart = 1
         icount = nband
      endif

      ! Write spectrum to file
      if (mype==cputerm) then

        ! Create file
        call nf(nf_create(fspectrum,nf_clobber,id_file))

        ! Define nfreq dimension
        ! Line-by-line
        if (.not. corr_k) then
          call nf(nf_def_dim(id_file,'nfreq',nfreq,id_nfreq))
        ! Correlated-k
        else
          call nf(nf_def_dim(id_file,'nfreq',nband,id_nfreq))
        endif

        ! Define wavenumber variable
        call nf(nf_def_var(id_file,'nu',nf_double,1,id_nfreq,id_var))
        call nf(nf_put_att_text(id_file,id_var,'units',4,'cm-1'))

        ! Define min and max wavenumber of bands and TOA flux
        if (corr_k) then
          ! ATMO fluxes
            call nf(nf_def_var(id_file,'nubandmin',nf_double,1,id_nfreq,id_var))
            call nf(nf_put_att_text(id_file,id_var,'units',4,'cm-1'))

            call nf(nf_def_var(id_file,'nubandmax',nf_double,1,id_nfreq,id_var))
            call nf(nf_put_att_text(id_file,id_var,'units',4,'cm-1'))

            call nf(nf_def_var(id_file,'fnu',nf_double,1,id_nfreq,id_var))
            call nf(nf_put_att_text(id_file,id_var,'units',15,'erg s-1 cm-2 cm'))

            ! Planet radius
            call nf(nf_def_var(id_file,'R_planet_TOA',nf_double,0,0,id_var))
            call nf(nf_put_att_text(id_file,id_var,'units',2,'cm'))

        else

          call nf(nf_def_var(id_file,'fnu',nf_double,1,id_nfreq,id_var))
          call nf(nf_put_att_text(id_file,id_var,'units',15,'erg s-1 cm-2 cm'))

          ! Planet radius
          call nf(nf_def_var(id_file,'R_planet_TOA',nf_double,0,0,id_var))
          call nf(nf_put_att_text(id_file,id_var,'units',2,'cm'))

        end if

        ! Stellar flux
        if (irrad .or. um_input) then

          call nf(nf_def_var(id_file,'fnu_star',nf_double,1,id_nfreq,id_var))
          call nf(nf_put_att_text(id_file,id_var,'units',15,'erg s-1 cm-2 cm'))

          call nf(nf_def_var(id_file,'R_star',nf_double,0,0,id_var))
          call nf(nf_put_att_text(id_file,id_var,'units',2,'cm'))
        endif

        ! Close file
        call nf(nf_close(id_file))

       endif

#if PARA == 1
      call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif

!      call nf(nf_open_par(fspectrum,nf_write,MPI_COMM_WORLD,MPI_INFO_NULL,id_file))
      iwrite = 0
      if(mype==0) iwrite = 1

#if PARA==1
      if (mype>0) call MPI_RECV(iwrite,1,MPI_INTEGER,mype-1,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,info)
#endif

      ! Write variables to file
      if (iwrite==1) then

        ! open file
        call nf(nf_open(fspectrum,nf_write,id_file))


          call nf(nf_inq_varid(id_file,'nu' ,id_var))
          call nf(nf_put_vara_double(id_file,id_var,istart,icount,nu_surf*(1/unitL)))

          if (corr_k) then
            call nf(nf_inq_varid(id_file,'nubandmin' ,id_var))
            call nf(nf_put_vara_double(id_file,id_var,istart,icount, &
                 nubandmin(1:nband)*(1/unitL)))

            call nf(nf_inq_varid(id_file,'nubandmax' ,id_var))
            call nf(nf_put_vara_double(id_file,id_var,istart,icount, &
                 nubandmax(1:nband)*(1/unitL)))
          end if

          call nf(nf_inq_varid(id_file,'fnu' ,id_var))
          call nf(nf_put_vara_double(id_file,id_var,istart,icount,fnu_surf_spec*unitSF))

          ! Planet radius
          call nf(nf_inq_varid(id_file,'R_planet_TOA',id_var))
          call nf(nf_put_vara_double(id_file,id_var,1,1,rrf(1)*unitL))

        ! Stellar flux
        if (irrad .or. um_input) then
          call nf(nf_inq_varid(id_file,'fnu_star' ,id_var))
            if (.not. corr_k) then
              call nf(nf_put_vara_double(id_file,id_var,istart,icount, &
                   (4*pi*irad)*unitSF))
            else
              call nf(nf_put_vara_double(id_file,id_var,istart,icount, &
                   (4*pi*irad(1:nfreqkap:nkmix)/ &
                   (nubandmax(1:nband)-nubandmin(1:nband)))*unitSF))
            endif
          call nf(nf_inq_varid(id_file,'R_star',id_var))
          call nf(nf_put_vara_double(id_file,id_var,1,1,rstar*unitL))
        endif

        ! Close file
        call nf(nf_close(id_file))

      endif

#if PARA==1
      if (mype<ncpu-1) call MPI_SEND(iwrite,1,MPI_INTEGER,mype+1,tag,MPI_COMM_WORLD,info)
#endif

    end subroutine compute_surf_spec


    !> Wrapper routine for calculation of phase curves. Assumes planet is
    !! tidally locked.
    !<
    subroutine compute_phase_curve

      use mod_grid
      use mod_opacity
      use mod_param
      use mod_cst
      use mod_eos
      use mod_util
      use mod_chem
      use mod_radtrans

      implicit none
      include 'netcdf.inc'

#if PARA==1
      include 'mpif.h'
      integer :: info,tag=0,stat
#endif

      integer :: id_file,id_var,id_nfreq,id_nphase,iwrite,istart,icount,i,id_nfreq_sw,id_nfreq_lw, icount_sw,icount_lw
      ! Surface spectrum
      real,dimension(:),allocatable :: fnu_surf_spec,nu_surf
      real,dimension(:),allocatable    :: nu_sw,nu_lw,fnu_sw_int_um,fnu_lw_int_um
      ! Array of phase angles
      real,dimension(:),allocatable :: phase_angle
      ! Longitude and latitude of observer
      real :: long_obs_loc,lat_obs_loc
      ! Resolution of phase curve
      integer :: n_points
      ! Loop index
      integer :: iband

      n_points = np_phase_curve

      allocate(phase_angle(n_points))

      if (.not. corr_k) then
         allocate(fnu_surf_spec(nfreqkap))
         allocate(nu_surf(nfreqkap))
      else
         allocate(fnu_surf_spec(nband))
         allocate(nu_surf(nband))
      endif

      if (.not. corr_k) then
        istart = inukap_start
        icount = nfreqkap
      else
        istart = 1
        icount = nband
       endif

      if (.not. corr_k) then
        nu_surf = nu(:)
      else
        do iband = 1,nband
          nu_surf(iband) = 0.5*(nubandmin(iband)+nubandmax(iband))
        enddo
      endif

      ! Construct array of phase angles
      call linspace(0.,2*pi,n_points,phase_angle)

      if (mype==cputerm) then

        call nf(nf_create(fspectrum,nf_clobber,id_file))
        if (.not. corr_k) then
          call nf(nf_def_dim(id_file,'nfreq',nfreq,id_nfreq))
        else
          call nf(nf_def_dim(id_file,'nfreq',nband,id_nfreq))
        endif

        call nf(nf_def_dim(id_file,'phase_angle',n_points,id_nphase))

        ! Define wavenumber variable
        call nf(nf_def_var(id_file,'nu',nf_double,1,id_nfreq,id_var))
        call nf(nf_put_att_text(id_file,id_var,'units',4,'cm-1'))

        if (corr_k) then
            call nf(nf_def_var(id_file,'nubandmin',nf_double,1,id_nfreq,id_var))
            call nf(nf_put_att_text(id_file,id_var,'units',4,'cm-1'))
            call nf(nf_def_var(id_file,'nubandmax',nf_double,1,id_nfreq,id_var))
            call nf(nf_put_att_text(id_file,id_var,'units',4,'cm-1'))
        end if

        call nf(nf_def_var(id_file,'phase_angle',nf_double,1,id_nphase,id_var))
        call nf(nf_put_att_text(id_file,id_var,'units',7,'degrees'))

          call nf(nf_def_var(id_file,'fnu',nf_double,2, &
               (/id_nphase,id_nfreq/),id_var))
          call nf(nf_put_att_text(id_file,id_var,'units',2,'erg s-1 cm-2 cm'))
          call nf(nf_def_var(id_file,'fnu_star',nf_double,1,id_nfreq,id_var))
          call nf(nf_put_att_text(id_file,id_var,'units',15,'erg s-1 cm-2 cm'))




        call nf(nf_def_var(id_file,'R_star',nf_double,0,0,id_var))
        call nf(nf_put_att_text(id_file,id_var,'units',2,'cm'))
        call nf(nf_def_var(id_file,'R_planet_TOA',nf_double,0,0,id_var))
        call nf(nf_put_att_text(id_file,id_var,'units',2,'cm'))

        call nf(nf_close(id_file))

        call nf(nf_open(fspectrum,nf_write,id_file))

          call nf(nf_inq_varid(id_file,'nu',id_var))
          call nf(nf_put_vara_double(id_file,id_var,istart,icount,nu_surf*(1/unitL)))

        !Define wavenumber limits of bands
        if (corr_k) then
            call nf(nf_inq_varid(id_file,'nubandmin',id_var))
            call nf(nf_put_vara_double(id_file,id_var,istart,icount,nubandmin*(1/unitL)))
            call nf(nf_inq_varid(id_file,'nubandmax',id_var))
            call nf(nf_put_vara_double(id_file,id_var,istart,icount,nubandmax*(1/unitL)))
            call nf(nf_inq_varid(id_file,'R_planet_TOA',id_var))
            call nf(nf_put_vara_double(id_file,id_var,1,1,rrf(1)*unitL))
        end if

        call nf(nf_inq_varid(id_file,'phase_angle',id_var))
        call nf(nf_put_vara_double(id_file,id_var,1,n_points,phase_angle*180./pi))

        call nf(nf_inq_varid(id_file,'fnu_star' ,id_var))
        if (.not. corr_k) then
          call nf(nf_put_vara_double(id_file,id_var,istart,icount, &
              (4*pi*irad)*unitSF))
        else
          call nf(nf_put_vara_double(id_file,id_var,istart,icount, &
              (4*pi*irad(1:nfreqkap:nkmix)/ &
              (nubandmax(1:nband)-nubandmin(1:nband)))*unitSF))
        endif

        call nf(nf_inq_varid(id_file,'R_star',id_var))
        call nf(nf_put_vara_double(id_file,id_var,1,1,rstar*unitL))

        call nf(nf_close(id_file))
      endif

      ! Loop over all phase angles and calculate surface spectrum
      do i=1,n_points

        ! Set direction of observer. Assume planet is tidally locked.
        long_obs_loc = 2*pi-phase_angle(i)
        lat_obs_loc = lat_obs!pi/2.

        if (debug>0 .and. mype==cputerm) write(*,*) 'Computing phase curve fluxes'
        if (debug>0 .and. mype==cputerm) write(*,*) 'phase_curve ',i,' of ',n_points

        if (.not. corr_k) then
          call calc_surf_spec_3d(long_obs_loc,lat_obs_loc,nfreqkap,fnu(1,:))
          fnu_surf_spec = fnu(1,:)
        else
          call calc_surf_spec_3d(long_obs_loc,lat_obs_loc,nfreqkap,fnu(1,:))
          do iband = 1,nband
            fnu_surf_spec(iband) = sum(fnu(1,(iband-1)*nkmix+1:iband*nkmix)*weightsfnu(1,(iband-1)*nkmix+1:iband*nkmix))/(nubandmax(iband)-nubandmin(iband))
          enddo
        endif

        ! Write current emission spectrum to file
#if PARA == 1
        call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif

!         call nf(nf_open_par(fspectrum,nf_write,MPI_COMM_WORLD,MPI_INFO_NULL,id_file))
        iwrite = 0
        if(mype==0) iwrite = 1

#if PARA==1
        if (mype>0) call MPI_RECV(iwrite,1,MPI_INTEGER,mype-1,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,info)
#endif

        !Write flux at current phase point to file
        if (iwrite==1) then

          call nf(nf_open(fspectrum,nf_write,id_file))

            call nf(nf_inq_varid(id_file,'fnu',id_var))
            call nf(nf_put_vara_double(id_file,id_var,(/i,istart/),(/1,icount/), &
                fnu_surf_spec*unitSF))

          call nf(nf_close(id_file))

        endif

#if PARA==1
        if (mype<ncpu-1) call MPI_SEND(iwrite,1,MPI_INTEGER,mype+1,tag,MPI_COMM_WORLD,info)
#endif
      end do

    end subroutine compute_phase_curve


    !> Wrapper routine for calculation of transmission spectrum
    !<
    subroutine compute_trans_spec

      use mod_grid
      use mod_opacity
      use mod_param
      use mod_cst
      use mod_eos
      use mod_util
      use mod_chem
      use mod_radtrans

      implicit none
      include 'netcdf.inc'

#if PARA==1
      include 'mpif.h'
      integer :: info,tag=0,stat
#endif


      integer :: id_file,id_var,id_dim,iwrite,istart,icount,iband
      ! Transit radius
      real,dimension(:),allocatable :: Rp_transit, Rp_transit_corrk,nu_corrk

      allocate(Rp_transit(nfreqkap))
      if (corr_k) then
      	allocate(Rp_transit_corrk(nband))
      	allocate(nu_corrk(nband))
      endif

      ! Calculate transit radius
      if (transmission_spectrum_3d) then
        call calc_trans_radius_3d(nfreqkap,Rp_transit)
      else
        call calc_trans_radius_1d(nfreqkap,Rp_transit)
      end if

      ! Set writing thread
      iwrite = 0
      if(mype==0) iwrite = 1

      if (iwrite==1) then
        call nf(nf_create(ftrans_spec,nf_clobber,id_file))
        if (corr_k) then
        	call nf(nf_def_dim(id_file,'nu',nband,id_dim))
        else
        	call nf(nf_def_dim(id_file,'nu',nfreq,id_dim))
        endif

        call nf(nf_def_var(id_file,'nu',nf_double,1,id_dim,id_var))
        call nf(nf_put_att_text(id_file,id_var,'units',4,'cm-1'))

        call nf(nf_def_var(id_file,'transit_radius',nf_double,1,id_dim,id_var))
        call nf(nf_put_att_text(id_file,id_var,'units',6,'R_star'))

        call nf(nf_def_var(id_file,'R_star',nf_double,0,0,id_var))
        call nf(nf_put_att_text(id_file,id_var,'units',2,'cm'))

        call nf(nf_close(id_file))
      endif

#if PARA == 1

      call MPI_BARRIER(MPI_COMM_WORLD,info)

#endif

      !call nf(nf_open_par(ftrans_spec,nf_write,MPI_COMM_WORLD,MPI_INFO_NULL,id_file))

#if PARA==1
      if (mype>0) call MPI_RECV(iwrite,1,MPI_INTEGER,mype-1,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,info)
#endif

      ! Set where in file to write
      if (corr_k) then
         istart = 1
         icount = nband
      else
         istart = inukap_start
         icount = nfreqkap
      endif

      if (corr_K) then
         Rp_transit = Rp_transit**2-rrf(ndepth+1)**2

         do iband = 1,nband
            nu_corrk(iband) = 0.5*(nubandmin(iband)+nubandmax(iband))
            Rp_transit_corrk(iband) = sum(Rp_transit((iband-1)*nkmix+1:iband*nkmix)*weightsfnu(1,(iband-1)*nkmix+1:iband*nkmix))
         enddo

         Rp_transit = sqrt(rrf(ndepth+1)**2 + Rp_transit)
         Rp_transit_corrk = sqrt(rrf(ndepth+1)**2 + Rp_transit_corrk)
      endif

      ! Write transmission spectrum to file
      if (iwrite==1) then
        call nf(nf_open(ftrans_spec,nf_write,id_file))

        call nf(nf_inq_varid(id_file,'nu',id_var))
!      !  call nf(nf_var_par_access(id_file,id_var,NF_COLLECTIVE))
		if (corr_k) then
			call nf(nf_put_vara_double(id_file,id_var,istart,icount,nu_corrk/unitL))
		else
        	call nf(nf_put_vara_double(id_file,id_var,istart,icount,nu/unitL))
        endif

        call nf(nf_inq_varid(id_file,'transit_radius',id_var))
!       ! call nf(nf_var_par_access(id_file,id_var,NF_COLLECTIVE))
		if (corr_k) then
       		call nf(nf_put_vara_double(id_file,id_var,istart,icount,Rp_transit_corrk/rstar))
       	else
       		call nf(nf_put_vara_double(id_file,id_var,istart,icount,Rp_transit/rstar))
        endif
        call nf(nf_inq_varid(id_file,'R_star',id_var))
!       ! call nf(nf_var_par_access(id_file,id_var,NF_COLLECTIVE))
        call nf(nf_put_vara_double(id_file,id_var,1,1,rstar*unitL))

        call nf(nf_close(id_file))
      endif

#if PARA==1
      if (mype<ncpu-1) call MPI_SEND(iwrite,1,MPI_INTEGER,mype+1,tag,MPI_COMM_WORLD,info)
#endif

    end subroutine compute_trans_spec


    !> Find the pressure/temp profile using a newton-raphson method
    !<
    subroutine find_profile(chem_loc)

      use mod_grid
      use mod_param
      use mod_opacity
      use mod_cst
      use mod_chem
      use mod_util
      use mod_eos
      use mod_radtrans

      implicit none

      integer :: i,nit,idepth,max_hydro_level,max_energy_level,max_dp_level,max_dt_level,cnt,nit_cnt
      integer :: varerrsteady
      real,dimension(2*(ndepth+1)) :: dpt
      real,dimension(2*(ndepth+1)) :: dpt_mean
      character(4)  :: fnit
      real :: err_hydro,err_energy,err,errold_hydro,errveryold_hydro,errold_energy,errveryold_energy,step,new_tau
      real :: varerr
      real,dimension(:),allocatable :: errold
      character(3) :: chem_loc

      real,dimension(ndepth+1)::ttf_coupled

      integer :: imin,imax,iker,kertt_coupled

      if(.not. (solve_energy .or. solve_hydro)) then
         return
      endif

      err        = 1.
      varerr     = 1.
      varerrsteady = 0
      allocate(errold(nitmax))
      errold(:) = 1.
      errold_hydro     = 1.
      errveryold_hydro = 1.
      errold_energy = 1.
      errveryold_energy = 1.
      nit  = 0
      cnt = 1

      if (solve_hydro) then
         ppf(1) = psurf
      endif

      !do while (err>accuracy .and. varerrsteady<3 .and. nit<nitmax)
      do while (err>accuracy .and. nit<nitmax)

         jacobi = 0.
         rhs    = 0.

         do i = 1,2*(ndepth+1)
            jacobi(i,i) = 1.
         enddo

         err_hydro  = 0.
         err_energy = 0.

         if (any(chem_loc==(/'eq ','ana','man'/))) then
            if (debug>0 .and. mype==cputerm) write(*,*) 'start get_eq_chemistry'
            call get_eq_chemistry
            if (debug>0 .and. mype==cputerm) write(*,*) 'end get_eq_chemistry'

         else if(chem_loc == 'neq') then

               if (mype==cputerm) then

               write(*,headrfmt) 'Error find_profile'
               write(*,'(A)') 'neq cannot be used inside find_profile'
               write(*,chrpmfmt) 'chem',chem
               write(*,*)
               write(*,separfmt)
               stop

            endif

         endif

         nn  = get_nn(ndepth,pp,tt)
         nnf = get_nn(ndepth+1,ppf,ttf)

         rho = get_rho(ndepth,pp,tt,amol)
         rhof = get_rho(ndepth+1,ppf,ttf,amolf)

         call get_rrgr

         if(debug>0 .and. mype==cputerm) write(*,*) 'start get_opacity'
         if (corr_k) then
            call get_corr_k_opacity
         else
            call get_lbl_opacity
         endif
         if (hminus) then
            call get_hm_opacity
         endif
         if (scatter) then
            call get_scat_opacity
         endif
         if(debug>0 .and. mype==cputerm) write(*,*) 'end get_opacity'

         !if (chem == 'cst') then

         !  open(unit=1,file='kapstd.dat')
         !  read(1,*) kapstd
         !  close(1)
         !  open(unit=1,file='kapfstd.dat')
         !  read(1,*) kapfstd
         !  close(1)

         !endif

         if (.not. include_dkapdpt) then

            dkapstd = 0.
            dkapfstd = 0.
            dkapnu = 0.
            dkapfnu = 0.

         endif

         if(solve_hydro) then
            if (debug>0 .and. mype==cputerm) write(*,*) 'start perturb_hydro'
            call perturb_hydro (err_hydro, err_hydro_prf)
            if (debug>0 .and. mype==cputerm) write(*,*) 'end perturb_hydro'
         endif

         if(solve_energy) then
            if (debug>0 .and. mype==cputerm) write(*,*) 'start perturb_energy'
            call perturb_energy(err_energy,err_energy_prf,nit)
            if (debug>0 .and. mype==cputerm) write(*,*) 'end perturb_energy'
         endif

         err = max(err_hydro,err_energy)
         varerr=minval(abs(err-errold(:))/err)

         if(varerr.ge.accuracy)varerrsteady=0
         if(varerr.lt.accuracy)varerrsteady=varerrsteady+1

         if (err>accuracy) then
         !if (err>accuracy .and. varerrsteady<3) then

            nit = nit + 1

            ! After nit_minstep iterations change the minstep
            if (change_minstep) then
              if (nit==nit_minstep) then
                minstep = minstep*minstep_factor
                if (mype==cputerm) write(*,*) 'minstep change to ', minstep
              endif
            endif

            ! After nit_chem_smooth iterations turn on chem_smooth
            if (change_chem_smooth) then
              if (nit==nit_chem_smooth) then
                chem_smooth=.True.
                if (mype==cputerm) write(*,*) 'chem smooth turned on'
                if (mype==cputerm) write(*,*) 'kerchem_smooth = ', kerchem_smooth
              endif
              if (nit==nit_chem_smooth*1.5) then
                kerchem_smooth = kerchem_smooth + 1
                if (mype==cputerm) write(*,*) 'kerchem_smooth = ', kerchem_smooth
              endif
            endif

            ! After nit_gradad_smooth iterations turn on gradad_smooth
            if (nit==nit_gradad_smooth) then 
              gradad_smooth=.True.
              init_gradad=.True.
              if (mype==cputerm) write(*,*) 'gradad_smooth turned on'  
              if (mype==cputerm) write(*,*) 'kereos_smooth = ', kereos_smooth
            endif

            if (mod(nit,mod_tau) == 0) then
              if (check_tau) then
                call update_tau
              endif
            endif

            ! Set correction step size
            if(nit .gt. 1)then
               ! Increase step size if residuum falls monotonically
               if( (nit .eq. 2) .and. (err .lt. errold(nit-1)) )then
                  step=min(step*1.1,maxstep)
               endif
               if ((nit .gt. 2) .and. (err .lt. errold(nit-1))) then
                  if (errold(nit-1) .lt. errold(nit-2)) then
                     step=min(step*1.1,maxstep)
                  endif
               endif
               ! Decrease step size if residuum increases
               if(err .gt. errold(nit-1))then
                  step=max(step*0.5,minstep)
               endif
            else
               ! Set a small step size for the beginning
               step=minstep
            endif

           ! Hand over maximum
            errold(nit)=err
            errveryold_hydro = errold_hydro
            errold_hydro     = err_hydro
            errveryold_energy= errold_energy
            errold_energy = err_energy

            write(fnit,'(i4)') nit

            dpt = 0.
            if (debug>0 .and. mype==cputerm) write(*,*) 'start linsolve'
            call linsolve(2*(ndepth+1),jacobi,rhs,dpt)
            if (debug>0 .and. mype==cputerm) write(*,*) 'end linsolve'
        
            ! Set small step size if relative changes in P, T become large
            if ((maxval(abs(dpt(1:ndepth+1)/ppf(:))) .ge. dpt_limit) .or. (maxval(abs(dpt(ndepth+2:2*(ndepth+1))/ttf(:))) .ge. dpt_limit)) then
                step = minstep
            endif


              ! Get the model levels that have the maximum error
			  do i = 1, ndepth+1
			    if (err_hydro_prf(i) == err_hydro) then
			      max_hydro_level = i
			    endif

			    if (err_energy_prf(i) == err_energy) then
			    max_energy_level = i
			    endif

			    if (abs(dpt(i)/ppf(i)) == maxval(abs(dpt(1:ndepth+1)/ppf(:)))) then
			     max_dp_level = i
			    endif

			    if (abs(dpt(ndepth+1+i)/ttf(i)) ==  maxval(abs(dpt(ndepth+2:2*(ndepth+1))/ttf(:)))) then
			      max_dt_level = i
			    endif

			  enddo

            if (mype == cputerm) then
              write(*,headrfmt) 'PT Iter.: '//fnit
              write(*,1000) maxval(abs(dpt(1:ndepth+1)/ppf(:))),max_dp_level,err_hydro,max_hydro_level,step
              write(*,1001) maxval(abs(dpt(ndepth+2:2*(ndepth+1))/ttf(:))),max_dt_level,err_energy,max_energy_level
              write(*,*)

              !write(*,900) max_hydro_level,ppf(max_hydro_level)*unitP/1E6,dpt(max_hydro_level)/ppf(max_hydro_level),err_hydro_save(max_hydro_level),step
              !write(*,901) dpt(ndepth+1+max_hydro_level)/ttf(max_hydro_level),err_energy_save(max_hydro_level)
              !write(*,*)

              !write(*,900) max_energy_level,ppf(max_energy_level)*unitP/1E6,dpt(max_energy_level)/ppf(max_energy_level),err_hydro_save(max_energy_level),step
              !write(*,901) dpt(ndepth+1+max_energy_level)/ttf(max_energy_level),err_energy_save(max_energy_level)
              !write(*,*)

              !do i = 1,ndepth+1
              !  write(*,900) i,ppf(i)*unitP/1E6,dpt(i)/ppf(i),err_hydro_save(i),step
              !  write(*,901) dpt(ndepth+1+i)/ttf(i),err_energy_save(i)
              !  write(*,*)
              !enddo

              !where (abs(dpt(ndepth+2:2*(ndepth+1))/ttf(:))>0.1) dpt(ndepth+2:2*(ndepth+1)) = 0.

900         format(i3,5x,'p[bar]:',es10.1,5x,' dp/p: ',es10.1,5x,' err. hydro : ',es10.2,5x,' step : ',es10.1)
901	    format(30x,' dt/t: ',es10.1,5x,' err. energy: ',es10.2)
1000        format(5x,' dp/p: ',es10.1,' (',i3,')',5x,' err. hydro : ',es10.2,' (',i3,')',5x,' step : ',es10.1)
1001        format(5x,' dt/t: ',es10.1,' (',i3,')',5x,' err. energy: ',es10.2,' (',i3,')')

            ppf = ppf + step*dpt(1:ndepth+1)
            ttf = ttf + step*dpt(ndepth+2:2*(ndepth+1))

            kertt_coupled = 1

            imin = kertt_coupled+1
            imax = ndepth+1-kertt_coupled

            do idepth = imin,imax

               ttf_coupled(idepth) = ttf(idepth)

               do iker = 1,kertt_coupled

                  ttf_coupled(idepth) = ttf_coupled(idepth) +  temp_coupling*ttf(idepth-iker)+ temp_coupling*ttf(idepth+iker)

               enddo

               ttf_coupled(idepth) = ttf_coupled(idepth)/(2*kertt_coupled*temp_coupling+1)

            enddo

            do idepth = imin,imax
               if (ppf(idepth)>ppmin_temp_coupling .and. ppf(idepth)<ppmax_temp_coupling) then
                  ttf(idepth) = ttf_coupled(idepth)
               endif
            enddo

            endif

            if (rainout_freeze) then

              if (nit > rainout_freeze_nit .and. maxval(abs(dpt(1:ndepth+1)/ppf(:))) < rainout_freeze_dpt .and. maxval(abs(dpt(ndepth+2:2*(ndepth+1))/ttf(:))) < rainout_freeze_dpt) then
                rainout_level = maxval(rainout_levels)
                freeze_rainout = .True.
              else
                ! Get the rainout level
		        rainout_levels(cnt) = rainout_level
			    cnt = cnt + 1

			    if (cnt == rainout_freeze_past+1) then
			     cnt = 1
			    endif
                freeze_rainout = .False.
              endif

            endif

         else

            if (mype == cputerm) then
               write(fnit,'(i4)') (nit+1)
               write(*,headrfmt) 'Iteration: '//fnit
               write(*,1002) err_hydro,err_energy
               write(*,*)
               write(*,separfmt)
            endif

1002        format('err. hydro : ',es10.1,' err. energy: ',es10.1)

         endif

         if(mype == cputerm .and. debug>0) then
            do i = 1,ndepth+1
               write(*,1003) i,taufstd(i),ppf(i)*unitP,ttf(i)*unitK
            enddo

         endif

1003        format(5x,' idepth: ',i4,' tauf: ',es10.2,' ppf: ',es10.2' ttf: ',f10.0)

         if (minval(ttf)<0. .or. minval(ppf)<0.) then

           if(mype==cputerm) then

             write(*,headrfmt) 'Error find_profile'
             write(*,'(A)') 'Reached a negative pressure or temperature'
             write(*,exppmfmt) 'min pressure',minval(ppf)
             write(*,exppmfmt) 'min temperature',minval(ttf)
             write(*,*)

             do i = 1,ndepth+1
               write(*,1003) i,taufstd(i),ppf(i)*unitP,ttf(i)*unitK
             enddo

             write(*,separfmt)

             stop

           endif
         endif

         !check that pressure temperature is positive
         if(mype == cputerm) then
            if (debug>0) write(*,*) 'start write_pt_profile'
            call write_pt_profile(ftemp)
            if (debug>0) write(*,*) 'end write_pt_profile'
         endif

            ! do idepth =2,ndepth

            !    if((ttf(idepth)-ttf(idepth+1))> 25*unitK) then


            !       new_tau = 10.**(log10(taufstd(idepth))-20.*unitK*(log10(taufstd(idepth))-log10(taufstd(idepth+1)))/(ttf(idepth)-ttf(idepth+1)))

            !     !  if (abs(new_tau-taufstd(idepth+1))/taufstd(idepth+1)>0.01) then

            !          taufstd(idepth+1) = new_tau

            !          ppf(idepth+1) = 10.**(log10(ppf(idepth))-20.*unitK*(log10(ppf(idepth))-log10(ppf(idepth+1)))/(ttf(idepth)-ttf(idepth+1)))

            !          ttf(idepth+1) = ttf(idepth)-20*unitK
            !          write(*,*) '****',idepth,taufstd(idepth+1)
            !      ! endif


            !    else if ((ttf(idepth)-ttf(idepth+1))< -25*unitK) then

            !       new_tau = 10.**(log10(taufstd(idepth))+20.*unitK*(log10(taufstd(idepth))-log10(taufstd(idepth+1)))/(ttf(idepth)-ttf(idepth+1)))

            !       !if (abs(new_tau-taufstd(idepth+1))/taufstd(idepth+1)>0.01) then

            !          taufstd(idepth+1) = new_tau

            !          ppf(idepth+1) = 10.**(log10(ppf(idepth))+20.*unitK*(log10(ppf(idepth))-log10(ppf(idepth+1)))/(ttf(idepth)-ttf(idepth+1)))

            !          ttf(idepth+1) = ttf(idepth)+20*unitK
            !          write(*,*) '****',idepth,taufstd(idepth+1)
            !       !endif

            !    endif


            ! enddo


            ! interpolation at cell centres
            do i =1,ndepth

               pp(i) = sqrt(ppf(i)*ppf(i+1))
               tt(i) = sqrt(ttf(i)*ttf(i+1))

            enddo

            nn  = get_nn(ndepth,pp,tt)
            nnf = get_nn(ndepth+1,ppf,ttf)
            rho = get_rho(ndepth,pp,tt,amol)
            rhof = get_rho(ndepth+1,ppf,ttf,amolf)

            call get_rrgr

         enddo

         if (calc_radtime) then
            call compute_radtime
         endif

! Add code here to stop it writing to the file at the maximum number of iterations.
!      if (err<accuracy .and. nit<nitmax) then
         if(mype == cputerm) then
            if (debug>0) write(*,*) 'start write_pt_profile'
            if (nit == nitmax)then
              call write_pt_profile(ftemp)
            else
              call write_pt_profile(fout)
            endif
            if (debug>0) write(*,*) 'end write_pt_profile'
         endif
!      endif

       end subroutine find_profile

 subroutine uvtrans_pt

   use mod_grid
   use mod_radtrans
   use mod_chem
   use mod_eos
   use mod_cst

   implicit none

   integer :: idepth,ise,ilambda,i,nse,iray
   integer,parameter::nlam_uv = 315
   integer :: nrays_uv
   logical :: scatter_uv
   real :: ilam,dlambda
   real,dimension(nraymax) :: alphap,sig

   real,dimension(ndepth+1,nlam_uv) :: fstar_uv,snu_uv,bnu_uv,hnup_uv,jnup_uv,hnud_uv,jnud_uv,snu_scat_loc
   real,dimension(ndepth  ,nlam_uv) :: kapnu_uv,kscatnu_uv,epsnu_uv,fuv_loc
   real,dimension(ndepth+1,nlam_uv) :: kapfnu_uv,kscatfnu_uv,epsfnu_uv,fuv_nu
   real,dimension(ndepth           ) :: kapstd_uv,mu
   real,dimension(ndepth+1         ) :: muf,frad_uv
   real,dimension(         nlam_uv) :: Irad_uv_loc

   real,dimension(:,:),allocatable :: se			!< cross sections
   real,dimension(:,:),allocatable :: selambda	!< wavelength arrays of cross section files
   integer,dimension(:),allocatable :: semol
   character(30)     :: dir_neq,fhnu_uv  !< name of the directory containing the non-eq network
   character(12), dimension(:), allocatable :: se_name
   logical :: first
   character(2) :: rd

   character(2) :: slname
   character(lname) :: semolname

   first = .true.
   dir_neq = '../../chem/venot2012/'
   fhnu_uv = 'hnu_uv_hd209.dat'

   nse=0
   open(1,file=trim(dir_neq)//'filenames_se.txt')
   do
      read(1,'(a2)',end=400) rd
      nse = nse + 1
   enddo
400 close(1)

  allocate(se_name(nse))
  allocate(semol(nse))
  allocate(selambda(nse,nlam_uv),se(nse,nlam_uv))
  se = 0.; selambda = 0.

  open(1,file=trim(dir_neq)//'filenames_se.txt')
2016 format(a13)
  do ise = 1, nse
     read(1,2016,end=16) se_name(ise)
  enddo
16 close(1)

  !Read in cross sections and branching ratios
  do ise = 1, nse
     ilambda = 0
     open(1,file=trim(dir_neq)//trim(se_name(ise))//'.dat')
2017 format(e11.5,4x,e11.5)
     do
        ilambda = ilambda + 1

        if(ilambda .le. nlam_uv) then

!!! selambda is in nm and the cross section 'se' in cm**2
           read(1,*,end=17) selambda(ise,ilambda),se(ise,ilambda)
        else
           read(1,*,end=17)
        endif
     enddo
17   close(1)
  enddo

   open(1,file=trim(dir_neq)//'filenames_se.txt')
   do ise = 1, nse
     write(slname,'(i2)') lname
     read(1,'(2x,a'//slname//')') semolname
     semol(ise) = imol(trim(semolname))
   end do
   close(1)

   dlambda = 1.

   mu  = get_mu(ndepth,amol)
   muf = get_mu(ndepth+1,amolf)

   nrays_uv = 16
   scatter_uv = .true.

   fstar_uv = 0. ; snu_uv = 0. ; bnu_uv = 0. ; hnup_uv = 0. ; jnup_uv = 0.
   hnud_uv = 0. ; jnud_uv = 0. ; fuv = 0.
   Irad_uv_loc = 0. ; kapstd_uv = 0.
   kapnu_uv  = 0. ; kscatnu_uv  = 0. ; epsnu_uv  = 1.
   kapfnu_uv = 0. ; kscatfnu_uv = 0. ; epsfnu_uv = 1.

   !compute kapstd_uv from hydrostatic equilibrium
   do idepth = 1,ndepth

      kapstd_uv(idepth) = gr(idepth)*(taufstd(idepth+1)-taufstd(idepth))/(ppf(idepth+1)-ppf(idepth))

   enddo




   !compute the absorption opacities

   do ise = 1, nse

      i = 0

      do ilambda = 1, nlam_uv

         if (ilambda == selambda(ise,i+1)) then
            i = i + 1
            kapnu_uv(:,ilambda) = kapnu_uv(:,ilambda) + se(ise,i)/unitL**2*amol(:,semol(ise))/(mu(:)*amu)

            kapfnu_uv(:,ilambda) = kapfnu_uv(:,ilambda) + se(ise,i)/unitL**2*amolf(:,semol(ise))/(muf(:)*amu)

         endif

      enddo
   enddo


    !compute the scattering opacities

   do ilambda = 1, nlam_uv

      ilam = ilambda*1E-7/unitL

     do iray = 1,num_ray_scatt

      alphap(iray) = (nr(iray)**2-1.)/(4.*pi*nl)
      sig(iray)  = 8.*pi/3.*alphap(iray)**2*(6.+3.*p_fac(iray))/(6.-7.*p_fac(iray))

      kscatnu_uv(:,ilambda) = (sig(iray)*(2.*pi/ilam)**4)*amol(:,imol(rayname(iray)))/(mu(:)*amu)+(sig(iray)*(2.*pi/ilam)**4)*amol(:,imol(rayname(iray)))/(mu(:)*amu)

      kscatfnu_uv(:,ilambda) = (sig(iray)*(2.*pi/ilam)**4)*amolf(:,imol(rayname(iray)))/(muf(:)*amu)+(sig(iray)*(2.*pi/ilam)**4)*amolf(:,imol(rayname(iray)))/(muf(:)*amu)
    enddo

   enddo

   !compute epsnu and total opacities
   where (kapnu_uv + kscatnu_uv /=0.)
      epsnu_uv = kapnu_uv / (kapnu_uv + kscatnu_uv)
      kapnu_uv = kapnu_uv + kscatnu_uv
   endwhere

   !compute epsfnu and total opacities
   where (kapfnu_uv + kscatfnu_uv /=0.)
      epsfnu_uv = kapfnu_uv / (kapfnu_uv + kscatfnu_uv)
      kapfnu_uv = kapfnu_uv + kscatfnu_uv
   endwhere

   !get Irad
   open(unit=1,file=fhnu_uv)


   do ilambda = 1, nlam_uv

      read(1,*) ilam,Irad_uv_loc(ilambda)

   enddo
   close(1)

   call startrans(nlam_uv,kapnu_uv,kapstd_uv,fstar_uv,Irad_uv_loc)

   !compute the source function
   if (scatter_uv) then
      snu_uv = (1.-epsfnu_uv)*abs(fstar_uv/murad)/(4.*pi)
   else
      !compute planck function
      snu_uv = 0.
   endif

   ! boundary condition at the bottom layer
   bnu_uv = snu_uv

   snu_scat_loc = snu_uv


   call radtrans(nrays_uv,scatter_uv,nlam_uv,kapnu_uv,kapstd_uv,epsfnu_uv,snu_uv,bnu_uv,hnup_uv,hnud_uv,jnup_uv,jnud_uv,snu_scat_loc)

   do idepth=1,ndepth

      fuv_nu(idepth,1:nlam_uv) = -(fstar_uv(idepth,1:nlam_uv)+4.*pi*(hnup_uv(idepth,1:nlam_uv)+hnud_uv(idepth,1:nlam_uv)))

   enddo

   fuv = 0.

   do idepth = 1,ndepth+1

      do ilambda = 1,nlam_uv
      fuv(idepth) = fuv(idepth)+(fstar_uv(idepth,ilambda)+4.*pi*(hnup_uv(idepth,ilambda)+hnud_uv(idepth,ilambda)))/(1./unitT/unitL**2)*hplanck*cvel/(ilambda*1E-7/unitL)*dlambda
   enddo
enddo



 end subroutine uvtrans_pt


 subroutine compute_radtime

   use mod_grid
   use mod_param
   use mod_opacity
   use mod_cst
   use mod_chem
   use mod_util
   use mod_eos
   use mod_radtrans

   implicit none

   integer :: i,idepth
   real    :: eps
   real,dimension(:,:),allocatable :: jnu,jnup,jnud,hnu,hnup,hnud,bnu,snu
   real,dimension(ndepth+1)        :: cpspec,gradad_loc
   character(len=6),dimension(ndepth+1)  :: eos_flag

   allocate(hnu(ndepth+1,nfreqkap),jnu(ndepth+1,nfreqkap),snu(ndepth+1,nfreqkap))
   allocate(hnup(ndepth+1,nfreqkap),jnup(ndepth+1,nfreqkap),bnu(ndepth+1,nfreqkap))
   allocate(hnud(ndepth+1,nfreqkap),jnud(ndepth+1,nfreqkap))

   eps = 1.0E-2

   do idepth = 2,ndepth
      write(*,*) 'radiative timescale computation ideth = ',idepth

      hnu = 0.; jnu = 0.; snu = 0.; hnup = 0.; jnup = 0.; bnu = 0.; hnud = 0.; jnud = 0.

      ttf(idepth) = ttf(idepth)*(1+eps)

      if (any(chem==(/'eq ','ana','man'/))) then
         if (debug>0 .and. mype==cputerm) write(*,*) 'start get_eq_chemistry'
         call get_eq_chemistry
         if (debug>0 .and. mype==cputerm) write(*,*) 'end get_eq_chemistry'

      else if(chem == 'neq') then

         if (mype==cputerm) then

            write(*,headrfmt) 'Error calc_radtime'
            write(*,'(A)') 'neq cannot be used inside calc_radtime'
            write(*,chrpmfmt) 'chem',chem
            write(*,*)
            write(*,separfmt)
            stop

         endif

      endif

      nn  = get_nn(ndepth,pp,tt)
      nnf = get_nn(ndepth+1,ppf,ttf)

      rho = get_rho(ndepth,pp,tt,amol)
      rhof = get_rho(ndepth+1,ppf,ttf,amolf)

      call compute_eos(ndepth+1,ppf,ttf,rhof,cpspec,gradad_loc,eos_flag)

      call get_rrgr

      if(debug>0 .and. mype==cputerm) write(*,*) 'start get_opacity'
      if (corr_k) then
         call get_corr_k_opacity
      else
         call get_lbl_opacity
      endif
      if (hminus) then
         call get_hm_opacity
      endif
      if (scatter) then
         call get_scat_opacity
      endif
      if(debug>0 .and. mype==cputerm) write(*,*) 'end get_opacity'

      !compute planck function
      bnu      = get_bnu (ndepth+1,nfreqkap,ttf)

      fstar_nu = 0.

      !compute star flux
      if (irrad) then
         call startrans(nfreqkap,kapnu,kapstd,fstar_nu,Irad)
      endif

      !source term for RT
      if (scatter) then
         snu = epsfnu*bnu+(1.-epsfnu)*abs(fstar_nu/murad)/(4.*pi)
      else
         !compute planck function
         snu      = bnu
      endif

      if (corr_k) then
         where(weightsfnu == 0.) snu = 0.
      endif

      !compute radiative flux at interfaces
      if(mype==cputerm .and. debug>0) write(*,*) 'start radtrans'
      call radtrans(nrays,scatter,nfreqkap,kapnu,kapstd,epsfnu,snu,bnu,hnup,hnud,jnup,jnud,snu_scat)
      if(mype==cputerm .and. debug>0) write(*,*) 'end radtrans'

      hnu = hnup+hnud
      jnu = jnup+jnud

      ! compute total irradiation and add uv
      if (irrad) then

         fstar    = int_freq(ndepth+1,nfreqkap,nuf,fstar_nu)

         if (include_uv) then

            fuv = 0.
            call uvtrans_pt
            fstar = fstar + fuv

         endif

      endif

      !compute radiative flux and heating rate
      do i=1,ndepth+1
         qnu(i,:) = 4.*pi*epsfnu(i,:)*(jnu(i,:)-snu(i,:))*kapfnu(i,:)*rhof(i)
         fnu(i,:) = 4.*pi*hnu(i,:)
      enddo

      !integrate the flux and heating rate
      frad = int_freq(ndepth+1,nfreqkap,nuf,fnu)
      qrad = int_freq(ndepth+1,nfreqkap,nuf,qnu)

      ttf(idepth) = ttf(idepth)/(1.+eps)

      radtime(idepth) = eps*ttf(idepth)*rhof(idepth)*cpspec(idepth)/abs((frad(idepth+1)-frad(idepth))/(rrf(idepth+1)-rrf(idepth)))

   enddo

   radtime(1) = radtime(2)**2/radtime(3)
   radtime(ndepth+1) = radtime(ndepth)**2/radtime(ndepth-1)

 end subroutine compute_radtime



 subroutine fill_matrices_2d(err)

   use mod_grid
   use mod_radtrans
   use mod_opacity
   use mod_chem
   use mod_cst
   use mod_eos
   use mod_util
   use mod_conv
   use mod_param

   implicit none

   real,dimension(ndepth+1,nfreqkap) :: jnup,jnud,hnup,hnud,bnu,dbnu,snu,dsnu,dhnup,dhnud,djnup,djnud,dfnu
   real,dimension(ndepth+1,0:ntheta+2) :: rhs_momr,rhs_momt,rhs_ener,dmomr,dmomt,dener,dppf,dttf,drhof,dfrad,dutf,dHf,Hf_2d,rhs_mass,dmass,ddupf,rhs_momp,dmomp,durf,csf_2d
    real,dimension(ndepth+1,2*(ndepth+1),0:ntheta+2) :: dfconv_2d
   real :: err,fint,epsilon, coseps,epsmu,flux_norm
   real,dimension(1,nfreqkap) :: bnusurf
   real,dimension(1) :: int_fsurf
   integer :: it,idepth,irhs,ijac,ndepth_um
   real,dimension(ndepth+1) :: Fmom


   bnusurf = get_bnu(1,nfreqkap,(/ Teff /))
   int_fsurf = int_freq(1,nfreqkap,nuf,bnusurf)
   fint = int_fsurf(1)*pi

   err = 0.

   ! compute fluxes
   do it = 1,ntheta

      ppf(:) = ppf_2d(:,it)
      ttf(:) = ttf_2d(:,it)

      do idepth =1,ndepth

         pp(idepth) = sqrt(ppf(idepth)*ppf(idepth+1))
         tt(idepth) = sqrt(ttf(idepth)*ttf(idepth+1))

      enddo

      nn  = get_nn(ndepth,pp,tt)
      nnf = get_nn(ndepth+1,ppf,ttf)

      call get_eq_chemistry

      rho = get_rho(ndepth,pp,tt,amol)
      rhof = get_rho(ndepth+1,ppf,ttf,amolf)
      muf_2d (:,it) = get_mu(ndepth+1,Amolf)
      rhof_2d(:,it) = get_rho(ndepth+1,ppf,ttf,amolf)
      cpf_2d (:,it) = get_cpspec(ndepth+1,ppf,ttf,Amolf)/(kb/(muf_2d(:,it)*amu))
      Hf_2d(:,it) = kb*ttf_2d(:,it)/(amu*muf_2d(:,it)*grf(:))
      csf_2d(:,it) = sqrt(ppf(:)/rhof(:))


      if(cloud_night_side) then

         if(cos(thetaf(it))<0.) then

            cloud = .true.

         endif

      endif

      if (corr_k) then
        call get_corr_k_opacity
      else
         write(*,*) 'err 2d needs corr_k'
         stop
      endif

      cloud = .false.

      if (hminus) then
         call get_hm_opacity
      endif
      if (scatter) then
         call get_scat_opacity
      endif

      kapnu_2d(:,it,:) = kapnu(:,:)
      kapstd_2d(:,it) = kapstd(:)
      epsfnu_2d(:,it,:) = epsfnu(:,:)

      taufstd(1) = taumin
      taufstd_2d(1,it) = taufstd(1)

      do idepth = 2,ndepth+1

         taufstd(idepth ) = taufstd(idepth-1)-kapstd(idepth-1)*rho(idepth-1)*(rrf(idepth)-rrf(idepth-1))
         taufstd_2d(idepth,it) = taufstd(idepth)

      enddo

      bnu      = get_bnu (ndepth+1,nfreqkap,ttf)
      dbnu     = get_dbnu(ndepth+1,nfreqkap,ttf)

      where(weightsfnu == 0.) dbnu = 0.

      dsnu_2d(:,it,:) = dbnu(:,:)*epsfnu(:,:)

      fstar_nu = 0.

      !compute stellar irradiation angle
      murad = cos(thetaf(it))

      frad_2d(:,it) = 0.

      !compute star flux
      if (irrad) then
         if (murad>0.) then
            call startrans(nfreqkap,kapnu,kapstd,fstar_nu,Irad)
            frad_2d(:,it) = int_freq(ndepth+1,nfreqkap,nuf,fstar_nu)
         endif
      endif

      flux_norm = abs(frad_2d(1,1))

      !source term for RT
      if (scatter) then
         if (murad >0.) then
            snu = epsfnu*bnu+(1.-epsfnu)*abs(fstar_nu)/(murad*4.*pi)
         else
            snu = epsfnu*bnu
         endif

      else
         !compute planck function
         snu      = bnu
      endif

      where(weightsfnu == 0.) snu = 0.

      hnup = 0. ; hnud = 0. ; jnup = 0. ; jnud = 0. ; fnu=0.

      !compute planet/diffuse flux
      call radtrans(nrays,scatter,nfreqkap,kapnu,kapstd,epsfnu,snu,bnu,hnup,hnud,jnup,jnud,snu_scat_2d(:,it,:))

      !compute radiative flux
      fnu = 4.*pi*(hnup+hnud)
      frad_2d(:,it) = frad_2d(:,it) + int_freq(ndepth+1,nfreqkap,nuf,fnu)


      call get_convection

      frad_2d(:,it) = frad_2d(:,it) + fconv(:)
      dfconv_2d(:,:,it) = dfconv(:,:)

   enddo



   !initialize winds
   if (first) then

      do idepth=2,ndepth

         call interp_1d(log10(p_2din),u_2din,log10(ppf_2d(idepth,1)),utf_2d(idepth,1))

      enddo


      do idepth =1,ndepth+1
           do it=1,ntheta

              utf_2d(idepth,it)  = utf_2d(idepth,1)
              dupf_2d(idepth,it) = 0.
              urf_2d(idepth,it)  = 0.

         enddo
      enddo

      first = .false.

   endif


   !boundary conditions (com in parallel to be done!)
   ppf_2d(:,0) = ppf_2d(:,ntheta) ; ttf_2d(:,0) = ttf_2d(:,ntheta) ; rhof_2d(:,0) = rhof_2d(:,ntheta)
   utf_2d(:,0) =utf_2d(:,ntheta) ; cpf_2d(:,0) = cpf_2d(:,ntheta) ; frad_2d(:,0) = frad_2d(:,ntheta)
   Hf_2d(:,0) = Hf_2d(:,ntheta) ; dupf_2d(:,0) = dupf_2d(:,ntheta) ; urf_2d(:,0) = urf_2d(:,ntheta)

   ppf_2d(:,ntheta+1) = ppf_2d(:,1) ; ttf_2d(:,ntheta+1) = ttf_2d(:,1) ; rhof_2d(:,ntheta+1) = rhof_2d(:,1)
   utf_2d(:,ntheta+1) =utf_2d(:,1) ; cpf_2d(:,ntheta+1) = cpf_2d(:,1) ; frad_2d(:,ntheta+1) = frad_2d(:,1)
   Hf_2d(:,ntheta+1) = Hf_2d(:,1) ; dupf_2d(:,ntheta+1) = dupf_2d(:,1) ; urf_2d(:,ntheta+1) = urf_2d(:,1)


   !fill the rhs terms
   rhs_momr = mom_r_2d(ppf_2d,rhof_2d,Hf_2d,utf_2d,urf_2d)
   rhs_momt = mom_t_2d(ppf_2d,rhof_2d,utf_2d,dupf_2d,urf_2d,cpf_2d)
   rhs_ener = energy_2d(ppf_2d,rhof_2d,ttf_2d,utf_2d,dupf_2d,urf_2d,frad_2d,fint,cpf_2d)
   rhs_mass = mass_2d(rhof_2d,ppf_2d,utf_2d,dupf_2d,urf_2d)
   rhs_momp = mom_p_2d(ppf_2d,rhof_2d,ttf_2d,utf_2d,dupf_2d,urf_2d,frad_2d,fint,cpf_2d)

   !compute error
   do it = 1,ntheta

      do idepth=2,ndepth
         if(ppf_2d(idepth,it)>pp2d_decoupling_err) then

            err = max(err,(abs(rhs_momr(idepth  ,it))/(grf(idepth)*rmax)))
            err = max(err,(abs(rhs_momt(idepth  ,it))/(grf(idepth)*rmax)))
            err = max(err,(abs(rhs_momp(idepth  ,it))/sqrt(grf(idepth)*rmax)   ))
            err = max(err,(abs(rhs_ener(idepth  ,it))/max(abs(cpf_2d(idepth,it)*ppf_2d(idepth,it)*sqrt(grf(idepth)*rmax)),abs(flux_norm*rmax/(rrf(1)-rrf(2))))))
            err = max(err,(abs(rhs_mass(idepth  ,it))/sqrt(grf(idepth)*rmax)   ))
         endif
      enddo
      err = max(err,abs(rhs_ener(ndepth+1  ,it))/abs(flux_norm*rmax/(rrf(1)-rrf(2))))

   enddo

    ! terminal output max err function of idepth
   if (mype == cputerm) then

      write(*,*) 'idepth    momr_err    momt_err    momp_err    ener_err    mass_err'

      do idepth=1,ndepth

         write(*,3000) idepth,                                       &
              maxval(abs(rhs_momr(idepth,1:ntheta))/(grf(idepth)*rmax)), &
              maxval(abs(rhs_momt(idepth,1:ntheta))/(grf(idepth)*rmax)), &
              maxval(abs(rhs_momp(idepth,1:ntheta))/sqrt(grf(idepth)*rmax)   ), &
              maxval(abs(rhs_ener(idepth,1:ntheta))/max(abs(cpf_2d(idepth,1:ntheta)*ppf_2d(idepth,1:ntheta)*sqrt(grf(idepth)*rmax)),abs(flux_norm*rmax/(rrf(1)-rrf(2))))), &
              maxval(abs(rhs_mass(idepth,1:ntheta))/sqrt(grf(idepth)*rmax)   )

      enddo
   endif

    write(*,3000) ndepth+1,                                       &
              maxval(abs(rhs_momr(ndepth+1,1:ntheta))/(grf(ndepth+1)*rmax)), &
              maxval(abs(rhs_momt(ndepth+1,1:ntheta))/(grf(ndepth+1)*rmax)), &
              maxval(abs(rhs_momp(ndepth+1,1:ntheta))/sqrt(grf(ndepth+1)*rmax)   ), &
              maxval(abs(rhs_ener(ndepth+1,1:ntheta))/abs((flux_norm*rmax/(rrf(1)-rrf(2))))), &
              maxval(abs(rhs_mass(ndepth+1,1:ntheta))/sqrt(grf(ndepth+1)*rmax)   )



3000 format(3x,i3,5(2x,es10.1))


   do idepth=1,ndepth
      do it = 1,ntheta

         irhs = idepth+(it-1)*(ndepth+1)
         rhs(irhs) = -rhs_momr(idepth,it)

         irhs = 1*(ndepth+1)*(ntheta+1)+idepth+(it-1)*(ndepth+1)
         rhs(irhs) = -rhs_momt(idepth,it)

         irhs = 2*(ndepth+1)*(ntheta+1)+idepth+(it-1)*(ndepth+1)
         rhs(irhs) = -rhs_ener(idepth,it)

         irhs = 3*(ndepth+1)*(ntheta+1)+idepth+(it-1)*(ndepth+1)
         rhs(irhs) = -rhs_mass(idepth,it)

         irhs = 4*(ndepth+1)*(ntheta+1)+idepth+(it-1)*(ndepth+1)
         rhs(irhs) = -rhs_momp(idepth,it)

      enddo
   enddo

   do it = 1,ntheta

      irhs = 2*(ndepth+1)*(ntheta+1)+ndepth+1+(it-1)*(ndepth+1)
      rhs(irhs) = -rhs_ener(ndepth+1,it)

   enddo

   !perturb the equations and fill jacobi matrix
   dppf = 0. ; dttf = 0. ; drhof = 0. ; dutf = 0. ; dfrad = 0. ; dHf = 0. ; ddupf = 0. ; durf = 0.

   do it = 1,ntheta

      taufstd = taufstd_2d(:,it)

      do idepth=2,ndepth

         dppf(idepth,it)  = 1.

         dppf(1,1:ntheta) = 2.*dppf(2,1:ntheta)*ppf_2d(2,1:ntheta)/ppf_2d(3,1:ntheta)-dppf(3,1:ntheta)*ppf_2d (2,1:ntheta)**2/ppf_2d(3,1:ntheta)**2

         drhof(:,1:ntheta) = dppf(:,1:ntheta)*(muf_2d(:,1:ntheta)*amu)/(kb*ttf_2d(:,1:ntheta))

         dfrad(:,it) = dfconv_2d(:,idepth,it)

         dppf(:,0)         = dppf(:,ntheta)
         dppf(:,ntheta+1)  = dppf(:,1)


         drhof(:,0)         = drhof(:,ntheta)
         drhof(:,ntheta+1)  = drhof(:,1)

         dfrad(:,0) = dfrad(:,ntheta)
         dfrad(:,ntheta+1) = dfrad(:,1)

         dmomr = dmom_r_2d (dppf,drhof,dHf,dutf,durf,ppf_2d,rhof_2d,Hf_2d,utf_2d,urf_2d)
         dmomt = dmom_t_2d (dppf,drhof,dutf,ddupf,durf,ppf_2d,rhof_2d,utf_2d,dupf_2d,urf_2d,cpf_2d)
         dener = denergy_2d(dppf,drhof,dttf,dutf,ddupf,durf,dfrad,ppf_2d,rhof_2d,ttf_2d,utf_2d,dupf_2d,urf_2d,frad_2d,fint,cpf_2d)
         dmass = dmass_2d(drhof,dppf,dutf,ddupf,durf,rhof_2d,ppf_2d,utf_2d,dupf_2d,urf_2d)
         dmomp = dmom_p_2d(dppf,drhof,dttf,dutf,ddupf,durf,dfrad,ppf_2d,rhof_2d,ttf_2d,utf_2d,dupf_2d,urf_2d,frad_2d,fint,cpf_2d)

         ijac = idepth+(it-1)*(ndepth+1)

         call write_jacobi(ijac,dmomr,dmomt,dener,dmomp,dmass)

         dppf  = 0. ; drhof = 0. ; dfrad = 0.

         dutf(idepth,it) = 1.
         dutf(1,:) = 2.*dutf(2,:)-dutf(3,:)

         dutf(idepth,0) = dutf(idepth,ntheta)
         dutf(idepth,ntheta+1) = dutf(idepth,1)


         dmomr = dmom_r_2d (dppf,drhof,dHf,dutf,durf,ppf_2d,rhof_2d,Hf_2d,utf_2d,urf_2d)
         dmomt = dmom_t_2d (dppf,drhof,dutf,ddupf,durf,ppf_2d,rhof_2d,utf_2d,dupf_2d,urf_2d,cpf_2d)
         dener = denergy_2d(dppf,drhof,dttf,dutf,ddupf,durf,dfrad,ppf_2d,rhof_2d,ttf_2d,utf_2d,dupf_2d,urf_2d,frad_2d,fint,cpf_2d)
         dmass = dmass_2d(drhof,dppf,dutf,ddupf,durf,rhof_2d,ppf_2d,utf_2d,dupf_2d,urf_2d)
         dmomp = dmom_p_2d(dppf,drhof,dttf,dutf,ddupf,durf,dfrad,ppf_2d,rhof_2d,ttf_2d,utf_2d,dupf_2d,urf_2d,frad_2d,fint,cpf_2d)

         ijac = (ndepth+1)*(ntheta+1)+idepth+(it-1)*(ndepth+1)
         call write_jacobi(ijac,dmomr,dmomt,dener,dmomp,dmass)

         dutf = 0.

         dttf(idepth,it)  = 1.

         dttf(1,:) = 2.*dttf(2,:)-dttf(3,:)


         drhof(:,1:ntheta) = -ppf_2d(:,1:ntheta)*muf_2d(:,1:ntheta)*amu/(kb*ttf_2d(:,1:ntheta)**2)*dttf(:,1:ntheta)
         dHf(idepth,it) = kb/(amu*muf_2d(idepth,it)*grf(idepth))

         !reset the perturbation
         dsnu = 0. ; dfnu = 0. ; dhnup = 0. ; dhnud = 0. ; djnup = 0. ; djnud = 0.

         dsnu(idepth,:) = dsnu_2d(idepth,it,:)

         call radtrans(1,.false.,nfreqkap,kapnu_2d(:,it,:),kapstd_2d(:,it),epsfnu_2d(:,it,:),dsnu,dsnu,dhnup,dhnud,djnup,djnud,dsnu)

         dfnu = 4.*pi*(dhnup+dhnud)
         dfrad(:,it) =  int_freq(ndepth+1,nfreqkap,nuf,dfnu)

         dfrad(:,it) = dfrad(:,it) + dfconv_2d(:,ndepth+1+idepth,it)

         dttf(:,0) = dttf(:,ntheta)
         drhof(:,0) = drhof(:,ntheta)
         dHf(idepth,0) = dHf(idepth,ntheta)
         dfrad(:,0) = dfrad(:,ntheta)

        dttf(:,ntheta+1) = dttf(:,1)
        drhof(:,ntheta+1) = drhof(:,1)
        dHf(idepth,ntheta+1) = dHf(idepth,1)
        dfrad(:,ntheta+1) = dfrad(:,1)

         dmomr = dmom_r_2d (dppf,drhof,dHf,dutf,durf,ppf_2d,rhof_2d,Hf_2d,utf_2d,urf_2d)
         dmomt = dmom_t_2d (dppf,drhof,dutf,ddupf,durf,ppf_2d,rhof_2d,utf_2d,dupf_2d,urf_2d,cpf_2d)
         dener = denergy_2d(dppf,drhof,dttf,dutf,ddupf,durf,dfrad,ppf_2d,rhof_2d,ttf_2d,utf_2d,dupf_2d,urf_2d,frad_2d,fint,cpf_2d)
         dmass = dmass_2d(drhof,dppf,dutf,ddupf,durf,rhof_2d,ppf_2d,utf_2d,dupf_2d,urf_2d)
         dmomp = dmom_p_2d(dppf,drhof,dttf,dutf,ddupf,durf,dfrad,ppf_2d,rhof_2d,ttf_2d,utf_2d,dupf_2d,urf_2d,frad_2d,fint,cpf_2d)


         ijac = 2*(ndepth+1)*(ntheta+1)+idepth+(it-1)*(ndepth+1)

         call write_jacobi(ijac,dmomr,dmomt,dener,dmomp,dmass)

         dttf = 0. ; drhof = 0. ; dfrad = 0. ; dhf = 0.


         ddupf(idepth,it)  = 1.
         ddupf(1,:) = 2.*ddupf(2,:)-ddupf(3,:)

         ddupf(idepth,0)         = ddupf(idepth,ntheta)
         ddupf(idepth,ntheta+1)  = ddupf(idepth,1)


         dmomr = dmom_r_2d (dppf,drhof,dHf,dutf,durf,ppf_2d,rhof_2d,Hf_2d,utf_2d,urf_2d)
         dmomt = dmom_t_2d (dppf,drhof,dutf,ddupf,durf,ppf_2d,rhof_2d,utf_2d,dupf_2d,urf_2d,cpf_2d)
         dener = denergy_2d(dppf,drhof,dttf,dutf,ddupf,durf,dfrad,ppf_2d,rhof_2d,ttf_2d,utf_2d,dupf_2d,urf_2d,frad_2d,fint,cpf_2d)
         dmass = dmass_2d(drhof,dppf,dutf,ddupf,durf,rhof_2d,ppf_2d,utf_2d,dupf_2d,urf_2d)
         dmomp = dmom_p_2d(dppf,drhof,dttf,dutf,ddupf,durf,dfrad,ppf_2d,rhof_2d,ttf_2d,utf_2d,dupf_2d,urf_2d,frad_2d,fint,cpf_2d)

         ijac = 3*(ndepth+1)*(ntheta+1)+idepth+(it-1)*(ndepth+1)

         call write_jacobi(ijac,dmomr,dmomt,dener,dmomp,dmass)

         ddupf  = 0.


         durf(idepth,it)  = 1.

         durf(idepth,0)         = durf(idepth,ntheta)
         durf(idepth,ntheta+1)  = durf(idepth,1)


         dmomr = dmom_r_2d (dppf,drhof,dHf,dutf,durf,ppf_2d,rhof_2d,Hf_2d,utf_2d,urf_2d)
         dmomt = dmom_t_2d (dppf,drhof,dutf,ddupf,durf,ppf_2d,rhof_2d,utf_2d,dupf_2d,urf_2d,cpf_2d)
         dener = denergy_2d(dppf,drhof,dttf,dutf,ddupf,durf,dfrad,ppf_2d,rhof_2d,ttf_2d,utf_2d,dupf_2d,urf_2d,frad_2d,fint,cpf_2d)
         dmass = dmass_2d(drhof,dppf,dutf,ddupf,durf,rhof_2d,ppf_2d,utf_2d,dupf_2d,urf_2d)
         dmomp = dmom_p_2d(dppf,drhof,dttf,dutf,ddupf,durf,dfrad,ppf_2d,rhof_2d,ttf_2d,utf_2d,dupf_2d,urf_2d,frad_2d,fint,cpf_2d)

         ijac = 4*(ndepth+1)*(ntheta+1)+idepth+(it-1)*(ndepth+1)

         call write_jacobi(ijac,dmomr,dmomt,dener,dmomp,dmass)

         durf  = 0.


      enddo

      dttf(ndepth+1,it)  = 1.
      drhof(ndepth+1,it) = -ppf_2d(ndepth+1,it)*muf_2d(ndepth+1,it)*amu/(kb*ttf_2d(ndepth+1,it)**2)
      dHf(ndepth+1,it) = kb/(amu*muf_2d(ndepth+1,it)*grf(ndepth+1))

      !reset the perturbation
      dsnu = 0. ; dfnu = 0. ; dhnup = 0. ; dhnud = 0. ; djnup = 0. ; djnud = 0.

      dsnu(ndepth+1,:) = dsnu_2d(ndepth+1,it,:)

      call radtrans(1,.false.,nfreqkap,kapnu_2d(:,it,:),kapstd_2d(:,it),epsfnu_2d(:,it,:),dsnu,dsnu,dhnup,dhnud,djnup,djnud,dsnu)

      dfnu = 4.*pi*(dhnup+dhnud)
      dfrad(:,it) =  int_freq(ndepth+1,nfreqkap,nuf,dfnu)

      dfrad(:,it) = dfrad(:,it) + dfconv_2d(:,2*(ndepth+1),it)

      dttf(ndepth+1,0) = dttf(ndepth+1,ntheta)
      drhof(ndepth+1,0) = drhof(ndepth+1,ntheta)
      dHf(ndepth+1,0) = dHf(ndepth+1,ntheta)
      dfrad(:,0) = dfrad(:,ntheta)

      dttf(ndepth+1,ntheta+1) = dttf(ndepth+1,1)
      drhof(ndepth+1,ntheta+1) = drhof(ndepth+1,1)
      dHf(ndepth+1,ntheta+1) = dHf(ndepth+1,1)
      dfrad(:,ntheta+1) = dfrad(:,1)


      dmomr = dmom_r_2d (dppf,drhof,dHf,dutf,durf,ppf_2d,rhof_2d,Hf_2d,utf_2d,urf_2d)
      dmomt = dmom_t_2d (dppf,drhof,dutf,ddupf,durf,ppf_2d,rhof_2d,utf_2d,dupf_2d,urf_2d,cpf_2d)
      dener = denergy_2d(dppf,drhof,dttf,dutf,ddupf,durf,dfrad,ppf_2d,rhof_2d,ttf_2d,utf_2d,dupf_2d,urf_2d,frad_2d,fint,cpf_2d)
      dmass = dmass_2d(drhof,dppf,dutf,ddupf,durf,rhof_2d,ppf_2d,utf_2d,dupf_2d,urf_2d)
      dmomp = dmom_p_2d(dppf,drhof,dttf,dutf,ddupf,durf,dfrad,ppf_2d,rhof_2d,ttf_2d,utf_2d,dupf_2d,urf_2d,frad_2d,fint,cpf_2d)

      ijac = 2*(ndepth+1)*(ntheta+1)+ndepth+1+(it-1)*(ndepth+1)

      call write_jacobi(ijac,dmomr,dmomt,dener,dmomp,dmass)

      dttf = 0. ; drhof= 0. ; dfrad = 0. ; dHf = 0.


   enddo




 end subroutine fill_matrices_2d

 subroutine write_jacobi(ijac2,dmomr,dmomt,dener,dmomp,dmass)

   use mod_grid

   implicit none

   integer :: ijac1,ijac2,it,idepth
   real,dimension(ndepth+1,0:ntheta+2) :: dmomr,dmomt,dener,dmass,dmomp

   do it=1,ntheta
      do idepth = 1,ndepth

         ijac1 = idepth+(it-1)*(ndepth+1)
         jacobi(ijac1,ijac2) = dmomr(idepth,it)

         ijac1 = (ndepth+1)*(ntheta+1)+idepth+(it-1)*(ndepth+1)
         jacobi(ijac1,ijac2) = dmomt(idepth,it)

         ijac1 = 3*(ndepth+1)*(ntheta+1)+idepth+(it-1)*(ndepth+1)
         jacobi(ijac1,ijac2) = dmass(idepth,it)

         ijac1 = 4*(ndepth+1)*(ntheta+1)+idepth+(it-1)*(ndepth+1)
         jacobi(ijac1,ijac2) = dmomp(idepth,it)


      enddo
   enddo

   do it=1,ntheta
      do idepth = 1,ndepth+1

         ijac1 = 2*(ndepth+1)*(ntheta+1)+idepth+(it-1)*(ndepth+1)
         jacobi(ijac1,ijac2) = dener(idepth,it)

      enddo
   enddo


 end subroutine write_jacobi



 function mom_r_2d(ploc,rholoc,Hloc,utloc,urloc)

   use mod_grid

   implicit none

   real,dimension(ndepth+1,0:ntheta+2) :: ploc,utloc,Hloc, mom_r_2d,rholoc,urloc
   integer :: it,idepth,itp,itm,idp,idm
   real    :: fur,fut,fpp

   mom_r_2d = 0.

   do it = 1,ntheta
      do idepth=2,ndepth

         if (urloc(idepth,it) == 0.) then
            idp = idepth+1
            idm = idepth-1
         else if (urloc(idepth,it)>0.) then
            idp = idepth+1
            idm = idepth
         else
            idp = idepth
            idm = idepth-1
         endif

         if (utloc(idepth,it) == 0.) then
            itp = it+1
            itm = it-1
         else if (utloc(idepth,it)>0.) then
            itp = it
            itm = it-1
         else
            itp = it+1
            itm = it
         endif

         fur = urloc(idepth,it)*(urloc(idp,it)-urloc(idm,it))/(rrf(idp)-rrf(idm))*rrf(idepth)

         fut =  - utloc(idepth,it)**2 -2.*omega*utloc(idepth,it)*rrf(idepth) + utloc(idepth,it)*(urloc(idepth,itp)-urloc(idepth,itm))/(thetaf(itp)-thetaf(itm))

         fpp = (ploc(idepth+1,it)-ploc(idepth-1,it))/(rrf(idepth+1)-rrf(idepth-1))/rholoc(idepth,it)*rrf(idepth)+grf(idepth)*rrf(idepth)

         if(ploc(idepth,it)<pp2d_decoupling_momr) then
            fut = fut * (ploc(idepth,it)/pp2d_decoupling_momr)
            fur = fur * (ploc(idepth,it)/pp2d_decoupling_momr)
         endif

         mom_r_2d(idepth,it) = fur + fut + fpp

      enddo
   enddo

 end function mom_r_2d

 function dmom_r_2d(dploc,drholoc,dHloc,dutloc,durloc,ploc,rholoc,Hloc,utloc,urloc)

   use mod_grid

   implicit none

   real,dimension(ndepth+1,0:ntheta+2) :: dploc,dutloc,dHloc,dmom_r_2d,ploc,Hloc,utloc,drholoc,rholoc,durloc,urloc

   integer :: it,idepth,itp,itm,idp,idm
   real    :: dfur,dfut,dfpp

   dmom_r_2d = 0.

   do it = 1,ntheta
      do idepth=2,ndepth

         if (urloc(idepth,it) == 0.) then
            idp = idepth+1
            idm = idepth-1
         else if (urloc(idepth,it)>0.) then
            idp = idepth+1
            idm = idepth
         else
            idp = idepth
            idm = idepth-1
         endif

         if (utloc(idepth,it) == 0.) then
            itp = it+1
            itm = it-1
         else if (utloc(idepth,it)>0.) then
            itp = it
            itm = it-1
         else
            itp = it+1
            itm = it
         endif

         dfur = durloc(idepth,it)*(urloc(idp,it)-urloc(idm,it))/(rrf(idp)-rrf(idm))*rrf(idepth)+urloc(idepth,it)*(durloc(idp,it)-durloc(idm,it))/(rrf(idp)-rrf(idm))*rrf(idepth)

         dfut =  - 2.*utloc(idepth,it)*dutloc(idepth,it) -2.*omega*dutloc(idepth,it)*rrf(idepth) + dutloc(idepth,it)*(urloc(idepth,itp)-urloc(idepth,itm))/(thetaf(itp)-thetaf(itm))+utloc(idepth,it)*(durloc(idepth,itp)-durloc(idepth,itm))/(thetaf(itp)-thetaf(itm))

         dfpp = (dploc(idepth+1,it)-dploc(idepth-1,it))/(rrf(idepth+1)-rrf(idepth-1))/rholoc(idepth,it)*rrf(idepth)-(ploc(idepth+1,it)-ploc(idepth-1,it))/(rrf(idepth+1)-rrf(idepth-1))/rholoc(idepth,it)**2*drholoc(idepth,it)*rrf(idepth)

         if(ploc(idepth,it)<pp2d_decoupling_momr) then
            dfut = dfut * (ploc(idepth,it)/pp2d_decoupling_momr)
            dfur = dfur * (ploc(idepth,it)/pp2d_decoupling_momr)
         endif

         dmom_r_2d(idepth,it) = dfur + dfut + dfpp

      enddo
   enddo

 end function dmom_r_2d

  function mom_t_2d(ploc,rholoc,utloc,duploc,urloc,cploc)

   use mod_grid
   use mod_cst
   use mod_util

   implicit none

   real,dimension(ndepth+1,0:ntheta+2) :: ploc,utloc,rholoc,mom_t_2d,cploc,duploc,urloc
   integer :: it,idepth,itp,itm,idp,idm
   real    :: fur,fut,fpp,uloc

   mom_t_2d = 0.

   do it = 2,ntheta
      do idepth=2,ndepth

         if (urloc(idepth,it) == 0.) then
            idp = idepth+1
            idm = idepth-1
         else if (urloc(idepth,it)>0.) then
            idp = idepth+1
            idm = idepth
         else
            idp = idepth
            idm = idepth-1
         endif

         if (abs(utloc(idepth,it)) <1E-10) then
            itp = it+1
            itm = it-1
         else if (utloc(idepth,it)>0.) then
            itp = it
            itm = it-1
         else
            itp = it+1
            itm = it
         endif

        fur = urloc(idepth,it)*(rrf(idp)*utloc(idp,it)-rrf(idm)*utloc(idm,it))/(rrf(idp)-rrf(idm))+2.*omega*rrf(idepth)*urloc(idepth,it)

        fut = utloc(idepth,it)*(utloc(idepth,itp)-utloc(idepth,itm))/(thetaf(itp)-thetaf(itm))

        fpp = (ploc(idepth,it+1)-ploc(idepth,it-1))/(thetaf(it+1)-thetaf(it-1))/rholoc(idepth,it)

        if(ploc(idepth,it)<pp2d_decoupling_momt) then
           fpp = fpp * ploc(idepth,it)/pp2d_decoupling_momt
        endif

        mom_t_2d(idepth,it) = fur + fut + fpp

      enddo
    enddo


    do idepth=2,ndepth

       call interp_1d(log10(p_2din),u_2din,log10(ploc(idepth,1)),uloc)

       mom_t_2d(idepth,1) = utloc(idepth,1)**2-(uloc)**2

    enddo


 end function mom_t_2d


 function dmom_t_2d(dploc,drholoc,dutloc,dduploc,durloc,ploc,rholoc,utloc,duploc,urloc,cploc)

   use mod_grid
   use mod_cst

   implicit none

   real,dimension(ndepth+1,0:ntheta+2) :: ploc,utloc,rholoc,dmom_t_2d,dploc,drholoc,dutloc,cploc,duploc,dduploc,urloc,durloc
   integer :: it,idepth,itp,itm,idp,idm
   real    :: dfur,dfut,dfpp

   dmom_t_2d = 0.

   do it = 2,ntheta
      do idepth=2,ndepth

          if (urloc(idepth,it) == 0.) then
            idp = idepth+1
            idm = idepth-1
         else if (urloc(idepth,it)>0.) then
            idp = idepth+1
            idm = idepth
         else
            idp = idepth
            idm = idepth-1
         endif

         if (abs(utloc(idepth,it)) < 1E-10) then
            itp = it+1
            itm = it-1
         else if (utloc(idepth,it)>0.) then
            itp = it
            itm = it-1
         else
            itp = it+1
            itm = it
         endif

         dfur = durloc(idepth,it)*(rrf(idp)*utloc(idp,it)-rrf(idm)*utloc(idm,it))/(rrf(idp)-rrf(idm))+urloc(idepth,it)*(rrf(idp)*dutloc(idp,it)-rrf(idm)*dutloc(idm,it))/(rrf(idp)-rrf(idm))+2.*omega*rrf(idepth)*durloc(idepth,it)

         dfut = dutloc(idepth,it)*(utloc(idepth,itp)-utloc(idepth,itm))/(thetaf(itp)-thetaf(itm))+utloc(idepth,it)*(dutloc(idepth,itp)-dutloc(idepth,itm))/(thetaf(itp)-thetaf(itm))

        dfpp = (dploc(idepth,it+1)-dploc(idepth,it-1))/(thetaf(it+1)-thetaf(it-1))/rholoc(idepth,it) - (ploc(idepth,it+1)-ploc(idepth,it-1))/(thetaf(it+1)-thetaf(it-1))/rholoc(idepth,it)**2*drholoc(idepth,it)

        if(ploc(idepth,it)<pp2d_decoupling_momt) then
           dfpp = dfpp * ploc(idepth,it)/pp2d_decoupling_momt
        endif

         dmom_t_2d(idepth,it) = dfur + dfut + dfpp


      enddo
   enddo

   do idepth=2,ndepth

       dmom_t_2d(idepth,1) = 2.*utloc(idepth,1)*dutloc(idepth,1)

    enddo


 end function dmom_t_2d

 function energy_2d(ploc,rholoc,tloc,utloc,duploc,urloc,Floc,fint,cploc)

   use mod_grid
   use mod_cst

   implicit none

   real,dimension(ndepth+1,0:ntheta+2) :: ploc,utloc,rholoc,energy_2d,Floc,cploc,tloc,duploc,urloc,gam_loc
   integer :: it,idepth,itp,itm,idp,idm
   real    :: fur,fut,fup,ff,fint

   energy_2d = 0.

   do it = 1,ntheta
      do idepth = 2,ndepth

         if (urloc(idepth,it) == 0.) then
            idp = idepth+1
            idm = idepth-1
         else if (urloc(idepth,it)>0.) then
            idp = idepth+1
            idm = idepth
         else
            idp = idepth
            idm = idepth-1
         endif

         if (abs(utloc(idepth,it)) < 1E-10) then
            itp = it+1
            itm = it-1
         else if (utloc(idepth,it)>0.) then
            itp = it
            itm = it-1
         else
            itp = it+1
            itm = it
         endif

         fur = cploc(idepth,it)*ploc(idepth,it)*(rrf(idepth+1)**2*urloc(idepth+1,it)-rrf(idepth-1)**2*urloc(idepth-1,it))/(rrf(idepth+1)-rrf(idepth-1))/rrf(idepth)+(cploc(idepth,it)-1.)*urloc(idepth,it)*(ploc(idp,it)-ploc(idm,it))/(rrf(idp)-rrf(idm))*rrf(idepth)

         fup = cploc(idepth,it)*duploc(idepth,it)*ploc(idepth,it)

         fut = cploc(idepth,it)*ploc(idepth,it)*(utloc(idepth,it+1)-utloc(idepth,it-1))/(thetaf(it+1)-thetaf(it-1))+(cploc(idepth,it)-1.)*utloc(idepth,it)*(ploc(idepth,itp)-ploc(idepth,itm))/(thetaf(itp)-thetaf(itm))

         ff = (Floc(idepth+1,it)-Floc(idepth,it))/(rrf(idepth+1)-rrf(idepth))*rrf(idepth)

         energy_2d(idepth,it) = fur + fup + fut + ff

      enddo

      energy_2d(ndepth+1,it) = fint - Floc(ndepth+1,it)

   enddo

 end function energy_2d

 function denergy_2d(dploc,drholoc,dtloc,dutloc,dduploc,durloc,dFloc,ploc,rholoc,tloc,utloc,duploc,urloc,Floc,fint,cploc)

   use mod_grid
   use mod_cst

   implicit none

   real,dimension(ndepth+1,0:ntheta+2) :: ploc,utloc,rholoc,denergy_2d,dFloc,cploc,dploc,drholoc,dutloc,tloc,dtloc,dduploc,duploc,urloc,durloc,floc
   integer :: it,idepth,itp,itm,idp,idm
   real    :: dfur,dfut,dfup,dff,fint


   denergy_2d = 0.

   do it = 1,ntheta
      do idepth=2,ndepth

         if (urloc(idepth,it) == 0.) then
            idp = idepth+1
            idm = idepth-1
         else if (urloc(idepth,it)>0.) then
            idp = idepth+1
            idm = idepth
         else
            idp = idepth
            idm = idepth-1
         endif

         if (abs(utloc(idepth,it)) < 1E-10) then
            itp = it+1
            itm = it-1
         else if (utloc(idepth,it)>0.) then
            itp = it
            itm = it-1
         else
            itp = it+1
            itm = it
         endif

         dfur = cploc(idepth,it)*dploc(idepth,it)*(rrf(idepth+1)**2*urloc(idepth+1,it)-rrf(idepth-1)**2*urloc(idepth-1,it))/(rrf(idepth+1)-rrf(idepth-1))/rrf(idepth)+(cploc(idepth,it)-1.)*durloc(idepth,it)*(ploc(idp,it)-ploc(idm,it))/(rrf(idp)-rrf(idm))*rrf(idepth)+cploc(idepth,it)*ploc(idepth,it)*(rrf(idepth+1)**2*durloc(idepth+1,it)-rrf(idepth-1)**2*durloc(idepth-1,it))/(rrf(idepth+1)-rrf(idepth-1))/rrf(idepth)+(cploc(idepth,it)-1.)*urloc(idepth,it)*(dploc(idp,it)-dploc(idm,it))/(rrf(idp)-rrf(idm))*rrf(idepth)

         dfup = cploc(idepth,it)*dduploc(idepth,it)*ploc(idepth,it) + cploc(idepth,it)*duploc(idepth,it)*dploc(idepth,it)

         dfut = cploc(idepth,it)*dploc(idepth,it)*(utloc(idepth,it+1)-utloc(idepth,it-1))/(thetaf(it+1)-thetaf(it-1))+(cploc(idepth,it)-1.)*dutloc(idepth,it)*(ploc(idepth,itp)-ploc(idepth,itm))/(thetaf(itp)-thetaf(itm)) + cploc(idepth,it)*ploc(idepth,it)*(dutloc(idepth,it+1)-dutloc(idepth,it-1))/(thetaf(it+1)-thetaf(it-1))+(cploc(idepth,it)-1.)*utloc(idepth,it)*(dploc(idepth,itp)-dploc(idepth,itm))/(thetaf(itp)-thetaf(itm))

         dff = (dFloc(idepth+1,it)-dFloc(idepth,it))/(rrf(idepth+1)-rrf(idepth))*rrf(idepth)

         denergy_2d(idepth,it) = dfur + dfup + dfut + dff

      enddo

       denergy_2d(ndepth+1,it) = -dFloc(ndepth+1,it)

    enddo

 end function denergy_2d

 function mass_2d(rholoc,ploc,utloc,duploc,urloc)

   use mod_grid

   implicit none

   real,dimension(ndepth+1,0:ntheta+2) :: rholoc,utloc,duploc,mass_2d,ploc,urloc
   integer :: it,idepth,itp,itm,idp,idm
   real    :: fur,frr,fut,frt,fup

   mass_2d = 0.

   do it = 1,ntheta
      do idepth=2,ndepth

         if (urloc(idepth,it) == 0.) then
            idp = idepth+1
            idm = idepth-1
         else if (urloc(idepth,it)>0.) then
            idp = idepth+1
            idm = idepth
         else
            idp = idepth
            idm = idepth-1
         endif

         if (utloc(idepth,it) == 0.) then
            itp = it+1
            itm = it-1
         else if (utloc(idepth,it)>0.) then
            itp = it
            itm = it-1
         else
            itp = it+1
            itm = it
         endif

         fur  = (rrf(idepth+1)**2*urloc(idepth+1,it)-rrf(idepth-1)**2*urloc(idepth-1,it))/(rrf(idepth+1)-rrf(idepth-1))/rrf(idepth)

         frr = rrf(idepth)*urloc(idepth,it)*(rholoc(idp,it)-rholoc(idm,it))/(rrf(idp)-rrf(idm))/rholoc(idepth,it)

         fup  = duploc(idepth,it)

         fut  = (utloc(idepth,it+1)-utloc(idepth,it-1))/(thetaf(it+1)-thetaf(it-1))

         frt  = utloc(idepth,it)*(rholoc(idepth,itp)-rholoc(idepth,itm))/(thetaf(itp)-thetaf(itm))/rholoc(idepth,it)

         mass_2d(idepth,it) = fur + frr + fup + fut + frt

      enddo
   enddo

 end function mass_2d

 function dmass_2d(drholoc,dploc,dutloc,dduploc,durloc,rholoc,ploc,utloc,duploc,urloc)

   use mod_grid

   implicit none

   real,dimension(ndepth+1,0:ntheta+2) :: rholoc,utloc,duploc,dmass_2d,drholoc,dutloc,dduploc,ploc,dploc,urloc,durloc
   integer :: it,idepth,itp,itm,idp,idm
   real    :: dfur,dfrr,dfup,dfut,dfrt

   dmass_2d = 0.

   do it = 1,ntheta
      do idepth=2,ndepth

         if (urloc(idepth,it) == 0.) then
            idp = idepth+1
            idm = idepth-1
         else if (urloc(idepth,it)>0.) then
            idp = idepth+1
            idm = idepth
         else
            idp = idepth
            idm = idepth-1
         endif

         if (utloc(idepth,it) == 0.) then
            itp = it+1
            itm = it-1
         else if (utloc(idepth,it)>0.) then
            itp = it
            itm = it-1
         else
            itp = it+1
            itm = it
         endif

         dfur  = (rrf(idepth+1)**2*durloc(idepth+1,it)-rrf(idepth-1)**2*durloc(idepth-1,it))/(rrf(idepth+1)-rrf(idepth-1))/rrf(idepth)

         dfrr = rrf(idepth)*durloc(idepth,it)*(rholoc(idp,it)-rholoc(idm,it))/(rrf(idp)-rrf(idm))/rholoc(idepth,it) + rrf(idepth)*urloc(idepth,it)*(drholoc(idp,it)-drholoc(idm,it))/(rrf(idp)-rrf(idm))/rholoc(idepth,it) - rrf(idepth)*urloc(idepth,it)*(rholoc(idp,it)-rholoc(idm,it))/(rrf(idp)-rrf(idm))/rholoc(idepth,it)**2*drholoc(idepth,it)

         dfup  = dduploc(idepth,it)

         dfut  = (dutloc(idepth,it+1)-dutloc(idepth,it-1))/(thetaf(it+1)-thetaf(it-1))

         dfrt  = dutloc(idepth,it)*(rholoc(idepth,itp)-rholoc(idepth,itm))/(thetaf(itp)-thetaf(itm))/rholoc(idepth,it) + utloc(idepth,it)*(drholoc(idepth,itp)-drholoc(idepth,itm))/(thetaf(itp)-thetaf(itm))/rholoc(idepth,it) - utloc(idepth,it)*(rholoc(idepth,itp)-rholoc(idepth,itm))/(thetaf(itp)-thetaf(itm))/rholoc(idepth,it)**2*drholoc(idepth,it)

         dmass_2d(idepth,it) = dfur + dfrr + dfup + dfut + dfrt

      enddo
   enddo

 end function dmass_2d

 function mom_p_2d(ploc,rholoc,tloc,utloc,duploc,urloc,Floc,fint,cploc)

   use mod_grid
   use mod_cst
   use mod_util

   implicit none

   real,dimension(ndepth+1,0:ntheta+2) :: ploc,utloc,rholoc,mom_p_2d,Floc,cploc,tloc,duploc,urloc,gam_loc
   integer :: it,idepth,itp,itm,idp,idm
   real    :: fur,fut,fup,ff,fint,frt,frr,alpha_loc

   mom_p_2d = 0.

   do it = 1,ntheta
      do idepth=2,ndepth


         if (urloc(idepth,it) == 0.) then
            idp = idepth+1
            idm = idepth-1
         else if (urloc(idepth,it)>0.) then
            idp = idepth+1
            idm = idepth
         else
            idp = idepth
            idm = idepth-1
         endif

         if (utloc(idepth,it) == 0.) then
            itp = it+1
            itm = it-1
         else if (utloc(idepth,it)>0.) then
            itp = it
            itm = it-1
         else
            itp = it+1
            itm = it
         endif

         fur  = (rrf(idepth+1)**2*urloc(idepth+1,it)-rrf(idepth-1)**2*urloc(idepth-1,it))/(rrf(idepth+1)-rrf(idepth-1))/rrf(idepth)

         frr = rrf(idepth)*urloc(idepth,it)*(rholoc(idp,it)-rholoc(idm,it))/(rrf(idp)-rrf(idm))/rholoc(idepth,it)

         fup  = duploc(idepth,it)

         fut  = (utloc(idepth,it+1)-utloc(idepth,it-1))/(thetaf(it+1)-thetaf(it-1))

         frt  = utloc(idepth,it)*(rholoc(idepth,itp)-rholoc(idepth,itm))/(thetaf(itp)-thetaf(itm))/rholoc(idepth,it)



         if (cos(thetaf(it))>0.) then
            call interp_1d(log10(p_2din),alpha_dayside_2din,log10(ploc(idepth,it)),alpha_loc)
         else
            call interp_1d(log10(p_2din),alpha_nightside_2din,log10(ploc(idepth,it)),alpha_loc)
         endif

         mom_p_2d(idepth,it) = fur + frr - 1./alpha_loc*fup

      enddo
   enddo

 end function mom_p_2d

 function dmom_p_2d(dploc,drholoc,dtloc,dutloc,dduploc,durloc,dFloc,ploc,rholoc,tloc,utloc,duploc,urloc,Floc,fint,cploc)

   use mod_grid
   use mod_cst
   use mod_util

   implicit none

   real,dimension(ndepth+1,0:ntheta+2) :: ploc,utloc,rholoc,dmom_p_2d,dFloc,cploc,dploc,drholoc,dutloc,tloc,dtloc,dduploc,duploc,urloc,durloc,floc
   integer :: it,idepth,itp,itm,idp,idm
   real    :: dfur,dfut,dfup,dff,fint,dfrt,dfrr,alpha_loc

   dmom_p_2d = 0.

   do it = 1,ntheta
      do idepth=2,ndepth

          if (urloc(idepth,it) == 0.) then
            idp = idepth+1
            idm = idepth-1
         else if (urloc(idepth,it)>0.) then
            idp = idepth+1
            idm = idepth
         else
            idp = idepth
            idm = idepth-1
         endif

         if (utloc(idepth,it) == 0.) then
            itp = it+1
            itm = it-1
         else if (utloc(idepth,it)>0.) then
            itp = it
            itm = it-1
         else
            itp = it+1
            itm = it
         endif

         dfur  = (rrf(idepth+1)**2*durloc(idepth+1,it)-rrf(idepth-1)**2*durloc(idepth-1,it))/(rrf(idepth+1)-rrf(idepth-1))/rrf(idepth)

         dfrr = rrf(idepth)*durloc(idepth,it)*(rholoc(idp,it)-rholoc(idm,it))/(rrf(idp)-rrf(idm))/rholoc(idepth,it) + rrf(idepth)*urloc(idepth,it)*(drholoc(idp,it)-drholoc(idm,it))/(rrf(idp)-rrf(idm))/rholoc(idepth,it) - rrf(idepth)*urloc(idepth,it)*(rholoc(idp,it)-rholoc(idm,it))/(rrf(idp)-rrf(idm))/rholoc(idepth,it)**2*drholoc(idepth,it)

         dfup  = dduploc(idepth,it)

         dfut  = (dutloc(idepth,it+1)-dutloc(idepth,it-1))/(thetaf(it+1)-thetaf(it-1))

         dfrt  = dutloc(idepth,it)*(rholoc(idepth,itp)-rholoc(idepth,itm))/(thetaf(itp)-thetaf(itm))/rholoc(idepth,it) + utloc(idepth,it)*(drholoc(idepth,itp)-drholoc(idepth,itm))/(thetaf(itp)-thetaf(itm))/rholoc(idepth,it) - utloc(idepth,it)*(rholoc(idepth,itp)-rholoc(idepth,itm))/(thetaf(itp)-thetaf(itm))/rholoc(idepth,it)**2*drholoc(idepth,it)

         if (cos(thetaf(it))>0.) then
            call interp_1d(log10(p_2din),alpha_dayside_2din,log10(ploc(idepth,it)),alpha_loc)
         else
            call interp_1d(log10(p_2din),alpha_nightside_2din,log10(ploc(idepth,it)),alpha_loc)
         endif

         dmom_p_2d(idepth,it) = dfur + dfrr - 1./alpha_loc*dfup

      enddo
   enddo

 end function dmom_p_2d

   !> Find the pressure/temp profile using a newton-raphson method in 2d
    !<
    subroutine find_profile_2d

      use mod_grid
      use mod_param
      use mod_opacity
      use mod_cst
      use mod_chem
      use mod_util
      use mod_eos
      use mod_radtrans

      implicit none
      include 'netcdf.inc'

      integer:: i,nit,idepth,it,idptu,nlevel,ncol,id_file,id_var,ndepth_um
      real,dimension(5*(ndepth+1)*(ntheta+1)) :: dptu
      real,dimension(ndepth+1,0:ntheta+2) :: dpp,dtt,dut,ppf_new,ddup,dur
      character(4)  :: fnit
      real :: err_hydro,err_energy,err,errold,errveryold,epsilon,coseps,coseps2,eta



      step = minstep
      err        = 1.
      errold     = 1.
      errveryold = 1.
      nit  = 0

      do while (err>accuracy)

         jacobi = 0.
         rhs    = 0.

         do i = 1,5*(ndepth+1)*(ntheta+1)

            jacobi(i,i) = 1.

         enddo

         err = 1.

         call fill_matrices_2d(err)

         if (err>accuracy) then

            nit = nit + 1

            ! Set correction step size
            if(nit .gt. 1 )then
               ! Increase step size if residuum falls monotonically
               if( (err.lt.errold) .and. (errold.lt.errveryold) )then
                  step=min(step*1.1,maxstep)
               endif
               ! Decrease step size if residuum increases
               if(err .gt. errold )then
                  step=max(step*0.5,minstep)
               endif
            else
               ! Set a small step size for the beginning
               step=minstep
            endif

            ! Hand over maximum
            errveryold = errold
            errold     = err


            do idptu = 1, 5*(ndepth+1)*(ntheta+1)

               if ((sum(abs(jacobi(idptu,:))))==0.) jacobi(idptu,idptu) = 1.

            enddo

            write(fnit,'(i4)') nit

            dptu = 0. ; dpp = 0. ; dtt = 0. ; dut = 0. ; ddup = 0. ; dur = 0.
            if (debug>0 .and. mype==cputerm) write(*,*) 'start linsolve'
            call linsolve(5*(ndepth+1)*(ntheta+1),jacobi,rhs,dptu)
            if (debug>0 .and. mype==cputerm) write(*,*) 'end linsolve'

            do it=1,ntheta
               do idepth = 1,ndepth+1

                  idptu = idepth+(it-1)*(ndepth+1)
                  dpp(idepth,it) = dptu(idptu)

                  idptu = (ndepth+1)*(ntheta+1)+idepth+(it-1)*(ndepth+1)
                  dut(idepth,it) = dptu(idptu)

                  idptu = 2*(ndepth+1)*(ntheta+1)+idepth+(it-1)*(ndepth+1)
                  dtt(idepth,it) = dptu(idptu)

                  idptu = 3*(ndepth+1)*(ntheta+1)+idepth+(it-1)*(ndepth+1)
                  ddup(idepth,it) = dptu(idptu)

                  idptu = 4*(ndepth+1)*(ntheta+1)+idepth+(it-1)*(ndepth+1)
                  dur(idepth,it) = dptu(idptu)


               enddo
            enddo

            if (mype == cputerm) then
               write(*,headrfmt) 'PT Iter.: '//fnit
               write(*,1010) maxval(abs(dpp(1:ndepth+1,1:ntheta)/ppf_2d(1:ndepth+1,1:ntheta))),step
               write(*,1011) maxval(abs(dtt(1:ndepth+1,1:ntheta)/ttf_2d(1:ndepth+1,1:ntheta))),err
               write(*,1012) maxval(abs(dut(1:ndepth+1,1:ntheta))),maxval(abs(utf_2d(1:ndepth+1,1:ntheta)))
               write(*,1013) maxval(abs(ddup(1:ndepth+1,1:ntheta))),maxval(abs(dupf_2d(1:ndepth+1,1:ntheta)))
               write(*,1014) maxval(abs(dur(1:ndepth+1,1:ntheta))),maxval(abs(urf_2d(1:ndepth+1,1:ntheta)))

               write(*,*)
            endif

1010        format(5x,' dp/p: ',es10.1,5x,' step : ',es10.1)
1011        format(5x,' dt/t: ',es10.1,5x,' err. : ',es10.2)
1012        format(5x,'  dut: ',es10.1,5x,' ut   : ',es10.2)
1013        format(5x,' ddup: ',es10.1,5x,' dup  : ',es10.2)
1014        format(5x,'  dur: ',es10.1,5x,' ur   : ',es10.2)

            ppf_2d(1:ndepth+1,1:ntheta ) = ppf_2d(1:ndepth+1,1:ntheta ) + step*dpp(1:ndepth+1,1:ntheta)
            ttf_2d(1:ndepth+1,1:ntheta ) = ttf_2d(1:ndepth+1,1:ntheta ) + step*dtt(1:ndepth+1,1:ntheta)
            utf_2d(1:ndepth+1,1:ntheta ) = utf_2d(1:ndepth+1,1:ntheta ) + step*dut(1:ndepth+1,1:ntheta)
            dupf_2d(1:ndepth+1,1:ntheta) = dupf_2d(1:ndepth+1,1:ntheta) + step*ddup(1:ndepth+1,1:ntheta)
            urf_2d(1:ndepth+1,1:ntheta ) = urf_2d(1:ndepth+1,1:ntheta ) + step*dur(1:ndepth+1,1:ntheta)


            urf_2d (ndepth+1,:) = 0.
            utf_2d (ndepth+1,:) = 0.
            dupf_2d(ndepth+1,:) = 0.
            ppf_2d (ndepth+1,:) = pmax

            urf_2d (1,:) = 2.*urf_2d (2,:) - urf_2d (3,:)
            utf_2d (1,:) = 2.*utf_2d (2,:) - utf_2d (3,:)
            dupf_2d(1,:) = 2.*dupf_2d(2,:) - dupf_2d(3,:)
            ttf_2d (1,:) = 2.*ttf_2d (2,:) - ttf_2d (3,:)
            ppf_2d (1,:) = ppf_2d(2,:)**2  / ppf_2d (3,:)


            ! do idepth=1,ndepth
            !    write(*,*) urf_2d(idepth,1),urf_2d(idepth,3),urf_2d(idepth,5),urf_2d(idepth,7),urf_2d(idepth,9),rrf(idepth)/(rrf(idepth+1)-rrf(idepth))

            ! enddo



         else


            if (mype == cputerm) then
               write(fnit,'(i4)') (nit+1)
               write(*,headrfmt) 'Iteration: '//fnit
               write(*,1015) err
               write(*,*)
               write(*,separfmt)
            endif
1015        format('err. : ',es10.1)


         endif

         !check that pressure temperature is positive
         if (minval(ttf_2d(1:ndepth+1,1:ntheta+1))<0. .or. minval(ppf_2d(1:ndepth+1,1:ntheta+1))<0.) then

            if(mype==cputerm) then

               write(*,headrfmt) 'Error find_profile'
               write(*,'(A)') 'Reached a negative pressure or temperature'
               write(*,exppmfmt) 'min pressure',minval(ppf_2d(1:ndepth+1,1:ntheta+1))
               write(*,exppmfmt) 'min temperature',minval(ttf_2d(1:ndepth+1,1:ntheta+1))
               write(*,*)



               stop

            endif
         endif


            !save profiles
            if(mype == cputerm) then
               if (debug>0) write(*,*) 'start write_pt_profile'
               call write_pt_profile(ftemp)
               if (debug>0) write(*,*) 'end write_pt_profile'
            endif


         enddo

         !save profiles
         if (err<accuracy) then
         if(mype == cputerm) then
            if (debug>0) write(*,*) 'start write_pt_profile'
            call write_pt_profile(fout)
            if (debug>0) write(*,*) 'end write_pt_profile'
         endif
      endif


    end subroutine find_profile_2d

!! Calling and writing contribution function ro netcdf file

    subroutine compute_con_fun

     use mod_radtrans
     use mod_util
     use mod_param
     use mod_cst
     use mod_opacity
     use mod_grid
     use mod_chem

    implicit none
    include 'netcdf.inc'

    integer :: id_file1,id_var1,id_var2,id_var3,id_var4,id_var5,id_var6,id_var7,id_var8,id_var9,id_lev,id_lev_int,id_freq, istart, icount, idepth, iband
    real,dimension(:,:),allocatable :: bnu, fnu_htr, qnu_htr, dfnu_htr_drrf
    real,dimension(:),allocatable :: nu_cf
    real,dimension(ndepth) :: dfrad_drrf, cp_chem_center


     if (.not. corr_k) then
      allocate(bnu(ndepth+1,nfreqkap))
      allocate(fnu_htr(ndepth+1,nfreqkap))
      allocate(qnu_htr(ndepth+1,nfreqkap))
      allocate(dfnu_htr_drrf(ndepth+1,nfreqkap))
      allocate(nu_cf(nfreqkap))
      nu_cf = nu(:)
     else
      allocate(bnu(ndepth+1,nband))
      allocate(fnu_htr(ndepth+1,nband))
      allocate(qnu_htr(ndepth+1,nband))
      allocate(dfnu_htr_drrf(ndepth+1,nband))
      allocate(nu_cf(nband))
      do iband = 1,nband
      nu_cf(iband) = 0.5*(nubandmin(iband)+nubandmax(iband))
      enddo
     endif

     if (.not. corr_k) then
      call calc_con_fun(nfreq,taub,bnu)
     else
      call calc_con_fun(nband,taub,bnu)
     endif



   if (.not. corr_k) then
    istart = 1
    icount = nfreqkap
   else
    istart = 1
    icount = nband
   endif

   if (.not. corr_k) then
      fnu_htr = fnu(:,:)
      qnu_htr = fnu(:,:)

   else

    do idepth =1,ndepth+1
      do iband = 1,nband
          fnu_htr(idepth,iband) = sum(fnu(idepth,(iband-1)*nkmix+1:iband*nkmix)*weightsfnu(1,(iband-1)*nkmix+1:iband*nkmix))/(nubandmax(iband)-nubandmin(iband))
          qnu_htr(idepth,iband) = sum(qnu(idepth,(iband-1)*nkmix+1:iband*nkmix)*weightsfnu(1,(iband-1)*nkmix+1:iband*nkmix))/(nubandmax(iband)-nubandmin(iband))
      enddo
    enddo
   endif

        !interpolate to get cp_chem at cell center, using geometric mean
          do idepth = 1,ndepth
             cp_chem_center(idepth) = sqrt(cp_chem(idepth)*cp_chem(idepth+1))
          end do

   !! Compute_Heating rates

      do idepth = 1,ndepth
        dfnu_htr_drrf(idepth,:) = ((fnu_htr(idepth,:) - fnu_htr(idepth+1,:))/(rrf(idepth) - rrf(idepth+1)))   !! This gives spectral heating rate
        dfrad_drrf(idepth) = ((frad(idepth) - frad(idepth+1))/(rrf(idepth) - rrf(idepth+1)))              !! This gives integrated heating rate
        ! dfrad_drrf(idepth) = ((frad(idepth) - frad(idepth+1))/(rrf(idepth) - rrf(idepth+1)))*(1/(cp_chem_center(idepth)*rho(idepth)))  !! If needed in K/S integrated heating rate
      end do



      !! if (mype==cputerm) then
         call nf(nf_create(fcfout,nf_clobber,id_file1))
         call nf(nf_def_dim(id_file1,'nlevel',ndepth,id_lev))
         call nf(nf_def_dim(id_file1,'nlevel_interface',ndepth+1,id_lev_int))

         if (.not. corr_k) then
         call nf(nf_def_dim(id_file1,'nfreq',nfreqkap,id_freq))
         else
          call nf(nf_def_dim(id_file1,'nband',nband,id_freq))
         endif





            call nf(nf_def_var(id_file1,'nu',nf_double,1,id_freq,id_var1))
            call nf(nf_put_att_text(id_file1,id_var1,'units',4,'cm-1'))

         if (corr_k) then

            call nf(nf_def_var(id_file1,'nubandmin',nf_double,1,id_freq,id_var1))
            call nf(nf_put_att_text(id_file1,id_var1,'units',4,'cm-1'))

            call nf(nf_def_var(id_file1,'nubandmax',nf_double,1,id_freq,id_var1))
            call nf(nf_put_att_text(id_file1,id_var1,'units',4,'cm-1'))

          end if

          call nf(nf_def_var(id_file1,'ncfnu',nf_double,2,(/id_lev,id_freq/),id_var2)) !!! normalised contribution function
          call nf(nf_put_att_text(id_file1,id_var2,'units',4,'none'))

          call nf(nf_def_var(id_file1,'pressure',nf_double,1,id_lev,id_var3))
          call nf(nf_put_att_text(id_file1,id_var3,'units',9,'dyne cm-2'))

          call nf(nf_def_var(id_file1,'pressure_interface',nf_double,1,id_lev_int,id_var4))
          call nf(nf_put_att_text(id_file1,id_var4,'units',9,'dyne cm-2'))

          call nf(nf_def_var(id_file1,'spectral_htr',nf_double,2,(/id_lev,id_freq/),id_var5)) !!! Wavelength dependent heating rate
          call nf(nf_put_att_text(id_file1,id_var5,'units',12,'erg s-1 cm-3'))

          call nf(nf_def_var(id_file1,'int_htr',nf_double,1,id_lev,id_var6)) !!! Integrated heating rate
          call nf(nf_put_att_text(id_file1,id_var6,'units',12,'erg s-1 cm-3'))

          call nf(nf_def_var(id_file1,'spectral_depth_fnu',nf_double,2,(/id_lev_int,id_freq/),id_var7)) !!! Wavelength and height dependent flux
          call nf(nf_put_att_text(id_file1,id_var7,'units',15,'erg s-1 cm-2 cm'))

          call nf(nf_def_var(id_file1,'rrf',nf_double,1,id_lev_int,id_var8))
          call nf(nf_put_att_text(id_file1,id_var8,'units',2,'cm'))

          call nf(nf_close(id_file1))

       !! endif

!#if PARA == 1
!      call MPI_BARRIER(MPI_COMM_WORLD,info)
!#endif

!      call nf(nf_open_par(fspectrum,nf_write,MPI_COMM_WORLD,MPI_INFO_NULL,id_file))
      !! iwrite = 0
      !! if(mype==0) iwrite = 1

!#if PARA==1
!      if (mype>0) call MPI_RECV(iwrite,1,MPI_INTEGER,mype-1,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,info)
!#endif

      !! if (iwrite==1) then

            call nf(nf_open(fcfout,nf_write,id_file1))

            if (.not.corr_k) then
             call nf(nf_inq_varid(id_file1,'nu' ,id_var1))
             call nf(nf_put_vara_double(id_file1,id_var1,1,nfreqkap,nu_cf*(1/unitL)))
            else

             call nf(nf_inq_varid(id_file1,'nu' ,id_var1))
             call nf(nf_put_vara_double(id_file1,id_var1,1,nband,nu_cf*(1/unitL)))

              call nf(nf_inq_varid(id_file1,'nubandmin' ,id_var1))
              call nf(nf_put_vara_double(id_file1,id_var1,1,nband, &
                  nubandmin(1:nband)*(1/unitL)))

              call nf(nf_inq_varid(id_file1,'nubandmax' ,id_var1))
              call nf(nf_put_vara_double(id_file1,id_var1,1,nband, &
                  nubandmax(1:nband)*(1/unitL)))

            end if

            call nf(nf_inq_varid(id_file1,'ncfnu' ,id_var2))
            call nf(nf_put_vara_double(id_file1,id_var2,(/1,1/),(/ndepth,icount/),ncf(:,1:icount)))

            call nf(nf_inq_varid(id_file1,'pressure',id_var3))
            call nf(nf_put_vara_double(id_file1,id_var3,1,ndepth,pp*unitP))

            call nf(nf_inq_varid(id_file1,'pressure_interface',id_var4))
            call nf(nf_put_vara_double(id_file1,id_var4,1,ndepth+1,ppf*unitP))

            call nf(nf_inq_varid(id_file1,'spectral_htr' ,id_var5))
            call nf(nf_put_vara_double(id_file1,id_var5,(/1,1/),(/ndepth,icount/),dfnu_htr_drrf(:,1:icount)*(unitE/(unitT*(unitL**3)))))

            call nf(nf_inq_varid(id_file1,'int_htr' ,id_var6))
            call nf(nf_put_vara_double(id_file1,id_var6,1,ndepth,dfrad_drrf*(unitE/(unitT*(unitL**3)))))

            call nf(nf_inq_varid(id_file1,'spectral_depth_fnu' ,id_var7))
            call nf(nf_put_vara_double(id_file1,id_var7,(/1,1/),(/ndepth+1,icount/),fnu_htr(:,1:icount)*unitSF))

            call nf(nf_inq_varid(id_file1,'rrf',id_var8))
            call nf(nf_put_vara_double(id_file1,id_var8,1,ndepth+1,rrf*unitL))



         call nf(nf_close(id_file1))

      !! endif

    end subroutine compute_con_fun

  !--------------------------------------------
  !Initialise the LSODE solver
  !--------------------------------------------
  SUBROUTINE lsode_init(ndim)

  IMPLICIT NONE

  ! Number of ODEs
  INTEGER,INTENT(IN) :: ndim

  ! Length of integer and real work arrays
  liw = 30 + 2*ndim
  lrw = 20 + 11*ndim + 2*ndim**2

  ! Allocate integer and real work arrays
  ALLOCATE(iwork(liw))
  ALLOCATE(rwork(lrw))

  ! Initial values of solver parameters
  lsode_istate   = 1
  lsode_itask    = 1
  lsode_itol     = 1
  lsode_iopt     = 1

  iwork    = 0
  rwork    = 0.

  iwork(5) = lsode_order
  iwork(6) = lsode_nit
  iwork(7) = lsode_print_warning

  END SUBROUTINE lsode_init

  !--------------------------------------------
  !Call DLSODES solver
  !--------------------------------------------
  SUBROUTINE solve_dlsodes(ndim,y,t,t_plus_dt,f,jac)

  IMPLICIT NONE

  EXTERNAL f, jac

  ! Number of ODEs
  INTEGER,INTENT(IN) :: ndim
  ! Initial time
  REAL,INTENT(INOUT)    :: t
  ! Requested time
  REAL,INTENT(IN)    :: t_plus_dt
  ! Initial values
  REAL,INTENT(INOUT) :: y(ndim)

  ! Call DLSODES
  ! y are the initial values
  ! and returned values are also stored in y
  CALL dlsodes(                                                              &
          f,ndim, y, t, t_plus_dt, lsode_itol, lsode_rtol, lsode_atol,       &
          lsode_itask, lsode_istate, lsode_iopt, rwork, lrw,    &
          iwork, liw, jac, lsode_imf)

  END SUBROUTINE solve_dlsodes

end module mod_solver
