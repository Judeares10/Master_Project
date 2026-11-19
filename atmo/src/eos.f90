!> \file
!! Contains module mod_eos, subroutine mod_eos::init_eos(), mod_eos::get_eos()\n
!! Initialization of the equation of state, ideal gas per default, the user can change this routine for other eos
!<

!> hold the parameters of the equation of state
!<
module mod_eos

  real :: gamma       !< Adiabatic exponent of a monatomic ideal gas
  character(len=30) :: feostab !< EOS table name
  real :: mucst  !< if 0. compute mean molecular weight from the abundances
  real :: pmin_conv !< indicate the pressure for convection and change of gamma
  real :: pmax_conv !< indicate the pressure for convection and change of gamma
  logical :: gradad_smooth    !< if gradad is smoothed
  integer :: kereos_smooth !< number of neighbors used for smoothing (1 means smoothing on 3 cells)
  integer :: nb_eos_smooth
  integer :: nit_gradad_smooth !< Turn on gradad_smooth after nit_gradad_smooth
                               !< iterations of find_profile solver
  real(kind=8) :: xHe !< Effective mass fraction of Helium


  integer, parameter :: ntt=101, npp=321, ntab=2, nq=6
  character(len=100) :: eos_dir !< directory for eos tables
  character(len=100) :: eos_file_Y25 !< Filename for helium mass fraction 0.25
  character(len=100) :: eos_file_Y30 !< Filename for helium mass fraction 0.30
  character(6) :: eos_type !< if ='ideal' then ideal gas if ='table' use feostab

  ! Global variables
  real(kind=8), dimension(ntab)      , save, private :: yt
  real(kind=8), dimension(ntt)       , save, private :: tlgt
  real(kind=8), dimension(npp)      , save, private :: plgt
  real(kind=8), dimension(6,npp,ntt), save, private :: data

  logical :: ocean          !< logical to include oceann layer or not
  real    :: opacity_ocean  !< Ocea absorption
  real    :: ppmin_ocean    !< minimum pressure above which ocean opacity is used
  real    :: ppmax_ocean    !< maximum pressure above which ocean opacity is used
  real    :: density_ocean  !< density of the ocean
  real    :: gamma_ocean    !< adiabatic index of the ocean
  real    :: mu_ocean       !< mean molecular weight ocean
  real    :: unitPocean
  real    :: unitKocean
  real    :: unitDocean

  integer, save, private :: np, nt

contains

  subroutine init_eos

    use mod_param
    use mod_cst
    use mod_grid
    use mod_chem

    implicit none

    real(kind=8), dimension(ntab,npp,ntt) :: rolgt,elgt,slgt,cplgt
    real(kind=8), dimension(ntab,npp,ntt) :: abadtab,gamtab
    integer :: n,iunit,ip,jt,imin,imax,i,stat
    real(kind=8)    :: alpha, xHe
    real(kind=8)    :: tiny

    namelist /eos/gamma,feostab,eos_type,eos_dir,eos_file_Y25,eos_file_Y30,xHe,mucst,pmin_conv,Gradad_smooth,kereos_smooth,pmax_conv,nb_eos_smooth,nit_gradad_smooth,opacity_ocean,unitKocean,ppmax_ocean,ppmin_ocean,ocean,unitPocean,density_ocean,unitDocean, gamma_ocean, mu_ocean

    ! Set standard values
    gradad_smooth = .false.
    kereos_smooth = 1
    gamma     = 5.0/3.0
    feostab   = 'None'
    eos_type  = 'ideal'
    eos_dir   = '../../eos/'
    eos_file_Y25 = 'eostpY25_atmo'
    eos_file_Y30 = 'eostpY30_atmo'
    xHe       = 0.28566
    mucst     = 0.
    pmax_conv = 1E-40
    pmin_conv = 1E40
    nb_eos_smooth = 1
    nit_gradad_smooth = -1

    ocean = .false.
    opacity_ocean = 1.
    ppmin_ocean = 1E-10
    ppmax_ocean = 1E6
    density_ocean = 1.
    gamma_ocean = 1.0001
    mu_ocean = 20.
    unitPocean  = 1.0E6 !bar
    unitKocean = 1. !cm2/g
    unitDocean = 1. !g/cm**3

    ! Read parameters
    open(unit=1,file=fparam)
    read(unit=1,nml=eos,iostat=stat)
     IF (stat/=0) THEN
      BACKSPACE(1)
      READ(1,fmt='(A)') line
      WRITE(*,'(A)') 'Invalid line in namelist: '//trim(line)
      STOP
     END IF
     close(1)

     ppmin_ocean  = ppmin_ocean*unitPocean/unitP
     ppmax_ocean  = ppmax_ocean*unitPocean/unitP

     opacity_ocean = opacity_ocean*unitKocean/(unitL**2/unitM)
     density_ocean = density_ocean*unitDocean/(unitM/unitL**3)

    if (eos_type=='scvh95') then

      ! loading the tables
      tiny=1.0d-14

      open(unit=71,file=trim(eos_dir)//trim(eos_file_Y25),status='old' &
            ,form='unformatted')
      open(unit=72,file=trim(eos_dir)//trim(eos_file_Y30),status='old' &
            ,form='unformatted')
      do n=1,ntab
        iunit=n+70
        read(iunit) nt,np,yt(n)
!       write(*,*)nt,np,yt(n)
        do ip=1,np
          do jt=1,nt
            read(iunit) tlgt(jt),plgt(ip),rolgt(n,ip,jt),elgt(n,ip,jt) &
                   ,slgt(n,ip,jt),cplgt(n,ip,jt),gamtab(n,ip,jt),abadtab(n,ip,jt)
          enddo
        enddo
      enddo

      close(71)
      close(72)

      ! Interpolate for the metallicity here in data
      imin=0
      imax=0
      do i=1,ntab
        if(xHe.ge.yt(i))imin=i
      enddo
      imax=imin+1
      if(xHe.eq.yt(ntab))then
        imin=ntab-1
        imax=ntab
      endif

      !check that xHe is within table
      if(imin.eq.0.or.imin.eq.ntab)then
        write(6,*)"helium mass fraction",xHe &
                   ,"outside range of table values of 0.275--0.30"
        stop
      endif

      alpha = (xHe-yt(imin))/(yt(imax)-yt(imin))

      data(1,:,:) = (1d0-alpha)*rolgt(imin,:,:) + alpha*rolgt(imax,:,:)
      data(2,:,:) = (1d0-alpha)*elgt(imin,:,:) + alpha*elgt(imax,:,:)
      data(3,:,:) = (1d0-alpha)*slgt(imin,:,:) + alpha*slgt(imax,:,:)
      data(4,:,:) = (1d0-alpha)*cplgt(imin,:,:) + alpha*cplgt(imax,:,:)
      data(5,:,:) = (1d0-alpha)*gamtab(imin,:,:) + alpha*gamtab(imax,:,:)
      data(6,:,:) = (1d0-alpha)*abadtab(imin,:,:) + alpha*abadtab(imax,:,:)

    else if (eos_type/='ideal') then
      if (mype==cputerm) then
        write(*,'(A)') 'This eos is not implemented, eos should be ideal or scvh95'
        write(*,chrpmfmt)'eos_type',eos_type
        write(*,*)
        write(*,separfmt)
        stop
      endif
    endif

    ! units
    pmin_conv = pmin_conv*1E6/unitP
    pmax_conv = pmax_conv*1E6/unitP

if (mype==cputerm) then
    write(*,headrfmt)'Equation of state'
    write(*,chrpmfmt)'eos_type',eos_type,'Saumon Chabrier Van Horn 1995 or ideal'
    if (eos_type=='scvh95') then
        write(*,chrpmfmt)'eos_dir',eos_dir, 'directory with eos table'
        write(*,fltpmfmt)'xHe',xHe, 'effective helium mass fraction'
    else if (eos_type=='ideal') then
       write(*,fltpmfmt) 'mucst',mucst,'If = 0 then compute mu from abundances'
    else
       write(*,'(A)') 'The requested eos is not implemented, should be ideal or scvh95'
       write(*,*)
       write(*,separfmt)
       stop
    endif
    write(*,exppmfmt)'gamma',gamma,'Adiabatic exponent, if = 0 use thermodynamic data'
    if (gamma/=0.0) then
       write(*,exppmfmt) 'pmin_conv [bar]',pmin_conv*unitP/1E6,'Pressure above which gradad is not constant'
       write(*,exppmfmt) 'pmax_conv [bar]',pmax_conv*unitP/1E6,'Pressure below which gradad is not constant'
    endif
    write(*,*)
    write(*,separfmt)
endif

end subroutine init_eos

  ! Compute either ideal or Saumon, Chabrier and van Horn 1995 EOS
  subroutine compute_eos(ndepth_loc,pp_loc,tt_loc,rho_loc,cp_loc,abad_loc,eos_flag)

    use mod_cst
    use mod_param
    use mod_chem
    use mod_grid

    implicit none

    integer, intent(in) :: ndepth_loc
    real, dimension(ndepth_loc), intent(inout)  :: pp_loc,tt_loc
    real, dimension(ndepth_loc)  :: grad_smooth
    real, dimension(ndepth_loc), intent(out) :: rho_loc,cp_loc,abad_loc
    character(len=6), dimension(ndepth_loc), intent(out) :: eos_flag
    real              :: eint, sint, gam
    real              :: tcalc,interp
    integer   :: idepth,iker,ismooth

    ! Calculate density, specific heat and adiabatic gradient in the ideal case
    ! This will be overidden if the eos_type is scvh95 and we are in the P, T limits of the tables
    rho_loc = get_rho(ndepth_loc,pp_loc,tt_loc,Amolf)
    cp_loc = get_cpspec(ndepth_loc,pp_loc,tt_loc,Amolf)
    abad_loc = get_gradad(ndepth_loc,pp_loc,tt_loc,Amolf)
    eos_flag = 'ideal'

    if (eos_type=='scvh95') then

      ! Convert input variables pem and tem to cgs units for eossc call
      pp_loc = pp_loc * unitP
      tt_loc = tt_loc * unitK

      do idepth=1,ndepth_loc
        ! If we are in the table limits, treat low T limit seperately
        if ((pp_loc(idepth)>=10.0**plgt(2)) .and. (pp_loc(idepth)<=10.0**plgt(np-1)) .and. (tt_loc(idepth)<=10.0**tlgt(nt-1))) then

          ! Below low T limit, keep EOS frozen at last valid value
          if (tt_loc(idepth)<10.0**tlgt(2)) then
            tcalc = (10.**tlgt(2))*1.001
          else
            tcalc = tt_loc(idepth)
          endif

          ! Call eossc
          call eossc(pp_loc(idepth),tcalc,rho_loc(idepth),eint,sint,cp_loc(idepth),gam,abad_loc(idepth))
          eos_flag(idepth) = eos_type
          ! Convert from cgs to code units
          rho_loc(idepth) = rho_loc(idepth) / (unitM / unitL**3)
          cp_loc(idepth)  = cp_loc(idepth)  / (unitE / (unitM*unitT))
        endif
      enddo

      ! Convert from cgs to code units
      pp_loc   = pp_loc  / unitP
      tt_loc   = tt_loc  / unitK

      if (.not. gamma==0.) then

         do idepth = 1, ndepth_loc

            if (idepth>1 .and. pp_loc(idepth)>pmin_conv .and. pp_loc(idepth-1)<pmin_conv) then
               interp = abs((log(pmin_conv)-log(pp_loc(idepth)))/(log(pp_loc(idepth))-log(pp_loc(idepth-1))))
               abad_loc(idepth) =  abad_loc(idepth)*(1.-interp)+ (gamma-1.0)/gamma*interp

            else if (idepth<ndepth+1 .and. pp_loc(idepth)<pmax_conv .and. pp_loc(idepth+1)>pmax_conv) then
               interp = abs((log(pmax_conv)-log(pp_loc(idepth)))/(log(pp_loc(idepth))-log(pp_loc(idepth+1))))
               abad_loc(idepth) =  abad_loc(idepth)*(1.-interp)+ (gamma-1.0)/gamma*interp

            else if (pp_loc(idepth)>pmin_conv .and. pp_loc(idepth)<pmax_conv)then
               abad_loc(idepth) = (gamma-1.0)/gamma

          endif

       enddo

    endif

    if (ocean) then

         do idepth = 1, ndepth_loc

            if (idepth>1 .and. pp_loc(idepth)>ppmin_ocean .and. pp_loc(idepth-1)<ppmin_ocean) then
               interp = abs((log(ppmin_ocean)-log(pp_loc(idepth)))/(log(pp_loc(idepth))-log(pp_loc(idepth-1))))
               abad_loc(idepth) = abad_loc(idepth) * (1.-interp) + (gamma_ocean-1.0)/gamma_ocean*interp
               rho_loc(idepth)  = rho_loc(idepth)  * (1.-interp) + density_ocean*interp
               cp_loc(idepth)   = cp_loc(idepth)   * (1.-interp) + gamma_ocean/(gamma_ocean-1.0)*kb/(mu_ocean*amu)*interp
            else if (idepth<ndepth+1 .and. pp_loc(idepth)<ppmax_ocean .and. pp_loc(idepth+1)>ppmax_ocean) then
               interp = abs((log(ppmax_ocean)-log(pp_loc(idepth)))/(log(pp_loc(idepth))-log(pp_loc(idepth+1))))
               abad_loc(idepth) =  abad_loc(idepth)*(1.-interp)+ (gamma_ocean-1.0)/gamma_ocean*interp
               rho_loc(idepth)  = rho_loc(idepth)  * (1.-interp) + density_ocean*interp
               cp_loc(idepth)   = cp_loc(idepth)   * (1.-interp) + gamma_ocean/(gamma_ocean-1.0)*kb/(mu_ocean*amu)*interp
            else if (pp_loc(idepth)>ppmin_ocean .and. pp_loc(idepth)<ppmax_ocean)then
               abad_loc(idepth) = (gamma_ocean-1.0)/gamma_ocean
               rho_loc(idepth)  = density_ocean*interp
               cp_loc(idepth)   = gamma_ocean/(gamma_ocean-1.0)*kb/(mu_ocean*amu)
            endif

       enddo

    endif


       if (gradad_smooth) then
          do ismooth = 1,nb_eos_smooth
             grad_smooth = abad_loc

             do idepth = 1+kereos_smooth,ndepth_loc-kereos_smooth

                do iker=1,kereos_smooth
                   grad_smooth(idepth) = grad_smooth(idepth) + abad_loc(idepth-iker) + abad_loc(idepth+iker)
                enddo

                grad_smooth(idepth) = grad_smooth(idepth)/(2*kereos_smooth+1)

             enddo

             do idepth = 1, ndepth_loc
                abad_loc(idepth) = grad_smooth(idepth)
                if (idepth<kereos_smooth+1) then
                   abad_loc(idepth) = grad_smooth(kereos_smooth+1)
                else if (idepth>ndepth_loc-kereos_smooth) then
                   abad_loc(idepth) = grad_smooth(ndepth_loc-kereos_smooth)
                endif

             enddo
          enddo
    endif

    else
      if (eos_type/='ideal') then
        if (mype==cputerm) then
          write(*,'(A)') 'This eos is not implemented, eos should be ideal or scvh95'
          write(*,chrpmfmt) 'eos_type', eos_type
          write(*,*)
          write(*,separfmt)
          stop
        endif
      endif
    endif

  end subroutine compute_eos

  ! Compute Saumon, Chabrier and van Horn 1995 EOS
  subroutine eossc(pem,tem,roint,eint,sint,cp,gamma,abad)
    real(kind=8), intent(in) :: pem, tem
    real(kind=8), intent(out) :: roint,eint,sint
    real(kind=8), intent(out) :: cp,gamma,abad
    real(kind=8), dimension(6) :: q, qp, qt
    real(kind=8), dimension(4) :: yrate, y1rate, y2rate, y12rate
    real(kind=8) :: plg,tlg,dplgt,dtlgt,dp,dt,ff1,ff2,ff3,ff4,ff5,ff6,ff7,ff8,ff9,fff,ffp,fft
    integer :: ip,jt,iq,ii,iim,iin

    plg=dlog10(pem)
    tlg=dlog10(tem)

    call mylocate(plgt,np,plg,ip)

    if(ip.lt.2.or.ip.gt.(np-2))then
       write(6,*) 'stop in eossc: check iro=',ip,plg
       stop
    endif

    call mylocate(tlgt,nt,tlg,jt)

    if(jt.lt.2.or.jt.gt.(nt-2))then
       write(6,*) 'stop in eossc: check jt=',jt,tlg
       stop
    endif

    dplgt = plgt(ip+1)-plgt(ip)
    dtlgt = tlgt(jt+1)-tlgt(jt)
    dp = (plg-plgt(ip))/dplgt
    dt = (tlg-tlgt(jt))/dtlgt

    !interpolation of the quantities required
    do iq=1,6
       do ii=1,4
          if(ii.eq.1) then
             iim = ip
             iin = jt
          end if
          if(ii.eq.2) iim = ip+1
          if(ii.eq.3) iin = jt+1
          if(ii.eq.4) iim = ip

          ff1 = data(iq,iim,iin)

          ff2 = data(iq,iim+1,iin)
          ff3 = data(iq,iim-1,iin)

          ff4 = data(iq,iim,iin+1)
          ff5 = data(iq,iim,iin-1)

          ff6 = data(iq,iim+1,iin+1)
          ff7 = data(iq,iim+1,iin-1)
          ff8 = data(iq,iim-1,iin+1)
          ff9 = data(iq,iim-1,iin-1)

          yrate(ii) = ff1
          y1rate(ii) = (ff2-ff3)/(2.d0*dplgt)
          y2rate(ii) = (ff4-ff5)/(2.d0*dtlgt)
          y12rate(ii) = (ff6-ff7-ff8+ff9)/(4.d0*dplgt*dtlgt)
       enddo
       call intercub(yrate,y1rate,y2rate,y12rate,dplgt,dtlgt,dp,dt,fff,ffp,fft)
       q(iq)=fff
       qp(iq)=ffp/dplgt
       qt(iq)=fft/dtlgt
    enddo

    roint=10.d0**q(1)
    eint=10.d0**q(2)
    sint=10.d0**q(3)
    cp=10.d0**q(4)
    gamma=q(5)
    abad=q(6)

    return
  end subroutine eossc

  !> Compute the number density
  !<
  function get_nn(ndepth_loc,pp_loc,tt_loc)

    use mod_grid
    use mod_cst

    implicit none

    integer :: ndepth_loc
    real,dimension(ndepth_loc)      :: pp_loc,tt_loc,get_nn

    get_nn = pp_loc/(kb*tt_loc)

    if(ocean) then
       where(pp_loc>ppmin_ocean .and. pp_loc<ppmax_ocean) get_nn = density_ocean/mu_ocean
    endif

  end function get_nn


  !> Compute the density in function of the mean molecular weights computed from the chemistry at cell centres or interfaces
  !<
  function get_rho(ndepth_loc,pp_loc,tt_loc,Amol_loc)

    use mod_grid
    use mod_cst
    use mod_chem

    implicit none

    integer :: ndepth_loc
    real,dimension(ndepth_loc)      :: pp_loc,tt_loc,get_rho
    real,dimension(ndepth_loc,nmol) :: Amol_loc

    get_rho = get_nn(ndepth_loc,pp_loc,tt_loc)*(get_mu(ndepth_loc,Amol_loc)*amu)

    if(ocean) then
       where(pp_loc>ppmin_ocean .and. pp_loc<ppmax_ocean) get_rho = density_ocean
    endif

  end function get_rho


  !> Compute number density gradient with respect to pp
  !<
  function get_dndp(ndepth_loc,pp_loc,tt_loc)

    use mod_cst
    use mod_grid

    implicit none

    integer :: ndepth_loc
    real,dimension(ndepth_loc) :: get_dndp,tt_loc,pp_loc


    get_dndp = 1./(kb*tt_loc)
    if(ocean) then
       where(pp_loc>ppmin_ocean .and. pp_loc<ppmax_ocean) get_dndp = 0.
    endif

  end function get_dndp

  !> Compute density gradient with respect to tt
  !<
  function get_dndt(ndepth_loc,pp_loc,tt_loc)

    use mod_cst
    use mod_grid

    implicit none

    integer :: ndepth_loc
    real,dimension(ndepth_loc) :: get_dndt,tt_loc,pp_loc

    get_dndt = -pp_loc/(kb*tt_loc**2)
    if(ocean) then
       where(pp_loc>ppmin_ocean .and. pp_loc<ppmax_ocean) get_dndt = 0.
    endif

  end function get_dndt


  !> Compute density gradient with respect to pp
  !<
  function get_drhodp(ndepth_loc,pp_loc,tt_loc,Amol_loc)

    use mod_cst
    use mod_grid
    use mod_chem

    implicit none

    integer :: ndepth_loc
    real,dimension(ndepth_loc) :: get_drhodp,tt_loc,pp_loc
    real,dimension(ndepth_loc,nmol) :: Amol_loc

    get_drhodp = get_dndp(ndepth_loc,pp_loc,tt_loc)*(get_mu(ndepth_loc,Amol_loc)*amu)

    if(ocean) then
       where(pp_loc>ppmin_ocean .and. pp_loc<ppmax_ocean) get_drhodp = 0.
     endif

  end function get_drhodp

  !> Compute density gradient with respect to tt
  !<
  function get_drhodt(ndepth_loc,pp_loc,tt_loc,Amol_loc)

    use mod_cst
    use mod_grid
    use mod_chem

    implicit none

    integer :: ndepth_loc
    real,dimension(ndepth_loc) :: get_drhodt,tt_loc,pp_loc
    real,dimension(ndepth_loc,nmol) :: Amol_loc

    get_drhodt = get_dndt(ndepth_loc,pp_loc,tt_loc)*(get_mu(ndepth_loc,Amol_loc)*amu)

    if(ocean) then
       where(pp_loc>ppmin_ocean .and. pp_loc<ppmax_ocean) get_drhodt = 0.
    endif

  end function get_drhodt

  !> Compute specific cp
  !<
  function get_cpspec(ndepth_loc,pp_loc,tt_loc,Amol_loc)

    use mod_cst
    use mod_grid
    use mod_chem
    use mod_param

    implicit none

    integer :: ndepth_loc,idepth,i
    real,dimension(ndepth_loc) :: get_cpspec,tt_loc,pp_loc
    real,dimension(ndepth_loc,nmol) :: Amol_loc
    real,dimension(nmol) :: mu, h, s, cpmoldR

    if (.not. gamma == 0. ) then

       get_cpspec = gamma/(gamma-1.0)
       get_cpspec = get_cpspec*kb/(get_mu(ndepth_loc,Amol_loc)*amu)

    else if (calc_cp_chem) then

      if (ndepth_loc == ndepth+1) then
          get_cpspec = cp_chem
      else
          write(*,headrfmt) 'Error get_cpspec'
          write(*,'(A)') 'calc_cp_chem selected, but number of levels do not match'
          write(*,*)
          write(*,separfmt)
          stop
      end if
    else

       get_cpspec = 0.

       do idepth=1,ndepth_loc

         ! calculate thermodynamic quantities
         CALL get_thermodynamics(nmol,tt_loc(idepth),molParam,mu,h,s,cpmoldR)

         ! calculate heat capacity of mixture
          do i = 1,nmol
             if (molphase(i) == 'g') then
               get_cpspec(idepth) = get_cpspec(idepth) + Amol_loc(idepth,i)*cpmoldR(i)
             endif
          enddo
       enddo

       get_cpspec = get_cpspec*kb/(get_mu(ndepth_loc,Amol_loc)*amu)!3.655*kb/(get_mu(ndepth_loc,Amol_loc)*amu)!get_cpspec*kb/(get_mu(ndepth_loc,Amol_loc)*amu)!3.655*kb/(get_mu(ndepth_loc,Amol_loc)*amu)

    endif

    if(ocean) then
       where(pp_loc>ppmin_ocean .and. pp_loc<ppmax_ocean) get_cpspec = gamma_ocean/(gamma_ocean-1.0)*kb/(mu_ocean*amu)
    endif

  end function get_cpspec


  !> Compute adiabatic gradient
  !<
  function get_gradad(ndepth_loc,pp_loc,tt_loc,Amol_loc)

    use mod_cst
    use mod_grid
    use mod_chem

    implicit none

    integer :: ndepth_loc,i,idepth,iker
    real,dimension(ndepth_loc) :: get_gradad,tt_loc,pp_loc,grad_smooth
    real,dimension(ndepth_loc,nmol) :: Amol_loc
    real,dimension(nmol) :: mu,h,s,cpmoldR
    real :: cpmoldR_tot



        get_gradad = 0.

        do idepth=1,ndepth_loc


           if (.not. gamma == 0. .and. ppf(idepth)>pmin_conv .and. ppf(idepth)<pmax_conv)then

              get_gradad(idepth) = (gamma-1.0)/gamma

           else

          cpmoldR_tot = 0.

          ! Calculate thermodynamic quantities
          CALL get_thermodynamics(nmol,tt_loc(idepth),molParam,mu,h,s,cpmoldR)

          do i = 1,nmol
            if (molphase(i) == 'g') then
              cpmoldR_tot = cpmoldR_tot + Amol_loc(idepth,i)*cpmoldR(i)
            endif
          enddo


          get_gradad(idepth) = 1./cpmoldR_tot !1./3.655 0.25

          ! if (ttf(idepth)>3000./unitK) get_gradad(idepth) = get_gradad(idepth-1)

       endif
    enddo

    if(ocean) then
       where(pp_loc>ppmin_ocean .and. pp_loc<ppmax_ocean) get_gradad = (gamma_ocean-1.0)/gamma_ocean
    endif

    if (gradad_smooth) then

       grad_smooth = get_gradad

       do idepth = 1+kereos_smooth,ndepth+1-kereos_smooth

          do iker=1,kereos_smooth
             grad_smooth(idepth) = grad_smooth(idepth) + get_gradad(idepth-iker) + get_gradad(idepth+iker)
          enddo


          grad_smooth(idepth) = grad_smooth(idepth)/(2*kereos_smooth+1)

          get_gradad(idepth) = grad_smooth(idepth)

       enddo


    endif


          ! if (chem=='eq') then
          !    !get_gradad(idepth) =1./cp_chem(idepth)
          ! endif


  end function get_gradad


  !> function to compute the mean molecular weight along the profile
  !<
  function get_mu(ndepth_loc,Amol_loc)

    use mod_chem

    implicit none

    integer:: ndepth_loc,idepth

    real,dimension(ndepth_loc,nmol) :: Amol_loc
    real,dimension(ndepth_loc) :: get_mu

    do idepth =1,ndepth_loc

       get_mu(idepth) = sum(Amol_loc(idepth,1:nmolg)*molmass(1:nmolg))

    enddo

    if(.not. mucst == 0.) then

       get_mu = mucst

    endif

  end function get_mu


  !> compute the vertical scale assuming hydrostatic equilibrium\n
  !! assume that rho is computed before
  !<
  subroutine get_rrgr

    use mod_grid
    use mod_param
    use mod_cst

    implicit none

    integer :: idepth,idepth_Rp


    if (Rp /=0.) then

       if(.not. pp_Rp == 0.) then

          do idepth=1,ndepth+1

             if (ppf(idepth)<pp_Rp) then
                idepth_Rp=idepth
             endif

          enddo

          rrf(idepth_Rp) = Rp

          do idepth = idepth_Rp-1,1,-1

             rrf(idepth) = 1./(1./rrf(idepth+1)+(ppf(idepth)-ppf(idepth+1))/(rho(idepth)*grav*Rp**2))

             rr(idepth) = 2./(1./rrf(idepth)+1./rrf(idepth+1))

          enddo

          do idepth = idepth_Rp+1,ndepth+1

             rrf(idepth) = 1./(1./rrf(idepth-1)-(ppf(idepth-1)-ppf(idepth))/(rho(idepth-1)*grav*Rp**2))

             if (idepth<ndepth+1) then
                rr(idepth) = 2./(1./rrf(idepth)+1./rrf(idepth-1))
             endif

          enddo


       else

          rrf(ndepth+1) = Rp

          do idepth = ndepth,1,-1

             rrf(idepth) = 1./(1./rrf(idepth+1)+(ppf(idepth)-ppf(idepth+1))/(rho(idepth)*grav*Rp**2))

             rr(idepth) = 2./(1./rrf(idepth)+1./rrf(idepth+1))

          enddo

       endif
       gr  = grav*Rp**2/rr**2
       grf = grav*Rp**2/rrf**2

    else

       do idepth = ndepth,1,-1

          rrf(idepth) = rrf(idepth+1)-(ppf(idepth)-ppf(idepth+1))/(rho(idepth)*grav)
          rr(idepth) = (rrf(idepth)+rrf(idepth+1))/2.

       enddo



       gr  = grav
       grf = grav

    endif


  end subroutine get_rrgr

  ! Calculate thickness of layers
  subroutine get_drdrf

  use mod_grid

  implicit none

  integer :: idepth

  DO idepth = 2,ndepth
    drf(idepth) = (rr(idepth)-rr(idepth-1))
  END DO

  DO idepth = 1,ndepth
    dr(idepth) = (rrf(idepth+1)-rrf(idepth))
  END DO

  ! Boundaries
  drf(1) = (rr(2)-rr(1))
  drf(ndepth+1) = (rr(ndepth)-rr(ndepth-1))

  end subroutine

  !> Copy vertical scale and gravity from UM output
  !<
  subroutine get_um_rrgr

    use mod_grid
    use mod_param

    if (Rp == 0) then
      write(*,headrfmt) 'Error get_um_rrgr'
      write(*,'(A)') 'Rp == 0. Please use radius consistent with UM output.'
      write(*,*)
      write(*,separfmt)
      stop
    endif

    ! Set vertical level heights
    rr = Rp + height_um
    rrf = Rp + heightf_um

    ! Calculate gravity
    gr  = grav*Rp**2/rr**2
    grf = grav*Rp**2/rrf**2

  end subroutine get_um_rrgr

  !>  bicubic interpolation with derivatives
  subroutine intercub(y,y1,y2,y12,d1,d2,t,u,rate,ratero,ratet)
    real(kind=8), dimension(4), INTENT(IN) :: y,y1,y2,y12
    real(kind=8), INTENT(IN) :: d1,d2,t,u
    REAL(kind=8) :: rate,ratero,ratet
    integer :: i
    real(kind=8), dimension(4,4) :: c
    real(kind=8), dimension(16) :: x,temp

    temp(1:4) = y(:)
    temp(5:8) = y1(:)*d1
    temp(9:12) = y2(:)*d2
    temp(13:16) = y12*d1*d2


    !MATRIX VECTOR PRODUCT for solving the linear problem
    !as part of cubic interpolation. See further notes at end of file
    x(1) = temp(1)
    x(2) = temp(9)
    x(3) = -3.d0*temp(1) + 3.d0*temp(4) -2.d0*temp(9) - temp(12)
    x(4) = 2.d0*temp(1) - 2.d0*temp(4) + temp(9) + temp(12)
    x(5) = temp(5)
    x(6) = temp(13)
    x(7) = -3.d0*temp(5) + 3.d0*temp(8) -2.d0*temp(13) -temp(16)
    x(8) = 2.d0*temp(5) - 2.d0*temp(8) + temp(13) + temp(16)
    x(9) = -3.d0*temp(1) + 3.d0*temp(2) -2.d0*temp(5)  - temp(6)
    x(10) = -3.d0*temp(9) + 3.d0*temp(10)  - 2.d0*temp(13) - temp(14)
    x(11) = 9.d0*temp(1) - 9.d0*temp(2) + 9.d0*temp(3) -9.d0*temp(4) &
         + 6.d0*temp(5) +3.d0*temp(6) -3.d0*temp(7) -6.d0*temp(8)  &
         + 6.d0*temp(9) -6.d0*temp(10) -3.d0*temp(11) +3.d0*temp(12) &
         +4.d0*temp(13)  +2.d0*temp(14) +temp(15) +2.d0*temp(16)
    x(12) = -6.d0*temp(1) +6.d0*temp(2) -6.d0*temp(3) + 6.d0*temp(4)  &
         -4.d0*temp(5) -2.d0*temp(6) + 2.d0*temp(7) + 4.d0*temp(8)  &
         -3.d0*temp(9) +3.d0*temp(10) + 3.d0*temp(11)  -3.d0*temp(12)    &
         -2.d0*temp(13) - temp(14) - temp(15)  - 2.d0*temp(16)
    x(13) = 2.d0*temp(1) -2.d0*temp(2) + temp(5) + temp(6)
    x(14) = 2.d0*temp(9) -2.d0*temp(10) + temp(13) + temp(14)
    X(15) = -6.d0*temp(1) + 6.d0*temp(2) -6.d0*temp(3)+ 6.d0*temp(4) &
         - 3.d0*temp(5) -3.d0*temp(6) +3.d0*temp(7) +3.d0*temp(8) &
         -4.d0*temp(9) +4.d0*temp(10) +2.d0*temp(11) -2.d0*temp(12) &
         -2.d0*temp(13) -2.d0*temp(14) - temp(15) - temp(16)
    x(16) = 4.d0*temp(1) - 4.d0*temp(2) + 4.d0*temp(3) - 4.d0*temp(4) &
         + 2.d0*temp(5) + 2.d0*temp(6) -2.d0*temp(7)  -2.d0*temp(8) &
         +2.d0*temp(9) -2.d0*temp(10) - 2.d0*temp(11)  + 2.d0*temp(12) &
         +temp(13) + temp(14) + temp(15) + temp(16)

    c = reshape(x, (/4,4/), order=(/2,1/))

    rate = 0.d0
    ratet = 0.d0
    ratero = 0.d0
    do i=4,1,-1
       rate = rate*t + ((c(i,4)*u+c(i,3))*u+c(i,2))*u + c(i,1)
       ratet = ratet*t + (3.d0*c(i,4)*u+2.d0*c(i,3))*u + c(i,2)
       ratero = ratero*u + (3.d0*c(4,i)*t+2.d0*c(3,i))*t + c(2,i)
    enddo
    return
  end subroutine intercub

  subroutine mylocate(xx,n,x,jl)
    integer ::j,n
    real(kind=8) ::x,xx(n)
    integer ::jl,jm,ju

    jl=0
    ju=n+1
10  if(ju-jl.gt.1)then
       jm=(ju+jl)/2
       if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
       else
          ju=jm
       endif
       goto 10
    endif
    j=jl
    return
  end subroutine mylocate

end module mod_eos
