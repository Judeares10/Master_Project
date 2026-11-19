!> \file
!! Contains module mod_conv  \n
!! Routines used for the convection
!<

!> hold the parameters for the convection
!<
module mod_conv


  real :: nuconv                !< Turbulent viscosity
  real :: alpha                 !< Mixing length
  real :: yparam                !< Optical depth scaling
  real :: pstop_conv            !< Pressure [bar] below which convection is removed
  real :: eps_gradad            !< Error allowed on gradad
  real :: eps_vconv             !< Factor for vconv decrease in jacobian
  real :: pp1d_decoupling_alpha !< Pressure below which mixing length is decreased as a power law with exponent decoupling_power
  integer :: decoupling_power
  real :: pp1d_neg_fconv        !< Pressure above which negative convective fluxes are allowed
  real :: pp1d_force_vconv      !< Pressure above which convective velocities are forced to 10% of local sounds speed
  real :: unitP_pp1d            !< Pressure units for pressure parameters

contains

  !> Initialize the variables for the convection
  !<
  subroutine init_conv

    use mod_grid
    use mod_param
    use mod_cst

    implicit none

    integer :: stat

    namelist/convection/ alpha,nuconv,yparam,pstop_conv,eps_gradad,eps_vconv,pp1d_decoupling_alpha,unitP_pp1d,decoupling_power,pp1d_neg_fconv,pp1d_force_vconv

    ! Set standard values
    alpha      = 1.5
    nuconv     = 8.0
    yparam     = 0.076
    pstop_conv = 0.0
    eps_gradad = 1.0E-2
    eps_vconv  = 0.0
    pp1d_decoupling_alpha  = 1E-10
    unitP_pp1d = 1.0E6  !bars
    decoupling_power = 0
    pp1d_neg_fconv = 1E10
    pp1d_force_vconv = 1E10

    ! Read parameters
    open(unit=1,file=fparam)
    read(unit=1,nml=convection,iostat=stat)
     IF (stat/=0) THEN
      BACKSPACE(1)
      READ(1,fmt='(A)') line
      WRITE(*,'(A)') 'Invalid line in namelist: '//trim(line)
      STOP
     END IF
    close(1)

    ! Switch from bar to cgs to code units
    pstop_conv = pstop_conv*unitPgrid/unitP
    pp1d_decoupling_alpha = pp1d_decoupling_alpha*unitP_pp1d/unitP
    pp1d_neg_fconv = pp1d_neg_fconv*unitP_pp1d/unitP
    pp1d_force_vconv = pp1d_force_vconv*unitP_pp1d/unitP


    if (mype==cputerm) then
       write(*,headrfmt)'Convection'
       write(*,exppmfmt)'pstop_conv',pstop_conv,'Pressure [bar] below which convection is removed'
       write(*,exppmfmt)'alpha',alpha,'Mixing length [dimless]'
       write(*,exppmfmt)'nuconv',nuconv,'Turbulent viscosity [dimless]'
       write(*,exppmfmt)'yparam',yparam,'Optical depth scaling [dimless]'
       write(*,exppmfmt)'eps_gradad',eps_gradad,'Error allowed on gradad'
       write(*,*)
       write(*,separfmt)
    endif


  end subroutine init_conv


  !> Compute the convective flux at cell interfaces
  !<
  subroutine get_convection

    use mod_grid
    use mod_eos
    use mod_cst
    use mod_param
    use mod_chem
    use mod_opacity

    implicit none

    !local variable

     integer :: idepth

    real,dimension(ndepth+1) :: hp,lconv,qconv,theta_loc,coeff,delta,drhofdt,drhofdp
    real,dimension(ndepth+1) :: gam,gradT,grad,vconv,kapfross,cpspec,gradad_test,gradadeps,dgradaddp,dgradaddt,ppfeps,ttfeps
    real,dimension(ndepth+1,2*(ndepth+1)) :: dgradT
    character(len=6),dimension(ndepth+1)  :: eos_flag
    logical :: recompute_gradad
    real :: ppp,ppm,ttp,ttm,convdpres,eps

    dfconv = 0.
    fconv  = 0.
    dgradT = 0.
    kconv  = 0.

    eps=1E-3

    !call compute eos
    ttfeps = ttf*(1+eps)
    ppfeps = ppf*(1+eps)

    call compute_eos(ndepth+1,ppf,ttf,rhof,cpspec,gradad_test,eos_flag)
    call compute_eos(ndepth+1,ppf,ttfeps,rhof,cpspec,gradadeps,eos_flag)

    dgradaddt = (gradadeps - gradad_test)/(eps*ttf)

    call compute_eos(ndepth+1,ppfeps,ttf,rhof,cpspec,gradadeps,eos_flag)

    dgradaddp = (gradadeps - gradad_test)/(eps*ppf)


    call compute_eos(ndepth+1,ppf,ttf,rhof,cpspec,gradad_test,eos_flag)

    !compute density
    drhofdt = get_drhodt(ndepth+1,ppf,ttf,Amolf)
    drhofdp = get_drhodp(ndepth+1,ppf,ttf,Amolf)

    !compute rosseland opacities
    kapfross = get_kapross(ndepth+1,ttf,kapfnu)

    !the local pressure scale height
    hp = ppf/(grav*rhof)

    ! the mixing length
    lconv = alpha*hp

    do idepth=1,ndepth+1

       if (ppf(idepth)<pp1d_decoupling_alpha) then
          lconv(idepth) = lconv(idepth)*(ppf(idepth)/pp1d_decoupling_alpha)**decoupling_power
       endif

    enddo

    ! thermal expansion coefficient
    qconv = -ttf/rhof*get_drhodt(ndepth+1,ppf,ttf,Amolf)

    if (init_gradad) then
       gradad = gradad_test
       write(*,*) 'Recomputed gradad'
       init_gradad = .false.
    endif

    recompute_gradad = .false.

    do idepth=1,ndepth+1
       if((abs(gradad(idepth) - gradad_test(idepth))/gradad_test(idepth))>eps_gradad) then
          recompute_gradad = .true.
       endif

    enddo


    if (recompute_gradad) then
      gradad = gradad_test
      write(*,*) 'Recomputed gradad'
    endif

    do idepth = 2,ndepth

       !linear interpolation for tt at cell centre
       ttp = ttf(idepth+1)
       ttm = ttf(idepth-1)

       !log-linear interpolation for pp at cell centre
       ppp = ppf(idepth+1)
       ppm = ppf(idepth-1)

       gradT(idepth) = (log10(ttp)-log10(ttm))/(log10(ppp)-log10(ppm))
       dgradT(idepth,ndepth+2+idepth) =  1./ttp/(log10(ppp)-log10(ppm))
       dgradT(idepth,ndepth+idepth) = -1./ttm/(log10(ppp)-log10(ppm))
       dgradT(idepth,idepth+1) = -1./ppp*(log10(ttp)-log10(ttm))/(log10(ppp)-log10(ppm))**2
       dgradT(idepth,idepth-1) =  1./ppm*(log10(ttp)-log10(ttm))/(log10(ppp)-log10(ppm))**2

       ! if (ppf(idepth)<pp1d_decoupling_alpha) then
       !    dgradT(idepth,idepth) = dgradT(idepth,idepth) - dgradaddp(idepth)
       !    dgradT(idepth,ndepth+1+idepth) = dgradT(idepth,ndepth+1+idepth) - dgradaddt(idepth)
       ! else
       !    dgradT(idepth,idepth) = (dgradT(idepth,idepth) - dgradaddp(idepth))!*(gradT(idepth)-gradad(idepth))
       !    dgradT(idepth,ndepth+1+idepth) = (dgradT(idepth,ndepth+1+idepth) - dgradaddt(idepth))!*(gradT(idepth)-gradad(idepth))
       ! endif

    enddo

    !linear interpolation for tt at cell centre
    ttp = ttf(2)
    ttm = ttf(1)**2/ttf(2)

    !log-linear interpolation for pp at cell centre
    ppp = ppf(2)
    ppm = ppf(1)**2/ppf(2)

    !extrapolate at boundaries
    gradT(1)        = (log10(ttp)-log10(ttm))/(log10(ppp)-log10(ppm))

    dgradT(1,ndepth+1+2) =  2./ttf(2)/(log10(ppp)-log10(ppm))
    dgradT(1,ndepth+1+1) = -2./ttf(1)/(log10(ppp)-log10(ppm))

    dgradT(1,2) =  -2./ppf(2)*(log10(ttp)-log10(ttm))/(log10(ppp)-log10(ppm))**2
    dgradT(1,1) =  2./ppf(1)*(log10(ttp)-log10(ttm))/(log10(ppp)-log10(ppm))**2

    !linear interpolation for tt at cell centre
    ttp = ttf(ndepth+1)**2/ttf(ndepth)
    ttm = ttf(ndepth)

    !log-linear interpolation for pp at cell centre
    ppp = ppf(ndepth+1)**2/ppf(ndepth)
    ppm = ppf(ndepth)

    gradT(ndepth+1) = (log10(ttp)-log10(ttm))/(log10(ppp)-log10(ppm))

    dgradT(ndepth+1,2*ndepth+1)   = -2./ttf(ndepth)/(log10(ppp)-log10(ppm))
    dgradT(ndepth+1,2*(ndepth+1)) =  2./ttf(ndepth+1)/(log10(ppp)-log10(ppm))

    dgradT(ndepth+1,ndepth)   =  2./ppf(ndepth)*(log10(ttp)-log10(ttm))/(log10(ppp)-log10(ppm))**2
    dgradT(ndepth+1,ndepth+1) =  -2./ppf(ndepth+1)*(log10(ttp)-log10(ttm))/(log10(ppp)-log10(ppm))**2

    ! dgradT(1,1) = dgradT(1,1) - dgradaddp(1)
    ! dgradT(1,ndepth+1+1) = dgradT(1,ndepth+1+1) - dgradaddt(1)

    ! dgradT(ndepth+1,ndepth+1) = (dgradT(ndepth+1,ndepth+1) - dgradaddp(ndepth+1))
    ! dgradT(ndepth+1,ndepth+1+ndepth+1) =( dgradT(ndepth+1,ndepth+1+ndepth+1) - dgradaddt(ndepth+1))

    !compute the difference between the T-P gradient ans adiabatic gradient
    do idepth = 1, ndepth+1
       if (ppf(idepth)<pp1d_neg_fconv) then
          grad(idepth) = max(gradT(idepth)-gradad(idepth),0.0)
       else
          grad(idepth) = abs(gradT(idepth)-gradad(idepth))
       endif
    enddo

    where( kapfross /= 0. .and. alpha /=0.)

       ! interpolation factor for transparent and opaque bubbles
       theta_loc = (1.0 + yparam*(kapfross*rhof*lconv)**2)/(kapfross*rhof*lconv)

       !compute the coefficient appearing in the computation of Gamma
       coeff = (alpha*rhof*cpspec*theta_loc/&
         & (8.0*(sigsb)*ttf**3))**2*grav*hp*qconv/nuconv*grad

       !compute Gamma
       gam = 0.5*(sqrt(1.0+4.0*coeff)-1.0)

       !compute Delta
       delta = gam/(1.0+gam)*grad

    elsewhere

        !compute Delta
       delta = grad

    end where

    !compute the convective velocity
    do idepth=1,ndepth+1
       if (ppf(idepth)<pp1d_force_vconv) then
          vconv(idepth) = lconv(idepth)/hp(idepth)*sqrt(grav*hp(idepth)*qconv(idepth)*delta(idepth)/nuconv)
       else
          vconv(idepth) = 0.1*sqrt(ppf(idepth)/rhof(idepth))
          delta(idepth) = (gradT(idepth)-gradad(idepth))
       endif
    enddo

    kconv = vconv*lconv

    !compute the convective flux
    fconv = 0.5*rhof*cpspec*ttf*lconv/hp*vconv*delta

    ! !Recalculate convective velocities for Jacobian
    grad = abs(gradT-gradad)

    do idepth=1,ndepth+1
       if (ppf(idepth)<pp1d_force_vconv) then
          vconv(idepth) = lconv(idepth)/hp(idepth)*sqrt(grav*hp(idepth)*qconv(idepth)*grad(idepth)/nuconv)
       else
          vconv(idepth) = 0.1*sqrt(ppf(idepth)/rhof(idepth))
          delta(idepth) = (gradT(idepth)-gradad(idepth))
       endif
    enddo

    do idepth=1,ndepth+1
       if (fconv(idepth)==0.0) vconv(idepth)=eps_vconv*vconv(idepth)
    enddo

    do idepth = 1,ndepth+1
       dfconv(:,ndepth+1+idepth) =  0.5*rhof*cpspec*ttf*lconv/hp*dgradT(:,ndepth+1+idepth)*vconv*3./2.

       dfconv(idepth,ndepth+1+idepth) = dfconv(idepth,ndepth+1+idepth) + 0.5*rhof(idepth)*cpspec(idepth)*lconv(idepth)/hp(idepth)*vconv(idepth)*delta(idepth)+0.25*drhofdt(idepth)*cpspec(idepth)*ttf(idepth)*lconv(idepth)/hp(idepth)*vconv(idepth)*delta(idepth)

       dfconv(:,idepth) =  0.5*rhof*cpspec*ttf*lconv/hp*dgradT(:,idepth)*vconv*3./2.
       dfconv(idepth,idepth) = dfconv(idepth,idepth) +0.25*drhofdp(idepth)*cpspec(idepth)*ttf(idepth)*lconv(idepth)/hp(idepth)*vconv(idepth)*delta(idepth)+0.25*rhof(idepth)*cpspec(idepth)*ttf(idepth)*lconv(idepth)/hp(idepth)*vconv(idepth)*delta(idepth)/ppf(idepth)

    enddo

    do idepth = 1,ndepth+1

       if (ppf(idepth) .le. pstop_conv) then

          fconv(idepth) = 0.
          dfconv(idepth,:) = 0.

       endif

    enddo

    do idepth = 2,ndepth

       if(fconv(idepth-1)==0. .and. fconv(idepth+1)==0.) then

          fconv(idepth) = 0.
          dfconv(idepth,:) = 0.
       endif

    enddo

    ! If eos_type=scvh95 check that none of the convective layers have been calculated using
    ! eos_type=ideal, which occurs when the P, T is outside the eos table boundaries
    if (eos_type=='scvh95') then
      do idepth=1,ndepth+1
        !write(*,*) eos_flag(idepth)
        if ((fconv(idepth)/=0.0) .and. (eos_flag(idepth)=='ideal')) then
           write(*,headrfmt) 'Stop get_convection'
           write(*,intpmfmt) 'idepth',idepth
           write(*,'(A)') 'Convective model level outside P,T range of eos tables'
           write(*,*)
           write(*,separfmt)
          !stop
        endif
      enddo
    endif

  end subroutine get_convection


end module mod_conv
