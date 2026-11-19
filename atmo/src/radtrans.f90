!> \file
!! Contains module mod_radtrans,  \n
!! Routines to compute the radiative transfer
!<

!> hold the parameters of radiative transfer scheme
!<
module mod_radtrans

  integer :: nrays       !< Number of rays used for the radiative transfer


  logical :: rayleigh    !< include rayleigh scattering 
  logical :: scatter     !< Switch to include scattering
  logical :: haze_scatt  !< Switch to include haze scattering
  logical :: print_scat  !< Print the error of scatter iteration
  real    :: eps_scat    !< Accuracy of the scattering solver
  integer :: nitmax_scat !< Maximum number of iteration for the scattering solver
  integer :: nit_scat    !< number of iteration for scattering

  character(256)  :: fspectrum     !< name of the output spectrum
  character(256)  :: ftrans_spec   !< name of the output transmission spectrum
  character(256)  :: firad         !< name of the input irradiation spectrum
  character(256) :: fcfout    !< name of output contribution function file

  real :: rorbit   !< radius of the orbit of the planet
  real :: rstar    !< radius of the star
  real :: rplanet  !< radius of the planet, in km

  integer :: nfreqrad  !< number of frequency points of the irradiation file 
  real    :: nuradmin  !< Minimum frequency of the irradiation file
  real    :: nuradmax  !< Maximum frequency of the irradiation file
  real    :: dnurad    !< resolution of the irradiation file
  real    :: unitnurad !< Frequency unit of the irradiation file
  real    :: unitirad  !< Irradiation unit of the irradiation file

  logical :: irrad       !< Switch to include the external irradiation
  real,dimension(:),allocatable :: Irad !< Irradiation intensity
  real    :: murad       !< Angle of the incoming irradiation
  real    :: ftau_top    !< factor by which the optical depth is multiplied for irrad absorption before entering the domain
  
  real,dimension(:,:),allocatable     :: Inu_surf !< Surface intensity
  real,dimension(:,:,:,:),allocatable :: Inu_surf_um !< Surface intensity for UM input
  real,dimension(:,:,:),allocatable   :: fnu_toa_um !< Top of atmosphere flux from UM input
  real,dimension(:),allocatable       :: xmu_surf !< Gaussian quadrature points used when interpolating in Inu_surf and Inu_surf_um
  
  logical :: trans_spec_emission  !< Switch to take into account the planetary emission when calculating transmission spectra

  real,dimension(2) :: wgl1,xgl1 !< Gauss-Legendre quadrature weights and angles, 1 ray
  real,dimension(:),allocatable :: wgln,xgln !< Gauss-Legendre quadrature weights and angles,n rays
 
  real,dimension(:,:),allocatable :: kr !< Rayleigh Scattering Opacity
  real,dimension(:),allocatable :: kh  !< Haze Scattering Opacity
  
  real :: art_haze		    !< Artificial haze scattering
  real :: art_haze_alpha    !< Exponent to control slope of the haze parameter

  real,dimension(:,:),allocatable :: snu_scat
  real,dimension(:,:,:),allocatable :: snu_scat_2d
  real :: fred !< heat redistribution on the planet
  
  integer :: num_ray_scatt       !< namelist variable that will control number of rayleigh scatterers in the run
 
  contains

    !> Initialize the variables for the radiative transfer scheme
    !<
    subroutine init_radtrans(nukapmin_loc,nukapmax_loc,nfreq_loc)

      use mod_param
      use mod_util
      use mod_grid
      use mod_cst

      implicit none
      include 'netcdf.inc'

      real :: nukapmin_loc,nukapmax_loc,nuloc,dnu_loc
      integer :: stat,idepth,istart,iend,icount,id_file,id_var,irad_start,irad_end,nbin,i,i_loc,ii_loc,nfreq_loc,ifreq,ib,iray
      real,dimension(:),allocatable :: irad_loc
      real,dimension(nraymax) :: alphap,sig
      real,dimension(1,nfreq_loc) :: bnusurf

      namelist /radtrans/nrays,scatter,eps_scat,nitmax_scat,irrad,murad,fspectrum,ftrans_spec,firad,rorbit,rstar,rplanet,nuradmin,nuradmax,dnurad,unitnurad,unitirad,nfreqrad,trans_spec_emission,print_scat,rayleigh,ftau_top,fred,fcfout,num_ray_scatt,&
      haze_scatt,art_haze,art_haze_alpha

      ! default values
      nrays       = 16
      scatter     = .false.
      rayleigh    = .true.
      eps_scat    = 1.0E-4
      nitmax_scat = 1000
      print_scat  = .true.
      ftau_top    = 0.
      fred = 1.

      nuradmin = 0.    ! cm-1
      nuradmax = 5.0E4 ! cm-1
      dnurad   = 0.1   ! cm-1
      nfreqrad = 500000

      num_ray_scatt = 10

      unitnurad = 1. ! file in cgs
      unitirad  = 1. ! file in cgs

      rorbit  = 1.
      rstar   = 1.
      rplanet = 1.

      irrad       = .false.
      murad       = 1.
      trans_spec_emission = .false.

      fspectrum = 'spectrum.ncdf'
      ftrans_spec = 'transmision_spectrum.ncdf'
      fcfout = 'contrib_func.ncdf'
      
      haze_scatt = .false.
      art_haze = 1
      art_haze_alpha = 4
      
      
      ! read namelist in the file fparam
      open(1,file=fparam)
      rewind(1) ; read(1,NML=radtrans,iostat=stat)
       IF (stat/=0) THEN
        BACKSPACE(1)
        READ(1,fmt='(A)') line
        WRITE(*,'(A)') 'Invalid line in namelist: '//trim(line)
        STOP
       END IF
      close(1)
      

      nuradmin = nuradmin*unitnurad/(1./unitL)
      nuradmax = nuradmax*unitnurad/(1./unitL)
      dnurad   = dnurad  *unitnurad/(1./unitL)

      nfreqrad = (nuradmax-nuradmin)/dnurad

      ! convert to code units
      rplanet = rplanet * Rjup
      rstar   = rstar   * rsun 
      rorbit  = rorbit  * au
      


      ! Compute Gauss-Legendre quadrature nodes and weights
      call gauleg(2,-1.0,1.0,wgl1,xgl1)
      allocate(wgln(2*nrays),xgln(2*nrays))
      call gauleg(2*nrays,-1.0,1.0,wgln,xgln)
      
      ! Allocate array to store surface intensity and Gaussian integration
      ! points. Only rays going up.
      allocate(Inu_surf(nrays,nfreq_loc))
      allocate(xmu_surf(nrays))
      
      ! Allocate UM variables if needed
      if (um_input) then
        allocate(Inu_surf_um(nlong,nlat,nrays,nfreq_loc))
      endif


      if (nfreq_loc>0) then
      dnu_loc = (nukapmax_loc - nukapmin_loc)/nfreq_loc
      endif

      if (scatter .and. rayleigh) then
      
         allocate(kr(nfreq_loc,num_ray_scatt)) 
         allocate(kh(nfreq_loc)) 
         kr = 0.
         kh = 0.
         
         do iray=1,num_ray_scatt
         alphap(iray) = (nr(iray)**2-1.)/(4.*pi*nl)
         sig(iray)  = 8.*pi/3.*alphap(iray)**2*(6.+3.*p_fac(iray))/(6.-7.*p_fac(iray))
         end do
      endif


      if (irrad .or. trans_spec_emission) then
         allocate(Irad(nfreq_loc))		
   
         Irad = 0.

         if (corr_k) then

             
            !recall: nurad = nuradmin +(i-1)*dnurad
            istart = floor((nuradmin-nuradmin)/dnurad)+1
            iend   = ceiling((nuradmax-nuradmin)/dnurad)
            
            !if there is only an overlap then take everything you can from the table
            if (iend.gt. nfreqrad) then
               iend = nfreqrad
            else if (istart.lt.1) then
               istart = 1
            endif
            
            icount = iend-istart+1
            
            allocate(Irad_loc(icount))
            Irad_loc = 0.
            call nf(nf_open(firad,nf_nowrite,id_file))
            call nf(nf_inq_varid(id_file,'hnu',id_var))
            call nf(nf_get_vara_double(id_file,id_var,istart,icount,Irad_loc))
            call nf(nf_close(id_file))

            Irad_loc = Irad_loc*unitirad/unitSF
            
!            if (minval(abs(Irad_loc)) == 0.0 ) then
!            
!            	print*, 'Problem with Irad_loc'
!            	
!            endif
            
            do ib = 1,nband

               nuloc = nuradmin+0.5*dnurad
               do ifreq = 1,nfreqrad

                  if(nuloc .lt. nubandmax(ib) .and. nuloc .ge. nubandmin(ib)) then

                     Irad((ib-1)*nkmix+1:ib*nkmix) = Irad((ib-1)*nkmix+1:ib*nkmix)+Irad_loc(ifreq)*dnurad

                     if (scatter .and. rayleigh) then
                        do iray = 1,num_ray_scatt
                           kr((ib-1)*nkmix+1:ib*nkmix,iray) = kr((ib-1)*nkmix+1:ib*nkmix,iray)+Irad_loc(ifreq)*dnurad*sig(iray)*(2.*pi*nuloc)**4
                        enddo

                        if (haze_scatt) then
                           kh((ib-1)*nkmix+1:ib*nkmix) = kh((ib-1)*nkmix+1:ib*nkmix)+Irad_loc(ifreq)*dnurad*(((5.31e-27)/(unitL**2))*(1/(2*amu)))*art_haze*(nuloc/(28571.42857*unitL))**art_haze_alpha !sigma_h2_at350*art_haze*(lambda/350nm)**art_haze_alpha
                        endif                                                                                                                                                                        !350nm = 28571.42857 cm-1 in wavenumber

                     endif

                  endif

                  nuloc = nuloc+dnurad

               enddo

               if(scatter .and. rayleigh) then
                  do iray = 1,num_ray_scatt
                     kr((ib-1)*nkmix+1:ib*nkmix,iray) = kr((ib-1)*nkmix+1:ib*nkmix,iray)/Irad((ib-1)*nkmix+1:ib*nkmix)
                  enddo
                  if (haze_scatt) then
                     kh((ib-1)*nkmix+1:ib*nkmix) = kh((ib-1)*nkmix+1:ib*nkmix)/Irad((ib-1)*nkmix+1:ib*nkmix)
                  endif
               endif

            enddo

         else

            !recall: nurad = nuradmin +(i-1)*dnurad
         istart = floor((nukapmin_loc-nuradmin)/dnurad)+1
         iend   = ceiling((nukapmax_loc-nuradmin)/dnurad)
              
         !if there is only an overlap then take everything you can from the table
         if (iend.gt. nfreqrad) then
            iend = nfreqrad
         else if (istart.lt.1) then
            istart = 1
         endif

         icount = iend-istart+1

         allocate(Irad_loc(icount))

         call nf(nf_open(firad,nf_nowrite,id_file))
         call nf(nf_inq_varid(id_file,'hnu',id_var))
         call nf(nf_get_vara_double(id_file,id_var,istart,icount,Irad_loc))
         call nf(nf_close(id_file))

         !frequency indices to fill with the irradiation file
         irad_start = 1
         irad_end   = nfreq_loc
               
         !deal with overlaps                
         if (nukapmax_loc .ge. nuradmax) then

            irad_end   = nfreq_loc - (nukapmax_loc-nuradmax)/dnu_loc   

         endif

         if (nukapmin_loc .le.nuradmin) then

            irad_start = (nuradmin-nukapmin_loc)/dnu_loc+1

         endif


         if (dnu_loc == dnurad) then
            
            Irad(irad_start:irad_end) = Irad_loc
            
         else if (dnu_loc<dnurad) then
            !increase the resolution of the file by direct injection
            nbin = int(dnurad/dnu_loc)
            
            do i = irad_start,irad_end
               i_loc = int((i-irad_start)/nbin)+1
               Irad(i) = Irad_loc(i_loc)
            enddo
            

         else if(dnu_loc>dnurad) then
            !decrease the resolution of the file by mean
            nbin = int(dnu_loc/dnurad)

            i_loc = 0
            do i = irad_start,irad_end
                                               
               Irad(i) = 0.
               do ii_loc = i_loc+1,i_loc+nbin
                            
                  Irad(i) = Irad(i) +  Irad_loc(ii_loc)
               enddo
                         
               Irad(i) = Irad(i)/nbin
   
               i_loc = i_loc + nbin
                         
            enddo

            
         endif

         Irad = Irad*unitirad/unitSF

         if (scatter .and. rayleigh) then

            nuloc = nukapmin_loc+0.5*dnu_loc
            do ifreq = 1,nfreq_loc

               do iray = 1,num_ray_scatt
                  kr(ifreq,iray) = sig(iray)*(2.*pi*nuloc)**4
               enddo
               if (haze_scatt) then
                  kh(ifreq) = (((5.31e-27)/(unitL**2))*(1/(2*amu)))*art_haze*(nuloc/(28571.42857*unitL))**art_haze_alpha
               endif

               nuloc = nuloc + dnu_loc

            enddo

         endif
      endif

   endif

   if ((.not. irrad) .and. scatter .and. rayleigh) then

     kr = 0.
     kh = 0.

      if (corr_k) then

         do ib = 1,nband

            dnu_loc = (nubandmax(ib)- nubandmin(ib))/nfreq_loc
            nuloc = nubandmin(ib)+0.5*dnu_loc

            do ifreq = 1,nfreq_loc

               do iray = 1,num_ray_scatt
                  kr ((ib-1)*nkmix+1:ib*nkmix,iray) = kr((ib-1)*nkmix+1:ib*nkmix,iray)+dnu_loc*sig(iray)*(2.*pi*nuloc)**4
               enddo

               if (haze_scatt) then
                  kh ((ib-1)*nkmix+1:ib*nkmix) = kh((ib-1)*nkmix+1:ib*nkmix)+dnu_loc*(((5.31e-27)/(unitL**2))*(1/(2*amu)))*art_haze*(nuloc/(28571.42857*unitL))**art_haze_alpha
               endif
               nuloc = nuloc+dnu_loc

            enddo

            do iray = 1,num_ray_scatt
               kr ((ib-1)*nkmix+1:ib*nkmix,iray) = kr((ib-1)*nkmix+1:ib*nkmix,iray)/(nubandmax(ib)- nubandmin(ib))
            enddo

            if (haze_scatt) then
             kh ((ib-1)*nkmix+1:ib*nkmix) = kh((ib-1)*nkmix+1:ib*nkmix)/(nubandmax(ib)- nubandmin(ib))
            endif

        enddo

         else

            nuloc = nukapmin_loc+0.5*dnu_loc

            do ifreq = 1,nfreq_loc

               do iray = 1,num_ray_scatt
                  kr (ifreq,iray) = sig(iray)*(2.*pi*nuloc)**4
               enddo

               if (haze_scatt) then
                  kh (ifreq) = (((5.31e-27)/(unitL**2))*(1/(2*amu)))*art_haze*(nuloc/(28571.42857*unitL))**art_haze_alpha
               endif
               nuloc = nuloc+dnu_loc
            enddo

      endif
   endif

      if (mype == cputerm) then
         write(*,headrfmt) 'Radiative transfer'
         write(*,intpmfmt) 'nrays',nrays,'Number of rays for Gauss-Legendre quadrature'
         write(*,logpmfmt) 'scatter',scatter,'if True includes scattering'
         if (scatter) then
            write(*,exppmfmt) 'eps_scat',eps_scat,'Accuracy of the scattering solver'
            write(*,intpmfmt) 'nitmax_scat',nitmax_scat,'Maximum number of iterations for the scattering solver'
         endif
         
         write(*,logpmfmt) 'irrad',irrad,'if True includes irradiation'
         if (irrad) then
            write(*,exppmfmt) 'murad',murad,'Angle for incoming irradiation'
            write(*,exppmfmt) 'rstar [Rsun]',rstar/Rsun,'star radius'
            write(*,exppmfmt) 'rorbit  [au]',rorbit/au,'planet orbit'
            write(*,exppmfmt) 'ftau_top',ftau_top,'scaling factor between tau=0 and tau=taumin'
         endif
         write(*,*) 
         write(*,separfmt)
      endif

      
    end subroutine init_radtrans

    !> Compute the star flux
    !< 
    subroutine startrans(nfreq_loc,kapnu_loc,kapstd_loc,fstar_loc,Irad_loc)

      use mod_grid, only : ndepth,taufstd
      use mod_cst

      implicit none

      integer :: nfreq_loc,idepth
      real,dimension(ndepth) :: kapstd_loc
      real,dimension(ndepth,nfreq_loc) :: kapnu_loc
      real,dimension(ndepth+1,nfreq_loc) :: fstar_loc
      real,dimension(nfreq_loc) :: Irad_loc,Iloop

      Iloop = Irad_loc*rstar**2/rorbit**2*4.*pi*fred

      Iloop(:) = Iloop(:)*exp(-ftau_top*(taufstd(1))*kapnu_loc(1,:)/kapstd_loc(1)/murad)

      fstar_loc(1,:) =  - Iloop(:)*murad

      do idepth = 2,ndepth+1
			
         Iloop(:) = Iloop(:)*exp(-(taufstd(idepth)-taufstd(idepth-1))*kapnu_loc(idepth-1,:)/kapstd_loc(idepth-1)/murad)
               
         fstar_loc(idepth,:) =  - Iloop(:)*murad     
        
      enddo
      
    end subroutine startrans
    


    !> Compute the radiative transfer. The number of rays and the source term are\n
    !! in argument because we will call the routine with others 1 ray and \n
    !! other source terms to compute the Jacobi matrix
    !<
    subroutine radtrans(nrays_loc,scatter_loc,nfreq_loc,kapnu_loc,kapstd_loc,epsfnu_loc,snu,bnu,hnup,hnud,jnup,jnud,snu_scat_loc)

      use mod_grid
      use mod_cst
      use mod_param

      implicit none

      character :: CR = CHAR(13)    ! carriage return character
      character(4)  :: fnit
      real :: diff_scat
      logical :: scatter_loc
      integer :: nrays_loc,idepth,l,mmu,nfreq_loc,ifreq
      real,dimension(ndepth+1,nfreq_loc) :: Snu,hnup,jnup,hnud,jnud,snu_tot,epsfnu_loc,bnu,snu_scat_loc !source term and moments
      real,dimension(:,:),allocatable   :: Inu !spectral intensity 
      real,dimension(nfreq_loc)          :: psi1,psi2,psi3,expmtau,dtaum,dtaup,tau !integral csts
            
      real,dimension(:),allocatable     :: xmu,wmu,dsnu
      real,dimension(ndepth) :: kapstd_loc,ntot
      real,dimension(ndepth,nfreq_loc) :: kapnu_loc 
      real,dimension(:,:  ),allocatable :: lambda
      real,dimension(:,:,:),allocatable :: psi1_scat,psi2_scat,psi3_scat,expmtau_scat

      !source function initialisation at planck (or planck + solar scatter if scatter)
      snu_tot = snu
      
      mmu = 2*nrays_loc
      allocate(xmu(mmu),wmu(mmu))

      if (scatter_loc) then
         allocate(psi1_scat(mmu,ndepth+1,nfreq_loc),psi2_scat(mmu,ndepth+1,nfreq_loc),psi3_scat(mmu,ndepth+1,nfreq_loc),expmtau_scat(mmu,ndepth+1,nfreq_loc))
         allocate(lambda(ndepth+1,nfreq_loc))
         allocate(dsnu(nfreq_loc))

         psi1_scat = 0. ; psi2_scat = 0. ; psi3_scat = 0.
         lambda = 0.
         dsnu   = 0.
        

      endif

       nit_scat  = 0
		
      ! set quadrature nodes and weights
      if (nrays_loc ==1) then
         xmu = xgl1
         wmu = wgl1*0.5
      else
         xmu = xgln
         wmu = wgln*0.5
      endif
   		
      !allocate spectra intensity
      allocate(Inu(mmu,nfreq_loc))
      
      diff_scat = 1.
      
      

      if (scatter_loc) then
         snu_tot = snu_scat_loc
      else
         snu_tot = snu
      endif

       
      if (mype == cputerm .and. print_scat .and. scatter_loc) then
         write(*,headrfmt) 'Scattering '
      endif

   do while(((diff_scat .gt. eps_scat) .and. (nit_scat<nitmax_scat)) .or. (nit_scat==1)) 

      if (scatter_loc) then
         nit_Scat = nit_scat + 1
         
      endif

      Hnup  = 0.
      Hnud  = 0.
      Jnup  = 0.
      Jnud  = 0.
      Inu   = 0.
      

      diff_scat = 0.

     

      !1. perform an integration from top to bottom for rays going gown 
      !mu from 1 to mmu/2, compute spectral intensity at interfaces

      !boundary condition at the top of the atmosphere
      do l = 1,mmu/2 

         Inu(l,:) = 0.

         !compute the mean intensity
         Jnud(1,:) = Jnud(1,:) + Inu(l,:)*wmu(l)
         
         !compute the flux
         Hnud(1,:) = Hnud(1,:) + Inu(l,:)*wmu(l)*xmu(l)
         
      enddo

        
      !interior points
      do idepth = 2, ndepth+1
        
         do l = 1,mmu/2 
            
            dtaum = 0. ; dtaup = 0.
             
            if (kapstd_loc(idepth-1) /= 0.) then
               dtaum(:) = (taufstd(idepth)-taufstd(idepth-1))* &
                    & kapnu_loc(idepth-1,:)/kapstd_loc(idepth-1)/abs(xmu(l))
            endif
            if (idepth == ndepth + 1) then
               dtaup(:) = 0.
            else
           if (kapstd_loc(idepth) /= 0.) then
                  dtaup(:) = (taufstd(idepth+1)-taufstd(idepth))* &
                       & kapnu_loc(idepth,:)/kapstd_loc(idepth)/abs(xmu(l))
               endif
            endif
            
            if (scatter_loc) then

               if (nit_scat == 1) then
                  call integral_cst(nfreq_loc,dtaum,dtaup,psi1,psi2,psi3,expmtau)          
                  psi1_scat(l,idepth,:)    = psi1(:)
                  psi2_scat(l,idepth,:)    = psi2(:)
                  psi3_scat(l,idepth,:)    = psi3(:)
                  expmtau_scat(l,idepth,:) = expmtau(:)
                
                  lambda(idepth,:) = lambda(idepth,:) + psi2_scat(l,idepth,:)*wmu(l)

                  if (idepth>2) then
                     lambda(idepth,:) = lambda(idepth,:) + psi3_scat(l,idepth-1,:)*expmtau_scat(l,idepth,:)*wmu(l)
                  endif

               else

                  psi1(:)    = psi1_scat(l,idepth,:)
                  psi2(:)    = psi2_scat(l,idepth,:)
                  psi3(:)    = psi3_scat(l,idepth,:)
                  expmtau(:) = expmtau_scat(l,idepth,:) 
               endif
            else
               
               call integral_cst(nfreq_loc,dtaum,dtaup,psi1,psi2,psi3,expmtau)

            endif
            
            !compute spectral intensity
            if (idepth == ndepth+1) then
               Inu(l,:) = Inu(l,:)*expmtau(:) + psi1(:)*snu_tot(idepth-1,:) &
                    & + psi2(:)*snu_tot(idepth,:)
            else
              
               Inu(l,:) = Inu(l,:)*expmtau(:) + psi1(:)*snu_tot(idepth-1,:) &
                    & + psi2(:)*snu_tot(idepth,:) + psi3(:)*snu_tot(idepth+1,:)
                       
            endif

            !compute the mean intensity
            Jnud(idepth,:) = Jnud(idepth,:) + Inu(l,:)*wmu(l)

            !compute the flux
            Hnud(idepth,:) = Hnud(idepth,:) + Inu(l,:)*wmu(l)*xmu(l)
            
         enddo
                   
      enddo
      
      !2. perform an integration from bottom to top for rays going up 
      !mu from mmu/2+1 to mmu, compute spectral intensity at interfaces
      
      !boundary condition at the bottom, assume diffusion approximation
      do l = mmu/2+1,mmu 

         where (kapnu_loc(ndepth,:)/=0.)

            Inu(l,:) = bnu(ndepth+1,:)+(bnu(ndepth+1,:)-bnu(ndepth,:))/(taufstd(ndepth+1)-taufstd(ndepth))*abs(xmu(l))/kapnu_loc(ndepth,:)*kapstd_loc(ndepth)
!            Inu(l,:) = snu_tot(ndepth+1,:)+(snu_tot(ndepth+1,:)-snu_tot(ndepth,:))/(taufstd(ndepth+1)-taufstd(ndepth))*abs(xmu(l))/kapnu_loc(ndepth,:)*kapstd_loc(ndepth)
            
         elsewhere

            Inu(l,:) = bnu(ndepth+1,:)
!            Inu(l,:) = snu_tot(ndepth+1,:)

         end where
     
         !compute the mean intensity
         Jnup(ndepth+1,:) = Jnup(ndepth+1,:) + Inu(l,:)*wmu(l)
         
         !compute the flux
         Hnup(ndepth+1,:) = Hnup(ndepth+1,:) + Inu(l,:)*wmu(l)*xmu(l)
         
        
      enddo
      
      !correct the first point of the source function
      if (scatter_loc) then

         !compute the source function correction
         dsnu(:) = (1.0-epsfnu_loc(ndepth+1,:))*(jnup(ndepth+1,:)+jnud(ndepth+1,:))+snu(ndepth+1,:)-snu_tot(ndepth+1,:)

         !apply convergence acceleration
         dsnu(:) = dsnu(:)/(1.0-(1.0-epsfnu_loc(ndepth+1,:))*lambda(ndepth+1,:))
         
         !apply source function correction
         snu_tot(ndepth+1,:) = snu_tot(ndepth+1,:) + dsnu(:)

         
         do ifreq = 1,nfreq_loc
            if(snu_tot(ndepth+1,ifreq)>0. .and. snu(ndepth+1,ifreq)>0.) then
               diff_scat = max(diff_scat,abs(dsnu(ifreq)/snu_tot(ndepth+1,ifreq)))
            endif
         enddo
       

      endif
      
      
      !interior points
      do idepth = ndepth,1,-1
        
         do l = mmu/2+1,mmu 
            
            dtaum = 0. ; dtaup = 0.

            if (kapstd_loc(idepth) /= 0.) then
               dtaum(:) = (taufstd(idepth+1)-taufstd(idepth))* &
                    & kapnu_loc(idepth,:)/kapstd_loc(idepth)/abs(xmu(l))
            endif
            if (idepth == 1) then
               dtaup(:) = 0.
            else
              
               if (kapstd_loc(idepth-1) /= 0.) then
                  dtaup(:) = (taufstd(idepth)-taufstd(idepth-1))* &
                       & kapnu_loc(idepth-1,:)/kapstd_loc(idepth-1)/abs(xmu(l))
               endif

            endif
            
            if(scatter_loc) then
            
               if (nit_scat == 1) then
                  call integral_cst(nfreq_loc,dtaum,dtaup,psi1,psi2,psi3,expmtau)
                   
                  psi1_scat(l,idepth,:) = psi1(:)
                  psi2_scat(l,idepth,:) = psi2(:)
                  psi3_scat(l,idepth,:) = psi3(:)
                  expmtau_scat(l,idepth,:) = expmtau(:)

                  lambda(idepth,:) = lambda(idepth,:) + psi2_scat(l,idepth,:)*wmu(l)

                  if (idepth<ndepth) then
                     lambda(idepth,:) = lambda(idepth,:) + psi3_scat(l,idepth+1,:)*expmtau_scat(l,idepth,:)*wmu(l)
                  endif

               else

                  psi1(:) = psi1_scat(l,idepth,:)
                  psi2(:) = psi2_scat(l,idepth,:)
                  psi3(:) = psi3_scat(l,idepth,:)
                  expmtau(:) = expmtau_scat(l,idepth,:)
                  
               endif
               
            else
               
               call integral_cst(nfreq_loc,dtaum,dtaup,psi1,psi2,psi3,expmtau)

            endif

            

            !compute spectral intensity
            if (idepth == 1) then
               Inu(l,:) = Inu(l,:)*expmtau(:) + psi1(:)*snu_tot(idepth+1,:) &
                    & + psi2(:)*snu_tot(idepth,:)
            else
               Inu(l,:) = Inu(l,:)*expmtau(:) + psi1(:)*snu_tot(idepth+1,:) &
                    & + psi2(:)*snu_tot(idepth,:) + psi3(:)*snu_tot(idepth-1,:)
            endif
            
            !compute the mean intensity
            Jnup(idepth,:) = Jnup(idepth,:) + Inu(l,:)*wmu(l)

            !compute the flux
            Hnup(idepth,:) = Hnup(idepth,:) + Inu(l,:)*wmu(l)*xmu(l)
          
         enddo
         
         !scattering corrections
         if(scatter_loc) then
            
            !post correct the ingoing mean intensity
            do l = 1,mmu/2

               jnud(idepth,:) = jnud(idepth,:) + dsnu(:)*psi3_scat(l,idepth,:)*wmu(l)

            enddo
            
            ! !compute the source function correction
             dsnu(:) = (1.0-epsfnu_loc(idepth,:))*(jnup(idepth,:)+jnud(idepth,:))+snu(idepth,:)-snu_tot(idepth,:) 

           
            !apply convergence acceleration
            dsnu(:) = dsnu(:)/(1.0-(1.0-epsfnu_loc(idepth,:))*lambda(idepth,:))

            !apply source function correction
            snu_tot(idepth,:) = snu_tot(idepth,:) + dsnu(:)

 
            do l = mmu/2+1,mmu
               Inu(l,:)  = Inu(l,:) + dsnu(:)*psi2(:)
            enddo
            
            do ifreq = 1,nfreq_loc
            
               if(snu_tot(idepth,ifreq)>0. .and. snu(idepth,ifreq)>0.) then
        
                  diff_scat = max(diff_scat,abs(dsnu(ifreq)/snu_tot(idepth,ifreq)))
                  
               endif
            enddo
          
         endif

  

      enddo
     
      
      !write(*,headrfmt,advance='no') 'Scattering '
      if (mype == cputerm .and. print_scat .and. scatter_loc) then
         write(*,'(A)',advance='no') CR
         write(*,1000,advance='no') nit_scat,diff_scat
         
      endif

1000        format(5x,'Iter. scat.: ',i5,5x,' dSnu/Snu: ',es10.1,5x)
           

   enddo
   
   ! Save intensity and Gaussian quadrature points at surface
   if (nfreq_loc == size(Inu_surf,2)) then
     Inu_surf = Inu(mmu/2+1:mmu,:)
     xmu_surf = xmu(mmu/2+1:mmu)
   endif

   if (scatter_loc .and. nit_scat>1) then

      snu_scat_loc = snu_tot

   endif

   if (mype == cputerm .and. print_scat .and. scatter_loc) then
      write(*,*)
      write(*,*)
      write(*,separfmt)
   endif
      
   
    
      
    end subroutine radtrans

    !> compute the integral constant of the short characteristic method to compute the radiative transfer
    !<
    subroutine integral_cst(nfreq_loc,dtaum,dtaup,psi1,psi2,psi3,expmtau)

      use mod_param
      use mod_grid

      
      implicit none

      integer      :: i,nfreq_loc
      integer      :: idepth !indice of the interface

      real,dimension(nfreq_loc) :: U0,U1,U2,dtaum,dtaup
      real,dimension(nfreq_loc) :: psi1,psi2,psi3,expmtau
      real :: dtau

      !prevent an overflow of the exponential function
      expmtau(:) = exp(-dtaum(:))

    !  dtau = 1.0E-6

    !  write(*,*) abs((dtau)**2-2.0*(dtau-(1.-exp(-dtau)))-(dtau)**3/3.)/((dtau)**3/3.),abs((dtau-(1.-exp(-dtau)))-(dtau)**2/2.)/((dtau)**2/2.),abs((1.-exp(-dtau))-(dtau))/((dtau))

      where (dtaum < 1.0E-9)
         U0 = dtaum
      elsewhere
         U0 = 1.0 - expmtau
      endwhere

      where (dtaum < 1.0E-6)
         U1 = dtaum**2/2.
      elsewhere
         U1 = dtaum - U0
      endwhere

      where (dtaum < 1.0E-3)
         U2 = dtaum**3/3.
      elsewhere
         U2 = dtaum**2 - 2.0*U1
      endwhere
    
      !if dtaum == 0., the medium is transparent
      where(dtaum == 0.)
 
         psi1 = 0.
         psi2 = 0.
         psi3 = 0.
         
      !if dtaup == 0., use linear interpolation
      elsewhere(dtaup == 0.) 

         psi1 = U0-U1/dtaum
         psi2 = U0 - psi1
         psi3 = 0.
        
      !else use parabolic interpolation

      elsewhere(U2-dtaum*U1 < 0.)

         psi1 = U0-U1/dtaum
         psi2 = U0 - psi1
         psi3 = 0.
         
      elsewhere

         psi1 = U0 + (U2-(2.0*dtaum+dtaup)*U1)/(dtaum*(dtaum+dtaup))
         psi2 = ((dtaup+dtaum)*U1-U2)/(dtaum*dtaup)
         psi3 = (U2-dtaum*U1)/(dtaup*(dtaum+dtaup))

       
            
      end where

      ! Check NaN values in psi1, psi2 and psi3
      if ((isnan(sum(psi1))) .or. (isnan(sum(psi2))) .or. (isnan(sum(psi3)))) then 
        write(*,*) 'error radtrans'
        write(*,*) 'NaN values in integral_cst subroutine'
        stop 
      endif


         if(sum(psi1+psi2+psi3+expmtau)<nfreq_loc*0.97 .or. sum(psi1+psi2+psi3+expmtau)>nfreq_loc*1.3) then
         
          !  write(*,*)mype, sum(psi1+psi2+psi3+expmtau)/nfreq_loc
            if (mype==0)then
               write(*,*) 'error radtrans'
               do i =1,nfreq_loc
                  write(*,*) i,dtaum(i),dtaup(i),psi1(i),psi2(i),psi3(i),expmtau(i),psi1(i)+psi2(i)+psi3(i)+expmtau(i)
        
               enddo
               stop
            endif
      endif
    

    end subroutine integral_cst
      

    !> compute the planck function at cell center or interfaces\n
    !! suppose that bnu is allocated accordingly !
    !<
    function get_bnu(ndepth_loc,nfreq_loc,tt_loc)

      use mod_grid
      use mod_param
      use mod_cst
      use mod_util

      implicit none
      
      integer :: idepth,ndepth_loc,nfreq_loc,ib,ifreq,n_points,i
      real,dimension(ndepth_loc,nfreq_loc) :: get_bnu
      real,dimension(ndepth_loc,nband) :: int_varmax, int_varmin
      real,dimension(ndepth_loc)           :: tt_loc
      real,dimension(nfreq_loc)            :: hcnu_kt
      real,dimension(nfreq_loc)            :: nu_int
      real,dimension(nfreq_loc)            :: bp_loc
      real :: dnu_int, d_var, d_var1, int_bnu
      real, allocatable :: integrand(:), int_var(:)

      if (.not. corr_k) then
      do idepth = 1,ndepth_loc

         hcnu_kt(:) = hplanck*cvel*nu(:)/(kb*tt_loc(idepth))
         
         where (hcnu_kt < 1.0E-4) !use rayleigh-jeans
            
            get_bnu(idepth,:) = 2.0*hplanck*cvel**2*nu(:)**3/hcnu_kt(:)
            
         elsewhere

            get_bnu(idepth,:) = 2.0*hplanck*cvel**2*nu(:)**3 &
                 & / (exp(hplanck*cvel*nu(:)/(kb*tt_loc(idepth)))-1.)

         endwhere

      enddo
   else
         
        
         
!       do ib = 1,nband
!          
!             dnu_int = (nubandmax(ib)-nubandmin(ib))/nfreq_loc
!             nu_int(1) = nubandmin(ib)+0.5*dnu_int
!             do ifreq = 2,nfreq_loc
!                nu_int(ifreq) = nu_int(ifreq-1)+dnu_int
!             enddo
! 
!             do idepth = 1,ndepth_loc
!                
!                hcnu_kt(:) = hplanck*cvel*nu_int(:)/(kb*tt_loc(idepth))
! 
!                where (hcnu_kt < 1.0E-4) !use rayleigh-jeans
!                   
!                   bp_loc = 2.0*hplanck*cvel**2*nu_int(:)**3/hcnu_kt(:)
!                   
!                elsewhere
!                   
!                   bp_loc = 2.0*hplanck*cvel**2*nu_int(:)**3 &
!                         / (exp(hplanck*cvel*nu_int(:)/(kb*tt_loc(idepth)))-1.)
!                   
!                endwhere
!                
!                get_bnu(idepth,(ib-1)*nkmix+1:ib*nkmix) = sum(bp_loc)*dnu_int
!                
!             enddo
!             
!          enddo
! 
  do ib = 1,nband  
        
        d_var = 1e-3
        
      do idepth = 1,ndepth_loc        
      
           int_varmax(idepth,ib) = (hplanck*cvel*nubandmax(ib))/(kb*tt_loc(idepth))
           int_varmin(idepth,ib) = (hplanck*cvel*nubandmin(ib))/(kb*tt_loc(idepth))
           
           n_points = nint((int_varmax(idepth,ib) - int_varmin(idepth,ib))/d_var) + 1 
           d_var1 = (int_varmax(idepth,ib) - int_varmin(idepth,ib))/(n_points-1)
           
       !! Calculating integrand
       allocate(int_var(1:n_points))
       allocate(integrand(1:n_points))
           
       do i = 1,n_points
           int_var(i) = int_varmin(idepth,ib) + d_var1*(i-1)
           integrand(i) = (int_var(i)**3)/(exp(int_var(i))-1.0)
       enddo
   
       !integrate using trapezoidal rule
       int_bnu = sum((integrand(1:n_points-1) + integrand(2:n_points))* &
  		d_var1/2.0)
  		
       get_bnu(idepth,(ib-1)*nkmix+1:ib*nkmix) = ((2.0*(kb**4)*(tt_loc(idepth)**4))/((hplanck**3)*(cvel**2)))*int_bnu
       
       deallocate(int_var)
       deallocate(integrand)
       
       enddo
    enddo
         
         
      endif

    end function get_bnu


    !> compute the derivative of planck function at cell center or interfaces\n
    !! suppose that dbnu is allocated accordingly !
    !<
    function get_dbnu(ndepth_loc,nfreq_loc,tt_loc)

      use mod_grid
      use mod_cst

      implicit none
     
      integer :: ndepth_loc,idepth,nfreq_loc
      real,dimension(ndepth_loc,nfreq_loc) :: get_dbnu,bnu_loc1,bnu_loc2
      real,dimension(ndepth_loc)          :: tt_loc
     
      real :: deps 

      deps = 1.0E-6

      bnu_loc1 = get_bnu(ndepth_loc,nfreq_loc,tt_loc*(1.0+deps))
      bnu_loc2 = get_bnu(ndepth_loc,nfreq_loc,tt_loc)

      do idepth = 1,ndepth_loc
         get_dbnu(idepth,:) = (bnu_loc1(idepth,:) - bnu_loc2(idepth,:))/(deps*tt_loc(idepth))
      enddo

     end function get_dbnu

    !> compute the second derivative of planck function at cell center or interfaces\n
    !! suppose that receiving array is allocated accordingly !
    !<
    function get_d2bnu(ndepth_loc,nfreq_loc,tt_loc)

      use mod_grid

      implicit none
     
      integer :: ndepth_loc,idepth,nfreq_loc
      real,dimension(ndepth_loc,nfreq_loc) :: get_d2bnu,dbnu_loc1,dbnu_loc2
      real,dimension(ndepth_loc)          :: tt_loc
     
      real :: deps 

      deps = 1.0E-6

      dbnu_loc1 = get_dbnu(ndepth_loc,nfreq_loc,tt_loc*(1.0+deps))
      dbnu_loc2 = get_dbnu(ndepth_loc,nfreq_loc,tt_loc)
   
       
      do idepth = 1,ndepth_loc
         get_d2bnu(idepth,:) = (dbnu_loc1(idepth,:) - dbnu_loc2(idepth,:))/(deps*tt_loc(idepth))

      enddo


    end function get_d2bnu
    
    !> Compute the transmission spectrum of the planet using a 1D P-T profile
    !<
    subroutine calc_trans_radius_1d(nfreq_loc,Rp_transit)
    
      use mod_cst
      use mod_param
      use mod_util
      use mod_grid
      
      implicit none
      
      ! Number of frequencies
      integer,intent(in) :: nfreq_loc
      ! Apparent planetary radius
      real,dimension(:),intent(out) :: Rp_transit

      double precision :: impact_param_1,impact_param_2,step_impact_param_1,step_impact_param_2
      ! Current interpolated densities
      real :: rho_cur_1, rho_cur_2
      ! Current calculated optical depth along path
      real,dimension(:),allocatable :: tau_path_1, tau_path_2
      ! Current interpolated opacities
      real,dimension(:),allocatable :: kapnu_cur_1, kapnu_cur_2
      ! Path along ray towards the observer and step for this path
      real :: path, step_path
      ! Equivalent radius for current path
      real :: rr_cur
      ! Loop indices
      integer :: idepth, idepth_ip, i
      ! True if top of atmosphere has been reached
      logical :: exit_path_loop
      
      allocate(tau_path_1(nfreq_loc),tau_path_2(nfreq_loc))
      allocate(kapnu_cur_1(nfreq_loc), kapnu_cur_2(nfreq_loc))
      
!      kapfnu = kapfnu + 0.012/(unitL**2/unitM)
!      do idepth = 1,ndepth+1
!        kapfnu(idepth,:) = kapfnu(idepth,:) + 0.1*(1./(nu/unitL*100)/(0.8e-6))**-4
!      end do
!      do idepth = 1,ndepth+1
!        kapfnu(idepth,:) = 0.1*(1./(nu/unitL*100)/(0.8e-6))**-4
!      end do
      
      Rp_transit = 0.
      impact_param_2 = rrf(ndepth+1)
      step_impact_param_2 = 0.
      idepth = ndepth+1
      i = 0

      do
        
        ! Update impact parameter
        impact_param_2 = impact_param_2 + step_impact_param_2
        
        ! If top of the atmosphere has been reached exit loop
        if (impact_param_2 >= rrf(1)) then
          exit
          
        ! Otherwise increment radius index if necessary
        else if (impact_param_2 >= rrf(idepth-1)) then
          idepth = idepth-1
        endif
        
        ! Update step in impact parameter
        step_impact_param_2 = rrf(idepth-1) - rrf(idepth)

        path = 0.
        step_path = rrf(idepth-1) - rrf(idepth)
        
        ! Get opacity and density at path = 0
        kapnu_cur_1(:) = kapfnu(idepth,:)
        rho_cur_1 = rhof(idepth)
        
        ! Reset lower bound height index for interpolation and logical for
        ! exiting loop
        idepth_ip = idepth
        exit_path_loop = .false.
        
        ! Reset optical depth along path
        tau_path_2 = 0.
        
        do
          ! Get radius corresponding to this value of the path
          path = path + step_path
          rr_cur = sqrt(path**2 + impact_param_2**2)
          
          ! If top of the atmosphere has been reached, add remaining bit and
          ! exit loop.
          if (rr_cur >= rrf(1)) then
            rr_cur = rrf(1)
            step_path = sqrt(rr_cur**2 - impact_param_2**2) - (path - step_path)
            path = sqrt(rr_cur**2 - impact_param_2**2)
            exit_path_loop = .true.
            
          ! Otherwise increment radius index if necessary
          else if (rr_cur > rrf(idepth_ip-1)) then
            idepth_ip = idepth_ip-1
          endif
          
          ! Interpolate to get value for opacity and density at this radius
          kapnu_cur_2(:) = kapfnu(idepth_ip-1,:) + &
              (kapfnu(idepth_ip,:) - kapfnu(idepth_ip-1,:))/ &
              (rrf(idepth_ip) - rrf(idepth_ip-1))* &
              (rr_cur - rrf(idepth_ip-1))
          rho_cur_2 = rhof(idepth_ip-1) + &
              (rhof(idepth_ip) - rhof(idepth_ip-1))/ &
              (rrf(idepth_ip) - rrf(idepth_ip-1))* &
              (rr_cur - rrf(idepth_ip-1))
          
          ! Add the contribution to the integral using the trapezoidal rule
          tau_path_2 = tau_path_2 + (kapnu_cur_1*rho_cur_1 + kapnu_cur_2*rho_cur_2)*step_path
          
          ! Exit loop if done
          if (exit_path_loop) exit
          
          ! Prepare for next iteration by saving newly calculated opacities and densities
          kapnu_cur_1 = kapnu_cur_2
          rho_cur_1 = rho_cur_2
          
        enddo
        
        ! Using symmetry and take into account -path_min to 0 in addition to
        ! 0 to path_max we must multiply tau_path_2 by 2. Using the trapezoidal
        ! rule we must divide tau_path_2 by 2. All in all, we need not modify
        ! tau_path_2.
          
        ! If not first iteration, add contribution to transit radius integral
        if (i /= 0) then
          Rp_transit = Rp_transit + (impact_param_1*(1 - exp(-tau_path_1)) + &
              impact_param_2*(1 - exp(-tau_path_2)))*step_impact_param_1
        endif
        
        ! Store old quantities
        tau_path_1 = tau_path_2
        impact_param_1 = impact_param_2
        step_impact_param_1 = step_impact_param_2
        
        ! Update loop index
        i = i+1

      enddo
      
      ! For idepth = 1, the optical depth is zero. Add this contribution to
      ! Rp_transit:
      tau_path_2 = 0
      impact_param_2 = impact_param_2 + step_impact_param_2
      Rp_transit = Rp_transit + (impact_param_1*(1 - exp(-tau_path_1)) + &
          impact_param_2*(1 - exp(-tau_path_2)))*step_impact_param_1
      
      ! Again, the integral calculating Rp_transit should be multiplied by 2,
      ! but due to the trapezoidal rule Rp_transit should also be divided by 2.
      ! Add contribution from opaque atmosphere below lower boundary and take
      ! the square root.
      Rp_transit = sqrt(rrf(ndepth+1)**2 + Rp_transit)
      
      ! Add contribution from thermal emission from planet. Ignore planetary
      ! and stellar limb darkening/brightening here for now. Since observers
      ! neglect thermal planet emission and remove stellar limb darkening, it
      ! is not going to be correct anyway.
      if (trans_spec_emission) then
        Rp_transit = Rp_transit/sqrt(1 + (2*pi*Rp**2*fnu(1,:))/ &
            (4*pi*rstar**2*irad(:)))
      endif

      deallocate(tau_path_1,tau_path_2)
      deallocate(kapnu_cur_1, kapnu_cur_2)
      
    end subroutine calc_trans_radius_1d
    
    !> Calculate transit radius from UM output.
    !<
    subroutine calc_trans_radius_3d(nfreq_loc,Rp_transit)
    
      use mod_grid
      use mod_cst
      use mod_param
      use mod_util
      
      implicit none
      
      ! Number of frequencies
      integer,intent(in) :: nfreq_loc
      ! Apparent planetary radius
      real,dimension(:),intent(out) :: Rp_transit
      ! Cylindrical coordinates. rho is equivalent to the impact parameter.
      real,dimension(:),allocatable :: rho_cyl, phi_cyl
      ! Optical depth along path
      real,dimension(:,:,:),allocatable :: tau_path
      ! Number of cylindrical points
      integer :: nrho_cyl, nphi_cyl
      ! Integral over phi
      real,dimension(:,:),allocatable :: int_phi
      ! Planet emission spectrum
      real,dimension(:),allocatable :: fnu_surf_spec
      ! Loop indices
      integer :: irho,iphi,inu
      
      ! Construct cylindrical rho grid for integration
      nrho_cyl = ndepth+1
      allocate(rho_cyl(nrho_cyl))
      call linspace(rrf(ndepth+1), rrf(1), nrho_cyl, rho_cyl)
      
      ! Construct cylindrical phi grid
      nphi_cyl = nlong
      allocate(phi_cyl(nphi_cyl))
      call linspace(longitude_um(1), longitude_um(nlong), nphi_cyl, &
          phi_cyl)
      
      ! Loop over cylindrical grid and calculate optical depths
      allocate(tau_path(nrho_cyl, nphi_cyl, nfreq_loc))
      allocate(int_phi(nrho_cyl, nfreq_loc))
      do irho=1,nrho_cyl
        do iphi=1,nphi_cyl
          call calc_tau_path(nfreq_loc, rho_cyl(irho), phi_cyl(iphi), tau_path(irho,iphi,:))
        end do
        
        ! Integrate over phi
        do inu=1,nfreq_loc
          call integrate_trapezoid_loop(phi_cyl, 1 - exp(-tau_path(irho,:,inu)), &
              int_phi(irho,inu))
        end do
      end do
      
      ! Integrate over rho (impact factor)
      do inu=1,nfreq_loc
        call integrate_trapezoid(rho_cyl, rho_cyl*int_phi(:,inu), Rp_transit(inu))
      end do
      
      ! Calculate final transit radius by what is below the lower boundary
      Rp_transit = sqrt(rrf(ndepth+1)**2 + Rp_transit/pi)
      
      ! Add contribution from thermal emission from planet. Calculate emission
      ! spectrum from planet from UM output.
      if (trans_spec_emission) then
        allocate(fnu_surf_spec(nfreq_loc))
        call calc_surf_spec_3d(long_obs,lat_obs,nfreq_loc, &
            fnu_surf_spec)
        Rp_transit = Rp_transit/sqrt(1 + (rrf(1)**2*fnu_surf_spec)/ &
            (rstar**2*(4*pi*irad)))
        deallocate(fnu_surf_spec)
      endif

      
      deallocate(tau_path, int_phi, rho_cyl, phi_cyl)

    end subroutine calc_trans_radius_3d
    
    
    !> Calculate optical depth through atmosphere along path between star and
    !! observer.
    !<
    subroutine calc_tau_path(nfreq_loc, rho_cyl, phi_cyl, tau_path)
    
      use mod_grid
      use mod_cst
      use mod_param
      use mod_util
      
      implicit none
      
      ! Number of frequencies
      integer,intent(in) :: nfreq_loc
      ! Cylindrical coordinate of path
      real,intent(in) :: rho_cyl, phi_cyl
      ! Calculated path optical depth
      real,dimension(:),intent(out) :: tau_path
      
      ! Desired step in path
      real :: z_cyl_step
      ! Start point and end point for path integration
      real :: z_cyl_min, z_cyl_max
      ! Inclination angle in a coordinate system where z points towards the
      ! observer
      real :: theta_cyl_cur
      ! Current path point in cartesian coordinates where z points towards the
      ! observer
      real :: x_cyl_cur, y_cyl_cur
      real,dimension(:),allocatable :: z_cyl_cur
      ! Current path point in cartesian coordinates where z points towards north
      ! pole
      real :: x_sph_cur, y_sph_cur, z_sph_cur
      ! Current path point in spherical coordinates where z points towards north
      ! pole
      real :: rrf_sph_cur, theta_sph_cur, phi_sph_cur
      ! Arrays to store quantities along path
      real,dimension(:),allocatable :: rho_cur
      real,dimension(:,:),allocatable :: kapnu_cur
      ! True if path loop should exit
      logical :: exit_path_loop
      ! Number of integration steps, loop index
      integer :: nsteps, i
      
      ! Calculate step in integration path
      z_cyl_step = (longitude_um(2) - longitude_um(1))*rrf(ndepth+1)/10.
      
      ! Calculate minimum z_cyl for current rho, phi
      theta_cyl_cur = asin(rho_cyl/rrf(1))
      z_cyl_min = -rrf(1)*cos(theta_cyl_cur)
      z_cyl_max = rrf(1)*cos(theta_cyl_cur)
      
      ! Calculate number of required integration steps
      nsteps = ceiling((z_cyl_max - z_cyl_min)/z_cyl_step + 1)
      
      ! Allocate arrays to store integration parameters
      allocate(z_cyl_cur(nsteps),rho_cur(nsteps),kapnu_cur(nsteps,nfreq_loc))
      
      ! Calculate the corresponding value for x,y in the cylindrical
      ! coordinate system
      x_cyl_cur = rho_cyl*cos(phi_cyl)
      y_cyl_cur = rho_cyl*sin(phi_cyl)
      
      ! Loop along z_cyl and calculate optical depth along path
      i = 0
      exit_path_loop = .false.
      do while (.not.exit_path_loop)
        i = i+1
        
        ! Calculate current value for z_cyl
        z_cyl_cur(i) = z_cyl_min + (i-1)*z_cyl_step
        
        ! Check if the current z_cyl_cur is out of range
        if (z_cyl_cur(i) >= z_cyl_max) then
          z_cyl_cur(i) = z_cyl_max
          exit_path_loop = .true.
        end if
        
        ! Transform x,y,z in the cylindrical coordinate system where z faces
        ! the observer to a coordinate system where z is towards the north pole.
        ! In this coordinate system x is towards the observer.
        x_sph_cur = z_cyl_cur(i)
        y_sph_cur = y_cyl_cur
        z_sph_cur = -x_cyl_cur
        
        ! Convert to spherical coordinates
        rrf_sph_cur = sqrt(x_sph_cur**2 + y_sph_cur**2 + z_sph_cur**2)
        theta_sph_cur = acos(z_sph_cur/rrf_sph_cur)
        phi_sph_cur = atan2(y_sph_cur, x_sph_cur) + pi
        
        ! Observer is towards long_obs, lat_obs, rotate to get the correct path
        phi_sph_cur = phi_sph_cur + long_obs
        theta_sph_cur = theta_sph_cur + (lat_obs - pi/2.)
        
        ! Make sure that phi is between 0 and 2*pi and theta is between 0 and pi
        x_sph_cur = rrf_sph_cur*sin(theta_sph_cur)*cos(phi_sph_cur)
        y_sph_cur = rrf_sph_cur*sin(theta_sph_cur)*sin(phi_sph_cur)
        z_sph_cur = rrf_sph_cur*cos(theta_sph_cur)
        rrf_sph_cur = sqrt(x_sph_cur**2 + y_sph_cur**2 + z_sph_cur**2)
        theta_sph_cur = acos(z_sph_cur/rrf_sph_cur)
        phi_sph_cur = atan2(y_sph_cur, x_sph_cur) + pi
        
        ! Interpolate UM output to get the density and opacity at this point
        call interp_3d(longitude_um, latitude_um, rrf(ndepth+1:1:-1), &
            rhof_um(:,:,ndepth+1:1:-1), phi_sph_cur, theta_sph_cur, &
            rrf_sph_cur, rho_cur(i))
        call interp_3d_array(longitude_um, latitude_um, rrf(ndepth+1:1:-1), &
            kapfnu_um(:,:,ndepth+1:1:-1,:), phi_sph_cur, theta_sph_cur, &
            rrf_sph_cur, kapnu_cur(i,:))
      end do
      
      ! Integrate to get the optical depth
      do i=1,nfreq_loc
        call integrate_trapezoid(z_cyl_cur, rho_cur*kapnu_cur(:,i), tau_path(i))
      end do
      
      deallocate(z_cyl_cur,rho_cur,kapnu_cur)

    end subroutine calc_tau_path
    
    
    !> Calculate surface spectrum observed at Earth from UM output.
    !<
    subroutine calc_surf_spec_3d(long_obs_loc,lat_obs_loc,nfreq_loc, &
        fnu_surf_spec)
    
      use mod_grid
      use mod_cst
      use mod_param
      use mod_util
      
      implicit none
      
      ! Longitude and latitude of observer
      real,intent(in) :: long_obs_loc,lat_obs_loc
      ! Number of frequencies
      integer,intent(in) :: nfreq_loc
      ! Planet emission spectrum
      real,dimension(:),intent(out) :: fnu_surf_spec
      
      ! Latitude and longitude integration coordinates
      real,dimension(:),allocatable :: theta_loc,phi_loc
      ! Length of theta and phi
      integer :: ntheta_loc,nphi_loc
      ! Maximum and minimum longitude for integration
      real :: long_min,long_max
      ! Maximum and minimum latitude for integration
      real :: lat_min,lat_max
      ! Interpolated intensity in direction of observer
      real,dimension(:,:),allocatable :: Inu_surf_obs
      ! Integral over phi
      real,dimension(:,:),allocatable :: int_phi
      ! Cartesian coordinates in the integral coordinate system
      real :: x_int,y_int,z_int
      ! Cartesian coordinates in the coordinate system used by the UM
      real :: x_um,y_um,z_um
      ! Spherical coordinates in the coordinate system of the UM
      real :: theta_um,phi_um
      ! Temporary variable
      real :: tmp
      ! Loop indices
      integer :: itheta,iphi,inu,imu
            
      ! Construct longitude and latitude grid for integration
      ntheta_loc = nlat/2
      allocate(theta_loc(ntheta_loc))
      call linspace(0., pi/2., ntheta_loc, theta_loc)
      
      nphi_loc = nlong
      allocate(phi_loc(nphi_loc))
      call linspacep(0., 2*pi, nphi_loc, phi_loc)
      
      ! Loop over theta and phi
      allocate(Inu_surf_obs(nphi_loc, nfreq_loc))
      allocate(int_phi(ntheta_loc, nfreq_loc))

      
      do itheta=1,ntheta_loc
        do iphi=1,nphi_loc
          
          ! Convert current theta and phi to the coordinate system of the UM.
          ! First convert to Cartesian coordinates. Use radius = 1.
          x_int = sin(theta_loc(itheta))*cos(phi_loc(iphi))
          y_int = sin(theta_loc(itheta))*sin(phi_loc(iphi))
          z_int = cos(theta_loc(itheta))
          
          ! Transform x,y,z in the coordinate system where z faces the observer
          ! to a coordinate system where z is towards the north pole.
          ! In this coordinate system x is towards the observer.
          x_um = z_int
          y_um = y_int
          z_um = -x_int
        
          ! Convert to spherical coordinates
          theta_um = acos(z_um)
          phi_um = atan2(y_um, x_um) + pi
        
          ! Observer is towards long_obs_loc, lat_obs_loc, rotate to get the
          ! correct point
          phi_um = phi_um + long_obs_loc
          theta_um = theta_um + (lat_obs_loc - pi/2.)
        
          ! Make sure that phi is between 0 and 2*pi and theta is between 0 and pi
          x_um = sin(theta_um)*cos(phi_um)
          y_um = sin(theta_um)*sin(phi_um)
          z_um = cos(theta_um)
          theta_um = acos(z_um)
          phi_um = atan2(y_um, x_um) + pi
          
          ! Get outgoing intensity in the direction of the observer at current
          ! latitude and longitude
          ! Get outgoing intensity in the direction of the observer at current
          ! latitude and longitude
          call interp_3d_array(longitude_um, latitude_um, xmu_surf, &
          Inu_surf_um(:,:,:,:), phi_um, theta_um, cos(theta_loc(itheta)), &
          Inu_surf_obs(iphi,:))
        end do
        
        ! Integrate over phi
        do inu=1,nfreq_loc
          call integrate_trapezoid_loop(phi_loc, Inu_surf_obs(:,inu), &
              int_phi(itheta,inu))
        end do
      
      end do
      
      ! Integrate over theta
      do inu=1,nfreq_loc
        call integrate_trapezoid(theta_loc, int_phi(:,inu)*cos(theta_loc)*sin(theta_loc), &
            fnu_surf_spec(inu))
      end do
      
      deallocate(Inu_surf_obs)
      deallocate(int_phi)
  
    end subroutine calc_surf_spec_3d
    
    
    !> Copy 1D profiles to 3D arrays for storage.
    !<
    subroutine save_um_radtrans_vars(ilong, ilat)
    
      implicit none
    
      ! Longitude and latitude indices
      integer,intent(in) :: ilong, ilat
      
      integer :: i, j
    
      Inu_surf_um(ilong,ilat,:,:) = Inu_surf(:,:)
      
    end subroutine save_um_radtrans_vars
    
    
    !> Set the star zenith angle from the longitude and latitude. Assumes planet
    !! is tidally locked with substellar point at phi_um = pi, theta_um = pi/2.
    !<
    subroutine set_murad(ilong, ilat)
      
      use mod_grid
      implicit none
    
      ! Longitude and latitude indices
      integer,intent(in) :: ilong, ilat
      ! Cartesian coordinates in the coordinate system of the UM
      real :: x_um,y_um,z_um
      ! Cartesian coordinates in the coordinate system where z points towards
      ! the star
      real :: x_mu,y_mu,z_mu
      
      ! Convert current coordinate to cartesian coordinates. Use radius = 1.
      x_um = cos(latitude_um(ilat))*cos(longitude_um(ilong))
      y_um = cos(latitude_um(ilat))*sin(longitude_um(ilong))
      z_um = sin(latitude_um(ilat))
      
      ! Convert to a coordinate system where the z axis points towards the star
      x_mu = z_um
      y_mu = y_um
      z_mu = -x_um
      
      ! The star zenith angle is now simply the zenith angle in this coordinate
      ! system.
      murad = z_mu
      
      ! The star zenith angle is now simply the zenith angle. Check if angle is
      ! smaller than 0. If this is the case, irradiation is turned off for this
      ! column. Otherwise, turn it on with the correct zenith angle.
      if (z_mu < 0) then
        irrad = .false.
        murad = 0.
      else
        irrad = .true.
        murad = z_mu
      end if
    
    end subroutine set_murad
    
    
!!! Subroutine for Calculating Contribution Function
    
subroutine calc_con_fun(nfreq_cf,taub_loc,bnu_loc)


use mod_grid, only : ndepth,cf,ncf,ppf,ttf,nuf,corr_k,nkmix,nband
use mod_cst

implicit none

integer :: idepth, ifreq, nfreq_cf

real,dimension((nfreq_cf)) :: mcf
real,dimension(ndepth+1,(nfreq_cf)) :: bnu_loc,taub_loc
real,dimension(ndepth+1,nfreq_cf*nkmix) :: bnu1



 if (.not. corr_k) then
  bnu_loc = get_bnu(ndepth+1,nfreq_cf,ttf)
  else
  
  bnu1 = get_bnu(ndepth+1,nfreq_cf*nkmix,ttf)
  
  do idepth = 1,ndepth+1
   do ifreq = 1,(nfreq_cf*nkmix),nkmix  !!! nfreq_cf is nband in corr-k
    bnu_loc(idepth,((ifreq-1)/nkmix)+1) = bnu1(idepth,ifreq)
   end do
  end do
 end if

 
  do idepth = 2,ndepth+1
    do ifreq = 1,(nfreq_cf)
  
     cf(idepth-1,ifreq) = ((exp(-taub_loc(idepth,ifreq))-exp(-taub_loc(idepth-1,ifreq)))/(log(ppf(idepth-1))-log(ppf(idepth))))*(bnu_loc(idepth,ifreq))
    end do 
  end do


  do ifreq = 1,(nfreq_cf)
   mcf(ifreq) = maxval(cf(:,ifreq), 1)
  end do

 
  do idepth = 1,ndepth
   do ifreq = 1,(nfreq_cf)
     ncf(idepth,ifreq) = cf(idepth,ifreq)/mcf(ifreq)  
   end do
  end do


 end subroutine calc_con_fun
 
    
end module mod_radtrans

