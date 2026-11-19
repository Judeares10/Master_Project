!> \file
!! Contains main program Atmo
!<

!> \mainpage
!! <CENTER> Atmo: 1D stellar and planetary atmosphere code\n
!! Version: V2.0, 11/2013\n
!! Contact: P. Tremblin, tremblin@astro.ex.ac.uk\n
!! </CENTER>
!!
!! \n
!! The documentation and call graph of the main program is here: \link atmo </A>
!! \n
!! - General architecture: P. Tremblin, I. Baraffe\n
!! - Solver for the hydrostatic and energy equations to compute the pressure/temperature profiles: W. Hayek, P. Tremblin\n
!!  \image html pt_t700.png
!! - Line by line radiative transfer including scattering and irradiation: D. Amundsen, W. Hayek, P. Tremblin\n
!! \image html gj436b_spectrums.png 
!! - Opacities for H2O, H2-H2, H2-He, CH4, NH3, CO, N2, TiO, VO: D. Amundsen, I. Baraffe\n
!! - Analytical abundances (ref Burrows & Sharp 1999): D. Amundsen, I. Baraffe\n
!! - Equilibrium chemistry by minimization of Gibbs free energy: P. Tremblin, G. Chabrier\n
!! \image html abundances_gibbs.png
!! - Out of equilibrium chemistry (network venot et al 2012): B. Drummond, P. Tremblin \n
!! - Photo-chemistry (in progress): B. Drummond
!! - Condensation: P. Mourier, P. Tremblin, I. Baraffe, G. Chabrier
!! - GCM chemistry/radiative transfer coupling: D. Amundsen, B. Drummond, P. Tremblin, N. Mayne\n 
!<
program atmo

  use mod_param
  use mod_cst
  use mod_grid
  use mod_chem
  use mod_chem_neq
  use mod_eos
  use mod_opacity
  use mod_radtrans
  use mod_conv
  use mod_solver
  use mod_chem_relax

  implicit none

#if PARA==1 

  include 'mpif.h'
  integer :: info

  call MPI_INIT(info)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mype,info)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,ncpu,info)
  cputerm = 0

#else

  ncpu    = 1
  mype    = 0
  cputerm = 0

#endif

 if (mype == cputerm) then
     write(*,*)
     write(*,separfmt)
     write(*,*) 
     write(*,'(a)') '                             _/    _/_/_/  _/  _/  _/_/_/                    '
     write(*,'(a)') '                           _/  _/    _/    _/_/_/  _/  _/                    '
     write(*,'(a)') '                           _/_/_/    _/    _/  _/  _/  _/                    '
     write(*,'(a)') '                           _/  _/    _/    _/  _/  _/  _/                    '
     write(*,'(a)') '                           _/  _/    _/    _/  _/  _/_/_/                    '
     write(*,'(a)') '                                                                             '
     write(*,'(a)') '                   Contact: P. Tremblin, tremblin@astro.ex.ac.uk             '
     write(*,'(a)') '                              Version: 2.0   (11/2013)                       '
     write(*,*) 
     write(*,separfmt) 
     write(*,*)
  endif

  ! Initialize modules
  if (debug>0 .and. mype==cputerm) write(*,*) 'start init_param'
  call init_param
  if (debug>0 .and. mype==cputerm) write(*,*) 'end init_param'

  if (debug>0 .and. mype==cputerm) write(*,*) 'start init_cst'
  call init_cst
  if (debug>0 .and. mype==cputerm) write(*,*) 'end init_cst'

  if (debug>0 .and. mype==cputerm) write(*,*) 'start init_grid'
  call init_grid
  if (debug>0 .and. mype==cputerm) write(*,*) 'end init_grid'

  if (debug>0 .and. mype==cputerm) write(*,*) 'start init_chem'
  call init_chem
  if (debug>0 .and. mype==cputerm) write(*,*) 'end init_chem'

  if (debug>0 .and. mype==cputerm) write(*,*) 'start init_eos'
  call init_eos
  if (debug>0 .and. mype==cputerm) write(*,*) 'end init_eos'
  
  if (debug>0 .and. mype==cputerm) write(*,*) 'start init_neq_chem'
  if(chem=='neq') then
     call init_neq_chem
  endif
  if (debug>0 .and. mype==cputerm) write(*,*) 'end init_neq_chem'

  if (debug>0 .and. mype==cputerm) write(*,*) 'start init_opacity'
  call init_opacity
  if (debug>0 .and. mype==cputerm) write(*,*) 'end init_opacity'

  if (debug>0 .and. mype==cputerm) write(*,*) 'start init_radtrans'
  call init_radtrans(nukapmin,nukapmax,nfreqkap)
  if (debug>0 .and. mype==cputerm) write(*,*) 'end init_radtrans'
  
  if (debug>0 .and. mype==cputerm) write(*,*) 'start init_conv'
  call init_conv
  if (debug>0 .and. mype==cputerm) write(*,*) 'end init_conv'

  if (debug>0 .and. mype==cputerm) write(*,*) 'start init_solver'
  call init_solver
  if (debug>0 .and. mype==cputerm) write(*,*) 'end init_solver'

  if(chem == 'neq') then
     if (debug>0 .and. mype==cputerm) write(*,*) 'start get_neq_chemistry'
     if (uncertainty) then
        if (debug>0 .and. mype==cputerm) write(*,*) 'start Monte Carlo Runs'
        do iuncert = 1, nuncert
           call get_neq_chemistry
           call init_chem
           call init_neq_chem
        enddo
        write(*,*) 'End of Monte Carlo Runs'
        stop
     endif
     ! Initialise abundances with equilibrium abundances 
     if (init_eq) then
       call get_eq_chemistry
     end if
     call get_neq_chemistry
     if (debug>0 .and. mype==cputerm) write(*,*) 'end get_neq_chemistry'
  else if (chem == 'relax') then
    if (debug>0 .and. mype==cputerm) write(*,*) 'start get_relax_chemistry'
      call init_relax_chemistry
      call get_relax_chemistry
    if (debug>0 .and. mype==cputerm) write(*,*) 'end get_relax_chemistry'
  end if
  
  if (ndim==1) then
     if (solve_hydro .or. solve_energy) then
        if (um_input .and. mype==cputerm) then
           write(*,headrfmt) 'Error atmo'
           write(*,'(A)') 'solve_hydro and solve_energy cannot be used with ' // &
          'um_input. Please turn solve_hydro and solve_energy off.'
           write(*,*)
           write(*,separfmt)
           stop
        end if
        if (debug>0 .and. mype==cputerm) write(*,*) 'start find_profile'
        if (any(chem==(/'eq ','ana','cst','man'/))) then
           call find_profile(chem)
        endif
        if (debug>0 .and. mype==cputerm) write(*,*) 'end find_profile'
     endif
  else if (ndim==2) then
     if (solve_hydro .and. solve_energy) then
        call find_profile_2d
     endif
  endif
     
     if(surface_spectrum .or. transmission_spectrum .or. phase_curve) then
 
    if (um_input) then
    
      if (debug>0 .and. mype==cputerm) write(*,*) 'start compute_spectrum_um'
      if (any(chem==(/'eq ','ana','man'/))) then
          call compute_spectrum_um(chem)
      else if(any(chem == (/'neq','cst'/))) then      
          call compute_spectrum_um('cst')
      endif
      if (debug>0 .and. mype==cputerm) write(*,*) 'end compute_spectrum_um'
      
    else
    
      if (debug>0 .and. mype==cputerm) write(*,*) 'start compute_spectrum'
      if (any(chem==(/'eq ','ana','man'/))) then
          call compute_spectrum(chem,0)
       else if(any(chem == (/'neq','cst'/))) then      
          call compute_spectrum('cst',0)
       endif
      if (debug>0 .and. mype==cputerm) write(*,*) 'end compute_spectrum'
      
    endif
 endif
 
 if (transmission_spectrum) then
   if (.not. um_input .and. transmission_spectrum_3d) then
     write(*,headrfmt) 'Error atmo'
     write(*,'(A)') 'Cannot calculate 3D transmission spectrum without ' // &
        'UM input. Please select 1D transmission spectrum or provide ' // &
        'UM input.'
     write(*,*)
     write(*,separfmt)
     stop
   else
     if (debug>0 .and. mype==cputerm) write(*,*) 'start compute_trans_spec'
     call compute_trans_spec
     if (debug>0 .and. mype==cputerm) write(*,*) 'end compute_trans_spec'
   endif
 endif
 
 if (phase_curve) then
   if (debug>0 .and. mype==cputerm) write(*,*) 'start compute_phase_curve'
   call compute_phase_curve
   if (debug>0 .and. mype==cputerm) write(*,*) 'end compute_phase_curve'
 end if

 if (.not. (solve_hydro .or. solve_energy .or. surface_spectrum .or. transmission_spectrum .or. phase_curve) .and. chem=='eq') then
    
    if (debug>0 .and. mype==cputerm) write(*,*) 'start get_eq_chemistry'
    call get_eq_chemistry
    if (debug>0 .and. mype==cputerm) write(*,*) 'end get_eq_chemistry'
 
 endif


# if PARA == 1

  call MPI_FINALIZE(info)

#endif

  stop

end program atmo



