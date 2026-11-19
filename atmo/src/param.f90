!> \file
!! Contains module mod_param, subroutine mod_param::init_param()\n
!! Initialization of the parameters, input parameter file, debug level, mpi cst, output file names
!<

!> Hold the parameters of the simulations, input file parameters, debug levels, number of cpus, output files names
module mod_param

  integer        :: ncpu     !< number of cpu used
  integer        :: mype     !< number of the current processor
  integer        :: cputerm  !< indice of the cpu used for terminal print
  integer        :: debug    !< level of debug for dev
  character(256) :: fparam   !< name of the input file for namelists
  character(256) :: fin      !< name of the input profile
  character(128) :: f2din    !< name of the input u/alpha profile in 2d
  character(256) :: fsp_sw_um !< name of the shortwave spectral file
  character(256) :: fsp_lw_um !< name of the shortwave spectral file
  character(256) :: fout     !< name of the output profile
  character(256) :: ftemp    !< name of the temporary output profile
  character(256) :: line    !< Variable to store each line in the namelist file for checking it is valid.
  character(128) :: um_p_name !< name of the pressure input variable from um
  character(128) :: um_theta_name !< name of the potential temperature input variable from um
  logical        :: uncertainty
  logical        :: neq_feedback
  integer        :: nuncert,iuncert

  character(len=100) :: headrfmt,separfmt,strpmfmt,chrpmfmt,intpmfmt,fltpmfmt,exppmfmt,logpmfmt
  data separfmt/'("****************************************************************************************")'/
  data headrfmt/'("*******************************",A24," ******************************")'/
  data strpmfmt/'(A18," = ",A50      )'/
  data chrpmfmt/'(A18," = ",A9  ,3x,A)'/
  data intpmfmt/'(A18," = ",I9  ,3x,A)'/
  data exppmfmt/'(A18," = ",ES9.2,3x,A)'/
  data fltpmfmt/'(A18," = ",F9.3,3x,A)'/
  data logpmfmt/'(A18," = ",L9  ,3x,A)'/

contains

  subroutine init_param

    implicit none

    integer :: stat
    logical :: exists

    namelist /param/debug,fin,fout,cputerm,uncertainty,nuncert,um_p_name,um_theta_name,f2din, &
    fsp_sw_um, fsp_lw_um,ftemp

    ! Check if a name for the input parameter file has been provided and
    ! if it exists, otherwise look for the standard file name 'param.in'
    fparam=' '
    call getarg(1,fparam)
     if(fparam.ne.' ')then
       inquire(file=fparam,exist=exists)
       if(.not.exists)then
          write(*,'(A,1x,A,1x,A)')'The file',trim(fparam),'does not exist!'
          stop
       endif
    else
       fparam='param.in'
       inquire(file=fparam,exist=exists)
       if(.not.exists)then
          write(*,'(A)')'No input parameter file found!'
          stop
       endif
    endif

    !default param

    f2din   ='None'
    fin     ='None'
    fout    ='fout.ncdf'
    ftemp   = 'temp.ncdf'
    um_p_name = 'p'
    um_theta_name = 'theta'
    fsp_sw_um = 'None'
    fsp_lw_um = 'None'
    debug   = 0
    cputerm = 0
    uncertainty = .false.
    nuncert = 1
    
    open(1,file=fparam)
    rewind(1) ; read(1,NML=param,iostat=stat)
     IF (stat/=0) THEN
      BACKSPACE(1)
      READ(1,fmt='(A)') line
      WRITE(*,'(A)') 'Invalid line in namelist: '//trim(line)
      STOP
     END IF    
    close(1)

if (mype==cputerm) then
    write(*,headrfmt)'Simulation Parameters'
    write(*,intpmfmt) 'ncpu',ncpu,'Number of cpus'
    write(*,intpmfmt) 'cputerm',cputerm,'Id of the cpu used for terminal print'
    write(*,intpmfmt) 'debug',debug,'Level of debug'
    write(*,strpmfmt) 'fparam','Parameter file: '//fparam
    write(*,strpmfmt) 'fin',   'Input file    : '//fin
    write(*,strpmfmt) 'fout',  'Output file   : '//fout
    write(*,strpmfmt) 'ftemp', 'Temp file     : '//ftemp
    write(*,*) 
    write(*,separfmt)
endif

  end subroutine init_param

end module mod_param
