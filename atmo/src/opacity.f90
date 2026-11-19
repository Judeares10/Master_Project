!> \file
!! Contains module mod_opacity, subroutine mod_opacity::init_opacity() \n
!! Initialization of the opacity table or constant
!<


!> hold the parameters of the opacity table
!<
module mod_opacity

  integer :: nkap        !< number of molecules used for opacities
  integer,parameter :: nkapmax = 22 !< maximum number of molecules for arrays
  integer,parameter :: lname   = 9  !< Maximum lenght of the molecule names
  character(lname),dimension(nkapmax) :: kapname !< names of the molecules used for opacities
  real,dimension(nkapmax) :: kapcst     !< Constant opacities, if 0. use a table
  
  character(128),dimension(nkapmax) :: fktab  !< name of the files containing the opacity table
  real,   dimension(nkapmax) :: dnuktab   !< frequency resolution of the tables
  real,   dimension(nkapmax) :: nuktabmin !< minimum frequency of the tables
  real,   dimension(nkapmax) :: nuktabmax !< maximum frequency of the tables
  integer,dimension(nkapmax) :: nfreqktab !< number of frequencies in the tables

  integer :: npktab      !< number of pressure point in the opacity table
  real    :: pktabmin    !< Minimum pressure in the opacity table
  real    :: pktabmax    !< Maximum pressure in the opacity table

  integer :: ntktab      !< number of temperature point in the opacity table
  real    :: tktabmin    !< Minimum temperature in the opacity table
  real    :: tktabmax    !< Maximum temperature in the opacity table
  
  real    :: unitPktab   !< Table pressure unit to cgs
  real    :: unitKktab   !< Table temperature unit to cgs
  real    :: unitnuktab  !< Table frequency unit to cgs
  real    :: unitLktab   !< Table length unit to cgs
  real    :: unitMktab   !< Table mass unit to cgs

  real    :: unitLkcst   !< kcst length unit to cgs
  real    :: unitMkcst   !< kcst mass unit to cgs

  integer :: ntab        !< number of opacity table used

  real   ,dimension(:,:,:,:),allocatable :: ktab,kstdtab  !< Opacity table

  real   ,dimension(:      ),allocatable :: pktab !< Pressure    for the opacity table
  real   ,dimension(:      ),allocatable :: tktab !< Temperature for the opacity table
  
  integer :: nfreqkap_glob  !< global number of frequencies with an opacity
  real    :: nukapmin_glob  !< global minimum frequency with an opacity 
  real    :: nukapmax_glob  !< global maximum frequency with an opacity

  integer :: nfreqkap       !< local-cpu number of frequency with an opacity
  real    :: nukapmin       !< local-cpu minimum frequency with an opacity 
  real    :: nukapmax       !< local-cpu maximum frequency with an opacity
  
  integer :: inukap_start   !< local-cpu start indice in nuf_glob
  integer :: inukap_end     !< local-cpu end indice in nuf_glob


  integer :: itauref !< reference layer for weights of the mixture
 
  logical :: kap_smooth    !< if T opacities are vertically smoothed
  real    :: ppmin_smooth  !< minimum pressure above which kap is smoothed
  real    :: ppmax_smooth  !< maximum pressure above which kap is smoothed
  integer :: kerkap_smooth !< number of neighbors used for smoothing (1 means smoothing on 3 cells)
  real    :: unitPsmooth
  
  !! Absorption Cloud deck 
  logical :: cloud          !< logical to include cloud deck or not
  real :: cloud_strength	  !< Cloud absorption
  real    :: ppmin_cloud  !< minimum pressure above which cloud opacity is added
  real    :: ppmax_cloud  !< maximum pressure above which cloud opacity is added
  real    :: unitPcloud

  ! Writing corr k opacity
  logical :: write_opacity !< Flag for writing mixed corr k opacity to file
  character(128) :: fkabsout !< filename for mixed corr k opacity
  integer :: ikapout !< integer for opacity to be written out. 0 for total opacity


contains

  !> initialize the opacity table
  !<
  subroutine init_opacity

    use mod_param
    use mod_cst
    use mod_chem
    use mod_grid
    use mod_util


    implicit none
 

    !local variable
    integer :: stat

    namelist /opacity/ nkap,kapname,fktab,nfreqktab,dnuktab,nuktabmin,nuktabmax,kapcst,npktab,ntktab,pktabmin,pktabmax,tktabmin,tktabmax,unitPktab,unitKktab,unitnuktab,unitLktab,unitMktab,unitLkcst,unitMkcst,kap_smooth,ppmin_smooth,ppmax_smooth,kerkap_smooth,unitPsmooth,cloud_strength,ppmax_cloud,ppmin_cloud,cloud,unitPcloud,write_opacity,fkabsout


    nkap    = 0
       
    cloud = .false.
    unitPcloud  = 1.0E6
    ppmin_cloud = 1E-10
    ppmax_cloud = 1E6
    cloud_strength = 1

    ! pressure, temperature table 800 points per default in cgs
    npktab   = 40
    ntktab   = 20
    pktabmin = 7.8804628E-4
    pktabmax = 1.0E9
    tktabmin = 70.
    tktabmax = 3000.



    !units pressure/temperature in cgs frequencies in m-1, convert to cm-1
    unitPktab  = 1.
    unitKktab  = 1.
    unitnuktab = 0.01
    unitLktab  = 100.   !m  -> cm
    unitMktab  = 1000.  !kg -> g

    unitLkcst  = 100.   !m  -> cm
    unitMkcst  = 1000.  !kg -> g

    !smoothing parameters
    kap_smooth    = .false.
    ppmin_smooth  = 1E-9
    ppmax_smooth  = 1000.
    unitPsmooth   = 1.0E6  !bars
    kerkap_smooth = 1

    ! Writing corr k opacity parameters
    write_opacity = .False.
    fkabsout = 'kabs_out.ncdf'
    
    if (.not. corr_k) then

! Reduced LbL stuff
!       kapname(1) = 'H2'  ; fktab(1) = '../../kabs/rabs_coeff_H2-H2_HITRAN_pt800.nc'
!       kapname(2) = 'He'  ; fktab(2) = '../../kabs/rabs_coeff_H2-He_HITRAN_pt800.nc'
!       kapname(3) = 'H2O' ; fktab(3) = '../../kabs/rabs_coeff_H2O_BT2_pt800.nc'
!       kapname(4) = 'CO'  ; fktab(4) = '../../kabs/rabs_coeff_CO_HITEMP_pt800.nc'
!       kapname(5) = 'CO2' ; fktab(5) = '../../kabs/rabs_coeff_CO2_CDSD-4000_pt800.nc'
!       kapname(6) = 'CH4' ; fktab(6) = '../../kabs/rabs_coeff_CH4_YT10to10_pt800.nc'
!       kapname(6) = 'CH4' ; fktab(6) = '../../kabs/rabs_coeff_CH4_STDS_pt800.nc'
!       kapname(7) = 'NH3' ; fktab(7) = '../../kabs/rabs_coeff_NH3_BYTe_pt800.nc'
!       kapname(8) = 'Na'  ; fktab(8) = '../../kabs/rabs_coeff_Na_VALD_pt800.nc'
!       kapname(9) = 'K'   ; fktab(9) = '../../kabs/rabs_coeff_K_VALD_pt800.nc'
!       kapname(10)= 'C2H2'; fktab(10)= '../../kabs/rabs_coeff_C2H2_HITRAN_pt800.nc'
!       kapname(11)= 'TiO' ; fktab(11)= '../../kabs/rabs_coeff_TiO_Plez_pt800.nc'
!       kapname(12)= 'VO'  ; fktab(12)= '../../kabs/rabs_coeff_VO_Plez_pt800.nc'
       

       
       !frequency range of the tables units: m-1
       nfreqktab(:) = 50000 ; dnuktab(:) = 100. 
       nuktabmin(:) = 0.; nuktabmax(:) = 5.0E6
      
! Standard LbL stuff
     kapname(1)  = 'H2'   ; fktab(1)  = '/scratch/dp015/dc-trem1/kabs/abs_coeff_H2-H2_HITRAN_pt800.nc'
     kapname(2)  = 'He'   ; fktab(2)  = '/scratch/dp015/dc-trem1/kabs/abs_coeff_H2-He_HITRAN_pt800.nc'
     kapname(3)  = 'H2O'  ; fktab(3)  = '/scratch/dp015/dc-trem1/kabs/abs_coeff_H2O_BT2_pt800.nc'
     kapname(4)  = 'CO2'  ; fktab(4)  = '/scratch/dp015/dc-trem1/kabs/abs_coeff_CO2_CDSD-4000_pt800.nc'
     kapname(5)  = 'CO'   ; fktab(5)  = '/scratch/dp015/dc-trem1/kabs/abs_coeff_CO_HITEMP_pt800.nc'
     kapname(6)  = 'CH4'  ; fktab(6)  = '/scratch/dp015/dc-trem1/kabs/abs_coeff_CH4_YT10to10_pt800.nc'
     kapname(7)  = 'NH3'  ; fktab(7)  = '/scratch/dp015/dc-trem1/kabs/abs_coeff_NH3_BYTe_pt800.nc'
     kapname(8)  = 'Na'   ; fktab(8)  = '/scratch/dp015/dc-trem1/kabs/abs_coeff_Na_VALD_pt800.nc'
     kapname(9)  = 'K'    ; fktab(9)  = '/scratch/dp015/dc-trem1/kabs/abs_coeff_K_VALD_pt800_A16_H2He.nc'
     kapname(10) = 'C2H2' ; fktab(10) = '/scratch/dp015/dc-trem1/kabs/abs_coeff_C2H2_HITRAN_pt800.nc'
     kapname(11) = 'TiO'  ; fktab(11) = '/scratch/dp015/dc-trem1/kabs/abs_coeff_TiO_Plez_pt800.nc'
     kapname(12) = 'VO'   ; fktab(12) = '/scratch/dp015/dc-trem1/kabs/abs_coeff_VO_VOMYT_pt800.nc'
     kapname(13) = 'Li'   ; fktab(13) = '/scratch/dp015/dc-trem1/kabs/abs_coeff_Li_VALD_pt800.nc'
     kapname(14) = 'Rb'   ; fktab(14) = '/scratch/dp015/dc-trem1/kabs/abs_coeff_Rb_VALD_pt800.nc'
     kapname(15) = 'Cs'   ; fktab(15) = '/scratch/dp015/dc-trem1/kabs/abs_coeff_Cs_VALD_pt800.nc'
     kapname(16) = 'FeH'  ; fktab(16) = '/scratch/dp015/dc-trem1/kabs/abs_coeff_FeH_Yueqi_pt800.nc'
     kapname(17) = 'PH3'  ; fktab(17) = '/scratch/dp015/dc-trem1/kabs/abs_coeff_PH3_SAlTY_pt800.nc'
     kapname(18) = 'H2S'  ; fktab(18) = '/scratch/dp015/dc-trem1/kabs/abs_coeff_H2S_HITRAN_pt800.nc'
     kapname(19) = 'HCN'  ; fktab(19) = '/scratch/dp015/dc-trem1/kabs/abs_coeff_HCN_Harris_pt800.nc'
     kapname(20) = 'SO2'  ; fktab(20) = '/scratch/dp015/dc-trem1/kabs/abs_coeff_SO2_ExoAmes_pt800.nc'


 !kapname(13) = 'Fe'   ; fktab(13) = '/scratch/dp015/dc-trem1/kabs/abs_coeff_Fe_VALD_pt800.nc'
! 
!     
!     !frequency range of the tables units: m-1
     nfreqktab(1)  = 10000     ; dnuktab(1)  = 100. ; nuktabmin(1)  = 0. ; nuktabmax(1)  = 1.0E6
     nfreqktab(2)  = 20000     ; dnuktab(2)  = 100. ; nuktabmin(2)  = 0. ; nuktabmax(2)  = 2.0E6
     nfreqktab(3)  = 30000000  ; dnuktab(3)  = 0.1  ; nuktabmin(3)  = 0. ; nuktabmax(3)  = 3.0E6
     nfreqktab(4)  = 12000000  ; dnuktab(4)  = 0.1  ; nuktabmin(4)  = 0. ; nuktabmax(4)  = 1.2E6
     nfreqktab(5)  = 10000000  ; dnuktab(5)  = 0.1  ; nuktabmin(5)  = 0. ; nuktabmax(5)  = 1.0E6
     nfreqktab(6)  = 15000000  ; dnuktab(6)  = 0.1  ; nuktabmin(6)  = 0. ; nuktabmax(6)  = 1.5E6
     nfreqktab(7)  = 12000000  ; dnuktab(7)  = 0.1  ; nuktabmin(7)  = 0. ; nuktabmax(7)  = 1.2E6
     nfreqktab(8)  = 50000000  ; dnuktab(8)  = 0.1  ; nuktabmin(8)  = 0. ; nuktabmax(8)  = 5.0E6
     nfreqktab(9)  = 500000    ; dnuktab(9)  = 10.0 ; nuktabmin(9)  = 0. ; nuktabmax(9)  = 5.0E6
     nfreqktab(10) = 10000000  ; dnuktab(10) = 0.1  ; nuktabmin(10) = 0. ; nuktabmax(10) = 1.0E6
     nfreqktab(11) = 35000000  ; dnuktab(11) = 0.1  ; nuktabmin(11) = 0. ; nuktabmax(11) = 3.5E6
     nfreqktab(12) = 40000000  ; dnuktab(12) = 0.1  ; nuktabmin(12) = 0. ; nuktabmax(12) = 4.0E6
     nfreqktab(13) = 50000000  ; dnuktab(13) = 0.1  ; nuktabmin(13) = 0. ; nuktabmax(13) = 5.0E6
     nfreqktab(14) = 50000000  ; dnuktab(14) = 0.1  ; nuktabmin(14) = 0. ; nuktabmax(14) = 5.0E6
     nfreqktab(15) = 50000000  ; dnuktab(15) = 0.1  ; nuktabmin(15) = 0. ; nuktabmax(15) = 5.0E6
     nfreqktab(16) = 20000000  ; dnuktab(16) = 0.1  ; nuktabmin(16) = 0. ; nuktabmax(16) = 2.0E6
     nfreqktab(17) = 12000000  ; dnuktab(17) = 0.1  ; nuktabmin(17) = 0. ; nuktabmax(17) = 1.2E6
     nfreqktab(18) = 15000000  ; dnuktab(18) = 0.1  ; nuktabmin(18) = 0. ; nuktabmax(18) = 1.5E6
     nfreqktab(19) = 20000000  ; dnuktab(19) = 0.1  ; nuktabmin(19) = 0. ; nuktabmax(19) = 2.0E6
     nfreqktab(20) = 10000000  ; dnuktab(20) = 0.1  ; nuktabmin(20) = 0. ; nuktabmax(20) = 1.0E6



      !nfreqktab(13) = 40000000  ; dnuktab(13) = 0.1  ; nuktabmin(13) = 0. ; nuktabmax(13) = 5.0E6
    
    else

       if (nband .le. 32) then
          
          kapname(1)  = 'H2'  ; fktab(1)  = '../../kabs/h2-h2_32_t5e-3_uw1116.nc'
          kapname(2)  = 'He'  ; fktab(2)  = '../../kabs/h2-he_32_t5e-3_uw1116.nc'
          kapname(3)  = 'H2O' ; fktab(3)  = '../../kabs/h2o_32_t5e-3_uw1116.nc'
          kapname(4)  = 'CO2' ; fktab(4)  = '../../kabs/co2_32_t5e-3_uw1116.nc'
          kapname(5)  = 'CO'  ; fktab(5)  = '../../kabs/co_32_t5e-3_uw1116.nc'
          kapname(6)  = 'CH4' ; fktab(6)  = '../../kabs/ch4_32_t5e-3_uw1116.nc'
          !kapname(6)  = 'CH4' ; fktab(6)  = '../../kabs/ch4_32_t5e-3_uw1116_STDS.nc'
          kapname(7)  = 'NH3' ; fktab(7)  = '../../kabs/nh3_32_t5e-3_uw1116.nc'
          kapname(8)  = 'Na'  ; fktab(8)  = '../../kabs/na_32_t5e-3_uw1116.nc'
          kapname(9)  = 'K'   ; fktab(9)  = '../../kabs/k_32_t5e-3_uw1116.nc'
          kapname(10) = 'Li'  ; fktab(10) = '../../kabs/li_32_t5e-3_uw1116.nc'
          kapname(11) = 'Rb'  ; fktab(11) = '../../kabs/rb_32_t5e-3_uw1116.nc'
          kapname(12) = 'Cs'  ; fktab(12) = '../../kabs/cs_32_t5e-3_uw1116.nc'
          kapname(13) = 'TiO' ; fktab(13) = '../../kabs/tio_32_t5e-3_uw1116.nc'
          kapname(14) = 'VO'  ; fktab(14) = '../../kabs/vo_32_t5e-3_uw1116_vomyt.nc'
          !kapname(14) = 'VO'  ; fktab(14) = '../../kabs/vo_32_t5e-3_uw1116.nc'
          kapname(15) = 'FeH' ; fktab(15) = '../../kabs/feh_32_t5e-3_uw1116.nc'
          kapname(16) = 'PH3' ; fktab(16) = '../../kabs/ph3_32_t5e-3_uw1116.nc'
          kapname(17) = 'H2S' ; fktab(17) = '../../kabs/h2s_32_t5e-3_uw1116.nc'
          kapname(18) = 'HCN' ; fktab(18) = '../../kabs/hcn_32_t5e-3_uw1116.nc'
          kapname(19) = 'C2H2'; fktab(19) = '../../kabs/c2h2_32_t5e-3_uw1116.nc'
          kapname(20) = 'SO2' ; fktab(20) = '../../kabs/so2_32_t5e-3_uw1116.nc'
          kapname(21) = 'Fe'  ; fktab(21) = '../../kabs/fe_32_t5e-3_uw1116.nc'
          !kapname(21) = 'CrH'  ; fktab(21) = '../../kabs/crh_32_t5e-3_uw1116.nc'
          

       else if (nband .le. 500) then

          kapname(1)  = 'H2'  ; fktab(1)  = '../../kabs/h2-h2_500_t5e-3_uw1116.nc'
          kapname(2)  = 'He'  ; fktab(2)  = '../../kabs/h2-he_500_t5e-3_uw1116.nc'
          kapname(3)  = 'H2O' ; fktab(3)  = '../../kabs/h2o_500_t5e-3_uw1116.nc'
          kapname(4)  = 'CO2' ; fktab(4)  = '../../kabs/co2_500_t5e-3_uw1116.nc'
          kapname(5)  = 'CO'  ; fktab(5)  = '../../kabs/co_500_t5e-3_uw1116.nc'
          kapname(6)  = 'CH4' ; fktab(6)  = '../../kabs/ch4_500_t5e-3_uw1116.nc'
          !kapname(6)  = 'CH4' ; fktab(6)  = '../../kabs/ch4_500_t5e-3_uw1116_STDS.nc'
          kapname(7)  = 'NH3' ; fktab(7)  = '../../kabs/nh3_500_t5e-3_uw1116.nc'
          kapname(8)  = 'Na'  ; fktab(8)  = '../../kabs/na_500_t5e-3_uw1116.nc'
          kapname(9)  = 'K'   ; fktab(9)  = '../../kabs/k_500_t5e-3_uw1116.nc'
          kapname(10) = 'Li'  ; fktab(10) = '../../kabs/li_500_t5e-3_uw1116.nc'
          kapname(11) = 'Rb'  ; fktab(11) = '../../kabs/rb_500_t5e-3_uw1116.nc'
          kapname(12) = 'Cs'  ; fktab(12) = '../../kabs/cs_500_t5e-3_uw1116.nc'
          kapname(13) = 'TiO' ; fktab(13) = '../../kabs/tio_500_t5e-3_uw1116.nc'
          kapname(14) = 'VO'  ; fktab(14) = '../../kabs/vo_500_t5e-3_uw1116_vomyt.nc'
          !kapname(14) = 'VO'  ; fktab(14) = '../../kabs/vo_500_t5e-3_uw1116.nc'
          kapname(15) = 'FeH' ; fktab(15) = '../../kabs/feh_500_t5e-3_uw1116.nc'
          kapname(16) = 'PH3' ; fktab(16) = '../../kabs/ph3_500_t5e-3_uw1116.nc'
          kapname(17) = 'H2S' ; fktab(17) = '../../kabs/h2s_500_t5e-3_uw1116.nc'
          kapname(18) = 'HCN' ; fktab(18) = '../../kabs/hcn_500_t5e-3_uw1116.nc'
          kapname(19) = 'C2H2'; fktab(19) = '../../kabs/c2h2_500_t5e-3_uw1116.nc'
          kapname(20) = 'SO2'; fktab(20) = '../../kabs/so2_500_t5e-3_uw1116.nc'
          kapname(21) = 'Fe';  fktab(21) = '../../kabs/fe_500_t5e-3_uw1116.nc'
          !kapname(21) = 'CrH'  ; fktab(21) = '../../kabs/crh_500_t5e-3_uw1116.nc'
         

       else if (nband .le. 5000) then

          kapname(1)  = 'H2'  ; fktab(1)  = '../../kabs/h2-h2_5000_t5e-3_uw1116.nc'
          kapname(2)  = 'He'  ; fktab(2)  = '../../kabs/h2-he_5000_t5e-3_uw1116.nc'
          kapname(3)  = 'H2O' ; fktab(3)  = '../../kabs/h2o_5000_t5e-3_uw1116.nc'
          kapname(4)  = 'CO2' ; fktab(4)  = '../../kabs/co2_5000_t5e-3_uw1116.nc'
          kapname(5)  = 'CO'  ; fktab(5)  = '../../kabs/co_5000_t5e-3_uw1116.nc'
          kapname(6)  = 'CH4' ; fktab(6)  = '../../kabs/ch4_5000_t5e-3_uw1116.nc'
          !kapname(6)  = 'CH4' ; fktab(6)  = '../../kabs/ch4_5000_t5e-3_uw1116_STDS.nc'
          kapname(7)  = 'NH3' ; fktab(7)  = '../../kabs/nh3_5000_t5e-3_uw1116.nc'
          kapname(8)  = 'Na'  ; fktab(8)  = '../../kabs/na_5000_t5e-3_uw1116.nc'
          kapname(9)  = 'K'   ; fktab(9)  = '../../kabs/k_5000_t5e-3_uw1116.nc'
          kapname(10) = 'Li'  ; fktab(10) = '../../kabs/li_5000_t5e-3_uw1116.nc'
          kapname(11) = 'Rb'  ; fktab(11) = '../../kabs/rb_5000_t5e-3_uw1116.nc'
          kapname(12) = 'Cs'  ; fktab(12) = '../../kabs/cs_5000_t5e-3_uw1116.nc'
          kapname(13) = 'TiO' ; fktab(13) = '../../kabs/tio_5000_t5e-3_uw1116.nc'
          kapname(14) = 'VO'  ; fktab(14) = '../../kabs/vo_5000_t5e-3_uw1116_vomyt.nc'
          !kapname(14) = 'VO'  ; fktab(14) = '../../kabs/vo_5000_t5e-3_uw1116.nc'
          kapname(15) = 'FeH'  ; fktab(15) = '../../kabs/feh_5000_t5e-3_uw1116.nc'
          kapname(16) = 'PH3' ; fktab(16) = '../../kabs/ph3_5000_t5e-3_uw1116.nc'
          kapname(17) = 'H2S' ; fktab(17) = '../../kabs/h2s_5000_t5e-3_uw1116.nc'
          kapname(18) = 'HCN' ; fktab(18) = '../../kabs/hcn_5000_t5e-3_uw1116.nc'
          kapname(19) = 'C2H2'; fktab(19) = '../../kabs/c2h2_5000_t5e-3_uw1116.nc'
          kapname(20) = 'SO2'; fktab(20) = '../../kabs/so2_5000_t5e-3_uw1116.nc'
          kapname(21) = 'Fe';  fktab(21) = '../../kabs/fe_5000_t5e-3_uw1116.nc'
          !kapname(21) = 'CrH'  ; fktab(21) = '../../kabs/crh_5000_t5e-3_uw1116.nc'
          
       else if (uniform_R .and. (nband .le. 5010)) then
       
        kapname(1)  = 'H2'    ; fktab(1)  = '../../kabs/h2-h2_5010_-t5e-3_uw_R1000.nc'
		  kapname(2)  = 'He'    ; fktab(2)  = '../../kabs/h2-he_5010_-t5e-3_uw_R1000.nc'
		  kapname(3)  = 'H2O'   ; fktab(3)  = '../../kabs/h2o_5010_-t5e-3_uw_R1000.nc'
		  kapname(4)  = 'CO2'   ; fktab(4)  = '../../kabs/co2_5010_-t5e-3_uw_R1000.nc'
		  kapname(5)  = 'CO'    ; fktab(5)  = '../../kabs/co_5010_-t5e-3_uw_R1000.nc'
		  kapname(6)  = 'CH4'   ; fktab(6)  = '../../kabs/ch4_5010_-t5e-3_uw_R1000.nc'
		  kapname(7)  = 'NH3'   ; fktab(7)  = '../../kabs/nh3_5010_-t5e-3_uw_R1000.nc'
		  kapname(8)  = 'Na'    ; fktab(8)  = '../../kabs/na_5010_-t5e-3_uw_R1000.nc'
		  kapname(9)  = 'K'     ; fktab(9)  = '../../kabs/k_5010_-t5e-3_uw_R1000.nc'
		  kapname(10)  = 'Li'   ; fktab(10)  = '../../kabs/li_5010_-t5e-3_uw_R1000.nc'
		  kapname(11)  = 'Rb'   ; fktab(11)  = '../../kabs/rb_5010_-t5e-3_uw_R1000.nc'
		  kapname(12)  = 'Cs'   ; fktab(12)  = '../../kabs/cs_5010_-t5e-3_uw_R1000.nc'
		  kapname(13)  = 'TiO'  ; fktab(13)  = '../../kabs/tio_5010_-t5e-3_uw_R1000.nc'
		  kapname(14)  = 'VO'   ; fktab(14)  = '../../kabs/vo_5010_-t5e-3_uw_R1000.nc'
		  kapname(15)  = 'FeH'  ; fktab(15)  = '../../kabs/feh_5010_-t5e-3_uw_R1000.nc'
		  kapname(16)  = 'PH3'  ; fktab(16)  = '../../kabs/ph3_5010_-t5e-3_uw_R1000.nc'
		  kapname(17)  = 'H2S'  ; fktab(17)  = '../../kabs/h2s_5010_-t5e-3_uw_R1000.nc'
		  kapname(18)  = 'HCN'  ; fktab(18)  = '../../kabs/hcn_5010_-t5e-3_uw_R1000.nc'
		  kapname(19)  = 'C2H2' ; fktab(19)  = '../../kabs/c2h2_5010_-t5e-3_uw_R1000.nc'
		  kapname(20)  = 'SO2'  ; fktab(20)  = '../../kabs/so2_5010_-t5e-3_uw_R1000.nc'
		  kapname(21)  = 'Fe'   ; fktab(21)  = '../../kabs/fe_5010_-t5e-3_uw_R1000.nc'
		  !kapname(22) = 'CrH'  ; fktab(22) = '../../kabs/crh_5010_-t5e-3_uw_R1000.nc'
       
      else if (uniform_R .and. (nband .le. 25256)) then

         kapname(1)  = 'H2'    ; fktab(1)  = '../../kabs/h2-h2_25256_-t5e-3_uw_R10000.nc'
         kapname(2)  = 'He'    ; fktab(2)  = '../../kabs/h2-he_25256_-t5e-3_uw_R10000.nc'
         kapname(3)  = 'H2O'   ; fktab(3)  = '../../kabs/h2o_25256_-t5e-3_uw_R10000.nc'
         kapname(4)  = 'CO2'   ; fktab(4)  = '../../kabs/co2_25256_-t5e-3_uw_R10000.nc'
         kapname(5)  = 'CO'    ; fktab(5)  = '../../kabs/co_25256_-t5e-3_uw_R10000.nc'
         kapname(6)  = 'CH4'   ; fktab(6)  = '../../kabs/ch4_25256_-t5e-3_uw_R10000.nc'
         kapname(7)  = 'NH3'   ; fktab(7)  = '../../kabs/nh3_25256_-t5e-3_uw_R10000.nc'
         kapname(8)  = 'Na'    ; fktab(8)  = '../../kabs/na_25256_-t5e-3_uw_R10000.nc'
         kapname(9)  = 'K'     ; fktab(9)  = '../../kabs/k_25256_-t5e-3_uw_R10000.nc'
         kapname(10)  = 'Li'   ; fktab(10)  = '../../kabs/li_25256_-t5e-3_uw_R10000.nc'
         kapname(11)  = 'Rb'   ; fktab(11)  = '../../kabs/rb_25256_-t5e-3_uw_R10000.nc'
         kapname(12)  = 'Cs'   ; fktab(12)  = '../../kabs/cs_25256_-t5e-3_uw_R10000.nc'
         kapname(13)  = 'TiO'  ; fktab(13)  = '../../kabs/tio_25256_-t5e-3_uw_R10000.nc'
         kapname(14)  = 'VO'   ; fktab(14)  = '../../kabs/vo_25256_-t5e-3_uw_R10000.nc'
         kapname(15)  = 'FeH'  ; fktab(15)  = '../../kabs/feh_25256_-t5e-3_uw_R10000.nc'
         kapname(16)  = 'PH3'  ; fktab(16)  = '../../kabs/ph3_25256_-t5e-3_uw_R10000.nc'
         kapname(17)  = 'H2S'  ; fktab(17)  = '../../kabs/h2s_25256_-t5e-3_uw_R10000.nc'
         kapname(18)  = 'HCN'  ; fktab(18)  = '../../kabs/hcn_25256_-t5e-3_uw_R10000.nc'
         kapname(19)  = 'C2H2' ; fktab(19)  = '../../kabs/c2h2_25256_-t5e-3_uw_R10000.nc'
         kapname(20)  = 'SO2'  ; fktab(20)  = '../../kabs/so2_25256_-t5e-3_uw_R10000.nc'
         kapname(21)  = 'Fe'   ; fktab(21)  = '../../kabs/fe_25256_-t5e-3_uw_R10000.nc'
         !kapname(22) = 'CrH'  ; fktab(22) = '../../kabs/crh_25256_-t5e-3_uw_R10000.nc'
                    
       endif
          
    endif
    

    kapcst  = 0.0

    ! read grid namelist in the file fparam
    open(1,file=fparam)
    rewind(1) ; read(1,NML=opacity,iostat=stat)
     IF (stat/=0) THEN
      BACKSPACE(1)
      READ(1,fmt='(A)') line
      WRITE(*,'(A)') 'Invalid line in namelist: '//trim(line)
      STOP
     END IF    
    close(1)
    
    ! convert in code units
    dnuktab   = dnuktab   * unitnuktab/(1./unitL)
    nuktabmin = nuktabmin * unitnuktab/(1./unitL)
    nuktabmax = nuktabmax * unitnuktab/(1./unitL)

    Tktabmin = Tktabmin*unitKktab/unitK
    Pktabmin = Pktabmin*unitPktab/unitP

    Tktabmax = Tktabmax*unitKktab/unitK
    Pktabmax = Pktabmax*unitPktab/unitP

    ppmin_smooth = ppmin_smooth*unitPsmooth/unitP
    ppmax_smooth = ppmax_smooth*unitPsmooth/unitP

    ppmin_cloud  = ppmin_cloud*unitPcloud/unitP
    ppmax_cloud  = ppmax_cloud*unitPcloud/unitP

    if (corr_k) then
       
       call init_corr_k

    else
       
       call init_lbl

    endif
    
  end subroutine init_opacity


  subroutine init_lbl

    use mod_param
    use mod_cst
    use mod_chem
    use mod_grid
    use mod_util

   implicit none
   include 'netcdf.inc'
   
#if PARA==1 

  include 'mpif.h'
  integer :: info

#endif


    !local variable
    integer :: stat,ikap,im,ip,it,i,iiloc,iktab,iktab_loc,nbin,iktab_start,iktab_end,pt_pair
    real    :: dlogpktab,dlogtktab,ktabloc_inv,kplanck,kross,tmin_calc_table,tmax_calc_table,pmin_calc_table,pmax_calc_table

    !ncdf variable
    integer,dimension(2) :: istart,iend,icount
    integer :: id_var,id_file

    !local opacity table
    real,dimension(:,:),allocatable  :: ktab_loc,kstdtab_loc

    ! PT grid in opacity tables
    real,dimension(:), allocatable :: t_calc_table, p_calc_table

    ! looks at the number of tables used
    ntab = 0
    do ikap = 1,nkap

       if (kapcst(ikap) == 0. .and. .not. fktab(ikap)=='') then
          ntab = ntab + 1
       endif

    enddo
    
    ! test if the frequencies are compatible with the opacity tables
    ! for "simplicity" it is assumed that numin and numax are always at
    ! interfaces of the different opacity tables. 
    ! The resolution can be different if dnu/dnutab or dnutab/dnu is an integer
   
    if(.not. ntab == 0) then

       do ikap = 1,nkap
          
          if (kapcst(ikap) == 0. .and. .not. fktab(ikap)=='') then

             !test if dnu/dnutab or dnutab/dnu is an integer
             if (dnu .gt. dnuktab(ikap)) then
                if (abs(dnu/dnuktab(ikap)-int(dnu/dnuktab(ikap))) >1.0E-13 ) then
                   if (mype==cputerm) then
                      
                      write(*,headrfmt) 'Error init_lbl'
                      write(*,'(A)') 'Frequency resolution incompatible with tables'
                      write(*,'(A)') 'needs dnu/dnuktab integer'
                      write(*,exppmfmt)'dnu',dnu
                      write(*,intpmfmt)'ikap',ikap
                      write(*,exppmfmt)'dnuktab',dnuktab(ikap)
                      write(*,*) 
                      write(*,separfmt)
                      stop
                   endif
                endif
             else
                if (abs(dnuktab(ikap)/dnu-int(dnuktab(ikap)/dnu)) >1.0E-13 ) then
                   if (mype==cputerm) then
                      
                      write(*,headrfmt) 'Error init_lbl'
                      write(*,'(A)') 'Frequency resolution incompatible with tables'
                      write(*,'(A)') 'needs dnuktab/dnu integer'
                      write(*,exppmfmt)'dnu',dnu
                      write(*,intpmfmt)'ikap',ikap
                      write(*,exppmfmt)'dnuktab',dnuktab(ikap)
                      write(*,*) 
                      write(*,separfmt)
                      stop
                   endif
                endif
             endif
                

          
             !test if numin or nuktabmin is on interfaces
             if (numin .lt. nuktabmin(ikap)) then
                
                if (abs((nuktabmin(ikap)-numin)/dnu-int((nuktabmin(ikap)-numin)/dnu)) >1.0E-13 ) then
                   
                   if (mype==cputerm) then
                      
                      write(*,headrfmt) 'Error init_lbl'
                      write(*,'(A)') 'numin incompatible with tables'
                      write(*,'(A)') 'needs (nuktabmin-numin)/dnu integer'
                      write(*,exppmfmt)'numin',numin
                      write(*,exppmfmt)'dnu',dnu
                      write(*,intpmfmt)'ikap',ikap
                      write(*,exppmfmt)'nuktabmin',nuktabmin(ikap)
                      write(*,*) 
                      write(*,separfmt)
                      stop
                   endif
                endif

             else if ((numin .gt. nuktabmin(ikap)).and. (numin .gt. nuktabmax(ikap))) then
                if(abs((numin-nuktabmin(ikap))/dnuktab(ikap)-int((numin-nuktabmin(ikap))/dnuktab(ikap))) >1.0E-13 ) then
                  
                   if (mype==cputerm) then
                      
                      write(*,headrfmt) 'Error init_lbl'
                      write(*,'(A)') 'numin incompatible with tables'
                      write(*,'(A)') 'needs (numin-nuktabmin)/dnuktab integer'
                      write(*,exppmfmt)'numin',numin
                      write(*,intpmfmt)'ikap',ikap
                      write(*,exppmfmt)'dnuktab',dnuktab(ikap)
                      write(*,exppmfmt)'nuktabmin',nuktabmin(ikap)
                      write(*,*) 
                      write(*,separfmt)
                      stop
                   endif
                  
                   
                endif
             endif

             !test if numax or nuktabmax is on interfaces
             if (numax .gt. nuktabmax(ikap)) then
                
                if (abs((numax-nuktabmax(ikap))/dnu-int((numax-nuktabmax(ikap))/dnu)) >1.0E-13 ) then
                   
                   if (mype==cputerm) then
                      
                      write(*,headrfmt) 'Error init_lbl'
                      write(*,'(A)') 'numax incompatible with tables'
                      write(*,'(A)') 'needs (numax-nuktabmax)/dnu integer'
                      write(*,exppmfmt)'numax',numax
                      write(*,exppmfmt)'dnu',dnu
                      write(*,intpmfmt)'ikap',ikap
                      write(*,exppmfmt)'nuktabmax',nuktabmax(ikap)
                      write(*,*) 
                      write(*,separfmt)
                      stop
                   endif
                endif

             else if ((numax .lt. nuktabmax(ikap)).and. (numax.gt. nuktabmin(ikap))) then
                if(abs((nuktabmax(ikap)-numax)/dnuktab(ikap)-int((nuktabmax(ikap)-numax)/dnuktab(ikap))) >1.0E-13 ) then
                  
                   if (mype==cputerm) then
                      
                      write(*,headrfmt) 'Error init_lbl'
                      write(*,'(A)') 'numax incompatible with tables'
                      write(*,'(A)') 'needs (nuktabmax-numax)/dnuktab integer'
                      write(*,exppmfmt)'numax',numax
                      write(*,intpmfmt)'ikap',ikap
                      write(*,exppmfmt)'dnuktab',dnuktab(ikap)
                      write(*,exppmfmt)'nuktabmax',nuktabmax(ikap)
                      write(*,*) 
                      write(*,separfmt)
                      stop
                   endif
                  
                   
                endif
             endif
             
          endif
       enddo   
       
    endif




    if (mype == cputerm) then
       write(*,headrfmt) 'Opacity'
       write(*,intpmfmt) 'Nkap',nkap,'Number of molecules used for opacities'
    endif
   
    ! if table are used:
    if (.not. ntab == 0) then
      
       ! create pressure and temperature points for tables
       dlogTktab = (log10(tktabmax)-log10(tktabmin))/(ntktab-1)
       dlogPktab = (log10(pktabmax)-log10(pktabmin))/(npktab-1)
       
       !allocate arrays
       allocate(pktab(npktab),tktab(ntktab))
       
       Tktab(1) = tktabmin
       Pktab(1) = pktabmin
       
       do iT = 2,ntktab
          Tktab(iT) = Tktab(iT-1)*10.**dlogTktab
       enddo
       
       do iP = 2,npktab
          Pktab(iP) = Pktab(iP-1)*10.**dlogPktab
       enddo
      
       
       !create global frequency range for which table opacity exists
       nukapmin_glob = numax
       nukapmax_glob = numin
       
       do ikap = 1,nkap
          if (kapcst(ikap) == 0.) then
             
             !look for the maximum range in the opacity tables
             nukapmin_glob = min(nukapmin_glob,nuktabmin(ikap))
             nukapmax_glob = max(nukapmax_glob,nuktabmax(ikap))
             
          endif
       enddo
       
       !do the intersection with the actual band required
       nukapmin_glob = max(numin,nukapmin_glob)
       nukapmax_glob = min(numax,nukapmax_glob)
       
       !global number of frequency points with a table opacity
       nfreqkap_glob = int((nukapmax_glob-nukapmin_glob)/dnu)
       
       !create the cpu local frequency range for the opacity table
       inukap_start = mype*ceiling(real(nfreqkap_glob)/ncpu)+1
       inukap_end   = (mype+1)*ceiling(real(nfreqkap_glob)/ncpu)
       
       if (inukap_end > nfreqkap_glob) inukap_end = Nfreqkap_glob
       nfreqkap = inukap_end - inukap_start + 1
      
       nukapmin = nukapmin_glob + (inukap_start-1)*dnu
       nukapmax = nukapmin_glob + (inukap_end)*dnu

       allocate(nuf(nfreqkap+1),nu(nfreqkap))
       allocate(fnu(ndepth+1,nfreqkap),qnu(ndepth+1,nfreqkap),fstar_nu(ndepth+1,nfreqkap))
       allocate(kapnu(ndepth,nfreqkap))
       allocate(kapfnu(ndepth+1,nfreqkap))
       allocate(dkapnu(2*ndepth,nfreqkap))
       allocate(dkapfnu(2*(ndepth+1),nfreqkap))
       allocate(kscatnu(ndepth,nfreqkap),epsnu(ndepth,nfreqkap))
       allocate(kscatfnu(ndepth+1,nfreqkap),epsfnu(ndepth+1,nfreqkap))
       if (um_input) allocate(kapfnu_um(nlong,nlat,ndepth+1,nfreqkap))

       kscatnu  = 0.
       kscatfnu = 0.
       epsnu   = 1. 
       epsfnu  = 1. 

       nuf(1) = nukapmin
       nuf(nfreqkap+1) = nukapmax
       nu(1)  = nukapmin + 0.5*dnu
       
       
       do i=2,nfreqkap

          nuf(i) = nuf(i-1) + dnu 
          nu(i)  = nuf(i) + 0.5*dnu
       
       enddo
       
       !allocate the local opacity table for the local frequency range
       allocate(ktab(ntktab,npktab,nfreqkap,nkap))
       allocate(kstdtab(ntktab,npktab,1,nkap))

       ktab = 0.
       
       do ikap = 1,nkap

          if (mype==cputerm) then
             !verify if the molecule is in the chemistry
             im = imol(kapname(ikap))
 
             if (kapcst(ikap)==0.) then
                if (fktab(ikap) == '') then
                   write(*,'(i6,5x,A9,A13,es10.3,A12)') ikap,kapname(ikap),' cst opacity ',kapcst(ikap),' [m^2 kg-1]'
                else
                   write(*,'(i6,5x,A9,A50)') ikap,kapname(ikap),fktab(ikap)
                endif
             else
                write(*,'(i6,5x,A9,A13,es10.3,A12)') ikap,kapname(ikap),' cst opacity ',kapcst(ikap),' [m^2 kg-1]'
             endif
          endif

          if(.not. kapcst(ikap) == 0. .or. (kapcst(ikap)==0. .and. fktab(ikap)=='')) then
             !cst opacity
             ktab(:,:,:,ikap) = kapcst(ikap)
             kstdtab(:,:,1,ikap) = kapcst(ikap)

          else
             !look at the table
             
             !1. if the local range is outside the table the opacity is 0.
             if(nukapmax .le. nuktabmin(ikap)) then
                ktab(:,:,:,ikap) = 0.
             else if (nukapmin.ge. nuktabmax(ikap)) then
                ktab(:,:,:,ikap) = 0.             
             else
                !2. assumes first the local range is fully inside the opacity table range
                !recall: nutab = nutabmin +(i-1)*dnutab
                istart(1) = floor((nukapmin-nuktabmin(ikap))/dnuktab(ikap))+1
                iend(1)   = ceiling((nukapmax-nuktabmin(ikap))/dnuktab(ikap))
              
                !3. if there is only an overlap then take everything you can from the table
                if (iend(1).gt. nfreqktab(ikap)) then
                   iend(1) = nfreqktab(ikap)
                else if (istart(1).lt.1) then
                   istart(1) = 1
                endif
               
                istart(2) = 1
                iend(2)   = npktab*ntktab
                
                icount(:) = iend(:) - istart(:) + 1
                
                !allocate and read the opacity table
                
                allocate(ktab_loc(icount(1),icount(2)))
                
                call nf(nf_open(fktab(ikap),nf_nowrite,id_file))

                 call nf(nf_inq_dimid(id_file,'pt_pair',id_var))
                 call nf(nf_inq_dimlen(id_file,id_var,pt_pair))

                 ! Check number of pt pairs in k-table matches
                 if (pt_pair/=(npktab*ntktab)) then 
                  write(*,headrfmt) 'Error init_corr_k'
                  write(*,*) 'Incorrect number of pt pairs in k-table'
                  write(*,intpmfmt)'pt_pair_table',pt_pair
                  write(*,intpmfmt)'pt_pair_code',npktab*ntktab
                  write(*,*) 
                  write(*,separfmt)
                  stop
                 endif

                 ! Check PT table 
                 allocate(t_calc_table(pt_pair), p_calc_table(pt_pair))
                 call nf(nf_inq_varid(id_file,'p_calc',id_var))
                 call nf(nf_get_vara_double(id_file,id_var,1,pt_pair,p_calc_table))
                 call nf(nf_inq_varid(id_file,'t_calc',id_var))
                 call nf(nf_get_vara_double(id_file,id_var,1,pt_pair,t_calc_table))

                 tmin_calc_table = minval(t_calc_table)*unitKktab/unitK
                 tmax_calc_table = maxval(t_calc_table)*unitKktab/unitK
                 pmin_calc_table = minval(p_calc_table)*10.0*unitPktab/unitP
                 pmax_calc_table = maxval(p_calc_table)*10.0*unitPktab/unitP

                 ! Check temperature grid
                 if ((tmin_calc_table/=tktabmin) .or. (tmax_calc_table/=tktabmax)) then 
                  write(*,headrfmt) 'Error init_corr_k'
                  write(*,*) 'Incorrect temperature grid in k-table'
                  write(*,fltpmfmt) 'Tmin k-table',tmin_calc_table
                  write(*,fltpmfmt) 'Tmax k-table',tmax_calc_table
                  write(*,fltpmfmt) 'Tmin code',tktabmin
                  write(*,fltpmfmt) 'Tmax code',tktabmax
                  write(*,*)
                  write(*,separfmt)
                  stop
                 endif

                 ! Check pressure grid - note pmin is slightly different between the two, so we shall check 
                 ! the percentage difference. If pmin or pmax differ by more than 0.1% then stop
                if ((ABS(((pmin_calc_table-pktabmin)/pmin_calc_table)*100.0)>=0.1) .or. &
                  & (ABS(((pmax_calc_table-pktabmax)/pmax_calc_table)*100.0)>=0.1)) then 
                  write(*,headrfmt) 'Error init_corr_k'
                  write(*,*) 'Incorrect pressure grid in k-table'
                  write(*,exppmfmt) 'Pmin k-table',pmin_calc_table
                  write(*,exppmfmt) 'Pmax k-table',pmax_calc_table
                  write(*,exppmfmt) 'Pmin code',pktabmin
                  write(*,exppmfmt) 'Pmax code',pktabmax
                  write(*,*)
                  write(*,separfmt)
                  stop
                 endif

                 deallocate(p_calc_table, t_calc_table)



                call nf(nf_inq_varid(id_file,'kabs',id_var))
                call nf(nf_get_vara_double(id_file,id_var,istart,icount,ktab_loc))
              
                if (mype == cputerm) then
                   istart(1) = int(nustd/dnuktab(ikap))
                   iend(1)   = int(nustd/dnuktab(ikap))
                   icount(1) = 1
                   allocate(kstdtab_loc(icount(1),icount(2)))
               
                   
                   call nf(nf_get_vara_double(id_file,id_var,istart,icount,kstdtab_loc))
                endif
                call nf(nf_close(id_file))
             
                !frequency indices to fill with the opacity table
                iktab_start = 1
                iktab_end   = nfreqkap
               
                !deal with overlaps                
                if (nukapmax .ge. nuktabmax(ikap)) then

                   iktab_end   = nfreqkap - (nukapmax-nuktabmax(ikap))/dnu   

                endif

                if (nukapmin.le.nuktabmin(ikap)) then

                    iktab_start = (nuktabmin(ikap)-nukapmin)/dnu+1

                 endif
                  
                do i = istart(2),iend(2)
                   iP = (i-1)/ntktab + 1
                   iT = i - (iP-1)*ntktab
                   
                   ktab(iT,iP,:,ikap) = 0.
       
                   if (mype==cputerm) kstdtab(iT,iP,1,ikap) = kstdtab_loc(1,i)
                   
                   if (dnu == dnuktab(ikap)) then
                      ktab(iT,iP,iktab_start:iktab_end,ikap) = ktab_loc(:,i)
                      
                   else if (dnu<dnuktab(ikap)) then
                      !increase the resolution of the file by direct injection
                      nbin = int(dnuktab(ikap)/dnu)

                      do iktab = iktab_start,iktab_end
                         iktab_loc = int((iktab-iktab_start)/nbin)+1
                         ktab(iT,iP,iktab,ikap) = ktab_loc(iktab_loc,i)
                      enddo

                   else if(dnu>dnuktab(ikap)) then
                      !decrease the resolution of the file by rosseland or planck mean (assumes that dBnu(T)/dT cst on the bin)
                      nbin = int(dnu/dnuktab(ikap))

                      iktab_loc = 0
                      do iktab = iktab_start,iktab_end
                         ktabloc_inv = 0.
                         do iiloc = iktab_loc+1,iktab_loc+nbin
                            if (.not. ktab_loc(iiloc,i) == 0.) then
                               ktabloc_inv = ktabloc_inv + 1./ktab_loc(iiloc,i)
                            endif
                         enddo
                         
                         if( .not. ktabloc_inv == 0.) then
                            kross = nbin/ktabloc_inv
                         else
                            kross = 0.
                         endif
                      
                         kplanck = 0.
                         do iiloc = iktab_loc+1,iktab_loc+nbin
                            
                            kplanck = kplanck+ ktab_loc(iiloc,i)
                         enddo
                         
                         kplanck = kplanck/nbin


                         ktab(it,ip,iktab,ikap) = max(kplanck,kross)
   
                         iktab_loc = iktab_loc + nbin
                         
                      enddo

                      
                   endif
                   
                enddo
                
                deallocate(ktab_loc)   
                if (mype==cputerm) deallocate(kstdtab_loc)   
                
             endif
             
             
             
          endif
       enddo
       
       !convert ktab units to cgs to code unit

       ktab = ktab*(unitLktab**2/unitMktab)/(unitL**2/unitM)
       kstdtab = kstdtab*(unitLktab**2/unitMktab)/(unitL**2/unitM)
       

#if PARA==1
       !share the standard opacity table among cpu
       call MPI_bcast(kstdtab,npktab*ntktab*nkap,mpi_double_precision,cputerm,mpi_comm_world,info)
#endif
       


   endif

    ! if no table is used, the opacity is cst for all the species
    if (ntab == 0 .and. nkap>1) then
    
       npktab   = 1
       ntktab   = 1
       nfreqkap = 1
       
       allocate(ktab(ntktab,npktab,nfreqkap,nkap))
       allocate(kstdtab(ntktab,npktab,1,nkap))
       
       ktab = 0.
       
       do ikap = 1,nkap
          
          ktab(1,1,1,ikap) = kapcst(ikap)
          kstdtab(1,1,1,ikap) = kapcst(ikap)
          
       enddo

       !create the cpu local frequency range for the opacity table
       nfreqkap_glob = nfreq
       
       nukapmin_glob = numin
       nukapmax_glob = numax

       inukap_start = mype*(nfreqkap_glob/ncpu)+1
       inukap_end   = (mype+1)*(nfreqkap_glob/ncpu)
       
       if (inukap_end > nfreqkap_glob) inukap_end = Nfreqkap_glob
       nfreqkap = inukap_end - inukap_start + 1
       
       nukapmin = nukapmin_glob + (inukap_start-1)*dnu
       nukapmax = nukapmin_glob + (inukap_end)*dnu
       
       allocate(nuf(nfreqkap+1),nu(nfreqkap))
       allocate(fnu(ndepth+1,nfreqkap),qnu(ndepth+1,nfreqkap),fstar_nu(ndepth+1,nfreqkap))
       allocate(kapnu(ndepth,nfreqkap))
       allocate(kapfnu(ndepth+1,nfreqkap))
       allocate(dkapnu(2*ndepth,nfreqkap))
       allocate(dkapfnu(2*(ndepth+1),nfreqkap))

       nuf(1) = nukapmin
       nuf(nfreqkap+1) = nukapmax
       nu(1)  = nuf(1) + 0.5*dnu
       
       do i=2,nfreqkap
          
          nuf(i) = nuf(i-1) + dnu 
          nu(i)  = nuf(i) + 0.5*dnu

       enddo

       !convert ktab units to cgs to code unit
       
       ktab = ktab*(unitLkcst**2/unitMkcst)/(unitL**2/unitM)
       kstdtab = kstdtab*(unitLkcst**2/unitMkcst)/(unitL**2/unitM)
    

    endif
    
    
    if (mype==cputerm) then
       write(*,*)
       write(*,logpmfmt) 'kap_smooth',kap_smooth,'if T opacities are vertically smoothed'
       if (kap_smooth) then
          write(*,exppmfmt) 'ppmin_smooth [bar]',ppmin_smooth*unitP/1.0E6, 'pressure above which kappa is smoothed'
          write(*,exppmfmt) 'ppmax_smooth [bar]',ppmax_smooth*unitP/1.0E6, 'pressure below which kappa is smoothed'
          write(*,intpmfmt)'kerkap_smooth',kerkap_smooth,'number of neighbors used for smoothing'
       endif
       write(*,*)  
       write(*,separfmt)
    endif

  end subroutine init_lbl

  subroutine init_corr_k

    use mod_param
    use mod_cst
    use mod_chem
    use mod_grid
    use mod_util

   implicit none
   include 'netcdf.inc'
   
#if PARA==1 

  include 'mpif.h'
  integer :: info

#endif


    !local variable
    integer :: stat,ikap,im,ip,it,i,iiloc,iktab,iktab_loc,nbin,iktab_start,iktab_end,ib,ib_loc,imax,pt_pair
    real    :: dlogpktab,dlogtktab,ktabloc_inv,kplanck,kross,tmin_calc_table,tmax_calc_table,pmin_calc_table,pmax_calc_table

    !ncdf variable
    integer,dimension(3) :: istart,iend,icount
    integer :: id_var,id_file, nband_loc

    !local opacity table
    real,dimension(:,:,:),allocatable  :: ktab_loc,kstdtab_loc
    real,dimension(:,:),allocatable  :: weights_loc
    integer,dimension(:    ),allocatable  :: iband,n_k
    real,dimension(:), allocatable :: t_calc_table, p_calc_table

    itauref = 0
    do i=1,ndepth+1

       if (taufstd(i)<1.) then
          itauref = i
       endif

    enddo
   

    ! Corr-k has to be used for all species, no cst opacity allowed
    do ikap = 1,nkap

       if (kapcst(ikap) /= 0. .or. fktab(ikap)=='') then
          
          if (mype==cputerm) then
             
             write(*,headrfmt) 'Error init_corr_k'
             write(*,'(A)') 'kapcst /= 0 or fktab = None'
             write(*,intpmfmt)'ikap',ikap
             write(*,*) 
             write(*,separfmt)
             stop
          
          endif

       endif

    enddo
    

    if (mype == cputerm) then
       write(*,headrfmt) 'Opacity'
       write(*,intpmfmt) 'Nkap',nkap,'Number of molecules used for opacities'
    endif
    
    ! create pressure and temperature points for tables
    dlogTktab = (log10(tktabmax)-log10(tktabmin))/(ntktab-1)
    dlogPktab = (log10(pktabmax)-log10(pktabmin))/(npktab-1)
    
    !allocate arrays
    allocate(pktab(npktab),tktab(ntktab))
    
    Tktab(1) = tktabmin
    Pktab(1) = pktabmin
       
    do iT = 2,ntktab
       Tktab(iT) = Tktab(iT-1)*10.**dlogTktab
    enddo
       
    do iP = 2,npktab
       Pktab(iP) = Pktab(iP-1)*10.**dlogPktab
    enddo
      
    nfreqkap = nkmix*nband

    allocate(nuf(nfreqkap+1),nu(nfreqkap))
    allocate(fnu(ndepth+1,nfreqkap),qnu(ndepth+1,nfreqkap),fstar_nu(ndepth+1,nfreqkap))
    allocate(kapnu(ndepth,nfreqkap))
    allocate(kapfnu(ndepth+1,nfreqkap))
    allocate(weightsnu(ndepth,nfreqkap))
    allocate(weightsfnu(ndepth+1,nfreqkap))
    allocate(dkapnu(2*ndepth,nfreqkap))
    allocate(dkapfnu(2*(ndepth+1),nfreqkap))
    allocate(kscatnu(ndepth,nfreqkap),epsnu(ndepth,nfreqkap))
    allocate(kscatfnu(ndepth+1,nfreqkap),epsfnu(ndepth+1,nfreqkap))
    if (um_input) allocate(kapfnu_um(nlong,nlat,ndepth+1,nfreqkap))

    
    if(ndim==2) then
         
       allocate(kapnu_2d(ndepth,0:ntheta+2,nfreqkap))   !< opacities in 2d
       allocate(epsfnu_2d(ndepth+1,0:ntheta+2,nfreqkap))  !< epsfnu in 2d
       allocate(dsnu_2d(ndepth+1,0:ntheta+2,nfreqkap))    !<  source term in 2d
       kapnu_2d = 0.
       epsfnu_2d = 1.
       dsnu_2d = 0.
       
    endif

    kscatnu  = 0.
    kscatfnu = 0.
    epsnu    = 1. 
    epsfnu   = 1.

    nuf = 0. ; nu = 0.
       
    !allocate the local opacity table for the local frequency range
    allocate(ktab(ntktab,npktab,nk*nband,nkap))
    allocate(weights(nk*nband,nkap))
    allocate(nkmax(nband,nkap))
    !allocate(kstdtab(ntktab,npktab,1,nkap))

    
    weights = 0.
    ktab = 0.
    nkmax = 0
    
       
    do ikap = 1,nkap

      
       
       if (mype==cputerm) then
          !verify if the molecule is in the chemistry
          im = imol(kapname(ikap))
             
          if (kapcst(ikap)==0.) then
             if (fktab(ikap) == '') then
                write(*,'(i6,5x,A9,A13,es10.3,A12)') ikap,kapname(ikap),' cst opacity ',kapcst(ikap),' [m^2 kg-1]'
             else
                write(*,'(i6,5x,A9,A40)') ikap,kapname(ikap),fktab(ikap)
             endif
          else
             write(*,'(i6,5x,A9,A13,es10.3,A12)') ikap,kapname(ikap),' cst opacity ',kapcst(ikap),' [m^2 kg-1]'
          endif
       endif
       
       nband_loc = 0
       call nf(nf_open(fktab(ikap),nf_nowrite,id_file))
       call nf(nf_inq_dimid (id_file,'band',id_var))
       call nf(nf_inq_dimlen(id_file,id_var,nband_loc))
       call nf(nf_inq_dimid(id_file,'pt_pair',id_var))
       call nf(nf_inq_dimlen(id_file,id_var,pt_pair))

       ! Check number of pt pairs in k-table matches
       if (pt_pair/=(npktab*ntktab)) then 
        write(*,headrfmt) 'Error init_corr_k'
        write(*,*) 'Incorrect number of pt pairs in k-table'
        write(*,intpmfmt)'pt_pair_table',pt_pair
        write(*,intpmfmt)'pt_pair_code',npktab*ntktab
        write(*,*) 
        write(*,separfmt)
        stop
       endif

       istart(1) = 1
       iend(1)   = nk

       istart(2) = 1
       iend(2)   = npktab*ntktab
       
       istart(3) = 1
       iend(3)   = nband_loc
          
       icount(:) = iend(:) - istart(:) + 1
       
       !allocate and read the opacity table
          
       allocate(ktab_loc(icount(1),icount(2),icount(3)))
     
       call nf(nf_inq_varid(id_file,'kopt',id_var))
       call nf(nf_get_vara_double(id_file,id_var,istart,icount,ktab_loc))
     
       allocate(iband(icount(3)))
       allocate(n_k(icount(3)))
       call nf(nf_inq_varid(id_file,'band',id_var))
       call nf(nf_get_vara_int(id_file,id_var,istart(3),icount(3),iband))
     
       allocate(weights_loc(icount(1),icount(3)))

       call nf(nf_inq_varid(id_file,'w_k',id_var))
       call nf(nf_get_vara_double(id_file,id_var,(/istart(1),istart(3)/),(/icount(1),icount(3)/),weights_loc))

       call nf(nf_inq_varid(id_file,'n_k',id_var))
       call nf(nf_get_vara_int(id_file,id_var,istart(3),icount(3),n_k))

       ! Check PT table 
       allocate(t_calc_table(pt_pair), p_calc_table(pt_pair))
       call nf(nf_inq_varid(id_file,'p_calc',id_var))
       call nf(nf_get_vara_double(id_file,id_var,1,pt_pair,p_calc_table))
       call nf(nf_inq_varid(id_file,'t_calc',id_var))
       call nf(nf_get_vara_double(id_file,id_var,1,pt_pair,t_calc_table))

       tmin_calc_table = minval(t_calc_table)*unitKktab/unitK
       tmax_calc_table = maxval(t_calc_table)*unitKktab/unitK
       pmin_calc_table = minval(p_calc_table)*10.0*unitPktab/unitP
       pmax_calc_table = maxval(p_calc_table)*10.0*unitPktab/unitP

       ! Check temperature grid
       if ((tmin_calc_table/=tktabmin) .or. (tmax_calc_table/=tktabmax)) then 
        write(*,headrfmt) 'Error init_corr_k'
        write(*,*) 'Incorrect temperature grid in k-table'
        write(*,fltpmfmt) 'Tmin k-table',tmin_calc_table
        write(*,fltpmfmt) 'Tmax k-table',tmax_calc_table
        write(*,fltpmfmt) 'Tmin code',tktabmin
        write(*,fltpmfmt) 'Tmax code',tktabmax
        write(*,*)
        write(*,separfmt)
        stop
       endif

       ! Check pressure grid - note pmin is slightly different between the two, so we shall check 
       ! the percentage difference. If pmin or pmax differ by more than 0.1% then stop
      if ((ABS(((pmin_calc_table-pktabmin)/pmin_calc_table)*100.0)>=0.1) .or. &
        & (ABS(((pmax_calc_table-pktabmax)/pmax_calc_table)*100.0)>=0.1)) then 
        write(*,headrfmt) 'Error init_corr_k'
        write(*,*) 'Incorrect pressure grid in k-table'
        write(*,exppmfmt) 'Pmin k-table',pmin_calc_table
        write(*,exppmfmt) 'Pmax k-table',pmax_calc_table
        write(*,exppmfmt) 'Pmin code',pktabmin
        write(*,exppmfmt) 'Pmax code',pktabmax
        write(*,*)
        write(*,separfmt)
        stop
       endif

       deallocate(p_calc_table, t_calc_table)
       
     
       call nf(nf_close(id_file))
      
       do i = istart(2),iend(2)
          iP = (i-1)/ntktab + 1
          iT = i - (iP-1)*ntktab

          

          do ib =1,nband_loc

             ib_loc = iband(ib)
             
             if (ib_loc.le.nband) then
                iktab_start = (ib_loc-1)*nk+1
                iktab_end   = ib_loc*nk
                ktab(iT,iP,iktab_start:iktab_end,ikap) = ktab_loc(:,i,ib)
             endif
             
          enddo
          
       enddo
       
       do ib =1,nband_loc
          
          ib_loc = iband(ib)
          
          if (ib_loc.le.nband) then
             iktab_start = (ib_loc-1)*nk+1
             iktab_end   = ib_loc*nk
             
             weights(iktab_start:iktab_end,ikap) = weights_loc(:,ib)
             
          if (n_k(ib)>nk) then

              if (mype==cputerm) then          
               write(*,headrfmt) 'Error init_corr_k'
               write(*,'(A)') 'weight number > nk'
               write(*,intpmfmt)'n_k',n_k(ib)
               write(*,intpmfmt)'nk',nk
               write(*,intpmfmt)'ib',ib
               write(*,intpmfmt)'ikap',ikap
               write(*,*) 
               write(*,separfmt)
               stop
            endif

             

          endif

          
          nkmax(iband(ib),ikap) =n_k(ib)

        
          endif
       enddo

                
       deallocate(ktab_loc,iband,weights_loc,n_k)   
        
    enddo
       
       !convert ktab units to cgs to code unit

       ktab = ktab*(unitLktab**2/unitMktab)/(unitL**2/unitM)
       !kstdtab = kstdtab*(unitLktab**2/unitMktab)/(unitL**2/unitM)
       
       if (mype==cputerm) then
          write(*,*)
          write(*,logpmfmt) 'kap_smooth',kap_smooth,'if T opacities are vertically smoothed'
       if (kap_smooth) then
          write(*,exppmfmt) 'ppmin_smooth [bar]',ppmin_smooth*unitP/1.0E6, 'pressure above which kappa is smoothed'
          write(*,exppmfmt) 'ppmax_smooth [bar]',ppmax_smooth*unitP/1.0E6, 'pressure below which kappa is smoothed'
          write(*,intpmfmt)'kerkap_smooth',kerkap_smooth,'number of neighbors used for smoothing'
       endif
       write(*,*)  
       write(*,separfmt)
    endif


 

        


     end subroutine init_corr_k





  !> Fill the opacity profile of the atmosphere in function of P,T
  !<
  subroutine get_lbl_opacity

    use mod_grid
    use mod_util
    use mod_chem
    use mod_cst
    use mod_param
    use mod_eos
    
    implicit none

    integer :: ikap,idepth,ip,it,ifreq,ier

    real,dimension(:),allocatable :: logTatm,logPatm
    real,dimension(nPktab) :: logPktab
    real,dimension(nTktab) :: logTktab
    real,dimension(:,:,:),allocatable  :: kapnumol,kapstdmol
    real,dimension(ntktab,npktab,nfreqkap,nkap) :: dktabdp,dktabdt
    real,dimension(ntktab,npktab,1,nkap) :: dkstddp,dkstddt
    real,dimension(ndepth)   :: mu, cfdtaufstd
    real,dimension(ndepth+1) :: muf
    real,dimension(ndepth,nfreqkap) :: cfdtaufnu
    
    kapnu   = 0.
    kapfnu  = 0.
    dkapnu  = 0.
    dkapfnu = 0.
    kapstd = 0.
    kapfstd = 0.
    dkapstd = 0.
    dkapfstd = 0.
     
    mu  = get_mu(ndepth,Amol)
    muf = get_mu(ndepth+1,Amolf)
     
    if (.not. ntab == 0) then

       logPktab = log10(Pktab)
       logTktab = log10(Tktab)

       !compute table of derivatives
       call dfdxdy(ntktab,npktab,nfreqkap,nkap,ktab,logTktab,logPktab,dktabdt,dktabdp)

       call dfdxdy(ntktab,npktab,1,nkap,kstdtab,logTktab,logPktab,dkstddt,dkstddp)

       do iP = 1,npktab
          do iT = 1,ntktab

             dktabdt(it,ip,:,:) = dktabdt(it,ip,:,:)/(log(10.)*Tktab(it))
             dktabdp(it,ip,:,:) = dktabdp(it,ip,:,:)/(log(10.)*Pktab(ip))

             dkstddt(it,ip,:,:) = dkstddt(it,ip,:,:)/(log(10.)*Tktab(it))
             dkstddp(it,ip,:,:) = dkstddp(it,ip,:,:)/(log(10.)*Pktab(ip))

          enddo
       enddo
          


       allocate(logTatm(ndepth),logPatm(ndepth))
       allocate(kapnumol(ndepth,nfreqkap,nkap))
       allocate(kapstdmol(ndepth,1,nkap))

       

       where (pp<pktabmin)
          logPatm = log10(pktabmin)
       elsewhere
          logPatm = log10(pp)
       endwhere
       
       where (tt<tktabmin)
          logTatm = log10(tktabmin)
       elsewhere (tt>tktabmax)
          logTatm = log10(tktabmax)
       elsewhere
          logTatm = log10(tt)
       endwhere       
       
       kapstd = 0.
       kapstdmol = 0.
    
       !cubic spline interpolation
       call cubint(ntktab,npktab,nfreqkap,nkap,logTktab,logPktab,ktab,ndepth,logTatm,logPatm,kapnumol)
       call cubint(ntktab,npktab,1,nkap,logTktab,logPktab,kstdtab,ndepth,logTatm,logPatm,kapstdmol)

       ! Set floor on opacity so it cannot be negative
       where (kapnumol .lt. 0.0)
         kapnumol = 0.0
       endwhere
       where (kapstdmol .lt. 0.0)
        kapstdmol = 0.0
       endwhere

       do ikap=1,nkap

          do idepth = 1,ndepth
             kapnu(idepth,:) = kapnu(idepth,:) + kapnumol(idepth,:,ikap)*Amol(idepth,imol(kapname(ikap)))*molmass(imol(kapname(ikap)))/mu(idepth)
             kapstd(idepth) = kapstd(idepth) + kapstdmol(idepth,1,ikap)*Amol(idepth,imol(kapname(ikap)))*molmass(imol(kapname(ikap)))/mu(idepth)

          enddo
       enddo

       !cubic spline interpolation derivatives
       kapnumol = 0.
       kapstdmol = 0.

       call cubint(ntktab,npktab,nfreqkap,nkap,logTktab,logPktab,dktabdp,ndepth,logTatm,logPatm,kapnumol)
       call cubint(ntktab,npktab,1,nkap,logTktab,logPktab,dkstddp,ndepth,logTatm,logPatm,kapstdmol)

       ! Set floor on opacity so it cannot be negative
       where (kapnumol .lt. 0.0)
         kapnumol = 0.0
       endwhere
       where (kapstdmol .lt. 0.0)
        kapstdmol = 0.0
       endwhere

       do ikap=1,nkap
          do idepth = 1,ndepth
             dkapnu(idepth,:) = dkapnu(idepth,:) + kapnumol(idepth,:,ikap)*Amol(idepth,imol(kapname(ikap)))*molmass(imol(kapname(ikap)))/mu(idepth)
             
             dkapstd(idepth) = dkapstd(idepth) + kapstdmol(idepth,1,ikap)*Amol(idepth,imol(kapname(ikap)))*molmass(imol(kapname(ikap)))/mu(idepth)

          enddo
       enddo

       kapnumol = 0.
       kapstdmol = 0.

       call cubint(ntktab,npktab,nfreqkap,nkap,logTktab,logPktab,dktabdt,ndepth,logTatm,logPatm,kapnumol)
       call cubint(ntktab,npktab,1,nkap,logTktab,logPktab,dkstddt,ndepth,logTatm,logPatm,kapstdmol)

       ! Set floor on opacity so it cannot be negative
       where (kapnumol .lt. 0.0)
         kapnumol = 0.0
       endwhere
       where (kapstdmol .lt. 0.0)
        kapstdmol = 0.0
       endwhere

       do ikap=1,nkap
          do idepth = 1,ndepth
             dkapnu(ndepth+idepth,:) = dkapnu(ndepth+idepth,:) + kapnumol(idepth,:,ikap)*Amol(idepth,imol(kapname(ikap)))*molmass(imol(kapname(ikap)))/mu(idepth)

             dkapstd(ndepth+idepth) = dkapstd(ndepth+idepth) + kapstdmol(idepth,1,ikap)*Amol(idepth,imol(kapname(ikap)))*molmass(imol(kapname(ikap)))/mu(idepth)

          enddo
       enddo
       
       deallocate(kapnumol)
       deallocate(logTatm,logPatm)
       deallocate(kapstdmol)

       !compute kapfnu at cell interfaces
       
       allocate(logTatm(ndepth+1),logPatm(ndepth+1))
       allocate(kapnumol(ndepth+1,nfreqkap,nkap))
       allocate(kapstdmol(ndepth+1,1,nkap))

       
       where (ppf<pktabmin)
          logPatm = log10(pktabmin)
       elsewhere
          logPatm = log10(ppf)
       endwhere

        where (ttf<tktabmin)
          logTatm = log10(tktabmin)
       elsewhere (ttf>tktabmax)
          logTatm = log10(tktabmax)
       elsewhere
          logTatm = log10(ttf)
       endwhere 


       kapnumol = 0.
       kapstdmol = 0.
    
       !cubic spline interpolation
       call cubint(ntktab,npktab,nfreqkap,nkap,logTktab,logPktab,ktab,ndepth+1,logTatm,logPatm,kapnumol)
       call cubint(ntktab,npktab,1,nkap,logTktab,logPktab,kstdtab,ndepth+1,logTatm,logPatm,kapstdmol)

       ! Set floor on opacity so it cannot be negative
       where (kapnumol .lt. 0.0)
         kapnumol = 0.0
       endwhere
       where (kapstdmol .lt. 0.0)
        kapstdmol = 0.0
       endwhere
      
       do ikap=1,nkap

          do idepth = 1,ndepth+1
             kapfnu(idepth,:) = kapfnu(idepth,:) + kapnumol(idepth,:,ikap)*Amolf(idepth,imol(kapname(ikap)))*molmass(imol(kapname(ikap)))/muf(idepth)

             kapfstd(idepth) = kapfstd(idepth) + kapstdmol(idepth,1,ikap)*Amolf(idepth,imol(kapname(ikap)))*molmass(imol(kapname(ikap)))/muf(idepth)
             
          enddo
         
       enddo
       
       kapnumol = 0.
       kapstdmol = 0.

       call cubint(ntktab,npktab,nfreqkap,nkap,logTktab,logPktab,dktabdp,ndepth+1,logTatm,logPatm,kapnumol)
       call cubint(ntktab,npktab,1,nkap,logTktab,logPktab,dkstddp,ndepth+1,logTatm,logPatm,kapstdmol)

       ! Set floor on opacity so it cannot be negative
       where (kapnumol .lt. 0.0)
         kapnumol = 0.0
       endwhere
       where (kapstdmol .lt. 0.0)
        kapstdmol = 0.0
       endwhere

       do ikap=1,nkap
          do idepth = 1,ndepth+1
             dkapfnu(idepth,:) = dkapfnu(idepth,:) + kapnumol(idepth,:,ikap)*Amolf(idepth,imol(kapname(ikap)))*molmass(imol(kapname(ikap)))/muf(idepth)
             dkapfstd(idepth) = dkapfstd(idepth) + kapstdmol(idepth,1,ikap)*Amolf(idepth,imol(kapname(ikap)))*molmass(imol(kapname(ikap)))/muf(idepth)
          enddo
       enddo

       kapnumol = 0.
       kapstdmol = 0.

       call cubint(ntktab,npktab,nfreqkap,nkap,logTktab,logPktab,dktabdt,ndepth+1,logTatm,logPatm,kapnumol)
       call cubint(ntktab,npktab,1,nkap,logTktab,logPktab,dkstddt,ndepth+1,logTatm,logPatm,kapstdmol)

       ! Set floor on opacity so it cannot be negative
       where (kapnumol .lt. 0.0)
         kapnumol = 0.0
       endwhere
       where (kapstdmol .lt. 0.0)
        kapstdmol = 0.0
       endwhere

       do ikap=1,nkap
          do idepth = 1,ndepth+1
             dkapfnu(ndepth+1+idepth,:) = dkapfnu(ndepth+1+idepth,:) + kapnumol(idepth,:,ikap)*Amolf(idepth,imol(kapname(ikap)))*molmass(imol(kapname(ikap)))/muf(idepth)
             dkapfstd(ndepth+1+idepth) = dkapfstd(ndepth+1+idepth) + kapstdmol(idepth,1,ikap)*Amolf(idepth,imol(kapname(ikap)))*molmass(imol(kapname(ikap)))/muf(idepth)
          enddo
       enddo
       
       deallocate(kapnumol)
       deallocate(logTatm,logPatm)
       deallocate(kapstdmol)

    else

       !uses only constant opacities
       kapnu = 0.
       do ikap=1,nkap
          do idepth = 1,ndepth
             kapnu(idepth,:) = kapnu(idepth,:) + ktab(1,1,1,ikap)*Amol(idepth,imol(kapname(ikap)))*molmass(imol(kapname(ikap)))/mu(idepth)

             kapstd(idepth) = kapstd(idepth) + kstdtab(1,1,1,ikap)*Amol(idepth,imol(kapname(ikap)))*molmass(imol(kapname(ikap)))/mu(idepth)

          enddo

          do idepth = 1,ndepth+1
             kapfnu(idepth,:) = kapfnu(idepth,:) + ktab(1,1,1,ikap)*Amolf(idepth,imol(kapname(ikap)))*molmass(imol(kapname(ikap)))/muf(idepth)

             kapfstd(idepth) = kapfstd(idepth) + kstdtab(1,1,1,ikap)*Amolf(idepth,imol(kapname(ikap)))*molmass(imol(kapname(ikap)))/muf(idepth)

          enddo

       enddo

    
       
    endif
    
    kapnu = max(kapnu,0.)
    kapfnu = max(kapfnu,0.)

    do idepth=1,ndepth+1
       
       if(kapfstd(idepth) .lt. 0.) then
          
            if (mype==cputerm) then          
               write(*,headrfmt) 'Error get_lbl_opacity'
               write(*,'(A)') 'kapfstd .lt. 0.'
               write(*,intpmfmt)'idepth',idepth
               write(*,exppmfmt)'kapfstd',kapfstd(idepth)
               write(*,*) 
               write(*,separfmt)
               stop
            endif

         endif

         if (idepth .le. ndepth) then
            
            if(kapstd(idepth) .lt. 0.) then
          
               if (mype==cputerm) then          
                  write(*,headrfmt) 'Error get_lbl_opacity'
                  write(*,'(A)') 'kapstd .lt. 0.'
                  write(*,intpmfmt)'idepth',idepth
                  write(*,exppmfmt)'kapstd',kapstd(idepth)
                  write(*,*) 
                  write(*,separfmt)
                  stop
               endif
               
            endif
            
         endif
         
         
         
      enddo

        !!!! Artificial Cloud deck !!!!

    if (cloud) then

       do idepth = 1,ndepth
          if(pp(idepth)>ppmin_cloud .and. pp(idepth)<ppmax_cloud) then
             kapnu(idepth,:) = kapnu(idepth,:) + ((5.31e-27)/(unitL**2))*1/(2*amu)*cloud_strength
          endif
       enddo

       do idepth = 1,ndepth+1
          if(ppf(idepth)>ppmin_cloud .and. ppf(idepth)<ppmax_cloud) then
             kapfnu(idepth,:) = kapfnu(idepth,:)  + ((5.31e-27)/(unitL**2))*1/(2*amu)*cloud_strength
          endif
       enddo

    endif

    ! 8.48e-27 is the scattering crossection of molecular hydrogen at around 350 nm (band 286 in 500 band run)
    ! this can be checked using "print *,"jm = ", kr(286*nkmix+1,1)" which gives 8.48e-43 in code units
    ! sing 2016 nature has value of 5.31e-27

    ! Amol that is hydrogen abundance for cloud taken as 1 and mean molecular weight with just molecular hydrogen as 2

!!!! artificial cloud deck ends !!!!
      
      if (kap_smooth) call smooth_opacity
      
!!! Calculation of Optical Depth in LBL Mode

  do idepth = 2,ndepth+1
   cfdtaufstd(idepth-1) = (ppf(idepth)-ppf(idepth-1))*kapstd(idepth-1)/gr(idepth-1) !! Difference in taufstd in hydrostatic equilibrium
  end do

  do idepth = 1,ndepth
    do ifreq = 1,nfreqkap

     cfdtaufnu(idepth,ifreq) = (cfdtaufstd(idepth)*kapnu(idepth,ifreq))/kapstd(idepth) !! Calculating optical depth at each levels 

    end do

  end do

do idepth = 1,ndepth+1
   do ifreq = 1,nfreqkap
    
    if (idepth == 1) then
    taub(idepth,ifreq) = 0
    else
    taub(idepth,ifreq) = sum(cfdtaufnu(1:(idepth-1),ifreq), 1)   !! Summing optical depth from top of the atmosphere upto each level for each frequency
    end if
    
   end do
end do


  end subroutine get_lbl_opacity

   !> Fill the opacity profile of the atmosphere in function of P,T with correlated k opacities
  !<
  subroutine get_corr_k_opacity

    use mod_grid
    use mod_util
    use mod_chem
    use mod_cst
    use mod_param
    use mod_eos
    
    implicit none
    include 'netcdf.inc'

    integer :: ikap,idepth,ip,it,ifreq,ier,ib,m,ik,iloc,ikapmin,ikap_start,ikap_end,bin,ik2,ikbin,imax,i
       
    real,dimension(:),allocatable :: logTatm,logPatm
    real,dimension(nPktab) :: logPktab
    real,dimension(nTktab) :: logTktab
    real,dimension(:,:,:),allocatable  :: kapnumol,kapstdmol
    real,dimension(ntktab,npktab,1,nkap) :: dkstddp,dkstddt
    real,dimension(ndepth)   :: mu, cfdtaufstd
    real,dimension(ndepth+1) :: muf
    real,dimension(:,:),allocatable :: kapmix,weightsmix
    real,dimension(:,:),allocatable :: kapmix_prod,weightsmix_prod
    real :: min_weights,min_kap,kaploc,dlkap,tauloc,wloc,kapstd_min
    ! real,dimension(ndepth+1,nband) :: taub
    real,dimension(ndepth,(nband*nkmix)) :: cfdtaufnu, taub1
    real,dimension(ndepth+1,nband) :: taub3, taub2
    
    kapnu   = 0.
    kapfnu  = 0.
    weightsnu = 0.
    weightsfnu = 0.
    dkapnu  = 0.
    dkapfnu = 0.
    kapstd = 0.
    kapfstd = 0.
    dkapstd = 0.
    dkapfstd = 0.

    mu  = get_mu(ndepth,Amol)
    muf = get_mu(ndepth+1,Amolf)

     
    logPktab = log10(Pktab)
    logTktab = log10(Tktab)
    
    allocate(logTatm(ndepth),logPatm(ndepth))
    allocate(kapnumol(ndepth,nband*nk,nkap))
    
!     logTatm = log10(tt)
!     logPatm = log10(pp)
! 
!     where (tt<tktabmin) logTatm = log10(tktabmin)

    where (pp<pktabmin)
       logPatm = log10(pktabmin)
    elsewhere
       logPatm = log10(pp)
    endwhere
    
    
    where (tt<tktabmin)
        logTatm = log10(tktabmin)
    elsewhere (tt>tktabmax)
        logTatm = log10(tktabmax)
    elsewhere
        logTatm = log10(tt)
    endwhere 

   ! ktab = max(ktab,1E-100
   
    !cubic spline interpolation
    call cubint(ntktab,npktab,nband*nk,nkap,logTktab,logPktab,ktab,ndepth,logTatm,logPatm,kapnumol)
    
    kapnumol = max(kapnumol,0.)

    !write(*,*) kapnumol(1:10,1,3)
    !stop

    

    do ib =1,nband
       
       m = 1
       
       allocate(kapmix(ndepth,m))
       allocate(weightsmix(ndepth,m))

       weightsmix = 1.
       kapmix = 0.
       
       do ikap = 1,nkap
          
          if (nkmax(ib,ikap)/=0) then
             allocate(kapmix_prod(ndepth,m*nkmax(ib,ikap)))
             allocate(weightsmix_prod(ndepth,m*nkmax(ib,ikap)))
          
             do ik = 1,nkmax(ib,ikap)
                
                iloc = nk*(ib-1)+ik
                
                do ik2 = 1,m


                   kapmix_prod(:,(ik-1)*m+ik2) = kapmix(:,ik2) + Amol(:,imol(kapname(ikap)))*molmass(imol(kapname(ikap)))/mu(:)*kapnumol(:,iloc,ikap)
                   
                   
                   weightsmix_prod(:,(ik-1)*m+ik2) = weightsmix(:,ik2)*weights(iloc,ikap)


                enddo
                
             enddo

             
             deallocate(kapmix,weightsmix)
             m = m*nkmax(ib,ikap)

             if (m.le.nkmix) then

                allocate(kapmix(ndepth,m))
                allocate(weightsmix(ndepth,m))

                kapmix = kapmix_prod
                weightsmix = weightsmix_prod
                
                deallocate(kapmix_prod,weightsmix_prod)

             else

                allocate(kapmix(ndepth,nkmix))
                allocate(weightsmix(ndepth,nkmix))
                kapmix = 0.
                weightsmix = 0.


                do idepth = 1,ndepth
                   
                   do ik = 1,m
                            
                      min_kap = kapmix_prod(idepth,ik)
                      min_weights = weightsmix_prod(idepth,ik)
                      ikapmin = ik
                            
                      do ik2 = ik,m
                               
                         if(kapmix_prod(idepth,ik2)<min_kap) then
                            min_kap = kapmix_prod(idepth,ik2)
                            min_weights = weightsmix_prod(idepth,ik2)
                            ikapmin = ik2
                         endif
                
                      enddo
                            
                      kapmix_prod(idepth,ikapmin) = kapmix_prod(idepth,ik)
                      weightsmix_prod(idepth,ikapmin) = weightsmix_prod(idepth,ik)
                      
                      kapmix_prod(idepth,ik) = min_kap
                      weightsmix_prod(idepth,ik) = min_weights
                
                   enddo
                enddo

                !do idepth=1,ndepth
                idepth=itauref
                               
                bin = m/nkmix
                         
                ikap_start = 0
                ikap_end = 0
                
                do ik=1,nkmix
                   
                   if (ik .le. m-bin*nkmix) then
                      
                      ikap_start = ikap_end +1
                      ikap_end   = ikap_end+ bin+1
                      
                   else
                      
                      ikap_start = ikap_end +1
                      ikap_end   = ikap_end+ bin
                      
                   endif
                   
                   if (ik==nkmix) ikap_end = m
                   
                   ikbin = ikap_end-ikap_start+1
                   
                   weightsmix(idepth,ik) = sum(weightsmix_prod(idepth,ikap_start:ikap_end))

                   !if (weightsmix(idepth,ik)>0.) then
                      kapmix(idepth,ik) = sum(kapmix_prod(idepth,ikap_start:ikap_end)*weightsmix_prod(idepth,ikap_start:ikap_end))/weightsmix(idepth,ik)
                   !else
                   !   kapmix(idepth,ik) = 0.
                   !endif
                   
                enddo

                      
                     ! enddo

                         do idepth=1,ndepth
                         
                      if (idepth /= itauref) then

                         weightsmix(idepth,:) = weightsmix(itauref,:)
                     
                         ik2 = 1
                         do ik = 1,nkmix
                            
                            wloc = 0.
                            
                      
                            do while(wloc+weightsmix_prod(idepth,ik2).lt. weightsmix(idepth,ik))
                              
                               kapmix(idepth,ik) = kapmix(idepth,ik)+kapmix_prod(idepth,ik2)*weightsmix_prod(idepth,ik2)
                      
                               if (ik2+1 .le. m) then
                               wloc = wloc + weightsmix_prod(idepth,ik2)
                               ik2 = ik2+1
                               else

                                  wloc =  weightsmix(idepth,ik) 
                                  
                               endif

                            enddo
                         
                            kapmix(idepth,ik) = kapmix(idepth,ik)+kapmix_prod(idepth,ik2)*(weightsmix(idepth,ik)-wloc)
                            weightsmix_prod(idepth,ik2) = weightsmix_prod(idepth,ik2)-(weightsmix(idepth,ik)-wloc)
                            
                            !if (weightsmix(idepth,ik)>0.) then
                               kapmix(idepth,ik) = kapmix(idepth,ik)/weightsmix(idepth,ik)
                           ! else
                           !    kapmix(idepth,ik) = 0.
                           ! endif
                               
                         enddo
                         
                      endif
                     
                      enddo

                   deallocate(kapmix_prod,weightsmix_prod)
                   m = nkmix

                endif
             
             

          endif
       
          
       enddo

       kapnu(:,nkmix*(ib-1)+1:nkmix*(ib-1)+m) = kapmix(:,1:m)
       weightsnu(:,nkmix*(ib-1)+1:nkmix*(ib-1)+m) = weightsmix(:,1:m)

     
       
       deallocate(kapmix,weightsmix)
    enddo
       
    deallocate(kapnumol)
    deallocate(logTatm,logPatm)
   
     
    !compute kapfnu at cell interfaces
    
    allocate(logTatm(ndepth+1),logPatm(ndepth+1))
    allocate(kapnumol(ndepth+1,nband*nk,nkap))

!     logTatm = log10(ttf)  
!     logPatm = log10(ppf)

!     where (ttf<tktabmin) logTatm = log10(tktabmin)

    where (ppf<pktabmin)
       logPatm = log10(pktabmin)
    elsewhere
       logPatm = log10(ppf)
    endwhere
    
    where (ttf<tktabmin)
        logTatm = log10(tktabmin)
    elsewhere (ttf>tktabmax)
       logTatm = log10(tktabmax)
    elsewhere
       logTatm = log10(ttf)
    endwhere 



    kapnumol = 0.

!    ktab = max(ktab,1E-100)
    
    !cubic spline interpolation
    call cubint(ntktab,npktab,nband*nk,nkap,logTktab,logPktab,ktab,ndepth+1,logTatm,logPatm,kapnumol)
   
    kapnumol = max(kapnumol,0.)
    
  !  write(*,*) kapnumol(1:15,(nband_std-1)*nk+1,3)
  !  stop

  do ib =1,nband
       
       m = 1
       
       allocate(kapmix(ndepth+1,m))
       allocate(weightsmix(ndepth+1,m))

       weightsmix = 1.
       kapmix = 0.
      
       do ikap = 1,nkap
     
          if (nkmax(ib,ikap)/=0) then
             allocate(kapmix_prod(ndepth+1,m*nkmax(ib,ikap)))
             allocate(weightsmix_prod(ndepth+1,m*nkmax(ib,ikap)))
          
             do ik = 1,nkmax(ib,ikap)
                
                iloc = nk*(ib-1)+ik
                
                do ik2 = 1,m

                   kapmix_prod(:,(ik-1)*m+ik2) = kapmix(:,ik2) + Amolf(:,imol(kapname(ikap)))*molmass(imol(kapname(ikap)))/muf(:)*kapnumol(:,iloc,ikap)

                  
                   
                   weightsmix_prod(:,(ik-1)*m+ik2) = weightsmix(:,ik2)*weights(iloc,ikap)


                   
                enddo

             enddo

             deallocate(kapmix,weightsmix)
             m = m*nkmax(ib,ikap)

             if (m.le.nkmix) then

                allocate(kapmix(ndepth+1,m))
                allocate(weightsmix(ndepth+1,m))

                kapmix = kapmix_prod
                weightsmix = weightsmix_prod
                
                deallocate(kapmix_prod,weightsmix_prod)

             else

                allocate(kapmix(ndepth+1,nkmix))
                allocate(weightsmix(ndepth+1,nkmix))
               
                kapmix = 0.
                weightsmix = 0.

                 
                do idepth = 1,ndepth+1
                  
                  
                   do ik = 1,m
                
                      min_kap = kapmix_prod(idepth,ik)
                      min_weights = weightsmix_prod(idepth,ik)
                      ikapmin = ik
                      
                      do ik2 = ik,m
                         
                         if(kapmix_prod(idepth,ik2)<min_kap) then
                            min_kap = kapmix_prod(idepth,ik2)
                            min_weights = weightsmix_prod(idepth,ik2)
                            ikapmin = ik2
                         endif
                
                      enddo

                      kapmix_prod(idepth,ikapmin) = kapmix_prod(idepth,ik)
                      weightsmix_prod(idepth,ikapmin) = weightsmix_prod(idepth,ik)

                      kapmix_prod(idepth,ik) = min_kap
                      weightsmix_prod(idepth,ik) = min_weights
                
                   enddo

                enddo
              
!                do idepth = 1,ndepth+1
                   idepth = itauref
                   bin = m/nkmix
                   ikap_start = 0
                   ikap_end = 0
                   
                   do ik=1,nkmix

                      if (ik .le. m-bin*nkmix) then
                         
                         ikap_start = ikap_end + 1
                         ikap_end = ikap_end + (bin+1)
                    

                      else

                          ikap_start = ikap_end + 1
                         ikap_end = ikap_end + bin

                                          
                      endif

                      if (ik==nkmix) ikap_end = m
                
                      ikbin = ikap_end-ikap_start+1
               
                      weightsmix(idepth,ik) = sum(weightsmix_prod(idepth,ikap_start:ikap_end))

                    !  if (weightsmix(idepth,ik)>0.) then
                         kapmix(idepth,ik) = sum(kapmix_prod(idepth,ikap_start:ikap_end)*weightsmix_prod(idepth,ikap_start:ikap_end))/weightsmix(idepth,ik)
                    !  else
                    !     kapmix(idepth,ik) = 0.
                    !  endif
                      
                      

                   enddo
                   
                  ! enddo

                   do idepth=1,ndepth+1

                      if (idepth /= itauref) then

                         weightsmix(idepth,:) = weightsmix(itauref,:)
                         ik2 = 1
                         do ik = 1,nkmix
                            
                            wloc = 0.
                            
                            
                            do while(wloc+weightsmix_prod(idepth,ik2).lt. weightsmix(idepth,ik))
                               kapmix(idepth,ik) = kapmix(idepth,ik)+kapmix_prod(idepth,ik2)*weightsmix_prod(idepth,ik2)

                                if (ik2+1 .le. m) then
                               wloc = wloc + weightsmix_prod(idepth,ik2)
                               ik2 = ik2+1
                               else

                                  wloc =  weightsmix(idepth,ik) 
                                  
                               endif


                            enddo

                            kapmix(idepth,ik) = kapmix(idepth,ik)+kapmix_prod(idepth,ik2)*(weightsmix(idepth,ik)-wloc)
                            weightsmix_prod(idepth,ik2) = weightsmix_prod(idepth,ik2)-(weightsmix(idepth,ik)-wloc)
                         !   if (weightsmix(idepth,ik)>0.) then
                               kapmix(idepth,ik) = kapmix(idepth,ik)/weightsmix(idepth,ik)
                         !   else
                         !      kapmix(idepth,ik) = 0.
                         !   endif
                            
                         enddo
                         endif

                   enddo

                deallocate(kapmix_prod,weightsmix_prod)
                m = nkmix
             
             endif
         
          endif
          
       enddo
       
       
       kapfnu(:,nkmix*(ib-1)+1:nkmix*(ib-1)+m) = kapmix(:,1:m)
       weightsfnu(:,nkmix*(ib-1)+1:nkmix*(ib-1)+m) = weightsmix(:,1:m)
             

       deallocate(kapmix,weightsmix)

    enddo

    ! weightsnu(1,:) = sum(weightsnu(:,:),1)/ndepth

    ! do idepth = 2,ndepth

    !    weightsnu(idepth,:) = weightsnu(1,:)

    ! enddo

    ! weightsfnu(1,:) = sum(weightsfnu(:,:),1)/(ndepth+1)

    ! do idepth = 1,ndepth+1

    !    weightsfnu(idepth,:) = weightsfnu(1,:)

    ! enddo

       
    deallocate(kapnumol)
    deallocate(logTatm,logPatm)

  
    

    ! do ib=1,nband

    ! do idepth=1,ndepth+1
       
    !    if(ppf(idepth)<pktabmin) then
          
    !       kapstd_min = sum(kapfnu(idepth+1,(ib-1)*nkmix+1:ib*nkmix)*weightsnu(idepth+1,(ib-1)*nkmix+1:ib*nkmix))

    !    endif
    
    ! enddo

    ! do idepth=1,ndepth+1

    !    if(ppf(idepth)<pktabmin) then

    !       if (ib==nband_std) then
    !          kapfnu(idepth,(ib-1)*nkmix+1:ib*nkmix) = kapstd_min
    !       else
    !          kapfnu(idepth,(ib-1)*nkmix+1:ib*nkmix) = 0.

    !       endif
    !    endif

    ! enddo

    ! do idepth=1,ndepth

    !    if(pp(idepth)<pktabmin) then
          
          
    !       if (ib==nband_std) then
    !          kapnu(idepth,(ib-1)*nkmix+1:ib*nkmix) = kapstd_min
    !       else
    !          kapnu(idepth,(ib-1)*nkmix+1:ib*nkmix) = 0.
             
    !       endif
          


    !    endif

    ! enddo

 !enddo

    imax = 0

    do idepth=1,ndepth+1

       if(ttf(idepth)>tktabmax .and. ppf(idepth)*unitP<1.) then

          imax = idepth

       endif

    enddo

    if (kap_smooth) call smooth_opacity

   kapstd(:) = sum(kapnu(:,(nband_std-1)*nkmix+1:nband_std*nkmix)*weightsnu(:,(nband_std-1)*nkmix+1:nband_std*nkmix),2)
   kapfstd(:) = sum(kapfnu(:,(nband_std-1)*nkmix+1:nband_std*nkmix)*weightsfnu(:,(nband_std-1)*nkmix+1:nband_std*nkmix),2)

     !!!! Artificial Cloud deck !!!!


    ! 8.48e-27 is the scattering crossection of molecular hydrogen at around 350 nm (band 286 in 500 band run)
    ! this can be checked using "print *,"jm = ", kr(286*nkmix+1,1)" which gives 8.48e-43 in code units
    ! sing 2016 nature has value of 5.31e-27

    ! Amol that is hydrogen abundance for cloud taken as 1 and mean molecular weight with just molecular hydrogen as 2

!!!! artificial cloud deck ends !!!!
    if (cloud) then

       do idepth = 1,ndepth
          if(pp(idepth)>ppmin_cloud .and. pp(idepth)<ppmax_cloud) then
             kapnu(idepth,:) = kapnu(idepth,:)+ ((8.48e-27)/(unitL**2))*1/(2*amu)*cloud_strength
          endif
       enddo

       do idepth = 1,ndepth+1
          if(ppf(idepth)>ppmin_cloud .and. ppf(idepth)<ppmax_cloud) then
             kapfnu(idepth,:) =  kapfnu(idepth,:)+ ((8.48e-27)/(unitL**2))*1/(2*amu)*cloud_strength
          endif
       enddo

    endif

     if (ocean) then

       do idepth = 1,ndepth
          if(pp(idepth)>ppmin_ocean .and. pp(idepth)<ppmax_ocean) then
             kapnu(idepth,:) =  opacity_ocean
          endif
       enddo

       do idepth = 1,ndepth+1
          if(ppf(idepth)>ppmin_ocean .and. ppf(idepth)<ppmax_ocean) then
             kapfnu(idepth,:) =  opacity_ocean
          endif
       enddo

    endif

    if (kap_smooth) call smooth_opacity

 !   do ib = 1,nband
 !   ib = 8
    
 !   do idepth = 1,ndepth+1
 !   idepth=1

 !   write(*,*)  ib,idepth,ppf(idepth)*unitP/1E6,sum(kapfnu(idepth,(ib-1)*nkmix+1:ib*nkmix)*weightsfnu(idepth,(ib-1)*nkmix+1:ib*nkmix))
 !   write(*,*) idepth,ttf(idepth),amolf(idepth,imol('Na')),amolf(idepth,imol('K'))
 
 !   enddo
 !   enddo
   
    do idepth=1,ndepth+1

       if(kapfstd(idepth) .lt. 0.) then
          
            if (mype==cputerm) then          
               write(*,headrfmt) 'Error get_corr_k_opacity'
               write(*,'(A)') 'kapfstd .lt. 0.'
               write(*,intpmfmt)'idepth',idepth
               write(*,exppmfmt)'kapfstd',kapfstd(idepth)
               write(*,*) 
               write(*,separfmt)
               stop
            endif

         endif

         if (idepth .le. ndepth) then
            
            if(kapstd(idepth) .lt. 0.) then
          
               if (mype==cputerm) then          
                  write(*,headrfmt) 'Error get_corr_k_opacity'
                  write(*,'(A)') 'kapstd .lt. 0.'
                  write(*,intpmfmt)'idepth',idepth
                  write(*,exppmfmt)'kapstd',kapstd(idepth)
                  write(*,*) 
                  write(*,separfmt)
                  stop
               endif
               
            endif
            
         endif
         
         
         
      enddo
      
 !! Calculation of Optical Depth in Correlated-k Mode

 do idepth = 2,ndepth+1
   cfdtaufstd(idepth-1) = (ppf(idepth)-ppf(idepth-1))*kapstd(idepth-1)/gr(idepth-1)
 end do

  do idepth = 1,ndepth
    do ifreq = 1,(nband*nkmix)

     cfdtaufnu(idepth,ifreq) = (cfdtaufstd(idepth)*(kapnu(idepth,ifreq)))/kapstd(idepth) !! Optical depth at each level for each g-space
     !taub1(idepth,ifreq) = exp(-cfdtaufnu(idepth,ifreq))*weightsnu(idepth,ifreq)         
     
    end do

  end do
  
  
   do idepth=1,ndepth+1
     do ib = 1,nband

     if (idepth == 1) then
      taub(idepth,ib) = 0
      !taub2(idepth,ib)= 0
      !taub3(idepth,ib)= 0
     else
     !taub2(idepth,ib) =  sum(taub1(idepth-1,(ib-1)*nkmix+1:ib*nkmix)) ! /(nubandmax(ib)-nubandmin(ib)) !! Summing transmittance in g-space for each band
     !taub3(idepth,ib) = -log(taub2(idepth,ib))  !! Converting back again from transmittance to optical depth 
     !taub(idepth,ib)= sum(taub3(1:idepth,ib))   !! Summing optical depth from top of the atmosphere up-to each level for each frequency.
     
     !! Summing optical depth in vertical then calculating transmittance by taking exponential and multiplying with weights in each g-space,
     !! then summing in g-space and finally taking log again to get optical depth
     taub(idepth,ib) =  -log(sum(exp(-sum(cfdtaufnu(1:idepth-1,(ib-1)*nkmix+1:ib*nkmix),1))*weightsfnu(idepth,(ib-1)*nkmix+1:ib*nkmix)))
      end if
    enddo
    
    
    enddo

    ! Write opacity to a file
    if (write_opacity) then
      call output_corr_k_opacity
    endif


   

      ! taub(1,:) = 0.

      ! do idepth=2,ndepth+1

      !    do ib = 1,nband

      !       taub(idepth,ib) = taub(idepth-1,ib)+sum(kapnu(idepth-1,(ib-1)*nkmix+1:ib*nkmix)*weightsnu(idepth-1,(ib-1)*nkmix+1:ib*nkmix))/kapstd(idepth-1)*(taufstd(idepth)-taufstd(idepth-1))

      !    enddo

      ! enddo
      
      ! do idepth=1,ndepth+1
       
      !    write(*,*)idepth,ppf(idepth)*unitP/1.0E6,ttf(idepth), taub(idepth,19),taub(idepth,17),taub(idepth,14),taub(idepth,11)
         
      ! enddo
      
! # Limits of spectral intervals (wavelengths in m.)
! # Band        Lower limit         Upper limit
! #     1        4.608294931E-05     3.225806452E-04
! #     2        2.000000000E-05     4.608294931E-05
! #     3        1.039501040E-05     2.000000000E-05
! #     4        6.451612903E-06     1.039501040E-05
! #     5        5.219206681E-06     6.451612903E-06
! #     6        4.399472063E-06     5.219206681E-06
! #     7        3.799392097E-06     4.399472063E-06
! #     8        3.288391976E-06     3.799392097E-06
! #     9        2.988643156E-06     3.288391976E-06
! #    10        2.505010020E-06     2.988643156E-06
! #    11        2.170138889E-06     2.505010020E-06
! #    12        2.020202020E-06     2.170138889E-06
! #    13        1.777145904E-06     2.020202020E-06
! #    14        1.593117731E-06     1.777145904E-06
! #    15        1.497005988E-06     1.593117731E-06
! #    16        1.329964091E-06     1.497005988E-06
! #    17        1.197031362E-06     1.329964091E-06
! #    18        1.099989000E-06     1.197031362E-06
! #    19        1.005025126E-06     1.099989000E-06
! #    20        9.599692810E-07     1.005025126E-06
! #    21        9.100009100E-07     9.599692810E-07
! #    22        8.599931201E-07     9.100009100E-07
! #    23        7.849909726E-07     8.599931201E-07
! #    24        7.449899426E-07     7.849909726E-07
! #    25        6.749915626E-07     7.449899426E-07
! #    26        6.119951040E-07     6.749915626E-07
! #    27        5.719842132E-07     6.119951040E-07
! #    28        4.950004950E-07     5.719842132E-07
! #    29        4.000000000E-07     4.950004950E-07
! #    30        3.571428571E-07     4.000000000E-07
! #    31        3.148515475E-07     3.571428571E-07
! #    32        2.000000000E-07     3.148515475E-07


      
    end subroutine get_corr_k_opacity

    ! Subroutine to output the mixed corr k opacity 
    subroutine output_corr_k_opacity

      USE mod_param
      USE mod_util
      USE mod_cst
      USE mod_grid

      IMPLICIT NONE 

      INCLUDE 'netcdf.inc'

      integer id_file,id_lev,id_nband,nvar,ib,idepth,id_int
      integer,allocatable :: id_var(:)

      real,dimension(:,:),allocatable :: kapfnu_out 
      real,dimension(:,:),allocatable :: kapfnumol_out
      real,dimension(:),allocatable :: nu_out

      nvar=6

      allocate(id_var(nvar)) 
      allocate(kapfnu_out(ndepth+1,nband))
      allocate(nu_out(nband))

      id_var=0

      ! Calculate opacity in each model level
      do idepth=1,ndepth+1
        do ib=1,nband
          kapfnu_out(idepth,ib) = sum(kapfnu(idepth,(ib-1)*nkmix+1:ib*nkmix)*weightsfnu(idepth,(ib-1)*nkmix+1:ib*nkmix))
        enddo
      enddo

      do ib=1,nband
        nu_out(ib) = 0.5*(nubandmin(ib)+nubandmax(ib))
      enddo

      ! Create netcdf file 
      call nf(nf_create(fkabsout,nf_clobber,id_file))
      
      ! Create dimensions
      call nf(nf_def_dim(id_file,'nlevel',ndepth+1,id_lev))
      call nf(nf_def_dim(id_file,'nband',nband,id_nband))
      
      ! Create variables 
      call nf(nf_def_var(id_file,'pressure',nf_double,1,id_lev,id_var(1)))
      call nf(nf_def_var(id_file,'temperature',nf_double,1,id_lev,id_var(2)))
      call nf(nf_def_var(id_file,'kapfnu',nf_double,2,(/id_lev,id_nband/),id_var(3)))
      call nf(nf_def_var(id_file,'tau',nf_double,2,(/id_lev,id_nband/),id_var(4)))
      call nf(nf_def_var(id_file,'nu',nf_double,1,id_nband,id_var(5)))

      ! Write to variables 
      call nf(nf_enddef(id_file))
      call nf(nf_put_vara_double(id_file,id_var(1),1,ndepth+1,ppf*unitP))
      call nf(nf_put_vara_double(id_file,id_var(2),1,ndepth+1,ttf*unitK))
      call nf(nf_put_vara_double(id_file,id_var(3),(/1,1/),(/ndepth+1,nband/),kapfnu_out*(unitL**2/unitM)))
      call nf(nf_put_vara_double(id_file,id_var(4),(/1,1/),(/ndepth+1,nband/),taub))
      call nf(nf_put_vara_double(id_file,id_var(5),1,nband,nu_out*(1.0/unitL)))

      ! Close file 
      call nf(nf_close(id_file))



    end subroutine output_corr_k_opacity


    !> Fill the opacity profile of the atmosphere in function of P,T with correlated k opacities
    !<
    subroutine smooth_opacity

      use mod_grid

      implicit none

      real,dimension(ndepth+1,nfreqkap)::kapfnu_smooth
      real,dimension(ndepth,nfreqkap)::kapnu_smooth

      integer :: imin,imax,ifmin,ifmax,idepth,iker

      kerkap_smooth = max(kerkap_smooth ,1)

      imin = kerkap_smooth+1
      imax = ndepth-kerkap_smooth
      
      ifmin = kerkap_smooth+1
      ifmax = ndepth+1-kerkap_smooth

      do idepth = ifmin,ifmax

         if (ppf(idepth)>ppmin_smooth .and. ppf(idepth)<ppmax_smooth) then
    
            kapfnu_smooth(idepth,:) = kapfnu(idepth,:)

            do iker = 1,kerkap_smooth

               kapfnu_smooth(idepth,:) = kapfnu_smooth(idepth,:) +  kapfnu(idepth-iker,:)+ kapfnu(idepth+iker,:)

            enddo

            kapfnu_smooth(idepth,:) = kapfnu_smooth(idepth,:)/(2*kerkap_smooth+1)

         endif

      enddo

      do idepth = ifmin,ifmax
    
         if (ppf(idepth)>ppmin_smooth .and. ppf(idepth)<ppmax_smooth) then
         
            kapfnu(idepth,:) = kapfnu_smooth(idepth,:)
         
         endif
   
      enddo

      do idepth = imin,imax
         
         if (pp(idepth)>ppmin_smooth .and. pp(idepth)<ppmax_smooth) then
    
            kapnu_smooth(idepth,:) = kapnu(idepth,:)

            do iker = 1,kerkap_smooth

               kapnu_smooth(idepth,:) = kapnu_smooth(idepth,:) +  kapnu(idepth-iker,:)+ kapnu(idepth+iker,:)

            enddo

            kapnu_smooth(idepth,:) = kapnu_smooth(idepth,:)/(2*kerkap_smooth+1)

         endif

      enddo

      do idepth = imin,imax
         
         if (pp(idepth)>ppmin_smooth .and. pp(idepth)<ppmax_smooth) then
            
            kapnu(idepth,:) = kapnu_smooth(idepth,:)
            
         endif
         
      enddo
      

    end subroutine smooth_opacity


      !> compute rosseland opacities at cell centre or cell interface
    !<
    function get_kapross(ndepth_loc,tt_loc,kapnu_loc)

      use mod_param
      use mod_grid
      use mod_radtrans

      implicit none

#if PARA==1
      include 'mpif.h'
      integer :: info
#endif

      integer :: idepth,ndepth_loc
      real,dimension(ndepth_loc,nfreqkap) :: dbnudkap,wross,wplanck,kapnu_loc,kapbnu
      real,dimension(ndepth_loc) :: get_kapross,sumwr,sumwrdk,sumwr_loc,sumwrdk_loc,tt_loc,sumwp,sumwpk,sumwp_loc,sumwpk_loc


      if (.not. corr_k) then

      dbnudkap = 0.
      kapbnu = 0.
      wross    = get_dbnu(ndepth_loc,nfreqkap,tt_loc)     
      wplanck  = get_bnu(ndepth_loc,nfreqkap,tt_loc)     
      
      where(kapnu_loc /= 0.) dbnudkap = wross/kapnu_loc
      kapbnu = wplanck*kapnu_loc

      sumwp_loc = sum(wplanck,2)
      sumwpk_loc = sum(kapbnu,2)

      sumwr_loc   = sum(wross,2)
      sumwrdk_loc = sum(dbnudkap,2)

#if PARA == 1

      call MPI_ALLREDUCE(sumwr_loc,sumwr,ndepth_loc,MPI_DOUBLE_PRECISION, &
           & MPI_SUM,MPI_COMM_WORLD,info)
      call MPI_ALLREDUCE(sumwrdk_loc,sumwrdk,ndepth_loc,MPI_DOUBLE_PRECISION, &
           & MPI_SUM,MPI_COMM_WORLD,info)
      call MPI_ALLREDUCE(sumwp_loc,sumwp,ndepth_loc,MPI_DOUBLE_PRECISION, &
           & MPI_SUM,MPI_COMM_WORLD,info)
      call MPI_ALLREDUCE(sumwpk_loc,sumwpk,ndepth_loc,MPI_DOUBLE_PRECISION, &
           & MPI_SUM,MPI_COMM_WORLD,info)

#else

      sumwr   = sumwr_loc
      sumwrdk = sumwrdk_loc
      sumwp   = sumwp_loc
      sumwpk = sumwp_loc

#endif
      
      get_kapross = 0.
      where (sumwrdk /=0.) get_kapross = sumwr/sumwrdk
     

      else
         dbnudkap = 0.
         wross    = get_dbnu(ndepth_loc,nfreqkap,tt_loc)   
         where (kapnu_loc /=0.)  dbnudkap= wross/kapnu_loc
         get_kapross = int_freq(ndepth_loc,nfreqkap,nuf,wross)/int_freq(ndepth_loc,nfreqkap,nuf,dbnudkap)
         
      endif

    end function get_kapross


    !> compute the rosseland opacities and derivatives
    !<
    function get_dkapross(ndepth_loc,tt_loc,kapnu_loc,dkapnu_loc)

      use mod_grid     
      use mod_radtrans

      implicit none

#if PARA==1
      include 'mpif.h'
      integer :: info
#endif

      integer :: idepth,ndepth_loc
      real,dimension(ndepth_loc) :: sum_w,sum_wok,sum_w_loc,sum_wok_loc,tt_loc,sum_dw,sum_dw_loc,sum_dwok,sum_dwok_loc
      real,dimension(2*ndepth_loc,nfreqkap) :: wdkok2,dkapnu_loc
      real,dimension(ndepth_loc,nfreqkap)   :: kapnu_loc,wok,dwok,w,dw
      real,dimension(2*ndepth_loc) :: get_dkapross,sum_wdkok2,sum_wdkok2_loc
      

      !kapross = sum(w)/sum(w/kapnu)
      !dkapross = sum(dw)/sum(w/kapnu)-sum(w)/sum(w/kapnu)**2*(sum(dw/kapnu)-sum(w/kapnu**2*dkapnu))
    
      w    = get_dbnu(ndepth_loc,nfreqkap,tt_loc)
      dw   = get_d2bnu(ndepth_loc,nfreqkap,tt_loc)
   
    
      wok    =  w/kapnu
      dwok   =  dw/kapnu
      wdkok2 =  w/kapnu**2*dkapnu

      

      where(kapnu_loc /= 0.) 

         wok  = w/kapnu_loc
         dwok = dw/kapnu_loc
         
         !pressure derivatives
         wdkok2(1:ndepth_loc,:) = w*dkapnu_loc(1:ndepth_loc,:) / kapnu_loc**2
            
         !temperature derivatives
         wdkok2((ndepth_loc+1):(2*ndepth_loc),:) = w &
              & *dkapnu_loc((ndepth_loc+1):(2*ndepth_loc),:) / kapnu_loc**2
      
      endwhere

      

      sum_w_loc      = sum(w     ,2)
      sum_dw_loc     = sum(dw    ,2)
      sum_wok_loc    = sum(wok   ,2)
      sum_dwok_loc   = sum(dwok  ,2)
      sum_wdkok2_loc = sum(wdkok2,2)
      

#if PARA == 1

      call MPI_ALLREDUCE(sum_w_loc,sum_w,ndepth_loc,MPI_DOUBLE_PRECISION, &
           & MPI_SUM,MPI_COMM_WORLD,info)
      call MPI_ALLREDUCE(sum_dw_loc,sum_dw,ndepth_loc,MPI_DOUBLE_PRECISION, &
           & MPI_SUM,MPI_COMM_WORLD,info)
      call MPI_ALLREDUCE(sum_wok_loc,sum_wok,ndepth_loc,MPI_DOUBLE_PRECISION, &
           & MPI_SUM,MPI_COMM_WORLD,info)
      call MPI_ALLREDUCE(sum_dwok_loc,sum_dwok,ndepth_loc,MPI_DOUBLE_PRECISION, &
            & MPI_SUM,MPI_COMM_WORLD,info)
      call MPI_ALLREDUCE(sum_wdkok2_loc,sum_wdkok2,2*ndepth_loc,MPI_DOUBLE_PRECISION, &
            & MPI_SUM,MPI_COMM_WORLD,info)

#else

      sum_w      = sum_w_loc
      sum_dw     = sum_dw_loc
      sum_wok    = sum_wok_loc
      sum_dwok   = sum_dwok_loc
      sum_wdkok2 = sum_wdkok2_loc 

#endif

      !dkapross = sum(dw)/sum(w/kapnu)-sum(w)/sum(w/kapnu)**2*(sum(dw/kapnu)-sum(w/kapnu**2*dkapnu))
      get_dkapross = 0.

      where(sum_wok /=0.)
      
         !pressure derivative
         get_dkapross(1:ndepth_loc) = sum_w/sum_wok**2*sum_wdkok2(1:ndepth_loc)
         
         !temperature derivatives
         get_dkapross(ndepth_loc+1:2*ndepth_loc) = sum_dw/sum_wok-sum_w/sum_wok**2*&
              & (sum_dwok-sum_wdkok2(ndepth_loc+1:2*ndepth_loc))
         
      endwhere
     
    end function get_dkapross


    function int_freq(ndepth_loc,nfreq_loc,nuf_loc,Anu)

      use mod_grid
      use mod_param
      use mod_cst
      
      implicit none

#if PARA==1
      include 'mpif.h'
      integer :: info
#endif


      integer :: ndepth_loc,ifreq,nfreq_loc,ib
      real,dimension(ndepth_loc,nfreq_loc) :: Anu
      real,dimension(ndepth_loc) :: Atot_loc,int_freq
      real,dimension(nfreq_loc+1) :: nuf_loc

      if(.not. corr_k) then

      Atot_loc = 0.
     
      do ifreq = 1,nfreq_loc
         Atot_loc(:) = Atot_loc(:) + Anu(:,ifreq)*(nuf_loc(ifreq+1)-nuf_loc(ifreq))
      enddo

     

#if PARA == 1

      call MPI_ALLREDUCE(atot_loc,int_freq,ndepth_loc,MPI_DOUBLE_PRECISION, &
           & MPI_SUM,MPI_COMM_WORLD,info)
      
#else

      int_freq = Atot_loc

#endif         
      
   else

      Atot_loc = 0.

      if(ndepth_loc == ndepth) then
         Atot_loc = sum(Anu*weightsnu,2)

      else if (ndepth_loc == ndepth+1) then
         Atot_loc = sum(Anu*weightsfnu,2)
         
         ! do ib =1,nband
         !    if (ib ==24) then
         !    write(*,*)ib, sum(Anu(:,(ib-1)*nkmix+1:ib*nkmix)*weightsfnu(:,(ib-1)*nkmix+1:ib*nkmix),2)/Atot_loc(:)
         ! stop   
         ! endif
            
            
        ! enddo
      else if (ndepth_loc==1) then
         ! integrate planck function at surface
         Atot_loc = sum(Anu)/nkmix

      else
         

         write(*,*) 'error int_freq with corr k ndepth_loc should be ndepth or ndepth+1'
         stop

      endif

      int_freq = Atot_loc
         

   endif
     


    end function int_freq


    !> Compute the scattering opacity for radiative transfer
    !<
    subroutine get_scat_opacity

      use mod_grid
      use mod_chem
      use mod_radtrans
      use mod_eos
      use mod_cst

      implicit none

      integer :: idepth,ib,iray
      real,dimension(ndepth+1) :: muf
      real,dimension(ndepth  ) :: mu

      kscatnu  = 0.
      kscatfnu = 0.

      mu  = get_mu(ndepth,Amol)
      muf = get_mu(ndepth+1,Amolf)
      
      if (rayleigh) then
      
       do iray = 1,num_ray_scatt
      
         if (imol(rayname(iray))>0) then
            do idepth = 1,ndepth
               kscatnu(idepth,:) = kscatnu(idepth,:) + kr(:,iray)*Amol(idepth,imol(rayname(iray)))/(mu(idepth)*amu)
               kscatfnu(idepth,:) = kscatfnu(idepth,:) + kr(:,iray)*Amolf(idepth,imol(rayname(iray)))/(muf(idepth)*amu)
            enddo
            kscatfnu(ndepth+1,:) = kscatfnu(ndepth+1,:) + kr(:,iray)*Amolf(ndepth+1,imol(rayname(iray)))/(muf(ndepth+1)*amu)
         endif

       enddo
       
       if (haze_scatt) then
            do idepth = 1,ndepth
               kscatnu(idepth,:) = kscatnu(idepth,:) + kh(:)
               kscatfnu(idepth,:) = kscatfnu(idepth,:) + kh(:)
            enddo
            kscatfnu(ndepth+1,:) = kscatfnu(ndepth+1,:) + kh(:)
       endif
      
      endif   
      
      where (kapnu+kscatnu /=0.) 
         epsnu = kapnu/(kapnu+kscatnu)
         kapnu = kapnu + kscatnu
      endwhere
      
      where (kapfnu+kscatfnu /=0.) 
         epsfnu = kapfnu/(kapfnu+kscatfnu)
         kapfnu = kapfnu + kscatfnu
      endwhere
      
      if (corr_k) then
      
        kapstd(:) = sum(kapnu(:,(nband_std-1)*nkmix+1:nband_std*nkmix)* &
            weightsnu(:,(nband_std-1)*nkmix+1:nband_std*nkmix),2)
        kapfstd(:) = sum(kapfnu(:,(nband_std-1)*nkmix+1:nband_std*nkmix)* &
            weightsfnu(:,(nband_std-1)*nkmix+1:nband_std*nkmix),2)
      
      else
      
        kapstd(:) = kapnu(:,1)
        kapfstd(:) = kapfnu(:,1)
      
      end if
      
 ! !  ib = 16
 !         idepth=1
 !  !  do idepth = 1,ndepth+1
 !    write(*,*)  sum(kapfnu(idepth,(ib-1)*nkmix+1:ib*nkmix)*weightsfnu(idepth,(ib-1)*nkmix+1:ib*nkmix))
 !    enddo
 !   ! write(*,*)
   ! do idepth = 1,ndepth+1
   !    write(*,*) kapfstd(idepth)
   !enddo
  !  stop
    end subroutine get_scat_opacity


          !> compute rosseland opacities at cell centre or cell interface
    !<
    function get_kapplanck(ndepth_loc,tt_loc,kapnu_loc)

      use mod_param
      use mod_grid
      use mod_radtrans
      

      implicit none

#if PARA==1
      include 'mpif.h'
      integer :: info
#endif

      integer :: idepth,ndepth_loc
      real,dimension(ndepth_loc,nfreqkap) :: dbnudkap,wross,wplanck,kapnu_loc,kapbnu
      real,dimension(ndepth_loc) :: get_kapplanck,sumwr,sumwrdk,sumwr_loc,sumwrdk_loc,tt_loc,sumwp,sumwpk,sumwp_loc,sumwpk_loc


      if (.not. corr_k) then

      dbnudkap = 0.
      kapbnu = 0.
      wross    = get_dbnu(ndepth_loc,nfreqkap,tt_loc)     
      wplanck  = get_bnu(ndepth_loc,nfreqkap,tt_loc)     
      
      where(kapnu_loc /= 0.) dbnudkap = wross/kapnu_loc
      kapbnu = wplanck*kapnu_loc

      sumwp_loc = sum(wplanck,2)
      sumwpk_loc = sum(kapbnu,2)

      sumwr_loc   = sum(wross,2)
      sumwrdk_loc = sum(dbnudkap,2)

#if PARA == 1

      call MPI_ALLREDUCE(sumwr_loc,sumwr,ndepth_loc,MPI_DOUBLE_PRECISION, &
           & MPI_SUM,MPI_COMM_WORLD,info)
      call MPI_ALLREDUCE(sumwrdk_loc,sumwrdk,ndepth_loc,MPI_DOUBLE_PRECISION, &
           & MPI_SUM,MPI_COMM_WORLD,info)
      call MPI_ALLREDUCE(sumwp_loc,sumwp,ndepth_loc,MPI_DOUBLE_PRECISION, &
           & MPI_SUM,MPI_COMM_WORLD,info)
      call MPI_ALLREDUCE(sumwpk_loc,sumwpk,ndepth_loc,MPI_DOUBLE_PRECISION, &
           & MPI_SUM,MPI_COMM_WORLD,info)

#else

      sumwr   = sumwr_loc
      sumwrdk = sumwrdk_loc
      sumwp   = sumwp_loc
      sumwpk = sumwpk_loc

#endif
      
      get_kapplanck = sumwpk/sumwp
     

      else
         dbnudkap = 0.
         wplanck    = get_bnu(ndepth_loc,nfreqkap,tt_loc)   
         dbnudkap= wplanck*kapnu_loc
         get_kapplanck = int_freq(ndepth_loc,nfreqkap,nuf,dbnudkap)/int_freq(ndepth_loc,nfreqkap,nuf,wplanck)
         
      endif

    end function get_kapplanck
    
  
      !> Compute the opacity of H- 
    !<
    subroutine get_hm_opacity
    
      use mod_param
      use mod_grid
      use mod_eos
      use mod_cst
      use mod_chem
      
      implicit none
      
      integer :: ifreq, j, idepth, iband, nfreq_loc, k
      real,dimension(:),allocatable :: nu_hm
      real,dimension(ndepth) :: nd, mu, rho_hm
      real,dimension(ndepth+1) :: ndf, muf, rhof_hm
      real,dimension(:),allocatable :: lam_loc, f_var, shmbf, var_test   !! nfreqkap = nband*nkmix
      real,dimension(:,:),allocatable :: kap_hmbf, kapf_hmbf, var_1, shmff, shmfff, kap_hmff, kapf_hmff, var_test2 !! kapf_hmbf is at cell edges
      real,dimension(6) :: coeff_hmbf, A1, B1, C1, D1, E1, F1, A2, B2, C2, D2, E2, F2 
      real, parameter :: lam_thresh = 1.6419   !! This is in microns as equation requires
      
    
    !!! Bound Free H- opacity using equation from John 1988
      
      if (corr_k) then
          allocate(nu_hm(nband))
          do iband = 1,nband
           nu_hm(iband) = 0.5*(nubandmin(iband)+nubandmax(iband))
          enddo
          allocate(lam_loc(nband))
          allocate(f_var(nband))
          allocate(shmbf(nband))
          allocate(kap_hmbf(ndepth,nband))
          allocate(kapf_hmbf(ndepth+1,nband))
          allocate(var_1(6,nband))
          allocate(shmff(ndepth,nband))
          allocate(shmfff(ndepth+1,nband))
          allocate(kap_hmff(ndepth,nband))
          allocate(kapf_hmff(ndepth+1,nband))
          allocate(var_test(nband))
          allocate(var_test2(6,ndepth))
          nfreq_loc = nband
      else
          allocate(nu_hm(nfreqkap))
          nu_hm = nu
          allocate(lam_loc(nfreqkap))
          allocate(f_var(nfreqkap))
          allocate(shmbf(nfreqkap))
          allocate(kap_hmbf(ndepth,nfreqkap))
          allocate(kapf_hmbf(ndepth+1,nfreqkap))
          allocate(var_1(6,nfreqkap))
          allocate(shmff(ndepth,nfreqkap))
          allocate(shmfff(ndepth+1,nfreqkap))
          allocate(kap_hmff(ndepth,nfreqkap))
          allocate(kapf_hmff(ndepth+1,nfreqkap))
          nfreq_loc = nfreqkap
      endif
      
      lam_loc(:) = 0
      f_var(:) = 0
      shmbf(:) = 0
      kap_hmbf(:,:) = 0
      kapf_hmbf(:,:) = 0
      var_1(:,:) = 0
      shmff(:,:) = 0
      shmfff(:,:) = 0
      kap_hmff(:,:) = 0
      kapf_hmff(:,:) = 0
      var_test(:) = 0
      var_test2(:,:) = 0
      
      coeff_hmbf = (/ 152.519, 49.534, -118.858, 92.536, -34.194, 4.982 /)
      lam_loc(:) = 1e4/(nu_hm(:)*(1/unitL))
 
      do ifreq = 1,nfreq_loc
       if (lam_loc(ifreq) <= lam_thresh) then
        do j = 1,6
  
         f_var(ifreq) = f_var(ifreq) + coeff_hmbf(j)*((((1.0/lam_loc(ifreq)) - (1.0/lam_thresh)))**((j-1)/2.0))

        enddo
        shmbf(ifreq) = 1e-18*(lam_loc(ifreq)**3.0)*(((1.0/lam_loc(ifreq)) - (1.0/lam_thresh))**(3.0/2.0))*f_var(ifreq) !! in units of cm**2
       else
    
         f_var(ifreq) = 0.0
         shmbf(ifreq) = 0.0
       endif 
      enddo

      
      !! Converting to cm2/g by first multiplying by H- ions per cm3 and then dividing by density in g/cm3 (Remember code units)
       !kapmix_prod(:,(ik-1)*m+ik2) = kapmix(:,ik2) + Amolf(:,imol(kapname(ikap)))*molmass(imol(kapname(ikap)))/muf(:)*kapnumol(:,iloc,ikap)
      
      mu  = get_mu(ndepth,Amol)
      muf = get_mu(ndepth+1,Amolf)
      nd  = get_nn(ndepth,pp,tt)*(1/(unitL**3)) ! Converting from code units to CGS
      ndf  = get_nn(ndepth+1,ppf,ttf)*(1/(unitL**3))
      
      rho_hm(:) = rho(:)*(unitM/(unitL**3))
      rhof_hm(:) = rhof(:)*(unitM/(unitL**3))
      
      do idepth=1,ndepth
        kap_hmbf(idepth,:) = (shmbf(:)*Amol(idepth,imol('H_'))*nd(idepth)/rho_hm(idepth))*(unitM/(unitL**2)) !/mu(idepth)*amu ! CGS to code units 
        kapf_hmbf(idepth,:) = (shmbf(:)*Amolf(idepth,imol('H_'))*ndf(idepth)/rhof_hm(idepth))*(unitM/(unitL**2))
      enddo
        kapf_hmbf(ndepth+1,:) = (shmbf(:)*Amolf(ndepth+1,imol('H_'))*ndf(ndepth+1)/rhof_hm(ndepth+1))*(unitM/(unitL**2))

!       print *,"jm = ",  lam_loc(:)
!       print *,"jm1 = ", shmbf(:)
!       print *,"jm1 = ", f_var(:)
!       print *,"jm2 = ", coeff_hmbf(:)
!       stop
!       print *,"jm1 = ", size(rho)
        
  !! Free Free Equation starts
  
       A1 = (/ 518.1021, 473.2636, -482.2089, 115.5291, 0.0, 0.0 /)
       B1 = (/ -734.8666, 1443.4137, -737.1616, 169.6374, 0.0, 0.0 /)
       C1 = (/ 1021.1775, -1977.3395, 1096.8827, -245.6490, 0.0, 0.0 /)
       D1 = (/ -479.0721, 922.3575, -521.1341, 114.2430, 0.0, 0.0 /)
       E1 = (/ 93.1373, -178.9275, 101.7963, -21.9972, 0.0, 0.0 /)
       F1 = (/ -6.4285, 12.3600, -7.0571, 1.5097, 0.0, 0.0 /)

       A2 = (/ 0.0, 2483.3460, -3449.8890, 2200.0400, -696.2710, 88.2830 /)
       B2 = (/ 0.0, 285.8270, -1158.3820, 2427.7190, -1841.4000, 444.5170 /)
       C2 = (/ 0.0, -2054.29, 8746.5230, -13651.1050, 8624.9700, -1863.8640 /)
       D2 = (/ 0.0, 2827.7760, -11485.6320, 16755.5240, -10051.5300, 2095.2880 /)
       E2 = (/ 0.0, -1341.5370, 5303.6090, -7510.4940, 4400.0670, -901.7880 /)
       F2 = (/ 0.0, 208.9520, -812.9390, 1132.7380, -655.0200, 132.9850 /)
       
       do ifreq = 1,nfreq_loc
          do k = 1,6
           if ((0.1823 < lam_loc(ifreq)) .and. (lam_loc(ifreq) < 0.3645)) then
           
             var_1(k,ifreq) = (1e-29)*(((lam_loc(ifreq)**2)*A1(k)) + (B1(k)) + (C1(k)/lam_loc(ifreq)) + (D1(k)/(lam_loc(ifreq)**2)) + (E1(k)/(lam_loc(ifreq)**3)) + (F1(k)/(lam_loc(ifreq)**4)))
   
           else if (lam_loc(ifreq) > 0.3645) then

             var_1(k,ifreq) = (1e-29)*(((lam_loc(ifreq)**2)*A2(k)) + (B2(k)) + (C2(k)/lam_loc(ifreq)) + (D2(k)/(lam_loc(ifreq)**2)) + (E2(k)/(lam_loc(ifreq)**3)) + (F2(k)/(lam_loc(ifreq)**4)))
  
           endif
          enddo
       enddo

       do k = 1,6
        var_test(:) = var_test(:) + var_1(k,:)
       enddo

      
      do idepth = 1, ndepth
       do k = 1,6
        var_test2(k,idepth) = ((5040.0/tt(idepth))**((k+1.0)/2.0))
       enddo
      enddo
      
        
      do ifreq = 1, nfreq_loc 
       do idepth = 1,ndepth
        do k = 1,6
         shmff(idepth,ifreq) = shmff(idepth,ifreq) + var_1(k,ifreq)*(((5040.0/tt(idepth))**((k+1.0)/2.0)))
         shmfff(idepth,ifreq) = shmfff(idepth,ifreq) + var_1(k,ifreq)*(((5040.0/ttf(idepth))**((k+1.0)/2.0)))
        enddo
         kap_hmff(idepth,ifreq) = (shmff(idepth,ifreq)*(1.380620e-16)*tt(idepth)*(Amol(idepth,imol('H'))*nd(idepth))*(Amol(idepth,imol('e_'))*nd(idepth))/(rho_hm(idepth)))*(unitM/(unitL**2))
         kapf_hmff(idepth,ifreq) = (shmfff(idepth,ifreq)*(1.380620e-16)*ttf(idepth)*(Amolf(idepth,imol('H'))*ndf(idepth))*(Amolf(idepth,imol('e_'))*ndf(idepth))/(rhof_hm(idepth)))*(unitM/(unitL**2))
       enddo
      enddo
     
      do ifreq = 1, nfreq_loc
       do k =1,6
        shmfff(ndepth+1,ifreq) = shmfff(ndepth+1,ifreq) + var_1(k,ifreq)*(((5040.0/ttf(ndepth+1))**((k+1.0)/2.0)))
       enddo
       kapf_hmff(ndepth+1,ifreq) = (shmfff(ndepth+1,ifreq)*(1.380620e-16)*ttf(ndepth+1)*(Amolf(ndepth+1,imol('H'))*ndf(ndepth+1))*(Amolf(ndepth+1,imol('e_'))*ndf(ndepth+1))/(rhof_hm(ndepth+1)))*(unitM/(unitL**2))
      enddo
   
   !!! Adding to ATMO opacity starts     
      if (corr_k) then
       do idepth = 1,ndepth
        do ifreq = 1,nfreq_loc
         kapnu(idepth,nkmix*(ifreq-1)+1:ifreq*nkmix) = kapnu(idepth,nkmix*(ifreq-1)+1:ifreq*nkmix)  + kap_hmbf(idepth,ifreq) + kap_hmff(idepth,ifreq)
         kapfnu(idepth,nkmix*(ifreq-1)+1:ifreq*nkmix) = kapfnu(idepth,nkmix*(ifreq-1)+1:ifreq*nkmix)  + kapf_hmbf(idepth,ifreq) + kapf_hmff(idepth,ifreq)
        enddo
       enddo
       do ifreq = 1,nfreq_loc
        kapfnu(ndepth+1,nkmix*(ifreq-1)+1:ifreq*nkmix) = kapfnu(ndepth+1,nkmix*(ifreq-1)+1:ifreq*nkmix)  + kapf_hmbf(ndepth+1,ifreq) + kapf_hmff(ndepth+1,ifreq)
       enddo
      
      else
      
       do idepth = 1,ndepth
        do ifreq = 1,nfreq_loc
         kapnu(idepth,ifreq) = kapnu(idepth,ifreq)  + kap_hmbf(idepth,ifreq) + kap_hmff(idepth,ifreq)
         kapfnu(idepth,ifreq) = kapfnu(idepth,ifreq)  + kapf_hmbf(idepth,ifreq) + kapf_hmff(idepth,ifreq)
        enddo
       enddo
       do ifreq = 1,nfreq_loc
        kapfnu(ndepth+1,ifreq) = kapfnu(ndepth+1,ifreq)  + kapf_hmbf(ndepth+1,ifreq) + kapf_hmff(ndepth+1,ifreq)
       enddo
      endif
      
      
!       open (unit = 12, file = "total_h-_opacity.txt")    
!       write (12,*),  ((kap_hmbf(1,:) + kap_hmff(1,:))*((unitL**2)/unitM)), lam_loc(:)
!       close(12)
!       stop
!         print *,"j  = ",  kapnu(10,:)
!         print *,"jm = ",  kap_hmbf(10,:)
!         print *,"jm1 = ", kap_hmff(10,:)

      
  end subroutine get_hm_opacity

end module mod_opacity
