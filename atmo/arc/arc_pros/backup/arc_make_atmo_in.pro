pro arc_make_atmo_in,teff,MdH,COratio,Rp,Cloud,haze,fitabundances,cloudtop,kzz,fast=fast,contribution=contribution
  ; Generate ATMO input file
  ; ISOTHRMAL SINGLE Metalicity
  Common SystemParameters,rplanet,logg,rstar,rorbit,stellarfile,Tstar,Cloudtype,Hazetype,Chemistry,Isothermal,Opacity_switch,Abund_fit_switch,Abund_fit_index,Abund_priorLow,Abund_priorUP,rain_switch,Element_fit_switch,Element_fit_index,Ele_priorLow,Ele_priorUP
  Common cpcount,ChildProcessCount,dirs,allmodel,allmodelpoints,Emission,nlevel,pmin,pmax

 MdHin=alog10(MdH)
 ;ATMO Crashes for MdH > 4.50, Set an upper limit on Metalicity of 10^4.5=31622x solar
  if (MdHin gt 4.50) then MdHin = 4.50
 ;Doesn't handle negitive clouds
 ; if (Cloud le 0) then Cloud = 0

  ;Test for ATMO Namelist file and spectra, delete if present
  Result = FILE_TEST('ARC'+string(ChildProcessCount,format='(I02)')+'.in')  
  if (Result eq 1 ) then spawn,'rm ARC'+string(ChildProcessCount,format='(I02)')+'.in'

  ;Open Namelist file to writeto
  close,1 & openw,1,'ARC'+string(ChildProcessCount,format='(I02)')+'.in'

  ;PARAM 
  printf,1,'&PARAM'
  printf,1,'Debug  = 0'
  printf,1,"fout = 'pt_ARC"+string(ChildProcessCount,format='(I02)')+".out.ncdf'"
  if (Isothermal eq 1) then printf,1,"fin = 'None'"
  if (Isothermal eq 0) then printf,1,"fin = 'pt_ARC"+string(ChildProcessCount,format='(I02)')+".ncdf'"
  printf,1,'/'
  printf,1,''
  printf,1,'&EOS gamma = 0./'
  printf,1,''

  ;GRID 
  printf,1,'&GRID'
  printf,1,'logg   = '+string(logg,format='(F4.2)')
  printf,1,'teff   = '+string(teff,format='(F14.9)')
  if keyword_set(fast) then begin  ;low-resolution
    printf,1,'ndepth = '+string(nlevel,format='(I2)')
    printf,1,'pmin = '+string(pmin,format='(E7.1)')
    printf,1,'pmax = '+string(pmax,format='(E7.1)')
  endif else begin
    printf,1,'ndepth = 50' ;high-resolution, more depth points
    printf,1,'pmin = '+string(pmin,format='(E7.1)')
    printf,1,'pmax = '+string(pmax,format='(E7.1)')
  endelse
  printf,1,'Rp = '+string(Rp,format='(F15.9)')
  if (Emission eq 'no') then printf,1,'pp_Rp = 0.001' ;+string(pp_Rp,format='(F16.12)')  ;pressure in bar of Rp
  printf,1,''
  printf,1,'taumin = 1E-6'
  printf,1,'taumax = 8E3'
  printf,1,''
  if keyword_set(fast) then begin
    if (Emission eq 'no')  then printf,1,'nfreq  = 350'
    if (Emission eq 'yes')  then printf,1,'nfreq  = 150'
    printf,1,'corr_k = t'
    printf,1,'nkmix  = 19'     
    if (Emission eq 'no')  then printf,1,'nband  = 350'   
    if (Emission eq 'yes')  then printf,1,'nband  = 150'   
  endif else begin
    printf,1,'nfreq  = 3500'
    printf,1,'corr_k = t'
    printf,1,'nkmix  = 17'
    printf,1,'nband  = 3500'
  endelse
  printf,1,'nband_std = 32'   ;32
  printf,1,'/'

  ;CHEMISTRY chem.f90
  ;NAMELIST /chemistry/ chem, fAin, fAeqout, fAneqout, fcoeff, fcoeffnine, MdH, Tcrit, &
  ;                     Pcrit, tfreeze_eq, Acst, print_chem, cond_NH3, cond_H2O, chem_conv, rainout, &
  ;                     calc_cp_chem, cp_chem_reac, calc_hform, print_thermodynamics, init_eq, chem_smooth, kerchem_smooth, &
  ;                     elename_freeze, change_chem_smooth, nit_chem_smooth, rainout_min_ele, ele_xfactor
  printf,1,'&CHEMISTRY' 
  if (Chemistry eq 0 or Chemistry eq 3) then printf,1,"chem  = 'eq'" 
  if (Chemistry eq 1) then printf,1,"chem  = 'man'" 
  if (Chemistry eq 2) then printf,1,"chem  = 'neq'" 
  printf,1,"fAin = 'chem_dummy_arc.ncdf'"
  printf,1,"!fAin  = 'chem_eq.ncdf'"
  printf,1,"fAeqout = 'chem_eq_ARC"+string(ChildProcessCount,format='(I02)')+".ncdf'"
  printf,1,"fAneqout = 'chem_neq_ARC"+string(ChildProcessCount,format='(I02)')+".ncdf'" 
  printf,1,"fcoeff = '../../chem/coeff_NASA_sc.dat'"
  printf,1,'MdH =   '+string(MdHin,format='(F13.9)')
  if (Chemistry eq 1) then for i=0,n_elements(fitabundances)-1 do printf,1,'Acst('+string(Abund_fit_index(i),format='(I03)')+') = '+string(fitabundances(i),format='(E15.9)') ;set all abundances
  ;printf,1,'COratio = '+string(COratio,format='(F14.6)')
  printf,1,'print_chem = f'  
  if (rain_switch eq 1) then printf,1,'rainout = t'
  if keyword_set(contribution) then if (Emission eq 'yes') then  printf,1,'calc_cp_chem = t'
  if (Chemistry eq 2) then printf,1,'init_eq = t'
  ;specifiy metallity of individual elements when chemistry=3
  if (Chemistry eq 3) then begin
    for i=0,n_elements(Element_fit_index)-1 do printf,1,'ele_xfactor('+string(Element_fit_index(i),format='(I03)')+') = '+string(fitabundances(i),format='(E15.9)') ;set all abundances
  endif

    ;CHEM_NEQ
  printf,1,'&CHEM_NEQ'
  printf,1,'mixing    = f'
  printf,1,'photochem = f'
  printf,1,'nmol_neq = 107'
  printf,1,'kzzcst  = '+string(kzz,format='(E15.8)')
  if (Chemistry eq 2) then printf,1,'check_lbound = t'
  printf,1,''
  printf,1,'tmax = 1.0E12'
  printf,1,'dtmax = 1.0E10'
  printf,1,''
  printf,1,'rate_limiter = t'
  printf,1,'Nmin = 1E-100'
  printf,1,'atol  = 1.0E-10'
  printf,1,'/'
  printf,1,''

  ;RADTRANS 
  printf,1,'&RADTRANS'
  printf,1,'nrays = 16'
  printf,1,'scatter = t'
  printf,1,'irrad = t'
  printf,1,"firad = '"+stellarfile+"'"    ;Stellar file
  printf,1,'rstar = '+string(rstar,format='(F5.3)')
  printf,1,'rorbit = '+string(rorbit,format='(F7.5)')
  printf,1,'murad = 0.5'
  printf,1,"fspectrum = 'ARC.emissionspec"+string(ChildProcessCount,format='(I02)')+".ncdf'"
  printf,1,"ftrans_spec = 'ARC.transitspec"+string(ChildProcessCount,format='(I02)')+".ncdf'"
  printf,1,"fcfout = 'ARC.emission.contribution"+string(ChildProcessCount,format='(I02)')+".ncdf'"
  printf,1,'/'
  printf,1,''

  ;OPACITY 
  printf,1,'&OPACITY '
  if keyword_set(fast) then begin
          ;R=500 opacities
          kapname_include=n_elements(where(opacity_switch eq 1)) ;number of opacitites to include
          kapname_arr=(where(opacity_switch eq 1))               ;array of opacitites to include
          kapnameA=strarr(22)  &             kapnameB=strarr(22)
          ;remove h- as it's handled differently                                                                                                                                                            
          if (opacity_switch[21] eq 1) then begin
            a=where(kapname_arr ne 21) & kapname_arr=kapname_arr(a) & kapname_include=kapname_include-1
          endif
          ;List from ~/pacode/atmo/src/opacity.f90
          kapnameA(1)  ="= 'H2'  ; fktab(" & kapnameB(1) =")  = '../../kabs/h2-h2_500_t5e-3_uw1116.nc'"
          kapnameA(2)  ="= 'He'  ; fktab(" & kapnameB(2) =")  = '../../kabs/h2-he_500_t5e-3_uw1116.nc'"
          kapnameA(3)  ="= 'H2O' ; fktab(" & kapnameB(3) =")  = '../../kabs/h2o_500_t5e-3_uw1116.nc'"
          kapnameA(4)  ="= 'CO2' ; fktab(" & kapnameB(4) =")  = '../../kabs/co2_500_t5e-3_uw1116.nc'"
          kapnameA(5)  ="= 'CO'  ; fktab(" & kapnameB(5) =")  = '../../kabs/co_500_t5e-3_uw1116.nc'"
          kapnameA(6)  ="= 'CH4' ; fktab(" & kapnameB(6) =")  = '../../kabs/ch4_500_t5e-3_uw1116.nc'"
          kapnameA(7)  ="= 'NH3' ; fktab(" & kapnameB(7) =")  = '../../kabs/nh3_500_t5e-3_uw1116.nc'"
          kapnameA(8)  ="= 'Na'  ; fktab(" & kapnameB(8) =")  = '../../kabs/na_500_t5e-3_uw1116.nc'"
          kapnameA(9)  ="= 'K'   ; fktab(" & kapnameB(9) =")  = '../../kabs/k_500_t5e-3_uw1116.nc'"
          kapnameA(10) ="= 'Li'  ; fktab(" & kapnameB(10) =") = '../../kabs/li_500_t5e-3_uw1116.nc'"
          kapnameA(11) ="= 'Rb'  ; fktab(" & kapnameB(11) =") = '../../kabs/rb_500_t5e-3_uw1116.nc'"
          kapnameA(12) ="= 'Cs'  ; fktab(" & kapnameB(12) =") = '../../kabs/cs_500_t5e-3_uw1116.nc'"
          kapnameA(13) ="= 'TiO' ; fktab(" & kapnameB(13) =") = '../../kabs/tio_500_t5e-3_uw1116.nc'"
          kapnameA(14) ="= 'VO'  ; fktab(" & kapnameB(14) =") = '../../kabs/vo_500_t5e-3_uw1116_vomyt.nc'"
          kapnameA(15) ="= 'FeH' ; fktab(" & kapnameB(15) =") = '../../kabs/feh_500_t5e-3_uw1116.nc'"
          kapnameA(16) ="= 'PH3' ; fktab(" & kapnameB(16) =") = '../../kabs/ph3_500_t5e-3_uw1116.nc'"
          kapnameA(17) ="= 'H2S' ; fktab(" & kapnameB(17) =") = '../../kabs/h2s_500_t5e-3_uw1116.nc'"
          kapnameA(18) ="= 'HCN' ; fktab(" & kapnameB(18) =") = '../../kabs/hcn_500_t5e-3_uw1116.nc'"
          kapnameA(19) ="= 'C2H2'; fktab(" & kapnameB(19) =") = '../../kabs/c2h2_500_t5e-3_uw1116.nc'"
          kapnameA(20) ="= 'SO2' ; fktab(" & kapnameB(20) =") = '../../kabs/so2_500_t5e-3_uw1116.nc'"
          kapnameA(21) ="= 'Fe'  ; fktab(" & kapnameB(21) =") = '../../kabs/fe_500_t5e-3_uw1116.nc'"
          printf,1,'nkap = '+string(kapname_include,format='(I02)')
          ;printf,1,"kapname(1)  = 'H2'  ; fktab(1)  = '../../kabs/h2-h2_5000_t5e-3_uw1116.nc'"
          for k=0,kapname_include-1 do printf,1,'kapname('+string(k+1,format='(I02)')+')  '+kapnameA(kapname_arr(k)+1)+string(k+1,format='(I02)')+kapnameB(kapname_arr(k)+1)
  endif else begin
          ;R=5000 opacities
          kapname_include=n_elements(where(opacity_switch eq 1)) ;number of opacitites to include
          kapname_arr=(where(opacity_switch eq 1))               ;array of opacitites to include
          kapnameA=strarr(22)  &             kapnameB=strarr(22)
          ;remove h- as it's handled differently
            if (opacity_switch[21] eq 1) then begin
              a=where(kapname_arr ne 21) & kapname_arr=kapname_arr(a) & kapname_include=kapname_include-1
            endif
          ;List from ~/pacode/atmo/src/opacity.f90
          kapnameA(1)  ="= 'H2'  ; fktab(" & kapnameB(1) =")  = '../../kabs/h2-h2_5000_t5e-3_uw1116.nc'"
          kapnameA(2)  ="= 'He'  ; fktab(" & kapnameB(2) =")  = '../../kabs/h2-he_5000_t5e-3_uw1116.nc'"
          kapnameA(3)  ="= 'H2O' ; fktab(" & kapnameB(3) =")  = '../../kabs/h2o_5000_t5e-3_uw1116.nc'"
          kapnameA(4)  ="= 'CO2' ; fktab(" & kapnameB(4) =")  = '../../kabs/co2_5000_t5e-3_uw1116.nc'"
          kapnameA(5)  ="= 'CO'  ; fktab(" & kapnameB(5) =")  = '../../kabs/co_5000_t5e-3_uw1116.nc'"
          kapnameA(6)  ="= 'CH4' ; fktab(" & kapnameB(6) =")  = '../../kabs/ch4_5000_t5e-3_uw1116.nc'"
          kapnameA(7)  ="= 'NH3' ; fktab(" & kapnameB(7) =")  = '../../kabs/nh3_5000_t5e-3_uw1116.nc'"
          kapnameA(8)  ="= 'Na'  ; fktab(" & kapnameB(8) =")  = '../../kabs/na_5000_t5e-3_uw1116.nc'"
          kapnameA(9)  ="= 'K'   ; fktab(" & kapnameB(9) =")  = '../../kabs/k_5000_t5e-3_uw1116.nc'"
          kapnameA(10) ="= 'Li'  ; fktab(" & kapnameB(10) =") = '../../kabs/li_5000_t5e-3_uw1116.nc'"
          kapnameA(11) ="= 'Rb'  ; fktab(" & kapnameB(11) =") = '../../kabs/rb_5000_t5e-3_uw1116.nc'"
          kapnameA(12) ="= 'Cs'  ; fktab(" & kapnameB(12) =") = '../../kabs/cs_5000_t5e-3_uw1116.nc'"
          kapnameA(13) ="= 'TiO' ; fktab(" & kapnameB(13) =") = '../../kabs/tio_5000_t5e-3_uw1116.nc'"
          kapnameA(14) ="= 'VO'  ; fktab(" & kapnameB(14) =") = '../../kabs/vo_5000_t5e-3_uw1116_vomyt.nc'"
          kapnameA(15) ="= 'FeH' ; fktab(" & kapnameB(15) =") = '../../kabs/feh_5000_t5e-3_uw1116.nc'"
          kapnameA(16) ="= 'PH3' ; fktab(" & kapnameB(16) =") = '../../kabs/ph3_5000_t5e-3_uw1116.nc'"
          kapnameA(17) ="= 'H2S' ; fktab(" & kapnameB(17) =") = '../../kabs/h2s_5000_t5e-3_uw1116.nc'"
          kapnameA(18) ="= 'HCN' ; fktab(" & kapnameB(18) =") = '../../kabs/hcn_5000_t5e-3_uw1116.nc'"
          kapnameA(19) ="= 'C2H2'; fktab(" & kapnameB(19) =") = '../../kabs/c2h2_5000_t5e-3_uw1116.nc'"
          kapnameA(20) ="= 'SO2' ; fktab(" & kapnameB(20) =") = '../../kabs/so2_5000_t5e-3_uw1116.nc'"
          kapnameA(21) ="= 'Fe'  ; fktab(" & kapnameB(21) =") = '../../kabs/fe_5000_t5e-3_uw1116.nc'"  
   printf,1,'nkap = '+string(kapname_include,format='(I02)')
   ;printf,1,"kapname(1)  = 'H2'  ; fktab(1)  = '../../kabs/h2-h2_5000_t5e-3_uw1116.nc'"
   for k=0,kapname_include-1 do printf,1,'kapname('+string(k+1,format='(I02)')+')  '+kapnameA(kapname_arr(k)+1)+string(k+1,format='(I02)')+kapnameB(kapname_arr(k)+1)
   ;
  endelse
  ;include artificial hazes if option is selected
  if (Hazetype eq 1) then printf,1,'art_haze = '+string(haze,format='(E15.8)')
  ;include artificial Clouds if option is selected
  if (Cloudtype eq 0) then begin 
    printf,1,'cloud = f'
  endif 
  if (Cloudtype eq 1) then begin 
    printf,1,'cloud = t'
    printf,1,'ppmax_cloud = '+string(pmax,format='(E15.8)')
    printf,1,'ppmin_cloud = '+string(cloudtop,format='(E15.8)')
    printf,1,'cloud_strength = '+string(Cloud,format='(E15.8)')
  endif
  printf,1,'kap_smooth = t'
  printf,1,'kernel_smooth = 2'
  printf,1,'/'
  printf,1,''

  ;SOLVER 
  printf,1,'&SOLVER' 
  printf,1,'solve_hydro  = f'
  printf,1,'solve_energy = f'
  printf,1,'minstep      = 1.0E-5'
  printf,1,'maxstep      = 9.0E-1'
  printf,1,'accuracy     = 1.0E-1'
  printf,1,'psurf        = 1E-6'
  printf,1,'print_err = f'
  if (OPACITY_SWITCH[21] eq 1) then printf,1,'hminus = t' ;turn on H- if requested
  if (Emission eq 'no')  then printf,1,'transmission_spectrum=t'
  if (Emission eq 'no')  then printf,1,'surface_spectrum=f'
  if (Emission eq 'yes') then printf,1,'transmission_spectrum=f'
  if (Emission eq 'yes') then printf,1,'surface_spectrum=t'
  if keyword_set(contribution) then if (Emission eq 'yes') then printf,1,'calc_cf = t'
  printf,1,'/'
  printf,1,''
  printf,1,'&CONVECTION alpha = 0.0/'

  close,1
  ;spawn,'more ARC.in'
  if keyword_set(fast) then begin
    spawn,'./atmo.x ARC'+string(ChildProcessCount,format='(I02)')+'.in',list
  endif else begin
    spawn,'./atmo.x ARC'+string(ChildProcessCount,format='(I02)')+'.in'
  endelse

end

