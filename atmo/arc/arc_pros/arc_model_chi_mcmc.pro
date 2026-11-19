FUNCTION arc_model_chi_mcmc,pars, determinant=determinant
; Read in parameters, make ATMO model, and determine CHI^2
determinant = 1d0
p=pars
prior=0D0
  Common ARC_retrieval_data, x,y,err,waverr
  Common SystemParameters,rplanet,rstar,rorbit,stellarfile,Tstar,Cloudtype,Hazetype,Chemistry,Isothermal,Opacity_switch,Abund_fit_switch,Abund_fit_index,Abund_priorLow,Abund_priorUP,rain_switch,Element_fit_switch,Element_fit_index,Ele_priorLow,Ele_priorUP,Distance,nband
  Common cpcount,ChildProcessCount,dirs,allmodel,allmodelpoints,Emission,nlevel,pmin,pmax,allchemistry,allTP
  common passmodel,ws,model,modelpoints,LowPriorsParams,UpPriorsParams,abundances,pressure,temperature
  Compile_Opt DEFINT32

 ;convert fit-parameters where necessary
 teff          =   p(0)
 MdH           =exp(p(1))
 logg          =   p(2)
 Rp            =   p(3) ;Plantary Radius at 1 bar
 Cloud         =exp(p(4))
 haze          =exp(p(5))
 cloudtop      =10^p(6) ;placeholder p(6) Cloud TOP
 art_haze_alpha=abs(p(7)) ;force to be possitive slope
 kzz           =10^p(8)  ;KZZ
 fsed          =   p(9)
 sigma_g       =   p(10)
 Kir           =10^p(11) 
 gamma1        =10^p(12) 
 gamma2        =10^p(13) 
 beta          =   p(14)  
 alpha         =   p(15) 
 Tint          =   p(16)
 delta         =   p(17) ;P-T lapse rate
 Pdelta        =10^p(18) ;Pressure lapserate base

 ; feed ATMO the fit molecular abundances for Acst(i)
  if (chemistry eq 1) then fitabundances=exp(p(19:n_elements(p)-1))
  
 ; feed ATMO the fit element abundances for ele_xfactor(i)
 ; Removes MdH=p(1) value first, as all elements are initially multipled by MdH in ATMO chemistry.f90
 ; before multipling by ele_xfactor.
 ; This removes element parameters dependance on MdH
 if (chemistry eq 3) then fitabundances=exp(p(19:n_elements(p)-1)-p(1))

 offender=''
 ;==Metallcity Prior
 ;test physical bounding priors (IDL maximum value is 1D308)
 ;ATMO Crashes for MdH > 4.50, Set an upper limit on Metalicity of 10^4.5=31622x solar
 MdHin=alog10(MdH)
 if (MdHin gt 4.50) then begin & prior=1D308 & offender=offender+'MdH ' & endif
 if (p(1) lt LowPriorsParams(1)) then begin & prior=1D308 & offender=offender+'MdHLow at '+string(p(1)) & endif
 if (p(1) gt UpPriorsParams(1)) then begin & prior=1D308 & offender=offender+'MdHhigh at '+string(p(1)) & endif
 ;==Temperature Prior
 if ((isothermal eq 1) and (teff lt LowPriorsParams(0))) then begin & prior=1D308 & offender=offender+'Teff_low at '+string(p(0)) & endif ;set maximum value to chi2
 if ((isothermal eq 1) and (teff gt  UpPriorsParams(0))) then begin & prior=1D308 & offender=offender+'Teff_hi  at '+string(p(0)) & endif ;set maximum value to chi2

 ;==logg ratio Prior
 ;if (logg le alog(1D-9)) then prior=1D308 ;set lower limit on COratio 
 if (p(2) lt (LowPriorsParams(2))) then begin & prior=1D308 & offender=offender+'logglow at '+string(p(2)) & endif ;set lower limit on COratio 
 if (p(2) gt ( UpPriorsParams(2)))  then begin & prior=1D308 & offender=offender+'logghi at '+string(p(2)) & endif    ;set upper limit on COratio 

 ;== clouds and haze ratio Priors
 for i=4,7 do begin
   if (p(i) lt (LowPriorsParams(i))) then begin & prior=1D308 & offender=offender+'clouds/haze low at '+string(p(i)) & endif ;set lower limit on COratio 
   if (p(i) gt ( UpPriorsParams(i)))  then begin & prior=1D308 & offender=offender+'clouds/haze hi at '+string(p(i)) & endif    ;set upper limit on COratio 
   ;print,i,p(i),LowPriorsParams(i),UpPriorsParams(i)
 endfor
 ;== kzz Prior
 for i=8,8 do begin
   if (p(i) lt (LowPriorsParams(i))) then begin & prior=1D308 & offender=offender+'kzz low at '+string(p(i)) & endif ;set lower limit on COratio 
   if (p(i) gt ( UpPriorsParams(i)))  then begin & prior=1D308 & offender=offender+'kzz hi at '+string(p(i)) & endif    ;set upper limit on COratio 
 endfor

 ;==Abundance Prior
 if (chemistry eq 1) then  fitted=where(Abund_fit_switch eq 0)+19 ;  select abundances that are being fit 
 ;if (chemistry eq 1 and Isothermal eq 0) then  fitted=where(Abund_fit_switch eq 0)+15 ;  select abundances that are being fit
 if (chemistry eq 1) then  fit=where(Abund_fit_switch eq 0) ;  select abundances that are being fit
 if (chemistry eq 1) then  for i=0,n_elements(fitted)-1 do if (p(fitted(i)) lt Abund_priorLow(fit(i))) then begin & prior=1D308  & offender=offender+'AbundL '+string(p(fitted(i))) & endif;set lower limit on fitted abundances 
 if (chemistry eq 1) then  for i=0,n_elements(fitted)-1 do if (p(fitted(i)) gt Abund_priorUp(fit(i))) then begin & prior=1D308  & offender=offender+'AbundU '+string(p(fitted(i))) & endif;set upper limit on fitted abundances 

 if (chemistry eq 3) then  fitted=where(Element_fit_switch eq 0)+19                                                                               ;  select abundances that are being fit
 ;if (chemistry eq 3 and Isothermal eq 0) then  fitted=where(Element_fit_switch eq 0)+15 ;  select abundances that are being fit
 if (chemistry eq 3) then  fit=where(Element_fit_switch eq 0) ;  select abundances that are being fit
 if (chemistry eq 3) then  for i=0,n_elements(fitted)-1 do if (p(fitted(i)) lt Ele_priorLow(fit(i))) then begin & prior=1D308  & offender=offender+'AbundL '+string(p(fitted(i))) & endif;set lower limit on fitted abundances 
 if (chemistry eq 3) then  for i=0,n_elements(fitted)-1 do if (p(fitted(i)) gt Ele_priorUp(fit(i))) then begin & prior=1D308  & offender=offender+'AbundU '+string(p(fitted(i))) & endif;set upper limit on fitted abundances 

    
 ;== Temperature struture priors ====
 if (Isothermal eq 0) then if (p(11) gt  UpPriorsParams(11)) then begin & prior=1D308  & offender=offender+' KirH '+string(p(11)) & endif;upper limit of Kir
 if (Isothermal eq 0) then if (p(11) lt LowPriorsParams(11)) then begin & prior=1D308  & offender=offender+' KirL '+string(p(11)) & endif ;lower limit of Kir
 if (Isothermal eq 0) then if (p(12) gt  UpPriorsParams(12)) then begin & prior=1D308  & offender=offender+' Gam1H '+string(p(12)) & endif ;upper limit of gamma1
 if (Isothermal eq 0) then if (p(12) lt LowPriorsParams(12)) then begin & prior=1D308  & offender=offender+' Gam1L '+string(p(12)) & endif ;lower limit of gamma1
 if (Isothermal eq 0) then if (p(13) gt  UpPriorsParams(13)) then begin & prior=1D308  & offender=offender+' Gam2H '+string(p(13)) & endif ;upper limit of gamma2
 if (Isothermal eq 0) then if (p(13) lt LowPriorsParams(13)) then begin & prior=1D308  & offender=offender+' Gam2L '+string(p(13)) & endif ;lower limit of gamma2
 if (Isothermal eq 0) then if (p(14) gt  UpPriorsParams(14)) then begin & prior=1D308  & offender=offender+' BetaH '+string(p(14)) & endif ;upper limit of beta
 if (Isothermal eq 0) then if (p(14) lt LowPriorsParams(14)) then begin & prior=1D308  & offender=offender+' BetaL '+string(p(14)) & endif ;lower limit of beta
 if (Isothermal eq 0) then if (p(15) gt  UpPriorsParams(15)) then begin & prior=1D308  & offender=offender+' AlphH '+string(p(15)) & endif ;lower limit of alpha
 if (Isothermal eq 0) then if (p(15) lt LowPriorsParams(15)) then begin & prior=1D308  & offender=offender+' AlphL '+string(p(15)) & endif ;lower limit of alpha
 if (Isothermal eq 0) then if (p(16) gt  UpPriorsParams(16)) then begin & prior=1D308  & offender=offender+' TintH '+string(p(16)) & endif ;lower limit of alpha
 if (Isothermal eq 0) then if (p(16) lt LowPriorsParams(16)) then begin & prior=1D308  & offender=offender+' TintL '+string(p(16)) & endif ;lower limit of alpha
 if (Isothermal eq 0) then if (p(17) gt  UpPriorsParams(17)) then begin & prior=1D308  & offender=offender+' delH '+string(p(17)) & endif ;upper limit of Kir
 if (Isothermal eq 0) then if (p(17) lt LowPriorsParams(17)) then begin & prior=1D308  & offender=offender+' delL '+string(p(17)) & endif ;upper limit of Kir
 if (Isothermal eq 0) then if (p(18) gt  UpPriorsParams(18)) then begin & prior=1D308  & offender=offender+' PdeH '+string(p(18)) & endif ;upper limit of Kir
 if (Isothermal eq 0) then if (p(18) lt LowPriorsParams(18)) then begin & prior=1D308  & offender=offender+' PdeL '+string(p(18)) & endif   ;upper limit of Kir   

 ;== Reject Model if fit is out of prior bounds
 if (prior ge 1D308) then chi2=1D308 ;set chi2 to maximum value, so step will be rejected
 if (prior ge 1D308) then print,offender+'Prior out of bounds from from Child ',ChildProcessCount
 if (prior ge 1D308) then goto,skipmodelrun_and_reject

 ;===Cleanup previous models
 if (Emission eq 'no') then begin
   Result = FILE_TEST('ARC.transitspec'+string(ChildProcessCount,format='(I02)')+'.ncdf')  
   if (Result eq 1 ) then spawn,'rm ARC.transitspec'+string(ChildProcessCount,format='(I02)')+'.ncdf'
 endif
 if (Emission eq 'yes' or Emission eq 'DI') then begin
   Result = FILE_TEST('ARC.emissionspec'+string(ChildProcessCount,format='(I02)')+'.ncdf')  
   if (Result eq 1 ) then spawn,'rm ARC.emissionspec'+string(ChildProcessCount,format='(I02)')+'.ncdf'
  endif

 Result = FILE_TEST('chem_eq_ARC'+string(ChildProcessCount,format='(I02)')+'.ncdf')  
 if (Result eq 1 ) then spawn,'rm chem_eq_ARC'+string(ChildProcessCount,format='(I02)')+'.ncdf'

 Result = FILE_TEST('chem_neq_ARC'+string(ChildProcessCount,format='(I02)')+'.ncdf')  
 if (Result eq 1 ) then spawn,'rm chem_neq_ARC'+string(ChildProcessCount,format='(I02)')+'.ncdf'

 ;===
 ;make copy P-T profile for use with each daughter process
 if (isothermal eq 0) then spawn,'cp '+dirs+'pt_ARC_start.ncdf pt_ARC'+string(ChildProcessCount,format='(I02)')+'.ncdf'
 ptfile='pt_ARC'+string(ChildProcessCount,format='(I02)')+'.ncdf' ;P-T profile which has the needed pressure grid
 ;Run T-P profile if needed
  if (isothermal eq 0) then  arc_TP_delta,ptfile,logg,Kir,gamma1,gamma2,beta,alpha,T,Tint,Press,delta,Pdelta

;======***Run ATMO**=====
 arc_make_atmo_in,teff,MdH,logg,Rp,Cloud,haze,cloudtop,art_haze_alpha,kzz,fsed,sigma_g,fitabundances,/fast
 ;test to see if file is there... if not then it likley crashed... if so reject model and move on
 if (Emission eq 'no') then begin
   Result = FILE_TEST('ARC.transitspec'+string(ChildProcessCount,format='(I02)')+'.ncdf')  
   if (Result eq 0 ) then prior=1D308
 endif
 if (Emission eq 'yes' or Emission eq 'DI') then begin
   Result = FILE_TEST('ARC.emissionspec'+string(ChildProcessCount,format='(I02)')+'.ncdf')  
   if (Result eq 0 ) then prior=1D308
  endif
 ;== Reject Model if fit is out of prior bounds
 if (prior ge 1D308) then chi2=1D308 ;set chi2 to maximum value, so step will be rejected
 if (prior ge 1D308) then print,offender+'Prior out of bounds from from Crashed Child ',ChildProcessCount
 if (prior ge 1D308) then goto,skipmodelrun_and_reject

  ;===Transmission Spectra ==
  if (Emission eq 'no') then begin
    ;==LOAD ATMO MODEL=======
    fid=ncdf_open('ARC.transitspec'+string(ChildProcessCount,format='(I02)')+'.ncdf')
    ncdf_varget,fid,'nu',nu
    ncdf_varget,fid,'transit_radius',rprs
    ncdf_close,fid
    ws=10000./nu
    Model=rprs;+p(3)
    ;Cloud Deck (Cloudtype=2)
    if (Cloudtype eq 2) then begin
      cloudindex=where(Model lt p(4))
      Model(cloudindex)=p(4)
    endif
    linterp,ws,model,x,modelinterp
    ;average model over bins (wsd, werr)
    model_data=dblarr(n_elements(x))
    for i=0,n_elements(y)-1 do begin & $
      a=where((ws ge x(i)-waverr(i)) and (ws le x(i)+waverr(i))) & $
      if (a(0) ne -1) then model_data(i)=mean(Model(a)) & $
      if (a(0) eq -1) then model_data(i)=modelinterp(i) & $
    endfor
    modelpoints=model_data
  endif
  ;============================
  
  ;===Emission Spectra ==
  if (Emission eq 'yes') then begin
    ;==LOAD ATMO MODEL=======
    atmofilenam='ARC.emissionspec'+string(ChildProcessCount,format='(I02)')+'.ncdf'
    arc_load_atmo_emission_spec,atmofilenam,ws,rplanet,modelfpfs,fnu,fstar
    model=modelfpfs
    ;
    linterp,ws,model,x,modelinterp
    ;average model over bins (wsd, werr)
    model_data=dblarr(n_elements(x))
    for i=0,n_elements(y)-1 do begin & $
      a=where((ws ge x(i)-waverr(i)) and (ws le x(i)+waverr(i))) & $
      if (a(0) ne -1) then model_data(i)=mean(Model(a)) & $
      if (a(0) eq -1) then model_data(i)=modelinterp(i) & $
    endfor
    modelpoints=model_data
  endif
  ;===================================

  ;===Direct Imaging Spectra ==
  if (Emission eq 'DI') then begin
    ;==LOAD ATMO MODEL=======
    atmofilenam='ARC.emissionspec'+string(ChildProcessCount,format='(I02)')+'.ncdf'
    arc_load_atmo_directimaging_spec,atmofilenam,Rp,Distance,ws,modelfpfs
    model=modelfpfs
    ;
    linterp,ws,model,x,modelinterp
    ;average model over bins (wsd, werr)
    model_data=dblarr(n_elements(x))
    for i=0,n_elements(y)-1 do begin & $
      a=where((ws ge x(i)-waverr(i)) and (ws le x(i)+waverr(i))) & $ ; index of model data points in wavelength bin
      if (a(0) ne -1) then model_data(i)=mean(Model(a)) & $ ; take mean of model datapoints if there is any
      if (a(0) eq -1) then model_data(i)=modelinterp(i) & $ ; if no model points over wavebin, use interpolation instead
    endfor
    modelpoints=model_data
 endif
 ;==============

chi2=total((y-modelpoints)^2./err^2.)

;read all the abundances to save
filenam='chem_eq_ARC'+string(ChildProcessCount,format='(I02)')+'.ncdf'
Result = FILE_TEST('chem_eq_ARC'+string(ChildProcessCount,format='(I02)')+'.ncdf')  
if (chemistry eq 1 or chemistry eq 3 and Result eq 1) then begin & $
   arc_readchem,'',filenam,abundances,pressure,temperature & $
endif else begin & $
   abundances=0 & pressure=0 & temperature=0 & $
endelse
;read all the NEQ abundances to save
filenam='chem_neq_ARC'+string(ChildProcessCount,format='(I02)')+'.ncdf'
Result = FILE_TEST('chem_neq_ARC'+string(ChildProcessCount,format='(I02)')+'.ncdf')  
if (chemistry eq 2 and Result eq 1 ) then begin  & $
   arc_readchem,'',filenam,abundances,pressure,temperature  & $
endif else begin & $
   abundances=0 & pressure=0 & temperature=0 & $
endelse

;*************************************************
;; Gaussian Priors can be added (not yet implimented)
;; add priors (priors[1,*] should be infinity for no prior)
;npriors = n_elements(priors)
;if npriors ne 0 and npriors ne 1 then begin
;   if npriors ne n_elements(pars)*2 then $
;      message, 'ERROR: PRIORS array must be a 2xNPARS array'
;   chi2 += total(((pars-priors[0,*])/priors[1,*])^2)
;endif
;*************************************************


print,chi2,' from Child ',ChildProcessCount
skipmodelrun_and_reject:
return,chi2

end
