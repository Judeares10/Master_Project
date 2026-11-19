FUNCTION arc_model_chi_mpfit,p,X=x,Y=y,err=err,waverr=waverr,plotit=plotit
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


 Result = FILE_TEST('ARC.transitspec'+string(ChildProcessCount,format='(I02)')+'.ncdf')  
 if (Result eq 1 ) then spawn,'rm ARC.transitspec'+string(ChildProcessCount,format='(I02)')+'.ncdf'
 Result = FILE_TEST('ARC.emissionspec'+string(ChildProcessCount,format='(I02)')+'.ncdf')  
 if (Result eq 1 ) then spawn,'rm ARC.emissionspec'+string(ChildProcessCount,format='(I02)')+'.ncdf'
 Result = FILE_TEST('chem_eq_ARC'+string(ChildProcessCount,format='(I02)')+'.ncdf')  
 if (Result eq 1 ) then spawn,'rm chem_eq_ARC'+string(ChildProcessCount,format='(I02)')+'.ncdf'

 ;make copy P-T profile for use with each daughter process
 if (isothermal eq 0) then spawn,'cp '+dirs+'pt_ARC_start.ncdf pt_ARC'+string(ChildProcessCount,format='(I02)')+'.ncdf'
 ptfile='pt_ARC'+string(ChildProcessCount,format='(I02)')+'.ncdf' ;P-T profile which has the needed pressure grid
 ;Run T-P profile
  if (isothermal eq 0) then print,'TP ',Kir,gamma1,gamma2,beta,alpha,Tint,delta,Pdelta 
  if (isothermal eq 0) then arc_TP_delta,ptfile,logg,Kir,gamma1,gamma2,beta,alpha,T,Tint,Press,delta,Pdelta

  ;RUN ATMO
  if (Emission eq 'yes' or Emission eq 'DI') then arc_make_atmo_in,teff,MdH,logg,Rp,Cloud,haze,cloudtop,art_haze_alpha,kzz,fsed,sigma_g,fitabundances,/fast;,/contribution
  if (Emission eq 'no') then arc_make_atmo_in,teff,MdH,logg,Rp,Cloud,haze,cloudtop,art_haze_alpha,kzz,fsed,sigma_g,fitabundances,/fast

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
  ;=======
  endif

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
      a=where((ws ge x(i)-waverr(i)) and (ws le x(i)+waverr(i))) & $ ; index of model data points in wavelength bin
      if (a(0) ne -1) then model_data(i)=mean(Model(a)) & $ ; take mean of model datapoints if there is any
      if (a(0) eq -1) then model_data(i)=modelinterp(i) & $ ; if no model points over wavebin, use interpolation instead
    endfor
    modelpoints=model_data
 endif
;==============
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

print,'model ave Rp/Rs=',mean(modelpoints)
print,'Chi2',total((y-modelpoints)^2./err^2.),'  Pars:',p

;Plot to screen if requested
if (plotit eq 'yes') then begin
  set_plot,'x'
  if (isothermal eq 0) then wset,1
  if (isothermal eq 0) then plot,T,press,/ylog,yrange=[max(press),min(press)],xrange=[min(T)-200,max(T)+200],/ystyle


  xmin=min(x-2*waverr)
  xmax=max(x+2*waverr)
  ymin=min(y)-0.2*min(y)
  if (ymin lt 0) then ymin=0
  ymax=max(y)+0.2*max(y)
  wset,0
  ploterror,x,y,waverr,err,/xlog,/xstyle,xrange=[xmin,xmax],yrange=[ymin,ymax],/ystyle,title='MPFIT'
  oplot,ws,model,color=1000
  plotsym,0,1.2,/fill
  oplot,x,modelpoints,psym=8
  plotsym,0,1,/fill
  oplot,x,modelpoints,psym=8,color=cgcolor('red5', !D.Table_Size-2)
  ;xyouts,1.5,max(y),string(total((y-modelpoints)^2./err^2.))+' chi^2'
  ;xyouts,1.5,max(y)-0.0002,string(p(0))+' Temp'
  ;xyouts,1.5,max(y)-0.0004,string(p(1))+' ln metalicity'
  ;xyouts,1.5,max(y)-0.0006,string(p(2))+' logg'
  ;xyouts,1.5,max(y)-0.0008,string(p(3))+' Rp'
  ;xyouts,1.5,max(y)-0.0010,string(p(4))+' Cloud'
  ;xyouts,1.5,max(y)-0.0012,string(p(5))+' Haze'
  ;xyouts,1.5,max(y)-0.0014,string(p(8))+' Kzz'
  ;if (Isothermal eq 0) then  xyouts,1.5,max(y)-0.0016,string(p(11))+' ln (Kir)'
  ;if (Isothermal eq 0) then  xyouts,1.5,max(y)-0.0018,string(p(12))+' ln (Gam1)'
  ;if (Isothermal eq 0) then  xyouts,1.5,max(y)-0.0020,string(p(13))+' ln (Gam2)'
  ;if (Isothermal eq 0) then  xyouts,1.5,max(y)-0.0022,string(p(14))+' ln (Beta)'
  ;if (Isothermal eq 0) then  xyouts,1.5,max(y)-0.0024,string(p(15))+' ln (Alpha)'
endif

save,filename=dirs+'MpfitModel.sav',modelpoints,ws,model,x,y,waverr,err
return,(y-modelpoints)/err

end
