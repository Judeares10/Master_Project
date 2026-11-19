pro arc_plotspec,dirs,filenam,MPFITfilenameprefix,params,wsd,r,werr,rerr,modelpoints,ws,model,fixit,Emission
  Common SystemParameters,rplanet,rstar,rorbit,stellarfile,Tstar,Cloudtype,Hazetype,Chemistry,Isothermal,Opacity_switch,Abund_fit_switch,Abund_fit_index,Abund_priorLow,Abund_priorUP,rain_switch,Element_fit_switch,Element_fit_index,Ele_priorLow,Ele_priorUP,Distance,nband
  set_plot,'ps'
  loadct,0

  ;===Transmission Spectra ==
  if (Emission eq 'no') then begin
    ;==LOAD ATMO MODEL=======
    fid=ncdf_open(dirs+MPFITfilenameprefix+'.highres.mcmc_bestfit.transitspec.ncdf')
    print,'loaded '+dirs+MPFITfilenameprefix+'.highres.mcmc_bestfit.transitspec.ncdf'
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
    linterp,ws,model,wsd,modelinterp
    ;average model over bins (wsd, werr)
    model_data=dblarr(n_elements(wsd))
    for i=0,n_elements(wsd)-1 do begin & $
      a=where((ws ge wsd(i)-werr(i)) and (ws le wsd(i)+werr(i))) & $
      if (a(0) ne -1) then model_data(i)=mean(Model(a)) & $
      if (a(0) eq -1) then model_data(i)=modelinterp(i) & $
    endfor
    modelpoints=model_data
  ;=======
  endif
  ;===Emission Spectra ==
  if (Emission eq 'yes') then begin
    ;==LOAD ATMO MODEL=======
    atmofilenam=dirs+MPFITfilenameprefix+'.highres.mcmc_bestfit.emissionspec.ncdf'
    arc_load_atmo_emission_spec,atmofilenam,ws,rplanet,modelfpfs,fnu,fstar
    model=modelfpfs
    ;
    linterp,ws,model,wsd,modelinterp
    ;average model over bins (wsd, werr)
    model_data=dblarr(n_elements(wsd))
    for i=0,n_elements(wsd)-1 do begin & $
      a=where((ws ge wsd(i)-werr(i)) and (ws le wsd(i)+werr(i))) & $
      if (a(0) ne -1) then model_data(i)=mean(Model(a)) & $
      if (a(0) eq -1) then model_data(i)=modelinterp(i) & $
    endfor
    modelpoints=model_data
  ;=======
  endif

  nfit=n_elements(where(fixit eq 0)) 
  
  chi,r,modelpoints,rerr,n_elements(r)-nfit,chiout

if (Emission eq 'no') then ytitle=textoidl('Transit Depth (R_{pl}/R_{s})')
if (Emission eq 'yes') then ytitle=textoidl('Planet-to-Star Flux Ratio (F_{pl}/F_{s})')
;===PLOT====
  set_plot,'ps'
  !p.font=0
  device,filename=dirs+filenam+'.log.ps',/helvetica,/landscape,/color
if (Emission eq 'no') then   plot,ws,Model,xrange=[0.3,5],yrange=[min(r)-0.005,max(r)+0.005],/xsty,/xlog,thick=2.0,/ysty,/nodata,$
    ytitle=ytitle,xtitle=textoidl('Wavelength (\mum)'),$
    title='',charthick=4,xthick=4,ythick=4,$
    Xtickv=[.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.5,3,3.5,4,5],xticks=22,xtickname=['.3','.4','.5','.6','.7','.8','.9','1',' ',' ',' ',' ','1.5',' ',' ',' ',' ','2','2.5','3','3.5','4','5'],xticklen=0.01

if (Emission eq 'yes') then   plot,ws,Model,xrange=[0.3,5],yrange=[1E-5,max(r)+0.005],/xsty,/xlog,/ylog,thick=2.0,/ysty,/nodata,$
    ytitle=ytitle,xtitle=textoidl('Wavelength (\mum)'),$
    title='',charthick=4,xthick=4,ythick=4,$
    Xtickv=[.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.5,3,3.5,4,5],xticks=22,xtickname=['.3','.4','.5','.6','.7','.8','.9','1',' ',' ',' ',' ','1.5',' ',' ',' ',' ','2','2.5','3','3.5','4','5'],xticklen=0.01

  oplot,ws,Model,thick=3,color=cgcolor('red', !D.Table_Size-2)
    plotsym,0,1.2,/fill
  oplot,wsd,modelpoints,psym=8
   plotsym,0,1,/fill
  oplot,wsd,modelpoints,psym=8,color=cgcolor('red5', !D.Table_Size-2)
   plotsym,0,1,/fill
  oploterror,wsd,r,werr,rerr,thick=4,psym=8
  device,/close
;=============

  spawn,'ps2eps -B -C -f '+filenam+'.log.ps'
;  spawn,'ps2pdf -dEPSCrop -dAutoRotatePages=/None -dAutoFilterColorImages=false -dColorImageFilter=/FlateEncode -dUseFlateCompression=true '+filenam+'.eps temporary.pdf'
;  spawn,'pdftk temporary.pdf cat 1L output '+filenam+'.pdf'
;  spawn,'open '+filenam+'.pdf'
  spawn,'ps2pdf -dEPSCrop -dAutoRotatePages=/None -dAutoFilterColorImages=false -dColorImageFilter=/FlateEncode -dUseFlateCompression=true '+dirs+filenam+'.log.ps '+dirs+filenam+'.log.pdf'
;  spawn,'rm '+dirs+filenam+'.ps'


;===PLOT====
  set_plot,'ps'
  !p.font=0
  device,filename=dirs+filenam+'.linear.ps',/helvetica,/landscape,/color
ymin=1E-4
if (min(r) lt 1E-4) then ymin=1E-5
if (Emission eq 'yes') then   plot,ws,Model,xrange=[min(wsd)-0.05,max(wsd)+0.5],yrange=[ymin,max(r)+0.005],/xsty,/xlog,thick=2.0,/ysty,/nodata,$
    ytitle=ytitle,xtitle=textoidl('Wavelength (\mum)'),$
    title='',charthick=4,xthick=4,ythick=4;,$
    ;Xtickv=[.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.5,3,3.5,4,5],xticks=22,xtickname=['.3','.4','.5','.6','.7','.8','.9','1',' ',' ',' ',' ','1.5',' ',' ',' ',' ','2','2.5','3','3.5','4','5'],xticklen=0.01
if (Emission eq 'no') then   plot,ws,Model,xrange=[min(wsd)-0.05,max(wsd)+0.5],yrange=[min(r)-0.005,max(r)+0.005],/xsty,/xlog,thick=2.0,/ysty,/nodata,$
    ytitle=ytitle,xtitle=textoidl('Wavelength (\mum)'),$
    title='',charthick=4,xthick=4,ythick=4;,$
    ;Xtickv=[.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.5,3,3.5,4,5],xticks=22,xtickname=['.3','.4','.5','.6','.7','.8','.9','1',' ',' ',' ',' ','1.5',' ',' ',' ',' ','2','2.5','3','3.5','4','5'],xticklen=0.01

  oplot,ws,Model,thick=3,color=cgcolor('red', !D.Table_Size-2)
    plotsym,0,1.2,/fill
  oplot,wsd,modelpoints,psym=8
   plotsym,0,1,/fill
  oplot,wsd,modelpoints,psym=8,color=cgcolor('red5', !D.Table_Size-2)
   plotsym,0,1,/fill
  oploterror,wsd,r,werr,rerr,thick=4,psym=8
  device,/close
;=============

  spawn,'ps2eps -B -C -f '+filenam+'.linear.ps'
;  spawn,'ps2pdf -dEPSCrop -dAutoRotatePages=/None -dAutoFilterColorImages=false -dColorImageFilter=/FlateEncode -dUseFlateCompression=true '+filenam+'.eps temporary.pdf'
;  spawn,'pdftk temporary.pdf cat 1L output '+filenam+'.pdf'
;  spawn,'open '+filenam+'.pdf'
  spawn,'ps2pdf -dEPSCrop -dAutoRotatePages=/None -dAutoFilterColorImages=false -dColorImageFilter=/FlateEncode -dUseFlateCompression=true '+dirs+filenam+'.linear.ps '+dirs+filenam+'.linear.pdf'
;  spawn,'rm '+dirs+filenam+'.ps'


end
