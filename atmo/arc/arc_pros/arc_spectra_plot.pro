pro arc_spectra_plot,dirs,filenam,MPFITfilenameprefix,params,wsd,r,werr,rerr,modelpoints,ws,cloudtype,Emission,Rpl;,planet  ;,model,fixit
  ;Common SystemParameters,logg,rstar,rorbit,stellarfile,Tstar,Cloudtype,Hazetype,Chemistry,Isothermal,Opacity_switch,Abund_fit_switch,Abund_fit_index

  ;arc_load_system_parameters,planet,Rp,Rjup,Rpl,Teff,logg,g,Rstar,Rs,rorbit,stellarfile,Tstar
     rplanet=Rpl ;Radius of the planet in cm.  Used for Emission Flux calculation, fnu/fstar*(rplanet/rs)^2.


!p.font=0 
filenam='Spectra.plot'
loadct,0

;set_plot,'x'
restore,dirs+'arc_inputs.sav'
restore,dirs+filenam+'arc_retrieval.setup.sav'
restore,dirs+'MPFIT.bestfit.sav' ;params Array[11], Kir=params(6)
restore,dirs+'MCMC.monitor.sav' ;allchains,allchainschi2,i,j,nchains,nthin,tofit
restore,dirs+'MCMC.Spectra.sav' ;allmodel
nfit=n_elements(tofit) 
numwavepoints=n_elements(allmodel(*,0,0));


 ;check if chains are full, reduce index by 1 if so to conform to idl array rules             
   if (n_elements(allchainschi2(*,0)) eq i) then i=i-1


  ;Chains may not find CHI^2 minimum, go rouge, ect.  Find them and destroy them
   fin=i
   strt=0 & if (i gt 4000) then strt=2000
   ;calculate average chi^2 of each chain to compare
   chiave=dblarr(nchains) & for w=0,nchains-1 do chiave(w)=mean(allchainschi2(strt:fin,w)-min(allchainschi2))
   sigma_cut=2  ;sigma cut on average chi^2 of the chains
   print,'sigma cut:',sigma_cut
   resistant_mean,chiave,sigma_cut,mean,sigma,num_rejected,goodvec=good & print,'rejected ',num_rejected,'       good:',good  
   print,'Good chains have ave CHI^2 values within delta Chi^2 :',sigma
   bad=where((chiave-mean)/(sigma_cut*sigma) gt sigma_cut)  
   print,'Bad chains have ave delta CHI^2 values of :',(chiave(bad)-mean)


  ;print out stats of best-fit model
  close,1
  openw,1,dirs+'Spectra.plot.bestfitmcmc.txt'
  printf,1,'min CHI^2 = ',min(allchainschi2(0:i,*))
  printf,1,'   CHIv^2 = ',min(allchainschi2(0:i,*))/(n_elements(wsd)-nfit)
  printf,1,'N         = ',n_elements(wsd)
  printf,1,'Nfit      = ',nfit
  printf,1,'DOF       = ',n_elements(wsd)-nfit
  printf,1,''
  for k=0,nchains-1 do if min(allchainschi2(0:i,k)) eq min(allchainschi2(0:i,*)) then begin print,k,min(allchainschi2(0:i,k),minindex),minindex & minchain=k & endif
  for j=0,nfit-1 do printf,1,parname(tofit(j))+' = '+string(allchains(j,minindex,minchain))
  close,1
  ;


  ;just keep good chains
   include=good ;[0,1,2,3,4,5,6,7,8,9,10,11] ;ALL
   trimchain=nchains-n_elements(include)
   allmodel=allmodel[*,*,include]
   allchains=allchains[*,*,include]
   allchainschi2=allchainschi2[*,include]
   j=j-trimchain
   nchains=nchains-trimchain
 
   print,'min CHI^2 = ',min(allchainschi2(0:i,*))
   for k=0,nchains-1 do if min(allchainschi2(0:i,k)) eq min(allchainschi2(0:i,*)) then begin print,k,min(allchainschi2(0:i,k),minindex),minindex & minchain=k & endif
   print,'min Parameters = ',allchains(*,minindex,minchain)
   print,'    from chain :',minchain,' index ',minindex

   x=wsd & y=r & waverr=werr 
   bestfitmodelmcmc=allmodel(*,minindex,minchain)

 


      ;; discard the burn-in (the first point at which all chains
      ;; have had at least one chi^2 lower than the median chi^2)
      medchi2 = median(allchainschi2[0:i,*])
      burnndx = 0L
      for j=0L, nchains-1 do begin & $
         tmpndx = (where(allchainschi2[0:i,j] lt medchi2))(0) & $
         if tmpndx gt burnndx then burnndx = tmpndx & $
      endfor
      ;; allows Gelman-Rubin calculation if one chain is being problematic
      burnndx = burnndx < (i-3) 


b=allmodel(*,burnndx:i,*)
a=(n_elements(b(0,*,0)))*1D0*nchains & a=long(a) ;chain length
c=reform(b,numwavepoints,a)
chi=reform(ALLCHAINSCHI2(burnndx:i,*),a)

SpecArr=c

SpecMean=dblarr(numwavepoints)
SpecHigh=dblarr(numwavepoints)
SpecLow=dblarr(numwavepoints)
SpecHigh2=dblarr(numwavepoints)
SpecLow2=dblarr(numwavepoints)
SpecHigh3=dblarr(numwavepoints)
SpecLow3=dblarr(numwavepoints)


for i=0,numwavepoints-1 do begin  & $
    e=percentiles(SpecArr(i,*)) & $;-mean(SpecArr(i,*))  & $
    SpecMean(i)=median(SpecArr(i,*)) & $
    SpecLow(i)=e(0) & $
    SpecHigh(i)=e(1) & $
    e=percentiles(SpecArr(i,*),CONFLIMIT=0.95449) & $;-mean(SpecArr(i,*))  & $
    SpecLow2(i)=e(0) & $
    SpecHigh2(i)=e(1) & $
    e=percentiles(SpecArr(i,*),CONFLIMIT=0.9973002) & $;-mean(SpecArr(i,*))  & $
    SpecLow3(i)=e(0) & $
    SpecHigh3(i)=e(1) & $
endfor

  ;===Transmission Spectra ==
  if (Emission eq 'no') then begin
    ;==LOAD ATMO MODEL=======
    fid=ncdf_open(dirs+MPFITfilenameprefix+'.lowres.mcmc_bestfit.transitspec.ncdf')
    ncdf_varget,fid,'nu',nu
    ncdf_varget,fid,'transit_radius',rprs
    ncdf_close,fid
    ws=10000./nu
    Model=SpecMean;plot mean spectra, not best-fit
    ;Cloud Deck (Cloudtype=2)
    if (Cloudtype eq 2) then begin
      cloudindex=where(Model lt p(4))
      Model(cloudindex)=p(4)
    endif
    linterp,ws,model,wsd,modelinterp
    ;average model over bins (wsd, werr)
    model_data=dblarr(n_elements(wsd))
    for i=0,n_elements(r)-1 do begin & $
      a=where((ws ge wsd(i)-werr(i)) and (ws le wsd(i)+werr(i))) & $
      if (a(0) ne -1) then model_data(i)=mean(Model(a)) & $
      if (a(0) eq -1) then model_data(i)=modelinterp(i) & $
    endfor
    modelpoints=model_data
    ;=========================
    ;nfit=n_elements(where(fixit eq 0)) 
    ;chi,r,modelpoints,rerr,n_elements(r)-nfit,chiout
    ;=======
  endif
  ;===Emission Spectra ==
  if (Emission eq 'yes') then begin
    ;==LOAD ATMO MODEL=======
    atmofilenam=dirs+MPFITfilenameprefix+'.lowres.mcmc_bestfit.emissionspec.ncdf'
    ;arc_load_atmo_emission_spec,atmofilenam,ws,modelfpfs,fnu,fstar
    arc_load_atmo_emission_spec,atmofilenam,ws,rplanet,modelfpfs,fnu,fstar
    model=SpecMean;modelfpfs
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

if (Emission eq 'no') then  a=where(ws ge 0.1 and ws le 30)
if (Emission eq 'yes') then a=where(ws ge 0.1 and ws le 30)
  ws2=ws(a) & spechigh=spechigh(a) & spechigh2=spechigh2(a) & spechigh3=spechigh3(a) & SpecLow=SpecLow(a) & SpecLow2=SpecLow2(a) & SpecLow3=SpecLow3(a)
;===PLOT====
 ;--------------
 ;OPTIAL TRASMISSION PLOT
if (Emission eq 'no') then  begin   
  set_plot,'ps'
  !p.font=0
  device,filename=dirs+filenam+'.optical.linear.ps',/helvetica,/landscape,/color
  plot,ws,specmean,xrange=[0.3,1],yrange=[min(r)-0.0025,max(r)+0.0035],/xsty,thick=2.0,/ysty,/nodata,$
    ytitle=textoidl('Planet-to-star radius ratio (R_{pl}/R_{s})'),xtitle=textoidl('Wavelength (\mum)'),$
    title='',charthick=4,xthick=4,ythick=4,ycharsize=1.5,xcharsize=1.5;,$
    ;Xtickv=[0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,3,4,5],xticks=15,xtickname=['','','1',' ',' ',' ',' ','1.5',' ',' ',' ',' ','2','3','4','5'],xticklen=0.01
    ;Xtickv=[.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.5,3,3.5,4,5], xticks=22,xtickname=['.3','.4','.5','.6','.7','.8','.9','1',' ','1.2',' ','1.4',' ','1.6',' ','1.8',' ','2','2.5','3','3.5','4','5'],xticklen=0.01

  cgcolorfill,[ws2,reverse(ws2),ws2(0)],[SpecHigh3,reverse(SpecLow3),SpecHigh3[0]],color='BLU2',noclip=0
  cgcolorfill,[ws2,reverse(ws2),ws2(0)],[SpecHigh2,reverse(SpecLow2),SpecHigh2[0]],color='BLU4',noclip=0
  cgcolorfill,[ws2,reverse(ws2),ws2(0)],[SpecHigh,reverse(SpecLow),SpecHigh[0]],color='BLU6',noclip=0
  oplot,ws,specmean,thick=6,color=cgcolor('red', !D.Table_Size-2)

    plotsym,0,1.2,/fill
  oplot,wsd,modelpoints,psym=8
   plotsym,0,1,/fill
  oplot,wsd,modelpoints,psym=8,color=cgcolor('red5', !D.Table_Size-2)
   plotsym,0,1,/fill
  oploterror,wsd,r,werr,rerr,thick=4,psym=8
  device,/close

  spawn,'ps2eps -B -C -f '+dirs+filenam+'.optical.linear.ps'
  spawn,'ps2pdf -dEPSCrop -dAutoRotatePages=/None -dAutoFilterColorImages=false -dColorImageFilter=/FlateEncode -dUseFlateCompression=true '+dirs+filenam+'.optical.linear.ps '+dirs+filenam+'.optical.linear.pdf'
endif

;--------
 ;Full TRASMISSION PLOT
if (Emission eq 'no') then begin   
  set_plot,'ps'
  !p.font=0
  device,filename=dirs+filenam+'.full.log.ps',/helvetica,/landscape,/color

  plot,ws,specmean,xrange=[0.3,5],yrange=[min(r)-0.0025,max(r)+0.0035],/xsty,/xlog,thick=2.0,/ysty,/nodata,$
    ytitle=textoidl('Planet-to-star radius ratio (R_{pl}/R_{s})'),xtitle=textoidl('Wavelength (\mum)'),$
    title='',charthick=4,xthick=4,ythick=4,ycharsize=1.5,xcharsize=1.5,$
    ;Xtickv=[0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,3,4,5],xticks=15,xtickname=['','','1',' ',' ',' ',' ','1.5',' ',' ',' ',' ','2','3','4','5'],xticklen=0.01
    Xtickv=[.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.5,3,3.5,4,5], xticks=22,xtickname=['.3','.4','.5','.6','.7','.8','.9','1',' ',' ',' ',' ','1.5',' ',' ',' ',' ','2','2.5','3','3.5','4','5'],xticklen=0.01

  cgcolorfill,[ws2,reverse(ws2),ws2(0)],[SpecHigh3,reverse(SpecLow3),SpecHigh3[0]],color='BLU2',noclip=0
  cgcolorfill,[ws2,reverse(ws2),ws2(0)],[SpecHigh2,reverse(SpecLow2),SpecHigh2[0]],color='BLU4',noclip=0
  cgcolorfill,[ws2,reverse(ws2),ws2(0)],[SpecHigh,reverse(SpecLow),SpecHigh[0]],color='BLU6',noclip=0
  oplot,ws,specmean,thick=6,color=cgcolor('red', !D.Table_Size-2)

    plotsym,0,1.2,/fill
  oplot,wsd,modelpoints,psym=8
   plotsym,0,1,/fill
  oplot,wsd,modelpoints,psym=8,color=cgcolor('red5', !D.Table_Size-2)
   plotsym,0,1,/fill
  oploterror,wsd,r,werr,rerr,thick=4,psym=8
  device,/close
  spawn,'ps2eps -B -C -f '+dirs+filenam+'.full.log.ps'
  spawn,'ps2pdf -dEPSCrop -dAutoRotatePages=/None -dAutoFilterColorImages=false -dColorImageFilter=/FlateEncode -dUseFlateCompression=true '+dirs+filenam+'.full.log.ps '+dirs+filenam+'.full.log.pdf'
endif

;--------------
if (Emission eq 'yes') then begin  
  set_plot,'ps'
  !p.font=0
  device,filename=dirs+filenam+'.ps',/helvetica,/landscape,/color

  plot,ws,Model*1E6,xrange=[0.6,5],yrange=[0,max(r)*1E6+2*max(rerr)*1E6],/xlog,/xstyle,thick=2.0,/ysty,/nodata,$
    ytitle=textoidl('Occultation Depth (ppm)'),xtitle=textoidl('Wavelength (\mum)'),$
    title='',charthick=4,xthick=4,ythick=4,ycharsize=1.5,xcharsize=1.5,$
    Xtickv=[0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,3,4,5],xticks=15,xtickname=['','','1',' ',' ',' ',' ','1.5',' ',' ',' ',' ','2','3','4','5'],xticklen=0.01
;plot 1-sigma range
  cgcolorfill,[ws2,reverse(ws2),ws2(0)],[SpecHigh3,reverse(SpecLow3),SpecHigh3[0]]*1E6,color='BLU2',noclip=0
  cgcolorfill,[ws2,reverse(ws2),ws2(0)],[SpecHigh2,reverse(SpecLow2),SpecHigh2[0]]*1E6,color='BLU4',noclip=0
  cgcolorfill,[ws2,reverse(ws2),ws2(0)],[SpecHigh,reverse(SpecLow),SpecHigh[0]]*1E6,color='BLU6',noclip=0
  oplot,ws,specmean*1E6,thick=6,color=cgcolor('red', !D.Table_Size-2)

    plotsym,0,1.2,/fill
  oplot,wsd,modelpoints*1E6,psym=8
   plotsym,0,1,/fill
  oplot,wsd,modelpoints*1E6,psym=8,color=cgcolor('red5', !D.Table_Size-2)
   plotsym,0,1,/fill
  oploterror,wsd,r*1E6,werr,rerr*1E6,thick=4,psym=8
  device,/close
  spawn,'ps2eps -B -C -f '+dirs+filenam+'.ps'
  spawn,'ps2pdf -dEPSCrop -dAutoRotatePages=/None -dAutoFilterColorImages=false -dColorImageFilter=/FlateEncode -dUseFlateCompression=true '+dirs+filenam+'.ps '+dirs+filenam+'.pdf'

endif
;--------------
longwave=where(ws gt 5)
if (Emission eq 'yes' and longwave[0] ne -1) then begin  
  set_plot,'ps'
  !p.font=0
  device,filename=dirs+filenam+'.wide.ps',/helvetica,/landscape,/color

  plot,ws,Model*1E6,xrange=[0.6,10],yrange=[0,max(r)*1E6+2*max(rerr)*1E6],/xlog,/xstyle,thick=2.0,/ysty,/nodata,$
    ytitle=textoidl('Occultation Depth (ppm)'),xtitle=textoidl('Wavelength (\mum)'),$
    title='',charthick=4,xthick=4,ythick=4,ycharsize=1.5,xcharsize=1.5,$
    Xtickv=[0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,3,4,5,6,7,8,9,10],xticks=20,xtickname=['','','1',' ',' ',' ',' ','1.5',' ',' ',' ',' ','2','3','4','5','6','7','8','9','10'],xticklen=0.01
;plot 1-sigma range
  cgcolorfill,[ws2,reverse(ws2),ws2(0)],[SpecHigh3,reverse(SpecLow3),SpecHigh3[0]]*1E6,color='BLU2',noclip=0
  cgcolorfill,[ws2,reverse(ws2),ws2(0)],[SpecHigh2,reverse(SpecLow2),SpecHigh2[0]]*1E6,color='BLU4',noclip=0
  cgcolorfill,[ws2,reverse(ws2),ws2(0)],[SpecHigh,reverse(SpecLow),SpecHigh[0]]*1E6,color='BLU6',noclip=0
  oplot,ws,specmean*1E6,thick=6,color=cgcolor('red', !D.Table_Size-2)

    plotsym,0,1.2,/fill
  oplot,wsd,modelpoints*1E6,psym=8
   plotsym,0,1,/fill
  oplot,wsd,modelpoints*1E6,psym=8,color=cgcolor('red5', !D.Table_Size-2)
   plotsym,0,1,/fill
  oploterror,wsd,r*1E6,werr,rerr*1E6,thick=4,psym=8
  device,/close
    spawn,'ps2eps -B -C -f '+dirs+filenam+'.wide.ps'
  spawn,'ps2pdf -dEPSCrop -dAutoRotatePages=/None -dAutoFilterColorImages=false -dColorImageFilter=/FlateEncode -dUseFlateCompression=true '+dirs+filenam+'.wide.ps '+dirs+filenam+'.wide.pdf'
endif

;=============

  ;spawn,'ps2eps -B -C -f '+dirs+filenam+'.wide.ps'
  ;spawn,'ps2pdf -dEPSCrop -dAutoRotatePages=/None -dAutoFilterColorImages=false -dColorImageFilter=/FlateEncode -dUseFlateCompression=true '+dirs+filenam+'.wide.ps '+dirs+filenam+'.wide.pdf'
  ;spawn,'ps2pdf -dEPSCrop -dAutoRotatePages=/None -dAutoFilterColorImages=false -dColorImageFilter=/FlateEncode -dUseFlateCompression=true '+dirs+filenam+'.wide.eps temporary.pdf'
  ;spawn,'pdftk temporary.pdf cat 1L output '+dirs+filenam+'.pdf'
  ;spawn,'open '+dirs+filenam+'.pdf'

;  spawn,'rm '+dirs+filenam+'.ps'

;==================================================================================================================================
;==================================================================================================================================
;==================================================================================================================================
;WFC3 PLOT
c=where(wsd lt 1.8 and wsd gt 1.1)
maxy=max(r(c)+err(c)*3)
miny=min(r(c)-err(c)*2)
b=where(wsd gt 0.1 and wsd le 30)
a=where(ws2 ge 0.1 and ws2 le 30)
  ws2=ws2(a) & spechigh=spechigh(a) & spechigh2=spechigh2(a) & spechigh3=spechigh3(a) & SpecLow=SpecLow(a) & SpecLow2=SpecLow2(a) & SpecLow3=SpecLow3(a)
  set_plot,'ps'
  !p.font=0
  device,filename=dirs+filenam+'.WFC3.ps',/helvetica,/landscape,/color
;--------------

if (Emission eq 'no') then  begin   plot,ws,specmean,xrange=[1,1.7],yrange=[miny,maxy],/xsty,thick=2.0,/ysty,/nodata,$
    ytitle=textoidl('Planet-to-star radius ratio (R_{pl}/R_{s})'),xtitle=textoidl('Wavelength (\mum)'),$
    title='',charthick=4,xthick=4,ythick=4,ycharsize=1.5,xcharsize=1.5
  cgcolorfill,[ws2,reverse(ws2),ws2(0)],[SpecHigh3,reverse(SpecLow3),SpecHigh3[0]],color='BLU2',noclip=0
  cgcolorfill,[ws2,reverse(ws2),ws2(0)],[SpecHigh2,reverse(SpecLow2),SpecHigh2[0]],color='BLU4',noclip=0
  cgcolorfill,[ws2,reverse(ws2),ws2(0)],[SpecHigh,reverse(SpecLow),SpecHigh[0]],color='BLU6',noclip=0
  oplot,ws,specmean,thick=6,color=cgcolor('red', !D.Table_Size-2)

    plotsym,0,1.2,/fill
  oplot,wsd,modelpoints,psym=8
   plotsym,0,1,/fill
  oplot,wsd,modelpoints,psym=8,color=cgcolor('red5', !D.Table_Size-2)
   plotsym,0,1,/fill
  oploterror,wsd,r,werr,rerr,thick=4,psym=8
endif
;--------------
if (Emission eq 'yes') then begin  plot,ws,Model*1E6,xrange=[1,1.8],yrange=[0,maxy*1E6],/xstyle,thick=2.0,/ysty,/nodata,$
    ytitle=textoidl('Occultation Depth (ppm)'),xtitle=textoidl('Wavelength (\mum)'),$
    title='',charthick=4,xthick=4,ythick=4,ycharsize=1.5,xcharsize=1.5
;plot 1-sigma range
  cgcolorfill,[ws2,reverse(ws2),ws2(0)],[SpecHigh3,reverse(SpecLow3),SpecHigh3[0]]*1E6,color='BLU2',noclip=0
  cgcolorfill,[ws2,reverse(ws2),ws2(0)],[SpecHigh2,reverse(SpecLow2),SpecHigh2[0]]*1E6,color='BLU4',noclip=0
  cgcolorfill,[ws2,reverse(ws2),ws2(0)],[SpecHigh,reverse(SpecLow),SpecHigh[0]]*1E6,color='BLU6',noclip=0
  oplot,ws,specmean*1E6,thick=6,color=cgcolor('red', !D.Table_Size-2)

    plotsym,0,1.2,/fill
  oplot,wsd,modelpoints*1E6,psym=8
   plotsym,0,1,/fill
  oplot,wsd,modelpoints*1E6,psym=8,color=cgcolor('red5', !D.Table_Size-2)
   plotsym,0,1,/fill
  oploterror,wsd,r*1E6,werr,rerr*1E6,thick=4,psym=8
endif

  device,/close
;=============

  spawn,'ps2eps -B -C -f '+dirs+filenam+'.WFC3.ps'
  spawn,'ps2pdf -dEPSCrop -dAutoRotatePages=/None -dAutoFilterColorImages=false -dColorImageFilter=/FlateEncode -dUseFlateCompression=true '+dirs+filenam+'.WFC3.ps '+dirs+filenam+'.WFC3.pdf'
  ;spawn,'ps2pdf -dEPSCrop -dAutoRotatePages=/None -dAutoFilterColorImages=false -dColorImageFilter=/FlateEncode -dUseFlateCompression=true '+dirs+filenam+'.WFC3.eps temporary.pdf'
  ;spawn,'pdftk temporary.pdf cat 1L output '+dirs+filenam+'.pdf'
  ;spawn,'open '+dirs+filenam+'.pdf'


;  spawn,'rm '+dirs+filenam+'.ps'
;==================================================================================================================================
;==================================================================================================================================
;Best-Fit MCMC PLOT

   linterp,ws,bestfitmodelmcmc,x,modelinterp
    ;average model over bins (wsd, werr)
    model_data=dblarr(n_elements(x))
    for j=0,n_elements(y)-1 do begin & $
      a=where((ws ge x(j)-waverr(j)) and (ws le x(j)+waverr(j))) & $
      if (a(0) ne -1) then model_data(j)=mean(bestfitmodelmcmc(a)) & $
      if (a(0) eq -1) then model_data(j)=modelinterp(j) & $
    endfor
    bestfitmodelmcmcmodelpoints=model_data


b=where(wsd gt 0.1 and wsd le 30)
a=where(ws2 ge 0.1 and ws2 le 30)
  ws2=ws2(a) & spechigh=spechigh(a) & spechigh2=spechigh2(a) & spechigh3=spechigh3(a) & SpecLow=SpecLow(a) & SpecLow2=SpecLow2(a) & SpecLow3=SpecLow3(a)
  set_plot,'ps'
  !p.multi = [0,1,2]  
  !p.font=0
  device,filename=dirs+filenam+'.bestfitmcmc.ps',/helvetica,/landscape,/color
;--------------

if (Emission eq 'no') then  begin   plot,ws,bestfitmodelmcmc,xrange=[0.3,5],yrange=[min(r[b])-0.002,max(r[b])+0.002],/xlog,/xsty,thick=2.0,/ysty,/nodata,$
    ytitle=textoidl('Planet-to-star radius ratio (R_{pl}/R_{s})'),xtitle=textoidl('Wavelength (\mum)'),$
    title='',charthick=4,xthick=4,ythick=4,POS=[0.1,0.21,  1,1],xcharsize=0.00001,ycharsize=1.5,$
Xtickv=[.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.5,3,3.5,4,5], xticks=22,xtickname=['.3','.4','.5','.6','.7','.8','.9','1',' ',' ',' ',' ','1.5',' ',' ',' ',' ','2','2.5','3','3.5','4','5'],xticklen=0.01
  
  oplot,ws,bestfitmodelmcmc,thick=6,color=cgcolor('red', !D.Table_Size-2)

    plotsym,0,1.2,/fill
  oplot,wsd,bestfitmodelmcmcmodelpoints,psym=8
   plotsym,0,1,/fill
  oplot,wsd,bestfitmodelmcmcmodelpoints,psym=8,color=cgcolor('red5', !D.Table_Size-2)
   plotsym,0,1,/fill
  oploterror,wsd,r,werr,rerr,thick=4,psym=8

 ploterror,wsd,(y-bestfitmodelmcmcmodelpoints),werr,rerr,POS=[0.1,0,  1,0.2],xrange=[0.3,5],/xlog,/xstyle,thick=2.0,/nohat,/ysty,/nodata,$
    ytitle=textoidl('Residuals'),xtitle=textoidl('Wavelength (\mum)'),$
    title='',charthick=4,xthick=4,ythick=4,xcharsize=1,ycharsize=1,charsize=1.5,$
Xtickv=[.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.5,3,3.5,4,5], xticks=22,xtickname=['.3','.4','.5','.6','.7','.8','.9','1',' ',' ',' ',' ','1.5',' ',' ',' ',' ','2','2.5','3','3.5','4','5'],xticklen=0.03
  plots,0.3,0 & plots,5,0,/continue,linestyle=2,thick=2
  oploterror,wsd,y-bestfitmodelmcmcmodelpoints,werr,rerr,thick=4,color=cgcolor('black', !D.Table_Size-2),psym=8,/nohat

endif
;--------------
if (Emission eq 'yes') then begin  plot,ws,bestfitmodelmcmc*1E6,xrange=[0.3,5],yrange=[0,max(r)*1E6+2*max(rerr)*1E6],/xlog,/xstyle,thick=2.0,/ysty,/nodata,$
    ytitle=textoidl('Occultation Depth (ppm)'),xtitle=textoidl('Wavelength (\mum)'),$
    title='',charthick=4,xthick=4,ythick=4,POS=[0.1,0.21,  1,1],xcharsize=0.00001,ycharsize=1.5,$
Xtickv=[0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.5,3,3.5,4,5], xticks=19,xtickname=['.6','.7','.8','.9','1',' ',' ',' ',' ','1.5',' ',' ',' ',' ','2','2.5','3','3.5','4','5'],xticklen=0.01
  oplot,ws,bestfitmodelmcmc*1E6,thick=6,color=cgcolor('red', !D.Table_Size-2)

    plotsym,0,1.2,/fill
  oplot,wsd,bestfitmodelmcmcmodelpoints*1E6,psym=8
   plotsym,0,1,/fill
  oplot,wsd,bestfitmodelmcmcmodelpoints*1E6,psym=8,color=cgcolor('red5', !D.Table_Size-2)
   plotsym,0,1,/fill
  oploterror,wsd,r*1E6,werr,rerr*1E6,thick=4,psym=8

 ploterror,wsd,(y-bestfitmodelmcmcmodelpoints)*1E6,werr,rerr*1E6,POS=[0.1,0,  1,0.2],xrange=[0.3,5],/xlog,/xstyle,thick=2.0,/nohat,/nodata,/ysty,$
    ytitle=textoidl('Residuals (ppm)'),xtitle=textoidl('Wavelength (\mum)'),$
    title='',charthick=4,xthick=4,ythick=4,xcharsize=1,ycharsize=1,charsize=1.5,$
Xtickv=[.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.5,3,3.5,4,5], xticks=22,xtickname=['.3','.4','.5','.6','.7','.8','.9','1',' ',' ',' ',' ','1.5',' ',' ',' ',' ','2','2.5','3','3.5','4','5'],xticklen=0.03
  plots,0.3,0 & plots,5,0,/continue,linestyle=2,thick=2
  oploterror,wsd,(y-bestfitmodelmcmcmodelpoints)*1E6,werr,rerr*1E6,thick=4,color=cgcolor('black', !D.Table_Size-2),psym=8,/nohat

endif

device,/close
;ps_2_pdf,filenam+'.bestfitmcmc',dirs
  spawn,'ps2eps -B -C -f '+dirs+filenam+'.bestfitmcmc.ps'
  spawn,'ps2pdf -dEPSCrop -dAutoRotatePages=/None -dAutoFilterColorImages=false -dColorImageFilter=/FlateEncode -dUseFlateCompression=true '+dirs+filenam+'.bestfitmcmc.ps '+dirs+filenam+'.bestfitmcmc.pdf'

;=============
!p.multi=0

  save,filename=dirs+'arc_spectra_plot.sav',/var,/COMPRESS


end
