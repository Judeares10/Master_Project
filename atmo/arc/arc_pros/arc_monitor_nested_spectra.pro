pro arc_monitor_nested_spectra,dirs
; Purpose:  Generate Posterior Distribution from Nested Sampling
;
; Resolution of Histograms/Controurss
;  nhistobins=15
;  nhistobinsCO=25
;  NContourlevels=15   ; Number of Contours for 2D plot10
;  colortable=0       ; Color of 2D plot  1==blue
;  bottomct= 1;  
;  res=40                         ; resolution of 2D plot
;===================
  result=isa(dirs) & if (result eq 0) then dirs=''

  dirsnow=dirs
  Result = FILE_SEARCH(dirs+'*arc_retrieval.setup.sav') 
  restore,Result[0]
  dirs=dirsnow
  low=LOWPRIORSPARAMS
  up=UPPRIORSPARAMS
  restore,dirs+'arc_inputs.sav'
  dirs=dirsnow
print,'EMISSION = ',emission
  ;Result = FILE_TEST(dirs+'*.MULTINEST.sav')  
  Result = FILE_SEARCH(dirs+'*.MULTINEST.sav') 
  restore,Result[0]
  dirs=dirsnow
  UPPRIORSPARAMS=up(tofit)
  LOWPRIORSPARAMS=low(tofit)
;LP              DOUBLE    = Array[390]
;MULTINEST_STR   STRUCT    = -> <Anonymous> Array[1]
;NESTPARS        DOUBLE    = Array[9, 390]
;NUM             LONG      =           30
;NUM_PARAMS      LONG      =            9
;PRIOR_FUNC      STRING    = 'arc_multinest_prior'
;SAMPLE_NUM      LONG      =           30
;TOFIT           LONG      = Array[9]
;TOL             DOUBLE    =        1.0000000

;restore,'hat41_transit_feeq_28Apr2019.MULTINEST.sav'

nfit=n_elements(tofit) 
filenam='arc_monitor_nested_spectra'



 ;get spectra data
 x=MULTINEST_STR.(10)
 y=MULTINEST_STR.(11)
 err=MULTINEST_STR.(12)
 waverr=MULTINEST_STR.(13)
 ;get models
 ws=MULTINEST_STR.(8)
 samples_model=MULTINEST_STR.(5)
 N=n_elements(samples_model (0,*))

 wsd=x & r=y & werr=waverr & rerr=err

 
 specmean=dblarr(n_elements(ws))
 SpecHigh=dblarr(n_elements(ws))
 SpecHigh2=dblarr(n_elements(ws))
 SpecHigh3=dblarr(n_elements(ws))
 SpecLow=dblarr(n_elements(ws))
 SpecLow2=dblarr(n_elements(ws))
 SpecLow3=dblarr(n_elements(ws))

 
;mean model points
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
 ;
;======================================================================================================================================================
; Calculate Distrubution Function using probability of each sampled
; parameter, use it to calculate the best-fit and 1,2,3-sigma ranges

Prob=(exp(multinest_str.(1)))   ;probability for each sample
parameter_prob=prob

print,'======================================='
for p=0,n_elements(ws)-1 do begin
   ;go through spectra channel-by-channel
   planetflux=transpose(samples_model(p,*)) ;all models in a given channel
   srt=sort(planetflux)
   parameter=planetflux[srt]
   parameter_prob=Prob[srt]

   ;weighted average parmaeter value
   w_ave=total(parameter*parameter_prob) 
   ;print,w_ave
   ;middle from cumulative distribution
   medpar=parameter[value_locate(total(parameter_prob,/cumulative),[0.5d0])]
   ;print,m
   e=parameter[value_locate(total(parameter_prob,/cumulative),[0.5d0-0.68268d0/2.,0.5d0+0.68268d0/2.])];-medpar[0]
   e2=parameter[value_locate(total(parameter_prob,/cumulative),[0.5d0-0.95449d0/2.,0.5d0+0.95449d0/2.])];-medpar[0]
   e3=parameter[value_locate(total(parameter_prob,/cumulative),[0.5d0-0.99730d0/2.,0.5d0+0.99730d0/2.])];-medpar[0]

   specmean(p)=medpar
   SpecLow(p)=e(0)
   SpecHigh(p)=e(1)
   SpecLow2(p)=e2(0)
   SpecHigh2(p)=e2(1)
   SpecLow3(p)=e3(0)
   SpecHigh3(p)=e3(1)
endfor


;======================================================================================================================================================


if (Emission eq 'no') then  a=where(ws ge 0.1 and ws le 30)
if (Emission eq 'yes') then a=where(ws ge 0.1 and ws le 30)

if (Emission eq 'DI') then begin
  factor=double(floor(alog10(max(y))))
  xrange=[min(wsd)-0.1,max(wsd)+0.1]
  yrange=[min(y-2*rerr),max(y+2*rerr)]/10^(factor)
  if yrange(0) lt 0 then yrange(0)=0
  a=where(ws ge 0.1 and ws le 30)
endif

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
if (Emission eq 'no') then begin 

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


  device,/close
;=============

  spawn,'ps2eps -B -C -f '+dirs+filenam+'.WFC3.ps'
  spawn,'ps2pdf -dEPSCrop -dAutoRotatePages=/None -dAutoFilterColorImages=false -dColorImageFilter=/FlateEncode -dUseFlateCompression=true '+dirs+filenam+'.WFC3.ps '+dirs+filenam+'.WFC3.pdf'
  ;spawn,'ps2pdf -dEPSCrop -dAutoRotatePages=/None -dAutoFilterColorImages=false -dColorImageFilter=/FlateEncode -dUseFlateCompression=true '+dirs+filenam+'.WFC3.eps temporary.pdf'
  ;spawn,'pdftk temporary.pdf cat 1L output '+dirs+filenam+'.pdf'
  ;spawn,'open '+dirs+filenam+'.pdf
endif

;  spawn,'rm '+dirs+filenam+'.ps'
;==================================================================================================================================

;DIRECT IMAGING
;--------------
if (Emission eq 'DI') then begin  
  set_plot,'ps'
  !p.font=0
  device,filename=dirs+filenam+'.ps',/helvetica,/landscape,/color

  plot,ws,Model/10^(factor),xrange=xrange,yrange=yrange,/xstyle,thick=2.0,/ysty,/nodata,$
    ytitle=textoidl('Flux (\times 10^{'+string(factor,format='(I3)')+'} W/m^2/\mum )'),xtitle=textoidl('Wavelength (\mum)'),$
    title='',charthick=4,xthick=4,ythick=4,ycharsize=1.5,xcharsize=1.5,xtickinterval=0.2;,$
    ;Xtickv=[0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,3,4,5],xticks=15,xtickname=['','','1',' ',' ',' ',' ','1.5',' ',' ',' ',' ','2','3','4','5'],xticklen=0.01
;plot 1-sigma range
  cgcolorfill,[ws2,reverse(ws2),ws2(0)],[SpecHigh3,reverse(SpecLow3),SpecHigh3[0]]/10^(factor),color='BLU2',noclip=0
  cgcolorfill,[ws2,reverse(ws2),ws2(0)],[SpecHigh2,reverse(SpecLow2),SpecHigh2[0]]/10^(factor),color='BLU4',noclip=0
  cgcolorfill,[ws2,reverse(ws2),ws2(0)],[SpecHigh,reverse(SpecLow),SpecHigh[0]]/10^(factor),color='BLU6',noclip=0
  oplot,ws,specmean/10^(factor),thick=6,color=cgcolor('red', !D.Table_Size-2)
  ;data
  plotsym,0,0.8,/fill
  oploterror,wsd,r/10^(factor),werr,rerr/10^(factor),thick=3,psym=8,/nohat

 ;   plotsym,0,1.2,/fill
 ; oplot,wsd,modelpoints/10^(factor),psym=8
 ;  plotsym,0,1,/fill
 ; oplot,wsd,modelpoints/10^(factor),psym=8,color=cgcolor('red5', !D.Table_Size-2)
 
  device,/close
  spawn,'ps2eps -B -C -f '+dirs+filenam+'.ps'
  spawn,'ps2pdf -dEPSCrop -dAutoRotatePages=/None -dAutoFilterColorImages=false -dColorImageFilter=/FlateEncode -dUseFlateCompression=true '+dirs+filenam+'.ps '+dirs+filenam+'.pdf'

endif
if (Emission eq 'DI') then begin  
  set_plot,'ps'
  !p.font=0
  device,filename=dirs+filenam+'Kband.ps',/helvetica,/landscape,/color
  xrange=[1.95,2.5]
  b=where(wsd ge 1.95 and wsd le 2.5)
  yrange=[min(y(b)-2*rerr(b)),max(y(b)+2*rerr(b))]/10^(factor)
  if yrange(0) lt 0 then yrange(0)=0

  plot,ws,Model/10^(factor),xrange=xrange,yrange=yrange,/xstyle,thick=2.0,/ysty,/nodata,$
    ytitle=textoidl('Flux (\times 10^{'+string(factor,format='(I3)')+'} W/m^2/\mum )'),xtitle=textoidl('Wavelength (\mum)'),$
    title='',charthick=4,xthick=4,ythick=4,ycharsize=1.5,xcharsize=1.5,xtickinterval=0.05;,$
    ;Xtickv=[0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,3,4,5],xticks=15,xtickname=['','','1',' ',' ',' ',' ','1.5',' ',' ',' ',' ','2','3','4','5'],xticklen=0.01
;plot 1-sigma range
  cgcolorfill,[ws2,reverse(ws2),ws2(0)],[SpecHigh3,reverse(SpecLow3),SpecHigh3[0]]/10^(factor),color='BLU2',noclip=0
  cgcolorfill,[ws2,reverse(ws2),ws2(0)],[SpecHigh2,reverse(SpecLow2),SpecHigh2[0]]/10^(factor),color='BLU4',noclip=0
  cgcolorfill,[ws2,reverse(ws2),ws2(0)],[SpecHigh,reverse(SpecLow),SpecHigh[0]]/10^(factor),color='BLU6',noclip=0
  oplot,ws,specmean/10^(factor),thick=6,color=cgcolor('red', !D.Table_Size-2)
  ;data
  plotsym,0,0.8,/fill
  oploterror,wsd,r/10^(factor),werr,rerr/10^(factor),thick=3,psym=8,/nohat

 ;   plotsym,0,1.2,/fill
 ; oplot,wsd,modelpoints/10^(factor),psym=8
 ;  plotsym,0,1,/fill
 ; oplot,wsd,modelpoints/10^(factor),psym=8,color=cgcolor('red5', !D.Table_Size-2)
 
  device,/close
  spawn,'ps2eps -B -C -f '+dirs+filenam+'Kband.ps'
  spawn,'ps2pdf -dEPSCrop -dAutoRotatePages=/None -dAutoFilterColorImages=false -dColorImageFilter=/FlateEncode -dUseFlateCompression=true '+dirs+filenam+'Kband.ps '+dirs+filenam+'Kband.pdf'

endif

;--------------

  


end
