pro arc_monitor_nested_plot,dirs,chi2min,chi=chi
; Purpose:  Plot The best-fitting spectra from the Nested sampling run
;
result=isa(dirs) & if (result eq 0) then dirs=''
  
restore,dirs+'multinest.monitor.sav'
MULTINEST_STR=DATA_STR_TEMP
 ;get spectra data
 x=MULTINEST_STR.(10)
 y=MULTINEST_STR.(11)
 err=MULTINEST_STR.(12)
 waverr=MULTINEST_STR.(13)
 ;get models
 ws=MULTINEST_STR.(8)
 samples_model=MULTINEST_STR.(5)
 N=n_elements(samples_model (0,*))

 ;calculate chi^2 for all
 chi2=dblarr(N)
for m=0,N-1 do begin & $
model=samples_model(*,m) & $
   linterp,ws,model,x,modelinterp & $
   model_data=dblarr(n_elements(x)) & $
    for i=0,n_elements(y)-1 do begin & $
      a=where((ws ge x(i)-waverr(i)) and (ws le x(i)+waverr(i))) & $ ; index of model data points in wavelength bin
      if (a(0) ne -1) then model_data(i)=mean(Model(a)) & $ ; take mean of model datapoints if there is any
      if (a(0) eq -1) then model_data(i)=modelinterp(i) & $ ; if no model points over wavebin, use interpolation instead
    endfor & $
    modelpoints=model_data & $
chi2(m)=total((y-modelpoints)^2./err^2.) & $
print,m,' ',chi2(m) & $
endfor

if keyword_set(chi) then begin
 set_plot,'x'

   ;plot data and best-fit model
   plotsym,0,1,/fill
   cgplot,x,y,/ystyle,/xstyle,xrange=[min(x),max(x)],/nodata,title='chi2: '+string(chi2(N-1))+' sample number:'+string(N)
   oploterror,x,y,waverr,err,color=cgcolor('black'),errcolor=cgcolor('black'),psym=8
   oplot,MULTINEST_STR.(8),MULTINEST_STR.(5)(*,N-1),color=cgcolor('red'),thick=2
endif

   set_plot,'ps'
   device,/landscape,filename='arc_monitor_nested_plot.ps'
   plotsym,0,1,/fill
   cgplot,x,y,/ystyle,/xstyle,xrange=[min(x),max(x)],/nodata,title='chi2: '+string(chi2(N-1))+' sample number:'+string(N)
   oploterror,x,y,waverr,err,color=cgcolor('black'),errcolor=cgcolor('black'),psym=8
   oplot,MULTINEST_STR.(8),MULTINEST_STR.(5)(*,N-1),color=cgcolor('red'),thick=2
   device,/close
   

   chi2min=chi2(N-1)

dirsn=dirs
restore,dirs+'arc_inputs.sav'
dirs=dirsn
;==================================================================================================================================
;Best-Fit MCMC PLOT
   wsd=x & r=y & werr=waverr & rerr=err
   bestfitmodelmcmcmodelpoints=model_data
   ws=MULTINEST_STR.(8)
   bestfitmodelmcmc=MULTINEST_STR.(5)(*,N-1)

  set_plot,'ps'
  !p.multi = [0,1,2]  
  !p.font=0
  filenam='arc_monitor_nested_plot_bestfit'
  device,filename=dirs+filenam+'.ps',/helvetica,/landscape,/color
;--------------

if (Emission eq 'no') then  begin   plot,ws,bestfitmodelmcmc,xrange=[0.3,5],yrange=[min(r)-0.002,max(r)+0.002],/xlog,/xsty,thick=2.0,/ysty,/nodata,$
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
if (Emission eq 'yes') then begin
   plot,ws,bestfitmodelmcmc*1E6,xrange=[0.3,5],yrange=[0,max(r)*1E6+2*max(rerr)*1E6],/xlog,/xstyle,thick=2.0,/ysty,/nodata,$
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
;--------------

factor=double(floor(alog10(max(y))))
xrange=[min(wsd)-0.1,max(wsd)+0.1]
yrange=[min(y-2*rerr),max(y+2*rerr)]/10^(factor)
if yrange(0) lt 0 then yrange(0)=0
if (Emission eq 'DI') then begin
   plot,ws,bestfitmodelmcmc/10^(factor),xrange=xrange,yrange=yrange,/xstyle,thick=2.0,/ysty,/nodata,$
    ytitle=textoidl('Flux (\times 10^{'+string(factor,format='(I3)')+'} W/m^2/\mum )'),xtitle=textoidl('Wavelength (\mum)'),$
        title='',charthick=4,xthick=4,ythick=4,POS=[0.1,0.21,  1,1],xcharsize=0.00001,ycharsize=1,xtickinterval=0.2
   ;,Xtickv=[0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.5,3,3.5,4,5], xticks=19,xtickname=['.6','.7','.8','.9','1',' ',' ',' ',' ','1.5',' ',' ',' ',' ','2','2.5','3','3.5','4','5'],xticklen=0.01
  oplot,ws,bestfitmodelmcmc/10^(factor),thick=6,color=cgcolor('red', !D.Table_Size-2)

  oploterror,wsd,r/10^(factor),werr,rerr/10^(factor),thick=4,psym=8
    plotsym,0,1.2,/fill
  oplot,wsd,bestfitmodelmcmcmodelpoints/10^(factor),psym=8
   plotsym,0,1,/fill
  oplot,wsd,bestfitmodelmcmcmodelpoints/10^(factor),psym=8,color=cgcolor('red5', !D.Table_Size-2)
   plotsym,0,1,/fill
  oplot,ws,bestfitmodelmcmc/10^(factor),thick=6,color=cgcolor('red', !D.Table_Size-2)

 ploterror,wsd,(y-bestfitmodelmcmcmodelpoints)/10^(factor),werr,rerr/10^(factor),POS=[0.1,0,  1,0.2],xrange=xrange,/xstyle,thick=2.0,/nohat,/nodata,/ysty,$
    ytitle=textoidl('Residuals'),xtitle=textoidl('Wavelength (\mum)'),$
    title='',charthick=4,xthick=4,ythick=4,xcharsize=1,ycharsize=1,charsize=1.5,xtickinterval=0.2;,$
;Xtickv=[.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.5,3,3.5,4,5], xticks=22,xtickname=['.3','.4','.5','.6','.7','.8','.9','1',' ',' ',' ',' ','1.5',' ',' ',' ',' ','2','2.5','3','3.5','4','5'],xticklen=0.03
  plots,xrange(0),0 & plots,xrange(1),0,/continue,linestyle=2,thick=2
  oploterror,wsd,(y-bestfitmodelmcmcmodelpoints)/10^(factor),werr,rerr/10^(factor),thick=4,color=cgcolor('black', !D.Table_Size-2),psym=8,/nohat

endif

device,/close
;ps_2_pdf,filenam+'.bestfitmcmc',dirs
  spawn,'ps2eps -B -C -f '+dirs+filenam+'.ps'
  spawn,'ps2pdf -dEPSCrop -dAutoRotatePages=/None -dAutoFilterColorImages=false -dColorImageFilter=/FlateEncode -dUseFlateCompression=true '+dirs+filenam+'.ps '+dirs+filenam+'.pdf'

;=============
!p.multi=0



   
end
