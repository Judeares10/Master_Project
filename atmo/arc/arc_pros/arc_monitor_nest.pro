pro arc_monitor_nest,dirs
; Monitor MCMC convergence
; To be run in a separte IDL session to ARC
; Loads current MCMC chain, calculates Gelman-Rubin statistics
; Requires: Coyote Library http://www.idlcoyote.com/documents/programs.php
!p.font=0 
loadct,0
set_plot,'x'





restore,dirs+'MPFIT.bestfit.sav' & tofit=where(parinfo.fixed eq 0) 
nfit=n_elements(tofit) 

result=file_search(dirs,'multinest.monitor.sav')
print,'restoring ',result(0)
restore,result(0) ;,multinest_str,NESTpars,LP,LOG_LIKELIHOOD_FUNC,PRIOR_FUNC,NUM_PARAMS,NUM,SAMPLE_NUM,EXPAND,TOL,tofit	
        NESTpars=data_str_temp.(0)
        LOGZ=data_str_temp.(2)
        LOGZ_ERR=data_str_temp.(3)
        LP=data_str_temp.(1)

 
;Number of bins in Histogram
histobins=25;100;28 ;45
res=12;100;50;25
denres=90;30;45;90



;plot progression
result=ascii_template(dirs+'multinest.monitor.log.txt')
data=read_ascii(dirs+'multinest.monitor.log.txt',template=result)
t=data.(0)
deltaz=data.(1)
if n_elements(t eq 1) then goto,skipthis ;check if there are enough data points
;                                 fac=1  & xtitle='Time (seconds)' 
;if (max(t-t(0))/60. lt 60 and max(t-t(0))/60. gt 1) then begin & fac=60  & xtitle='Time (minutes)' & endif 
;if (max(t-t(0))/3600. lt 10 and max(t-t(0))/3600. gt 1) then begin & fac=3600  & xtitle='Time (hours)' & endif
;if (max(t-t(0))/86400. gt 1.5) then begin & fac=86400  & xtitle='Time (days)' & endif

ymin=0.01
if (ymin gt tol) then ymin=tol
set_plot,'ps'
  !p.font=0
device,/landscape,filename=dirs+'multinest.monitor.ps',/helvetica
;fac=1  & xtitle='Time (seconds)' 
;  plot,(t-t(0))/fac,(deltaz) ,psym=-1,/ylog,yrange=[0.1,max(deltaz)],xrange=[0.1,max((t-t(0))/fac)*20],xtitle=xtitle,ytitle='Delta Z',thick=2,charthick=2,/xlog
;  oplot,(t-t(0))/fac,(deltaz),psym=-1,thick=2
;  x=(t-t(0))/fac & y=deltaz
;  a=poly_fit(alog10(x[1:n_elements(x)-1]),alog10(y[1:n_elements(x)-1]),1,/double,yfit=yfit) 
;  XX=findgen(100)*0.1+x(1) 
;  oplot,xx,10^(a(0)+alog10(XX)*a(1))
;  plots,0,tol & plots,max((t-t(0))/fac)*2,tol,linestyle=1,/continue,thick=2     

;fac=60  & xtitle='Time (minutes)' 
;  plot,(t-t(0))/fac,(deltaz) ,psym=-1,/ylog,yrange=[0.1,max(deltaz)],xrange=[0.1,max((t-t(0))/fac)*20],xtitle=xtitle,ytitle='Delta Z',thick=2,charthick=2,/xlog
;  oplot,(t-t(0))/fac,(deltaz),psym=-1,thick=2
;  x=(t-t(0))/fac & y=deltaz
;  a=poly_fit(alog10(x[1:n_elements(x)-1]),alog10(y[1:n_elements(x)-1]),1,/double,yfit=yfit) 
;  XX=findgen(100)*0.1+x(1) 
;  oplot,xx,10^(a(0)+alog10(XX)*a(1))
;  plots,0,tol & plots,max((t-t(0))/fac)*2,tol,linestyle=1,/continue,thick=2     
fac=3600  & xtitle='Time (hours)'
  plot,(t-t(0))/fac,(deltaz) ,psym=-1,/ylog,yrange=[0.1*Tol,max(deltaz)],xrange=[0.1,1000],xtitle=xtitle,ytitle='Delta Z',thick=2,charthick=2,/xlog
  oplot,(t-t(0))/fac,(deltaz),psym=-1,thick=2
  x=(t-t(0))/fac & y=deltaz
  a=poly_fit(alog10(x[1:n_elements(x)-1]),alog10(y[1:n_elements(x)-1]),1,/double,yfit=yfit) 
  XX=findgen(2000)*0.1+x(1) 
  oplot,xx,10^(a(0)+alog10(XX)*a(1))
  plots,0.1,tol & plots,1000,tol,linestyle=1,/continue,thick=2     
;fac=86400  & xtitle='Time (days)'
;  plot,(t-t(0))/fac,(deltaz) ,psym=-1,/ylog,yrange=[0.1,max(deltaz)],xrange=[0.1,max((t-t(0))/fac)*20],xtitle=xtitle,ytitle='Delta Z',thick=2,charthick=2,/xlog
;  oplot,(t-t(0))/fac,(deltaz),psym=-1,thick=2
;  x=(t-t(0))/fac & y=deltaz
;  a=poly_fit(alog10(x[1:n_elements(x)-1]),alog10(y[1:n_elements(x)-1]),1,/double,yfit=yfit) 
;  XX=findgen(100)*0.1+x(1) 
;  oplot,xx,10^(a(0)+alog10(XX)*a(1))
;  plots,0,tol & plots,max((t-t(0))/fac)*2,tol,linestyle=1,/continue,thick=2     
device,/close
goto,skiphisto
skipthis:



; Construct Marginalized Posterior Distribution from NESTPAR samples
; see 3.3 of Feroz & Hobson 2008, MNRAS, 384, 449
; log-weight of each sample is saved in array LP
; Nested sampling parameters are saved in nestpars
; Makes posterior copying each parameter proportonal to
; the weight.  Very low weighted values will not show up, while the
; highest weighted value will be duplicated the most times.

expandfactor=10000                 ;factor to duplicate
dupnum=round(exp(lp)*expandfactor) ;size of posterior with duplicates

a=where(dupnum gt 0)          ;non-zero weighted parameters
nonzeronestpars=nestpars(*,a) ;array of significant parameters
nonzeroLP=nestpars(a) ;array of significant weights

posterior=dblarr(nfit,total(round(exp(lp)*expandfactor))) 
LParray=dblarr(total(round(exp(lp)*expandfactor)))
;populate posteorior, duplicating 
count=0
for i=0,n_elements(a)-1 do begin & $
  for j=0,dupnum(a(i))-1 do begin & $
    posterior(*,count)=nonzeronestpars(*,i) & $
    LParray(count)=nonzeroLP(i) & $
    count=count+1  & $
  endfor & $
endfor



if (tofit(0) ne 0) then parameternamessimple=['Temp ','lnZ  ','lnC/O','pp_Rp','Cloud','CloudTop','CloudBott','CloudWave','Kzz','Haze','Kappa','gamma1','gamma2','beta','alpha','H2','He','H2O','CO2','CO','CH4','NH3','Na','K','Li','Rb','Cs','TiO','VO','FeH','PH3','H2S','HCN','C2H2']
if (tofit(0) eq 0) then parameternamessimple=['Temp ','lnZ  ','lnC/O','pp_Rp','Cloud','CloudTop','CloudBott','CloudWave','Kzz','H2','He','H2O','CO2','CO','CH4','NH3','Na','K','Li','Rb','Cs','TiO','VO','FeH','PH3','H2S','HCN','C2H2']

if (tofit(0) ne 0) then parameternames=['Temperature (K)','Metallicity ln(Z/Z$\tex_\odot$)','C/O Ratio ln(C/O)','R$\tex_{pl}$(R$\tex_{Jup}$,1 mbar)',$
'cloud opacity ln($\sigma$/$\sigma$$\tex_0$)','Haze Opacity ln($\sigma$/$\sigma$$\tex_0$-1)','Cloud Top (R$\tex_{pl}$)','Cloud Bottom (R$\tex_{pl}$)','Cloud Wavelength','K$\tex_{ZZ}$ (cm$\tex^2$ s)',$ 
'Infrared opacity $\kappa$$\tex_{IR}$ (cm$\tex^2$/g)','Optical/IR opacity ratio $\gamma$$\tex_1$','Optical/IR opacity ratio $\gamma$$\tex_2$','$\beta$','$\alpha$',$
'VMR ln(H$\tex_2$)','VMR ln(He)','VMR ln(H$\tex_2$O)','VMR ln(CO$\tex_2$)','VMR ln(CO)','VMR ln(CH$\tex_4$)','VMR ln(NH$\tex_3)','VMR ln(Na)','VMR ln(K)','VMR ln(Li)'$
,'VMR ln(Rb)','VMR ln(Cs)','VMR ln(TiO)','VMR ln(VO)','VMR ln(FeH)','VMR ln(PH$\tex_3)','VMR ln(H$\tex_2$S)','VMR ln(HCN)','VMR ln(C$\tex_2$H$\tex_2$)']

if (tofit(0) eq 0) then parameternames=['Temperature (K)','Metallicity ln(Z/Z$\tex_\odot$)','C/O Ratio ln(C/O)','R$\tex_{pl}$(R$\tex_{Jup}$,1 mbar)',$
'cloud opacity ln($\sigma$/$\sigma$$\tex_0$)','Haze Opacity ln($\sigma$/$\sigma$$\tex_0$-1)','Cloud Top (R$\tex_{pl}$)','Cloud Bottom (R$\tex_{pl}$)','Cloud Wavelength','K$\tex_{ZZ}$ (cm$\tex^2$ s)',$ 
'VMR ln(H$\tex_2$)','VMR ln(He)','VMR ln(H$\tex_2$O)','VMR ln(CO$\tex_2$)','VMR ln(CO)','VMR ln(CH$\tex_4$)','VMR ln(NH$\tex_3)','VMR ln(Na)','VMR ln(K)','VMR ln(Li)',$
'VMR ln(Rb)','VMR ln(Cs)','VMR ln(TiO)','VMR ln(VO)','VMR ln(FeH)','VMR ln(PH$\tex_3)','VMR ln(H$\tex_2$S)','VMR ln(HCN)','VMR ln(C$\tex_2$H$\tex_2$)']


  

MCMCpars=dblarr(nfit)
MCMCPCERRORlow=dblarr(nfit)
MCMCPCERRORhigh=dblarr(nfit)


burnndx=0

c=posterior;NESTPARS
p0=transpose(c(0,*))
if (nfit ge 2) then p1=transpose(c(1,*))
if (nfit ge 3) then p2=transpose(c(2,*))
if (nfit ge 4) then p3=transpose(c(3,*))
if (nfit ge 5) then p4=transpose(c(4,*))
if (nfit ge 6) then p5=transpose(c(5,*))
if (nfit ge 7) then p6=transpose(c(6,*))
if (nfit ge 8) then p7=transpose(c(7,*))
if (nfit ge 9) then p8=transpose(c(8,*))
if (nfit ge 10) then p9=transpose(c(9,*))
if (nfit ge 11) then p10=transpose(c(10,*))
if (nfit ge 12) then p11=transpose(c(11,*))
if (nfit ge 13) then p12=transpose(c(12,*))
if (nfit ge 14) then p13=transpose(c(13,*))
if (nfit ge 15) then p14=transpose(c(14,*))

if (nfit eq 1) then pall=[[p0]]
if (nfit eq 2) then pall=[[p0],[p1]]
if (nfit eq 3) then pall=[[p0],[p1],[p2]]
if (nfit eq 4) then pall=[[p0],[p1],[p2],[p3]]
if (nfit eq 5) then pall=[[p0],[p1],[p2],[p3],[p4]]
if (nfit eq 6) then pall=[[p0],[p1],[p2],[p3],[p4],[p5]]
if (nfit eq 7) then pall=[[p0],[p1],[p2],[p3],[p4],[p5],[p6]]
if (nfit eq 8) then pall=[[p0],[p1],[p2],[p3],[p4],[p5],[p6],[p7]]
if (nfit eq 9) then pall=[[p0],[p1],[p2],[p3],[p4],[p5],[p6],[p7],[p8]]
if (nfit eq 10) then pall=[[p0],[p1],[p2],[p3],[p4],[p5],[p6],[p7],[p8],[p9]]
if (nfit eq 11) then pall=[[p0],[p1],[p2],[p3],[p4],[p5],[p6],[p7],[p8],[p9],[p10]]
if (nfit eq 12) then pall=[[p0],[p1],[p2],[p3],[p4],[p5],[p6],[p7],[p8],[p9],[p10],p[11]]
if (nfit eq 13) then pall=[[p0],[p1],[p2],[p3],[p4],[p5],[p6],[p7],[p8],[p9],[p10],p[11],p[12]]
if (nfit eq 14) then pall=[[p0],[p1],[p2],[p3],[p4],[p5],[p6],[p7],[p8],[p9],[p10],p[11],p[12],p[13]]
if (nfit eq 15) then pall=[[p0],[p1],[p2],[p3],[p4],[p5],[p6],[p7],[p8],[p9],[p10],p[11],p[12],p[13],p[14]]


;parameternames=['Temperature (K)','Metallicity ln(Z/Z$\tex_\odot$)','Reference Radius ($\tex\DeltaR_{pl}$/R$\tex_{star}$)','cloud opacity ln($\sigma$/$\sigma$$\tex_0$)']


title='1-e^{-\Delta\chi^{2}}^{/2}  (%)'

;chi=reform(ALLCHAINSCHI2(burnndx:i,*),a)
allchains=posterior;NESTPARS

   levels   =[68,95,99]
   cc=floor((levels)/95.*256)
MCMCpars=dblarr(nfit)
MCMCPCERRORlow=dblarr(nfit)
MCMCPCERRORhigh=dblarr(nfit)

for p=0,nfit-1 do begin
    e=percentiles(allchains(p,*))-median(allchains(p,*)) 
    e2=percentiles(allchains(p,*),CONFLIMIT=0.95449)-median(allchains(p,*)) 
    e3=percentiles(allchains(p,*),CONFLIMIT=0.9973002)-median(allchains(p,*)) 
    medpar=median(allchains(p,*)) & $
    meanpar=median(allchains(p,*)) & $
  print,parameternames(tofit(p))+string(medpar)+' 1-sig:'+string(e(0))+string(e(1))+'  2-sig:'+string(e2(0))+string(e2(1))+'  3-sig:'+string(e3(0))+string(e3(1))
    MCMCpars(p)=medpar
    MCMCPCERRORlow(p)=e(0)
    MCMCPCERRORhigh(p)=e(1)

 endfor



;res=16
filenam='arc_monitor_density'
cgPS_Open, dirs+filenam+'.ps'
   cgDisplay
;device,filename=dirs+filenam+'.ps',/color,/landscape;,/encapsul
   pos = cgLayout([nfit,nfit],XGap=0.6,YGap=0.6)
   k=0
     xchrsz=0.001
     ychrsz=0.001
   FOR j=0,nfit-1 DO BEGIN
   FOR i=0,nfit-1 DO BEGIN
     if (j eq nfit-1) then xchrsz=0.3
     if (i eq 0) then ychrsz=0.3
     ;plot top corner histogram
     if (j eq 0 and i eq 0) then begin
        cghistoplot,pall(*,j),/fill,locations=locations,nbins=histobins,histdata=histdata,DATACOLORNAME='Dark Grey', POLYCOLOR='sky blue',xtitle=parameternames(tofit(j)),/frequenc,xthick=2,ythick=2,charthick=2,xcharsize=xchrsz,ycharsize=ychrsz,Position=pos[*,k],NoErase=j NE 0,xrange=[min(pall(*,0)),max(pall(*,0))],thick=2;,yrange=[0,max(histdata)]
        plots,MCMCpars(k),0 & plots,MCMCpars(k),max(histdata)*1.05,/continue,thick=2
        plots,MCMCpars(k)+MCMCPCERRORlow(k),0 & plots,MCMCpars(k)+MCMCPCERRORlow(k),max(histdata)*1.05,/continue,linestyle=2,thick=2
        plots,MCMCpars(k)+MCMCPCERRORhigh(k),0 & plots,MCMCpars(k)+MCMCPCERRORhigh(k),max(histdata)*1.05,/continue,linestyle=2,thick=2
     endif
     ;Plot diagnal histograms
     if (j eq i and j ne 0 and i ne 0) then begin 
       cghistoplot,pall(*,j),/fill,locations=locations,nbins=histobins,histdata=histdata,DATACOLORNAME='Dark Grey', POLYCOLOR='sky blue',xtitle=parameternames(tofit(j)),/frequenc,xthick=2,ythick=2,charthick=2,xcharsize=xchrsz,ycharsize=ychrsz,yticklen=0,yminor=1,ytickinterval=1,Position=pos[*,k],NoErase=j NE 0,xrange=[min(pall(*,i)),max(pall(*,i))],thick=2
       ;overplot 1-sigma ranges
       cgplot,[MCMCpars(i),MCMCpars(i)],[0,1],/overplot,thick=2,Position=pos[*,k], NoErase=j NE 0
       cgplot,[MCMCpars(i)+MCMCPCERRORlow(i),MCMCpars(i)+MCMCPCERRORlow(i)],[0,1],/overplot,thick=2,linestyle=2,Position=pos[*,k], NoErase=j NE 0
       cgplot,[MCMCpars(i)+MCMCPCERRORhigh(i),MCMCpars(i)+MCMCPCERRORhigh(i)],[0,1],/overplot,thick=2,linestyle=2,Position=pos[*,k], NoErase=j NE 0
     endif
     ;label Y-axis histogram axis on left side for all but top corner
     if (j eq i and j ne 0 and i ne 0) then cgaxis,yaxis=1,ytitle='Relative Frequency',ycharsize=0.25,xthick=2,ythick=2,charthick=2;,yrange=[0,max(histdata)]
     ;plot Distributions
     if (i lt j and j ne i) then begin
       red_band = pall(*,i)
       green_band = pall(*,j)
       xrange = [Min(red_band ), Max(red_band )]
       yrange = [Min(green_band ), Max(green_band )]
       xbinsize = (xrange(1)-xrange(0))/denres
       ybinsize = (yrange(1)-yrange(0))/denres
       h2d = HIST_2D(red_band, green_band, Min1=xrange[0], Max1=xrange[1], Bin1=xbinsize, $
                           Min2=yrange[0], Max2=yrange[1], Bin2=ybinsize) 
       ; Resize (with nearest-neighbor sampling) the 
       ; output array of bins to the original range of 
       ; intensities in the red and green bands, [0-255].
       ;h2d = CONGRID(h2d, MAX(red_band), MAX(green_band))
       ; Squash peaks by displaying logs of bins.
       ;h2d = BYTSCL(ALOG10(h2d > 1))
       maxDensity = Ceil(Max(h2d)/1e2) * 1e2

       ; Plot the 2D histogram.
       ct = COLORTABLE(0, /REVERSE)
       ;g0 = IMAGE(h2d, RGB_TABLE=ct, AXIS_STYLE=2, MARGIN=0.1, $
       
       cgLoadCT, 74,/reverse;33;65
          TVLCT, cgColor('white', /Triple), 0
          TVLCT, r, g, b, /Get
          palette = [ [r], [g], [b] ]     
          cgContour,LParray, pall(*,i), pall(*,j), NLevels=10, /Irregular, /Fill, /Outline,xtitle=parameternames(tofit(i)),ytitle=parameternames(tofit(j)),label=0,xthick=2,ythick=2,charthick=1,resolution=[res,res],xcharsize=xchrsz,ycharsize=ychrsz,Position=pos[*,k], NoErase=j NE 0,olevels=olevels,xrange=[min(pall(*,i)),max(pall(*,i))],/xstyle,cell_fill=1       
;          cgimage,h2d,XRange=xrange, YRange=yrange, /Axes, Palette=palette,  XTITLE=parameternames(tofit(i)), $
;          YTITLE=parameternames(tofit(j)),Position=pos[*,k],NoErase=j NE 0,charsize=ychrsz,/overplot
 
         cgimage,h2d,XRange=xrange, YRange=yrange, /Axes, Palette=palette,Position=pos[*,k],NoErase=j NE 0,charsize=0.001,/overplot


      ;cgContour,100-100*exp(-0.5*(chi-min(chi))), pall(*,i), pall(*,j), NLevels=10, /Irregular, /Fill, /Outline,xtitle=parameternames(tofit(i)),ytitle=parameternames(tofit(j)),label=0,xthick=2,ythick=2,charthick=1,resolution=[res,res],xcharsize=xchrsz,ycharsize=ychrsz,Position=pos[*,k], NoErase=j NE 0,olevels=olevels,xrange=[min(pall(*,i)),max(pall(*,i))],/xstyle,cell_fill=1
     
    endif
     xchrsz=0.001
     ychrsz=0.001
     k=k+1  
  ENDFOR
ENDFOR
;print,olevels
  ;cgLoadCT, 6,/brewer,/reverse;, /Brewer, /Reverse
   cgColorbar, /Vertical, Range=[0, maxDensity], NColors=254,charsize=1,charthick=2, Position=[0.96, 0.5, 0.98, 0.89],title='Density'
   cgPS_Close


  spawn,'ps2eps -B -C -f '+dirs+filenam+'.ps'
  spawn,'ps2pdf -dEPSCrop -dAutoRotatePages=/None -dAutoFilterColorImages=false -dColorImageFilter=/FlateEncode -dUseFlateCompression=true '+dirs+filenam+'.eps '+dirs+filenam+'.pdf'
  spawn,'pdftk temporary.pdf cat 1L output '+dirs+filenam+'.pdf'
;  spawn,'open '+dirs+filenam+'.ps'


!P.Multi = 0

;print,water abudnance if MdH and C/O is fit
a=where(tofit eq 2) 
if (a(0) ne -1) then begin
  if (tofit(0) eq 0) then h2o=exp(p1)/exp(p2) 
  if (tofit(0) eq 1) then h2o=exp(p0)/exp(p1)
  medpar=median(alog(h2O)) & $
  e=percentiles(alog(h2O))-median(alog(h2O)) 
  e2=percentiles(alog(h2O),CONFLIMI=0.95449)-median(alog(h2O)) 
  e3=percentiles(alog(h2O),CONFLIMI=0.9973)-median(alog(h2O)) 
  ev=percentiles(alog(h2O))
  e2v=percentiles(alog(h2O),CONFLIMI=0.95449)
  e3v=percentiles(alog(h2O),CONFLIMI=0.9973)
  print,''
  print,'H2O median value and error bars'
  print,'ln(H2O):'+string(medpar)+' 1-sig:'+string(e(0))+string(e(1))+'  2-sig:'+string(e2(0))+string(e2(1))+'  3-sig:'+string(e3(0))+string(e3(1))
  print,'H2O x Solar'
  print,'  (H2O):'+string(exp(medpar))+' 1-sig:'+string(-1*(-exp(ev(0))+exp(medpar)))+string(exp(ev(1))-exp(medpar))+'  2-sig:'+string(-1*(-exp(e2v(0))+exp(medpar)))+string(exp(e2v(1))-exp(medpar))+'  3-sig:'+string(-1*(-exp(e3v(0))+exp(medpar)))+string(exp(e3v(1))-exp(medpar))
  b=where((h2O) gt 67.0) & c=where((h2O) lt 67.) & print,string((1-double(n_elements(b))/double(n_elements(c)))*100),' % of distribution < 67x solar'
  b=where((h2O) gt 89.0) & c=where((h2O) lt 89.) & print,string((1-double(n_elements(b))/double(n_elements(c)))*100),' % of distribution < 89x solar'
  set_plot,'ps'                                                                                                                                       
 device,/landscape,filename=dirs+'H2Ohisto.ps'
  cghistoplot,alog10(h2o),/frequency,locations=locations,histdata=histdata,nbins=30,POLYCOLOR='sky blue',/fill,DATACOLORNAME='Dark Grey',xthick=4,ythick=4,charthick=4,xtitle='log$\tex_{10}$(H$\tex_{2}$O x Solar)'
 device,/close

skiphisto:
print,'Number of Points in Posterior: ',n_elements(LP)
print,'lnZ=',LOGZ,' +/-',LOGZ_ERR
;Nnest=n_elements(logz)
;M=100d0 ;Sample_number
;A=(M/(M+1))
;B=(M+1)/(M+2)
;subtotA=0
;subtotB=0
;subtotC=0
;first term eq 23
;for k=1,Nnest do begin & $
; subtotA=LP(k-1)*A^k & $
; for i=1,k do subtotB=subtotB+LP(i-1)*B^i & $
; subtotC=subtotC+subtotA*subtotB & $
; subtotB=0 & $
;endfor
;subtotD=0
;for i=1,Nnest do begin subtotD=subtotD+LP(i-1)*A^i
;sigz=sqrt(2/(M*(M+1))*subtotC-1/(M^2)*(subtotD)^2 ) ;problem is I
;have evidence saved, not L itself.

;print,water abudnance in freechemistry case
;a=where(tofit eq 13)
;solarH2Ovmr=0.000360745215794534;6E-4
;if (a(0) ne -1) then begin
;  h20=where(tofit eq 13)
;  b=where(exp(p3)/solarH2Ovmr ge 67.0) & c=where(exp(p3)/solarH2Ovmr lt 67) & print,string((1-double(n_elements(b))/double(n_elements(c)))*100),' % of distribution < 67x solar'
;  b=where(exp(p3)/solarH2Ovmr ge 89.0) & c=where(exp(p3)/solarH2Ovmr lt 89) & print,string((1-double(n_elements(b))/double(n_elements(c)))*100),' % of distribution < 89x solar'
;  H2O=exp(P3)   
;  CO2=exp(p4)
;  CO=exp(p5)
;  CH4=exp(p6)
;  O=(H2O+2*CO2+CO)
;  C=(CO2+CO+CH4)
;  set_plot,'ps'                                                                                                                                       
;  device,/landscape,filename=dirs+'COratio.ps'
;  cghistoplot,alog10(c/o),/frequency,locations=locations,histdata=histdata,nbins=30,POLYCOLOR='sky blue',/fill,DATACOLORNAME='Dark Grey',xthick=4,ythick=4,charthick=4,xtitle='log$\tex_{10}$(C/O ratio)'
;  e=percentiles(alog10(c/o))-median(alog10(c/o)) 
;  plots,median(alog10(c/o)),0 & plots,median(alog10(c/o)),max(histdata)*1.05,/continue,thick=7
;  plots,median(alog10(c/o))+e(0),0 & plots,median(alog10(c/o))+e(0),max(histdata)*1.05,/continue,linestyle=2,thick=7
;  plots,median(alog10(c/o))+e(1),0 & plots,median(alog10(c/o))+e(1),max(histdata)*1.05,/continue,linestyle=2,thick=7
;  print,'log10(C/O ratio): ',median(alog10(c/o)),e(0),e(1)
;  eln=percentiles(alog(c/o))-median(alog(c/o)) 
;  print,'ln(C/O ratio): ',median(alog(c/o)),eln(0),eln(1)
;endif

;device,/close



end



