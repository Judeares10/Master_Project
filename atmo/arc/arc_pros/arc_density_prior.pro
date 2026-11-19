pro arc_density_prior,prior,dirs
; Monitor MCMC convergence
; To be run in a separte IDL session to AIT
; Loads current MCMC chain, calculates Gelman-Rubin statistics
; Requires: Coyote Library http://www.idlcoyote.com/documents/programs.php
!p.font=0 
loadct,0
set_plot,'x'
restore,dirs+'MCMC.monitor.sav' ;allchains,allchainschi2,i,j,nchains,nthin,tofit

i=1.4E4
  ;Chains may not find CHI^2 minimum, go rouge, ect.  Find them and destroy them
   fin=i
   strt=0 & if (i gt 4000) then strt=2000
   ;calculate average chi^2 of each chain to compare
   chiave=dblarr(nchains) & for w=0,nchains-1 do chiave(w)=mean(allchainschi2(strt:fin,w)-min(allchainschi2(0:i,*)))
   sigma_cut=2  ;sigma cut on average chi^2 of the chains
   print,'min CHI^2 = ',min(allchainschi2(0:i,*))
   print,'sigma cut:',sigma_cut
   resistant_mean,chiave,sigma_cut,mean,sigma,num_rejected,goodvec=good & print,'rejected ',num_rejected,'       good:',good  
   print,'Good chains have ave CHI^2 values within delta Chi^2 :',sigma
   bad=where((chiave-mean)/(sigma_cut*sigma) gt sigma_cut)  
   print,'Bad chains have ave delta CHI^2 values of :',(chiave(bad)-mean)

  ;just keep good chains
   include=good ;[0,1,2,3,4,5,6,7,8,9,10,11] ;ALL
  ;include=[0,1,2,3,5,7,8,9,10,11] & i=25000; ; M3 frechemm with CO2
 
  trimchain=nchains-n_elements(include)
   allchains=allchains[*,*,include]
   allchainschi2=allchainschi2[*,include]
   j=j-trimchain
   nchains=nchains-trimchain
  

colr=[0,1000,100013,1000134,10001345,100013454,1000134542,10001345422,100013454228,1000134542285,10001345422853,100013454228535]
strt=0
if (i gt 1000) then strt=1000
fin=i
plot,allchainschi2(strt:fin,0)-min(allchainschi2),/ylog,title='chi2'
for w=1,n_elements(include)-1 do oplot,allchainschi2(strt:fin,w)-min(allchainschi2),color=colr(w)

;par=0 
;plot,allchains(par,*,0)
;for w=1,n_elements(include)-1 do oplot,allchains(par,*,w),color=colr(w)

;H2O
;plot,exp(allchains(1,*,0))/(exp(allchains(2,*,0))),/ylog
;for w=1,n_elements(include)-1 do oplot,exp(allchains(1,*,w))/(exp(allchains(2,*,w))),color=colr(w)

;C/O
;plot,allchains(2,*,0)
;for w=1,n_elements(include)-1 do oplot,allchains(2,*,w),color=colr(w)


;Number of bins in Histogram
histobins=25;100;28 ;45
res=12;100;50;25
denres=90;30;45;90

nfit=n_elements(tofit) 
if (tofit(0) ne 0) then parameternamessimple=['Temp ','lnZ  ','lnC/O','pp_Rp','Cloud','CloudTop','CloudBott','CloudWave','Kzz','Haze','Kappa','gamma1','gamma2','beta','alpha','H2','He','H2O','CO2','CO','CH4','NH3','Na','K','Li','Rb','Cs','TiO','VO','FeH','PH3','H2S','HCN','C2H2']
if (tofit(0) eq 0) then parameternamessimple=['Temp ','lnZ  ','lnC/O','pp_Rp','Cloud','CloudTop','CloudBott','CloudWave','Kzz','H2','He','H2O','CO2','CO','CH4','NH3','Na','K','Li','Rb','Cs','TiO','VO','FeH','PH3','H2S','HCN','C2H2']

if (tofit(0) ne 0) then parameternames=['Temperature (K)','Metallicity ln(Z/Z$\tex_\odot$)','C/O Ratio ln(C/O)','R$\tex_{pl}$(R$\tex_{Jup}$,1 mbar)',$
'cloud opacity ln($\sigma$/$\sigma$$\tex_0$)','Haze Opacity ln($\sigma$/$\sigma$$\tex_0$-1)','Cloud Top (R$\tex_{pl}$)','Cloud Bottom (R$\tex_{pl}$)','Cloud Wavelength','K$\tex_{ZZ}$ (cm$\tex^2$/s)',$ 
'Infrared opacity $\kappa$$\tex_{IR}$ (cm$\tex^2$/g)','Optical/IR opacity ratio $\gamma$$\tex_1$','Optical/IR opacity ratio $\gamma$$\tex_2$','$\beta$','$\alpha$',$
'VMR ln(H$\tex_2$)','VMR ln(He)','VMR ln(H$\tex_2$O)','VMR ln(CO$\tex_2$)','VMR ln(CO)','VMR ln(CH$\tex_4$)','VMR ln(NH$\tex_3$)','VMR ln(Na)','VMR ln(K)','VMR ln(Li)'$
,'VMR ln(Rb)','VMR ln(Cs)','VMR ln(TiO)','VMR ln(VO)','VMR ln(FeH)','VMR ln(PH$\tex_3)','VMR ln(H$\tex_2$S)','VMR ln(HCN)','VMR ln(C$\tex_2$H$\tex_2$)']

if (tofit(0) eq 0) then parameternames=['Temperature (K)','Metallicity ln(Z/Z$\tex_\odot$)','C/O Ratio ln(C/O)','R$\tex_{pl}$(R$\tex_{Jup}$,1 mbar)',$
'cloud opacity ln($\sigma$/$\sigma$$\tex_0$)','Haze Opacity ln($\sigma$/$\sigma$$\tex_0$-1)','Cloud Top (R$\tex_{pl}$)','Cloud Bottom (R$\tex_{pl}$)','Cloud Wavelength','K$\tex_{ZZ}$ (cm$\tex^2$/s)',$ 
'VMR ln(H$\tex_2$)','VMR ln(He)','VMR ln(H$\tex_2$O)','VMR ln(CO$\tex_2$)','VMR ln(CO)','VMR ln(CH$\tex_4$)','VMR ln(NH$\tex_3$)','VMR ln(Na)','VMR ln(K)','VMR ln(Li)',$
'VMR ln(Rb)','VMR ln(Cs)','VMR ln(TiO)','VMR ln(VO)','VMR ln(FeH)','VMR ln(PH$\tex_3)','VMR ln(H$\tex_2$S)','VMR ln(HCN)','VMR ln(C$\tex_2$H$\tex_2$)']



      ;; discard the burn-in (the first point at which all chains
      ;; have had at least one chi^2 lower than the median chi^2)
      medchi2 = median(allchainschi2[0:i,*])
      burnndx = 0L
      for j=0L, nchains-1 do begin
         tmpndx = (where(allchainschi2[0:i,j] lt medchi2))(0)
         if tmpndx gt burnndx then burnndx = tmpndx
      endfor
      ;; allows Gelman-Rubin calculation if one chain is being problematic
      burnndx = burnndx < (i-3) 
      
      ;; calculate the Gelman-Rubin statistic (remove burn-in)
      ;converged = exofast_gelmanrubin(pars[0:nfit-1,burnndx:i,*],$
      ;                                gelmanrubin,tz,angular=gelmanangular)

  converged = exofast_gelmanrubin(allchains[0:nfit-1,burnndx:i,*],gelmanrubin,tz,angular=gelmanangular) 
  print,'Burn index:',burnndx
  print,'Gelman-Rubin < 1.01 => MCMC is Converged, and Indep. Draws > 1000' 
  print,'Gelman-Rubin      =',gelmanrubin 
  print,'Num. Indep. Draws =',tz
  ;; estimate the number of steps it will take until it is converged
    mintz = 0d0 & $
    tzsteps = 0L & $
    mintz = [mintz,min(tz)] & $
    tzsteps = [tzsteps,i] & $
    ntz = n_elements(mintz) & $
    coeffs = poly_fit(tzsteps[0:ntz-1],mintz[0:ntz-1],1,yfit=fit) & $
    stepstoconvergence = (2000 - coeffs[0])/coeffs[1] & $
    print,'MCMC of '+string(j,format='(I3)'),' chains with '+string(i+1,format='(I9)')+' steps run out of ',N_ELEMENTS(allchainschi2(*,0))
    print,'Steps to Converge =   ',stepstoconvergence/j,' for each chain' & $ 
    print,'Taken ',i,' steps with each chain, which is ',(i*j)/stepstoconvergence*100,' % of needed for G-R<1.01'
    print,'Taken ',i,' steps with each chain, which is ',((min(tz))/1000)*100,' % of needed for 1000 draws'



MCMCpars=dblarr(nfit)
MCMCPCERRORlow=dblarr(nfit)
MCMCPCERRORhigh=dblarr(nfit)

for p=0,nfit-1 do begin
   e=percentiles(allchains(p,burnndx:i,*))-median(allchains(p,burnndx:i,*)) 
    e3=percentiles(allchains(p,burnndx:i,*),CONFLIMIT=0.9973002)-median(allchains(p,burnndx:i,*)) 
    medpar=median(allchains(p,burnndx:i,*)) & $

    MCMCpars(p)=medpar
    MCMCPCERRORlow(p)=e(0)
    MCMCPCERRORhigh(p)=e(1)
    
endfor
MCMCPCERROR=(abs(MCMCPCERRORhigh)+abs(MCMCPCERRORlow))/2.


b=allchains(*,burnndx:i,*)
a=(n_elements(b(0,*,0)))*1D0*nchains & a=long(a) ;chain length
c=reform(b,nfit,a)

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

chi=reform(ALLCHAINSCHI2(burnndx:i,*),a)

   levels   =[68,95,99]
   cc=floor((levels)/95.*256)

for p=0,nfit-1 do begin
    e=percentiles(allchains(p,burnndx:i,*))-median(allchains(p,burnndx:i,*)) 
    e2=percentiles(allchains(p,burnndx:i,*),CONFLIMIT=0.95449)-median(allchains(p,burnndx:i,*)) 
    e3=percentiles(allchains(p,burnndx:i,*),CONFLIMIT=0.9973002)-median(allchains(p,burnndx:i,*)) 
    medpar=median(allchains(p,burnndx:i,*)) & $
    meanpar=median(allchains(p,burnndx:i,*)) & $
  ;print,parameternames(tofit(p))+string(medpar)+' 1-sig:'+string(e(0))+string(e(1))+'  2-sig:'+string(e2(0))+string(e2(1))+'  3-sig:'+string(e3(0))+string(e3(1))
 endfor



;res=16
filenam='ait_monitor_density'
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
          cgContour,100-100*exp(-0.5*(chi-min(chi))), pall(*,i), pall(*,j), NLevels=10, /Irregular, /Fill, /Outline,xtitle=parameternames(tofit(i)),ytitle=parameternames(tofit(j)),label=0,xthick=2,ythick=2,charthick=1,resolution=[res,res],xcharsize=xchrsz,ycharsize=ychrsz,Position=pos[*,k], NoErase=j NE 0,olevels=olevels,xrange=[min(pall(*,i)),max(pall(*,i))],/xstyle,cell_fill=1       
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
  spawn,'ps2pdf -dEPSCrop -dAutoRotatePages=/None -dAutoFilterColorImages=false -dColorImageFilter=/FlateEncode -dUseFlateCompression=true '+dirs+filenam+'.eps temporary.pdf'
  spawn,'pdftk temporary.pdf cat 1L output '+dirs+filenam+'.pdf'
;  spawn,'open '+dirs+filenam+'.ps'


!P.Multi = 0

;Equilbrium Chemistry
;get water abudnance if MdH and C/O is fit 
a=where(tofit eq 2) 
if (a(0) ne -1) then begin
  
  ;isothermal
  if (tofit(0) eq 0) then h2o=exp(p1)/exp(p2) 
  if (tofit(0) eq 0) then   lowerCOindex=where(exp(p2) gt prior)
  ;non-isothermal
  if (tofit(0) eq 1) then h2o=exp(p0)/exp(p1)
  if (tofit(0) eq 1) then   lowerCOindex=where(exp(p1) gt prior)
  print,'-----'
  print,n_elements(lowerCOindex),' have the acceptable prior ',n_elements(lowerCOindex)/n_elements(p0)*100,'% of the chain'
  h2oOrig=h2O
  h2o=h2o(lowerCOindex)
  medpar=median(alog(h2O)) & $
  e=percentiles(alog(h2O))-median(alog(h2O)) 
  e2=percentiles(alog(h2O),CONFLIMI=0.95449)-median(alog(h2O)) 
  e3=percentiles(alog(h2O),CONFLIMI=0.9973)-median(alog(h2O)) 
  ev=percentiles(alog(h2O))
  e2v=percentiles(alog(h2O),CONFLIMI=0.95449)
  e3v=percentiles(alog(h2O),CONFLIMI=0.9973)
  print,''
  print,'H2O median value and error bars, with prior'
  print,'ln(H2O):'+string(medpar)+' 1-sig:'+string(e(0))+string(e(1))+'  2-sig:'+string(e2(0))+string(e2(1))+'  3-sig:'+string(e3(0))+string(e3(1))
  print,'H2O x Solar'
  print,'  (H2O):'+string(exp(medpar))+' 1-sig:'+string(-1*(-exp(ev(0))+exp(medpar)))+string(exp(ev(1))-exp(medpar))+'  2-sig:'+string(-1*(-exp(e2v(0))+exp(medpar)))+string(exp(e2v(1))-exp(medpar))+'  3-sig:'+string(-1*(-exp(e3v(0))+exp(medpar)))+string(exp(e3v(1))-exp(medpar))
  b=where((h2O) gt 67.0) & c=where((h2O) lt 67.) & print,string((1-double(n_elements(b))/double(n_elements(c)))*100),' % of distribution < 67x solar'
  b=where((h2O) gt 89.0) & c=where((h2O) lt 89.) & print,string((1-double(n_elements(b))/double(n_elements(c)))*100),' % of distribution < 89x solar'
  set_plot,'ps'                                                                                                                                       
 device,/landscape,filename=dirs+'H2Ohisto.ps'
  cghistoplot,alog10(h2o),/frequency,locations=locations,histdata=histdata,nbins=30,POLYCOLOR='sky blue',/fill,DATACOLORNAME='Dark Grey',xthick=4,ythick=4,charthick=4,xtitle='log$\tex_{10}$(H$\tex_{2}$O x Solar)'
 device,/close
 if (tofit(0) eq 0) then coratio=(p2)
 if (tofit(0) eq 1) then coratio=(p1)
 goto,skipnext
endif

;Free chemsitry
;print,water abudnance in isothermal-cloudfree-freechemistry case
a=where(tofit eq 6)
  if (tofit(0) eq 0) then a=where(tofit eq 8)
solarH2Ovmr=0.000360745215794534;6E-4
if (a(0) ne -1) then begin
  print,'tofit =',tofit
  a=where(tofit eq 3) 
  ;cloud-free case
  if (a(0) eq -1) then begin
     print,'**cloud-free free-chem-case**'
    H2O=exp(P2)  & h2oOrig=h2O
    CO2=exp(p3)
    CO=exp(p4)
    CH4=exp(p5)
    O=(H2O+2*CO2+CO)
    C=(CO2+CO+CH4)
   end
  ;cloudy-case
  if (a(0) ne -1) then begin
    print,'**cloudy free-chem-case**'
    print,'H2O=p3'
    H2O=exp(P3) & h2oOrig=h2O
    CO2=exp(p4)
    CO=exp(p5)
    CH4=exp(p6)
    O=(H2O+2*CO2+CO)
    C=(CO2+CO+CH4)
 endif 
 ;T-p-cloud-free free-chem-case
  a=where(tofit eq 13) 
  if (a(0) ne -1) then begin
    print,'**T-p cloudy free-chem case**'
    print,'H2O=p4'
    H2O=exp(P4) & h2oOrig=h2O
    CO2=exp(p5)
    CO=exp(p6)
    CH4=exp(p7)
    O=(H2O+2*CO2+CO)
    C=(CO2+CO+CH4)
 endif 
 ;T-p-cloudy-free-chem-case
  ;a=where(tofit eq 13) 
  ;if (a(0) ne -1) then begin
  ;  print,'**T-p cloudy free-chem case**'
  ;  print,'H2O=p5'
  ;  H2O=exp(P5) & h2oOrig=h2O
  ;  CO2=exp(p6)
  ;  CO=exp(p7)
  ;  CH4=exp(p8)
  ;  O=(H2O+2*CO2+CO)
  ;  C=(CO2+CO+CH4)
 ;endif 


  coratio=alog(c/o)
  lowerCOindex=where(c/o gt prior)
  print,'-----'
  print,n_elements(lowerCOindex),' have the acceptable prior ',n_elements(lowerCOindex)/n_elements(p0),'% of the chain'

  h2o=h2o(lowerCOindex)
  medpar=median(alog(h2O)) & $
  e=percentiles(alog(h2O))-median(alog(h2O)) 
  e2=percentiles(alog(h2O),CONFLIMI=0.95449)-median(alog(h2O)) 
  e3=percentiles(alog(h2O),CONFLIMI=0.9973)-median(alog(h2O)) 
  ev=percentiles(alog(h2O))
  e2v=percentiles(alog(h2O),CONFLIMI=0.95449)
  e3v=percentiles(alog(h2O),CONFLIMI=0.9973)
  print,''
  print,'H2O median value and error bars, with prior'
  print,'VMRln(H2O):'+string(medpar)+' 1-sig:'+string(e(0))+string(e(1))+'  2-sig:'+string(e2(0))+string(e2(1))+'  3-sig:'+string(e3(0))+string(e3(1))

  b=where(H2O/solarH2Ovmr ge 67.0) & c=where(H2O/solarH2Ovmr lt 67) & print,string((1-double(n_elements(b))/double(n_elements(c)))*100),' % of distribution < 67x solar'
  b=where(H2O/solarH2Ovmr ge 89.0) & c=where(H2O/solarH2Ovmr lt 89) & print,string((1-double(n_elements(b))/double(n_elements(c)))*100),' % of distribution < 89x solar'

  set_plot,'ps'                                                                                                                                       
  device,/landscape,filename=dirs+'COratio.ps'
  cghistoplot,alog10(c/o),/frequency,locations=locations,histdata=histdata,nbins=30,POLYCOLOR='sky blue',/fill,DATACOLORNAME='Dark Grey',xthick=4,ythick=4,charthick=4,xtitle='log$\tex_{10}$(C/O ratio)'
  e=percentiles(alog10(c/o))-median(alog10(c/o)) 
  plots,median(alog10(c/o)),0 & plots,median(alog10(c/o)),max(histdata)*1.05,/continue,thick=7
  plots,median(alog10(c/o))+e(0),0 & plots,median(alog10(c/o))+e(0),max(histdata)*1.05,/continue,linestyle=2,thick=7
  plots,median(alog10(c/o))+e(1),0 & plots,median(alog10(c/o))+e(1),max(histdata)*1.05,/continue,linestyle=2,thick=7
  print,'log10(C/O ratio): ',median(alog10(c/o)),e(0),e(1)
  eln=percentiles(alog(c/o))-median(alog(c/o)) 
  print,'ln(C/O ratio): ',median(alog(c/o)),eln(0),eln(1)
endif
skipnext:
;device,/close

;================================================================================
  h2o=h2oOrig
c=where(exp(coratio) gt prior) 
print,'Minimum Chi^2 with Prior:',min(chi(c))

colimit=findgen(104)*(0.11220324)-11.803794
h2oT=dblarr(104)
h2oM=dblarr(104)
h2o1=dblarr(2,104)
h2o2=dblarr(2,104)
h2o3=dblarr(2,104)

for i=10,100 do begin & $
  c=where(coratio gt colimit(i)) & $
    medpar=median(alog(h2O(c))) & $
    medparT=median(p0(c)) & $
    e=percentiles(alog(h2O(c)))& $
    e2=percentiles(alog(h2O(c)),CONFLIMI=0.95449) & $
    e3=percentiles(alog(h2O(c)),CONFLIMI=0.9973) & $
  h2oM(i)=medpar & $
  h2oT(i)=medparT & $
  h2o1(*,i)=e & $
  h2o2(*,i)=e2 & $
  h2o3(*,i)=e3 & $
endfor
;plot,colimit,h2om
;plot,exp(colimit),exp(h2om),/xlog

;
sm=6
loadct,0
set_plot,'ps'
filenam='COprior_H2O'
device, /color, /landscape,/encapsul, /helvetica,filename=dirs+'COprior_H2O.ps'
plot,exp(colimit(0:93)),smooth(exp(h2om),sm),/xlog,/ylog,yrange=[1,400],xrange=[3E-5,0.56],/xstyle,xtitle='C/O Ratio Lower Limit Prior Value',ytitle=textoidl('H_2O Median VMR (x Solar) & Upper Limits'),/ystyle,thick=5, position=[0.23,0.4,0.95,0.95], xthick=5, ythick=5, charthick=5, charsize=1.2
oplot,exp(colimit(0:93)),smooth(exp(h2om),sm),thick=5
oplot,exp(colimit(0:93)),smooth(exp(h2o1(1,*)),sm),thick=2,linestyle=2
oplot,exp(colimit(0:93)),smooth(exp(h2o2(1,*)),sm),thick=2,linestyle=3
;oplot,exp(colimit(0:93)),smooth(exp(h2o3(1,*)),sm),thick=2,linestyle=4
device,/close
  spawn,'ps2eps -B -C -f '+dirs+filenam+'.ps'
  spawn,'ps2pdf -dEPSCrop -dAutoRotatePages=/None -dAutoFilterColorImages=false -dColorImageFilter=/FlateEncode -dUseFlateCompression=true '+dirs+filenam+'.eps temporary.pdf'
  spawn,'pdftk temporary.pdf cat 1L output '+dirs+filenam+'.pdf'

;temperature vs C/O ratio prior
set_plot,'ps'
device, /color, /landscape,/encapsul, /helvetica,filename=dirs+'T_H2O.ps'
filenam='T_H2O'

device,/landscape
x=exp(colimit)
a=where(x lt 0.3)
plot,x(a),[smooth(h2oT,5)],xrange=[3E-5,0.56] ,/xlog,/xstyle,yrange=[500,1000],xtitle='C/O ratio lower-limit Prior Value',ytitle='Median Temperature (K)',thick=5, position=[0.23,0.4,0.95,0.95], xthick=5, ythick=5, charthick=5, charsize=1.2
oplot,x(a),[smooth(h2oT,5)],thick=5
device,/close
  spawn,'ps2eps -B -C -f '+dirs+filenam+'.ps'
  spawn,'ps2pdf -dEPSCrop -dAutoRotatePages=/None -dAutoFilterColorImages=false -dColorImageFilter=/FlateEncode -dUseFlateCompression=true '+dirs+filenam+'.eps '+dirs+filenam+'.pdf'
  spawn,'pdftk temporary.pdf cat 1L output '+dirs+filenam+'.pdf'





end



