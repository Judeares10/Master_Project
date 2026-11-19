pro arc_monitor,dirs
; Monitor MCMC convergence
; To be run in a separte IDL session to AIT
; Loads current MCMC chain, calculates Gelman-Rubin statistics
; Requires: Coyote Library http://www.idlcoyote.com/documents/programs.php
!p.font=0 
set_plot,'x'
restore,dirs+'MCMC.monitor.sav' ;allchains,allchainschi2,i,j,nchains,nthin,tofit

;Result = FILE_Search(dirs,'*.MPFIT.FitPars.sav')
;restore,result(0)

nfit=n_elements(tofit) 
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

;===
;running gelman 
gelmanrubinall=gelmanrubin
for k=burnndx+10,i-1 do begin & $
  converged = exofast_gelmanrubin(allchains[0:nfit-1,burnndx:k,*],gelmanrubin,tz,angular=gelmanangular)  & $
  gelmanrubinall=[[gelmanrubinall],[gelmanrubin]] & $
endfor
  window,5,xsize=800,ysize=500
  wset,5
  plot,gelmanrubinall(0,*),/xlog,xrange=[10,n_elements(allchains(0,*,0))],yrange=[0.8,2.5],/xstyle,ytitle='Gelman-Rubin Statistic',xtitle='MCMC Steps'
  if (nfit ge 2) then oplot,gelmanrubinall(1,*),color=cgcolor('red', !D.Table_Size-2)
  if (nfit ge 3) then oplot,gelmanrubinall(2,*),color=cgcolor('green', !D.Table_Size-2)
  if (nfit ge 4) then oplot,gelmanrubinall(3,*),color=cgcolor('dodger blue', !D.Table_Size-2)
  if (nfit ge 5) then oplot,gelmanrubinall(4,*),color=cgcolor('dodger blue', !D.Table_Size-2)
  if (nfit ge 6) then oplot,gelmanrubinall(5,*),color=cgcolor('dodger blue', !D.Table_Size-2)
  if (nfit ge 7) then oplot,gelmanrubinall(6,*),color=cgcolor('dodger blue', !D.Table_Size-2)
  if (nfit ge 8) then oplot,gelmanrubinall(7,*),color=cgcolor('dodger blue', !D.Table_Size-2)
  plots,10,1.01 & plots,n_elements(allchains(0,*,0)),1.01,/continue,linestyle=1 
  xyouts,20,0.9,parameternamessimple(tofit(0))
  if (nfit ge 2) then   xyouts,20,0.8,parameternamessimple(tofit(1)),color=cgcolor('red', !D.Table_Size-2)
  if (nfit ge 3) then   xyouts,20,0.7,parameternamessimple(tofit(2)),color=cgcolor('green', !D.Table_Size-2)
  if (nfit ge 4) then   xyouts,20,0.6,parameternamessimple(tofit(3)),color=cgcolor('dodger blue', !D.Table_Size-2)
  if (nfit ge 5) then   xyouts,20,0.5,parameternamessimple(tofit(4)),color=cgcolor('purple', !D.Table_Size-2)
  if (nfit ge 6) then   xyouts,20,0.5,parameternamessimple(tofit(5)),color=cgcolor('purple', !D.Table_Size-2)
  if (nfit ge 7) then   xyouts,20,0.5,parameternamessimple(tofit(6)),color=cgcolor('purple', !D.Table_Size-2)
  if (nfit ge 8) then   xyouts,20,0.5,parameternamessimple(tofit(7)),color=cgcolor('purple', !D.Table_Size-2)

;====


MCMCpars=dblarr(nfit)
MCMCPCERRORlow=dblarr(nfit)
MCMCPCERRORhigh=dblarr(nfit)

for p=0,nfit-1 do begin
   e=percentiles(allchains(p,burnndx:i,*))-mean(allchains(p,burnndx:i,*)) 
    e3=percentiles(allchains(p,burnndx:i,*),CONFLIMIT=0.9973002)-mean(allchains(p,burnndx:i,*)) 
    medpar=median(allchains(p,burnndx:i,*)) & $

    MCMCpars(p)=medpar
    MCMCPCERRORlow(p)=e(0)
    MCMCPCERRORhigh(p)=e(1)
    
endfor
MCMCPCERROR=(abs(MCMCPCERRORhigh)+abs(MCMCPCERRORlow))/2.
save,filename=dirs+'MCMC.monitor.FitPars.sav',MCMCpars,MCMCPCERRORlow,MCMCPCERRORhigh,MCMCPCERROR

for p=0,nfit-1 do begin
  window,p,xsize=800,ysize=500
  wset,p
  plot,allchains(p,burnndx:i,1),xrange=[1,i],/xstyle,yrange=[MCMCpars(p)+4*MCMCPCERRORlow(p),MCMCpars(p)+4*MCMCPCERRORhigh(p)],/ystyle,title='Parameter '+string(p)
  for k=0,nchains-1 do oplot,allchains(p,burnndx:i,k),color=123123*k
  plots,1,MCMCpars(p) & plots,i,MCMCpars(p),/continue,thick=2
  plots,1,MCMCpars(p)+MCMCPCERRORlow(p) & plots,i,MCMCpars(p)+MCMCPCERRORhigh(p),/continue,linestyle=1,thick=2
  plots,1,MCMCpars(p)-MCMCPCERRORlow(p) & plots,i,MCMCpars(p)-MCMCPCERRORhigh(p),/continue,linestyle=1,thick=2

    e=percentiles(allchains(p,burnndx:i,*))-mean(allchains(p,burnndx:i,*)) 
    e2=percentiles(allchains(p,burnndx:i,*),CONFLIMIT=0.95449)-mean(allchains(p,burnndx:i,*)) 
    e3=percentiles(allchains(p,burnndx:i,*),CONFLIMIT=0.9973002)-mean(allchains(p,burnndx:i,*)) 
    medpar=median(allchains(p,burnndx:i,*)) & $

    MCMCpars(p)=medpar
    MCMCPCERRORlow(p)=e(0)
    MCMCPCERRORhigh(p)=e(1)
  plots,1,medpar & plots,i,medpar,/continue,thick=3,color=1000
  print,parameternamessimple(p)+string(medpar)+' 1-sig:'+string(e(0))+string(e(1))+'  2-sig:'+string(e2(0))+string(e2(1))+'  3-sig:'+string(e3(0))+string(e3(1))
  plots,1,medpar+e(0) & plots,i,medpar+e(0),/continue,linestyle=2,thick=3,color=1000
  plots,1,medpar+e(1) & plots,i,medpar+e(1),/continue,linestyle=2,thick=3,color=1000
  ;cghistoplot,allchains(p,0:i,*),binsize=(e(1)/6.),/fill,DATACOLORNAME='black', POLYCOLOR='sky blue',/frequency,/rotate
endfor
;window,nfit
;wset,nfit
;plot,allchains(1,0:i,1),allchains(3,0:i,1),psym=1,xrange=[1,2.5],yrange=[0,4]
;  for k=0,nchains-1 do oplot,allchains(1,0:i,k),allchains(3,0:i,k),psym=1

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
   res=25
!P.Multi = 0
;!P.Multi = [0, 3, 3, 0, 1]
set_plot,'ps'
plotsym,0,0.25,/fill
filenam='ait_monitor_distributions'
device,filename=dirs+filenam+'.ps',/color,/landscape;,/encapsul
;Metalicty vs Clouds ;25


for k=0,nfit-1 do begin
for l=0,nfit-1 do begin
if ((k ne l) and (k gt l)) then begin
  cgLoadCT, 6,/brewer,/reverse;,/reverse;, /Brewer, /Reverse
    cgContour,100-100*exp(-0.5*(chi-min(chi))), pall(*,k), pall(*,l), NLevels=10, /Irregular, /Fill, /Outline,xtitle=parameternames(tofit(k)),ytitle=parameternames(tofit(l)),label=0,xthick=4,ythick=4,charthick=4,resolution=[res,res]
    ;cgContour,100-100*exp(-0.5*(chi-min(chi))), pall(*,k), pall(*,l), levels=levels, /Irregular, /Outline,label=1,/overplot,thick=3,resolution=[res,res]

    ;cgPlot, pall(*,k), pall(*,l),PSYM=8,/overplot,symcolor=cgcolor('dodger blue', !D.Table_Size-2)
;    oploterror,MPFITpars(1),MPFITpars(3),MPFITPCERROR(1),MPFITPCERROR(3),thick=2,errcolor=cgcolor('green', !D.Table_Size-2)
  ;  oploterror,MCMCpars(k),MCMCpars(l),MCMCPCERRORlow(k),MCMCPCERRORlow(l),/LOBAR,errcolor=cgcolor('blue', !D.Table_Size-2),errthick=2
  ;  oploterror,MCMCpars(k),MCMCpars(l),MCMCPCERRORhigh(k),MCMCPCERRORhigh(l),/HIBAR,errcolor=cgcolor('blue', !D.Table_Size-2),errthick=2
 endif

endfor
endfor
for k=0,nfit-1 do begin
cghistoplot,pall(*,k),/fill,locations=locations,histdata=histdata,DATACOLORNAME='Dark Grey', POLYCOLOR='sky blue',xtitle=parameternames(tofit(k)),/frequenc,xthick=4,ythick=4,charthick=4
;a=where((pall(*,k) gt MCMCpars(k)+MCMCPCERRORlow(k)) and (pall(*,k)
;lt MCMCpars(k)+MCMCPCERRORhigh(k)))
;cghistoplot,pall(*,k),locations=locations(a),POLYCOLOR='blue',/oplot,/fill,/frequenc,maxinput=MCMCpars(k)+MCMCPCERRORhigh(k)
plots,MCMCpars(k),0 & plots,MCMCpars(k),max(histdata)*1.05,/continue,thick=7
plots,MCMCpars(k)+MCMCPCERRORlow(k),0 & plots,MCMCpars(k)+MCMCPCERRORlow(k),max(histdata)*1.05,/continue,linestyle=2,thick=7
plots,MCMCpars(k)+MCMCPCERRORhigh(k),0 & plots,MCMCpars(k)+MCMCPCERRORhigh(k),max(histdata)*1.05,/continue,linestyle=2,thick=7

endfor


device,/close

  ;spawn,'ps2eps -B -C -f '+filenam+'.ps'
  ;spawn,'ps2pdf -dEPSCrop -dAutoRotatePages=/None -dAutoFilterColorImages=false -dColorImageFilter=/FlateEncode -dUseFlateCompression=true '+filenam+'.eps temporary.pdf'
  ;spawn,'pdftk temporary.pdf cat 1L output '+filenam+'.pdf'
  ;spawn,'open '+filenam+'.ps'


!P.Multi = 0
end



