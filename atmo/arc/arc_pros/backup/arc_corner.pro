pro arc_corner,dirs
; Monitor MCMC convergence
; To be run in a separte IDL session to AIT
; Loads current MCMC chain, calculates Gelman-Rubin statistics
; Requires: Coyote Library http://www.idlcoyote.com/documents/programs.php
!p.font=0 
set_plot,'x'
restore,dirs+'MCMC.monitor.sav' ;allchains,allchainschi2,i,j,nchains,nthin,tofit

;Result = FILE_Search(dirs,'*.MPFIT.FitPars.sav')
;restore,result(0)

histobins=26
res=10

nfit=n_elements(tofit) 
parameternamessimple=['Temp ','lnZ  ','lnC/O','pp_Rp','Cloud','Haze','H2','He','H2O','CO2','CO','CH4','NH3','Na','K','Li','Rb','Cs','TiO','VO','FeH','PH3','H2S','HCN','C2H2']
;parameternames=['Temperature (K)','Metallicity ln(Z/Z$\tex_\odot$)','C/O Ratio ln(C/O)','R$\tex_{pl}$(R$\tex_{Jup}$,1 bar)','cloud opacity ln($\sigma$/$\sigma$$\tex_0$)','Haze Opacity ln($\sigma$/$\sigma$$\tex_0$-1)','VMR ln(H$\tex_2$)','VMR ln(He)','VMR ln(H$\tex_2$O)','VMR ln(CO$\tex_2)','VMR ln(CO)','VMR ln(CH$\tex_4$)','VMR ln(NH$\tex_3)','VMR ln(Na)','VMR ln(K)','VMR ln(Li)','VMR ln(Rb)','VMR ln(Cs)','VMR ln(TiO)','VMR ln(VO)','VMR ln(FeH)','VMR ln(PH$\tex_3)','VMR ln(H$\tex_2$S)','VMR ln(HCN)','VMR ln(C$\tex_2$H$\tex_2$)']

parameternames=['Temperature (K)','Metallicity ln(Z/Z$\tex_\odot$)','C/O Ratio ln(C/O)','R$\tex_{pl}$(R$\tex_{Jup}$,1 mbar)','cloud opacity ln($\sigma$/$\sigma$$\tex_0$)','Haze Opacity ln($\sigma$/$\sigma$$\tex_0$-1)','Infrared opacity $\kappa$$\tex_{IR}$ (cm$\tex^2$/g)','Optical/IR opacity ratio $\gamma$$\tex_1$','Optical/IR opacity ratio $\gamma$$\tex_2$','$\beta$','$\alpha$','VMR ln(H$\tex_2$)','VMR ln(He)','VMR ln(H$\tex_2$O)','VMR ln(CO$\tex_2)','VMR ln(CO)','VMR ln(CH$\tex_4$)','VMR ln(NH$\tex_3$)','VMR ln(Na)','VMR ln(K)','VMR ln(Li)','VMR ln(Rb)','VMR ln(Cs)','VMR ln(TiO)','VMR ln(VO)','VMR ln(FeH)','VMR ln(PH$\tex_3)','VMR ln(H$\tex_2$S)','VMR ln(HCN)','VMR ln(C$\tex_2$H$\tex_2$)']

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
    cgContour,100-100*exp(-0.5*(chi-min(chi))), pall(*,k), pall(*,l), levels=levels, /Irregular, /Outline,label=1,/overplot,thick=3,resolution=[res,res]

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

;res=15
filenam='ait_monitor_corner'
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
        cghistoplot,pall(*,j),/fill,locations=locations,nbins=histobins,histdata=histdata,DATACOLORNAME='Dark Grey', POLYCOLOR='sky blue',xtitle=parameternames(tofit(j)),/frequenc,xthick=2,ythick=2,charthick=2,xcharsize=xchrsz,ycharsize=ychrsz,Position=pos[*,k],NoErase=j NE 0,xrange=[min(pall(*,0)),max(pall(*,0))],/xstyle;,yrange=[0,max(histdata)]
        plots,MCMCpars(k),0 & plots,MCMCpars(k),max(histdata)*1.05,/continue,thick=2
        plots,MCMCpars(k)+MCMCPCERRORlow(k),0 & plots,MCMCpars(k)+MCMCPCERRORlow(k),max(histdata)*1.05,/continue,linestyle=2,thick=2
        plots,MCMCpars(k)+MCMCPCERRORhigh(k),0 & plots,MCMCpars(k)+MCMCPCERRORhigh(k),max(histdata)*1.05,/continue,linestyle=2,thick=2
     endif
     ;Plot diagnal histograms
     if (j eq i and j ne 0 and i ne 0) then begin 
       cghistoplot,pall(*,j),/fill,locations=locations,nbins=histobins,histdata=histdata,DATACOLORNAME='Dark Grey', POLYCOLOR='sky blue',xtitle=parameternames(tofit(j)),/frequenc,xthick=2,ythick=2,charthick=2,xcharsize=xchrsz,ycharsize=ychrsz,yticklen=0,yminor=1,ytickinterval=1,Position=pos[*,k],NoErase=j NE 0,xrange=[min(pall(*,i)),max(pall(*,i))],/xstyle
       ;overplot 1-sigma ranges
       cgplot,[MCMCpars(i),MCMCpars(i)],[0,1],/overplot,thick=2,Position=pos[*,k], NoErase=j NE 0
       cgplot,[MCMCpars(i)+MCMCPCERRORlow(i),MCMCpars(i)+MCMCPCERRORlow(i)],[0,1],/overplot,thick=2,linestyle=2,Position=pos[*,k], NoErase=j NE 0
       cgplot,[MCMCpars(i)+MCMCPCERRORhigh(i),MCMCpars(i)+MCMCPCERRORhigh(i)],[0,1],/overplot,thick=2,linestyle=2,Position=pos[*,k], NoErase=j NE 0
     endif
     ;label Y-axis histogram axis on left side for all but top corner
     if (j eq i and j ne 0 and i ne 0) then cgaxis,yaxis=1,ytitle='Relative Frequency',ycharsize=0.25,xthick=2,ythick=2,charthick=2;,yrange=[0,max(histdata)]
     ;plot Distributions
     if (i lt j and j ne i) then cgContour,100-100*exp(-0.5*(chi-min(chi))), pall(*,i), pall(*,j), NLevels=10, /Irregular, /Fill, /Outline,xtitle=parameternames(tofit(i)),ytitle=parameternames(tofit(j)),label=0,xthick=1,ythick=1,charthick=1,resolution=[res,res],xcharsize=xchrsz,ycharsize=ychrsz,Position=pos[*,k], NoErase=j NE 0,olevels=olevels,xrange=[min(pall(*,i)),max(pall(*,i))],/xstyle
     ;if (i lt j and j ne i) then cgContour,100-100*exp(-0.5*(chi-min(chi))), pall(*,i), pall(*,j), levels=levels, /Irregular, /Outline,label=1,/overplot,thick=2,resolution=[res,res],charsize=0.2,Position=pos[*,k], NoErase=j NE 0
     xchrsz=0.001
     ychrsz=0.001
     k=k+1  
  ENDFOR
ENDFOR
print,olevels
  cgLoadCT, 6,/brewer,/reverse;, /Brewer, /Reverse
   cgColorbar, /Vertical, Range=[0,1],charsize=1,charthick=2, Position=[0.96, 0.5, 0.99, 0.90],title='1-e$\tex^{\Delta\chi^2}$$\tex^{/2}$'
   cgPS_Close


  spawn,'ps2eps -B -C -f '+dirs+filenam+'.ps'
  spawn,'ps2pdf -dEPSCrop -dAutoRotatePages=/None -dAutoFilterColorImages=false -dColorImageFilter=/FlateEncode -dUseFlateCompression=true '+dirs+filenam+'.eps temporary.pdf'
  spawn,'pdftk temporary.pdf cat 1L output '+dirs+filenam+'.pdf'
;  spawn,'open '+dirs+filenam+'.ps'


!P.Multi = 0
end



