pro arc_density,dirs,chi=chi
; Monitor MCMC convergence
; To be run in a separte IDL session to ARC
; Loads current MCMC chain, calculates Gelman-Rubin statistics
; Requires: Coyote Library
; http://www.idlcoyote.com/documents/programs.php
!path = '../../arc/coyote/:'+ !path                                        ;remember to change arc_multinest.pro threads setup

dirsnow=dirs
restore,dirs+'arc_inputs.sav'
dirs=dirsnow

;set_plot,'ps'
!p.font=0 
;loadct,0
;set_plot,'x'
restore,dirs+'MCMC.monitor.sav' ;allchains,allchainschi2,i,j,nchains,nthin,tofit,Chemistry

close,1
openw,1,'arc_monitor_density.statistics.txt'
printf,1,'UTC '+systime(/utc)
printf,1,''
printf,1,Planet
printf,1,'Data Source: ',DATAFILE
printf,1,FILENAM

if (Chemistry eq 0) then printf,1,'Chemical Equilibrium Retrieval'
if (Chemistry eq 1) then printf,1,'Free-Chemistry Retrieval'
if (Chemistry eq 2) then printf,1,'Non-Equilibrium Chemsitry Retrieval'
if (Chemistry eq 3) then printf,1,'Equilibrium Chemsitry Retrieval with Free Non-Solar Element Abundances'
printf,1,''

   ;Chains may not find CHI^2 minimum, go rouge, ect.  Find them and destroy them
   fin=i
   strt=0 & if (i gt 4000) then strt=2000
   ;calculate average chi^2 of each chain to compare
   chiave=dblarr(nchains) & for w=0,nchains-1 do chiave(w)=mean(allchainschi2(strt:fin,w)-min(allchainschi2(0:i,*)))
   sigma_cut=2  ;sigma cut on average chi^2 of the chains
   print,'min CHI^2 = ',min(allchainschi2(0:i,*))
   print,'Nfree = ',n_elements(tofit)
   printf,1,'min CHI^2 = ',min(allchainschi2(0:i,*))
   printf,1,'Nfree = ',n_elements(tofit)
      
   for k=0,nchains-1 do if min(allchainschi2(0:i,k)) eq min(allchainschi2(0:i,*)) then begin print,k,min(allchainschi2(0:i,k),minindex),minindex & minchain=k & endif
   print,'min Parameters = ',allchains(*,minindex,minchain)
   print,'    from chain :',minchain,' index ',minindex
   print,''
   print,'sigma cut:',sigma_cut
   resistant_mean,chiave,sigma_cut,mean,sigma,num_rejected,goodvec=good & print,'rejected ',num_rejected,'       good:',good  
   print,'Good chains have ave CHI^2 values within delta Chi^2 :',sigma
   bad=where((chiave-mean)/(sigma_cut*sigma) gt sigma_cut)  
   print,'Bad chains have ave delta CHI^2 values of :',(chiave(bad)-mean)

   printf,1,'min Parameters = ',allchains(*,minindex,minchain)
   printf,1,'    from chain :',minchain,' index ',minindex
   printf,1,''
   printf,1,'sigma cut:',sigma_cut
   printf,1,'Good chains have ave CHI^2 values within delta Chi^2 :',sigma
   printf,1,'Bad chains have ave delta CHI^2 values of :',(chiave(bad)-mean)
   
   
  ;just keep good chains
   include=good ;
   ;include=[0,1,2,3,4,5,6,7,8,9,10,11] ;ALL
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
if keyword_set(chi) then begin
 set_plot,'x'
 plot,allchainschi2(strt:fin,0)-min(allchainschi2),/ylog,title='chi2',yrange=[min(allchainschi2(strt:fin,*)),max(allchainschi2(strt:fin,*))]
 for w=1,n_elements(include)-1 do oplot,allchainschi2(strt:fin,w)-min(allchainschi2),color=1000*w
 set_plot,'ps'
endif


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
histobins=20;25;25;100;28 ;45
res=8;12;100;50;25
denres=90;75;90;30;45;90

nfit=n_elements(tofit) 
; Note: This indexing will not work in the case of a fixed-temperature retrieval

;redefine variable Names
;   parameternamessimple=['Temp ','lnZ  ','logg','pp_Rp','Cloud','haze','CloudTop','haze_alpha','Kzz','fsed','sigma_g','Kir','gamma1','gamma2','beta','alpha','Tint','delta ','Pdelta']
   parameternames=['Temperature (K)','Metallicity ln(Z/Z$\tex_{sun}$)','log(g)','R$\tex_{pl}$(R$\tex_{Jup}$,1 mbar)',$
                   'cloud opacity ln($\sigma$/$\sigma$$\tex_0$)','Haze Opacity ln($\sigma$/$\sigma$$\tex_0$)','Cloud Top (R$\tex_{pl}$)','Haze_$\alpha$','K$\tex_{ZZ}$ (cm$\tex^2$ s)','fsed','sigma_g',$ 
                   'IR opacity $\kappa$$\tex_{IR}$ (cm$\tex^2$/g)','Opt/IR opacity $\gamma$$\tex_1$','Opt/IR opacity $\gamma$$\tex_2$','$\beta$','$\alpha$', 'Temp Int (K)','$\delta$','P$\delta$']
parameternames_main=parameternames

if (Chemistry eq 1) then begin
     Molecule_names_ln    =  ['VMR ln(H$\tex_2$)','VMR ln(He)','VMR ln(H$\tex_2$O)','VMR ln(CO$\tex_2$)','VMR ln(CO)','VMR ln(CH$\tex_4$)','VMR ln(NH$\tex_3$)','VMR ln(Na)','VMR ln(K)','VMR ln(Li)',$
                              'VMR ln(Rb)','VMR ln(Cs)','VMR ln(TiO)','VMR ln(VO)','VMR ln(FeH)','VMR ln(PH$\tex_3$)','VMR ln(H$\tex_2$S)','VMR ln(HCN)','VMR ln(C$\tex_2$H$\tex_2$)','VMR ln(SO$\tex_2$)','VMR ln(Fe)','VMR ln(H_)']
     Molecule_names_log   =  ['VMR log(H$\tex_2$)','VMR log(He)','VMR log(H$\tex_2$O)','VMR log(CO$\tex_2$)','VMR log(CO)','VMR log(CH$\tex_4$)','VMR log(NH$\tex_3$)','VMR log(Na)','VMR log(K)','VMR log(Li)',$
                             'VMR log(Rb)','VMR log(Cs)','VMR log(TiO)','VMR log(VO)','VMR log(FeH)','VMR log(PH$\tex_3$)','VMR log(H$\tex_2$S)','VMR log(HCN)','VMR log(C$\tex_2$H$\tex_2$)','VMR ln(SO$\tex_2$)','VMR ln(Fe)','VMR log(H_)']
;     parameternamessimple=[parameternamessimple,Molecule_names]
     parameternames=[parameternames,Molecule_names_ln]
     parameternameslog=[parameternames_main,Molecule_names_log]
endif


if (Chemistry eq 3) then begin
     Element_names      =  ['H','He','C','N','O','Na','K','Si','Ar','Ti','V','S','Cl','Mg','Al','Ca','Fe','Cr','Li','Cs','Rb','F','P'] ;Spectrally active molecules
     Element_names_ln   =  ['ln(H/H$\tex_{sun}$)','ln(He/He$\tex_{sun}$)','ln(C/C$\tex_{sun}$)','ln(N/N$\tex_{sun}$)','ln(O/O$\tex_{sun}$)','ln(Na/Na$\tex_{sun}$)','ln(K/K$\tex_{sun}$)','ln(Si/Si$\tex_{sun}$)','ln(Ar/Ar$\tex_{sun}$)','ln(Ti/Ti$\tex_{sun}$)','ln(V/V$\tex_{sun}$)','ln(S/S$\tex_{sun}$)','ln(Cl/Cl$\tex_{sun}$)','ln(Mg/Mg$\tex_{sun}$)','ln(Al/Al$\tex_{sun}$)','ln(Ca/Ca$\tex_{sun}$)','ln(Fe/Fe$\tex_{sun}$)','ln(Cr/Cr$\tex_{sun}$)','ln(Li/Li$\tex_{sun}$)','ln(Cs/Cs$\tex_{sun}$)','ln(Rb/Rb$\tex_{sun}$)','ln(F/F$\tex_{sun}$)','ln(P/P$\tex_{sun}$)'] ;Spectrally active molecules
     Element_names_log  =  ['log(H/H$\tex_{sun}$)','log(He/He$\tex_{sun}$)','log(C/C$\tex_{sun}$)','log(N/N$\tex_{sun}$)','log(O/O$\tex_{sun}$)','log(Na/Na$\tex_{sun}$)','log(K/K$\tex_{sun}$)','log(Si/Si$\tex_{sun}$)','log(Ar/Ar$\tex_{sun}$)','log(Ti/Ti$\tex_{sun}$)','log(V/V$\tex_{sun}$)','log(S/S$\tex_{sun}$)','log(Cl/Cl$\tex_{sun}$)','log(Mg/Mg$\tex_{sun}$)','log(Al/Al$\tex_{sun}$)','log(Ca/Ca$\tex_{sun}$)','log(Fe/Fe$\tex_{sun}$)','log(Cr/Cr$\tex_{sun}$)','log(Li/Li$\tex_{sun}$)','log(Cs/Cs$\tex_{sun}$)','log(Rb/Rb$\tex_{sun}$)','log(F/F$\tex_{sun}$)','log(P/P$\tex_{sun}$)'] ;Spectrally active molecules
     numabunds=where(Element_fit_switch eq 0) ;index of activly fit elements
     Element_names = Element_names(numabunds)
     Element_names_ln   = Element_names_ln(numabunds)
     Element_names_log  = Element_names_log(numabunds)
;     parameternamessimple=[parameternamessimple,Element_names]
     parameternames=[parameternames,Element_names_ln]
     parameternameslog=[parameternames_main,Element_names_log]
     parameternameslog(1)='Metallicity log(Z/Z$\tex_{sun}$)'
endif







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
      
      ;; calculate the Gelman-Rubin statistic (remove burn-in)
      ;converged = exofast_gelmanrubin(pars[0:nfit-1,burnndx:i,*],$
      ;                                gelmanrubin,tz,angular=gelmanangular)

  converged = exofast_gelmanrubin(allchains[0:nfit-1,burnndx:i,*],gelmanrubin,tz,angular=gelmanangular) 
  print,'Burn index:',burnndx
  print,'Gelman-Rubin < 1.01 => MCMC is Converged, and Indep. Draws > 1000' 
  print,'Gelman-Rubin      =',gelmanrubin 
  print,'Num. Indep. Draws =',tz
  print,'Acceptance %      =',NACCEPT/(i*nchains)*100
  printf,1,'Burn index:',burnndx
  printf,1,'Gelman-Rubin < 1.01 => MCMC is Converged, and Indep. Draws > 1000 [Strict Convergence]'
  printf,1,'Gelman-Rubin      =',gelmanrubin
  printf,1,'Num. Indep. Draws =',tz
  printf,1,'Acceptance %      =',NACCEPT/(i*nchains)*100  
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
    printf,1,'MCMC of '+string(j,format='(I3)'),' chains with '+string(i+1,format='(I9)')+' steps run out of ',N_ELEMENTS(allchainschi2(*,0))
    printf,1,'Steps to Converge =   ',stepstoconvergence/j,' for each chain' & $
    printf,1,'Taken ',i,' steps with each chain, which is ',(i*j)/stepstoconvergence*100,' % of needed for G-R<1.01'
    printf,1,'Taken ',i,' steps with each chain, which is ',((min(tz))/1000)*100,' % of needed for 1000 draws'
    


;ln units
MCMCpars=dblarr(nfit)
MCMCPCERRORlow=dblarr(nfit)
MCMCPCERRORhigh=dblarr(nfit)
;log10 units
LMCMCpars=dblarr(nfit)
LMCMCPCERRORlow=dblarr(nfit)
LMCMCPCERRORhigh=dblarr(nfit)

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
if (nfit eq 12) then pall=[[p0],[p1],[p2],[p3],[p4],[p5],[p6],[p7],[p8],[p9],[p10],[p11]]
if (nfit eq 13) then pall=[[p0],[p1],[p2],[p3],[p4],[p5],[p6],[p7],[p8],[p9],[p10],[p11],[p12]]
if (nfit eq 14) then pall=[[p0],[p1],[p2],[p3],[p4],[p5],[p6],[p7],[p8],[p9],[p10],[p11],[p12],[p13]]
if (nfit eq 15) then pall=[[p0],[p1],[p2],[p3],[p4],[p5],[p6],[p7],[p8],[p9],[p10],[p11],[p12],[p13],[p14]]


;parameternames=['Temperature (K)','Metallicity ln(Z/Z$\tex_{sun}$)','Reference Radius ($\tex\DeltaR_{pl}$/R$\tex_{star}$)','cloud opacity ln($\sigma$/$\sigma$$\tex_0$)']


title='1-e^{-\Delta\chi^{2}}^{/2}  (%)'

chi=reform(ALLCHAINSCHI2(burnndx:i,*),a)

   levels   =[68,95,99]
   cc=floor((levels)/95.*256)

;Calculate Convidence intervals of all parameters and print
print,'--------'
printf,1,'-------'
for p=0,nfit-1 do begin
    e=percentiles(allchains(p,burnndx:i,*))-median(allchains(p,burnndx:i,*)) 
    e2=percentiles(allchains(p,burnndx:i,*),CONFLIMIT=0.95449)-median(allchains(p,burnndx:i,*)) 
    e3=percentiles(allchains(p,burnndx:i,*),CONFLIMIT=0.9973002)-median(allchains(p,burnndx:i,*)) 
    medpar=median(allchains(p,burnndx:i,*)) & $
    meanpar=median(allchains(p,burnndx:i,*)) & $
  print,parameternames(tofit(p))+string(medpar)+' 1-sig:'+string(e(0))+string(e(1))+'  2-sig:'+string(e2(0))+string(e2(1))+'  3-sig:'+string(e3(0))+string(e3(1))
  printf,1,parameternames(tofit(p))+string(medpar)+' 1-sig:'+string(e(0))+string(e(1))+'  2-sig:'+string(e2(0))+string(e2(1))+'  3-sig:'+string(e3(0))+string(e3(1))
 endfor




;printout posterior abudnances in Log_10
print,'-------- in Log_10 units -----'
printf,1,'-------- in Log_10 units -----'
if (tofit(0) eq 0) then freefit=(where(tofit gt 9)) ; array indicies of free fit molecules
if (tofit(0) ne 0) then freefit=(where(tofit gt 14))

;Make Chains with Log10
Lpall=pall
MdHindex=where(tofit eq 1) 
CORatindex=where(tofit eq 2) 
if (MdHindex(0) ne -1) then Lpall(*,MdHindex)=alog10(exp(pall(*,MdHindex)))
if (CORatindex(0) ne -1) then Lpall(*,CORatindex)=alog10(exp(pall(*,CORatindex)))
for p=0,n_elements(freefit)-1 do Lpall(*,freefit(p))=alog10(exp(pall(*,freefit(p))))


for p=0,nfit-1 do begin
    eL=percentiles(Lpall(*,p))-median(Lpall(*,p)) 
    Le3=percentiles(Lpall(*,p),CONFLIMIT=0.9973002)-median(Lpall(*,p)) 
    Lmedpar=median(Lpall(*,p)) & $
    LMCMCpars(p)=Lmedpar
    LMCMCPCERRORlow(p)=eL(0)
    LMCMCPCERRORhigh(p)=eL(1)
endfor


if (MdHindex(0) ne -1) then begin
    Lle=percentiles(Lpall(*,MdHindex))-median(Lpall(*,MdHindex)) 
    Le2=percentiles(Lpall(*,MdHindex),CONFLIMIT=0.95449)-median(Lpall(*,MdHindex) )
    Le3=percentiles(Lpall(*,MdHindex),CONFLIMIT=0.9973002)-median(Lpall(*,MdHindex))
    Lmedpar=median(Lpall(*,MdHindex)) & $
    Lmeanpar=median(Lpall(*,MdHindex)) & $
    print,parameternameslog(tofit(MdHindex))+string(Lmedpar)+' 1-sig:'+string(Lle(0))+string(Lle(1))+'  2-sig:'+string(Le2(0))+string(Le2(1))+'  3-sig:'+string(Le3(0))+string(Le3(1))
    printf,1,parameternameslog(tofit(MdHindex))+string(Lmedpar)+' 1-sig:'+string(Lle(0))+string(Lle(1))+'  2-sig:'+string(Le2(0))+string(Le2(1))+'  3-sig:'+string(Le3(0))+string(Le3(1))
endif
if (CORatindex(0) ne -1) then begin
    Lle=percentiles(Lpall(*,CORatindex))-median(Lpall(*,CORatindex)) 
    Le2=percentiles(Lpall(*,CORatindex),CONFLIMIT=0.95449)-median(Lpall(*,CORatindex) )
    Le3=percentiles(Lpall(*,CORatindex),CONFLIMIT=0.9973002)-median(Lpall(*,CORatindex))
    Lmedpar=median(Lpall(*,CORatindex)) & $
    Lmeanpar=median(Lpall(*,CORatindex)) & $
    print,parameternameslog(tofit(CORatindex))+string(Lmedpar)+' 1-sig:'+string(Lle(0))+string(Lle(1))+'  2-sig:'+string(Le2(0))+string(Le2(1))+'  3-sig:'+string(Le3(0))+string(Le3(1))
    printf,1,parameternameslog(tofit(CORatindex))+string(Lmedpar)+' 1-sig:'+string(Lle(0))+string(Lle(1))+'  2-sig:'+string(Le2(0))+string(Le2(1))+'  3-sig:'+string(Le3(0))+string(Le3(1))   
endif
if (freefit(0) gt -1) then begin
  for p=0,n_elements(freefit)-1 do begin
    Lle=percentiles(alog10(exp(allchains(freefit(p),burnndx:i,*))))-median(alog10(exp(allchains(freefit(p),burnndx:i,*)))) 
    Le2=percentiles(alog10(exp(allchains(freefit(p),burnndx:i,*))),CONFLIMIT=0.95449)-median(alog10(exp(allchains(freefit(p),burnndx:i,*)) ))
    Le3=percentiles(alog10(exp(allchains(freefit(p),burnndx:i,*))),CONFLIMIT=0.9973002)-median(alog10(exp(allchains(freefit(p),burnndx:i,*)) ))
    Lmedpar=median(alog10(exp(allchains(freefit(p),burnndx:i,*)))) & $
    Lmeanpar=median(alog10(exp(allchains(freefit(p),burnndx:i,*)))) & $
    print,parameternameslog(tofit(freefit(p)))+string(Lmedpar)+' 1-sig:'+string(Lle(0))+string(Lle(1))+'  2-sig:'+string(Le2(0))+string(Le2(1))+'  3-sig:'+string(Le3(0))+string(Le3(1))
    printf,1,parameternameslog(tofit(freefit(p)))+string(Lmedpar)+' 1-sig:'+string(Lle(0))+string(Lle(1))+'  2-sig:'+string(Le2(0))+string(Le2(1))+'  3-sig:'+string(Le3(0))+string(Le3(1))
 endfor  
endif
print,'--------'

;reduce text size with lots of free parameters
textsize=0.3
if (nfit gt 9) then textsize=0.2


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
     if (j eq nfit-1) then xchrsz=textsize
     if (i eq 0) then ychrsz=textsize
     ;plot top corner histogram
     if (j eq 0 and i eq 0) then begin
        cghistoplot,Lpall(*,j),/fill,locations=locations,nbins=histobins,histdata=histdata,DATACOLORNAME='Dark Grey', POLYCOLOR='Charcoal',xtitle=parameternameslog(tofit(j)),/frequenc,xthick=2,ythick=2,charthick=2,xcharsize=xchrsz,ycharsize=0.001,Position=pos[*,k],NoErase=j NE 0,xrange=[min(Lpall(*,0)),max(Lpall(*,0))],thick=2,PEAK_HEIGHT=1,yrange=[0,1.05],ystyle=1
        axis,yaxis=1,ytitle='Relative Freq.',xthick=2,ythick=2,charthick=2,ycharsize=textsize,yrange=[0,1.05],ystyle=1
       cgplot,[LMCMCpars(i),LMCMCpars(i)],[0,1.05],/overplot,thick=3,Position=pos[*,k], color=cgcolor("Orange",!D.Table_Size-2), NoErase=j NE 0
       cgplot,[LMCMCpars(i)+LMCMCPCERRORlow(i),LMCMCpars(i)+LMCMCPCERRORlow(i)],[0,1.05], color=cgcolor("Orange",!D.Table_Size-2),/overplot,thick=3,linestyle=1,Position=pos[*,k], NoErase=j NE 0
       cgplot,[LMCMCpars(i)+LMCMCPCERRORhigh(i),LMCMCpars(i)+LMCMCPCERRORhigh(i)],[0,1.05], color=cgcolor("Orange",!D.Table_Size-2),/overplot,thick=3,linestyle=1,Position=pos[*,k], NoErase=j NE 0
     endif
     ;Plot diagnal histograms
     if (j eq i and j ne 0 and i ne 0) then begin 
       cghistoplot,Lpall(*,j),/fill,locations=locations,nbins=histobins,histdata=histdata,DATACOLORNAME='Dark Grey', POLYCOLOR='Charcoal',xtitle=parameternameslog(tofit(j)),/frequenc,xthick=2,ythick=2,charthick=2,xcharsize=xchrsz,ycharsize=ychrsz,yticklen=0,yminor=1,ytickinterval=1,Position=pos[*,k],NoErase=j NE 0,xrange=[min(Lpall(*,i)),max(Lpall(*,i))],thick=2,PEAK_HEIGHT=1,yrange=[0,1.05],ystyle=1
       ;overplot 1-sigma ranges
       cgplot,[LMCMCpars(i),LMCMCpars(i)],[0,1.05],/overplot,thick=3,Position=pos[*,k], color=cgcolor("Orange",!D.Table_Size-2), NoErase=j NE 0
       cgplot,[LMCMCpars(i)+LMCMCPCERRORlow(i),LMCMCpars(i)+LMCMCPCERRORlow(i)],[0,1.05], color=cgcolor("Orange",!D.Table_Size-2),/overplot,thick=3,linestyle=1,Position=pos[*,k], NoErase=j NE 0
       cgplot,[LMCMCpars(i)+LMCMCPCERRORhigh(i),LMCMCpars(i)+LMCMCPCERRORhigh(i)],[0,1.05], color=cgcolor("Orange",!D.Table_Size-2),/overplot,thick=3,linestyle=1,Position=pos[*,k], NoErase=j NE 0
     endif
     ;label Y-axis histogram axis on left side for all but top corner
     if (j eq i and j ne 0 and i ne 0) then cgaxis,yaxis=1,ytitle='Relative Freq.',ycharsize=textsize,xthick=2,ythick=2,charthick=2,yrange=[0,1.05],ystyle=1;,yrange=[0,max(histdata)]
     ;plot Distributions
     if (i lt j and j ne i) then begin
       red_band = Lpall(*,i)
       green_band = Lpall(*,j)
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
          cgContour,100-100*exp(-0.5*(chi-min(chi))), Lpall(*,i), Lpall(*,j), NLevels=10, /Irregular, /Fill, /Outline,xtitle=parameternameslog(tofit(i)),ytitle=parameternameslog(tofit(j)),label=0,xthick=2,ythick=2,charthick=1,resolution=[res,res],xcharsize=xchrsz,ycharsize=ychrsz,Position=pos[*,k], NoErase=j NE 0,olevels=olevels,xrange=[min(Lpall(*,i)),max(Lpall(*,i))],/xstyle,cell_fill=1       
;          cgimage,h2d,XRange=xrange, YRange=yrange, /Axes, Palette=palette,  XTITLE=parameternames(tofit(i)), $
;          YTITLE=parameternames(tofit(j)),Position=pos[*,k],NoErase=j NE 0,charsize=ychrsz,/overplot
 
         cgimage,h2d,XRange=xrange, YRange=yrange, /Axes, Palette=palette,Position=pos[*,k],NoErase=j NE 0,charsize=0.001,/overplot

         cgplot,[LMCMCpars[i],LMCMCpars[i]],[-10000,10000],/overplot,thick=2,Position=pos[*,k], color=cgcolor("Orange",!D.Table_Size-2), NoErase=j NE 0
         cgplot,[-10000,10000],[LMCMCpars[j],LMCMCpars[j]],/overplot,thick=2,Position=pos[*,k], color=cgcolor("Orange",!D.Table_Size-2), NoErase=j NE 0
       
         plotsym,8, /fill
         plots, [LMCMCpars[i]],[LMCMCpars[j]], psym=8, symsize=0.4, color=cgcolor("Orange",!D.Table_Size-2)
         plotsym,8, /fill
         plots, [LMCMCpars[i]],[LMCMCpars[j]], psym=8, symsize=0.2, color=cgcolor("Red",!D.Table_Size-2)

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
;  spawn,'pdftk temporary.pdf cat 1L output '+dirs+filenam+'.pdf'
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
  ;b=where((h2O) gt 67.0) & c=where((h2O) lt 67.) & print,string((1-double(n_elements(b))/double(n_elements(c)))*100),' % of distribution < 67x solar'
  ;b=where((h2O) gt 89.0) & c=where((h2O) lt 89.) & print,string((1-double(n_elements(b))/double(n_elements(c)))*100),' % of distribution < 89x solar'
  set_plot,'ps'                                                                                                                                       
 device,/landscape,filename=dirs+'H2Ohisto.ps'
  cghistoplot,alog10(h2o),/frequency,locations=locations,histdata=histdata,nbins=30,POLYCOLOR='sky blue',/fill,DATACOLORNAME='Dark Grey',xthick=4,ythick=4,charthick=4,xtitle='log$\tex_{10}$(H$\tex_{2}$O x Solar)'
 device,/close
endif

;print,water abudnance in freechemistry case
   ; Isothermal       10   11   12    13    14   15    16    17  18  19   20   21   22   23     24   25   26     27     28                          
   ; T-P fit          15   16   17    18    19   20    21    22  23  24   25   26   27   28     29   30   31     32     33
  ;Opacity_names   =['H2','He','H2O','CO2','CO','CH4','NH3','Na','K','Li','Rb','Cs','TiO','VO','FeH','PH3','H2S','HCN','C2H2'] ;Spectrally active molecules

a=where(tofit eq 2) 
if (a(0) eq -1) then begin

if (tofit(0) eq 0) then begin ; isothermal T-P case
  h20indx=where(tofit eq 12)  
  CO2indx=where(tofit eq 13)
  COindx =where(tofit eq 14)
  CH4indx=where(tofit eq 15)
  TiOindx=where(tofit eq 22)
   VOindx=where(tofit eq 23)
endif
if (tofit(0) ne 0) then begin; fit T-P case
  h20indx=where(tofit eq 17)  
  CO2indx=where(tofit eq 18)
  COindx =where(tofit eq 19)
  CH4indx=where(tofit eq 20)
  TiOindx=where(tofit eq 27)
   VOindx=where(tofit eq 28)
endif

solarH2Ovmr=0.000360745215794534;6E-4
;if (a(0) ne -1) then begin
;   O=(exp(pall(*,5)+pall(*,0))*5.754E-04)  ;O/H
;  C=(exp(pall(*,4)+pall(*,0))*3.162E-04)   ;C/H
;  set_plot,'ps'                                                                                                                                       
;  device,/landscape,filename=dirs+'COratio.ps'
;  cghistoplot,(c/o),/frequency,PEAK_HEIGHT=1,locations=locations,histdata=histdata,nbins=30,POLYCOLOR='sky blue',/fill,DATACOLORNAME='Dark Grey',xthick=4,ythick=4,charthick=4,xtitle='C/O ratio'
;   e=percentiles(alog10(c/o))-median(alog10(c/o)) 
;  e2=percentiles(alog10(c/o),CONFLIMI=0.95449)-median(alog10(c/o))
;  e3=percentiles(alog10(c/o),CONFLIMI=0.9973)-median(alog10(c/o))
;  
;  plots,median(alog10(c/o)),0 & plots,median(alog10(c/o)),max(histdata)*1.05,/continue,thick=7
;  plots,median(alog10(c/o))+e(0),0 & plots,median(alog10(c/o))+e(0),max(histdata)*1.05,/continue,linestyle=2,thick=7
;  plots,median(alog10(c/o))+e(1),0 & plots,median(alog10(c/o))+e(1),max(histdata)*1.05,/continue,linestyle=2,thick=7
;  print,'log10(C/O ratio): '+string(median(alog10(c/o)))+' 1-sig:'+string(e(0))+string(e(1))+' 2-sig:'+string(e2(0))+string(e2(1))+' 3-sig:'+string(e3(0))+string(e3(1))
;  eln=percentiles(alog(c/o))-median(alog(c/o)) 
;  print,'ln(C/O ratio): ',median(alog(c/o)),eln(0),eln(1)

;  device,/close

endif

save,filename=dirs+'ARC_density.sav',/var,/comm
;device,/close
close,1


;=====
;plot best-fit mcmc
;restore,dirs+'arc_inputs.sav' 
;restore,dirs+filenam+'arc_retrieval.setup.sav'
;arc_spectra_plot,dirs,filenam,MPFITfilenameprefix,params,wsd,r,werr,rerr,modelpoints,ws,cloudtype,Emission,Rpl


end



