pro arc_monitor_nested,dirs
; Purpose:  Generate Posterior Distribution from Nested Sampling
;
; Resolution of Histograms/Controurss
  nhistobins=15
  nhistobinsCO=25
  NContourlevels=15   ; Number of Contours for 2D plot10
  colortable=1       ; Color of 2D plot  1==blue
  bottomct= 1;  
  res=35                         ; resolution of 2D plot
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


close,1
openw,1,dirs+'arc_monitor_nested.statistics.txt'
printf,1,'UTC '+systime(/utc)
printf,1,''
printf,1,Planet
printf,1,'Data Source: ',DATAFILE
printf,1,FILENAM
         printf,1,''
         printf,1,'Nested Sampling'
         printf,1,'Settings'
         printf,1,'  NUMBER OF LIVE POINTS: ',Num
         printf,1,'  Recluster Sample Number SAMPLE_NUM     : ',SAMPLE_NUM
         printf,1,'  TOLERANCE ON MEASURING the EVIDENCE LnZ: ',Tol
         printf,1,'  ln Evidence LnZ = '+string(multinest_str.(2)),'+/-',multinest_str.(3)
         printf,1,'  H = '+string(multinest_str.(4))+' information content (bits)'
         printf,1,''


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
;print,m,' ',chi2(m) & $
endfor
print,'Ndata points :',n_elements(x)
print,'min CHI^2 = ',min(chi2)
print,'    CHIv^2 = ',min(chi2)/(n_elements(x)-n_elements(tofit))
print,'DOF =',n_elements(x)-n_elements(tofit)
print,'Nfree = ',n_elements(tofit)
printf,1,'Ndata points :',n_elements(x)
printf,1,'min CHI^2 = ',min(chi2)
printf,1,'    CHIv^2 = ',min(chi2)/(n_elements(x)-n_elements(tofit))
printf,1,'DOF =',n_elements(x)-n_elements(tofit)
printf,1,'Nfree = ',n_elements(tofit)
print,''
printf,1,''
print,'Minimum model parameters:',NESTPARS(*,n_elements(LP)-1)
printf,1,'Minimum model parameters:',NESTPARS(*,n_elements(LP)-1)

print,'--------'
printf,1,'-------'


;===============================redefine variable Names==============================
   parameternames=['Temperature (K)','Metallicity ln(Z/Z$\tex_{sun}$)','log(g)','R$\tex_{pl}$(R$\tex_{Jup}$,1 mbar)',$
                   'cloud opacity ln($\sigma$/$\sigma$$\tex_0$)','Haze Opacity ln($\sigma$/$\sigma$$\tex_0$)','Cloud Top (R$\tex_{pl}$)','Haze_$\alpha$','K$\tex_{ZZ}$ (cm$\tex^2$ s)','fsed','sigma_g',$ 
                   'IR opacity $\kappa$$\tex_{IR}$ (cm$\tex^2$/g)','Opt/IR opacity $\gamma$$\tex_1$','Opt/IR opacity $\gamma$$\tex_2$','$\beta$','$\alpha$', 'Temp Int (K)','$\delta$','P$\delta$']
   parameternames_main=parameternames

   if (Chemistry eq 1) then begin
     Molecule_names_ln    =  ['VMR ln(H$\tex_2$)','VMR ln(He)','VMR ln(H$\tex_2$O)','VMR ln(CO$\tex_2$)','VMR ln(CO)','VMR ln(CH$\tex_4$)','VMR ln(NH$\tex_3$)','VMR ln(Na)','VMR ln(K)','VMR ln(Li)',$
                              'VMR ln(Rb)','VMR ln(Cs)','VMR ln(TiO)','VMR ln(VO)','VMR ln(FeH)','VMR ln(PH$\tex_3$)','VMR ln(H$\tex_2$S)','VMR ln(HCN)','VMR ln(C$\tex_2$H$\tex_2$)','VMR ln(SO$\tex_2$)','VMR ln(Fe)','VMR ln(H_)']
     Molecule_names_log   =  ['VMR log(H$\tex_2$)','VMR log(He)','VMR log(H$\tex_2$O)','VMR log(CO$\tex_2$)','VMR log(CO)','VMR log(CH$\tex_4$)','VMR log(NH$\tex_3$)','VMR log(Na)','VMR log(K)','VMR log(Li)',$
                             'VMR log(Rb)','VMR log(Cs)','VMR log(TiO)','VMR log(VO)','VMR log(FeH)','VMR log(PH$\tex_3$)','VMR log(H$\tex_2$S)','VMR log(HCN)','VMR log(C$\tex_2$H$\tex_2$)','VMR ln(SO$\tex_2$)','VMR ln(Fe)','VMR log(H_)']
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
     parameternames=[parameternames,Element_names_ln]
     parameternameslog=[parameternames_main,Element_names_log]
     parameternameslog(1)='Metallicity log(Z/Z$\tex_{sun}$)'
endif
;===========================================================================
;ln units
MCMCpars=dblarr(nfit)
MCMCPCERRORlow=dblarr(nfit)
MCMCPCERRORhigh=dblarr(nfit)
;log10 units
LMCMCpars=dblarr(nfit)
LMCMCPCERRORlow=dblarr(nfit)
LMCMCPCERRORhigh=dblarr(nfit)
LMCMCPCERROR2low=dblarr(nfit)
LMCMCPCERROR2high=dblarr(nfit)
LMCMCPCERROR3low=dblarr(nfit)
LMCMCPCERROR3high=dblarr(nfit)
;===========================================================================


;======================================================================================================================================================
; Calculate Distrubution Function using probability of each sampled
; parameter, use it to calculate the best-fit and 1,2,3-sigma ranges

Prob=(exp(multinest_str.(1)))   ;probability for each sample
parameter_prob=prob

print,'======================================='
;find median and 1,2,3 sigma of 1D marginalized posterior distriubtions
for p=0,nfit-1 do begin
   ;sort parameters
   srt=sort(NESTPARS[p,*])
   parameter=NESTPARS[p,srt]
   parameter_prob=Prob[srt]

   ;weighted average parmaeter value
   w_ave=total(parameter*parameter_prob) 
   ;print,w_ave
   ;middle from cumulative distribution
   medpar=parameter[value_locate(total(parameter_prob,/cumulative),[0.5d0])]
   ;print,m
   e=parameter[value_locate(total(parameter_prob,/cumulative),[0.5d0-0.68268d0/2.,0.5d0+0.68268d0/2.])]-medpar[0]
   e2=parameter[value_locate(total(parameter_prob,/cumulative),[0.5d0-0.95449d0/2.,0.5d0+0.95449d0/2.])]-medpar[0]
   e3=parameter[value_locate(total(parameter_prob,/cumulative),[0.5d0-0.99730d0/2.,0.5d0+0.99730d0/2.])]-medpar[0]

   print,parameternames(tofit(p))+string(medpar)+' 1-sig:'+string(e(0))+string(e(1))+'  2-sig:'+string(e2(0))+string(e2(1))+'  3-sig:'+string(e3(0))+string(e3(1))
   printf,1,parameternames(tofit(p))+string(medpar)+' 1-sig:'+string(e(0))+string(e(1))+'  2-sig:'+string(e2(0))+string(e2(1))+'  3-sig:'+string(e3(0))+string(e3(1))
   LMCMCpars(p)=medpar
   LMCMCPCERRORlow(p)=e(0)
   LMCMCPCERRORhigh(p)=e(1)
   LMCMCPCERROR2low(p)=e2(0)
   LMCMCPCERROR2high(p)=e2(1)
   LMCMCPCERROR3low(p)=e3(0)
   LMCMCPCERROR3high(p)=e3(1)

endfor
   ;get approximate probability of 1,2,3-sigma
   p=0
   srt=sort(NESTPARS[p,*])
   parameter=NESTPARS[p,srt]
   parameter_prob=Prob[srt]
   e=parameter[value_locate(total(parameter_prob,/cumulative),[0.5d0-0.68268d0/2.,0.5d0+0.68268d0/2.])]
   e2=parameter[value_locate(total(parameter_prob,/cumulative),[0.5d0-0.95449d0/2.,0.5d0+0.95449d0/2.])]
   e3=parameter[value_locate(total(parameter_prob,/cumulative),[0.5d0-0.99730d0/2.,0.5d0+0.99730d0/2.])]
   a1=where(transpose(parameter) ge e[0])  & p1=parameter_prob[a1[0]]
   a2=where(transpose(parameter) ge e2[0]) & p2=parameter_prob[a2[0]]
   a3=where(transpose(parameter) ge e3[0]) & p3=parameter_prob[a3[0]]
   ;print,'1,2,3 sig prob:',p1,p2,p3

;printout posterior abudnances in Log_10
print,'-------- in Log_10 units -----'
printf,1,'-------- in Log_10 units -----'
;if (tofit(0) eq 0) then freefit=(where(tofit gt 9)) ; array indicies of free fit molecules
if (tofit(0) ne 0) then freefit=(where(tofit ge 19))

;Make Chains with Log10
LNESTPARS=NESTPARS
MdHindex=where(tofit eq 1) 
if (MdHindex(0) ne -1) then LNESTPARS[MdHindex,*]=alog10(exp(NESTPARS[MdHindex,*]))
for p=0,n_elements(freefit)-1 do LNESTPARS[freefit(p),*]=alog10(exp(NESTPARS[freefit(p),*]))

if (MdHindex(0) ne -1) then begin
  UPPRIORSPARAMS(MdHindex)=alog10(exp(UPPRIORSPARAMS(MdHindex)))
  LOWPRIORSPARAMS(MdHindex)=alog10(exp(LOWPRIORSPARAMS(MdHindex)))

   p=MdHindex
   ;sort parameters
   srt=sort(LNESTPARS[p[0],*])
   parameter=LNESTPARS[p[0],srt]
   parameter_prob=Prob[srt]
   ;weighted average parmaeter value
   w_ave=total(parameter*parameter_prob) 
   ;middle from cumulative distribution
   medpar=parameter[value_locate(total(parameter_prob,/cumulative),[0.5d0])]
   e=parameter[value_locate(total(parameter_prob,/cumulative),[0.5d0-0.68268d0/2.,0.5d0+0.68268d0/2.])]-medpar[0]
   e2=parameter[value_locate(total(parameter_prob,/cumulative),[0.5d0-0.95449d0/2.,0.5d0+0.95449d0/2.])]-medpar[0]
   e3=parameter[value_locate(total(parameter_prob,/cumulative),[0.5d0-0.99730d0/2.,0.5d0+0.99730d0/2.])]-medpar[0]
   print,parameternameslog(tofit(p))+string(medpar)+' 1-sig:'+string(e(0))+string(e(1))+'  2-sig:'+string(e2(0))+string(e2(1))+'  3-sig:'+string(e3(0))+string(e3(1))
   printf,1,parameternameslog(tofit(p))+string(medpar)+' 1-sig:'+string(e(0))+string(e(1))+'  2-sig:'+string(e2(0))+string(e2(1))+'  3-sig:'+string(e3(0))+string(e3(1))
   LMCMCpars(MdHindex)=medpar
   LMCMCPCERRORlow(MdHindex)=e(0)
   LMCMCPCERRORhigh(MdHindex)=e(1)
   LMCMCPCERROR2low(MdHindex)=e2(0)
   LMCMCPCERROR2high(MdHindex)=e2(1)
   LMCMCPCERROR3low(MdHindex)=e3(0)
   LMCMCPCERROR3high(MdHindex)=e3(1)
endif
if (freefit(0) gt -1) then begin
  for p=0,n_elements(freefit)-1 do begin
                                ;sort parameters
   UPPRIORSPARAMS(freefit(p))=alog10(exp(UPPRIORSPARAMS(freefit(p))))
   LOWPRIORSPARAMS(freefit(p))=alog10(exp(LOWPRIORSPARAMS(freefit(p))))
   srt=sort(LNESTPARS[freefit(p),*])
   parameter=LNESTPARS[freefit(p),srt]
   parameter_prob=Prob[srt]
   ;weighted average parmaeter value
   w_ave=total(parameter*parameter_prob) 
   ;middle from cumulative distribution
   medpar=parameter[value_locate(total(parameter_prob,/cumulative),[0.5d0])]
   e=parameter[value_locate(total(parameter_prob,/cumulative),[0.5d0-0.68268d0/2.,0.5d0+0.68268d0/2.])]-medpar[0]
   e2=parameter[value_locate(total(parameter_prob,/cumulative),[0.5d0-0.95449d0/2.,0.5d0+0.95449d0/2.])]-medpar[0]
   e3=parameter[value_locate(total(parameter_prob,/cumulative),[0.5d0-0.99730d0/2.,0.5d0+0.99730d0/2.])]-medpar[0]
   print,parameternameslog(tofit(freefit(p)))+string(medpar)+' 1-sig:'+string(e(0))+string(e(1))+'  2-sig:'+string(e2(0))+string(e2(1))+'  3-sig:'+string(e3(0))+string(e3(1))
   printf,1,parameternameslog(tofit(freefit(p)))+string(medpar)+' 1-sig:'+string(e(0))+string(e(1))+'  2-sig:'+string(e2(0))+string(e2(1))+'  3-sig:'+string(e3(0))+string(e3(1))
   LMCMCpars(freefit(p))=medpar
   LMCMCPCERRORlow(freefit(p))=e(0)
   LMCMCPCERRORhigh(freefit(p))=e(1)
   LMCMCPCERROR2low(freefit(p))=e2(0)
   LMCMCPCERROR2high(freefit(p))=e2(1)
   LMCMCPCERROR3low(freefit(p))=e3(0)
   LMCMCPCERROR3high(freefit(p))=e3(1)
 endfor  
endif
Onesigma=(abs(LMCMCPCERRORlow)+abs(LMCMCPCERRORhigh))/2.
print,'--------'
printf,1,'--------'
printf,1,'Parameters Reconstructed From Posterior'
;======================================================================================================================================================
;save,filename='temp.sav',/all

;historgram
j=0
  srt=sort(NESTPARS[j,*])
  parameter=NESTPARS[j,srt]
  parameter_prob=Prob[srt]
  ;plot,parameter,total(parameter_prob,/cumulative)
  ;a=lin_rebin(parameter,parameter_prob,binwidth=(0.75*abs(min(Onesigma[j]))),bin_x=bin_x)
  ;cgplot,bin_x,a/max(a),psym=10,xrange=[LOWPRIORSPARAMS(j),UPPRIORSPARAMS(j)],yrange=[0,1.05],xtitle=parameternameslog(tofit(j))
  

       ;cgplot,[LMCMCpars[j],LMCMCpars[j]],[0,1.05],/overplot,thick=3, color=cgcolor("Orange",!D.Table_Size-2);,Position=pos[*,k], NoErase=j NE 0
       ;cgplot,[LMCMCPCERRORlow[j]+LMCMCpars[j],LMCMCPCERRORlow[j]+LMCMCpars[j]],[0,1.05], color=cgcolor("Orange",!D.Table_Size-2),/overplot,thick=3,linestyle=1;,Position=pos[*,k], NoErase=j NE 0
       ;cgplot,[LMCMCPCERRORhigh[j]+LMCMCpars[j],LMCMCPCERRORhigh[j]+LMCMCpars[j]],[0,1.05], color=cgcolor("Orange",!D.Table_Size-2),/overplot,thick=3,linestyle=1;,Position=pos[*,k], NoErase=j NE 0

;2D histogram
       cgLoadCT, colortable,/reverse,bottom=bottomct;33;65
          TVLCT, cgColor('white', /Triple), 0
          TVLCT, r, g, b, /Get
          palette = [ [r], [g], [b] ]     
;cgcontour,Prob,transpose(LNESTPARS[0,*]),transpose(LNESTPARS[1,*]), NLevels=5, /Irregular, /Fill, Palette=palette,xrange=[700,1800],yrange=[-3,3] 



; set X-range and Y-range
; Take Upper and Lower Prior limits, or 4.5-sigma range
expand=2
Lowrange=LOWPRIORSPARAMS
 Highrange=UPPRIORSPARAMS
 upper3range=LMCMCPCERROR3high+LMCMCpars
 lower3range=LMCMCPCERROR3low+LMCMCpars
 diffpriorup=UPPRIORSPARAMS-LMCMCpars
 diff3sigup=(LMCMCPCERROR3high*expand+LMCMCpars)-LMCMCpars
 a=where(diff3sigup lt diffpriorup)
 Highrange(a)=LMCMCPCERROR3high(a)*expand+LMCMCpars(a)
 diffpriorlow=LOWPRIORSPARAMS-LMCMCpars
 diff3siglow=(LMCMCPCERROR3low*expand+LMCMCpars)-LMCMCpars
 a=where(diff3siglow gt diffpriorlow)
 Lowrange(a)=LMCMCPCERROR3low(a)*expand+LMCMCpars(a)

;
 binwidth=(Highrange-Lowrange)/nhistobins  ;(0.75*abs(min(Onesigma[j])))
 
set_plot,'ps'
loadct,colortable,/silent,bottom=bottomct
plotsym,0,1,/fill

textsize=0.3
if (nfit gt 9) then textsize=0.2

filenam='arc_monitor_nested_corner_IDL'
;device,/landscape,filename=filenam+'.ps',set_font='Helvetica',/tt_font
   cgPS_Open, dirs+filenam+'.ps'
   cgDisplay

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
        srt=sort(LNESTPARS[j,*])
        parameter=LNESTPARS[j,srt]
        parameter_prob=Prob[srt]
        a=lin_rebin(parameter,parameter_prob,binwidth=binwidth(j),bin_x=bin_x) ;Rebin historgram 
       cgLoadCT, 0;33;65
        cgplot,bin_x,a/max(a),psym=10,xrange=[Lowrange(j),Highrange(j)],yrange=[0,1.05],Position=pos[*,k],NoErase=j NE 0,xthick=2,ythick=2,charthick=2,xcharsize=xchrsz,ycharsize=0.001,/xstyle,/ystyle
         for z=10,1000 do cgplot,bin_x,a/max(a)*(z*0.001),psym=10,color=cgcolor('light grey'),thick=0.1,Position=pos[*,k], NoErase=j NE 0,/overplot ;fill in the plot
          cgplot,bin_x,a/max(a),psym=10,color=cgcolor('black'),thick=3 ,Position=pos[*,k], NoErase=j NE 0,/overplot
          cgplot,[LMCMCpars[j],LMCMCpars[j]],[0,1.05],/overplot,thick=3, color=cgcolor("Orange",!D.Table_Size-2),Position=pos[*,k], NoErase=j NE 0
          cgplot,[LMCMCPCERRORlow[j]+LMCMCpars[j],LMCMCPCERRORlow[j]+LMCMCpars[j]],[0,1.05], color=cgcolor("Orange",!D.Table_Size-2),/overplot,thick=3,linestyle=1,Position=pos[*,k], NoErase=j NE 0
          cgplot,[LMCMCPCERRORhigh[j]+LMCMCpars[j],LMCMCPCERRORhigh[j]+LMCMCpars[j]],[0,1.05], color=cgcolor("Orange",!D.Table_Size-2),/overplot,thick=3,linestyle=1,Position=pos[*,k], NoErase=j NE 0
          cgaxis,/xaxis,xrange=[Lowrange(j),Highrange(j)],xcharsize=0.001,xthick=2,/xstyle
          cgplot,bin_x,a/max(a),psym=10,xrange=[Lowrange(j),Highrange(j)],yrange=[0,1.05],Position=pos[*,k],NoErase=j NE 0,xthick=2,ythick=2,charthick=2,xcharsize=xchrsz,ycharsize=0.001,/xstyle,/ystyle,/nodata,/overplot
          cgaxis,yaxis=1,ytitle='Relative Freq.',ycharsize=textsize,xthick=2,ythick=2,charthick=2,yrange=[0,1.05],ystyle=1
       endif
     ;Plot diagnal histograms
     if (j eq i and j ne 0 and i ne 0) then begin 
        srt=sort(LNESTPARS[j,*])
        parameter=LNESTPARS[j,srt]
        parameter_prob=Prob[srt]
        a=lin_rebin(parameter,parameter_prob,binwidth=binwidth(j),bin_x=bin_x) ;Rebin historgram 
        cgLoadCT, 0                                                            ;33;65
        ;print,'xcharsize',xchrsz
        cgplot,bin_x,a/max(a),psym=10,xrange=[Lowrange(j),Highrange(j)],yrange=[0,1.05],Position=pos[*,k],NoErase=j NE 0,xthick=2,ythick=2,charthick=2,xcharsize=xchrsz,ycharsize=0.001,color=cgcolor('black',!D.Table_Size-2),/xstyle,/ystyle,xtitle=parameternameslog(tofit(j))
          for z=10,1000 do cgplot,bin_x,a/max(a)*(z*0.001),psym=10,color=cgcolor('light grey'),thick=0.1,Position=pos[*,k], NoErase=j NE 0,/overplot ;fill in the plot
          cgplot,bin_x,a/max(a),psym=10,color=cgcolor('black'),thick=3 ,Position=pos[*,k], NoErase=j NE 0,/overplot
          ;plot 1-sig ranges
          cgplot,[LMCMCpars[j],LMCMCpars[j]],[0,1.05],/overplot,thick=3, color=cgcolor("Orange",!D.Table_Size-2),Position=pos[*,k], NoErase=j NE 0
          cgplot,[LMCMCPCERRORlow[j]+LMCMCpars[j],LMCMCPCERRORlow[j]+LMCMCpars[j]],[0,1.05], color=cgcolor("Orange",!D.Table_Size-2),/overplot,thick=3,linestyle=1,Position=pos[*,k], NoErase=j NE 0
          cgplot,[LMCMCPCERRORhigh[j]+LMCMCpars[j],LMCMCPCERRORhigh[j]+LMCMCpars[j]],[0,1.05], color=cgcolor("Orange",!D.Table_Size-2),/overplot,thick=3,linestyle=1,Position=pos[*,k], NoErase=j NE 0
 
     endif
     ;
     ;label Y-axis histogram axis on left side for all but top corner
     if (j eq i and j ne 0 and i ne 0) then cgaxis,yaxis=1,ytitle='Relative Freq.',ycharsize=textsize,xthick=2,ythick=2,charthick=2,yrange=[0,1.05],ystyle=1;,yrange=[0,max(histdata)]
     ;plot 2D Distributions
     if (i lt j and j ne i) then begin
       cgLoadCT, colortable,/reverse,bottom=bottomct;33;65
          TVLCT, cgColor('white', /Triple), 0
          TVLCT, r, g, b, /Get
          palette = [ [r], [g], [b] ]
          ; interpolates scattered data values on a plane
            nbins=50
            x=transpose(LNESTPARS[i,*])
            y=transpose(LNESTPARS[j,*]) 
            z=(Prob)/max(prob) 
            b=(Highrange(i)-Lowrange(i))/nbins
            xout=findgen(nbins)*b+Lowrange(i)
            b=(Highrange(j)-Lowrange(j))/nbins
            yout=findgen(nbins)*b+Lowrange(j)
            ;result=griddata(x,y,z,/GRID,XOUT=xout,YOUT=yout) ; interpolates scattered data values on a plane
            ;cgcontour,result,xout,yout, NLevels=4, /Fill, Palette=palette,xrange=[Lowrange(i),Highrange(i)],yrange=[Lowrange(j),Highrange(j)],xtitle=parameternameslog(tofit(i)),ytitle=parameternameslog(tofit(j)),label=0,xthick=2,ythick=2,charthick=1,xcharsize=xchrsz,ycharsize=ychrsz,Position=pos[*,k], NoErase=j NE 0,/xstyle,/ystyle;,resolution=[res,res]

       cgcontour,(Prob)/max(prob),transpose(LNESTPARS[i,*]),transpose(LNESTPARS[j,*]), NLevels=NContourlevels, /Irregular,/fill, Palette=palette,xrange=[Lowrange(i),Highrange(i)],yrange=[Lowrange(j),Highrange(j)],xtitle=parameternameslog(tofit(i)),ytitle=parameternameslog(tofit(j)),label=0,xthick=2,ythick=2,charthick=1,xcharsize=xchrsz,ycharsize=ychrsz,Position=pos[*,k], NoErase=j NE 0,/xstyle,/ystyle,resolution=[res,res]
       cgcontour,(Prob)/max(prob),transpose(LNESTPARS[i,*]),transpose(LNESTPARS[j,*]), Levels=[1E-4,1E-3], /Irregular, /overplot,thick=2,resolution=[res,res],color=cgcolor('grey')
       cgcontour,(Prob)/max(prob),transpose(LNESTPARS[i,*]),transpose(LNESTPARS[j,*]), Levels=[1E-2,1E-1], /Irregular, /overplot,thick=2,resolution=[res,res],color=cgcolor('dark grey')
    endif
     xchrsz=0.001
     ychrsz=0.001
     k=k+1  
  ENDFOR
ENDFOR
   cgLoadCT, colortable,/reverse,bottom=bottomct;33;65
   cgColorbar, /Vertical, Range=[0,1 ], NColors=254-bottomct,charsize=1,charthick=2, Position=[0.96, 0.5, 0.98, 0.89],title='Relative Probability',bottom=bottomct
;   cgPS_Close
device,/close
!P.Multi = 0
  spawn,'ps2eps -B -C -f '+filenam+'.ps'
  spawn,'ps2pdf -dEPSCrop -dAutoRotatePages=/None -dAutoFilterColorImages=false -dColorImageFilter=/FlateEncode -dUseFlateCompression=true '+filenam+'.eps '+dirs+filenam+'.pdf'


   ;====C/O ratio when C and O are both fit
   filenam='arc_monitor_nested_corner_CO'
   cgPS_Open, dirs+filenam+'.ps'
   cgDisplay
   Cindex=where(parameternameslog eq 'log(C/C$\tex_{sun}$)')
   Oindex=where(parameternameslog eq 'log(O/O$\tex_{sun}$)')
   if (Chemistry eq 3 and Cindex[0] ne -1 and Oindex[0] ne -1) then begin
      ci=where(tofit eq cindex[0])
      oi=where(tofit eq Oindex[0])
      C=(10^(LNESTPARS[ci[0],*])*3.162E-04)  ;C VMR
      O=(10^(LNESTPARS[oi[0],*])*5.754E-04)  ;O VMR
      C_O=C/O
      srt=sort(C_O)
      parameter=    C_O(srt)
      parameter_prob=Prob[srt]
      medpar=parameter[value_locate(total(parameter_prob,/cumulative),[0.5d0])]
      e=parameter[value_locate(total(parameter_prob,/cumulative),[0.5d0-0.68268d0/2.,0.5d0+0.68268d0/2.])]-medpar[0]
      e2=parameter[value_locate(total(parameter_prob,/cumulative),[0.5d0-0.95449d0/2.,0.5d0+0.95449d0/2.])]-medpar[0]
      e3=parameter[value_locate(total(parameter_prob,/cumulative),[0.5d0-0.99730d0/2.,0.5d0+0.99730d0/2.])]-medpar[0]
      print,'(C/O): '+string(medpar)+' 1-sig:'+string(e(0))+string(e(1))+'  2-sig:'+string(e2(0))+string(e2(1))+'  3-sig:'+string(e3(0))+string(e3(1))
      printf,1,'(C/O): '+string(medpar)+' 1-sig:'+string(e(0))+string(e(1))+'  2-sig:'+string(e2(0))+string(e2(1))+'  3-sig:'+string(e3(0))+string(e3(1))
     ; print,'log(C/O): '+string(exp(medpar))+' 1-sig:'+string(e(0))+string(e(1))+'  2-sig:'+string(e2(0))+string(e2(1))+'  3-sig:'+string(e3(0))+string(e3(1))
        xrange=[medpar+e3(0)*1.5,medpar+e3(1)*1.5]
        ;hw=where(parameter_prob gt 0.0001)
        a=lin_rebin(parameter,parameter_prob,binwidth=(xrange(1)-xrange(0))/nhistobinsCO,bin_x=bin_x) ;Rebin historgram 
        cgLoadCT, 0;33;65
        cgplot,bin_x,a/max(a),psym=10,xrange=xrange,yrange=[0,1.05],xthick=2,ythick=2,charthick=2,xcharsize=1,ycharsize=1,/xstyle,/ystyle,xtitle='C/O ratio',ytitle='Relative Frequency'
        for z=10,1000 do cgplot,bin_x,a/max(a)*(z*0.001),psym=10,color=cgcolor('light grey'),thick=0.1,/overplot ;fill in the plot
        cgplot,bin_x,a/max(a),psym=10,color=cgcolor('black'),thick=3 ,/overplot
        cgplot,[medpar,medpar],[0,1.05],/overplot,thick=3, color=cgcolor("Orange",!D.Table_Size-2)
        cgplot,[e[0]+medpar,e[0]+medpar],[0,1.05], color=cgcolor("Orange",!D.Table_Size-2),/overplot,thick=3,linestyle=1
        cgplot,[e[1]+medpar,e[1]+medpar],[0,1.05], color=cgcolor("Orange",!D.Table_Size-2),/overplot,thick=3,linestyle=1
        cgaxis,/xaxis,xrange=xrange,xcharsize=0.001,xthick=2,/xstyle
        cgplot,bin_x,a/max(a),psym=10,xrange=xrange,yrange=[0,1.05],xthick=2,ythick=2,charthick=2,xcharsize=xchrsz,ycharsize=0.001,/xstyle,/ystyle,/nodata,/overplot
        ;cgaxis,yaxis=1,ytitle='Relative Freq.',ycharsize=textsize,xthick=2,ythick=2,charthick=2,yrange=[0,1.05],ystyle=1
   device,/close
   spawn,'ps2eps -B -C -f '+filenam+'.ps'
   spawn,'ps2pdf -dEPSCrop -dAutoRotatePages=/None -dAutoFilterColorImages=false -dColorImageFilter=/FlateEncode -dUseFlateCompression=true '+filenam+'.eps '+dirs+filenam+'.pdf'
   endif
   !P.Multi = 0

   ;===Mass
   if (Emission eq 'DI') then begin
     filenam='arc_monitor_nested_corner_Mj'
     cgPS_Open, dirs+filenam+'.ps' & cgDisplay

     loggindex=where(parameternameslog eq 'log(g)')
     Rpindex=where(parameternameslog eq 'R$\tex_{pl}$(R$\tex_{Jup}$,1 mbar)')

      Gr=6.67428E-8 ;cgs
      gi=where(tofit eq loggindex[0])
      Ri=where(tofit eq Rpindex[0])
      grav= 10^(LNESTPARS[gi[0],*])  ; g
      R=    LNESTPARS[Ri[0],*]   ;O VMR
      Mass=(grav*(R*71492D5)^2./Gr)/1.898E30  ;MJup
      srt=sort(Mass)
      parameter=  Mass(srt)
      parameter_prob=Prob[srt]
      medpar=parameter[value_locate(total(parameter_prob,/cumulative),[0.5d0])]
      e=parameter[value_locate(total(parameter_prob,/cumulative),[0.5d0-0.68268d0/2.,0.5d0+0.68268d0/2.])]-medpar[0]
      e2=parameter[value_locate(total(parameter_prob,/cumulative),[0.5d0-0.95449d0/2.,0.5d0+0.95449d0/2.])]-medpar[0]
      e3=parameter[value_locate(total(parameter_prob,/cumulative),[0.5d0-0.99730d0/2.,0.5d0+0.99730d0/2.])]-medpar[0]
      print,'Mass (MJup) : '+string(medpar)+' 1-sig:'+string(e(0))+string(e(1))+'  2-sig:'+string(e2(0))+string(e2(1))+'  3-sig:'+string(e3(0))+string(e3(1))
      printf,1,'Mass (MJup) : '+string(medpar)+' 1-sig:'+string(e(0))+string(e(1))+'  2-sig:'+string(e2(0))+string(e2(1))+'  3-sig:'+string(e3(0))+string(e3(1))
        xrange=[medpar+e3(0)*1.5,medpar+e3(1)*1.5]
        ;hw=where(parameter_prob gt 0.0001)
        a=lin_rebin(parameter,parameter_prob,binwidth=(xrange(1)-xrange(0))/nhistobinsCO,bin_x=bin_x) ;Rebin historgram 
        cgLoadCT, 0;33;65
      cgplot,bin_x,a/max(a),psym=10,xrange=xrange,yrange=[0,1.05],xthick=2,ythick=2,charthick=2,xcharsize=1,ycharsize=1,/xstyle,/ystyle,xtitle='Mass (M$\tex_{Jup}$)',ytitle='Relative Frequency'
        for z=10,1000 do cgplot,bin_x,a/max(a)*(z*0.001),psym=10,color=cgcolor('light grey'),thick=0.1,/overplot ;fill in the plot
        cgplot,bin_x,a/max(a),psym=10,color=cgcolor('black'),thick=3 ,/overplot
        cgplot,[medpar,medpar],[0,1.05],/overplot,thick=3, color=cgcolor("Orange",!D.Table_Size-2)
        cgplot,[e[0]+medpar,e[0]+medpar],[0,1.05], color=cgcolor("Orange",!D.Table_Size-2),/overplot,thick=3,linestyle=1
        cgplot,[e[1]+medpar,e[1]+medpar],[0,1.05], color=cgcolor("Orange",!D.Table_Size-2),/overplot,thick=3,linestyle=1
        cgaxis,/xaxis,xrange=xrange,xcharsize=0.001,xthick=2,/xstyle
        cgplot,bin_x,a/max(a),psym=10,xrange=xrange,yrange=[0,1.05],xthick=2,ythick=2,charthick=2,xcharsize=xchrsz,ycharsize=0.001,/xstyle,/ystyle,/nodata,/overplot
        ;cgaxis,yaxis=1,ytitle='Relative Freq.',ycharsize=textsize,xthick=2,ythick=2,charthick=2,yrange=[0,1.05],ystyle=1
   
     device,/close
     spawn,'ps2eps -B -C -f '+filenam+'.ps'
     spawn,'ps2pdf -dEPSCrop -dAutoRotatePages=/None -dAutoFilterColorImages=false -dColorImageFilter=/FlateEncode -dUseFlateCompression=true '+filenam+'.eps '+dirs+filenam+'.pdf'
     !P.Multi = 0
   endif

   
close,1
end

 ; nbins=20
 ; x=transpose(LNESTPARS[i,*])
 ; y=transpose(LNESTPARS[j,*]) 
 ; z=(Prob)/max(prob) 
 ; b=(Highrange(i)-Lowrange(i))/nbins
 ; xout=findgen(nbins)*b+Lowrange(i)
 ; b=(Highrange(j)-Lowrange(j))/nbins
 ; yout=findgen(nbins)*b+Lowrange(j)

 ; cgLoadCT, 0,/reverse
 ; result=griddata(x,y,z,/GRID,XOUT=xout,YOUT=yout)                       
 ; cgcontour,result,xout,yout,label=0,/fill        
