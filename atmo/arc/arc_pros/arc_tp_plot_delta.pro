pro arc_tp_plot_delta,dirs,R,Rs
  Common SystemParameters,rplanet,rstar,rorbit,stellarfile,Tstar,Cloudtype,Hazetype,Chemistry,Isothermal,Opacity_switch,Abund_fit_switch,Abund_fit_index,Abund_priorLow,Abund_priorUP,rain_switch,Element_fit_switch,Element_fit_index,Ele_priorLow,Ele_priorUP,Distance,nband

  ;arc_load_system_parameters,planet,Rp,Rjup,Rpl,Teff,logg,g,Rstar,Rs,rorbit,stellarfile,Tstar
  dirsnow=dirs
  save,filename='tmp.sav',dirsnow,dirs ;
  
   Rjup  = 71492D5         ; cm
   Rstar=  Rs/6.960D10
   !p.font=0 
   filenam='TP.plot'
   Rmin=min(r)
   Rmax=max(r)  ; close to 1mbar in best-fit

       Result = FILE_TEST(dirsnow+'arc_tp_plot.sav')  
       if (Result eq 1) then begin restore,dirsnow+'arc_tp_plot.sav' & restore,'tmp.sav' & Result=1 & endif
       if (Result eq 1) then goto,skipload


set_plot,'x'
restore,dirs+'MPFIT.bestfit.sav' ;params Array[11], Kir=params(6)
restore,dirs+'MCMC.monitor.sav' ;allchains,allchainschi2,i,j,nchains,nthin,tofit
nfit=n_elements(tofit) 
dirs=dirsnow

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

  ;just keep good chains
   include=good ;[0,1,2,3,4,5,6,7,8,9,10,11] ;ALL
   trimchain=nchains-n_elements(include)
   allchains=allchains[*,*,include]
   allchainschi2=allchainschi2[*,include]
   j=j-trimchain
   nchains=nchains-trimchain



;Number of bins in Histogram
histobins=28;100;28 ;45

nfit=n_elements(tofit) 
parameternamessimple=['Temp ','lnZ  ','lnC/O','pp_Rp','Cloud','Haze','H2','He','H2O','CO2','CO','CH4','NH3','Na','K','Li','Rb','Cs','TiO','VO','FeH','PH3','H2S','HCN','C2H2']

parameternames=['Temperature (K)','Metallicity ln(Z/Z$\tex_\odot$)','C/O Ratio ln(C/O)','R$\tex_{pl}$(R$\tex_{Jup}$,1 bar)','cloud opacity ln($\sigma$/$\sigma$$\tex_0$)','Haze Opacity ln($\sigma$/$\sigma$$\tex_0$-1)','VMR ln(H$\tex_2$)','VMR ln(He)','VMR ln(H$\tex_2$O)','VMR ln(CO$\tex_2$)','VMR ln(CO)','VMR ln(CH$\tex_4$)','VMR ln(NH$\tex_3$)','VMR ln(Na)','VMR ln(K)','VMR ln(Li)','VMR ln(Rb)','VMR ln(Cs)','VMR ln(TiO)','VMR ln(VO)','VMR ln(FeH)','VMR ln(PH$\tex_3$)','VMR ln(H$\tex_2$S)','VMR ln(HCN)','VMR ln(C$\tex_2$H$\tex_2$)']


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

print,i,j
print,'Burnindex =',burnndx

b=allchains(*,burnndx:i,*)
a=(n_elements(b(0,*,0)))*1D0*nchains & a=long(a) ;chain length
c=reform(b,nfit,a)
chi=reform(ALLCHAINSCHI2(burnndx:i,*),a)

thinnum=floor(a/300.)
if thinnum lt 10 then thinnum=n_elements(a)-1
print,a,thinnum



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
if (nfit ge 16) then p14=transpose(c(15,*))
if (nfit ge 17) then p14=transpose(c(16,*))
if (nfit ge 18) then p14=transpose(c(17,*))

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
if (nfit eq 16) then pall=[[p0],[p1],[p2],[p3],[p4],[p5],[p6],[p7],[p8],[p9],[p10],[p11],[p12],[p13],[p14],[p15]]
if (nfit eq 17) then pall=[[p0],[p1],[p2],[p3],[p4],[p5],[p6],[p7],[p8],[p9],[p10],[p11],[p12],[p13],[p14],[p15],[p16]]
if (nfit eq 18) then pall=[[p0],[p1],[p2],[p3],[p4],[p5],[p6],[p7],[p8],[p9],[p10],[p11],[p12],[p13],[p14],[p15],[p16],[p17]]

loggArr=dblarr(a)   & loggArr(*)=params(2)
RplArr=dblarr(a)    & RplArr(*)=params(3)
KirArr=dblarr(a)    & kirarr(*)=params(11)
gamma1Arr=dblarr(a) & gamma1Arr(*)=params(12)
gamma2Arr=dblarr(a) & gamma2Arr(*)=params(13)
betaArr=dblarr(a)   & betaArr(*)=params(14)
alphaArr=dblarr(a)  & alphaArr(*)=params(15)
TintArr=dblarr(a)   & TintArr(*)=params(16)
deltaArr=dblarr(a)  & deltaArr(*)=params(17)
PdeltaArr=dblarr(a) & PdeltaArr(*)=params(18)

a=where(tofit eq 2)  & if (a ne -1) then loggArr=pall(*,a)
a=where(tofit eq 3)  & if (a ne -1) then RplArr=pall(*,a)
a=where(tofit eq 11) & if (a ne -1) then KirArr=pall(*,a)
a=where(tofit eq 12) & if (a ne -1) then gamma1Arr=pall(*,a)
a=where(tofit eq 13) & if (a ne -1) then gamma2Arr=pall(*,a)
a=where(tofit eq 14) & if (a ne -1) then betaArr=pall(*,a)
a=where(tofit eq 15) & if (a ne -1) then alphaArr=pall(*,a)
a=where(tofit eq 16) & if (a ne -1) then TintArr=pall(*,a)
a=where(tofit eq 17) & if (a ne -1) then deltaArr=pall(*,a)
a=where(tofit eq 18) & if (a ne -1) then PdeltaArr=pall(*,a)


KirArr=10^KirArr
gamma1Arr=10^gamma1Arr
gamma2Arr=10^gamma2Arr
PdeltaArr=10^PdeltaArr

;Open ATMO NCDF file of P-T profile
ptfile=dirs+'pt_ARC_start.ncdf'
fid=ncdf_open(ptfile)
ncdf_varget,fid,'pressure',pressure ;cgs Barye, bar is pressure*1D-6
p=alog10(pressure*1D-6)

arc_TP_delta,ptfile,params(2),10^params(11),10^params(12),10^params(13),params(14),params(15),T,params(16),press,params(17),10^params(18)
;arc_tp_delta,ptfile,logg,Kir,gamma1,gamma2,beta,alpha,T,Tint,Press,delta,Pdelta
            ;ptfile,Kir          ,gamma1       ,gamma2       ,beta      ,alpha     ,T,Tint     ,Press    ,delta,Pdelta
Tbest=T ;best-fit Temperature profile

;plot,T,pressure*1D-6,/ylog,yrange=[10,1E-6],xrange=[500,1500]
Tarr=T

for i=0,n_elements(KirArr)-1 do begin & $
  arc_TP_delta,ptfile,loggArr,KirArr(i),gamma1Arr(i),gamma2Arr(i),betaArr(i),alphaArr(i),T,TintArr(0),press,deltaArr(i),PdeltaArr(i) & $
  Tarr=[[Tarr],[T]] & $
endfor
  ;oplot,Tbest,pressure*1D-6,color=1000

T1siglow=Tbest
T1sighigh=Tbest
T2siglow=Tbest
T2sighigh=Tbest
T3siglow=Tbest
T3sighigh=Tbest
Tmean=Tbest
for i=0,n_elements(p)-1 do begin  & $
    e=percentiles(Tarr(i,*)) & $;-mean(Tarr(i,*))  & $
    e2=percentiles(Tarr(i,*),CONFLIMIT=0.95449) & $;-mean(Tarr(i,*))  & $
    e3=percentiles(Tarr(i,*),CONFLIMIT=0.9973002) & $;-mean(Tarr(i,*))  & $
    Tmean(i)=median(Tarr(i,*)) & $
    T1siglow(i)=e(0) & $
    T1sighigh(i)=e(1) & $
    T2siglow(i)=e2(0) & $
    T2sighigh(i)=e2(1) & $
    T3siglow(i)=e3(0) & $
    T3sighigh(i)=e3(1) & $
endfor
;help
skipload:

xmax=max(Tarr)
if (xmax gt 3500) then xmax=3500
set_plot,'ps'
  !p.font=0
  device,filename=dirs+filenam+'.ps',/helvetica,/landscape,/color
  cgplot,Tbest,p,xrange=[min(Tarr)-100,max(Tarr)+100],yrange=[2,-6],/xsty,thick=2.0,/ysty,/nodata,$
    ytitle=textoidl('log_{10}[Pressure (bar)]'),xtitle=textoidl('Temperature (K)'),$
    title='',charthick=4,xthick=4,ythick=4,ytickinterval=1,yminor=1
  for i=0,n_elements(KirArr)-1,thinnum do oplot,Tarr(*,i),p,thick=1  ,color=cgcolor('grey', !D.Table_Size-2)
  ;oplot,Tbest,p,thick=5,color=cgcolor('red', !D.Table_Size-2)



    K     = 1.380658D-16    ; gm*cm^2/s^2 * Kel
    mu    = 1.6726D-24*2.3  ; cgs
    Tau   = 0.56D0          ; Optical depth where Atmo Pressure corresponds to measured Transit Radius
 ;   result=min(abs(p+3),index) ;get 1mbar array index 
 ;     Rmbar=params(3)*rjup/rs
 ;     Tmbar=Tbest(index)
 ;     Hmbar=Tbest(index)*K/(mu*g)/rs ;units Ra
 ;     Pmin=exp(-(rmin-rmbar)/Hmbar)*0.001d0
 ;     Pmax=exp(-(rmax-rmbar)/Hmbar)*0.001d0
 ;     Pave=(alog10(Pmin)+alog10(Pmax))/2.
 ;     result=min(abs(p-Pave),index) ;get mid-pressure array index, update 
 ;     Tmbar=Tbest(index)
 ;     Hmbar=Tbest(index)*K/(mu*g)/rs ;units Ra
 ;     Pmin=exp(-(rmin-rmbar)/Hmbar)*0.001d0
 ;     Pmax=exp(-(rmax-rmbar)/Hmbar)*0.001d0
 ;     Contributionbest=where((pressure*1D-6 gt Pmax) and (pressure*1D-6 lt Pmin))
 ;     ;oplot,Tbest(Contributionbest),p(Contributionbest),thick=5,color=cgcolor('red6', !D.Table_Size-2)

 ;   for i=0,n_elements(KirArr)-1,thinnum do begin
 ;     Rmbar=RplArr(i)*rjup/rs
 ;     Tmbar=Tarr(index,i)
 ;     Hmbar=Tarr(index,i)*K/(mu*g)/rs ;units Ra
 ;     Pmin=exp(-(rmin-rmbar)/Hmbar)*0.001d0
 ;     Pmax=exp(-(rmax-rmbar)/Hmbar)*0.001d0
 ;     Pave=(alog10(Pmin)+alog10(Pmax))/2.
 ;     result=min(abs(p-Pave),index) ;get mid-pressure array index, update 
 ;     Tmbar=Tbest(index)
 ;     Hmbar=Tbest(index)*K/(mu*g)/rs ;units Ra
 ;     Pmin=exp(-(rmin-rmbar)/Hmbar)*0.001d0
 ;     Pmax=exp(-(rmax-rmbar)/Hmbar)*0.001d0
 ;     Contribution=where((pressure*1D-6 gt Pmax) and (pressure*1D-6 lt Pmin))
 ;     if (Contribution(0) ne -1) then oplot,Tarr(Contribution,i),p(Contribution),thick=6,color=cgcolor('sky blue', !D.Table_Size-2)
 ;   endfor


    

    press=p
    Fe_H=0.0                    ;log10[Fe/H]
    TempM_Hoffset=500
    condthick=3

                                ;References Vischer 2006,2010, Wakeford 2017
                                ;Ca
    Tcond_Ca3Ti2O7  = 1D4/(5.114 - 0.284*press - 0.568*Fe_H)
    Tcond_Ca4Ti3O10 = 1E4/(5.113 - 0.281*press - 0.563*Fe_H)
    Tcond_CaTiO3    = 1E4/(5.125 - 0.277*press - 0.554*Fe_H)
    oplot,Tcond_Ca3Ti2O7,press,linestyle=2,color=cgcolor('dodger blue', !D.Table_Size-2) ,thick=condthick
    oplot,Tcond_Ca4Ti3O10,press,linestyle=2,color=cgcolor('dodger blue', !D.Table_Size-2) ,thick=condthick
    oplot,Tcond_CaTiO3,press,linestyle=2,color=cgcolor('dodger blue', !D.Table_Size-2) ,thick=condthick
    if (min(Tarr)-100 lt 2000) then xyouts,2015,0.8,'Ca-Ti',ORIENTATION=-60,font=1,charsize=1,color=cgcolor('dodger blue', !D.Table_Size-2)

                                ;Al
    Tcond_Al2O3     = 1E4/(5.014 - 0.2179*press + 2.264D-3*press^2d0 - 0.580*Fe_H)
    oplot,Tcond_Al2O3,press,linestyle=2,color=cgcolor('purple', !D.Table_Size-2) ,thick=condthick
    if (min(Tarr)-100 lt 2100) then xyouts,2100,1,textoidl('Al_2O_3'),ORIENTATION=-60,font=1,charsize=1,color=cgcolor('purple', !D.Table_Size-2)

                                ;Fe
    Tcond_Fe         =1E4/(5.44 - 0.48*press - 0.48*Fe_H)
    oplot,Tcond_Fe,press,linestyle=2,color=cgcolor('dark red', !D.Table_Size-2) ,thick=condthick
    if (min(Tarr)-100 lt 1950) then xyouts,1980,1,'Fe',ORIENTATION=-40,font=1,charsize=1,color=cgcolor('dark red', !D.Table_Size-2)

                                ;Mg
    Tcond_Mg2SiO4   = 1E4/(5.89  - 0.37*press - 0.73*Fe_H)
    Tcond_MgSiO3    = 1E4/(6.26  - 0.35*press - 0.70*Fe_H)
    oplot,Tcond_Mg2SiO4,press,linestyle=2,color=cgcolor('red5', !D.Table_Size-2) ,thick=condthick
    if (min(Tarr)-100 lt 1800) then xyouts,1830,1,textoidl('Mg_2SiO_4'),ORIENTATION=-60,font=1,charsize=1,color=cgcolor('red5', !D.Table_Size-2)

    oplot,Tcond_MgSiO3,press,linestyle=2,color=cgcolor('red', !D.Table_Size-2) ,thick=condthick
    if (min(Tarr)-100 lt 1700) then xyouts,1700,1,textoidl('MgSiO_3'),ORIENTATION=-60,font=1,charsize=1,color=cgcolor('red', !D.Table_Size-2)

                                ;MgS
    Tcond_MgS       = 1D4/(7.45  - 0.42*press - 0.841*Fe_H)
    oplot,Tcond_MgS,press,linestyle=2,color=cgcolor('Tomato', !D.Table_Size-2) ,thick=condthick
    if (min(Tarr)-100 lt 1350) then xyouts,1395,1,textoidl('MgS'),ORIENTATION=-65,font=1,charsize=1,color=cgcolor('Tomato', !D.Table_Size-2)

                                ;Na2S
    Tcond_Na2S      = 1D4/(10.5-0.72*press-1.08*Fe_H)
    oplot,Tcond_Na2S,press,linestyle=2,color=cgcolor('Dark Goldenrod', !D.Table_Size-2) ,thick=condthick
    if (min(Tarr)-100 lt 1000) then xyouts,1030,1,textoidl('Na_2S'),ORIENTATION=-60,font=1,charsize=1,color=cgcolor('Dark Goldenrod', !D.Table_Size-2)
                                ;KCL
    Tcond_KCL       = 1D4/(12.479 -0.879*press -0.879*Fe_H)
    oplot,Tcond_KCL,press,linestyle=2,color=cgcolor('dark green', !D.Table_Size-2) ,thick=condthick
    if (min(Tarr)-100 lt 870) then xyouts,870,1,textoidl('KCL'),ORIENTATION=-60,font=1,charsize=1,color=cgcolor('dark green', !D.Table_Size-2)

                                ;ZnS
    Tcond_ZnS       = 1D4/(12.52 - 0.63*press - 1.26*Fe_H)
    oplot,Tcond_ZnS,press,linestyle=2,color=cgcolor('dark grey', !D.Table_Size-2) ,thick=condthick
    if (min(Tarr)-100 lt 800) then xyouts,810,1,textoidl('ZnS'),ORIENTATION=-70,font=1,charsize=1,color=cgcolor('dark grey', !D.Table_Size-2)

                                ;H2O
    Tcond_H2O       = 1D4/(38.84 - 3.93*Fe_H - 3.83*press - 0.20*Fe_H*press)
    oplot,Tcond_H2O,press,linestyle=2,color=cgcolor('blue', !D.Table_Size-2) ,thick=condthick
    if (min(Tarr)-100 lt 300) then xyouts,300,1,textoidl('H_2O'),ORIENTATION=-80,font=1,charsize=1,color=cgcolor('blue',!D.Table_Size-2)
    

    



    
   smoothparam=1
   TmeanSm=smooth(tmean,smoothparam)
   ; plot 1 sigma profiles
   oplot,smooth(t1siglow,smoothparam),p,color=cgcolor('blue', !D.Table_Size-2),linestyle=2,thick=9
   oplot,smooth(tmean,smoothparam),p,color=cgcolor('blue', !D.Table_Size-2),thick=9              
   ;oplot,TmeanSm(Contributionbest),p(Contributionbest),color=cgcolor('blue', !D.Table_Size-2),thick=25             
   oplot,smooth(t1sighigh,smoothparam),p,color=cgcolor('blue', !D.Table_Size-2),linestyle=2,thick=9
   oplot,smooth(t2sighigh,smoothparam),p,color=cgcolor('BLU4', !D.Table_Size-2),linestyle=2,thick=9
   oplot,smooth(t2siglow,smoothparam),p,color=cgcolor('BLU4', !D.Table_Size-2),linestyle=2,thick=9
   oplot,smooth(t3sighigh,smoothparam),p,color=cgcolor('BLU2', !D.Table_Size-2),linestyle=2,thick=9
   oplot,smooth(t3siglow,smoothparam),p,color=cgcolor('BLU2', !D.Table_Size-2),linestyle=2,thick=9


;plot best fit
  ;oplot,Tbest,p,thick=6,color=cgcolor('red', !D.Table_Size-2)
    ;oplot,Tbest(Contributionbest),p(Contributionbest),thick=25,color=cgcolor('red6', !D.Table_Size-2)

;plot H26 isothermal best fit
  ;plotsym,0,2,/fill
  ;oplot,[958.44223],[-2.24],psym=8
  ;oploterror,958.44223,-2.24,114.96869,0.6,/lobar,thick=5,/nohat
  ;oploterror,958.44223,-2.24,126.73626,0.5,/hibar,thick=5,/nohat


  device,/close

  spawn,'ps2pdf -dEPSCrop -dAutoRotatePages=/None -dAutoFilterColorImages=false -dColorImageFilter=/FlateEncode -dUseFlateCompression=true '+dirs+filenam+'.ps '+dirs+filenam+'.pdf'
  spawn,'rm '+dirs+filenam+'.ps'
  save,filename=dirs+'arc_tp_plot.sav',/var,/COMPRESS
  spawn,'rm tmp.sav'
end
