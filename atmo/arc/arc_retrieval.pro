pro arc_retrieval,input
; ATMO Retrieval Code
; PURPOSE:
;  Perform model fit retrieval techniques to ATMO model using DEMC,
;  Nested Sampling
;
; DESCRIPTION:
;  Uses ATMO model Fit Transmission Spectra for Abundances.
;  ATMO currently uses 19 opacity sources: 
;    1) H2   2) He   3) H2O  4) CO2  5) CO  6) CH4  7) NH3  
;    8) Na   9) K   10) Li  11) Rb  12) Cs  13) TiO  14) VO  
;   15) FeH 16) PH3 17) H2S 18) HCN 19) C2H2
;
;  ATMO equilbrium chemistry tracks 166 species 
;
;  Fit starts with Chi^2 minimzation with a Levenberg-Marqart routine (mpfit) 
;  Fit then does (optional) a differential evolution chain monte carlo fit (Eastman et al. 2014)
;  and/or a nested sampling routine to meaure the evidence
;
;  Temperature is either Isothermal, or can use parameterized T-P
;  profile from Gillot 2010, consisiting of 5 free parameters:
;     Kir    = Planck mean thermal IR opacity (cm^2 g^-1)
;     gamma1 = ratios of the Planck mean opacities in the visible to thermal
;     gamma2 = ratios of the Planck mean opacities in the visible to thermal
;     beta   = albedo/redistribution factor
;     alpha  = partitions opt/ir flux
;
;
; STEPS: 
;   0) EDIT arc_load_system_parameters.pro and enter Planet data if needed
;   1) INPUT PARAMETES into arc_retrieval_input_EXONAME.pro 
;      a) Sets Saving directory, Saving Filename, and DATA file 
;      b) Sets Starting Parameters 
;      c) Sets Switches 
;   4) Save arc_retrieval.pro with Desired Inputs to be run
;   5) Compile and Run  
;      IDL> .r arc_retrival
;      IDL> arc_retrieval,'arc_retrieval_input_EXONAME'
;
;   Parameters to set in arc_retrieval_input_EXONAME.pro :
;     dirs,filenam,datafile                                                      (File setup)
;     Planet,DoMCMC,nthin,nchains,maxsteps,wsdstart,Cloudtype,Hazetype,Chemistry (ARC Optional Switches)
;     FIXIT, StartingParams                                                      (Set fixed parameters)
;     Opacity_switch, Abund_fit_switch                                           (Setup Chemistry)
;
;   Note: ARC will skip MPFIT and MCMC_step if previouly ran and files exist. 
;         Delete MPFIT.bestfit.sav and/or MCMC.scale.sav if new clean run is desired.
;
;   Transmission Data File requires 4 column ascii file.
;      wavelength (microns) | wavelenth bin (microns) | radius ratio | radius ratio error
;      saved as IDL varables: ws, werr, r, rerr
;
;   Fit is done at low resolution (500 wavenumber bands; R~140 at 1um)
;   Plots are made at medium resolution (5000 wavenumber bands; R~1400 at 1um)
;
;   PRIORS:
;   Physical bounding box priors:  
;      Bounding box physically motivated priors.  In MCMC, if random step
;      is outside of these bounds, chi^2 is set to (near) infinity so
;      step will always be rejected.  Priors are passed as an array   
;         Temperature > 200    (temperature must be > 200 K)         [200,3500]
;         Metallicity < 10^4.5 = 31,622x solar (Max value for ATMO)  [-8,4.50]
;         Metallicity > lower limit may be needed, currently 1D-9 from formatting input into arc_make_atmo_in.pro
;         C/O ratio   > 10^-5  (must be positive number, lower practical bound)  [exp(1D-5),exp(1D+5)]
;         C/O ratio   < 10^+5  (must be positive number, upper practical bound)-> may need adjusting, Oxygen is derived from C/O and Metallicity
;         Abundances  > 10^-12 (lower practical bound for fit abundances)[-12,1.3]
;         Abundances  < 1000   (upper practical bound, must be number less than 1000)
;         Kir         < 10^2
;         Kir         > 10^-5
;         Gamma1,Gamma2 < 10^2
;         Gamma1,Gamma2 > 10^-4
;         alpha is between 0 and 1
;         beta is between 0 and2
;      Futher Gaussian priors can be added, but are not yet implimented
;
;    Outputs:
;       ARC/*.pro                        - Directory containing IDL programs used in the fit
;       ARC.transpec.*.ncdf              - Netcdf files containing each piece/molecule of the transit spectrum.  
;                                          Note Vertical the offset is NOT applied
;       chem_eq_ARC_start00.ncdf         - Netcdf file with the Starting Chemistry 
;       pt_ARC_start.ncdf                - Netcdf file with the Starting P-T profile 
;       arc_monitor_distribution.ps      - Postscript file showing MCMC posterior cornerplots and histograms
;       FILENAM.highres.ARC.in           - ATMO namelist file of the best-fit model
;       FILENAM.highres.transitspec.ncdf - Netcdf files containing best-fit model
;       FILENAM.MCMC.FitPars.txt         - ascii file of MCMC fitted parameters
;       FILENAM.MCMC.FitPars.sav         - IDL save file of MCMC fitted parameters
;       FILENAM.MCMC.sav                 - IDL save file of MCMC posterior distribution
;       FILENAM.MPFIT.FitPars.sav        - IDL save file of L-M fitted parameters
;       FILENAM.sav                      - IDL save file of Data and fitted parameters
;       FILENAM.compare_pieces_Press.ps  - Postscript file showing each piece/molecule of the transit spectrum in pressure units
;       FILENAM.compare_pieces_Spec.ps   - Postscript file showing each piece/molecule of the transit spectrum in radius units
;       FILENAM.pdf                      - Bestfit-model and Data plotted
;       MCMC.monitor.FitPars.sav         - IDL save file of MCMC fitted parameters
;       MCMC.scale.sav                   - IDL save file of MCMC stepsizes 
;       MCMC.monitor.sav                 - IDL save file of current MCMC Mchain (Updated after every step)
;       MPFIT.bestfit.sav                - IDL save file of L-M fitted parameters
;
; DKS (9/3/2016)
; Version 3.0   12 Dec  2010  DEMC adjustments, log g, haze adjustments
; Version 2.10  14 Jun  2019  DEMC adjustments, cloud-top pressure added
; Version 2.00  14 Dec  2018  Element Variable and rainout
; Version 1.00  08 Mar  2017  Added NonEQ Chem, Kzz, Priors, Nested Sampling
;                             Priors & Nested Sampling need further testing to verifty robustness
; Version 0.92  02 Feb  2017  Fixed Emission spectra to not use Top Of Atmosphere Rpl, but rather user input Rpl 
; Version 0.91  14 Dec  2016  Now calculates Emission Spectra Contribution Functions during MPFIT run
; Version 0.90  06 Oct  2016  Now can do Emission Spectra
; Version 0.75  19 Aug  2016  P-T profile fit, now ARC
; Version 0.70  20 June 2016  'Free Chemistry' now added.
; Version 0.60  01 June 2016  Box Physical Priors added
; Version 0.55  25  May 2016  Opacities can now be switched, artificial haze included
; Version 0.50  19  May 2016  Chemical Equilibrum Version Operational
;
;================================================================
  ;================================================================
  ;Set Common Blocks
  Common ARC_retrieval_data, x,y,err,waverr
  Common SystemParameters,rplanet,rstar,rorbit,stellarfile,Tstar,Cloudtype,Hazetype,Chemistry,Isothermal,Opacity_switch,Abund_fit_switch,Abund_fit_index,Abund_priorLow,Abund_priorUP,rain_switch,Element_fit_switch,Element_fit_index,Ele_priorLow,Ele_priorUP,Distance,nband
  Common cpcount,ChildProcessCount,dirs,allmodel,allmodelpoints,Emission,nlevel,pmin,pmax,allchemistry,allTP
  common passmodel,ws,model,modelpoints,LowPriorsParams,UpPriorsParams,abundances,pressure,temperature

   set_plot,'x'
   Compile_Opt DEFINT32
  ;================================================================
   ChildProcessCount=0
  ;====================================================================================================================
  ; INPUTS
  ; All inputs are stored in a seperate file 
  ;=================
   boostmcmc=1 ;default
   restore,input; EleBest,Element_fit_switch,Ele_priorLow,Ele_priorUP,AbunBest,ABUND_FIT_SWITCH,CHEMISTRY,CLOUDTYPE,DATAFILE,DIRS,DOMCMC,EMISSION,FILENAM,FIXIT,HAZETYPE,$
                ; ISOTHERMAL,MAXSTEPS,NCHAINS,NLEVEL,NTHIN,OPACITY_NAMES,OPACITY_SWITCH,PARNAMEHIST,PLANET,PMAX,PMIN,RESTARTMCMC,$ 
                ; STARTINGPARAMS,TINT,TP_NAMES,TP_PARAMS,TP_SWITCH,WSDSTART,LowPriorsParams,UpPriorsParams,Abund_priorLow, Abund_priorUP,TP_priorLow,TP_priorUP,$
                ; NUM,SAMPLE_NUM,TOL,NTHREADS,EXPAND,PLOTIT,Rp,Rjup,Rpl,Teff,logg,g,Rstar,Rs,rorbit,stellarfile,Tstar,boostmcmc,monitor,DoMPFIT,RestartNest,nband
   if (OPACITY_SWITCH[21] eq 1) then hminus='yes'                    ; add H- opacity
   if (EMISSION ne 'DI') then Distance=0                             ; ignore distance unless doing Direct Imaging Spectra
   result=isa(RestartNest) & if (result eq 0) then RestartNest='no'  ; ignore if not defined
   result=isa(nband) & if (result eq 0) then nband=350               ; Default to 350 bands if not defined
   ;====================================================================================================================


   ;====================================================================================================================
   ; dirctory of library files
   !path = '../../arc/arc_pros/:'+ !path                                      ;remember to change arc_multinest.pro threads setup
   !path = '../../arc/coyote/:'+ !path                                        ;remember to change arc_multinest.pro threads setup
   !path = '../../arc/arc_multinest_pros/:'+ !path                            ;remember to change arc_multinest.pro threads setup
 
    if (Monitor eq 'yes') or (Monitor eq 'YES') or (Monitor eq 'Yes') then print,'** SKIPING MCMC, GOING TO MAKE PRELIMINARY PLOTS**'
    if (Monitor eq 'yes') or (Monitor eq 'YES') or (Monitor eq 'Yes') then goto,skip_arc_fitforwardmodel

   ;======
   ;                0                1              2               3                   4                 5                6               7              8                   9                10
   parnameHist2=  ['Teff'         ,'lnMdH'        ,'logg'         ,'Rp'            ,'lnCloud'      ,'ln(Haze)'     ,'cloud_top log(bar)','haze alpha'    ,'log10(Kzz)'       ,'fsed'          ,'sigma_g']
   parname=       ['Teff(K)      ','ln(MdH)      ','log(g)       ','Rp(RJup 1mbar)','ln(Cloud)    ','ln(haze)     ','cloud_top log(bar)','haze alpha    ','log10(Kzz cm^2s)' ,'fsed'          ,'sigma_g']
   if (Isothermal eq 0) then FIXIT(0)=1 ;Fix Isothermal Temperature if Free T-P profile
   if (Cloudtype eq 0) then FIXIT(4)=1  ;Don't Vary Clouds when there are none there. StartingParam[4] will be ignored
   if (Cloudtype eq 0) then FIXIT(6)=1  ;Don't Vary Cloud top when there are none there. StartingParam[6] will be ignored
   if (Cloudtype ne 3) then FIXIT(9)=1  ;Don't Vary Eddysed clouds unless turned on
   if (Cloudtype ne 3) then FIXIT(10)=1  ;Don't Vary Eddysed clouds unless turned on
   if (FIXIT(6) eq 1)  then StartingParams(6)=alog10(1.01*pmin)  ;Set minimum pressure of cloud to almost the top of model by default if not a free parameter
   if (Hazetype  eq 0) then FIXIT(5)=1  ;Don't Vary Haze when there is none there. StartingParam[5] will be ignored
   if (Hazetype  eq 0) then FIXIT(7)=1  ;Don't Vary Haze alpha when there is none there. StartingParam[7] will be ignored
  ;if (FIXIT(7)  eq 1) then StartingParams(7)=4.0d0  ;Default to Rayleigh Scattering if not Fit
   if (Chemistry eq 1) then FIXIT(1)=1  ;Don't Vary Metalcity when Fitting for individual Abundances
   ;if (Chemistry eq 1) then FIXIT(2)=1  ;Don't Vary C/O when Fitting for individual Abundances
   if (Chemistry eq 2) then FIXIT(9)=0  ;Vary Kzz when fitting non-equibrium chemsitry
   if (Chemistry ne 2) then FIXIT(9)=1  ;Don't Vary Kzz unless fitting non-equilibrium chemistry
   if (Cloudtype eq 2) then parnameHist=   ['Temperature (K)','Abundance ln(Z/Zs)','log(g)','Rpl(1bar)','cloud deck height','haze opacity','cloud_top log(bar)','haze alpha' ,'log10(Kzz)','fsed' ,'sigma_g']
   if (Cloudtype eq 2) then parnameHist2=  ['Teff'           ,'lnMdH'             ,'log(g)','Rp'       ,'Cloud_z'          ,'ln(Haze)'    ,'CloudTop(Rp)'      ,'haze alpha' ,'log10(Kzz)','fsed' ,'sigma_g']
   FIXIT=[FIXIT,TP_switch]
   ;if (Isothermal eq 0) then FIXIT(0)=1 ;Don't Vary internal temperature if Fitting for T-P profile
   if (Chemistry eq 1) then FIXIT=[FIXIT,Abund_fit_switch]
   ;if (Chemistry eq 3) then FIXIT=[FIXIT,Element_fit_switch]
   ; Isothermal       10   11   12    13    14   15    16    17  18  19   20   21   22   23     24   25   26     27     28      29    30   31             
   ; T-P fit          15   16   17    18    19   20    21    22  23  24   25   26   27   28     29   30   31     32     33      34    35   36 
   Opacity_names   =['H2','He','H2O','CO2','CO','CH4','NH3','Na','K','Li','Rb','Cs','TiO','VO','FeH','PH3','H2S','HCN','C2H2','SO2','Fe','H_'] ;Spectrally active molecules
   Abund_fit_index =[   1,  62,    9,   15,   4,   12,   68, 108,114, 156, 164, 160,  128, 130,  147, 170,   137,   76,    17,  122, 146, 177] ;Atmo Index (IDL index is -1 of this value)
  ;Abund_fit_index =[   1,  62,    9,   15,   4,   12,   68, 108,113, 147, 155, 151,  126, 128,  140, 161,   133,   76,    17,] ;Atmo Index (IDL index is -1 of this value)
   numabunds=n_elements(Abund_fit_index)-1
  
   ;Chemistry = 3 ; Elements (rather than molecules) themselves are varied,
   ;Element indicies
   ; Isothermal       10   11  12  13  14  15   16   17   18  19   20  21   22   23  24   25   26   27   28  29    30   31 32                   
   ; T-P fit          15   16  17  18  19  20   21   22   23  24   25  26   27   28  29   30   31   32   33  34    35   36 37

   ;Trim elements to only those that are fit
   if (Chemistry eq 3) then begin
     Element_names   =  ['H','He','C','N','O','Na','K','Si','Ar','Ti','V','S','Cl','Mg','Al','Ca','Fe','Cr','Li','Cs','Rb','F','P'] ;Spectrally active molecules
     Element_fit_index =[   1,  2,  3,  4,  5,   6,  7,   8,   9,  10, 11, 12,  13,  14,  15,  16,  17,  18,  19,  20,  21, 22, 23] ;Atmo Index (IDL index is -1 of this value)
     numabunds=where(Element_fit_switch eq 0) ;index of activly fit elements
     Element_names=Element_names(numabunds)
     Element_fit_index=Element_fit_index(numabunds)
     EleBest=EleBest(numabunds)
     Element_fit_switch=Element_fit_switch(numabunds)
     Ele_priorLow=Ele_priorLow(numabunds)
     Ele_priorUP=Ele_priorUP(numabunds)
     FIXIT=[FIXIT,Element_fit_switch]
     numabunds=n_elements(numabunds) ;number of fit elements
  endif

   ;add in T-P parameters 
   ;if (Isothermal eq 0) then parname=[parname,TP_names]
   ;if (Isothermal eq 0) then parnameHist=[parnameHist,TP_names]
   ;if (Isothermal eq 0) then parnameHist2=[parnameHist2,TP_names]
   parname=[parname,TP_names]
   parnameHist=[parnameHist,TP_names]
   parnameHist2=[parnameHist2,TP_names]
   ;add in element/abundances if fit
   if (Chemistry eq 1) then parnameHist=[parnameHist,Opacity_names(0:numabunds)]
   if (Chemistry eq 1) then parnameHist2=[parnameHist2,Opacity_names(0:numabunds)]
   if (Chemistry eq 1) then parname=[parname,Opacity_names(0:numabunds)+'            ']
   if (Chemistry eq 3) then parnameHist=[parnameHist,Element_names]
   if (Chemistry eq 3) then parnameHist2=[parnameHist2,Element_names]
   if (Chemistry eq 3) then parname=[parname,Element_names+'            ']
   ;
   ;if (Isothermal eq 0) then StartingParams(0) = Tint ;Set temperature to internal temperature if fitting for T-P profile
   ;if (Isothermal eq 0) then parname(0)='Tint'
   ;if (Isothermal eq 0) then parnameHist2(0)='Tint'
   ;if (Isothermal eq 0) then StartingParams =  [StartingParams,TP_params]    ;Set temperature to internal temperature if fitting for T-P profile
   ;if (Isothermal eq 0) then UpPriorsParams  = [UpPriorsParams,TP_priorUP]   ;Add temperature profile priors
   ;if (Isothermal eq 0) then LowPriorsParams = [LowPriorsParams,TP_priorLow] ;Add temperature profile priors
   StartingParams =  [StartingParams,TP_params]    ;Set temperature to internal temperature if fitting for T-P profile
   UpPriorsParams  = [UpPriorsParams,TP_priorUP]   ;Add temperature profile priors
   LowPriorsParams = [LowPriorsParams,TP_priorLow] ;Add temperature profile priors
   ;
   fx=where(FIXIT eq 1)
   fre=where(FIXIT eq 0)
   print,'Fixed params :',parnameHist2(fx)
   print,'Fitting for  :',parnameHist2(fre)
   ;================================================================
   ;Constants  
     K     = 1.380658D-16    ; gm*cm^2/s^2 * Kel
     mu    = 1.6726D-24*2.3  ; cgs
     Tau   = 0.56D0          ; Optical depth where Atmo Pressure corresponds to measured Transit Radius
     Rjup  = 71492D5         ; cm
   ;================================================================

   ;================================================================
   ;Create Needed Directory Structure and load in data
     MPFITfilenameprefix=filenam
     ;create directory
       Result = FILE_TEST(dirs, /DIRECTORY)
       if (Result eq 0) then begin 
         spawn,'mkdir '+dirs
       endif
     spawn,'cp '+datafile+' '+dirs
     arc_loaddata,dirs+datafile,wsdstart,wsd,r,rerr,werr   
     ;create directory for AIT (keep copy of what's run)
       Result = FILE_TEST(dirs+'ARC/', /DIRECTORY) & if (Result eq 0) then spawn,'mkdir '+dirs+'ARC/'
       spawn,'cp arc_*.pro '+dirs+'ARC/'
     ;arc_load_system_parameters,planet,Rp,Rjup,Rpl,Teff,logg,g,Rstar,Rs,rorbit,stellarfile,Tstar
     rplanet=Rpl ;Radius of the planet in cm.  Used for Emission Flux calculation, fnu/fstar*(rplanet/rs)^2.
     logg=StartingParams(2)
     if (Emission eq 'yes') then StartingParams(3)=rplanet/71492D5 ;Force Fitted Planet Radius to be Same as Emission-used value in case it wasn't set the same
     ;initialize arrays to save models during MCMC

     allmodel = dblarr(nband,maxsteps,nchains)  ;350 wavelentgh points at lowres resolution 
     ;if (Emission eq 'yes')  then allmodel = dblarr(150,maxsteps,nchains)  ;350 wavelentgh points at lowres resolution 
     ;if (Emission eq 'DI')  then allmodel = dblarr(1100,maxsteps,nchains)  ;350 wavelentgh points at lowres resolution 
     ;allmodelpoints = dblarr(n_elements(wsd),maxsteps,nchains) 
     ;allchemistry   = dblarr(nlevel,n_elements(Opacity_switch),maxsteps,nchains)
     ;allTP          = dblarr(nlevel,nlevel,maxsteps,nchains)
   ;================================================================

   ;================================================================
   ;================================================================
   ;DO FIT 
   ;Starting Paratmeters for fit
     ChildProcessCount=0
     A = StartingParams     ; e.g. z0, Teq, Tio,Vo,FeH 
     ;
     if (Isothermal eq 1) then arc_ptmaker,nlevel,pmin,pmax
     if (Isothermal eq 0) then begin
         ;run T-P profile
         arc_ptmaker,nlevel,pmin,pmax
         spawn,'cp input_TP_arc.ncdf pt_ARC'+string(ChildProcessCount,format='(I02)')+'.ncdf'
         ptfile='pt_ARC'+string(ChildProcessCount,format='(I02)')+'.ncdf'
         Kir=10^A(11) & gamma1=10^A(12) & gamma2=10^A(13) & beta=A(14)  & alpha=A(15) & Tint=A(16) & delta=A(17) & Pdelta=10^A(18) 
         arc_tp_delta,ptfile,logg,Kir,gamma1,gamma2,beta,alpha,T,Tint,Press,delta,Pdelta
         spawn,'cp pt_ARC00.ncdf '+dirs+'pt_ARC_start.ncdf'
      endif
     if (plotit eq 'yes') then set_plot,'x'
     if (plotit eq 'yes') then      window,0 
     if (plotit eq 'yes') then      window,1
     if (plotit eq 'yes') then      wset,1
     if (plotit eq 'yes' and Isothermal eq 0) then plot,T,press,/ylog,yrange=[pmax,pmin],xrange=[min(T)-200,max(T)+200],/ystyle
          if (plotit eq 'yes') then wset,0
     ;------------------------------------------------------------
     ;--initialze non-fit abundances at starting metallicity value using Chemical Equibrium
     if (Chemistry eq 1) then begin
       chemistry=0 ;temporarly set chemistry=0 to start from Equ value
       Result = FILE_TEST(dirs+'chem_eq_ARC_start'+string(ChildProcessCount,format='(I02)')+'.ncdf')  
       if (Result eq 0) then begin 
         ;run model
         fitabundances=''
         cloudtop=10^A(6) ;placeholder p(6) Cloud TOP
         art_haze_alpha=abs(A(7))
         fsed=A(9)
         sigma_g=A(10)
         ;arc_make_atmo_in,teff,MdH,logg,Rp,Cloud,haze,cloudtop,art_haze_alpha,kzz,fsed,sigma_g,fitabundances,fast=fast,contribution=contribution
         arc_make_atmo_in,A(0),exp(A(1)),A(2),A(3),exp(A(4)),exp(A(5)),cloudtop,art_haze_alpha,10^A(9),fsed,sigma_g,fitabundances,/fast
         ;
         spawn,'mv chem_eq_ARC'+string(ChildProcessCount,format='(I02)')+'.ncdf '+dirs+'chem_eq_ARC_start'+string(ChildProcessCount,format='(I02)')+'.ncdf'
       endif
       chemistry=1 ; set chemistry back to 1
       fid=ncdf_open(dirs+'chem_eq_ARC_start'+string(ChildProcessCount,format='(I02)')+'.ncdf')
       ncdf_varget,fid,'abundances',abundances
       ncdf_varget,fid,'pressure',pressure
       ncdf_varget,fid,'molname',molname
       ncdf_varget,fid,'element_abundances',element_abundances
       ;
       pressIndex=where(abs(pressure*1D-6-0.001) eq min(abs(pressure*1D-6-0.001)))
       if (Emission eq 'yes') then pressIndex=where(abs(pressure*1D-6-1) eq min(abs(pressure*1D-6-1))) ;use 1 bar for emission
       ;fit=where(Abund_fit_switch eq 0)
       ;AbunStart=abundances(milibarIndex(0),Abund_fit_index(fit)-1) passes only fitted abundances
       if (Chemistry eq 1) then  AbunStart=abundances(pressIndex(0),Abund_fit_index(0:numabunds)-1) ;passes all abundances for spectrally active molecules
       AbunStart=transpose(AbunStart)
       setthese=where(Abund_fit_switch eq 0)
       if (n_elements(AbunBest) gt 0) then AbunStart[setthese]=AbunBest 
       ;set out of bounds abundances to boundry
       for k=0,n_elements(AbunStart)-1 do begin & $
          if (AbunStart[k] lt exp(Abund_priorLow[k])) then AbunStart[k] = exp(Abund_priorLow[k])  & $
          if (AbunStart[k] gt exp(Abund_priorUP[k])) then AbunStart[k] = exp(Abund_priorUP[k])  & $
       endfor
        
       print,'---'
         print,'starting molecular abudnances:'
         print,'Name     ln(VMR)           log10(VMR)         VMR'
         for i=0,n_elements(Opacity_names)-1 do print,Opacity_names(i),((alog(AbunStart(i)))),((alog10(AbunStart(i)))),(((AbunStart(i))))
         print,'starting element abudnances:'
         print,'Name     ln(VMR)           log10(VMR)         VMR'
         for i=0,n_elements(Element_names)-1 do print,Element_names(i),element_abundances(i)
       print,'---'
       if (Chemistry eq 1) then A = [StartingParams,(alog(AbunStart))] ;add ln(abundances) to starting parameters
       if (Chemistry eq 1) then UpPriorsParams  = [UpPriorsParams,Abund_priorUP]   ;add Abudnance Priors
       if (Chemistry eq 1) then LowPriorsParams = [LowPriorsParams,Abund_priorLow] ;add Abudnance Priors
     endif
     ;--------------------
     if (Chemistry eq 3) then begin 
       A = [StartingParams,[EleBest]]    ;add ln(abundances) to starting parameters
       UpPriorsParams  = [UpPriorsParams,Ele_priorUP]   ;add Abudnance Priors
       LowPriorsParams = [LowPriorsParams,Ele_priorLow] ;add Abudnance Priors
     endif
     ;------------------------------------------------------------

     ;== print starting parameters and priors
     print,'---'
     print,'PARAMETER         Starting-value     Upper-prior     Lower-Prior    FITTED?'
     for i=0,n_elements(a)-1 do begin
       if fixit(i) eq 0 then willfit='Yes'
       if fixit(i) eq 1 then willfit='No'
      print,parname(i),A(i),UpPriorsParams(i),LowPriorsParams(i),'        ',willfit
       if (A(i) lt LowPriorsParams(i) and fixit(i) eq 0) then print,'   **WARNING!! ',parname(i),'Fit Param Starting Value or Lower Prior is BAD, Fix Needed !!**'
       if (A(i) gt UpPriorsParams(i) and fixit(i) eq 0)  then print,'   **WARNING!! ',parname(i),'Fit Param Starting Value or Upper Prior is BAD,  Fix Needed !!**'
       if (A(i) lt LowPriorsParams(i) and fixit(i) eq 1) then print,'   **WARNING!! ',parname(i),'Non-fit Param Starting Value or Lower Prior is BAD, Fix Needed**'
       if (A(i) gt UpPriorsParams(i) and fixit(i) eq 1)  then print,'   **WARNING!! ',parname(i),'Non-fit Param Starting Value or Upper Prior is BAD,  Fix Needed**'
     endfor
     print,'---'
 
     x=wsd & y=r & err=rerr & waverr=werr ;data
     save,filename=dirs+filenam+'arc_retrieval.setup.sav',/COMM,/VARIABLES,/COMPRESS
     allmodelpoints = dblarr(n_elements(wsd),maxsteps,nchains)
     allchemistry   = dblarr(nlevel+1,n_elements(Opacity_switch),maxsteps,nchains)
     if (isothermal eq 0) then allTP = dblarr(nlevel+1,maxsteps,nchains)
     ;=========================================================== 
     ;Leverberg-Marquart Least Squares Minimization (MPFIT) & Nested Sampling &  MCMC

     arc_fitforwardmodel,A,wsd,r,rerr,werr,wmodel,fmodel,params,pcerror,chi,parname,MPFITfilenameprefix,FIXIT,parnameHist,parnameHist2,DoMCMC,nthin,nchains,maxsteps,RestartMCMC,DoNEST,NUM,SAMPLE_NUM,TOL,NTHREADS,EXPAND,plotit,boostmcmc,DoMPFIT,RestartNest

     skip_arc_fitforwardmodel:
     if (Monitor eq 'yes') or (Monitor eq 'YES') or (Monitor eq 'Yes') then begin
       restore,dirs+filenam+'arc_retrieval.setup.sav'
       restore,dirs+'MPFIT.bestfit.sav'
       chi=BESTNORM;total((y-modelpoints)^2./err^2.)
     endif

     ;================================================================
     ;================================================================

   ;================================================================
   ;Post Processing Plots
   ;==================
   ;Use MCMC parameters if they are lower Chi^2 (comment out to use chi^2 min from MPFIT)
     if (DoMCMC eq 'yes') then begin
         arc_density,dirs
         restore,dirs+'MCMC.monitor.sav' ;allchains,allchainschi2,i,j,nchains,nthin,tofit
         ;Get min MCMC params 
           min_mcmc_chi=min(allchainschi2(0:i,*))
             for k=0,nchains-1 do if min(allchainschi2(0:i,k)) eq min(allchainschi2(0:i,*)) then begin print,k,min(allchainschi2(0:i,k),minindex),minindex & minchain=k & endif
           mcmc_params=allchains(*,minindex,minchain)
         if (min_mcmc_chi lt chi) then params(tofit)=mcmc_params  ;REPLACE MPFIT minimum parameters with MCMC determained values
     endif
     
     ; set parameters index 0 to 9
     teff          =   params(0)
     MdH           =exp(params(1))
     logg          =   params(2)
     Rp            =   params(3) ;Plantary Radius at 1 bar
     Cloud         =exp(params(4))
     haze          =exp(params(5))
     cloudtop      =10^params(6) ;placeholder p(6) Cloud TOP
     art_haze_alpha=abs(params(7))
     kzz           =10^params(8) ;KZZ
     fsed          =   params(9)
     sigma_g       =   params(10)
     Kir           =10^params(11) 
     gamma1        =10^params(12) 
     gamma2        =10^params(13) 
     beta          =   params(14)  
     alpha         =   params(15) 
     Tint          =   params(16)
     delta         =   params(17) ;P-T lapse rate
     Pdelta        =10^params(18) ;Pressure lapserate base

     ;teff=params(0) & MdH=exp(params(1))   & COratio=exp(params(2))   & Rp=(params(3))   & Cloud=exp(params(4))   & Haze=exp(params(5))+1d0   & cloudtop=10^params(6)   & delta=params(7)   & Pdelta=10^params(8)   & Kzz=10^params(9)
     ;teffin=teff    & MdHin=exp(params(1)) & COratioin=exp(params(2)) & Rpin=(params(3)) & Cloudin=exp(params(4)) & Hazein=exp(params(5))+1d0 & cloudtopin=10^params(6) & deltain=params(7) & Pdeltain=10^params(8) & Kzzin=10^params(9)
     if (Isothermal eq 0 and chemistry eq 1) then fitabundances=exp(params(19:n_elements(params)-1))
     if (Isothermal eq 1 and chemistry eq 1) then fitabundances=exp(params(19:n_elements(params)-1))
     if (chemistry eq 3 and isothermal eq 1) then fitabundances=exp(params(19:n_elements(params)-1)-params(1))
     if (chemistry eq 3 and isothermal eq 0) then fitabundances=exp(params(19:n_elements(params)-1)-params(1))
     if (Cloudtype eq 2) then Cloud=params(4) ;height of hard cloud deck, not opacity
     ; delete old files if necessary
      ChildProcessCount=0
      Result = FILE_TEST('ARC.transitspec'+string(ChildProcessCount,format='(I02)')+'.ncdf')  
      if (Result eq 1 ) then spawn,'rm ARC.transitspec'+string(ChildProcessCount,format='(I02)')+'.ncdf'
      Result = FILE_TEST('chem_eq_ARC'+string(ChildProcessCount,format='(I02)')+'.ncdf')  
      if (Result eq 1 ) then spawn,'rm chem_eq_ARC'+string(ChildProcessCount,format='(I02)')+'.ncdf'
 
    ;Make Low and High Resolution Versions of fitted ATMO Model
    ;==  Low RES==================
     ;If T-P profile was fit, make a best-fit T-P profile
     ChildProcessCount=98
     arc_ptmaker,nlevel,pmin,pmax  ;makes dummy files chem_dummy_arc.ncdf & input_TP_arc.ncdf with appropriate levels
     if (Isothermal eq 0) then begin
       ptfile='input_TP_arc.ncdf'
       arc_tp_delta,ptfile,logg,Kir,gamma1,gamma2,beta,alpha,T,Tint,Press,delta,Pdelta
      ;make copy P-T profile for use with each daughter process
      spawn,'cp input_TP_arc.ncdf pt_ARC'+string(ChildProcessCount,format='(I02)')+'.ncdf'
     endif
     ;RUN ATMO
     ;arc_make_atmo_in,teffin,MdHin,COratioin,Rpin,Cloudin,Hazein,fitabundances,cloudtop,Kzzin,/fast,/contribution ;run ATMO
     arc_make_atmo_in,teff,MdH,logg,Rp,Cloud,haze,cloudtop,art_haze_alpha,kzz,fsed,sigma_g,fitabundances,/fast,/contribution
     ;save files
       if (Emission eq 'no') then spawn,'cp ARC.transitspec'+string(ChildProcessCount,format='(I02)')+'.ncdf '+dirs+MPFITfilenameprefix+'.lowres.mcmc_bestfit.transitspec.ncdf'
       if (Emission eq 'yes') then spawn,'cp ARC.emissionspec'+string(ChildProcessCount,format='(I02)')+'.ncdf '+dirs+MPFITfilenameprefix+'.lowres.mcmc_bestfit.emissionspec.ncdf'
       spawn,'cp ARC'+string(ChildProcessCount,format='(I02)')+'.in '+dirs+MPFITfilenameprefix+'.lowres.mcmc_bestfit.ARC.in'
       spawn,'cp chem_eq_ARC98.ncdf '+dirs+MPFITfilenameprefix+'.lowres.mcmc_bestfit.chemistry.ncdf'
       spawn,'cp chem_dummy_arc.ncdf '+dirs+'chem_dummy_arc_lowres.ncdf'
       if (Isothermal eq 0) then spawn,'cp pt_ARC98.ncdf '+dirs+MPFITfilenameprefix+'.lowres.mcmc_bestfit.PT.ncdf'
       if (Isothermal eq 0) then spawn,'cp pt_ARC98.ncdf '+dirs
       if (Emission eq 'yes') then spawn,'cp ARC.emission.contribution'+string(ChildProcessCount,format='(I02)')+'.ncdf '+dirs+MPFITfilenameprefix+'.lowres.mcmc_bestfit.contribution_function.ncdf'
   ;==  HIGH RES==================
      ChildProcessCount=99 ;save high-res files with 99 label
      arc_ptmaker,nlevel,pmin,pmax  ;makes dummy files chem_dummy_arc.ncdf & input_TP_arc.ncdf with appropriate levels
      if (Isothermal eq 0) then begin
          ;set first 5 T-P profile parameters (delta & pdelta are set in first 9)
          ;=>!!!!!need to re-generate pt_File and chem_dummy at new number of levels
          ;Kir=10^params(10) & gamma1=10^params(11) & gamma2=10^params(12) & beta=params(13)  & alpha=params(14)
          ptfile='input_TP_arc.ncdf'
          arc_tp_delta,ptfile,logg,Kir,gamma1,gamma2,beta,alpha,T,Tint,Press,delta,Pdelta
          spawn,'cp input_TP_arc.ncdf pt_ARC'+string(ChildProcessCount,format='(I02)')+'.ncdf'
      endif
      ;arc_make_atmo_in,teffin,MdHin,COratioin,Rpin,Cloudin,Hazein,fitabundances,cloudtop,Kzzin;,/contribution ;run ATMO
      arc_make_atmo_in,teff,MdH,logg,Rp,Cloud,haze,cloudtop,art_haze_alpha,kzz,fsed,sigma_g,fitabundances;,/fast,/contribution
      ;save files
       if (Emission eq 'no') then spawn,'cp ARC.transitspec'+string(ChildProcessCount,format='(I02)')+'.ncdf '+dirs+MPFITfilenameprefix+'.highres.mcmc_bestfit.transitspec.ncdf'
       if (Emission eq 'yes') then spawn,'cp ARC.emissionspec'+string(ChildProcessCount,format='(I02)')+'.ncdf '+dirs+MPFITfilenameprefix+'.highres.mcmc_bestfit.emissionspec.ncdf'
       spawn,'cp ARC'+string(ChildProcessCount,format='(I02)')+'.in '+dirs+MPFITfilenameprefix+'.highres.mcmc_bestfit.ARC.in'
       spawn,'cp chem_dummy_arc.ncdf '+dirs+'chem_dummy_arc_highres.ncdf'
       spawn,'cp chem_eq_ARC99.ncdf '+dirs+MPFITfilenameprefix+'.highres.mcmc_bestfit.chemistry.ncdf'
       if (Isothermal eq 0) then spawn,'cp pt_ARC99.ncdf '+dirs+MPFITfilenameprefix+'.highres.mcmc_bestfit.PT.ncdf'
       if (Isothermal eq 0) then spawn,'cp pt_ARC99.ncdf '+dirs
       if (Emission eq 'yes') then spawn,'cp ARC.emission.contribution'+string(ChildProcessCount,format='(I02)')+'.ncdf '+dirs+MPFITfilenameprefix+'.highres.mcmc_bestfit.contribution_function.ncdf'
  ;================================================================
   ;==================
   ;Save State of ARC
     save,filename=dirs+'ARC.postprocessing.sav',/COMM,/VARIABLES,/COMPRESS ;Save state such that this file can be restored and individual plots below ran seperatly
     ;To Run programs below individually after a MPFIT/MCMC, restore following file first
     ;restore,dirs+'ARC.postprocessing.sav'
   ;==================

   ;==================
   ;Plot Abundances, i.e. Chemistry Spagetii plot
     spawn,'cp chem_eq_ARC'+string(ChildProcessCount,format='(I02)')+'.ncdf '+dirs+filenam+'.chem_eq_ARC.ncdf'
     if (chemistry eq 0 or chemistry eq 2 or chemistry eq 3) then arc_plotchem,dirs,filenam
   ;==================

   ;==================
   ;Make Corner Plot, T-P Plot, Spectra range plot, for MCMC run
     if (DoMCMC eq 'yes' and Isothermal eq 0) then arc_tp_plot_delta,dirs,r,Rs  ;if running seperatley first compile IDL>.r arc_tp.pro
     if (DoMCMC eq 'yes') then arc_spectra_plot,dirs,filenam,MPFITfilenameprefix,params,wsd,r,werr,rerr,modelpoints,ws,cloudtype,Emission,Rpl;,planet
     if (DoMCMC eq 'yes') then arc_density,dirs
   ;==================
   ;Make and Plot Individual Spectral Pieces For Transmission Spectra
    ChildProcessCount=99 ;save high-res files with 99 label
    if (Emission eq 'no') then arc_make_atmo_pieces,dirs,filenam,teff,MdH,logg,Rp,Cloud,haze,cloudtop,art_haze_alpha,kzz,fsed,sigma_g,fitabundances
    if (Emission eq 'no') then arc_plot_atmo_pieces,dirs,filenam,teff,MdH,COratio,Rp,Cloud,MPFITfilenameprefix,Rpl,Rs,logg

   ;================================================================
   ;Plot Highres Spectra and Data
   arc_plotspec,dirs,filenam,MPFITfilenameprefix,params,wsd,r,werr,rerr,modelpoints,ws,model,fixit,Emission
   ;==================

   ;==================
   ;Save IDL variables
   save,filename=dirs+filenam+'.sav',wsd,r,werr,rerr,x,y,err,waverr,StartingParams,FIXIT,params,teff,MdH,COratio,z0,Cloud,modelpoints,ws,model,haze,pp_rp,fitabundances
   ;================================================================


end
