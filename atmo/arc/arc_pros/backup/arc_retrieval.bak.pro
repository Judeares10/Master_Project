pro arc_retrieval
; ATMO Retrieval Code
; PURPOSE:
;  Perform model fit retrieval techniques to ATMO model using DEMC
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
;  Fit then does a differential evolution chain monte carlo fit
;  (Eastman et al. 2014)
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
;   1) INPUT PARAMETES into arc_retrieval.pro (below)
;      a) Set Saving directory, Saving Filename, and DATA file below
;      b) Set Starting Parameters Below
;      c) Set Switches below
;   4) Save arc_retrieval.pro with Desired Inputs to be run
;   5) Compile and Run  
;      IDL> .r arc_retrival
;      IDL> arc_retrieval
;
;   Parameters to set:
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
;      step will always be rejected.  Hard coded in arc_model_chi_mcmc.pro   
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
;      Gaussian priors can be added, but are not yet implimented
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
; Version 0.95 To Do 1) add array of priors 2) add Kzz parameter & switch for NEQ
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
  Common SystemParameters,rplanet,logg,rstar,rorbit,stellarfile,Tstar,Cloudtype,Hazetype,Chemistry,Isothermal,Opacity_switch,Abund_fit_switch,Abund_fit_index,Abund_priorLow,Abund_priorUP
  Common cpcount,ChildProcessCount,dirs,allmodel,allmodelpoints,Emission,nlevel,pmin,pmax
  common passmodel,model,modelpoints,LowPriorsParams,UpPriorsParams

   set_plot,'x'
  ;================================================================
   ChildProcessCount=0
   ;====================================================================================================================
   ; INPUTS
   ;=================

   ;== Data & Filename Information ===
    dirs='/data/sing/ARC_Results/ARCv1Kzztest/'
    filenam='Hat26.kzztest'
    datafile='H26_transmission_ecc_10_8d_mean2_113-167_newwave.STISMarg.G096cust.txt'

   ;== ARC Switches ===
    Planet='Hat26'    ;Wasp39'   ; 'Hat18','Hat26','Wasp121','Wasp39'
    Emission   = 'no' ; 'yes','no' Fit Emission Spectra

    ;Nested Sampling
    DoNEST='no'      ; Perform Nested Sampling fit?  (Skips Nested Sampling if set to 'no')
    NUM        =30   ; OPTIONAL: THIS IS THE NUMBER OF LIVE POINTS 
    SAMPLE_NUM =30    ; OPTIONAL: THIS IS THE NUMBER OF SAMPLES THAT IT TAKES BEFORE RECLUSTERING
    EXPAND     =1.06d ; OPTIONAL: THIS EXPANSION FACTOR EXPANDS THE SIZE OF THE ELLIPSOIDS
    TOL        =500d0   ; OPTIONAL: THIS IS THE TOLERANCE/Accuracy set on measuring the Evidence LnZ
    NTHREADS   = 10    ; Number of parallel processors to sample with (should be < 1x number of processors)

    ;MCMC
    DoMCMC='no'       ; Perform MCMC fit?  (Skips MCMC if set to 'no')
    RestartMCMC='no'  ;'yes' or 'no' Restart MCMC based on last MCMC.monitor.sav?
    nthin=1           ; MCMC thinning parameter                                       (should work but still untested)
    nchains=12;       ; Number of MCMC chains and Number of processor to run on  (should be < 1x number of processors)
    maxsteps=50;30000    ; Number of MCMC steps for each chain  (MCMC runs until Gel.-Rub.<1.01 & 1000 independant draws)

    ;Data
    wsdstart=0.3      ; Perform Blue Cut in Wavelength if desired

    ;Atmosphere
    Cloudtype  = 1    ; Cloud is  0 = none                ; 1 = Uniform                 ; 2 = Hard Base   ; 3 = Wavelength Dependant
    Hazetype   = 0    ; Haze  is  0 = none                ; 1 = Rayleigh               
    Chemistry  = 2    ; Chemistry 0 = Chemical Equibrium  ; 1 = Manual 'free' Chemistry ; 2 =  'Non-Equilbrium' Chemistry
    Isothermal = 1    ; Temperature is 0 = fit            ; 1 = Isothermal

    ;Atmosphere Grid 
    nlevel=20;32         ; levels in the grid (grid is uniform in log-pressure)
    pmin=1E-5         ; Minimum pressure bar
    pmax=5E3          ; Maximum pressure bar

    ;
    Tint      = 100  ; Internal Temperature (used for T-P profile fits)

   ;=================
   ; Select Parameters to Fit and Their Starting Value
   ; SET FIXIT  0=free 1=fix 
   ;          index   0                 1                     2         3                    4                           5            6               7            8           9
    parnameHist=   ['Temperature (K)','Abundance ln(Z/Zs)','ln(C/O) Ratio','Rp(1 bar)' ,'cloud opacity ln(sig/sig0)','haze opacity','cloud top','cloud bottom','cloud wave','log10(Kzz cm^2 s)']
    FIXIT  =       [ 0               , 0                   , 0            , 0          , 0                          , 1            ,     1     ,       1      ,     1        ,        0         ]  
    StartingParams=[896.57271d0      , 0.29853348d0        , alog(0.56d0) ,0.57273258d0, -0.52471765d0              , 0.0d0        ,     0.0d0 ,       0.0d0  ,     4.0d0   ,        9.5d0    ]   

   ;Bounding Box Priors
    LowPriorsParams=[ 550d0          , -5                  , alog(1D-2)   ,0.561*0.9  , -2.5                       ,  -10         ,     0     ,       0     ,     0       ,         5         ]   
    UpPriorsParams =[1800d0          ,  4.25               , alog(1D5)    ,0.561*1.1  ,  6                         ,  10          ,     2.5   ,     2.5     ,     6       ,        11.0         ]   
 
   ;==============================
   ; Set Opacities to include for transmission spectra calculation 
    Opacity_names   =['H2','He','H2O','CO2','CO','CH4','NH3','Na','K','Li','Rb','Cs','TiO','VO','FeH','PH3','H2S','HCN','C2H2']
   ; Switch:      1=include  0=Omit
    Opacity_switch  =[   1,   1,    1,    1,   1,    1,    1,   0,  0,   0,   0,   0,    0,   0,    0,   0,     0,    0,     0] ;Hat26
   ; SET FIXIT    0=free 1=fix 
    Abund_fit_switch=[   1,   1,    0,    0,   0,    0,    1,   1,  1,   1,   1,   1,    1,   1,    1,   1,     1,    1,     1] ;
   ; Set Abundances to include to vary (reqiures Chemistry=1 & Opacity_switch=1)

   ;Option Start Abundances from Specific Values, order same as Abund_fit_switch,uncomment line below
       ;AbunBest=[exp(-6.4579836), exp(-4.6593652),exp(-20.636102),exp(-22.821792),exp(-13.396199),exp(-16.240493)] 

   ;Bounding Box Priors
    Abund_priorLow = [(dblarr(19)+1)*(alog(1D-12))]               ;minimum abundance VMR of 1x10^-12 for non-H/He 
    Abund_priorUP  = [[-0.05,-0.05],(dblarr(17)+1)*(alog(0.05))]  ;maximum abundance VMR of 0.05 for non-H/He 

   ;==============================
   ; Fit T-P profile parameters (0=free) 
   ; Set Isothermal = 0 to activate
   ;             10              11                 12                 13           14
   TP_names =  ['Kir'         ,'gamma1'          ,'gamma2'          ,'beta'      ,'alpha' ]
   ;TP_params=  [alog10(7.1D-2), alog10(0.03317D0), alog10(0.03317D0), 0.946d0    , 0.5d0 ]  ;Line et al. 
   TP_params=  [alog10(19.67) , alog10(0.001D0)  , alog10(0.03317D0), 0.5473d0   , 0.0d0 ]
   TP_params=  [-2.5014450    , -0.86841635      , alog10(0.03317D0), 1.0184612d0, 0.0d0 ]  ; HAT-P-26b M6
   TP_switch=  [          0   ,          0       ,        1         ,    0       ,   1   ]  
   ;Bounding Box Priors
   TP_priorLow=[  -5.0        , -4.0             ,  -4.0            , 0          ,  0    ] 
   TP_priorUP= [  -0.5        ,  1.5             ,   1.5            , 2          ,  1    ]

   ;==============================
   ; Save Inputs
   Result = FILE_TEST(dirs, /DIRECTORY) & if (Result eq 0) then spawn,'mkdir '+dirs
   save,filename=dirs+'arc_inputs.sav',ABUND_FIT_SWITCH,CHEMISTRY,CHILDPROCESSCOUNT,CLOUDTYPE,DATAFILE,DIRS,DOMCMC,EMISSION,FILENAM,FIXIT,HAZETYPE,$
                                       ISOTHERMAL,MAXSTEPS,NCHAINS,NLEVEL,NTHIN,OPACITY_NAMES,OPACITY_SWITCH,PARNAMEHIST,PLANET,PMAX,PMIN,RESTARTMCMC,$ 
                                       STARTINGPARAMS,TINT,TP_NAMES,TP_PARAMS,TP_SWITCH,WSDSTART,LowPriorsParams,UpPriorsParams,Abund_priorLow, Abund_priorUP,TP_priorLow,TP_priorUP,$
                                       NUM,SAMPLE_NUM,TOL,NTHREADS,EXPAND
   ;====================================================================================================================
   ; End of Inputs
   ;====================================================================================================================



   ;================================================================
   ;Should't need to change things below here
   ;====================================================================================================================
   !path = './arc_pros/:'+ !path                                                 ;remember to change arc_multinest.pro threads setup
   ;!PATH = Expand_Path('+/Users/sing/sing/IDL/Routines/coyote/') + ':' + !PATH

   ;======

   parnameHist2=  ['Teff'         ,'lnMdH'        ,'COratio'      ,'Rp'            ,'lnCloud'      ,'ln(Haze-1)'   ,'Cloud Top'     ,'Cloud Bottom'  ,'Cloud Wave'    ,'log10(Kzz)'       ]
   parname=       ['Teff(K)      ','ln(MdH)      ','ln(C/O)      ','Rp(RJup 1mbar)','ln(Cloud)    ','ln(haze)     ','Rp(RJup)       ','Rp(RJup)     ','Cloud Wave    ','log10(Kzz cm^2s)' ]
   if (Cloudtype eq 0) then FIXIT(4)=1  ;Don't Vary Clouds when there are none there. StartingParam[4] will be ignored
   if (Hazetype  eq 0) then FIXIT(5)=1  ;Don't Vary Haze when there is none there. StartingParam[5] will be ignored
   if (Chemistry eq 1) then FIXIT(1)=1  ;Don't Vary Metalcity when Fitting for individual Abundances
   if (Chemistry eq 1) then FIXIT(2)=1  ;Don't Vary C/O when Fitting for individual Abundances
   if (Chemistry eq 2) then FIXIT(9)=0  ;Vary Kzz when fitting non-equibrium chemsitry
   if (Chemistry ne 2) then FIXIT(9)=1  ;Don't Vary Kzz unless fitting non-equilibrium chemistry
   if (Cloudtype eq 2) then parnameHist=   ['Temperature (K)','Abundance ln(Z/Zs)','C/O Ratio','Rpl(1bar)','cloud deck height','haze opacity','Cloud Top'     ,'Cloud Bottom'  ,'Cloud Wave'    ,'log10(Kzz)'  ]
   if (Cloudtype eq 2) then parnameHist2=  ['Teff'           ,'lnMdH'             ,'COratio'  ,'Rp'              ,'Cloud'            ,'ln(Haze-1)','CloudTop(Rp)','CloudBot(Rp)','CloudWave','log10(Kzz)']
   if (Isothermal eq 0) then FIXIT=[FIXIT,TP_switch]
   if (Isothermal eq 0) then FIXIT(0)=1 ;Don't Vary internal temperature if Fitting for T-P profile
   if (Chemistry eq 1) then FIXIT=[FIXIT,Abund_fit_switch]
   ; Isothermal       10   11   12    13    14
   ; T-P fit          15   16   17    18    19   20    21    22  23  24   25   26   27   28     29   30   31     32     33
   Opacity_names   =['H2','He','H2O','CO2','CO','CH4','NH3','Na','K','Li','Rb','Cs','TiO','VO','FeH','PH3','H2S','HCN','C2H2'] ;Spectrally active molecules
   Abund_fit_index =[   1,  62,    9,   15,   4,   12,   68, 108,114, 156, 164, 160,  128, 130,  147, 170,   137,   76,    17] ;Atmo Index (IDL index is -1 of this value)
  ;Abund_fit_index =[   1,  62,    9,   15,   4,   12,   68, 108,113, 147, 155, 151,  126, 128,  140, 161,   133,   76,    17] ;Atmo Index (IDL index is -1 of this value)
   ;
   if (Isothermal eq 0) then parname=[parname,TP_names]
   if (Isothermal eq 0) then parnameHist=[parnameHist,TP_names]
   if (Isothermal eq 0) then parnameHist2=[parnameHist2,TP_names]
   if (Chemistry eq 1) then parnameHist=[parnameHist,Opacity_names(0:18)]
   if (Chemistry eq 1) then parnameHist2=[parnameHist2,Opacity_names(0:18)]
   if (Chemistry eq 1) then parname=[parname,Opacity_names(0:18)]
   ;
   if (Isothermal eq 0) then StartingParams(0) = Tint ;Set temperature to internal temperature if fitting for T-P profile
   if (Isothermal eq 0) then parname(0)='Tint'
   if (Isothermal eq 0) then parnameHist2(0)='Tint'
   if (Isothermal eq 0) then StartingParams =  [StartingParams,TP_params]    ;Set temperature to internal temperature if fitting for T-P profile
   if (Isothermal eq 0) then UpPriorsParams  = [UpPriorsParams,TP_priorUP]   ;Add temperature profile priors
   if (Isothermal eq 0) then LowPriorsParams = [LowPriorsParams,TP_priorLow] ;Add temperature profile priors
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
     arc_load_system_parameters,planet,Rp,Rjup,Rpl,Teff,logg,g,Rstar,Rs,rorbit,stellarfile,Tstar
     rplanet=Rpl ;Radius of the planet in cm.  Used for Emission Flux calculation, fnu/fstar*(rplanet/rs)^2.
     ;initialize arrays to save models during MCMC

     allmodel = dblarr(350,maxsteps,nchains)  ;350 wavelentgh points at lowres resolution 
     if (Emission eq 'yes')  then allmodel = dblarr(150,maxsteps,nchains)  ;350 wavelentgh points at lowres resolution 
     allmodelpoints = dblarr(n_elements(wsd),maxsteps,nchains) 
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
         Kir=10^A(10) & gamma1=10^A(11) & gamma2=10^A(12) & beta=A(13)  & alpha=A(14)
         arc_tp,ptfile,Kir,gamma1,gamma2,beta,alpha,T,Tint,Press
         spawn,'cp pt_ARC00.ncdf '+dirs+'pt_ARC_start.ncdf'
      endif
     set_plot,'x'
     window,0 
     window,1
     wset,1
     if (Isothermal eq 0) then plot,T,press,/ylog,yrange=[100,1E-6],xrange=[500,2500],/ystyle
     wset,0
     ;initialze abundances at starting value using Chemical Equibrium
     if (Chemistry eq 1) then begin
       chemistry=0 ;temporarly set chemistry=0 to start from Equ value
       Result = FILE_TEST(dirs+'chem_eq_ARC_start'+string(ChildProcessCount,format='(I02)')+'.ncdf')  
       if (Result eq 0) then begin 
         ;run model
         arc_make_atmo_in,A(0),exp(A(1)),exp(A(2)),A(3),exp(A(4)),exp(A(5))+1d0,10^A(9),/fast
         ;
         spawn,'mv chem_eq_ARC'+string(ChildProcessCount,format='(I02)')+'.ncdf '+dirs+'chem_eq_ARC_start'+string(ChildProcessCount,format='(I02)')+'.ncdf'
       endif
       chemistry=1 ; set chemistry back to 1
       fid=ncdf_open(dirs+'chem_eq_ARC_start'+string(ChildProcessCount,format='(I02)')+'.ncdf')
       ncdf_varget,fid,'abundances',abundances
       ncdf_varget,fid,'pressure',pressure
       ncdf_varget,fid,'molname',molname
       milibarIndex=where(abs(pressure*1D-6-0.001) eq min(abs(pressure*1D-6-0.001)))
       ;fit=where(Abund_fit_switch eq 0)
       ;AbunStart=abundances(milibarIndex(0),Abund_fit_index(fit)-1) passes only fitted abundances
       if (Chemistry eq 1) then AbunStart=abundances(milibarIndex(0),Abund_fit_index(0:18)-1) ;passes all abundances for spectrally active molecules
       setthese=where(Abund_fit_switch eq 0)
       if (n_elements(AbunBest) gt 0) then AbunStart(setthese)=AbunBest 
       if (Chemistry eq 1) then A = [StartingParams,transpose(alog(AbunStart))] ;add ln(abundances) to starting parameters
       if (Chemistry eq 1) then UpPriorsParams  = [UpPriorsParams,TP_priorUP]   ;add Abudnance Priors
       if (Chemistry eq 1) then LowPriorsParams = [LowPriorsParams,TP_priorLow] ;add Abudnance Priors
      endif

     x=wsd & y=r & err=rerr & waverr=werr ;data
     save,filename=dirs+filenam+'arc_retrieval.setup.sav',/COMM,/VARIABLES,/COMPRESS
     ;=========================================================== 
     ;Leverberg-Marquart Least Squares Minimization (MPFIT) & Nested Sampling &  MCMC
     arc_fitforwardmodel,A,wsd,r,rerr,werr,wmodel,fmodel,params,pcerror,chi,parname,MPFITfilenameprefix,FIXIT,parnameHist,parnameHist2,DoMCMC,nthin,nchains,maxsteps,RestartMCMC,DoNEST,NUM,SAMPLE_NUM,TOL,NTHREADS,EXPAND
 
     ;================================================================
     ;================================================================

   ;================================================================
   ;Post Processing Plots
   ;==================
   ;Use MCMC parameters (comment out to use chi^2 min from MPFIT)
     ;if (DoMCMC eq 'yes') then restore,dirs+'MCMC.monitor.FitPars.sav';MCMCpars,MCMCPCERRORlow,MCMCPCERRORhigh,MCMCPCERROR
     ;if (DoMCMC eq 'yes') then params=MCMCpars
     ; set parameters
     teff=params(0) & MdH=exp(params(1)) & COratio=exp(params(2)) & Rp=(params(3)) & Cloud=exp(params(4)) & Haze=exp(params(5))+1d0 & Kzz=params(9)
     teffin=teff    & MdHin=exp(params(1)) & COratioin=exp(params(2)) & Rpin=(params(3)) & Cloudin=exp(params(4)) & Hazein=exp(params(5))+1d0
     Kzzin=10^params(9)
     if (Isothermal eq 0 and chemistry eq 1) then fitabundances=exp(params(15:n_elements(params)-1))
     if (Isothermal eq 1 and chemistry eq 1) then fitabundances=exp(params(10:n_elements(params)-1))
     if (Cloudtype eq 2) then Cloud=params(4) ;height of hard cloud deck, not opacity
     ; delete old files if necessary
      ChildProcessCount=0
      Result = FILE_TEST('ARC.transitspec'+string(ChildProcessCount,format='(I02)')+'.ncdf')  
      if (Result eq 1 ) then spawn,'rm ARC.transitspec'+string(ChildProcessCount,format='(I02)')+'.ncdf'
      Result = FILE_TEST('chem_eq_ARC'+string(ChildProcessCount,format='(I02)')+'.ncdf')  
      if (Result eq 1 ) then spawn,'chem_eq_ARC'+string(ChildProcessCount,format='(I02)')+'.ncdf'
 
    ;Make High Resolution Version of fitted ATMO Model
     ;If T-P profile was fit, make a best-fit T-P profile
     if (Isothermal eq 0) then begin
       Kir=10^params(10) & gamma1=10^params(11) & gamma2=10^params(12) & beta=params(13)  & alpha=params(14)
       arc_tp,ptfile,Kir,gamma1,gamma2,beta,alpha,T,Tint,Press
       spawn,'cp pt_ARC00.ncdf '+dirs+'pt_ARC_MPFIT_bestfit.ncdf'
     endif
     arc_make_atmo_in,teffin,MdHin,COratioin,Rpin,Cloudin,Hazein,fitabundances,Kzzin
     if (Emission eq 'no') then spawn,'cp ARC.transitspec'+string(ChildProcessCount,format='(I02)')+'.ncdf '+dirs+MPFITfilenameprefix+'.highres.transitspec.ncdf'
     if (Emission eq 'yes') then spawn,'cp ARC.emissionspec'+string(ChildProcessCount,format='(I02)')+'.ncdf '+dirs+MPFITfilenameprefix+'.highres.emissionspec.ncdf'
     spawn,'cp ARC'+string(ChildProcessCount,format='(I02)')+'.in '+dirs+MPFITfilenameprefix+'.highres.ARC.in'
   ;================================================================
   ;==================
   ;Save State of ARC
     save,filename=dirs+'ARC.postprocessing.sav',/COMM,/VARIABLES,/COMPRESS ;Save state such that this file can be restored and individual plots below ran seperatly
     ;To Run programs below individually after a MPFIT/MCMC, restore following file first
     ;restore,dirs+'ARC.postprocessing.sav'
   ;==================
   ;================================================================
   ;Plot Spectra and Data
     arc_plotspec,dirs,filenam,MPFITfilenameprefix,params,wsd,r,werr,rerr,modelpoints,ws,model,fixit,Emission
   ;==================

   ;==================
   ;Plot Abundances, i.e. Chemistry Spagetii plot
     spawn,'cp chem_eq_ARC'+string(ChildProcessCount,format='(I02)')+'.ncdf '+dirs+filenam+'.chem_eq_ARC.ncdf'
     if (chemistry eq 0) then arc_plotchem,dirs,filenam
   ;==================

   ;==================
   ;Make and Plot Individual Spectral Pieces For Transmission Spectra
     if (Emission eq 'no') then arc_make_atmo_pieces,dirs,filenam,teffin,MdHin,COratioin,Rpin,Cloudin,Hazein,fitabundances,kzzin
     if (Emission eq 'no') then arc_plot_atmo_pieces,dirs,filenam,planet,teff,MdH,COratio,Rp,Cloud,MPFITfilenameprefix
   ;==================
   ;Make Corner Plot, T-P Plot, Spectra range plot, for MCMC run
     if (DoMCMC eq 'yes') then arc_density,dirs
     if (DoMCMC eq 'yes' and Isothermal eq 0) then arc_tp_plot,planet,dirs,r
     if (DoMCMC eq 'yes') then arc_spectra_plot,dirs,filenam,MPFITfilenameprefix,params,wsd,r,werr,rerr,modelpoints,ws,cloudtype,Emission
   ;==================
   ;Save IDL variables
   save,filename=dirs+filenam+'.sav',wsd,r,werr,rerr,x,y,err,waverr,StartingParams,FIXIT,params,teff,MdH,COratio,z0,Cloud,modelpoints,ws,model,haze,pp_rp,fitabundances
   ;================================================================


end
