pro arc_retrieval_input_Wasp121example
   ;====================================================================================================================
   ; INPUTS
   ; Four Sections to Enter:
   ; (1) Data and directory information
   ; (2) Switches controling ARC Data, fitting tools, basic ATMO setup
   ; (3) Starting Values and Priors of fitted Parameters
   ; (4) Planet and Star System Parameters
   ;====================================================================================================================
   Compile_Opt DEFINT32

   ;====================================================================================================================
   ; **Begining of Inputs to EDIT**
   ;====================================================================================================================

   ;======================================================================================================
   ;== (1) Data & Filename Information ===
    dirs='/data/sing/ARC_Results/ARC_Wasp121_Trans_Fe_25May2018/'             ; recomend in /data/.. on um-serve & post-zen
    filenam='Wasp121.Trans.Fe.25May2018'                      ; Name to give files
    datafile='w121_trspec.g430l_n30.g750l_n30.ldatlas_free.G430Lshifted.txt'  ; ALL DATA, look in  atmo/arc/data   
   ;======================================================================================================

   ;======================================================================================================
   ;== (2) ARC Switches === Controls Fitting Routines & Basic ATMO Setup 
   ;================================================================
    Planet    ='Wasp121'; 'Wasp39','Hat18','Hat26','Wasp121','Wasp39','Wasp43'
    Emission  = 'no'  ; 'yes','no' ; Fit Emission Spectra? yes=Emission ; no=Transmission

    ;Nested Sampling
    DoNEST     = 'no' ; Perform Nested Sampling fit?  (Skips Nested Sampling if set to 'no')
    NUM        = 30   ; OPTIONAL: THIS IS THE NUMBER OF LIVE POINTS 
    SAMPLE_NUM = 30   ; OPTIONAL: THIS IS THE NUMBER OF SAMPLES THAT IT TAKES BEFORE RECLUSTERING
    EXPAND     = 1.06d; OPTIONAL: THIS EXPANSION FACTOR EXPANDS THE SIZE OF THE ELLIPSOIDS
    TOL        = 500d0; OPTIONAL: THIS IS THE TOLERANCE/Accuracy set on measuring the Evidence LnZ
    NTHREADS   = 10   ; Number of parallel processors to sample with (should be < 1x number of processors)

    ;MCMC
    DoMCMC     = 'yes' ; Perform MCMC fit?  (Skips MCMC if set to 'no')
    RestartMCMC= 'no' ;'yes' or 'no' Restart MCMC based on last MCMC.monitor.sav?
    nthin      =   1  ; MCMC thinning parameter                                       (should work but still untested)
    nchains    =  12  ; Number of MCMC chains and Number of processor to run on  (should be < 1x number of processors)
    maxsteps   = 30000 ; Number of MCMC steps for each chain  (MCMC runs until Gel.-Rub.<1.01 & 1000 independant draws)

    ;Data
    wsdstart   = 0.3  ; Perform Blue Cut in Wavelength of data if desired

    ;Atmosphere
    Cloudtype  = 1    ; Cloud is  0 = none                ; 1 = Uniform                 ; 2 = Hard Base   ; 3 = Wavelength Dependant
    Hazetype   = 0    ; Haze  is  0 = none                ; 1 = Rayleigh               
    Chemistry  = 1    ; Chemistry 0 = Chemical Equibrium  ; 1 = Manual 'free' Chemistry ; 2 =  'Non-Equilbrium' Chemistry
    Isothermal = 1    ; Temperature is 0 = fit            ; 1 = Isothermal

    ;Atmosphere Grid 
    nlevel    = 50    ; levels in the grid (grid is uniform in log-pressure)
    pmin      = 1E-8  ; Minimum pressure bar
    pmax      = 5E2   ; Maximum pressure bar

    ;
    Tint      = 100  ; Internal Temperature (used for T-P profile fits)

    ;Screen Plotting Option
    plotit    = 'yes' ; 'yes' or 'no' Plot fit to screen during MPFIT
                     if (plotit eq 'no') then print,"RECOMEND starting XTERM IDL within 'screen' first"
                     if (plotit eq 'no') then print,''
   ;======================================================================================================

   ;======================================================================================================
   ;== (3) Select Parameters to Fit and Their Starting Value ==
   ;======================================================================================================
   ; SET FIXIT  0=free 1=fix 
   ;          index   0                 1                     2               3                     4                           5            6               7            8           9
    parnameHist=   ['Temperature (K)','Abundance ln(Z/Zs)','ln(C/O) Ratio','Rp(1 bar)'          ,'cloud opacity ln(sig/sig0)','haze opacity','cloud top','delta',        'Pdelta     ','log10(Kzz cm^2/s)']
    FIXIT  =       [ 0               , 1                   , 1            , 0                   , 0                          , 1            ,     1     ,       1      ,     1        ,        1         ]  
    StartingParams=[1815.2959d0      , alog(25d0)          , alog(0.56d0) ,1.7395418D0          , 0.25d0                     , 0.0d0        ,     0.0d0 ,       0.0d0  ,     4.0d0   ,        9.5d0    ]   
   ;Bounding Box Priors
    LowPriorsParams=[ 400d0          , -5                  , alog(1D-3)   ,StartingParams[3]*0.9, -10                        ,  -10         ,     -1     ,      -1.0     , alog10(pmin),         5         ]   
    UpPriorsParams =[3000d0          ,  4.25               , alog(1D5)    ,StartingParams[3]*1.1,  10                        ,  10          ,     2.5   ,        1.0     , alog10(pmax),        11.0         ] 
   ;==============================
   ; Set Opacities to include for transmission/emission spectra calculation 
    Opacity_names   =['H2','He','H2O','CO2','CO','CH4','NH3','Na','K','Li','Rb','Cs','TiO','VO','FeH','PH3','H2S','HCN','C2H2','SO2','Fe','H_'] ; Switch:      1=include  0=Omit
    Opacity_switch  =[   1,   1,    1,    1,   1,    1,    1,   1,  1,   0,   0,   0,    1,   1 ,    1,   0,     0,    0,    0,    0,   1,   0] ;Hat26
   ; SET FIXIT    0=free 1=fix 
    Abund_fit_switch=[   1,   1,    0,    1,   1,    1,    1,   1,  1,   1,   1,   1,    0,    0,    0,   1,     1,    1,     1,   1,   0,   1] ;
   ; Set Abundances to include to vary (reqiures Chemistry=1 & Opacity_switch=1)

   ;==Option Start Abundances from Specific Values 
   ;==order is same as Abund_fit_switch, number of elements set to num of freely fit abundances
   ;==Uncomment line below to activate
       ;AbunBest=[exp(-2.9957323), exp(-6.5921373)];,exp(-20.636102),exp(-22.821792),exp(-13.396199),exp(-16.240493)] ;Volume Mixing Ratio
       ;AbunBest=[exp(-4.2026960d0),exp(-7.9149090d0)];,exp(-20.636102),exp(-22.821792),exp(-13.396199),exp(-16.240493)] ;Volume Mixing Ratio                                                                                                                    
   ;Abundance Bounding Box Priors
    Abund_priorLow = [(dblarr(22)+1)*(alog(1D-15))]               ;minimum abundance VMR of 1x10^-12 for non-H/He 
    Abund_priorUP  = [[-0.05,-0.05],(dblarr(20)+1)*(alog(0.05))]  ;maximum abundance VMR of 0.05 for non-H/He 

   ;==============================
   ; Fit T-P profile parameters (0=free) 
   ; Set Isothermal = 0 to activate
   ;             10              11                 12                 13           14
   TP_names =  ['Kir'         ,'gamma1'          ,'gamma2'          ,'beta'      ,'alpha' ]
   ;TP_params=  [alog10(7.1D-2), alog10(0.03317D0), alog10(0.03317D0), 0.946d0    , 0.5d0 ]  ;Line et al. HD209
   ;TP_params=  [-1.5130983,      1.49999          ,1.4999    ,0.60413045  ,0.0d0]  ;W121newmin
   TP_params=  [-1.6565581,      1.1306451        ,1.4999    ,0.72919111  ,0.0d0]  ;W121newmin                                                                                                                                                              
   TP_switch =  [        0,            0,             1,         0,            1]  
   ;Bounding Box Priors
   TP_priorLow=[  -5.0        , -4.0             ,  -4.0            , 0          ,  0    ] 
   TP_priorUP= [  -0.5        ,  1.5             ,   1.5            , 2          ,  1    ]
   ;======================================================================================================

   ;================================================================
   ;== (4) Planet and Star SYETEM PARAMETERS ==
   ;================================================================
   Rp    = 1.6938182D0               ; Radius of Planet (minimum of transit spectra)
   Teff  = 1700D0               ; Original Equilbrium Temperature
   g     = 8.43D2                ; cm/s^2 gravity
   rstar = 1.458                ; Rsun  - Will be used in emission spectra calculation with (Fp/Fs)*(Rp/Rs)^2
   rorbit= 0.02544               ;AU 
   stellarfile='../../arc/stellar_spectra/WASP121A.Phoenix.ncdf'  ;location of stellar file
   Tstar = 6460D0                                               ;K 
   ;================================================================

   ;====================================================================================================================
   ; **End of Inputs to EDIT**
   ;====================================================================================================================



   ;================================================================
   ;unit convertion
   Rjup  = 71492D5         ; cm
   Rpl   = Rp*Rjup         ; Rjup (cm)      
   logg  = alog10(g)       ; logg
   Rs    = rstar*6.96D10   ;cm
   ;================================================================
   ;==============================
   ; Save Inputs
   Result = FILE_TEST(dirs, /DIRECTORY) & if (Result eq 0) then spawn,'mkdir '+dirs
   save,filename=dirs+'arc_inputs.sav',AbunBest,ABUND_FIT_SWITCH,CHEMISTRY,CLOUDTYPE,DATAFILE,DIRS,DOMCMC,EMISSION,FILENAM,FIXIT,HAZETYPE,$
                                       ISOTHERMAL,MAXSTEPS,NCHAINS,NLEVEL,NTHIN,OPACITY_NAMES,OPACITY_SWITCH,PARNAMEHIST,PLANET,PMAX,PMIN,RESTARTMCMC,$ 
                                       STARTINGPARAMS,TINT,TP_NAMES,TP_PARAMS,TP_SWITCH,WSDSTART,LowPriorsParams,UpPriorsParams,Abund_priorLow, Abund_priorUP,TP_priorLow,TP_priorUP,$
                                       NUM,SAMPLE_NUM,TOL,NTHREADS,EXPAND,PLOTIT,DONEST,Rp,Rjup,Rpl,Teff,logg,g,Rstar,Rs,rorbit,stellarfile,Tstar
   ;return,dirs
   ;==============================

   ;====================================================================================================================
   ; RUN RETRIEVAL
   ;====================================================================================================================
   !path = '../../arc/:'+ !path      ; adds path to so IDL finds atmo/arc/arc_retrieval.pro                                        
   input=dirs+'arc_inputs.sav'
   arc_retrieval,input
   ;====================================================================================================================

end

