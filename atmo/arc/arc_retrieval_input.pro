pro arc_retrieval_input
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
    dirs='/mnt/data/dsing3/ARC_Results/Wasp127_feeq_hazealpha2_TPfree_8July2020/'             ; recomend in /data/.. on um-serve & post-zen
    filenam='Wasp127.haze_alpha_TPfree'                      ; Name to give files
    datafile='wasp127_final_spec_common.txt'  ; ALL DATA, look in  atmo/arc/data   
   ;======================================================================================================

   ;======================================================================================================
   ;== (2) ARC Switches === Controls Fitting Routines & Basic ATMO Setup 
   ;================================================================
    Planet    ='Wasp127'; 'Wasp39','Hat18','Hat26','Wasp121','Wasp39','Wasp43'
    Emission  = 'no'  ; 'yes','no','DI' ; Fit Emission Spectra? yes=Emission ; no=Transmission; DI=Direct Imaging

    Monitor   = 'no'  ; 'yes',no'  ; Make preliminary post-processing plots (e.g. posterior) while seperate on-going MCMC is still running? 

    ;Levenberg-Marquardt Least Squares Fit
    DoMPFIT   = 'no'  ; 'yes','no'  ; Perform Levenberg-Marquardt Least Squares Fit? Recomended if doing MCMC, as it starts around best-fit value.

    ;Nested Sampling
    DoNEST     = 'yes'; Perform Nested Sampling fit?  (Skips Nested Sampling if set to 'no')
    RestartNest= 'no' ;'yes' or 'no' Restart MCMC based on last MCMC.monitor.sav?
    NUM        = 50   ; OPTIONAL: THIS IS THE NUMBER OF LIVE POINTS 
    SAMPLE_NUM = 10   ; OPTIONAL: THIS IS THE NUMBER OF SAMPLES THAT IT TAKES BEFORE RECLUSTERING
    EXPAND     = 1.06d; OPTIONAL: THIS EXPANSION FACTOR EXPANDS THE SIZE OF THE ELLIPSOIDS
    TOL        = 0.5d0; OPTIONAL: THIS IS THE TOLERANCE/Accuracy set on measuring the Evidence LnZ
    NTHREADS   = 5    ; Number of parallel processors to sample with (should be < 1x number of processors)

    ;MCMC
    DoMCMC     = 'no' ; Perform MCMC fit?  (Skips MCMC if set to 'no')
    RestartMCMC= 'no'  ;'yes' or 'no' Restart MCMC based on last MCMC.monitor.sav?
    nthin      =   1   ; MCMC thinning parameter                                       (should work but still untested)
    nchains    =  25   ; Number of MCMC chains and Number of processor to run on  (should be < 1x number of processors/threads)
    maxsteps   = 30000 ; Number of MCMC steps for each chain  (MCMC runs until Gel.-Rub.<1.01 & 1000 independant draws)
    boostmcmc  = 0.8   ; Adjusts default DEMCMC step size (gamma = boostmcmc * 2.38d0/sqrt(2d0*nfit)). Use with low acceptance rates (e.g. <10%), set this number below the default 1.0 (e.g. 0.5).
                       ; Trial run is usefull to adjust 'boostmcmc' such that Acceptance rate is near ~20%
                       ; Usefull for when Marginalizing over unconstrained paramters (not likely needed)
    ;Data
    wsdstart   = 0.25  ; Perform Blue Cut in Wavelength of data if desired

    ;Atmosphere
    Cloudtype  = 0    ; Cloud is  0 = none                ; 1 = Uniform                 ; 2 = Hard Base                   ; 3 = Eddysed
    Hazetype   = 1    ; Haze  is  0 = none                ; 1 = Rayleigh               
    Chemistry  = 3    ; Chemistry 0 = Chemical Equibrium  ; 1 = Manual 'free' Chemistry ; 2 =  'Non-Equilbrium' Chemistry ; Chemistry 3 = Free-Element Chemical Equibrium 
    Isothermal = 0    ; Temperature is 0 = fit            ; 1 = Isothermal
    
    ;Condensation
    rain_switch = 1   ; rainout is on = 1                 ; rainout is off = 0

    ;Atmosphere Grid 
    nlevel    = 35    ; levels in the grid (grid is uniform in log-pressure)
    pmin      = 1E-6  ; Minimum pressure bar
    pmax      = 5E0   ; Maximum pressure bar
    nband     = 350   ; Number of correlated K-frequencies to add.  > 500 will load 5000 band data.
                      ; nband = 150  -> Low Resolution   500 band, wavelegnth range = 0.664 to 198 microns
                      ; nband = 350  -> Low Resolution   500 band, wavelegnth range = 0.286 to 198 microns
                      ; nband = 1100 -> High Resolution 5000 band, wavelegnth range = 0.904 to 198 microns

    ;Screen Plotting Option
    plotit    = 'no' ; 'yes' or 'no' Plot fit to screen during MPFIT
                     if (plotit eq 'no') then print,"RECOMEND starting XTERM IDL within 'screen' first"
                     if (plotit eq 'no') then print,''
   ;======================================================================================================

   ;======================================================================================================
   ;== (3) Select Parameters to Fit and Their Starting Value ==
   ;======================================================================================================
   ; SET FIXIT  0=free 1=fix 
   ; parameter index  0                1                     2                       3                     4                           5               6                    7             8                    9           10   
    parnameHist=   ['Temperature (K)','Abundance ln(Z/Zs)','log(g) cgs'            ,'Rp(1 mbar)'         ,'cloud opacity ln(sig/sig0)','haze opacity','cloud_top log(bar)','haze alpha', 'log10(Kzz cm^2/s)', 'fsed' , 'sigma_g']
    FIXIT  =       [ 1               , 0                  , 1                      , 0                   , 1                          , 0            ,     1              ,     0      ,     1              ,  1     ,      1   ]  
    StartingParams=[952.73d0         , -0.3925d0          , alog10(2.14D2)         ,1.38448D0            , 1d0                        , 2.971d0      ,     -1.0d0         ,    1.83d0  ,     9.5d0          ,  2.5   ,      2.0 ]   

   ;Bounding Box Priors
    LowPriorsParams=[ 400d0          , -2.2               , StartingParams[2]*0.5  ,StartingParams[3]*0.9 , -10                       ,  -10         ,  alog10(pmin)     ,        0.0  ,      5.0           ,  0.1   ,   0.1    ]   
    UpPriorsParams =[3000d0          , +4.6               , StartingParams[2]*2    ,StartingParams[3]*1.1 ,  10                       ,   10         ,  alog10(pmax)     ,         5.  ,      11.0          ,  10    ,    10    ] 

   ;==============================
   ; Fit T-P profile parameters (0=free) 
   ; Set alpha=0 and gamma2=0 for Three Param T-P fit.
   ; Set Isothermal = 0 to activate
   ;             11              12                 13                 14           15         16    17       18
   TP_names =  ['Kir'         ,'gamma1'           ,'gamma2'          ,'beta'      ,'alpha' , 'Tint','delta ','Pdelta  '   ]
   TP_params= [alog10(7.1D-2), alog10(0.03317D0), alog10(0.03317D0), 0.946d0    , 0.5d0    ,   500  , 0      , -4          ] ;Line et al. HD209
   ;TP_params=  [-2.2256445d0   ,   0.23395873d0   ,1.1685362d0       ,1.17296864d0, 0.0d0  ,  100  , 0      , -4          ];W121newmin                           
   TP_switch=  [          0    ,          0       ,        1         ,    0       ,   1   ,     1  , 1      ,  1          ]  
   ;Bounding Box Priors
   TP_priorLow=[  -5.0         , -4.0             ,  -4.0            , 0          ,  0    ,    99  , -1     ,alog10(pmin) ] 
   TP_priorUP= [  -0.5         ,  1.5             ,   1.5            , 2          ,  1    ,   800  ,  1     ,alog10(pmax) ]
   ;======================================================================================================

 
   ;==============================================================
   ; Chemistry = 1 = Manual 'free' Chemistry
   ; Set Opacities to include for transmission/emission spectra calculation 
   ; SET  Abund_fit_switch    0=free 1=fix 
   ;                   19   20   21    22    23   24    25    26  27   28   29   30   31    32    33   34    35    36    37    38    39   40
    Opacity_names   =['H2','He','H2O','CO2','CO','CH4','NH3','Na','K','Li','Rb','Cs','TiO','VO','FeH','PH3','H2S','HCN','C2H2','SO2','Fe','H_'] ; Switch:      1=include  0=Omit
    Opacity_switch  =[   1,   1,    1,    1,   1,    1,    1,   1,  1,   1,   0,   0,    1,   1,    1,   1,     1,    1,     1,    1,   1,   1] ; 
    Abund_fit_switch=[   1,   1,    0,    0,   0,    0,    1,   0,  0,   1,   1,   1,    1,   1,    1,   1,     1,    1,     1,    1,   1,   1] ;
   ; Set Abundances to include to vary (reqiures Chemistry=1 & Opacity_switch=1)
       ;==Option Start Abundances from Specific Values 
       ;==order is same as Abund_fit_switch, number of elements set to num of freely fit abundances
       ;==Uncomment line below to activate
       ;AbunBest=[exp(-6.4579836), exp(-4.6593652),exp(-20.636102),exp(-22.821792),exp(-13.396199),exp(-16.240493)] ;Volume Mixing Ratio
   ;Abundance Bounding Box Priors
    Abund_priorLow = [(dblarr(22)+1)*(alog(1D-14))]               ;minimum abundance VMR of 1x10^-12 for non-H/He 
    Abund_priorUP  = [[-0.05,-0.05],(dblarr(20)+1)*(alog(0.05))]  ;maximum abundance VMR of 0.05 for non-H/He 
   ;==============================================================

   ;==============================================================
   ; Chemistry = 3 Chemical Equibrium Non-Solar Elements
   ; SET Abund_fit_switch(i)    0=free 1=fix 
   ;                   19   20   21  22  23  24   25  26   27   28   29  30  31   32   33   34   35   36   37   38   39  40   41
   Element_names     =['H','He','C','N','O','Na','K','Si','Ar','Ti','V','S','Cl','Mg','Al','Ca','Fe','Cr','Li','Cs','Rb','F','P'] ; Elements
   Element_fit_switch=[   1,  1,  0,  1,  0,   0,  1,   1,   1,   1,  1,  1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, 1] ; Set to Vary Element Abundances
  

   EleBest = (dblarr(n_elements(Element_fit_switch)) +1.d0) * StartingParams(1); start all elements at same StartingParams metallicity factor as Abundance ln(Z/Zs)
   ;uncomment as needed to refine individual starting values
   ;EleBest(0)  = alog(2.) ;  Set H      ln(Z/Zs) metallicity
   ;EleBest(1)  = alog(2.) ;  Set He     ln(Z/Zs) metallicity
   EleBest(2)  = 3.94197d0 ;  Set Carbon ln(Z/Zs) metallicity
   ;EleBest(3)  = alog(2.) ;  Set N      ln(Z/Zs) metallicity
   EleBest(4)  = 3.4605d0 ;  Set Oxygen ln(Z/Zs) metallicity
   EleBest(5)  = 3.80191d0 ;  Set Na     ln(Z/Zs) metallicity
   ;EleBest(6)  = alog(2.) ;  Set K      ln(Z/Zs) metallicity
   ;EleBest(7)  = alog(2.) ;  Set Si     ln(Z/Zs) metallicity
   ;EleBest(8)  = alog(2.) ;  Set Ar     ln(Z/Zs) metallicity
   ;EleBest(9)  = alog(2.89832d0) ;  Set Ti     ln(Z/Zs) metallicity
   ;EleBest(10) = alog(20.7892d0) ;  Set V      ln(Z/Zs) metallicity 
   ;EleBest(12) = alog(2.) ;  Set S      ln(Z/Zs) metallicity 
   ;EleBest(13) = alog(2.) ;  Set Cl     ln(Z/Zs) metallicity 
   ;EleBest(14) = alog(2.) ;  Set Al     ln(Z/Zs) metallicity 
   ;EleBest(15) = alog(2.) ;  Set Ca     ln(Z/Zs) metallicity 
   ;EleBest(16) = alog(2.) ;  Set Fe     ln(Z/Zs) metallicity 
   ;EleBest(17) = alog(2.) ;  Set Cr     ln(Z/Zs) metallicity 
   ;EleBest(18) = alog(2.) ;  Set Li     ln(Z/Zs) metallicity 
   ;EleBest(19) = alog(2.) ;  Set Cs     ln(Z/Zs) metallicity 
   ;EleBest(20) = alog(2.) ;  Set Rb     ln(Z/Zs) metallicity 
   ;EleBest(21) = alog(2.) ;  Set F      ln(Z/Zs) metallicity 
   ;EleBest(22) = alog(2.) ;  Set P      ln(Z/Zs) metallicity 

   Opacity_switch  = Opacity_switch         ;Set this using above array
   Ele_priorLow = [(dblarr(23)+1)*(-2.2)]     ;minimum metallicity  -4.6 = 0.01x
   Ele_priorUP  = [(dblarr(23)+1)*(+4.6)]  ;maximum abundance VMR of 0.05 for non-H/He 
   ;==============================================================


 
   ;================================================================
   ;== (4) Planet and Star SYETEM PARAMETERS ==
   ;================================================================
   Rp    = 1.37D0               ; Radius of Planet (minimum of transit spectra)
   rstar = 1.39d0               ; Rsun  - Will be used in emission spectra calculation with (Fp/Fs)*(Rp/Rs)^2
   rorbit= 0.052                ; AU 
   stellarfile='../../stellar_spectra/WASP-127_lte058-4.0-0.0a+0.0.BT-Settl.7.ncdf'  ;location of stellar file
   Tstar = 6600D0                                                                    ;K
   Dist  = 19.4288              ; pc  distance to the star (used for DI spectra Rpl^2/Distance^2 flux dilution)
   ;================================================================

   ;====================================================================================================================
   ; **End of Inputs to EDIT**
   ;====================================================================================================================


   ;================================================================
   ;unit convertion
   Rjup  = 71492D5         ; cm
   Rpl   = Rp*Rjup         ; Rjup (cm)      
   Rs    = rstar*6.96D10   ;cm
   Distance = Dist*3.08567818585E13 ; Distance to star in km
   Distance = Distance/71492d0      ; Rjup units
   ;================================================================
   ;==============================
   ; Save Inputs
   Result = FILE_TEST(dirs, /DIRECTORY) & if (Result eq 0) then spawn,'mkdir '+dirs
   save,filename=dirs+'arc_inputs.sav',EleBest,Element_fit_switch,Ele_priorLow,Ele_priorUP,AbunBest,ABUND_FIT_SWITCH,CHEMISTRY,CLOUDTYPE,DATAFILE,DIRS,DOMCMC,EMISSION,FILENAM,FIXIT,HAZETYPE,$
                                       ISOTHERMAL,MAXSTEPS,NCHAINS,NLEVEL,NTHIN,OPACITY_NAMES,OPACITY_SWITCH,PARNAMEHIST,PLANET,PMAX,PMIN,RESTARTMCMC,$ 
                                       STARTINGPARAMS,TINT,TP_NAMES,TP_PARAMS,TP_SWITCH,WSDSTART,LowPriorsParams,UpPriorsParams,Abund_priorLow, Abund_priorUP,TP_priorLow,TP_priorUP,$
                                       NUM,SAMPLE_NUM,TOL,NTHREADS,EXPAND,PLOTIT,DONEST,Rp,Rjup,Rpl,Teff,logg,g,Rstar,Rs,rorbit,stellarfile,Tstar,Element_names,rain_switch,boostmcmc,Monitor,Distance,DoMPFIT,RestartNest,nband
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

