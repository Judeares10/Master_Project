pro arc_fitforwardmodel,A,wsd,r,rerr,werr,wmodel,fmodel,params,pcerror,chi,parname,MPFITfilenameprefix,FIXIT,parnameHist,parnameHist2,DoMCMC,nthin,nchains,maxsteps,RestartMCMC,DoNEST,NUM,SAMPLE_NUM,TOL,NTHREADS,EXPAND,plotit,boostmcmc,DoMPFIT,RestartNest

  Common ARC_retrieval_data, x,y,err,waverr
  Common SystemParameters,rplanet,rstar,rorbit,stellarfile,Tstar,Cloudtype,Hazetype,Chemistry,Isothermal,Opacity_switch,Abund_fit_switch,Abund_fit_index,Abund_priorLow,Abund_priorUP,rain_switch,Element_fit_switch,Element_fit_index,Ele_priorLow,Ele_priorUP,Distance,nband
  Common cpcount,ChildProcessCount,dirs,allmodel,allmodelpoints,Emission,nlevel,pmin,pmax,allchemistry,allTP
  common passmodel,ws,model,modelpoints,LowPriorsParams,UpPriorsParams,abundances,pressure,temperature

 chi=1E100 ;set to high value
 ;Setup MPFIT
 fa = {X:x, Y:y, ERR:err,waverr:waverr,plotit:plotit}  ;create structure of data to be fitt
 p0=A   ;initial guess
 parinfo = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0], step:0.D ,relstep:0.D}, n_elements(p0))
 parinfo[*].value=p0
 for i=0,n_elements(fixit)-1 do if (FIXIT(i) eq 1) then parinfo[i].fixed = 1 ;fix parameters 
 ;if (FIXIT(0) eq 1) then parinfo[0].fixed = 1         ; fixes Teff parameter

 ; Step Size used calculate numerical derivatives default=1D-7
  ;parinfo[0].STEP = 1  ; Temp: K
  ;parinfo[1].STEP = 0.01  ; MdH

; Enable MPFIT boundping-box priors on all parameters
for i=0,n_elements(p0)-1 do begin & $
    parinfo[i].LIMITED[1] = 1 & $ ;set upper limit
    parinfo[i].LIMITED[0] = 1 & $ ;set lower limit
    parinfo[i].LIMITS[1]  =  UpPriorsParams[i]   & $;upper limit of 10^X
    parinfo[i].LIMITS[0]  = LowPriorsParams[i]   & $;lower limit of 10^-X
endfor

 
 ;======

;  PARINFO Example:
;  parinfo = replicate({value:0.D, fixed:0, limited:[0,0], $
;                       limits:[0.D,0]}, 5)
;  parinfo(0).fixed = 1
;  parinfo(4).limited(0) = 1
;  parinfo(4).limits(0)  = 50.D
;  parinfo(*).value = [5.7D, 2.2, 500., 1.5, 2000.]
;  
;  A total of 5 parameters, with starting values of 5.7,
;  2.2, 500, 1.5, and 2000 are given.  The first parameter
;  is fixed at a value of 5.7, and the last parameter is
;  constrained to be above 50.
;

 tofit=where(parinfo.fixed eq 0)     ; Fitted parameters from MFIT
 tonotfit=where(parinfo.fixed ne 0)  ; Fixed  parameters from MFIT
 nfree=n_elements(tofit)
 N=n_elements(x)
 Dof     = N_ELEMENTS(X) - n_elements(where(parinfo.fixed ne 1)) ; deg of freedom

 ;Already ran MPFIT? skip if so or if deliberatly skipped
 result=ISA(DoMPFIT)          ;check if definded
 if (result eq 1) then params=p0 ;set all parameters to intital value
 if (result eq 1) then if (DoMPFIT eq 'no' or DoMPFIT eq 'No' or DoMPFIT eq 'NO') then goto,skipmpfit
 
 Result = FILE_TEST(dirs+'MPFIT.bestfit.sav')  
 if (Result eq 0) then begin
   params = mpfit('arc_model_chi_mpfit',functargs=fa,BESTNORM=bestnorm,PERROR=perror,PARINFO=parinfo,niter=niter,maxiter=maxer,status=status,Nfree=nfree,COVAR=COV)
   print,'MPFIT status:'
   if (status eq 0) then print,'improper input parameters. Check priors vs. starting values.'

    ;Save MPFIT spectra and Parameters
    if (Emission eq 'no')  then spawn,'cp ARC.transitspec'+string(ChildProcessCount,format='(I02)')+'.ncdf '+dirs+MPFITfilenameprefix+'.lowres.mpfit.transitspec.ncdf'
    if (Emission eq 'no')  then spawn,'cp ARC'+string(ChildProcessCount,format='(I02)')+'.in '+dirs+MPFITfilenameprefix+'.lowres.mpfit.transitspec.in'
    if (Emission eq 'yes') then spawn,'cp ARC.emissionspec'+string(ChildProcessCount,format='(I02)')+'.ncdf '+dirs+MPFITfilenameprefix+'.lowres.mpfit.emissionspec.ncdf'
    if (Emission eq 'yes') then spawn,'cp ARC'+string(ChildProcessCount,format='(I02)')+'.in '+dirs+MPFITfilenameprefix+'.lowres.mpfit.emissionspec.in'
    if (Emission eq 'yes') then spawn,'cp pt_ARC'+string(ChildProcessCount,format='(I02)')+'.ncdf '+dirs+MPFITfilenameprefix+'.emissionspec.lowres.mpfit.PTprofile.ncdf'
    ;save,filename=dirs+MPFITfilenameprefix+'.MPFIT.FitPars.sav',params,perror,bestnorm,chi,bic,nfree,dof,N,tofit,tonotfit

 endif

 if (Result eq 1) then begin & restore,dirs+'MPFIT.bestfit.sav' & PERROR=PCERROR & endif
   PCERROR = PERROR ;* SQRT(BESTNORM / DOF)   ; scaled uncertainties
   chi=BESTNORM
   BIC=BESTNORM+nfree*alog(n_elements(x))
   
   ;check if best-fitting solution is within mcmc boundry
   if (chemistry eq 1) then  fitted=where(Abund_fit_switch eq 0)+19 ;  select abundances that are being fit
   if (chemistry eq 1) then  fit=where(Abund_fit_switch eq 0) ;  select abundances that are being fit
   if (chemistry eq 1) then  for i=0,n_elements(fitted)-1 do if (params(fitted(i)) lt Abund_priorLow(fit(i))) then begin & params(fitted(i))=Abund_priorLow(fit(i)) & print,'Adjusted AbundL ' & endif
   if (chemistry eq 1) then  for i=0,n_elements(fitted)-1 do if (params(fitted(i)) gt Abund_priorUp(fit(i))) then begin &  params(fitted(i))=Abund_priorUp(fit(i))  & print,'Adjusted AbundL ' & endif


  ;Save MPFIT state for re-running purposes
  if (Result eq 0) then save,filename=dirs+'MPFIT.bestfit.sav',params,DoF,PCERROR,COV,nfree,status,parinfo,bestnorm,chi,bic,N,tofit,tonotfit,x,y,err,waverr,modelpoints,ws,model

  ;printscreen bestfit params
  ;tofit=where(parinfo.fixed eq 0)     ; Fitted parameters from MFIT
  for i=0,n_elements(tofit)-1 do print,parname(tofit(i)),' = ',params(tofit(i)),'+/-',pcerror(tofit(i))

  ;Print MPFIT stats to file
  fitparsfile=dirs+MPFITfilenameprefix+'.MCMC.FitPars.txt'
  openw,1,fitparsfile
    printf,1,'MPFIT'
    printf,1,'Chi^2=',BESTNORM,' DOF =',DOF,' BIC=',BIC,' Nfree=',nfree,' N=',n_elements(x)
    for i=0,n_elements(parname)-1 do  printf,1,parname(i),',',params(i),',',pcerror(i)
  close,1

    MPFITpars=pars
    MPFITPCERROR=PCERROR



    
    skipmpfit:
  ;=============================================
    pars=PARAMS  ;from MPFIT
    best=PARAMS  ;from MPFIT feeds into MCMC
    bestpars=pars;from MPFIT feeds into Nested 



  ;; ======= Nested Sampling =============================================
  if ((DoNEST eq 'no') or (DoNEST eq 'NO') or (DoNEST eq 'No')) then goto,skipNEST

	LOG_LIKELIHOOD_FUNC = 'arc_model_logZ_multinest' ;;REQUIRED: THIS IS THE LOG LIKELIHOOD FUNCTION USED BY MULTINEST
	PRIOR_FUNC = 'arc_multinest_prior';'example_prior' ;;REQUIRED: THIS IS THE PRIOR FUNCTION USED BY MULTINEST
	NUM_PARAMS=n_elements(tofit);Nfree ;;REQUIRED: THIS IS THE NUMBER OF FREE PARAMETERS

	;if not keyword_set(NUM) then NUM=30;0 ;;OPTIONAL: THIS IS THE NUMBER OF LIVE POINTS 
	;if not keyword_set(SAMPLE_NUM) then SAMPLE_NUM=20;0 ;;OPTIONAL: THIS IS THE NUMBER OF SAMPLES THAT IT TAKES BEFORE RECLUSTERING
	;EXPAND=1.2d ;;OPTIONAL: THIS EXPANSION FACTOR EXPANDS THE SIZE OF THE ELLIPSOIDS
                    ;Need to further expand the ellipsoid along each principal axis by a factor f to ensure the entire prior 
                    ;mass X(Li) is encompassed; this is required for unbiased results. Allison & Dunkley (2014, MNRAS) adopt f = 1.06 following Shaw et al. (2007). 
	;TOL=5;0.5d ;.2d ;;OPTIONAL: THIS IS THE TOLERANCE
	
	if (NTHREADS eq 1) then multinest_str=arc_multinest(LOG_LIKELIHOOD_FUNC,PRIOR_FUNC,NUM_PARAMS,NUM=NUM,$
							SAMPLE_NUM=SAMPLE_NUM,TOL=TOL,EXPAND=EXPAND,bestpars=bestpars,tofit=tofit,RestartNest=RestartNest) ;,NTHREADS=10,/USE_THREADS);,USE_THREADS=USE_THREADS)
	if (NTHREADS gt 1) then multinest_str=arc_multinest(LOG_LIKELIHOOD_FUNC,PRIOR_FUNC,NUM_PARAMS,NUM=NUM,$
							SAMPLE_NUM=SAMPLE_NUM,TOL=TOL,EXPAND=EXPAND,NTHREADS=NTHREADS,/USE_THREADS,bestpars=bestpars,tofit=tofit,RestartNest=RestartNest);,USE_THREADS=USE_THREADS)
        NESTpars=MULTINEST_STR.(0) ;all sampled parameters from nested sampling
        LP=MULTINEST_STR.(1)       ;ln weight of each sample
        save,filename=dirs+MPFITfilenameprefix+'.MULTINEST.sav',multinest_str,NESTpars,LP,LOG_LIKELIHOOD_FUNC,PRIOR_FUNC,NUM_PARAMS,NUM,SAMPLE_NUM,EXPAND,TOL,tofit	

        ;================================Save into to .txt file==
        fitparsfile=dirs+MPFITfilenameprefix+'.MCMC.FitPars.txt'
         openw,1,fitparsfile,/append
         printf,1,''
         printf,1,'Nested Sampling'
         printf,1,'Settings'
         printf,1,'  Live Points NUM : ',NUM
         printf,1,'  Recluster Sample Number SAMPLE_NUM : ',SAMPLE_NUM
         printf,1,'  LOGZ Tolerance : ',TOL
         printf,1,'Results'
         printf,1,'  Likelihood L uses Normal Distribution' 
         printf,1,'  lnL=-0.5D0*n*alog(2*pi)-total(alog(err))-total((y-modelpoints)^2./(2.0d* err^2.d))'
         printf,1,'      -0.5D0*n*alog(2D0*!pi)-total(alog(err)) = ',-0.5D0*N*alog(2D0*!pi)-total(alog(err))
         printf,1,'  ln Evidence = lnZ'
         printf,1,'  LnZ = '+string(multinest_str.(2)),'+/-',multinest_str.(3)
         printf,1,'  H = '+string(multinest_str.(4))+' information content (bits)'
         printf,1,'Parameters (need to reconstruct Posterior)'
         close,1
        ;make corer plot of results & show samples
        arc_monitor_nested_plot,dirs,chi2min ; bestfit model
        arc_monitor_nested_spectra,dirs      ; reconstruct model spectra posterior
        arc_monitor_nested,dirs              ; reconstruct posterior
        spawn,'python ../../arc/arc_pros/arc_monitor_nested_plot_pt.py '+dirs ; reconstruct and plot model T-P posterior (python program)
        spawn,'python ../../arc/arc_pros/arc_monitor_nested_plot_chem.py '+dirs ; reconstruct and plot model T-P posterior (python program)
        spawn,'python ../../arc/arc_pros/arc_monitor_nested_posterior.py '+dirs ; reconstruct posterior
        bestfitnestpars=nestpars(*,n_elements(nestpars(0,*))-1)  ;best-fitting parameter from nested sampling
        if (chi2min lt chi) then PARAMS(tofit) = bestfitnestpars ;set to bestfitting params if lowest chi2 is better
        if (chi2min lt chi) then begin
           best=PARAMS          ;from MPFIT feeds into MCMC
           chi=chi2min
           BIC=chi2min+nfree*alog(n_elements(x))
        endif
        
        ;================================
  skipNEST: 
  ;=============================================

  ;; ======= MCMC =============================================
  chi2func = 'arc_model_chi_MCMC'    ; 4 parameter Limb-darkening

  ;; Parameters to fit for the Transit fit
  ;  parinfo[5].fixed = 1
  tofit=where(parinfo.fixed eq 0)     ; Fitted parameters from MFIT
    if (n_elements(tofit) eq 1) then nchain=10 ;Avoid DEMC getting stuck if not enough chains defined
  ;help,tofit,pars


  if ((DoMCMC eq 'no') or (DoMCMC eq 'NO') or (DoMCMC eq 'No')) then goto,skipMCMC

  ;restore scale to seed 
  Result = FILE_test(dirs+'MCMCscale_seed.sav')
  if (Result eq 1) then restore,dirs+'MCMCscale_seed.sav';mcmcscale
  ;restore scale if already found
  Result = FILE_test(dirs+'MCMC.scale.sav')
  if (Result eq 1) then restore,dirs+'MCMC.scale.sav'
  ;Too Computationally Intensive
  ;if (Result eq 0) then scale = arc_exofast_getmcmcscale(bestpars,chi2func,tofit=tofit,seedscale=mcmcscale,angular=angular) ;don't run if already done
  if (Result eq 0) then begin
     ;initialize MCMC around best-fit solution, scale 5% or 1-sigma
     scale=abs(params(tofit)*0.05)
     smaller_scale=where(abs(pcerror(tofit)) gt 0 and abs(pcerror(tofit)) lt abs(params(tofit)*0.05))
     scale(smaller_scale)=abs(pcerror(tofit(smaller_scale)))
     print,'Inititizing MCMC parameters arround best-fit using scale: ',scale
  endif
  zeroscale=where(scale eq 0)
  if (zeroscale[0] ne -1) then scale[zeroscale]=0.1 ;set minimum scale if value was zero
  if (Result eq 0) then save,filename=dirs+'MCMC.scale.sav',scale ;save Meausred Scale
  ;save,filename='tempmcmc.sav',/all


;  if keyword_set(quiet) then begin
;   exofast_demc_quiet, best, chi2func, pars, chi2=chi2, tofit=tofit,$
;                 scale=scale, nthin=nthin,maxsteps=maxsteps,$
;                 burnndx=burnndx, seed=seed, randomfunc=randomfunc0,/REMOVEBURN
;  endif else begin
    ;exofast_demc, best, chi2func, pars, chi2=chi2, tofit=tofit,$
    ;             scale=scale, nthin=nthin,maxsteps=maxsteps,nchains=nchains,$
    ;             burnndx=burnndx, seed=seed, randomfunc=randomfunc0,/REMOVEBURN

    arc_exofast_demc_multicore, best, chi2func, pars, chi2=chi2, tofit=tofit,$
                 scale=scale, nthin=nthin,maxsteps=maxsteps,nchains=nchains,$
                 burnndx=burnndx, seed=seed, randomfunc=randomfunc0,/REMOVEBURN,RestartMCMC=RestartMCMC,boostmcmc=boostmcmc


   print, 'Synthesizing results; for long chains and/or many fitted parameters, this may take up to 15 minutes'

   allchains=pars
   allchainschi2=chi2

   ;; combine all chains
   sz = size(pars)
   npars = sz[1]
   nsteps = sz[2]
   nchains = sz[3]
   pars = reform(pars,npars,nsteps*nchains)
   chi2 = reform(chi2,nsteps*nchains)
   minchi2 = min(chi2,bestndx)
  ;; =============================================================
  ;; Save MCMC Chain
  save,filename=dirs+MPFITfilenameprefix+'.MCMC.sav',pars,chi2,minchi2,best,allchains,npars,nsteps,nchains,allchainschi2

  ;; =============================================================
  ;; TEST Convergence see exofast_demc for example
  ;; Test for convergence as outlined in Ford 2006
  ;; must be converged for 6 consecutive passes
  ;; tz > 1000 and Gelman-Rubin < 1.01 => converged 
     ;i=maxsteps-1 & j=nchains
     ;save,filename=dirs+'MCMC.monitor.sav',allchains,allchainschi2,i,j,nchains,nthin,tofit
     ;arc_density,DIRS
  ;; =============================================================


  ;================================
  ;save,filename='temp.sav',/variables,/comm

  fitparsfile=dirs+MPFITfilenameprefix+'.MCMC.FitPars.txt'
  openw,1,fitparsfile,/append
  printf,1,''
  printf,1,'MCMC'
  ;output fitted MCMC parameters, median and 1-sigma error (68%)
  for i=0,n_elements(tofit)-1 do begin & $
    e=percentiles(pars(i,*))-mean(pars(i,*)) & $
    medpar=median(pars(i,*)) & $
    printf,1,parname(tofit(i)),',',medpar,',',e(0),',',e(1)  & $ 
    ;cghistoplot,pars(i,*),binsize=(e(1)/6.),/fill,output=dirs+MPFITfilenameprefix+'.HistogramLin.'+parnameHist2(i)+'.pdf',DATACOLORNAME='black', POLYCOLOR='sky blue',xtitle=parnameHist(tofit(i)),/frequency  & $ 
    ;cghistoplot,alog(pars(i,*)),binsize=0.1,/fill,output=dirs+MPFITfilenameprefix+'.HistogramLog.'+parnameHist2(i)+'.pdf',DATACOLORNAME='black', POLYCOLOR='sky blue',xtitle=parnameHist(tofit(i)),/frequency  & $ 
 endfor
  ;output fixed parameters
  for i=0,n_elements(tonotfit)-1 do begin & $
    printf,1,parname(tonotfit(i)),',',params(tonotfit(i))  & $
  endfor
  close,1
  ;================================
  save,filename=dirs+MPFITfilenameprefix+'.MCMC.FitPars.sav',MPFITpars,MPFITPCERROR,chi,bic,nfree,dof,N,tofit,tonotfit

skipMCMC: 


end


