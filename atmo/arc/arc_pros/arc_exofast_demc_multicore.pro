;+
; NAME:
;   exofast_demc
; PURPOSE:
;   Make NCHAINS independent Markov Chain Monte Carlo chains to fit data
;   and return the parameter distributions.
;
; DESCRIPTION:
;   Begins by determining the correct stepping scale by finding
;   differences in each parameters that yield a delta chi^2 = 1, then
;   uses 5*scale*randomn offsets in each parameter for the starting
;   points of each chain.
;
;   It then begins a Differential Evolution Markov Chain Monte Carlo
;   fit (ter Braak, 2006)
;   http://www.stat.columbia.edu/~gelman/stuff_for_blog/cajo.pdf.  The
;   only slight modification to this basic alogorithm is that we
;   determine the magnitude of the uniform random deviate for each
;   parameter. We have found the dynamic range of our parameters can
;   be large, and thus not well suited to the one-size fits all
;   approach recommended there.
;
;   After taking 5% of MAXSTEPS, the program will estimate how many steps
;   will be necessary to be well-mixed. If it is not expected to be
;   well-mixed before taking the maximum number of steps, it will
;   output a warning and a recommended setting for NTHIN. This should
;   be accurate at the factor of 2-3 level.
;
;   The program stops when the chains are well-mixed, as defined by
;   Ford 2006 (http://adsabs.harvard.edu/abs/2006ApJ...642..505F)
;   using the Gelman-Rubin Statistic and number of independent draws,
;   or when each chain has taken MAXSTEPS, whichever is first.
;
;   Every step in the chain before all chains have crossed below the
;   median chi^2 will be considered the (the "burn-in"), and not used
;   for convergence tests. BURNNDX will specify the index of the first
;   useable chain so as not to be biased by the starting criteria.
;
; CALLING SEQUENCE:
;   exofast_demc, bestpars, 'mychi2', pars [,CHI2=, TOFIT=,$
;                 SCALE=, SEED=,RANDOMFUNC=,NTHIN=, MAXSTEPS=,/DONTSTOP,$
;                 NCHAINS=, ANGULAR=,BURNNDX=,/REMOVEBURN]
;
; INPUTS:
;   BESTPARS   - Array of the best-fit parameters. An accurate initial
;                fit is required to find the correct step size. MCMC
;                is not ideal for fitting data, just characterizing
;                the errors.
;   CHI2FUNC   - A string that specifies the name of the user-defined
;                function that calculates the chi^2 of a given
;                parameter set. Required data should be put in a
;                COMMON block. Its calling sequence must be:
;                chi2 = chi2func(pars, determinant=determinant)
;
; OPTIONAL INPUTS:
;   TOFIT      - An array that indexes the parameters to fit. If not
;                specified, all parameters are fit.
;   ANGULAR    - An array that indexes the parameters that are
;                angular (must be in radians). This will enable
;                special ways to calculate the median and convergence
;                statistics, which may otherwise fail. Default is none.
;   SCALE      - An NPARS array containing the stepping scales for
;                each parameter. This is not recommended for normal
;                use. If not specified, it will be automatically
;                determined using EXOFAST_GETMCMCSCALE.
;   NTHIN      - Saves only every NTHINth link in each chain. This
;                is only recommended when the autocorrelation between
;                steps is large and memory management issues may
;                arise. Results will be more accurate the smaller this
;                is -- as long as the chain is well-mixed (pay
;                attention to the warnings).
;   MAXSTEPS   - The maximum number of steps (after thinning) for each
;                chain to take.  The default is 100,000. Take care when
;                increasing this number, as memory allocation issues
;                may arise (especially with 32 bit machines). If the
;                chain is not well-mixed, it is usually better to
;                increase NTHIN instead.
;   NCHAINS    - The number of independent chains to run. The
;                execution time scales ~linearly with this number, and
;                must be at least 3 to calculated convergence
;                statistics. The default is 10. 
;                NCHAINS = Number of processors to run on.
;   SEED       - A seed for the random number generator. Do not mix
;                seeds from different random number
;                generators, and do not use several different seeds
;                for the same generator. If you use a random number
;                generator elsewhere, you should call EXOFAST_DEMC with
;                SEED. Default is -systime(/seconds). 
;   RANDOMFUNC - A string for the name of the random number generator
;                function. You can use your own as long as it returns
;                a uniform double precision deviate between 0 and 1,
;                accepts the /NORMAL keyword to generate random
;                gaussian deviates, and can return at least a two
;                dimensional array of deviates. You may also choose
;                between "EXOFAST_RANDOM" (default) or "RANDOMU".
;              - EXOFAST_RANDOM is an IDL wrapper for pg_ran. This is
;                a state of the art, bullet-proof generator with no
;                known statistical flaws, based on the third edition
;                of numerical recipes, pg 342. It is ~120x slower than
;                RANDOMU, though for typical DEMC fits, this accounts
;                for ~1% of the total runtime.
;              - RANDOMU is IDL's built-in, ran1-based generator with a
;                periodicity of ~10^9. ran1 is formally admonished by
;                its author, Numerical Recipes, but shown to give
;                identical results as the state of the art generator
;                See Eastman et. al., 2013 for more details.
;                http://adsabs.harvard.edu/abs/2013PASP..125...83E
;
; OPTIONAL KEYWORDS:
;   DONTSTOP   - If set, the program will not stop when it is
;                converged, but will instead take MAXSTEPS steps in
;                each chain.
;   REMOVEBURN - If set, the burn-in will be removed in the returned
;                parameter array.
; OUTPUTS:
;   PARS       - NFIT x NSTEPS x NCHAINS array of parameters.
;                NOTE: NSTEPS will be the lesser of when it is
;                converged or MAXSTEPS. 
;   CHI2       - NSTEPS x NCHAINS array with the chi^2 for each set of
;                parameters
;
; REVISION HISTORY:
;   2016/05 - David Sing (Exeter) multicore version of exofast_demc
;   2012/06 - Public Release - Jason Eastman (LCOGT)
;   2012/12 - When parameters aren't mixed, display the index
;             to the (more) static PARS array, not the index to the
;             dynamic TOFIT array.
;           - Removed memory-hungry vestigial code
;-

pro arc_exofast_demc_multicore, bestpars,chi2func,pars,chi2=chi2, tofit=tofit,$
                  scale=scale,seed=seed,randomfunc=randomfunc,$
                  nthin=nthin, maxsteps=maxsteps, dontstop=dontstop,$
                  nchains=nchains, angular=angular,$
                  burnndx=burnndx, removeburn=removeburn,RestartMCMC=RestartMCMC,boostmcmc=boostmcmc
  ;Pass ARC data along
  Common ARC_retrieval_data, x,y,err,waverr
  Common SystemParameters,rplanet,rstar,rorbit,stellarfile,Tstar,Cloudtype,Hazetype,Chemistry,Isothermal,Opacity_switch,Abund_fit_switch,Abund_fit_index,Abund_priorLow,Abund_priorUP,rain_switch,Element_fit_switch,Element_fit_index,Ele_priorLow,Ele_priorUP,Distance,nband
  Common cpcount,ChildProcessCount,dirs,allmodel,allmodelpoints,Emission,nlevel,pmin,pmax,allchemistry,allTP
  common passmodel,ws,model,modelpoints,LowPriorsParams,UpPriorsParams,abundances,pressure,temperature
  Compile_Opt DEFINT32
  ;  !! Don't forget to change common in Child PRocess below if editing Common

;; default values
if n_elements(maxsteps) eq 0 then maxsteps = 1d5
if n_elements(nthin) eq 0 then nthin = 1L
if n_elements(tofit) eq 0 then tofit = indgen(n_elements(bestpars))
nfit = n_elements(tofit)
!except = 0 ;; don't display errors for NaN and infinity

;; defaults for Differential Evolution MCMC
;adjust demc step size to help acceptance rate
gamma = (boostmcmc)*2.38d0/sqrt(2d0*nfit)
if n_elements(nchains) eq 0 then nchains=2d0*nfit

;; use robust generator if estimated steps is large 
;; NOTE: necessarily ignores any random numbers used by chi^2 routine
if n_elements(randomfunc) eq 0 then begin
   nfit = n_elements(tofit)
   ;; number of random deviates required
   ;; Gaussian deviates are made in pairs, and 4/pi must be made for each
   nrand = ((ceil(nfit/2d0)*2d0)*(4d0/!dpi)+1d0)*maxsteps*nchains*nthin
   if nrand le 5d7 then randomfunc = 'randomu' $
   else randomfunc = 'exofast_random'
endif

;; angular must be indexed within tofit to use with gelman-rubin
nang = n_elements(angular)
if nang gt 0 then gelmanangular=intarr(nang)
for i=0, nang-1 do gelmanangular[i] = (where(tofit eq angular[i]))(0)

;; initialize arrays
pars = dblarr(nfit,maxsteps,nchains)
chi2 = dblarr(maxsteps,nchains)
niter = 1L
newpars = bestpars
oldpars = bestpars
olddet = dblarr(nchains)
if (RestartMCMC eq 'yes') then begin 
  ;restore,dirs+'MCMC.monitor.sav';,allmodel,allmodelpoints,allchains,allchainschi2,i,j,nchains,nthin,tofit,naccept
  restore,dirs+'MCMC.Spectra.sav';,allmodel,allmodelpoints,allchains,allchainschi2,i,j,nchains,nthin,tofit,naccept,allchemistry,allTP
  spawn,'cp '+dirs+'MCMC.Spectra.sav '+dirs+'MCMC.Spectra.backup.sav'
  ;dirs+'MCMC.monitor.sav' ;               allchains,allchainschi2,i,j,nchains,nthin,tofit,naccept
  pars=allchains & chi2=allchainschi2
  istart=i
  for j=0,nchains-1 do olddet[j] = 1d0
endif

;; get the MCMC step scale for each parameter
if n_elements(scale) eq 0 then $
   scale = arc_exofast_getmcmcscale(bestpars,chi2func,tofit=tofit,angular=angular)
if scale[0] eq -1 then message, 'No scale found'

;; initialize each chain
if (RestartMCMC eq 'no') then begin
for j=0, nchains-1 do begin
   ;; repeat until a finite chi^2 
   ;; i.e., don't start outside a boundary; you'll never get back
   repeat begin
      ;; start 4 steps from best value (dks edit 2)
      pars[0:nfit-1,0,j] = bestpars[tofit] + $
                           2d0*scale*call_function(randomfunc,seed,nfit,/normal)
      newpars[tofit] = pars[0:nfit-1,0,j]

      ;; find the chi^2
      chi2[0,j] = call_function(chi2func, newpars, determinant=det)
      pars[*,0,j]=newpars[tofit] ;; in case chi2 function changes pars
      olddet[j] = det
   endrep until ((finite(chi2[0,j])) and (chi2[0,j] lt 1D100)) ;priors set chi2 to 1D308 if out of bounds
   allmodel[*,0,j]=model
   allmodelpoints[*,0,j]=modelpoints
   allchemistry[*,*,0,j] = abundances[*,Abund_fit_index-1]
   if (Isothermal eq 0) then allTP[*,0,j] = temperature

endfor
endif

nextrecalc = 100L
npass = 1L 
nstop = 0L
if (RestartMCMC eq 'no') then naccept = 1d0
mintz = 0d0
tzsteps = 0L
alreadywarned = 0L
t0 = systime(/seconds)

;r1s = lonarr(maxsteps,nchains)
;r2s = lonarr(maxsteps,nchains)

;Initialize Child Process Bridges, array of bridges obridge
;Number of Bridges (i.e. parallel processors) equals number of chains, nchains
  obridge=obj_new("IDL_IDLBridge", output='')  ;creates new IDL bridge on seperate processor
  status=intarr(nchains)                       ;array to track status of IDL bridges
  CD, Current=theDirectory                     ;current directory
  ;naccept=intarr(nchains)
  for i=0,nchains-2 do obridge=[obridge,obj_new("IDL_IDLBridge", output='')] ;makes array of bridges
  for i=0,nchains-1 do begin & $
      ;oBridge(i)->EXECUTE, '@' + PREF_GET('IDL_STARTUP') & $ ;IDL Setup 
      oBridge(i)->Execute,"cd, '"+theDirectory+"'" & $ ;change to working directory
      obridge(i)->Execute,'Common ARC_retrieval_data, x,y,err,waverr' & $
      obridge(i)->Execute,'Common SystemParameters,rplanet,rstar,rorbit,stellarfile,Tstar,Cloudtype,Hazetype,Chemistry,Isothermal,Opacity_switch,Abund_fit_switch,Abund_fit_index,Abund_priorLow,Abund_priorUP,rain_switch,Element_fit_switch,Element_fit_index,Ele_priorLow,Ele_priorUP,Distance,nband' 
      obridge(i)->Execute,'Common cpcount,ChildProcessCount,dirs,allmodel,allmodelpoints,Emission,nlevel,pmin,pmax,allchemistry,allTP' 
      obridge(i)->Execute,'Common passmodel,ws,model,modelpoints,LowPriorsParams,UpPriorsParams,abundances,pressure,temperature' 
      obridge(i)->Execute,"!path = '../../arc/arc_pros/:'+ !path"
      obridge(i)->Execute,"!path = '../../arc/arc_multinest_pros/:'+ !path"
      obridge(i)->setvar,'x',x & $  ;Load data
      obridge(i)->setvar,'y',y & $
      obridge(i)->setvar,'err',err & $
      obridge(i)->setvar,'waverr',waverr & $
      obridge(i)->setvar,'rplanet',rplanet & $
      ;obridge(i)->setvar,'logg',logg & $
      ;obridge(i)->setvar,'Rp',Rp & $
      obridge(i)->setvar,'rstar',rstar & $
      obridge(i)->setvar,'rorbit',rorbit & $
      obridge(i)->setvar,'stellarfile',stellarfile & $
      obridge(i)->setvar,'Tstar',Tstar & $
      obridge(i)->setvar,'Cloudtype',Cloudtype & $
      obridge(i)->setvar,'Hazetype',Hazetype & $
      obridge(i)->setvar,'Chemistry',Chemistry & $
      obridge(i)->setvar,'Isothermal',Isothermal & $
      obridge(i)->setvar,'Opacity_switch',Opacity_switch & $
      obridge(i)->setvar,'Abund_fit_switch',Abund_fit_switch & $
      obridge(i)->setvar,'Abund_fit_index',Abund_fit_index & $
      obridge(i)->setvar,'Abund_priorLow',Abund_priorLow & $
      obridge(i)->setvar,'Abund_priorUP',Abund_priorUP & $
      obridge(i)->setvar,'LowPriorsParams',LowPriorsParams & $
      obridge(i)->setvar,'UpPriorsParams',UpPriorsParams & $
      obridge(i)->setvar,'rain_switch',rain_switch & $
      obridge(i)->setvar,'Element_fit_switch',Element_fit_switch & $
      if (Chemistry eq 3) then obridge(i)->setvar,'Element_fit_index',Element_fit_index & $
      if (Chemistry eq 3) then obridge(i)->setvar,'Ele_priorLow',Ele_priorLow & $
      if (Chemistry eq 3) then obridge(i)->setvar,'Ele_priorUP',Ele_priorUP & $
      obridge(i)->setvar,'Distance',Distance & $
      obridge(i)->setvar,'nband',nband & $
      obridge(i)->setvar,'Emission',Emission & $
      obridge(i)->setvar,'nlevel',nlevel & $
      obridge(i)->setvar,'pmin',pmin & $
      obridge(i)->setvar,'pmax',pmax & $
      oBridge(i)->Execute,'.r arc_exofast_demc_mc_step' & $ ;Precompile routines
      oBridge(i)->Execute,'.r arc_model_chi_mcmc' & $
      oBridge(i)->Execute,'.r arc_make_atmo_in' & $
      oBridge(i)->Execute,'.r arc_load_atmo_emission_spec' & $
      oBridge(i)->Execute,'.r arc_tp' & $
      oBridge(i)->Execute,'.r linterp' & $
      oBridge(i)->Execute,'.r zparcheck' & $
      oBridge(i)->Execute,'.r tabinv' & $
      oBridge(i)->Execute,'.r isarray' & $
      oBridge(i)->Execute,"print,call_function('randomu',seed,/normal)" & $ ;initilize random number seed sequence
    endfor

 

;; start MCMC chain
for i=1L,maxsteps-1L do begin
   print,''
   ChildProcessCount=0 ;count to ensure separate file structure for each chain
   ;for nchains = number of processors (hard coded for now)
   ;Launch all chains, don't wait untill the last chain

   ; Relaunch an exisiting chain if requested
   if (RestartMCMC eq 'yes') then begin
     RestartMCMC='no' 
     i=istart+1
   endif

   for j=0L, nchains-1L do begin
       ;send mcmc data to each chain/process
       obridge(j)->setvar,'ChildProcessCount',ChildProcessCount
       obridge(j)->setvar,'dirs',dirs  ;Load data
       obridge(j)->setvar,'i', i
       obridge(j)->setvar,'j', j
       obridge(j)->setvar,'pars', pars
       obridge(j)->setvar,'chi2', chi2 
       obridge(j)->setvar,'nthin', nthin 
       obridge(j)->setvar,'nchains', nchains
       obridge(j)->setvar,'randomfunc',randomfunc 
       obridge(j)->setvar,'tofit', tofit
       obridge(j)->setvar,'nfit', nfit
       obridge(j)->setvar,'scale', scale
       obridge(j)->setvar,'chi2func', chi2func 
       obridge(j)->setvar,'newpars', newpars 
       obridge(j)->setvar,'gamma', gamma 
       obridge(j)->setvar,'oldpars', oldpars 
       obridge(j)->setvar,'olddet', olddet 
       ;Run MCMC Step
       if (j lt nchains-1) then oBridge(j)->Execute,'arc_exofast_demc_mc_step,i,j,ChildProcessCount,pars,chi2,nthin,nchains,randomfunc,tofit,nfit,scale,chi2func,newpars,gamma,oldpars,oldchi2,olddet,naccept,seed=seed',/nowait
       ; ** last Processes/Chain ** Make it Wait for this Chain by commenting ;,/nowait
       if (j eq nchains-1) then oBridge(j)->Execute,'arc_exofast_demc_mc_step,i,j,ChildProcessCount,pars,chi2,nthin,nchains,randomfunc,tofit,nfit,scale,chi2func,newpars,gamma,oldpars,oldchi2,olddet,naccept,seed=seed'
       ChildProcessCount++ ;=ChildProcessCount+1
    endfor

    ;Wait untill all chains/Processors are Completed before continuing
    ; Status values:  0=idle, 1=Exectuting Command, 2=Completed, 3=Error,4 = Aborted 
    for st=0,nchains-1 do status(st)=obridge(st)->status() ;get status of each bridge
    test=where(status eq 1) 
    ;when all processes are finnished, test(0) will = -1
    while(test(0) ne -1) do begin
       wait,0.5 ;wait 0.5 second before re-testing
       for st=0,nchains-1 do status(st)=obridge(st)->status() ;get status of each bridge
       test=where(status eq 1) 
        print,'Status:',status,' ',i
    endwhile

   ;  retrieve parameters from Bridges/Child Processes
   ;; update the link in the chain
   for j=0,nchains-1 do begin
     naccept     = naccept + obridge(j)->GetVar('naccept')
     olddet[j]   = obridge(j)->GetVar('olddet[j]')  
     oldpars     = obridge(j)->GetVar('oldpars')
     oldchi2     = obridge(j)->GetVar('oldchi2')
     pars[*,i,j] = oldpars[tofit]
     chi2[i,j]   = oldchi2
     if (obridge(j)->GetVar('naccept') gt 0) then allmodel[*,i,j]       = obridge(j)->GetVar('model')        else allmodel[*,i,j] = allmodel[*,i-1,j]
     if (obridge(j)->GetVar('naccept') gt 0) then allmodelpoints[*,i,j] = obridge(j)->GetVar('modelpoints')  else allmodelpoints[*,i,j] = allmodelpoints[*,i-1,j]
     if (obridge(j)->GetVar('naccept') gt 0 and chemistry ne 1) then allchemistry[*,*,i,j] = obridge(j)->GetVar('abundances[*,Abund_fit_index-1]') else allchemistry[*,*,i,j] = allchemistry[*,*,i-1,j]
     if (obridge(j)->GetVar('naccept') gt 0 and isothermal eq 0) then allTP[*,i,j] = obridge(j)->GetVar('temperature') 

  endfor

   allchains=pars & allchainschi2=chi2
   save,filename=dirs+'MCMC.monitor.sav',allchains,allchainschi2,i,j,nchains,nthin,tofit,naccept,Chemistry
   ;save spectra every 1%
   if (i mod round(maxsteps/1000d0)) eq 0 then save,filename=dirs+'MCMC.Spectra.sav',allmodel,allmodelpoints,allchains,allchainschi2,i,j,nchains,nthin,tofit,naccept,Chemistry,allchemistry,allTP,/COMPRESS
   ;; Test for convergence as outlined in Ford 2006
   ;; must be converged for 6 consecutive passes
   ;; tz > 1000 and Gelman-Rubin < 1.01 => converged 
   if i eq nextrecalc then begin
      ;; discard the burn-in (the first point at which all chains
      ;; have had at least one chi^2 lower than the median chi^2)
      medchi2 = median(chi2[0:i,*])
      burnndx = 0L
      for j=0L, nchains-1 do begin
         tmpndx = (where(chi2[0:i,j] lt medchi2))(0)
         if tmpndx gt burnndx then burnndx = tmpndx
      endfor
      ;; allows Gelman-Rubin calculation if one chain is being problematic
      burnndx = burnndx < (i-3) 
      
      ;; calculate the Gelman-Rubin statistic (remove burn-in)
      converged = exofast_gelmanrubin(pars[0:nfit-1,burnndx:i,*],$
                                      gelmanrubin,tz,angular=gelmanangular)
      
      ;; estimate the number of steps it will take until it is converged
      mintz = [mintz,min(tz)]
      tzsteps = [tzsteps,i]
      ntz = n_elements(mintz)
      if not alreadywarned then begin
         if i gt maxsteps/20d0 and ntz gt 3 then begin
            ;; fit a line to the number of independent draws
            coeffs = poly_fit(tzsteps[1:ntz-1],mintz[1:ntz-1],1,yfit=fit)
            
            ;; extrapolate to 2000, when it'll be considered converged 
            ;; 1000 is considered converged, but this is far from exact...
            stepstoconvergence = (2000 - coeffs[0])/coeffs[1]
            
            ;; if it won't be converged, warn the user
            if stepstoconvergence gt maxsteps then begin
               bestnthin = round(stepstoconvergence*nthin/maxsteps)
               print,'WARNING: The chain is not expected to be well-mixed; '+$
                     'set NTHIN to >~ ' + strtrim(bestnthin,2) + ' to fix'
            endif else print, string(stepstoconvergence*50d0/maxsteps,$
                     format='("EXOFAST_DEMC: The chain is expected to be ' + $
                                   'well-mixed after ",i2,"% completed     ")')
            alreadywarned = 1 ;; only display message once
         endif
      endif
      
      ;; must pass 6 consecutive times
      if converged then begin
         if nstop eq 0 then nstop = i
         nextrecalc = long(nstop/(1.d0-npass/100.d0))
         npass++
         if npass eq 6 then begin
            if not keyword_set(dontstop) then break
            nextrecalc = maxsteps
         endif
      endif else begin
         nextrecalc = long(i/0.9d0)
         nstop = 0L
         npass = 1L
      endelse
   endif
   
   ;; print the progress, cumulative acceptance rate, and time remaining
   acceptancerate = strtrim(string(naccept/double(i*nchains*nthin)*100,$
                                   format='(f6.2)'),2)
   timeleft = (systime(/seconds)-t0)*(maxsteps/(i+1d)-1d)
   units = ' seconds '
   if timeleft gt 60 then begin
      timeleft /= 60
      units = ' minutes '
      if timeleft gt 60 then begin
         timeleft /= 60
         units = ' hours   '
         if timeleft gt 24 then begin
            timeleft /= 24
            units = ' days    '
         endif
      endif
   endif

   ;; Windows formatting is messy, only output every 1%
   if !version.os_family ne 'Windows' or $
      i mod round(maxsteps/1000d0) eq 0 then begin
      timeleft = strtrim(string(timeleft,format='(f255.2)'),2)
      format='("EXOFAST_DEMC:",f6.2,"% done; acceptance rate = ",a,"%; ' + $ 
             'time left: ",a,a,$,%"\r")'
      print, format=format,100.d0*(i+1)/maxsteps,acceptancerate,timeleft,units
   endif

endfor
save,filename=dirs+'MCMC.Spectra.sav',allmodel,allmodelpoints,allchains,allchainschi2,i,j,nchains,nthin,tofit,naccept,allchemistry,allTP,/COMPRESS

print, '' ;; don't overwrite the final line
print, ''

;; if didn't converge, or user didn't want to stop use all steps
if npass ne 6 or keyword_set(dontstop) then nstop = maxsteps-1

;; discard the burn-in (the first point at which all chains
;; have had at least one chi^2 lower than the median chi^2)
medchi2 = median(chi2[0:nstop,*])
burnndx = 0L
for j=0L, nchains-1 do begin
   tmpndx = (where(chi2[0:nstop,j] lt medchi2))(0)
   if tmpndx gt burnndx then burnndx = tmpndx
endfor
burnndx = burnndx < (maxsteps-3) 

;; if hit maxsteps and never converged
if npass ne 6 then begin

   converged = exofast_gelmanrubin(pars[0:nfit-1,burnndx:nstop,*],$
                                   gelmanrubin,tz,angular=gelmanangular)
   bad = where(tz lt 1000 or gelmanrubin gt 1.01)

   if bad[0] ne -1 then begin
      message, 'WARNING: The Gelman-Rubin statistic indicates ' + $
               'the following parameters are not well-mixed', /continue
      print, '       Parameter  Gelman-Rubin Independent Draws'
      forprint, tofit[bad], gelmanrubin[bad],tz[bad],textout=1
   endif else message, 'WARNING: The chain did not pass 6 consecutive ' + $
                       'tests and may be marginally well-mixed.',/continue
   print, ''
endif

;; remove the burn-in/uncalculated parameters
;; keep them by default so people can do their own analysis
if keyword_set(removeburn) then begin
   pars = pars[*,burnndx:nstop,*]
   chi2 = chi2[burnndx:nstop,*]
endif else begin
   pars = pars[*,0:nstop,*]
   chi2 = chi2[0:nstop,*]
endelse

;; calculate the runtime and units
runtime = systime(/seconds)-t0
units = ' seconds'
if runtime gt 60 then begin
   runtime /= 60
   units = ' minutes'
   if runtime gt 60 then begin
      runtime /= 60
      units = ' hours'
      if runtime gt 24 then begin
         runtime /= 24
         units = ' days'
      endif
   endif
endif

runtime = strtrim(string(runtime,format='(f100.2)'),2)

format = '("EXOFAST_DEMC: done in ",a,a,"; took ",a,"% of trial steps")'
print, format=format, runtime, units, acceptancerate

end

    
