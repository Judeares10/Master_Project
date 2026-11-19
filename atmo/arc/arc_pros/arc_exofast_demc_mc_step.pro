pro arc_exofast_demc_mc_step,i,j,ChildProcessCount,pars,chi2,nthin,nchains,randomfunc,tofit,nfit,scale,chi2func,newpars,gamma,oldpars,oldchi2,olddet,naccept,seed=seed
; perform one step of mcmc chain, direct copy of exofast_demc j loop
; over nchains
;
; Input:  pars,chi2,nthin,nchains,seed,randomfunc,tofit,nfit,scale,chi2func,newpars
; Output: 
      naccept = 0;  DKS modify naccept

      oldpars[tofit] = pars[*,i-1,j]
      oldchi2 = chi2[i-1,j]
      ;; automatically thin the chain (saves memory)
      for k=0L, nthin-1L do begin      
         ;; differential evolution mcmc step -- see Ter Braak 2006
         repeat r1=floor(call_function(randomfunc,seed)*(nchains)) $
           until r1 ne j ;; a random chain
         repeat r2=floor(call_function(randomfunc,seed)*(nchains)) $
           until r2 ne j and r2 ne r1 ;; another random chain
         newpars[tofit] = oldpars[tofit] + $ 
                          gamma*(pars[*,i-1,r1]-pars[*,i-1,r2] + $
                          (call_function(randomfunc,seed,nfit)-0.5d0)*scale/10d0)
         ;; calculate the chi^2 of the new step
         newchi2 = call_function(chi2func,newpars,determinant=det)
         C = olddet[j]/det*exp((oldchi2 - newchi2)/2d0)
         ;; accept the step; update values
         if call_function(randomfunc,seed) lt C then begin
;            naccept++
            naccept=1  ;have to separate acceptance 
            olddet[j] = det
            oldpars = newpars
            oldchi2 = newchi2
         endif ;; else keep previous values
      endfor
      ;; Don't Update these yet, 
      ;; but output oldpars,oldchi2
      ;; update the link in the chain
      ;pars[*,i,j] = oldpars[tofit]
      ;chi2[i,j] = oldchi2 
      


end
