pro arc_redo_chimin
;Re-run old ARC plots to make a lower MCMC chi^2 the focus of some outputs
 
  Common ARC_retrieval_data, x,y,err,waverr
  Common SystemParameters,rplanet,logg,rstar,rorbit,stellarfile,Tstar,Cloudtype,Hazetype,Chemistry,Isothermal,Opacity_switch,Abund_fit_switch,Abund_fit_index,Abund_priorLow,Abund_priorUP,rain_switch,Element_fit_switch,Element_fit_index,Ele_priorLow,Ele_priorUP
  Common cpcount,ChildProcessCount,dirs,allmodel,allmodelpoints,Emission,nlevel,pmin,pmax
  common passmodel,model,modelpoints,LowPriorsParams,UpPriorsParams
   !path = '../../arc/arc_pros/:'+ !path                                      ;remember to change arc_multinest.pro threads setup
   !path = '../../arc/coyote/:'+ !path                                        ;remember to change arc_multinest.pro threads setup
   !path = '../../arc/arc_multinest_pros/:'+ !path                            ;remember to change arc_multinest.pro threads setup
   dirs='/data/sing/ARC_Results/ARC_Wasp12_eq_Hminus_5_Sept2018/'
   
    ;-=-If MCMC not finnished 
    ;restore,dirs+'Wasp12.Emission.COm214Sept2018arc_retrieval.setup.sav'
    ;restore,dirs+'MPFIT.bestfit.sav' & chi=bestnorm
    ;ChildProcessCount=0
    
    ;-=-If already finnished MCMC
    restore,dirs+'ARC.postprocessing.sav'


  if (DoMCMC eq 'yes') then begin & $
         restore,dirs+'MCMC.monitor.sav'  & $ ;allchains,allchainschi2,i,j,nchains,nthin,tofit
         ;Get min MCMC params 
           min_mcmc_chi=min(allchainschi2(0:i,*))  & $
           for k=0,nchains-1 do if min(allchainschi2(0:i,k)) eq min(allchainschi2(0:i,*)) then begin print,k,min(allchainschi2(0:i,k),minindex),minindex & minchain=k & endif  & $
           mcmc_params=allchains(*,minindex,minchain)  & $
         if (min_mcmc_chi lt chi) then params(tofit)=mcmc_params  & $ ;REPLACE MPFIT minimum parameters with MCMC determained values if chi^2 is lower
     endif
      ; set parameters index 0 to 9
     teff=params(0) & MdH=exp(params(1))   & COratio=exp(params(2))   & Rp=(params(3))   & Cloud=exp(params(4))   & Haze=exp(params(5))+1d0   & cloudtop=params(6)   & delta=params(7)   & Pdelta=10^params(8)   & Kzz=10^params(9)
     teffin=teff    & MdHin=exp(params(1)) & COratioin=exp(params(2)) & Rpin=(params(3)) & Cloudin=exp(params(4)) & Hazein=exp(params(5))+1d0 & cloudtopin=params(6) & deltain=params(7) & Pdeltain=10^params(8) & Kzzin=10^params(9)
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
       ;set first 5 T-P profile parameters (delta & pdelta are set in first 9)
       Kir=10^params(10) & gamma1=10^params(11) & gamma2=10^params(12) & beta=params(13)  & alpha=params(14)
       arc_tp_delta,ptfile,Kir,gamma1,gamma2,beta,alpha,T,Tint,Press,delta,Pdelta
       spawn,'cp pt_ARC'+string(ChildProcessCount,format='(I02)')+'.ncdf '+dirs+'pt_ARC_bestfit.ncdf'
     endif
     arc_make_atmo_in,teffin,MdHin,COratioin,Rpin,Cloudin,Hazein,fitabundances,cloudtop,Kzzin
     if (Emission eq 'no') then spawn,'cp ARC.transitspec'+string(ChildProcessCount,format='(I02)')+'.ncdf '+dirs+MPFITfilenameprefix+'.highres.transitspec.ncdf'
     if (Emission eq 'yes') then spawn,'cp ARC.emissionspec'+string(ChildProcessCount,format='(I02)')+'.ncdf '+dirs+MPFITfilenameprefix+'.highres.emissionspec.ncdf'
     spawn,'cp ARC'+string(ChildProcessCount,format='(I02)')+'.in '+dirs+MPFITfilenameprefix+'.highres.ARC.in'
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
end
