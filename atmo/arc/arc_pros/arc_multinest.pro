;NAME:
;	MULTINEST
;AUTHOR:
;	LUKE STAGNER 
;	UNIVERSITY OF CALIFORNIA, IRVINE
;	LSTAGNER@UCI.EDU
;VERSION: 1.0 
;PURPOSE:
;	MULTINEST PERFORMS MULTIMODEL NESTED SAMPLING
;DESCRIPTION:
;	MULTINEST USES X-MEANS CLUSTERING TO PARTITION THE SET OF LIVE POINTS. EACH CLUSTER IS THEN SURROUNDED BY AN 
;	ELLIPSOID. AN ELLIPSOID IS SELECTED ACCORDING TO THE NUMBER OF LIVE POINTS ENCLOSED. A SAMPLE IS THEN DRAWN 
;	FROM THE ELLIPSOID SUBJECT TO THE LIKELIHOOD CONSTRAINT. THIS IS REPEATED UNTIL A CERTAIN TOLERANCE .  
;CALLING SEQUENCE:
;	S = MULTINEST( LOG_LIKELIHOOD_FUNC,		$
;				   PRIOR_FUNC,				$
;				   NUM_PARAMS,				$
;				   NUM = num,				$
;				   SAMPLE_NUM = sample_num, $
;				   TOL = tol				$
;				   EXPAND = expand,         $
;				   PLOT = plot, 			$	
;				   SILENT=silent			)
;PARAMETERS:
;	LOG_LIKELIHOOD_FUNC: NAME OF LOG LIKELIHOOD FUNCTION OF THE FORM: VALUE=LOG_LIKELIHOOD_FUNC(ARRAY[NUM_PARAMS]) 
;						 (LOG LIKELIHOOD FUNCTION RETURNS THE LOG OF THE LIKELIHOOD FUNCTION GIVEN A SET OF PARAMETERS)
;	PRIOR_FUNC: NAME OF PRIOR FUNCTION OF THE FORM: ARRAY[NUM_PARAMS,NUM] =  PRIOR_FUNC(ARRAY[NUM_PARAMS,NUM])
;				(PRIOR FUNCTION MAPS UNIFORM RANDOM SAMPLES FROM HYPERCUBE TO WANTED PRIOR PROBABILITY DISTRIBUTIONS)
;	NUM_PARAMS: THE NUMBER OF PARAMETERS BEING ESTIMATED
;KEYWORDS:
;	NUM: NUMBER OF LIVE POINTS DEFAULT SET TO 1000
;	SAMPLE_NUM: NUMBER OF SAMPLES DRAWN BETWEEN CLUSTERING DEFAULT SET TO 1
;	TOL: LOG EVIDENCE TOLERANCE DEFAULT SET TO .001 IN LOG EVIDENCE 
;	EXPAND: EXPANSION FACTOR FOR SAMPLING ELLIPSOIDS 
;	PLOT: SHOWS HISTOGRAMS FOR EACH PARAMETER BEING ESTIMATED
;RETURN VALUE:
;	A STRUCTURE WHICH CONTAINS THE SAMPLES DRAWN ACCORDING TO THE POSTERIOR DISTRIBUTION, THE INFORMATION, THE LOG_EVIDENCE, 
;	ASSOCIATED LOG_EVIDENCE ERROR, AND LOG_PROBABILITY FOR EACH SAMPLE
;REQUIRED PROGRAMS:
;	CLUSTER_DATA.PRO
;	SAMPLE_ELLIPSOID.PRO
;	IN_ELLIPSOIDS.PRO
;	RANDOMNUMBERGENERATOR__DEFINE.PRO
;	PLOT_2D_CLUSTERS.PRO
;COMMON BLOCKS:
;	NONE 
;REFERENCES:
;F. FEROZ, M.P. HOBSON, "MULTIMODAL NESTED SAMPLING: AN EFFICIENT AND ROBUST ALTERVATIVE TO MCMC METHODS FOR ASTRONOMICAL DATA ANALYSIS"
;	Mon. Not. Roy. Astron. Soc., 384, 2, 449-463 (2008)
;
;F. FEROZ, M.P. HOBSON, M. BRIDGES, "MULTINEST: AN EFFICIENT AND ROBUST BAYESIAN INFERENCE TOOL FOR COSMOLOGY AND PARTICLE PHYSICS" 
;	Mon. Not. Roy. Astron. Soc. 398: 1601-1614 (2009)
;DATE MODIFIED:
;	08/06/2012
;       01/23/2017 added the necessary setup for the threads

FUNCTION arc_multinest,LOG_LIKELIHOOD_FUNC,PRIOR_FUNC,NUM_PARAMS,$
				   NUM=NUM,SAMPLE_NUM=SAMPLE_NUM,TOL=TOL,EXPAND=EXPAND,$
				   PLOT=PLOT,SILENT=SILENT,NTHREADS=NTHREADS,USE_THREADS=USE_THREADS,bestpars=bestpars,tofit=tofit,RestartNest=RestartNest
  Common ARC_retrieval_data, x,y,err,waverr
  Common SystemParameters,rplanet,rstar,rorbit,stellarfile,Tstar,Cloudtype,Hazetype,Chemistry,Isothermal,Opacity_switch,Abund_fit_switch,Abund_fit_index,Abund_priorLow,Abund_priorUP,rain_switch,Element_fit_switch,Element_fit_index,Ele_priorLow,Ele_priorUP,Distance,nband
  Common cpcount,ChildProcessCount,dirs,allmodel,allmodelpoints,Emission,nlevel,pmin,pmax,allchemistry,allTP
  common passmodel,ws,model,modelpoints,LowPriorsParams,UpPriorsParams,abundances,pressure,temperature
	
	DefSysV, '!RNG', Obj_New('RandomNumberGenerator')
	if not keyword_set(NUM) then NUM=1000
	if not keyword_set(SAMPLE_NUM) then SAMPLE_NUM=1
	if not keyword_set(EXPAND) then EXPAND=1.0d
	if not keyword_set(TOL) then TOL=0.001d
	if not keyword_set(NTHREADS) then NTHREADS=4

        CD, Current=theDirectory ;current directory
	if keyword_set(USE_THREADS) then begin
		;;Create threads
                ;RESOLVE_ROUTINE,'thread_callback' ;sometimes (not always) crashes running create_threads... try pre compiling first
                ;RESOLVE_ROUTINE,'CREATE_THREADS'
                RESOLVE_ROUTINE,'WAIT_THREADS' 
                RESOLVE_ROUTINE,'UPDATE_STATUS'
		;create_threads,NTHREADS,threads,/idl_startup
         	stat_str=['idle','running','completed','error','aborted']
	        threads=objarr(nthreads)
	        stat=intarr(nthreads)
                threads=obj_new("IDL_IDLBridge", output='')                                 ;creates new IDL bridge on seperate processor
                for i=0,NTHREADS-2 do begin  & $
                   threads=[threads,obj_new("IDL_IDLBridge", output='')]  & $;makes array of bridges
                   stat[i]=threads[i]->status()  & $
                endfor
                for i=0,NTHREADS-1 do begin  & $
                   threads(i)->EXECUTE, '@' + PREF_GET('IDL_STARTUP') & $ ;IDL Setup 
                   threads(i)->Execute,"cd, '"+theDirectory+"'" & $    ;change to working directory
                   threads(i)->Execute,'Common ARC_retrieval_data, x,y,err,waverr' & $
                   threads(i)->Execute,'Common SystemParameters,rplanet,rstar,rorbit,stellarfile,Tstar,Cloudtype,Hazetype,Chemistry,Isothermal,Opacity_switch,Abund_fit_switch,Abund_fit_index,Abund_priorLow,Abund_priorUP,rain_switch,Element_fit_switch,Element_fit_index,Ele_priorLow,Ele_priorUP,Distance,nband' & $
                   threads(i)->Execute,'Common cpcount,ChildProcessCount,dirs,allmodel,allmodelpoints,Emission,nlevel,pmin,pmax,allchemistry,allTP' & $
                   threads(i)->Execute,'common passmodel,ws,model,modelpoints,LowPriorsParams,UpPriorsParams,abundances,pressure,temperature' & $
                   threads(i)->Execute,"!path = '../../arc/arc_pros/:'+ !path" & $ 
                   threads(i)->Execute,"!path = '../../arc/arc_multinest_pros/:'+ !path" & $ 
                   threads(i)->Execute,"!path = '../../arc/coyote/:'+ !path" & $ 
                   threads(i)->setvar,'initialSeed',Systime(1)+i*100 & $ ;need to give different seeds for parallel processors
                   threads(i)->Execute,"DefSysV, '!RNG', Obj_New('RandomNumberGenerator',initialSeed) & !EXCEPT=0" & $ 
                   state={id:i,status:stat[i],status_str:stat_str[stat[i]]} & $ 
		   threads[i]->setproperty,userdata=state & $ 
	           threads(i)->setvar,'x',x & $  ;Load data
	           threads(i)->setvar,'y',y & $
 	           threads(i)->setvar,'err',err & $
	           threads(i)->setvar,'waverr',waverr & $
	           threads(i)->setvar,'rplanet',rplanet & $
	           ;obridge(i)->setvar,'logg',logg & $
	           ;obridge(i)->setvar,'Rp',Rp & $
	           threads(i)->setvar,'rstar',rstar & $
	           threads(i)->setvar,'rorbit',rorbit & $
	           threads(i)->setvar,'stellarfile',stellarfile & $
 	           threads(i)->setvar,'Tstar',Tstar & $
 	           threads(i)->setvar,'Cloudtype',Cloudtype & $
 	           threads(i)->setvar,'Hazetype',Hazetype & $
	           threads(i)->setvar,'Chemistry',Chemistry & $
 	           threads(i)->setvar,'Isothermal',Isothermal & $
	           threads(i)->setvar,'Opacity_switch',Opacity_switch & $
	           threads(i)->setvar,'Abund_fit_switch',Abund_fit_switch & $
	           threads(i)->setvar,'Abund_fit_index',Abund_fit_index & $
 	           threads(i)->setvar,'Abund_priorLow',Abund_priorLow & $
	           threads(i)->setvar,'Abund_priorUP',Abund_priorUP & $
	           threads(i)->setvar,'LowPriorsParams',LowPriorsParams & $
	           threads(i)->setvar,'UpPriorsParams',UpPriorsParams & $
	           threads(i)->setvar,'rain_switch',rain_switch & $
	           threads(i)->setvar,'Element_fit_switch',Element_fit_switch & $
                   if (Chemistry eq 3) then threads(i)->setvar,'Element_fit_index',Element_fit_index & $
                   if (Chemistry eq 3) then threads(i)->setvar,'Ele_priorLow',Ele_priorLow & $
                   if (Chemistry eq 3) then threads(i)->setvar,'Ele_priorUP',Ele_priorUP & $
                   threads(i)->setvar,'Distance',Distance & $
                   threads(i)->setvar,'nband',nband & $
                   threads(i)->setvar,'Emission',Emission & $
                   threads(i)->setvar,'nlevel',nlevel & $
                   threads(i)->setvar,'pmin',pmin & $
                   threads(i)->setvar,'pmax',pmax & $
                   threads(i)->setvar,'bestpars',bestpars & $
                   threads(i)->setvar,'tofit',tofit & $
                   threads(i)->setvar,'dirs',dirs & $
                   threads(i)->setvar,'ChildProcessCount',i & $ ;fix childprocess count to thread number, keeps parallel ATMO runs separate
                   ;threads(i)->Execute,'.r arc_exofast_demc_mc_step' & $ ;Precompile routines
                   threads(i)->Execute,'.r arc_model_logz_multinest' & $
                   threads(i)->Execute,'.r arc_make_atmo_in' & $
                   threads(i)->Execute,'.r arc_load_atmo_emission_spec' & $
                   threads(i)->Execute,'.r arc_tp' & $
                   threads(i)->Execute,'.r linterp' & $
                   threads(i)->Execute,'.r zparcheck' & $
                   threads(i)->Execute,'.r tabinv' & $
                   threads(i)->Execute,'.r isarray' & $
                   threads(i)->Execute,'.r sample_ellipsoid' & $
                   threads(i)->Execute,'.r get_sample' & $
                   threads(i)->Execute,'.r arc_multinest_prior' & $
                   threads(i)->Execute,"print,call_function('randomu',seed,/normal)" & $ ;initilize random number seed sequence
                endfor
	endif

        openw,1,dirs+'multinest.monitor.log.txt',/append
           printf,1,' Julian Date UTC           DeltaZ           logZ'
           printf,1,string(systime(/Julian,/utc),format='(f19.7)')
        close,1
        
	COMMON mcmc_common,accept,tot,sigma
	accept=0.0d & tot=0.0d  & sigma=.01d
        threadtrack=0
        counttrack=0
        
	;;Get initial samples
        ;Explanation of Pointer Arrays http://www.idlcoyote.com/code_tips/ptrarray.html
        live_points=ptrarr(NUM,/allocate_heap)
        live_model=ptrarr(NUM,/allocate_heap)      ;dks edit
        live_abundances=ptrarr(NUM,/allocate_heap) ;;dks edit
        live_temperature=ptrarr(NUM,/allocate_heap);;dks edit
	for i=0,NUM-1 do begin  & $
		*live_points[i] = !RNG -> GetRandomNumbers(NUM_PARAMS,/double)  & $
	endfor

	live_samples=CALL_FUNCTION(PRIOR_FUNC,live_points,tofit)

	print,"GOT INITIAL LIVE SAMPLES"
	;;Get initial likelihoods and likelihood min
	live_likelihoods=dblarr(NUM)

        if (RestartNest eq 'yes' or RestartNest eq 'Yes' or RestartNest eq 'YES') then goto,restartnested        
        ; Load same intital Samples into all threads
        ; Can't pass pointer vectors to threads, do so by .sav file instead
        save,filename=dirs+'Multinest_initial_samples.sav',live_samples,live_likelihoods,live_points,live_model,live_abundances,live_temperature
	for j=0,NTHREADS-1 do begin  & $
           threads(j)->setvar,'LOG_LIKELIHOOD_FUNC',LOG_LIKELIHOOD_FUNC   & $
           cmd="restore,'"+dirs+"Multinest_initial_samples.sav'"   & $
           threads(j)->Execute,cmd   & $
        endfor
        status=intarr(NTHREADS)
        ;Calculate the initial Likelihoods for all points with parallel processing
	for i=0L, NUM-1 do begin
            ;launch likelihood calculation on the threads in parallel
              while (i/NTHREADS le num/nthreads-1) do begin   ;while (i+NUM/NTHREADS lt NUM-1) do begin
              for j=0,NTHREADS-1 do begin
                ;live_likelihoods[i]=CALL_FUNCTION(LOG_LIKELIHOOD_FUNC,live_samples[i])
	        ;execute_threads,threads,cmd
                threads(j)->setvar,'i',i ;pass index variable
                threads(j)->Execute,'live_likelihoods[i]=CALL_FUNCTION(LOG_LIKELIHOOD_FUNC,live_samples[i],bestpars,tofit)',/nowait
                print,'Execting model ',i,j
                ;if (j eq NTHREADS-1) then threads(j)->Execute,'live_likelihoods[i]=CALL_FUNCTION(LOG_LIKELIHOOD_FUNC,live_samples[i],bestpars,tofit)'
                i=i+1; 
              endfor
              ;Wait untill all chains/Processors are Completed before continuing
              ; Status values:  0=idle, 1=Exectuting Command, 2=Completed, 3=Error,4 = Aborted 
              for st=0,NTHREADS-1 do status(st)=threads(st)->status() ;get status of each bridge
              test=where(status eq 1) 
              ;when all processes are finnished, test(0) will = -1
              while(test(0) ne -1) do begin
                 wait,0.5                                     ;wait 0.5 second before re-testing
                 for st=0,NTHREADS-1 do status(st)=threads(st)->status() ;get status of each bridge
                 test=where(status eq 1) 
                 print,'Status:',status,' ',i
              endwhile
	      ;wait_threads,threads ;wait for parallel threads to all finnnish

              ;pull_variable,threads,out_vars,output
              for j=0,NTHREADS-1 do begin
                print,'Getting values for ',j
                live_likelihoods[i-NTHREADS+j]    = threads[j]->getvar('live_likelihoods[i]')
                *live_model[i-NTHREADS+j]         = threads[j]->getvar('model')
                *live_abundances[i-NTHREADS+j]    = threads[j]->getvar('abundances')
                *live_temperature[i-NTHREADS+j]   = threads[j]->getvar('temperature')
             endfor
            endwhile
            ;finnish off remainder not in parallel
            if (i lt NUM) then begin
               live_likelihoods[i]  = CALL_FUNCTION(LOG_LIKELIHOOD_FUNC,live_samples[i],bestpars,tofit)
               *live_model[i]       = model
               *live_abundances[i]  = abundances
               *live_temperature[i] = temperature
            endif
         endfor
             ;Get very first spectra, find number of points, define array to hold all the spectra
         bst=where(live_likelihoods gt -1E100)
         bst=reverse(bst)
         test=j-NTHREADS ;backup by number of threads
         while (live_likelihoods(test) lt -1E100) do test=test+1  ;pick a model that actually ran
            ws=threads[test]->getvar('ws')  ;get wavelength solution 
            pressure=threads[test]->getvar('pressure')           ; PRESSURE        DOUBLE    = Array[36]
            abundances=threads[test]->getvar('abundances')       ; ABUNDANCES      DOUBLE    = Array[36, 277]
            NUM_WS=n_elements(ws)    ;number of wavelength points
            NUM_PRESS=n_elements(pressure) ;number of pressure points
            NUM_SPECIES=n_elements(abundances(0,*))
         
         ;add molnames to output
            molname=''
            Result = FILE_TEST('chem_eq_ARC'+string(ChildProcessCount,format='(I02)')+'.ncdf')
            if (Result eq 1) then begin
               Fid=ncdf_open('chem_eq_ARC'+string(ChildProcessCount,format='(I02)')+'.ncdf')  ;test if file is ARC output, if not then load as is
               ncdf_varget,fid,'molname',molname
               ncdf_close,fid
               molname=string(molname)
            endif
            
            
        save,filename=dirs+'Multinest_initial_samples.sav',live_samples,live_likelihoods,live_points,live_model,live_abundances,live_temperature,ws,pressure,molname

	likelihood_min=MIN(live_likelihoods,likelihood_min_loc,MAX=likelihood_max)
	
	;;Initialize Varibles	
	H = 0.0
	logZ=-(machar(/double)).xmax
	k=0L
	logw=alog(1.0-exp(-1.0/double(NUM)))
	deltaZ = TOL

	samples=ptrarr(10000)
	log_prob=dblarr(10000,/nozero)
        n_samp=10000L
        samples_model=ptrarr(10000)
        samples_abundances=ptrarr(10000)
        samples_temperature=ptrarr(10000)

        ;-----load restarted MCMC values----
        restartnested:
        if (RestartNest eq 'yes' or RestartNest eq 'Yes' or RestartNest eq 'YES') then begin
                restore,dirs+'Multinest_current_samples.sav' ;,live_samples,live_likelihoods,live_points,live_model,live_abundances,live_temperature,ws,pressure,H,logZ,k,logw,deltaZ,samples,log_prob,n_samp,samples_model,samples_abundances,samples_temperature,likelihood_min,NUM_PARAMS,NUM_WS,NUM_PRESS,NUM_SPECIES
           	likelihood_min=MIN(live_likelihoods,likelihood_min_loc,MAX=likelihood_max)
   	        for j=0,NTHREADS-1 do begin
                  threads(j)->setvar,'LOG_LIKELIHOOD_FUNC',LOG_LIKELIHOOD_FUNC
                  cmd="restore,'"+dirs+"Multinest_current_samples_threads.sav'";live_samples,live_likelihoods,live_points,live_model,live_abundances,live_temperature
                  threads(j)->Execute,cmd
                endfor
                print,"**************************"
                print,"Restarting Nested Sampling"
                print,"Current deltaZ      : ",deltaZ
                print,"Current Z Tolerance : ",TOL
                print,"Samples Taken       : ",k         
                print,"**************************"
        endif
        ;--------------------------------

        ;===============================================================================================
	;;Begin Nested Sampling Loop
	while (deltaZ ge TOL) do begin
;		print,'Clustering Data'
		;;Cluster the live points 		
		result=cluster_data(live_points,/k_medoids)
		
		;;Determine the number of clusters
		clusters=n_elements(result)
		print,'Number of clusters:',clusters
		;;Determine the number of points in each cluster
		ellnum=dblarr(clusters)
		for i=0,clusters-1 do begin
			ellnum[i]=result[i].data_num
		endfor
		totnum=total(ellnum)
		cnt=0L

		if keyword_set(plot) then begin
			if NUM_PARAMS eq 2 then begin
				no_ellipse=0
				plot_2d_clusters,result,no_ellipse=no_ellipse
			endif else plot_3d_clusters,result
		endif
		index=0
		while cnt le SAMPLE_NUM-1 do begin
			hit=0.0 & miss=0.0
			;;Pick an ellipse to sample from based on the number of points in the ellipses
			ellprob=total(double(ellnum)/double(totnum),/cumulative)
			r=!RNG->GetRandomNumbers(1,/long)
			start_point=live_points[(r mod NUM)]
			w=in_ellipsoids(live_points[(r mod NUM)],result, expand=EXPAND, /location)
			if w[0] eq -1 then continue
			r=!RNG->GetRandomNumbers(1,/long) ;;pick new random number to be independent of the choice of ellipsoid
			if n_elements(w) gt 1 then w[0]=w[r mod n_elements(w)]
			cov=result[w[0]].cov
			mean=result[w[0]].mean
			sevecs=result[w[0]].sevecs ;Scaled Eigenvectors                        
			k_fac_max=result[w[0]].kmax ;Ellipsoid surface constant i.e. (x/a)^2 + (y/b)^2 + ... == kmax
			if keyword_set(USE_THREADS) then begin
				in_vars=['PRIOR_FUNC','LOG_LIKELIHOOD_FUNC','cov','mean','k_fac_max','EXPAND','sevecs']
				out_vars=['samp','tsamp','samp_likelihood','model','abundances','temperature']
				push_variable,threads,in_vars ;load needed variables into threads 
				wait_threads,threads
			endif
			;;Sample in Parallel
			samp_likelihood=likelihood_min-1
			count=0L
			while 1 do begin ;repeat indefintly until samp_likelihood gt likelihood_min
				if count lt 100 then begin ;do 100 random tries
                                        ;random sample in parallel to find samp_likelihood gt likelihood_min
                                        if keyword_set(USE_THREADS) then begin
                                                if index eq 0 then begin
                                                        print,'Current likelihood min = ',likelihood_min
							cmd='get_sample,PRIOR_FUNC,LOG_LIKELIHOOD_FUNC,cov,mean,k_fac_max,EXPAND,sevecs,samp,tsamp,samp_likelihood,bestpars,tofit,model,abundances,temperature'
							execute_threads,threads,cmd
							wait_threads,threads
							pull_variable,threads,out_vars,output
						endif
						tmp=index
						for i=tmp,NTHREADS-1 do begin
							if clusters eq 1 then index+=1
							samp_likelihood= (*(output[i,2]))
							if samp_likelihood gt likelihood_min then begin
								samp=output[i,0]   ; sampled points (0 to 1)
                                                                tsamp=output[i,1]  ; sampled parameters
                                                                modelsamp=output[i,3]     ;model which is a pointer
                                                                abundancessamp=output[i,4] ;abundances which is a pointer 
                                                                temperaturesamp=output[i,5] ;temperture which is a pointer
                                                                print,'Likelihood Live Point with min ',likelihood_min,' updated from thread: ',i
                                                                counttrack=[[counttrack],count]
                                                                threadtrack=[[threadtrack],i]
								if clusters eq 1 then print,index
								break ;end while 1 do 
							endif
						endfor
						if index ge NTHREADS then index=0
                                             endif else begin
                                          ;if not in parallel sample one at a time  
					  get_sample,PRIOR_FUNC,LOG_LIKELIHOOD_FUNC,cov,mean,k_fac_max,EXPAND,sevecs,samp,tsamp,samp_likelihood,bestpars,tofit,model,abundances,temperature,/ptr
					endelse
					count+=1
				endif else begin
;					print,'mcmc start'
                                        ;100 random tries were not enough, do short mcmc on liklihood_min point
                                        ;mcmc steps the minimum point to find a higher likelihood 
					mcmc_step,start_point,likelihood_min,LOG_LIKELIHOOD_FUNC,PRIOR_FUNC,samp,tsamp,samp_likelihood,error,bestpars,tofit,/ptr
					if error eq 1 then begin
						samp_likelihood=likelihood_min-1.
						print,'mcmc failed'
					endif else print,'mcmc succeeded'
					count=0L
				endelse
				if samp_likelihood gt likelihood_min then break else miss=miss+1
			endwhile		

			;;Determine which ellipsoids does the sample lie in. Ellipsoids could intersect so even though
			;;we sampled from one ellipsoid it could lie in more than one
			samp_num=in_ellipsoids(samp,result,expand=EXPAND)
			r=!RNG->GetRandomNumbers(1,/double)

			;;Accept point with probability 1.0/(number of ellipsoids the point lies in)
			if r le (1.0d/samp_num) then begin
				hit+=1.0

				;;Check if we need more space                            	    				
				if (k+cnt) ge n_samp then begin
					 print,'expand'
					 samples=[samples,ptrarr(NUM)]
					 log_prob=[log_prob,dblarr(NUM,/nozero)]
                                         samples_model=[samples_model,ptrarr(NUM)]
                                         samples_abundances=[samples_abundances,ptrarr(NUM)]
                                         samples_temperature=[samples_temperature,ptrarr(NUM)]
                                         n_samp+=NUM                                         
				endif
				
				;;Determine newest likelihood min and max
				new_likelihood_min=MIN([live_likelihoods,samp_likelihood])
	
				;;Update Log Evidence(logZ) and Information(H) 
				;;Formula in D.S. Sivia Data Analysis: A Bayesian Tutorial
;				print,'Updating Evidence'
				logwt=double(logw + log_plus(new_likelihood_min,likelihood_min) - alog(2.0d));;Trapazoid Integration
				logZnew=log_plus(logZ,logwt) 
				H = exp(logwt-logZnew)*likelihood_min + exp(logZ-logZnew)*(H+logZ)-logZnew
				logZ = logZnew
				logw = logw-1.0d/(double(NUM))

				;;add exiting point to samples and record prob 
				log_prob[k+cnt]=logwt
				samples[k+cnt]=live_samples[likelihood_min_loc]
                                samples_model[k+cnt]=live_model[likelihood_min_loc]
                                samples_abundances[k+cnt]=live_abundances[likelihood_min_loc]
                                samples_temperature[k+cnt]=live_temperature[likelihood_min_loc]

                                ;;Update/Add point to live points
				live_points[likelihood_min_loc]=samp[0]
				live_samples[likelihood_min_loc]=tsamp[0]
				live_likelihoods[likelihood_min_loc]=samp_likelihood
                                live_model[likelihood_min_loc]=modelsamp             ;
                                live_abundances[likelihood_min_loc]=abundancessamp ;
                                live_temperature[likelihood_min_loc]=temperaturesamp ;

                                likelihood_min=MIN(live_likelihoods,likelihood_min_loc,MAX=likelihood_max)

                                ;save,filename='tmp.sav',live_points,live_samples,live_likelihoods,live_model,live_abundances,live_temperature,samples, samples_model,samples_abundances,samples_temperature,log_prob,cnt,k,samp,tsamp,temperaturesamp

                                cnt=cnt+1
			endif else miss+=1
		endwhile
		k=k+cnt
		;;Determine the largest contribution to the log evidence from the existing live points in the sampling loop
		;;If this is less then some tolerance the sampling loop is terminated
		deltaZ = log_plus(logZ,likelihood_max-(k+cnt)/double(NUM))-logZ

		if k mod SAMPLE_NUM eq 0 and not keyword_set(SILENT) then begin
                   print,'DELTA LOGZ: '+strtrim(string(deltaZ))+string(hit/(hit+miss))
                   ;goto,skipsthis
                   ;=====save temporary data for monitoring purposes (DKS edit)
                   logw_temp=-double(k)/double(NUM)-alog(double(NUM))
                   log_prob_temp=log_prob
                   samples_temp=samples
                   samples_model_temp=samples_model
                   samples_abundances_temp=samples_abundances
                   samples_temperature_temp=samples_temperature
                   k_temp=k
                   for i=0, NUM-1 do begin
                      logwt_temp=double(logw_temp+live_likelihoods[i])		
                      logZnew_temp=log_plus(logZ,logwt_temp)
                      H_temp = exp(logwt_temp-logZnew_temp)*live_likelihoods[i] + exp(logZ-logZnew_temp)*(H+logZ)-logZnew_temp
                      logZ_temp=logZnew_temp
                      if (k_temp+cnt) ge n_samp then begin
                         samples_temp=[samples_temp,ptrarr(NUM)]
                         samples_model_temp=[samples_model_temp,ptrarr(NUM)]
                         samples_abundances_temp=[samples_abundances_temp,ptrarr(NUM)]
                         samples_temperature_temp=[samples_temperature_temp,ptrarr(NUM)]
                         log_prob_temp=[log_prob_temp,dblarr(NUM,/nozero)]
                         ;n_samp+=NUM
                      endif
                      log_prob_temp[k_temp]=logwt_temp
                      samples_temp[k_temp]=live_samples[i]
                      samples_model_temp[k_temp]=live_model[i]
                      samples_abundances_temp[k_temp]=live_abundances[i]
                      samples_temperature_temp[k_temp]=live_temperature[i]
                      k_temp+=1
                   endfor
                   samples_temp            =samples_temp[0:k_temp-1]
                   samples_model_temp      =samples_model_temp[0:k_temp-1]
                   samples_abundances_temp =samples_abundances_temp[0:k_temp-1]
                   samples_temperature_temp=samples_temperature_temp[0:k_temp-1]
                   log_prob_temp=log_prob_temp[0:k_temp-1]
                   ;;Normalize log probability
                   log_prob_temp=log_prob_temp-logZ_temp
                   ;;Sort the points for convienence 
                   sort_order=SORT(log_prob_temp)
                   samples_temp            =samples_temp[sort_order]
                   samples_model_temp      =samples_model_temp[sort_order]
                   samples_abundances_temp =samples_abundances_temp[sort_order]
                   samples_temperature_temp=samples_temperature_temp[sort_order]
                   log_prob_temp=log_prob_temp[sort_order]
                   ;Make regular arrays, feed in pointer
                   samps_temp            =dblarr(NUM_PARAMS,n_elements(log_prob_temp))
                   samps_model_temp      =dblarr(NUM_WS,n_elements(log_prob_temp))
                   samps_abundances_temp =dblarr(NUM_PRESS,NUM_SPECIES,n_elements(log_prob_temp))
                   samps_temperature_temp=dblarr(NUM_PRESS,n_elements(log_prob_temp))
                   for i=0,n_elements(log_prob_temp)-1 do begin & $
                      samps_temp[*,i]              = (*(samples_temp[i])) & $
                      samps_model_temp[*,i]        = (*(samples_model_temp[i])) & $
                      samps_abundances_temp[*,*,i] = (*(samples_abundances_temp[i])) & $
                      samps_temperature_temp[*,i]  = (*(samples_temperature_temp[i])) & $
                   endfor	
                   ;print,'H = ' , H 
                   ;print,'Log(Z) = ' + string(logZ) + ' +- '+ string(sqrt(H/(1.0*NUM)))
                   data_str_temp={samples:samps_temp,log_prob:log_prob_temp,logZ:logZ_temp,logZ_err:sqrt(H_temp/(1.0*NUM)),H:H_temp,samples_model:samps_model_temp,samples_abundances:samps_abundances_temp,samples_temperature:samps_temperature_temp,ws:ws,pressure:pressure,x:x,y:y,err:err,waverr:waverr,molname:molname}
                   openw,1,dirs+'multinest.monitor.log.txt',/append
                     printf,1,string(systime(/Julian,/utc),format='(f19.7)')+' '+strtrim(string(deltaZ))+' '+strtrim(string(logZ_temp))
                   close,1
                   ;save temporary monitoring data
                   ;LPtemp=(log_prob_temp-alog(total(exp(log_prob_temp))))  ;normalize so total probability=1
                   save,filename=dirs+'multinest.monitor.sav',data_str_temp,TOL,deltaZ,threadtrack,counttrack;,LPtemp
                   ;=======end save temporary data
                   ;skipsthis:
                endif
                             
                for i=0,clusters-1 do ptr_free,result[i].data_ptr
                ;save current iteration to be able to re-start Nested Run
                save,filename=dirs+'Multinest_current_samples.sav',live_samples,live_likelihoods,live_points,live_model,live_abundances,live_temperature,ws,pressure,molname,H,logZ,k,logw,deltaZ,samples,log_prob,n_samp,samples_model,samples_abundances,samples_temperature,likelihood_min,NUM_PARAMS,NUM_WS,NUM_PRESS,NUM_SPECIES,threadtrack,counttrack
                save,filename=dirs+'Multinest_current_samples_threads.sav',live_samples,live_likelihoods,live_points,live_model,live_abundances,live_temperature
        endwhile
        ;===============================================================================================

	
	;;Add the log evidence contribution from the remaining live points 
	logw=-double(k)/double(NUM)-alog(double(NUM))

	for i=0, NUM-1 do begin
		logwt=double(logw+live_likelihoods[i])		
		logZnew=log_plus(logZ,logwt)
		H = exp(logwt-logZnew)*live_likelihoods[i] + exp(logZ-logZnew)*(H+logZ)-logZnew
		logZ=logZnew
                if (k+cnt) ge n_samp then begin
                   samples=[samples,ptrarr(NUM)]
                   samples_model=[samples_model,ptrarr(NUM)]
                   samples_abundances=[samples_abundances,ptrarr(NUM)]
                   samples_temperature=[samples_temperature,ptrarr(NUM)]
		   log_prob=[log_prob,dblarr(NUM,/nozero)]
                   n_samp+=NUM
                endif
                log_prob[k]=logwt
                samples[k]=live_samples[i]
                samples_model[k]=live_model[i]
                samples_abundances[k]=live_abundances[i]
                samples_temperature[k]=live_temperature[i]
                k+=1
	endfor
        samples            =samples[0:k-1]
        samples_model      =samples_model[0:k-1]
        samples_abundances =samples_abundances[0:k-1]
        samples_temperature=samples_temperature[0:k-1]
	log_prob=log_prob[0:k-1]
	;;Normalize log probability
	log_prob=log_prob-logZ
	;;Sort the points for convienence 
	sort_order=SORT(log_prob)
        samples            =samples[sort_order]
        samples_model      =samples_model[sort_order]
        samples_abundances =samples_abundances[sort_order]
        samples_temperature=samples_temperature[sort_order]
	log_prob=log_prob[sort_order]
        ;Make regular arrays, feed in pointer
	samps            =dblarr(NUM_PARAMS,n_elements(log_prob))
        samps_model      =dblarr(NUM_WS,n_elements(log_prob))
        samps_abundances =dblarr(NUM_PRESS,NUM_SPECIES,n_elements(log_prob))
        samps_temperature=dblarr(NUM_PRESS,n_elements(log_prob))
	for i=0,n_elements(log_prob)-1 do begin
          samps[*,i]              = (*(samples[i]))
          samps_model[*,i]        = (*(samples_model[i]))
          samps_abundances[*,*,i] = (*(samples_abundances[i]))
          samps_temperature[*,i]  = (*(samples_temperature[i]))
	endfor	
	print,'H = ' , H 
	print,'Log(Z) = ' + string(logZ) + ' +- '+ string(sqrt(H/(1.0*NUM)))
	
	;;Return data structure that contains the samples,the log probability,
	;;the log evidence and error, and the information respectively
	data_str={samples:samps,log_prob:log_prob,logZ:logZ,logZ_err:sqrt(H/(1.0*NUM)),H:H,samples_model:samps_model,samples_abundances:samps_abundances,samples_temperature:samps_temperature,ws:ws,pressure:pressure,x:x,y:y,err:err,waverr:waverr,molname:molname}
	if keyword_set(USE_THREADS) then destroy_threads,threads
	heap_free,samples
	heap_free,live_points
	heap_free,live_samples
	return, data_str
	GET_OUT:
END
