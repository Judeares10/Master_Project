PRO mcmc_step,start,constraint,log_likelihood_func,prior_func,samp,tsamp,log_likelihood,error,bestpars,tofit,NUM_STEPS=NUM_STEPS,ptr=ptr

;	DefSysV, '!RNG', Obj_New('RandomNumberGenerator')
	;;SET NUMBER OF MCMC STEPS
	if not keyword_set(NUM_STEPS) then NUM_STEPS=20
	COMMON mcmc_common  ;accept,tot,sigma

	tstart = CALL_FUNCTION(prior_func,start,tofit)
	p = CALL_FUNCTION(LOG_LIKELIHOOD_FUNC,tstart,bestpars,tofit)  ;CALL_FUNCTION(log_likelihood_func,tstart)
	xx = *start[0]
	r = !RNG -> GetRandomNumbers(n_elements(xx),NUM_STEPS,/double,/NORMAL)
	i=0
	xt_prime=tstart
	while i lt NUM_STEPS do begin
;		print,'step',i
		tot+=1.0d
		x_prime = xx + sigma*r[*,i]
;		print,x_prime,xx
		w=where(x_prime gt 1.0D or x_prime lt 0.0D,nw)

		if nw ne 0 then begin
			i=i+1
			continue
		endif
		xt_prime=CALL_FUNCTION(prior_func,ptr_new(x_prime),tofit)
		p_prime=CALL_FUNCTION(log_likelihood_func,xt_prime,bestpars,tofit)
		a = p_prime/p
		if (r[i] le a) and (p_prime gt constraint) then begin
			xx = x_prime
			p = p_prime
			accept+=1.0d
		endif else tot+=1.0d
		i = i+1 
	endwhile	

	if accept/tot gt .9 then sigma = sigma * exp(1.0d / accept)
	if accept/tot lt .6 then sigma = sigma * exp(-1.0d / (tot-accept))

	samp=xx
	error=0
	if array_equal(samp,*start[0]) eq 1 then error=1
	if keyword_set(ptr) then begin
		samp=ptr_new(samp)
		tsamp=xt_prime
	endif else tsamp= *xt_prime
	log_likelihood=p
END 
