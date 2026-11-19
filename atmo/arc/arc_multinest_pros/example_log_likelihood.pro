FUNCTION example_log_likelihood,params

	ptmp=params[0]
    x=[-.35d0,-.4d0,-.2d0,.45d0,.1d0]
    y=[.2d0,-.4d0,.15d0,.1d0,-.15d0]
    sig=[.01d0,.01d0,.03d0,.05d0,.02d0]
    A=[1.0d0,.5d0,.8d0,.6d0,.5d0]
    d=A*exp(-((x- (*ptmp)[0])^2.0d + (y- (*ptmp)[1])^2.0d)/(2.0d * sig^2.0d))

    w=where(d eq 0.0d or d eq !values.f_nan,nw)  ;find where d is zero or NaN
    if nw ne 0 then d[w]=(MACHAR(/double)).xmin  ; replace zeros or NaN with minumum usable value
    ;help,(MACHAR(/double)).xmin 
    ;<Expression>    DOUBLE    =   2.2250739e-308
    ;alog((MACHAR(/double)).xmin) =       -708.39642

;print,d
    tmp=total(d)
;print,'alog(tmp)=',alog(tmp)
;print,''
	return,alog(tmp)
END
                                    
