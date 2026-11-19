FUNCTION example_prior_work,params
        ;params is a pointer array
        xmin=-2d0
        xmax=1d0
        ymin=-1d0
        ymax=2d0

        x0=(xmax+xmin)/2.
        xwidth=(xmax-xmin)/2.
        y0=(ymin+ymax)/2.
        ywidth=(ymax-ymin)/2.
	n_samps=n_elements(params)
	output=ptrarr(n_samps)
	for i=0,n_samps-1 do begin
		tmp=*params[i]
            ;transform random numbers to correspond to a uniform prior between xmin and xmax
	    tmp[0] = xmin + (xmax-xmin) * tmp[0]
            ;impose all the random numbers be contained within an ellipse
            ;  ellipse centered at x0,y0 with width defined by prior bounds
	    a= -( xwidth^2. - (tmp[0]-x0)^2.0d )^.5d
	    b=  ( xwidth^2 - (tmp[0]-x0)^2.0d )^.5d
	    tmp[1]=(a[0] + (b[0]-a[0])*tmp[1])*ywidth/xwidth +y0
            output[i]=ptr_new(tmp)
	endfor
    return,output
                                                	
END
