FUNCTION arc_multinest_prior,params,tofit
  common passmodel,ws,model,modelpoints,LowPriorsParams,UpPriorsParams,abundances,pressure,temperature

  ;reads in priors, calculate width & center of boundries
  Nparams=n_elements(tofit)
  MinD=dblarr(Nparams)
  MaxD=dblarr(Nparams)
  WidthD=dblarr(Nparams)
  CenD=dblarr(Nparams)
  for i=0,Nparams-1 do begin
    MinD(i)=LowPriorsParams(tofit(i))
    MaxD(i)=UpPriorsParams(tofit(i))
    CenD(i)=(LowPriorsParams(tofit(i))+UpPriorsParams(tofit(i)))/2.
    WidthD(i)=(UpPriorsParams(tofit(i))-LowPriorsParams(tofit(i)))/2.
   endfor

	n_samps=n_elements(params)
	output=ptrarr(n_samps)
	for i=0,n_samps-1 do begin
		tmp=*params[i]
            ;transform random numbers to correspond to a uniform prior between xmin and xmax
	    tmp[0] = minD(0) + (maxD(0)-minD(0)) * tmp[0]
            ;impose all the random numbers be contained within an n-dimentional ellipse
            ;  ellipse centered at x0,y0 with width defined by prior bounds
	    a= -( WidthD(0)^2. - (tmp[0]-CenD(0))^2.0d )^.5d
	    b=  ( widthD(0)^2. - (tmp[0]-CenD(0))^2.0d )^.5d

	    for j=1,Nparams-1 do begin
              tmp[j]=(a[0] + (b[0]-a[0])*tmp[j])*WidthD(j)/WidthD(0) + CenD(j)
	      ;tmp[2]=(a[0] + (b[0]-a[0])*tmp[2])*zwidth/xwidth +z0
	      ;tmp[3]=(a[0] + (b[0]-a[0])*tmp[3])*cwidth/xwidth +c0
            endfor
            output[i]=ptr_new(tmp)
	endfor
    return,output
                                                	
END
