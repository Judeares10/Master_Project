FUNCTION arc_box_prior,params
; turn random numbers generated between [0,1] to 
; pramameter values spanning the bounding box of the fit space
        ;params pointer array
	n_samps=n_elements(params)
	output=ptrarr(n_samps)
	for i=0,n_samps-1 do begin
		tmp=*params[i]
                print,'tmp:',tmp
	    tmp[0] = -1.0d + 2.0d * tmp[0]
	    a= -( 1.0d - (tmp[0])^2.0d )^.5d
	    b=  ( 1.0d - (tmp[0])^2.0d )^.5d
	    tmp[1]= a[0] + (b[0]-a[0])*tmp[1]
            print,'out:',tmp
		output[i]=ptr_new(tmp)
	endfor
    return,output
                                                	
                                                	                                                	
END
