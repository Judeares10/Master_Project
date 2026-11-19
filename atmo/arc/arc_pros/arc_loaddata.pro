pro arc_loaddata,datafile,wsdstart,wsd,r,rerr,werr   

   ;load in data
   data=read_ascii(datafile)  
   data=data.(0)
   wsd=transpose(data[0,*])
   werr=transpose(data[1,*])
   r=transpose(data[2,*])
   rerr=transpose(data[3,*])

   ;sort data by wavelength
   srt=sort(wsd)
   wsd=wsd[srt]
   werr=werr[srt]
   r=r[srt]
   rerr=rerr[srt]
   ;wavelengnth cutoff if requested
   if keyword_set(wsdstart) then begin & a=where(wsd gt wsdstart) & wsd=wsd(a) & r=r(a) & rerr=rerr(a) & werr=werr(a) & endif

end
