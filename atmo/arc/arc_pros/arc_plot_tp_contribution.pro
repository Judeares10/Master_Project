pro arc_plot_tp_contribution,dirs,tpfile,contribfile
;tpfile='pt_ARCnlevel36.ncdf'
;contribfile='ARC.emission.contribution.highrest.nlevel36.ncdf'

;Open ATMO NCDF file of P-T profile
  fid=ncdf_open(dirs+tpfile)
  ncdf_varget,fid,'pressure',pressure ;cgs Barye, bar is pressure*1D-6
  ncdf_varget,fid,'temperature',temperature ;cgs Barye, bar is pressure*1D-6
  p=alog10(pressure*1D-6)
  ncdf_close,fid


  fid=ncdf_open(dirs+contribfile)
  ncdf_varget,fid,'nu',nu ;cgs Barye, bar is pressure*1D-6
  ncdf_varget,fid,'pressure_si',pressure_si ;cgs Barye, bar is pressure*1D-6
  ncdf_varget,fid,'pressure_interface',pressure_interface    
  ncdf_varget,fid,'ncfnu',ncfnu ;cgs Barye, bar is pressure*1D-6
  pint=alog10(pressure_interface*1D-6)
  p_si=alog10(pressure_si*1D-6)
  ncdf_close,fid
  ws=10000./nu

  a=where(ws gt 1.4 and ws lt 1.44)
  a=where(ws gt 0.8 and ws lt 0.9)
  a=where(ws gt 1.6 and ws lt 1.7)

  ymin=2;max(p)
  ymax=-8;min(p)
  maxp=where(p ge ymin) & minp=where(p ge ymax)
  plot,temperature,p,psym=-1,yrange=[ymin,ymax],xrange=[1500,3500];temperature(maxp(0)),temperature(minp(0))]
  rescale=max(+temperature(minp(0)))-min(temperature(maxp(0)))
  for i=0,n_elements(a)-1,2 do oplot,ncfnu(*,a(i))*rescale+min(temperature),p_si,color=1000*i+1000
    for i=0,n_elements(a)-1,2 do oplot,ncfnu(*,a(i))*rescale+min(temperature),pint,color=1000*i+2000




end
