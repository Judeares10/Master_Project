pro arc_load_atmo_directimaging_spec,atmofilenam,Rp,Distance,ws,modelfpfs
  ;==LOAD ATMO MODEL=======
  fid=ncdf_open(atmofilenam)
  ncdf_varget,fid,'nu',nu             ;cm^-1
  ncdf_varget,fid,'fnu',fnu           ;ers s-1 cm-2 cm
 ; ncdf_varget,fid,'R_planet_TOA',Rpl ;cm  Top of atmosphere (BAD)
 ; ncdf_varget,fid,'R_star',Rs        ;cm
 ; ncdf_varget,fid,'fnu_star',fstar   ;ers s-1 cm-2 cm
  ncdf_close,fid
  ws=10000./nu
  ws_cm=1./nu  ;wavelength in cm
  flamb= fnu/ws_cm^2. * 1D-7          ; flambda in W/m^2/um (converting cm**2 -> m cancels factor of cm->um)
  modelfpfs=flamb*(Rp/Distance)^2.    
end
