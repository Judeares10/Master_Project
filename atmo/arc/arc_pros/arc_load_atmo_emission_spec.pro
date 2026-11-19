pro arc_load_atmo_emission_spec,atmofilenam,ws,rplanet,modelfpfs,fnu,fstar
  ;==LOAD ATMO MODEL=======
  fid=ncdf_open(atmofilenam)
  ncdf_varget,fid,'nu',nu            ;cm^-1
  ncdf_varget,fid,'fnu',fnu          ;ers s-1 cm-2 cm
 ;ncdf_varget,fid,'R_planet_TOA',Rpl ;cm  Top of atmosphere (BAD)
  ncdf_varget,fid,'R_star',Rs        ;cm
  ncdf_varget,fid,'fnu_star',fstar   ;ers s-1 cm-2 cm
  ncdf_close,fid
  ws=10000./nu
  modelfpfs=fnu/fstar*(rplanet/rs)^2.
end
