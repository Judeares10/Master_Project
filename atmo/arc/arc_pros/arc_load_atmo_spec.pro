pro arc_load_atmo_spec,atmofilenam,ws,model
  ;==LOAD ATMO MODEL=======
  fid=ncdf_open(atmofilenam)
  ncdf_varget,fid,'nu',nu
  ncdf_varget,fid,'transit_radius',rprs
  ncdf_close,fid
  ws=10000./nu
  Model=rprs


end
