pro arc_readchem,dirs,filenam,abundances,pressure,temperature

Fid=ncdf_open(dirs+filenam)  ;test if file is ARC output, if not then load as is
ncdf_varget,fid,'abundances',abundances
ncdf_varget,fid,'pressure',pressure
ncdf_varget,fid,'temperature',temperature
ncdf_varget,fid,'molname',molname
ncdf_close,fid

end

