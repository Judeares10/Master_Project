fid=ncdf_open('pt_hat1_v3.ncdf')
ncdf_varget,fid,'pressure_si',pressure_si
ncdf_varget,fid,'temperature',temperature
ncdf_close,fid

set_plot,'ps'
!p.font=0
device,filename='hat1_atmo_pt.eps',/encapsulated,/helvetica
plot,temperature,pressure_si*1e-5,/ylog,yrange=[100,1e-5],thick=2.0,/ysty,/nodata,$
  ytitle='Pressure (bars)',xtitle='Temperature (K)',xthi=4.0,ythi=4.0,charth=4.0,charsi=1.4,$
  title='HAT-P-1b'
oplot,temperature,pressure_si*1e-5,thick=2.0

device,/close
set_plot,'X'
!p.font=-1

end

