;read in transmission spectrum 
filename='transmision_spectrum.ncdf'

fid=ncdf_open(filename)
ncdf_varget,fid,'nu',nu
ncdf_varget,fid,'transit_radius',rprs
ncdf_close,fid

print,nu
print,rprs

nu = 1./nu
nu = nu/100.
nu = nu*1d6

;Convert to 'transit depth'
;rprs = rprs^2.

;;plot transmission spectrum
;set_plot,'ps'
;!p.font=0
;device,/encapsul,filename='hd209_trans_spec.eps',/helvetica
;plot,nu,rprs,/xlog,thick=2.0,xran=[0.3,30],/xsty,/ysty,/nodata,$
;  ytitle='Rp/Rs',xtitle='Wavelength (microns)',xthi=4.0,ythi=4.0,charth=4.0,charsi=1.4,$
;  title='HAT-P-1b'
;oplot,nu,rprs,thick=2.0
;
;device,/close

;now read in emission spectrum
filename2='spectrum.ncdf'

fid=ncdf_open(filename2)
ncdf_varget,fid,'nu',nu
ncdf_varget,fid,'fnu',fnup
ncdf_varget,fid,'fnu_star',fnus
ncdf_varget,fid,'R_planet_TOA',rp
ncdf_varget,fid,'R_star',rs
ncdf_close,fid

;Convert to wavelength and SI
nu = 1.d0/nu
nu = nu/100.d0

;Convert to um
nu = nu*1.d6

;Flux ratio
fpfs = fnup/fnus
fpfs = fpfs*(rp^2./rs^2.)

;plot transmission spectrum
set_plot,'ps'
!p.font=0
device,/encapsul,filename='hd209_trans_spec.eps',/helvetica
plot,nu,rprs,/xlog,thick=2.0,xran=[0.3,30],/xsty,/ysty,/nodata,$
  ytitle='Rp/Rs',xtitle='Wavelength (microns)',xthi=4.0,ythi=4.0,charth=4.0,charsi=1.4,$
  title='HD209458b'
oplot,nu,rprs,thick=2.0

device,/close

;plot emission spectrum
set_plot,'ps'
!p.font=0
device,/encapsul,filename='hd209_emspec.eps',/helvetica
plot,nu,fpfs,/xlog,thick=2.0,xran=[0.3,30],/xsty,/ysty,/nodata,$
  ytitle='Fp/Fs',xtitle='Wavelength (microns)',xthi=4.0,ythi=4.0,charth=4.0,charsi=1.4,$
  title='HD209458b'
oplot,nu,fpfs,thick=2.0

device,/close

set_plot,'X'
!p.font=-1

end
