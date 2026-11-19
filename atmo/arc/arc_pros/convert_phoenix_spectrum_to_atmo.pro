pro convert_phoenix_spectrum_to_atmo,filename
; read in stellar model from NextGen database
; Download model from here:
; http://svo2.cab.inta-csic.es/theory/newov/index.php?model=bt-nextgen-agss2009
; filename='HatP11A.Phoenix'

readcol,filename+'.txt',nu2,hnu2,format='D,D'
arr_ind=n_elements(hnu2)

;define arrays
nu_convert=dblarr(arr_ind)
hnu_convert=dblarr(arr_ind)

;convert nu and hnu to wavenumber units
for i=0,arr_ind-1 do begin & $
    nu_convert(i) = 1./(1.d-8*nu2(i)) & $
    hnu_convert(i) = ((hnu2(i)*1.d8)/(4.*!pi))/(nu_convert(i)^2) & $
endfor

;reverse order of arrays, required by ATMO
nu_convert=reverse(nu_convert)
hnu_convert=reverse(hnu_convert)

;Generate nu array with 500000 points and interpolate fluxes
;This is required for ATMO input
nu_final=dindgen(500000L,increment=0.1d0,start=0.05d0)
hnu_final=interpol(hnu_convert,nu_convert,nu_final)

;Create and write NetCDF file
id = NCDF_CREATE(filename+'.ncdf', /CLOBBER) ; Create the NetCDF file.
xid = NCDF_DIMDEF(id, 'nfreq',  n_elements(nu_final)) ; Define the dimension for the two variables.
vid0 = NCDF_VARDEF(id, 'nu', [xid], /DOUBLE) ;define variable nu
vid1 = NCDF_VARDEF(id, 'hnu', [xid], /DOUBLE) ;define variable hnu
NCDF_CONTROL, id, /ENDEF ; Put the file into data mode.
NCDF_VARPUT, id, vid0, nu_final  ; write converted nu to file
NCDF_VARPUT, id, vid1, hnu_final ; write converted hnu to file

;if you like, check the structure of the file to ensure it is correct.
;file_inq = NCDF_INQUIRE(id)
;HELP, file_inq, /STRUCTURE
NCDF_CLOSE, id ; Close the NetCDF file.


end
