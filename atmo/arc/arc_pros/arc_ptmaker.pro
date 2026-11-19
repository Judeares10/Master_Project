pro arc_ptmaker,nlevel,pmin,pmax
; Used to create the needed ncdf files with a user defined pressure
; array for ATMO.  Chem_dummy files is also crates (with rainout).
; IDL wrapper writes /pacode/atmo/visu/arc_ptmaker.py which uses
; /pacode/atmo/visu/write_atmo.py 
; to generate the needed ncdf files, saving the ncdf files to the current
; working directory
; 
; Inputs:  nlevel = number of levels (2 digits)
;          pmin   = minimum pressure (bars)
;          pmax   = maximum pressure (bars)
;
; Outputs: input_TP_arc.ncdf   = T-P profile ncdf file
;          chem_dummy_arc.ncdf = Input chemistry ncdf file
;
; Requires Python.
; 
; DKS 20/9/2016

CD, Current=theDirectory ;current working directory

close,1 & openw,1,'../../visu/arc_ptmaker.py'

  printf,1,'import numpy as np'
  printf,1,'import write_atmo as wa'
  printf,1,'nlevel='+string(nlevel,format='(I2)')  ; #levels in the grid'
  printf,1,'pmin='+string(pmin,format='(E8.1)')    ;1E-5 #bar'
  printf,1,'pmax='+string(pmax,format='(E8.1)')   ;1E3  #bar'

  printf,1,'pp=np.zeros(nlevel)'
  printf,1,'pp[0]      = pmin*1E6 #bar -> cgs'
  printf,1,'pp[nlevel-1] = pmax*1E6 #bar -> cgs'
  printf,1,'dlp = (np.log10(pmax)-np.log10(pmin))/(nlevel-1)'
  printf,1,'for i in range(1,nlevel,1):'
  printf,1,'  pp[i] = pp[i-1]*10.**dlp'
 
  printf,1,'tt = np.arange(1001,1001+nlevel,1) #dummy temperatures'
  printf,1,"wa.write_atmo(pp,tt,kzz='none',fout='"+theDirectory+"/input_TP_arc.ncdf')"
  printf,1,"wa.write_abundances(fin='../../visu/inidat_dummy_cond',fout='"+theDirectory+"/chem_dummy_arc.ncdf',nlevel=nlevel,pmin=pmin,pmax=pmax)"
  close,1
  spawn,'python ../../visu/arc_ptmaker.py'
  print,'== Pressure Array==='
  spawn,'ncdump input_TP_arc.ncdf'
  print,'===================='
end

