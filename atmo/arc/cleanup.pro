pro cleanup
; removes all ATMO files generated from parallel processing
; does not touch ARC output files in the save directory dirs
; Also removes anoying Mac Droppings .DS_Store .Apple_Double
;
; If connected to directory with OSX, AppleDouble files will re-appear 
;
spawn,'rm ARC*.in'
spawn,'rm ARC*spec*.ncdf'
spawn,'rm chem_*eq*.ncdf'
spawn,'rm pt_ARC*.ncdf'
spawn,'rm ARC.emission.contribution*.ncdf'
;
spawn,"find . -name '*.DS_Store' -type f -delete"
spawn,"find . -name '*.Trashes' -type d -delete"

spawn,"find . -name '*.AppleDouble' -type d "  ;lists all anyoing double files
spawn,'rm -rR arc_multinest_pros/backup/.AppleDouble/'
spawn,'rm -rR arc_multinest_pros/.AppleDouble/'
spawn,'rm -rR stellar_spectra/.AppleDouble/'
spawn,'rm -rR coyote/public/.AppleDouble/'
spawn,'rm -rR coyote/.AppleDouble/'
spawn,'rm -rR arc_pros/backup/.AppleDouble/'
spawn,'rm -rR arc_pros/.AppleDouble/'
spawn,'rm -rR .AppleDouble/'
spawn,"find . -name '*.AppleDouble' -type d "

end
