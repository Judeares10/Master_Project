pro arc_make_atmo_pieces,dirs,filenam,teff,MdH,logg,Rp,Cloud,haze,cloudtop,art_haze_alpha,kzz,fsed,sigma_g,fitabundances
  Common SystemParameters,rplanet,rstar,rorbit,stellarfile,Tstar,Cloudtype,Hazetype,Chemistry,Isothermal,Opacity_switch,Abund_fit_switch,Abund_fit_index,Abund_priorLow,Abund_priorUP,rain_switch,Element_fit_switch,Element_fit_index,Ele_priorLow,Ele_priorUP,Distance,nband
    Opacity_switch_original=Opacity_switch
    Cloudtype_original=Cloudtype
    Hazetype_original=Hazetype
    ChildProcessCount=0
;==================
; Purpose:
;  Makes spectra with individual opacity contribution


   ;Cloud
   if (Cloudtype_original eq 1) then begin
     Opacity_names =['H2','He','H2O','CO2','CO','CH4','NH3','Na','K','Li','Rb','Cs','TiO','VO','FeH','PH3','H2S','HCN','C2H2','SO2','Fe','H_']
     Opacity_switch=[   1,   1,    0,    0,   0,    0,    0,   0,  0,   0,   0,   0,    0,   0,   0,    0,     0,    0,     0,    0,   0,   0] ;Hat26
     Cloudtype=1
     Hazetype=0
     spawn,'rm ARC.transitspec.ncdf'
     arc_make_atmo_in,teff,MdH,logg,Rp,Cloud,haze,cloudtop,art_haze_alpha,kzz,fsed,sigma_g,fitabundances
     spawn,'cp ARC.transitspec99.ncdf '+dirs+'ARC.transitspec.Cloud.ncdf'  
   end   

   ;Haze
   if (Hazetype_original eq 1) then begin
     Opacity_names =['H2','He','H2O','CO2','CO','CH4','NH3','Na','K','Li','Rb','Cs','TiO','VO','FeH','PH3','H2S','HCN','C2H2','SO2','Fe','H_']
     Opacity_switch=[   1,   1,    0,    0,   0,    0,    0,   0,  0,   0,   0,   0,    0,   0,   0,    0,     0,    0,     0,    0,   0,   0] ;Hat26
     Cloudtype=0
     Hazetype=1
     spawn,'rm ARC.transitspec.ncdf'
     arc_make_atmo_in,teff,MdH,logg,Rp,Cloud,haze,cloudtop,art_haze_alpha,kzz,fsed,sigma_g,fitabundances
     spawn,'cp ARC.transitspec99.ncdf '+dirs+'ARC.transitspec.Haze.ncdf'  
   end   
   ;turn off clouds/haze for rest of molecules
   Cloudtype=0
   Hazetype=0

   ;H2
   if (Opacity_switch_original(0) eq 1) then begin
     Opacity_names =['H2','He','H2O','CO2','CO','CH4','NH3','Na','K','Li','Rb','Cs','TiO','VO','FeH','PH3','H2S','HCN','C2H2','SO2','Fe','H_']
     Opacity_switch=[   1,   0,    0,    0,   0,    0,    0,   0,  0,   0,   0,   0,    0,   0,   0,    0,     0,    0,     0,    0,   0,   0] ;Hat26
     spawn,'rm ARC.transitspec99.ncdf'
     arc_make_atmo_in,teff,MdH,logg,Rp,Cloud,haze,cloudtop,art_haze_alpha,kzz,fsed,sigma_g,fitabundances
     spawn,'cp ARC.transitspec99.ncdf '+dirs+'ARC.transitspec.H2.ncdf'  
   end

   ;He
   if (Opacity_switch_original(1) eq 1) then begin
     Opacity_names =['H2','He','H2O','CO2','CO','CH4','NH3','Na','K','Li','Rb','Cs','TiO','VO','FeH','PH3','H2S','HCN','C2H2','SO2','Fe','H_']
     Opacity_switch=[   0,   1,    0,    0,   0,    0,    0,   0,  0,   0,   0,   0,    0,   0,   0,    0,     0,    0,     0,    0,   0,   0] ;Hat26
     spawn,'rm ARC.transitspec99.ncdf'
     arc_make_atmo_in_molecules,teff,MdH,logg,Rp,Cloud,haze,cloudtop,art_haze_alpha,kzz,fsed,sigma_g,fitabundances
     spawn,'cp ARC.transitspec99.ncdf '+dirs+'ARC.transitspec.He.ncdf'  
   endif

   ;H2O
   if (Opacity_switch_original(2) eq 1) then begin
     Opacity_names =['H2','He','H2O','CO2','CO','CH4','NH3','Na','K','Li','Rb','Cs','TiO','VO','FeH','PH3','H2S','HCN','C2H2','SO2','Fe','H_']
     Opacity_switch=[   0,   0,    1,    0,   0,    0,    0,   0,  0,   0,   0,   0,    0,   0,   0,    0,     0,    0,     0,    0,   0,   0] ;Hat26
     spawn,'rm ARC.transitspec99.ncdf'
     arc_make_atmo_in_molecules,teff,MdH,logg,Rp,Cloud,haze,cloudtop,art_haze_alpha,kzz,fsed,sigma_g,fitabundances
     spawn,'cp ARC.transitspec99.ncdf '+dirs+'ARC.transitspec.H2O.ncdf'  
   endif

   ;CO2
   if (Opacity_switch_original(3) eq 1) then begin
     Opacity_names =['H2','He','H2O','CO2','CO','CH4','NH3','Na','K','Li','Rb','Cs','TiO','VO','FeH','PH3','H2S','HCN','C2H2','SO2','Fe','H_']
     Opacity_switch=[   0,   0,    0,    1,   0,    0,    0,   0,  0,   0,   0,   0,    0,   0,   0,    0,     0,    0,     0,    0,   0,   0] ;Hat26
     spawn,'rm ARC.transitspec99.ncdf'
     arc_make_atmo_in_molecules,teff,MdH,logg,Rp,Cloud,haze,cloudtop,art_haze_alpha,kzz,fsed,sigma_g,fitabundances
     spawn,'cp ARC.transitspec99.ncdf '+dirs+'ARC.transitspec.CO2.ncdf'  
   endif

   ;CO
   if (Opacity_switch_original(4) eq 1) then begin
     Opacity_names =['H2','He','H2O','CO2','CO','CH4','NH3','Na','K','Li','Rb','Cs','TiO','VO','FeH','PH3','H2S','HCN','C2H2','SO2','Fe','H_']
     Opacity_switch=[   0,   0,    0,    0,   1,    0,    0,   0,  0,   0,   0,   0,    0,   0,   0,    0,     0,    0,     0,    0,   0,   0] ;Hat26
     spawn,'rm ARC.transitspec99.ncdf'
     arc_make_atmo_in_molecules,teff,MdH,logg,Rp,Cloud,haze,cloudtop,art_haze_alpha,kzz,fsed,sigma_g,fitabundances
     spawn,'cp ARC.transitspec99.ncdf '+dirs+'ARC.transitspec.CO.ncdf'  
   endif

   ;CH4
   if (Opacity_switch_original(5) eq 1) then begin
     Opacity_names =['H2','He','H2O','CO2','CO','CH4','NH3','Na','K','Li','Rb','Cs','TiO','VO','FeH','PH3','H2S','HCN','C2H2','SO2','Fe','H_']
     Opacity_switch=[   0,   0,    0,    0,   0,    1,    0,   0,  0,   0,   0,   0,    0,   0,   0,    0,     0,    0,     0,    0,   0,   0] ;Hat26
     spawn,'rm ARC.transitspec99.ncdf'
     arc_make_atmo_in_molecules,teff,MdH,logg,Rp,Cloud,haze,cloudtop,art_haze_alpha,kzz,fsed,sigma_g,fitabundances
     spawn,'cp ARC.transitspec99.ncdf '+dirs+'ARC.transitspec.CH4.ncdf'  
   endif

   ;NH3
   if (Opacity_switch_original(6) eq 1) then begin
     Opacity_names =['H2','He','H2O','CO2','CO','CH4','NH3','Na','K','Li','Rb','Cs','TiO','VO','FeH','PH3','H2S','HCN','C2H2','SO2','Fe','H_']
     Opacity_switch=[   0,   0,    0,    0,   0,    0,    1,   0,  0,   0,   0,   0,    0,   0,   0,    0,     0,    0,     0,    0,   0,   0] ;Hat26
     spawn,'rm ARC.transitspec99.ncdf'
     arc_make_atmo_in_molecules,teff,MdH,logg,Rp,Cloud,haze,cloudtop,art_haze_alpha,kzz,fsed,sigma_g,fitabundances
     spawn,'cp ARC.transitspec99.ncdf '+dirs+'ARC.transitspec.NH3.ncdf'  
   endif

   ;Na
   if (Opacity_switch_original(7) eq 1) then begin
     Opacity_names =['H2','He','H2O','CO2','CO','CH4','NH3','Na','K','Li','Rb','Cs','TiO','VO','FeH','PH3','H2S','HCN','C2H2','SO2','Fe','H_']
     Opacity_switch=[   0,   0,    0,    0,   0,    0,    0,   1,  0,   0,   0,   0,    0,   0,   0,    0,     0,    0,     0,    0,   0,   0] ;Hat26
     spawn,'rm ARC.transitspec99.ncdf'
     arc_make_atmo_in_molecules,teff,MdH,logg,Rp,Cloud,haze,cloudtop,art_haze_alpha,kzz,fsed,sigma_g,fitabundances
     spawn,'cp ARC.transitspec99.ncdf '+dirs+'ARC.transitspec.Na.ncdf'  
   endif

   ;K
   if (Opacity_switch_original(8) eq 1) then begin
     Opacity_names =['H2','He','H2O','CO2','CO','CH4','NH3','Na','K','Li','Rb','Cs','TiO','VO','FeH','PH3','H2S','HCN','C2H2','SO2','Fe','H_']
     Opacity_switch=[   0,   0,    0,    0,   0,    0,    0,   0,  1,   0,   0,   0,    0,   0,   0,    0,     0,    0,     0,    0,   0,   0] ;Hat26
     spawn,'rm ARC.transitspec99.ncdf'
     arc_make_atmo_in_molecules,teff,MdH,logg,Rp,Cloud,haze,cloudtop,art_haze_alpha,kzz,fsed,sigma_g,fitabundances
     spawn,'cp ARC.transitspec99.ncdf '+dirs+'ARC.transitspec.K.ncdf'  
   endif

   ;Li
   if (Opacity_switch_original(9) eq 1) then begin
     Opacity_names =['H2','He','H2O','CO2','CO','CH4','NH3','Na','K','Li','Rb','Cs','TiO','VO','FeH','PH3','H2S','HCN','C2H2','SO2','Fe','H_']
     Opacity_switch=[   0,   0,    0,    0,   0,    0,    0,   0,  0,   1,   0,   0,    0,   0,   0,    0,     0,    0,     0,    0,   0,   0] ;Hat26
     spawn,'rm ARC.transitspec99.ncdf'
     arc_make_atmo_in_molecules,teff,MdH,logg,Rp,Cloud,haze,cloudtop,art_haze_alpha,kzz,fsed,sigma_g,fitabundances
     spawn,'cp ARC.transitspec99.ncdf '+dirs+'ARC.transitspec.Li.ncdf'  
   endif

   ;Rb
   if (Opacity_switch_original(10) eq 1) then begin
     Opacity_names =['H2','He','H2O','CO2','CO','CH4','NH3','Na','K','Li','Rb','Cs','TiO','VO','FeH','PH3','H2S','HCN','C2H2','SO2','Fe','H_']
     Opacity_switch=[   0,   0,    0,    0,   0,    0,    0,   0,  0,   0,   1,   0,    0,   0,   0,    0,     0,    0,     0,    0,   0,   0] ;Hat26
     spawn,'rm ARC.transitspec99.ncdf'
     arc_make_atmo_in_molecules,teff,MdH,logg,Rp,Cloud,haze,cloudtop,art_haze_alpha,kzz,fsed,sigma_g,fitabundances
     spawn,'cp ARC.transitspec99.ncdf '+dirs+'ARC.transitspec.Rb.ncdf'  
   endif

   ;Cs
   if (Opacity_switch_original(11) eq 1) then begin
     Opacity_names =['H2','He','H2O','CO2','CO','CH4','NH3','Na','K','Li','Rb','Cs','TiO','VO','FeH','PH3','H2S','HCN','C2H2','SO2','Fe','H_']
     Opacity_switch=[   0,   0,    0,    0,   0,    0,    0,   0,  0,   0,   0,   1,    0,   0,   0,    0,     0,    0,     0,    0,   0,   0] ;Hat26
     spawn,'rm ARC.transitspec99.ncdf'
     arc_make_atmo_in_molecules,teff,MdH,logg,Rp,Cloud,haze,cloudtop,art_haze_alpha,kzz,fsed,sigma_g,fitabundances
     spawn,'cp ARC.transitspec99.ncdf '+dirs+'ARC.transitspec.Cs.ncdf'  
   endif

   ;TiO
   if (Opacity_switch_original(12) eq 1) then begin
     Opacity_names =['H2','He','H2O','CO2','CO','CH4','NH3','Na','K','Li','Rb','Cs','TiO','VO','FeH','PH3','H2S','HCN','C2H2','SO2','Fe','H_']
     Opacity_switch=[   0,   0,    0,    0,   0,    0,    0,   0,  0,   0,   0,   0,    1,   0,   0,    0,     0,    0,     0,    0,   0,   0] ;Hat26
     spawn,'rm ARC.transitspec99.ncdf'
     arc_make_atmo_in_molecules,teff,MdH,logg,Rp,Cloud,haze,cloudtop,art_haze_alpha,kzz,fsed,sigma_g,fitabundances
     spawn,'cp ARC.transitspec99.ncdf '+dirs+'ARC.transitspec.TiO.ncdf'  
   endif

   ;VO
   if (Opacity_switch_original(13) eq 1) then begin
     Opacity_names =['H2','He','H2O','CO2','CO','CH4','NH3','Na','K','Li','Rb','Cs','TiO','VO','FeH','PH3','H2S','HCN','C2H2','SO2','Fe','H_']
     Opacity_switch=[   0,   0,    0,    0,   0,    0,    0,   0,  0,   0,   0,   0,    0,   1,   0,    0,     0,    0,     0,    0,   0,   0] ;Hat26
     spawn,'rm ARC.transitspec99.ncdf'
     arc_make_atmo_in_molecules,teff,MdH,logg,Rp,Cloud,haze,cloudtop,art_haze_alpha,kzz,fsed,sigma_g,fitabundances
     spawn,'cp ARC.transitspec99.ncdf '+dirs+'ARC.transitspec.VO.ncdf'  
   endif

   ;FeH
   if (Opacity_switch_original(14) eq 1) then begin
     Opacity_names =['H2','He','H2O','CO2','CO','CH4','NH3','Na','K','Li','Rb','Cs','TiO','VO','FeH','PH3','H2S','HCN','C2H2','SO2','Fe','H_']
     Opacity_switch=[   0,   0,    0,    0,   0,    0,    0,   0,  0,   0,   0,   0,    0,   0,   1,    0,     0,    0,     0,    0,   0,   0] ;Hat26
     spawn,'rm ARC.transitspec99.ncdf'
     arc_make_atmo_in_molecules,teff,MdH,logg,Rp,Cloud,haze,cloudtop,art_haze_alpha,kzz,fsed,sigma_g,fitabundances
     spawn,'cp ARC.transitspec99.ncdf '+dirs+'ARC.transitspec.FeH.ncdf'  
   endif

   ;PH3
   if (Opacity_switch_original(15) eq 1) then begin
     Opacity_names =['H2','He','H2O','CO2','CO','CH4','NH3','Na','K','Li','Rb','Cs','TiO','VO','FeH','PH3','H2S','HCN','C2H2','SO2','Fe','H_']
     Opacity_switch=[   0,   0,    0,    0,   0,    0,    0,   0,  0,   0,   0,   0,    0,   0,   0,    1,     0,    0,     0,    0,   0,   0] ;Hat26
     spawn,'rm ARC.transitspec99.ncdf'
     arc_make_atmo_in_molecules,teff,MdH,logg,Rp,Cloud,haze,cloudtop,art_haze_alpha,kzz,fsed,sigma_g,fitabundances
     spawn,'cp ARC.transitspec99.ncdf '+dirs+'ARC.transitspec.PH3.ncdf'  
   endif

   ;H2S
   if (Opacity_switch_original(16) eq 1) then begin
     Opacity_names =['H2','He','H2O','CO2','CO','CH4','NH3','Na','K','Li','Rb','Cs','TiO','VO','FeH','PH3','H2S','HCN','C2H2','SO2','Fe','H_']
     Opacity_switch=[   0,   0,    0,    0,   0,    0,    0,   0,  0,   0,   0,   0,    0,   0,   0,    0,     1,    0,     0,    0,   0,   0] ;Hat26
     spawn,'rm ARC.transitspec99.ncdf'
     arc_make_atmo_in_molecules,teff,MdH,logg,Rp,Cloud,haze,cloudtop,art_haze_alpha,kzz,fsed,sigma_g,fitabundances
     spawn,'cp ARC.transitspec99.ncdf '+dirs+'ARC.transitspec.H2S.ncdf'  
   endif

   ;HCN
   if (Opacity_switch_original(17) eq 1) then begin
     Opacity_names =['H2','He','H2O','CO2','CO','CH4','NH3','Na','K','Li','Rb','Cs','TiO','VO','FeH','PH3','H2S','HCN','C2H2','SO2','Fe','H_']
     Opacity_switch=[   0,   0,    0,    0,   0,    0,    0,   0,  0,   0,   0,   0,    0,   0,   0,    0,     0,    1,     0,    0,   0,   0] ;Hat26
     spawn,'rm ARC.transitspec99.ncdf'
     arc_make_atmo_in_molecules,teff,MdH,logg,Rp,Cloud,haze,cloudtop,art_haze_alpha,kzz,fsed,sigma_g,fitabundances
     spawn,'cp ARC.transitspec99.ncdf '+dirs+'ARC.transitspec.HCN.ncdf'  
   endif

   ;C2H2
   if (Opacity_switch_original(18) eq 1) then begin
     Opacity_names =['H2','He','H2O','CO2','CO','CH4','NH3','Na','K','Li','Rb','Cs','TiO','VO','FeH','PH3','H2S','HCN','C2H2','SO2','Fe','H_']
     Opacity_switch=[   0,   0,    0,    0,   0,    0,    0,   0,  0,   0,   0,   0,    0,   0,   0,    0,     0,    0,     1,    0,   0,   0] ;Hat26
     spawn,'rm ARC.transitspec99.ncdf'
     arc_make_atmo_in_molecules,teff,MdH,logg,Rp,Cloud,haze,cloudtop,art_haze_alpha,kzz,fsed,sigma_g,fitabundances
     spawn,'cp ARC.transitspec99.ncdf '+dirs+'ARC.transitspec.C2H2.ncdf'  
   endif

   ;SO2
   if (Opacity_switch_original(19) eq 1) then begin
     Opacity_names =['H2','He','H2O','CO2','CO','CH4','NH3','Na','K','Li','Rb','Cs','TiO','VO','FeH','PH3','H2S','HCN','C2H2','SO2','Fe','H_']
     Opacity_switch=[   0,   0,    0,    0,   0,    0,    0,   0,  0,   0,   0,   0,    0,   0,   0,    0,     0,    0,     0,    1,   0,   0] ;Hat26
     spawn,'rm ARC.transitspec99.ncdf'
     arc_make_atmo_in_molecules,teff,MdH,logg,Rp,Cloud,haze,cloudtop,art_haze_alpha,kzz,fsed,sigma_g,fitabundances
     spawn,'cp ARC.transitspec99.ncdf '+dirs+'ARC.transitspec.SO2.ncdf'  
   endif

   ;Fe
   if (Opacity_switch_original(20) eq 1) then begin
     Opacity_names =['H2','He','H2O','CO2','CO','CH4','NH3','Na','K','Li','Rb','Cs','TiO','VO','FeH','PH3','H2S','HCN','C2H2','SO2','Fe','H_']
     Opacity_switch=[   0,   0,    0,    0,   0,    0,    0,   0,  0,   0,   0,   0,    0,   0,   0,    0,     0,    0,     0,    0,   1,   0] ;Hat26
     spawn,'rm ARC.transitspec99.ncdf'
     arc_make_atmo_in_molecules,teff,MdH,logg,Rp,Cloud,haze,cloudtop,art_haze_alpha,kzz,fsed,sigma_g,fitabundances
     spawn,'cp ARC.transitspec99.ncdf '+dirs+'ARC.transitspec.Fe.ncdf'  
   endif

   ;H_
   if (Opacity_switch_original(21) eq 1) then begin
     Opacity_names =['H2','He','H2O','CO2','CO','CH4','NH3','Na','K','Li','Rb','Cs','TiO','VO','FeH','PH3','H2S','HCN','C2H2','SO2','Fe','H_']
     Opacity_switch=[   0,   0,    0,    0,   0,    0,    0,   0,  0,   0,   0,   0,    0,   0,   0,    0,     0,    0,     0,    0,   0,   1] ;Hat26
     spawn,'rm ARC.transitspec99.ncdf'
     arc_make_atmo_in_molecules,teff,MdH,logg,Rp,Cloud,haze,cloudtop,art_haze_alpha,kzz,fsed,sigma_g,fitabundances
     spawn,'cp ARC.transitspec99.ncdf '+dirs+'ARC.transitspec.Hminus.ncdf'  
   endif



   ;set back to orginally set
   Opacity_switch  = Opacity_switch_original
   Cloudtype = Cloudtype_original
   Hazetype = Hazetype_original

end
