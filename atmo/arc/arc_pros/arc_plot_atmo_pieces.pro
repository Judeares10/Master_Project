pro arc_plot_atmo_pieces,dirs,filenam,teff,MdH,COratio,Rp,Cloud,MPFITfilenameprefix,Rpl,Rs,logg
  Common SystemParameters,rplanet,rstar,rorbit,stellarfile,Tstar,Cloudtype,Hazetype,Chemistry,Isothermal,Opacity_switch,Abund_fit_switch,Abund_fit_index,Abund_priorLow,Abund_priorUP,rain_switch,Element_fit_switch,Element_fit_index,Ele_priorLow,Ele_priorUP,Distance,nband

  loadct,0
  g=10^logg
  FitTemp=teff
   ;arc_load_system_parameters,planet,Rp,Rjup,Rpl,Teff,logg,g,Rstar,Rs,rorbit,stellarfile
  ;====
     K     = 1.380658D-16    ; gm*cm^2/s^2 * Kel
     mu    = 1.6726D-24*2.3  ; cgs
     Tau   = 0.56D0          ; Optical depth where Atmo Pressure corresponds to measured Transit Radius
     Rjup  = 71492D5         ; cm
     alpha=-4
     zoffset=0;pp_Rp
     T=FitTemp
  ;====
   ;Opacity_names   =['H2','He','H2O','CO2','CO','CH4','NH3','Na','K','Li','Rb','Cs','TiO','VO','FeH','PH3','H2S','HCN','C2H2']
   ; Switch:      1=include  0=Omit
   ;Opacity_switch  =[   1,   1,    1,    0,   1,    1,    1,   0,  0,   0,   0,   0,    0,   0,    0,   0,     0,    0,     0] ;Hat26
  
  ;==Load Full model=======
  atmofilenam=dirs+MPFITfilenameprefix+'.highres.mcmc_bestfit.transitspec.ncdf'
  arc_load_atmo_spec,atmofilenam,ws,model
  fullmodel=model
  modelmin=min(fullmodel)

  ;==LOAD H2 Model=======
  if (Opacity_switch(0) eq 1) then begin
   atmofilenam=dirs+'ARC.transitspec.H2.ncdf'  
   arc_load_atmo_spec,atmofilenam,ws,model
   ModelH2=model
   if (min(model) lt modelmin) then modelmin=min(model)
  endif

  ;==LOAD Cloud Model=======
  if (Cloudtype eq 1) then begin
    atmofilenam=dirs+'ARC.transitspec.Cloud.ncdf'  
    arc_load_atmo_spec,atmofilenam,ws,model
    ModelCloud=model
    if (min(model) lt modelmin) then modelmin=min(model)
  endif

  ;==LOAD Cloud Model=======
  if (hazetype eq 1) then begin
    atmofilenam=dirs+'ARC.transitspec.Haze.ncdf'  
    arc_load_atmo_spec,atmofilenam,ws,model
    ModelHaze=model
    if (min(model) lt modelmin) then modelmin=min(model)
  endif

  ;==LOAD He Model=======
  if (Opacity_switch(1) eq 1) then begin
    atmofilenam=dirs+'ARC.transitspec.He.ncdf'  
    arc_load_atmo_spec,atmofilenam,ws,model
    ModelHe=model
    if (min(model) lt modelmin) then modelmin=min(model)
  endif

  ;==LOAD H2O Model=======
  if (Opacity_switch(2) eq 1) then begin
    atmofilenam=dirs+'ARC.transitspec.H2O.ncdf'  
    arc_load_atmo_spec,atmofilenam,ws,model
    ModelH2O=model
    if (min(model) lt modelmin) then modelmin=min(model)
  endif


  ;==LOAD CO2 Model======= 
  if (Opacity_switch(3) eq 1) then begin
    atmofilenam=dirs+'ARC.transitspec.CO2.ncdf'  
    arc_load_atmo_spec,atmofilenam,ws,model
    ModelCO2=model
    if (min(model) lt modelmin) then modelmin=min(model)
  endif

  ;==LOAD CO Model=======
  if (Opacity_switch(4) eq 1) then begin
    atmofilenam=dirs+'ARC.transitspec.CO.ncdf'  
    arc_load_atmo_spec,atmofilenam,ws,model
    ModelCO=model
    if (min(model) lt modelmin) then modelmin=min(model)
 endif

  ;==LOAD CH4 Model=======
  if (Opacity_switch(5) eq 1) then begin
    atmofilenam=dirs+'ARC.transitspec.CH4.ncdf'  
    arc_load_atmo_spec,atmofilenam,ws,model
    ModelCH4=model
    if (min(model) lt modelmin) then modelmin=min(model)
   endif


  ;==LOAD NH3 Model=======
  if (Opacity_switch(6) eq 1) then begin
    atmofilenam=dirs+'ARC.transitspec.NH3.ncdf'  
    arc_load_atmo_spec,atmofilenam,ws,model
    ModelNH3=model
    if (min(model) lt modelmin) then modelmin=min(model)
   endif

  ;==LOAD Na Model=======
  if (Opacity_switch(7) eq 1) then begin
    atmofilenam=dirs+'ARC.transitspec.Na.ncdf'  
    arc_load_atmo_spec,atmofilenam,ws,model
    ModelNa=model
    if (min(model) lt modelmin) then modelmin=min(model)
   endif

  ;==LOAD K Model=======
  if (Opacity_switch(8) eq 1) then begin
    atmofilenam=dirs+'ARC.transitspec.K.ncdf'  
    arc_load_atmo_spec,atmofilenam,ws,model
    ModelK=model
    if (min(model) lt modelmin) then modelmin=min(model)
 endif

  ;==LOAD Li Model=======
  if (Opacity_switch(9) eq 1) then begin
    atmofilenam=dirs+'ARC.transitspec.Li.ncdf'  
    arc_load_atmo_spec,atmofilenam,ws,model
    ModelLi=model
    if (min(model) lt modelmin) then modelmin=min(model)
 endif

  ;==LOAD Rb Model=======
  if (Opacity_switch(10) eq 1) then begin
    atmofilenam=dirs+'ARC.transitspec.Rb.ncdf'  
    arc_load_atmo_spec,atmofilenam,ws,model
    ModelRb=model
    if (min(model) lt modelmin) then modelmin=min(model)
 endif

  ;==LOAD Cs Model=======
  if (Opacity_switch(11) eq 1) then begin
    atmofilenam=dirs+'ARC.transitspec.Cs.ncdf'  
    arc_load_atmo_spec,atmofilenam,ws,model
    ModelCs=model
    if (min(model) lt modelmin) then modelmin=min(model)
 endif

  ;==LOAD TiO Model=======
  if (Opacity_switch(12) eq 1) then begin
    atmofilenam=dirs+'ARC.transitspec.TiO.ncdf'  
    arc_load_atmo_spec,atmofilenam,ws,model
    ModelTiO=model
    if (min(model) lt modelmin) then modelmin=min(model)
 endif

  ;==LOAD VO Model=======
  if (Opacity_switch(13) eq 1) then begin
    atmofilenam=dirs+'ARC.transitspec.VO.ncdf'  
    arc_load_atmo_spec,atmofilenam,ws,model
    ModelVO=model
    if (min(model) lt modelmin) then modelmin=min(model)
 endif
  ;==LOAD FeH  Model=======
  if (Opacity_switch(14) eq 1) then begin
    atmofilenam=dirs+'ARC.transitspec.FeH.ncdf'  
    arc_load_atmo_spec,atmofilenam,ws,model
    ModelFeH=model
    if (min(model) lt modelmin) then modelmin=min(model)
 endif

  ;==LOAD PH3 Model=======
  if (Opacity_switch(15) eq 1) then begin
    atmofilenam=dirs+'ARC.transitspec.PH3.ncdf'  
    arc_load_atmo_spec,atmofilenam,ws,model
    ModelPH3=model
    if (min(model) lt modelmin) then modelmin=min(model)
   endif

  ;==LOAD H2S Model=======
  if (Opacity_switch(16) eq 1) then begin
    atmofilenam=dirs+'ARC.transitspec.H2S.ncdf'  
    arc_load_atmo_spec,atmofilenam,ws,model
    ModelH2S=model
    if (min(model) lt modelmin) then modelmin=min(model)
   endif


  ;==LOAD HCN Model=======
  if (Opacity_switch(17) eq 1) then begin
    atmofilenam=dirs+'ARC.transitspec.HCN.ncdf'  
    arc_load_atmo_spec,atmofilenam,ws,model
    ModelHCN=model
    if (min(model) lt modelmin) then modelmin=min(model)
  endif

  ;==LOAD C2H2 Model=======
  if (Opacity_switch(18) eq 1) then begin
    atmofilenam=dirs+'ARC.transitspec.C2H2.ncdf'  
    arc_load_atmo_spec,atmofilenam,ws,model
    ModelC2H2=model
    if (min(model) lt modelmin) then modelmin=min(model)
  endif

    ;==LOAD SO2 Model=======
  if (Opacity_switch(19) eq 1) then begin
    atmofilenam=dirs+'ARC.transitspec.SO2.ncdf'  
    arc_load_atmo_spec,atmofilenam,ws,model
    ModelSO2=model
    if (min(model) lt modelmin) then modelmin=min(model)
  endif

    ;==LOAD Fe Model=======
  if (Opacity_switch(20) eq 1) then begin
    atmofilenam=dirs+'ARC.transitspec.Fe.ncdf'  
    arc_load_atmo_spec,atmofilenam,ws,model
    ModelFe=model
    if (min(model) lt modelmin) then modelmin=min(model)
  endif

    ;==LOAD H- Model=======
  if (Opacity_switch(21) eq 1) then begin
    atmofilenam=dirs+'ARC.transitspec.Hminus.ncdf'  
    arc_load_atmo_spec,atmofilenam,ws,model
    ModelHminus=model
    if (min(model) lt modelmin) then modelmin=min(model)
  endif


  ;==LOAD Cloud Deck Model=======
  ;ws defined from above, cloud deck heigh applied
  if (Cloudtype eq 2) then begin
     ModelCloud=ws/ws*Cloud; rprs+zoffset
     ;cloudindex=where(Model lt Cloud)
     ;Model(cloudindex)=Cloud
     ;ModelCloud=model
  endif


ymin=modelmin-0.001;Rpl/Rs*0.75
ymax=max(fullmodel)+0.001;Rpl/Rs*1.15

;===PLOT====
  set_plot,'ps'
  !p.font=0
  loadct,0
  device,filename=dirs+filenam+'compare_pieces_Spec.ps',/helvetica,/landscape,/color
  plot,ws,ModelH2,xrange=[0.7,6],yrange=[ymin,ymax],/xsty,/xlog,thick=2.0,/ysty,/nodata,$
    ytitle=textoidl('planet-to-star radius ratio (R_{pl}/R_{s})'),xtitle=textoidl('Wavelength (\mum)'),$
    title='',charthick=4,xthick=4,ythick=4,$
    Xtickv=[.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.5,3,3.5,4,5],xticks=22,xtickname=['.3','.4','.5','.6','.7','.8','.9','1',' ',' ',' ',' ','1.5',' ',' ',' ',' ','2','2.5','3','3.5','4','5'],xticklen=0.01
  if (Opacity_switch(0) eq 1) then oplot,ws,ModelH2,thick=3,color=cgcolor('black', !D.Table_Size-2)
  if (Cloudtype eq 1) then oplot,ws,ModelCloud,thick=3,color=cgcolor('Dark Grey', !D.Table_Size-2)
  if (Cloudtype eq 2) then oplot,ws,ModelCloud,thick=3,color=cgcolor('Dark Grey', !D.Table_Size-2)
  if (Hazetype eq 1) then oplot,ws,ModelHaze,thick=3,color=cgcolor('Blue Violet', !D.Table_Size-2)
  if (Opacity_switch(1) eq 1) then oplot,ws,ModelHe,thick=3,color=cgcolor('blue', !D.Table_Size-2)  
  if (Opacity_switch(2) eq 1) then oplot,ws,ModelH2O,thick=3,color=cgcolor('red', !D.Table_Size-2)
  if (Opacity_switch(3) eq 1) then oplot,ws,ModelCO2,thick=3,color=cgcolor('orange', !D.Table_Size-2)
  if (Opacity_switch(4) eq 1) then oplot,ws,ModelCO,thick=3,color=cgcolor('dark green', !D.Table_Size-2)
  if (Opacity_switch(5) eq 1) then oplot,ws,ModelCH4,thick=3,color=cgcolor('purple', !D.Table_Size-2)
  if (Opacity_switch(6) eq 1) then oplot,ws,ModelNH3,thick=3,color=cgcolor('Gold', !D.Table_Size-2)
  if (Opacity_switch(7) eq 1) then oplot,ws,ModelNa,thick=3,color=cgcolor('green', !D.Table_Size-2)
  if (Opacity_switch(8) eq 1) then oplot,ws,ModelK,thick=3,color=cgcolor('dodger blue', !D.Table_Size-2)
  if (Opacity_switch(9) eq 1) then oplot,ws,ModelLi,thick=3,color=cgcolor('red5', !D.Table_Size-2)
  if (Opacity_switch(10) eq 1) then oplot,ws,ModelRb,thick=3,color=cgcolor('red6', !D.Table_Size-2)
  if (Opacity_switch(11) eq 1) then oplot,ws,ModelCs,thick=3,color=cgcolor('red7', !D.Table_Size-2)
  if (Opacity_switch(12) eq 1) then oplot,ws,ModelTiO,thick=3,color=cgcolor('red8', !D.Table_Size-2)
  if (Opacity_switch(13) eq 1) then oplot,ws,ModelVO,thick=3,color=cgcolor('RYB8', !D.Table_Size-2)
  if (Opacity_switch(14) eq 1) then oplot,ws,ModelFeH,thick=3,color=cgcolor('PUR5', !D.Table_Size-2)
  if (Opacity_switch(15) eq 1) then oplot,ws,ModelPH3,thick=3,color=cgcolor('red5', !D.Table_Size-2)
  if (Opacity_switch(16) eq 1) then oplot,ws,ModelH2S,thick=3,color=cgcolor('violet', !D.Table_Size-2)
  if (Opacity_switch(17) eq 1) then oplot,ws,ModelHCN,thick=3,color=cgcolor('green', !D.Table_Size-2)
  if (Opacity_switch(18) eq 1) then oplot,ws,ModelC2H2,thick=3,color=cgcolor('yellow', !D.Table_Size-2)
  if (Opacity_switch(19) eq 1) then oplot,ws,ModelSO2,thick=3,color=cgcolor('Cyan', !D.Table_Size-2)
  if (Opacity_switch(20) eq 1) then oplot,ws,ModelFe,thick=3,color=cgcolor('CG9', !D.Table_Size-2)
  if (Opacity_switch(21) eq 1) then oplot,ws,ModelHminus,thick=3,color=cgcolor('Orchid', !D.Table_Size-2)
  oplot,ws,fullmodel,thick=2,color=cgcolor('black', !D.Table_Size-2)

  device,/close
;=============
   ;                    0   1    2     3     4    5      6     7
   ;Opacity_names   =['H2','He','H2O','CO2','CO','CH4','NH3','Na','K','Li','Rb','Cs','TiO','VO','FeH','PH3','H2S','HCN','C2H2']

;measure mean molecular weight from Rayleigh scatttering
mufac=1./(Rs*g*(ModelH2(3499)-ModelH2(1000))/alog(ws(3499)/ws(1000))/(alpha*T*k))/1.6726D-24
mu=mufac*1.6726D-24
H=K*T/(mu*g)
sig0=2.52E-28; cm^2 at 750nm
P0=0.4 ;400 mbar at 750nm
modelH=ModelH2/(h/rs) ;convert to scale height units

;===PLOT====
  set_plot,'ps'
  !p.font=0
  device,filename=dirs+filenam+'.compare_pieces_Press.ps',/helvetica,/landscape,/color
  plot,ws,P0*exp(-(modelH-modelH(1333))),xrange=[0.3,6],yrange=[100,0.000001],/xsty,/xlog,/ylog,thick=2.0,/ysty,/nodata,$
    ytitle=textoidl('Pressure (bar)'),xtitle=textoidl('Wavelength (\mum)'),$
    title='',charthick=4,xthick=4,ythick=4,$
    Xtickv=[.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.5,3,3.5,4,5],xticks=22,xtickname=['.3','.4','.5','.6','.7','.8','.9','1',' ',' ',' ',' ','1.5',' ',' ',' ',' ','2','2.5','3','3.5','4','5'],xticklen=0.01
  if (Opacity_switch(0) eq 1) then oplot,ws,P0*exp(-(ModelH2/(h/rs)-modelH(1333))),thick=3,color=cgcolor('black', !D.Table_Size-2)
  if (Cloudtype eq 1) then oplot,ws,P0*exp(-(ModelCloud/(h/rs)-modelH(1333))),thick=3,color=cgcolor('Dark Grey', !D.Table_Size-2)
  if (Cloudtype eq 2) then oplot,ws,P0*exp(-(ModelCloud/(h/rs)-modelH(1333))),thick=3,color=cgcolor('Dark Grey', !D.Table_Size-2)
  if (Hazetype eq 1) then oplot,ws,P0*exp(-(ModelHaze/(h/rs)-modelH(1333))),thick=3,color=cgcolor('Blue Violet', !D.Table_Size-2)
  if (Opacity_switch(1) eq 1) then oplot,ws,P0*exp(-(ModelHe/(h/rs)-modelH(1333))),thick=3,color=cgcolor('blue', !D.Table_Size-2)  
  if (Opacity_switch(2) eq 1) then oplot,ws,P0*exp(-(ModelH2O/(h/rs)-modelH(1333))),thick=3,color=cgcolor('red', !D.Table_Size-2)
  if (Opacity_switch(3) eq 1) then oplot,ws,P0*exp(-(ModelCO2/(h/rs)-modelH(1333))),thick=3,color=cgcolor('orange', !D.Table_Size-2)
  if (Opacity_switch(4) eq 1) then oplot,ws,P0*exp(-(ModelCO/(h/rs)-modelH(1333))),thick=3,color=cgcolor('dark green', !D.Table_Size-2)
  if (Opacity_switch(5) eq 1) then oplot,ws,P0*exp(-(ModelCH4/(h/rs)-modelH(1333))),thick=3,color=cgcolor('purple', !D.Table_Size-2)
  if (Opacity_switch(6) eq 1) then oplot,ws,P0*exp(-(ModelNH3/(h/rs)-modelH(1333))),thick=3,color=cgcolor('Gold', !D.Table_Size-2)
  if (Opacity_switch(7) eq 1) then oplot,ws,P0*exp(-(ModelNa/(h/rs)-modelH(1333))),thick=3,color=cgcolor('green', !D.Table_Size-2)
  if (Opacity_switch(8) eq 1) then oplot,ws,P0*exp(-(ModelK/(h/rs)-modelH(1333))),thick=3,color=cgcolor('dodger blue', !D.Table_Size-2)
  if (Opacity_switch(9) eq 1) then oplot,ws,P0*exp(-(ModelLi/(h/rs)-modelH(1333))),thick=3,color=cgcolor('red5', !D.Table_Size-2)
  if (Opacity_switch(10) eq 1) then oplot,ws,P0*exp(-(ModelRb/(h/rs)-modelH(1333))),thick=3,color=cgcolor('red6', !D.Table_Size-2)
  if (Opacity_switch(11) eq 1) then oplot,ws,P0*exp(-(ModelCs/(h/rs)-modelH(1333))),thick=3,color=cgcolor('red7', !D.Table_Size-2)
  if (Opacity_switch(12) eq 1) then oplot,ws,P0*exp(-(ModelTiO/(h/rs)-modelH(1333))),thick=3,color=cgcolor('red8', !D.Table_Size-2)
  if (Opacity_switch(13) eq 1) then oplot,ws,P0*exp(-(ModelVO/(h/rs)-modelH(1333))),thick=3,color=cgcolor('RYB8', !D.Table_Size-2)
  if (Opacity_switch(14) eq 1) then oplot,ws,P0*exp(-(ModelFeH/(h/rs)-modelH(1333))),thick=3,color=cgcolor('PUR5', !D.Table_Size-2)
  if (Opacity_switch(15) eq 1) then oplot,ws,P0*exp(-(ModelPH3/(h/rs)-modelH(1333))),thick=3,color=cgcolor('red5', !D.Table_Size-2)
  if (Opacity_switch(16) eq 1) then oplot,ws,P0*exp(-(ModelH2S/(h/rs)-modelH(1333))),thick=3,color=cgcolor('violet', !D.Table_Size-2)
  if (Opacity_switch(17) eq 1) then oplot,ws,P0*exp(-(ModelHCN/(h/rs)-modelH(1333))),thick=3,color=cgcolor('green', !D.Table_Size-2)
  if (Opacity_switch(18) eq 1) then oplot,ws,P0*exp(-(ModelC2H2/(h/rs)-modelH(1333))),thick=3,color=cgcolor('yellow', !D.Table_Size-2)
  if (Opacity_switch(19) eq 1) then oplot,ws,P0*exp(-(ModelSO2/(h/rs)-modelH(1333))),thick=3,color=cgcolor('Cyan', !D.Table_Size-2)
  if (Opacity_switch(20) eq 1) then oplot,ws,P0*exp(-(ModelFe/(h/rs)-modelH(1333))),thick=3,color=cgcolor('CG9', !D.Table_Size-2)
  if (Opacity_switch(21) eq 1) then oplot,ws,P0*exp(-(ModelHminus/(h/rs)-modelH(1333))),thick=3,color=cgcolor('Orchid', !D.Table_Size-2)
  oplot,ws,P0*exp(-(fullmodel/(h/rs)-modelH(1333))),thick=2,color=cgcolor('black', !D.Table_Size-2)

  device,/close
;=============

end
