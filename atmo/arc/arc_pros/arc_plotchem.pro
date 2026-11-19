pro arc_plotchem,dirs,filenam

Result = FILE_TEST(dirs+filenam+'.chem_eq_ARC.ncdf')  ;test if file is ARC output, if not then load as is
 if (Result eq 1) then fid=ncdf_open(dirs+filenam+'.chem_eq_ARC.ncdf') ;ARC output
 if (Result eq 0) then fid=ncdf_open(dirs+filenam+'.highres.mcmc_bestfit.chemistry.ncdf')             ;other output
ncdf_varget,fid,'abundances',abundances
ncdf_varget,fid,'pressure',pressure
ncdf_varget,fid,'molname',molname

;AtmoKapIndex=[0,61,8,14,3,11,67,107,112,130,138,134,125,127,142]+1
;           0 h2
;          61 he
;           8 h2o
;          15 CO2
;           3 CO        
;          11 CH4       
;          67 NH3       
;         107 Na        
;         112 K         
;        130 Li        
;        138 Rb        
;        134 Cs        
;        125 TiO       
;         127 VO        
;         142 FeH  
;Sept 2018
;       1: H2           2: O-3P         3: O-1D         4: CO           5: C          
;       6: CH           7: 3CH2         8: 1CH2         9: H2O         10: O2         
;      11: H2O2        12: CH4         13: H2CO        14: CH3OH       15: CO2        
;      16: CH3OOH      17: C2H2        18: C2H4        19: C2H6        20: CH2CO      
;      21: CH3CHO      22: C2H5OH      23: C2H5OOH     24: CH3COOOH    25: C3H8       
;      26: C4H8Y       27: C4H10       28: C2H5CHO     29: C3H7OH      30: C2H6CO     
;      31: C3H8CO      32: C2H3CHOZ    33: cC2H4O      34: H           35: C7H8       
;      36: OH          37: OOH         38: CH3         39: HCO         40: CH2OH      
;      41: CH3O        42: CH3OO       43: C2H         44: C2H3        45: C2H5       
;      46: CHCO        47: CH2CHO      48: CH3CO       49: C2H5O       50: C2H4OOH    
;      51: C2H5OO      52: CH3COOO     53: 1C3H7       54: 1C4H9       55: CH3OCO     
;      56: CO2H        57: 2C2H4OH     58: 1C2H4OH     59: 2C3H7       60: 2C4H9      
;      61: N2          62: He          63: Ar          64: N-4S        65: N-2D       
;      66: NH          67: NH2         68: NH3         69: NNH         70: NO         
;      71: NO2         72: N2O         73: NCN         74: HNO         75: CN         
;      76: HCN         77: H2CN        78: HCNN        79: HCNO        80: HOCN       
;      81: HNCO        82: HON         83: NCO         84: HNO2        85: HONO       
;      86: NO3         87: HONO2       88: CH3ONO      89: CH3NO2      90: CH3NO      
;      91: C3H7O       92: C4H9O       93: cC6H6       94: N2O3        95: NH2OH      
;      96: N2O4        97: N2H2        98: N2H3        99: N2H4       100: HNNO       
;     101: HNOH       102: HNO3       103: H2NO       104: CNN        105: H2CNO      
;     106: C2N2       107: HCNH       108: Na         109: NaH        110: NaO        
;     111: NaOH       112: NaCl       113: NaF        114: K          115: KH         
;     116: KO         117: KOH        118: KCl        119: KF         120: HO2        
;     121: SO         122: SO2        123: Cl         124: HCl        125: ClO        
;     126: Cl2        127: Ti         128: TiO        129: V          130: VO         
;     131: Si         132: SiH        133: SiO        134: SiS        135: S          
;     136: SH         137: H2S        138: Mg         139: MgH        140: MgS        
;     141: Al         142: AlH        143: AlOH       144: Al2O       145: Al2O3      
;     146: Fe         147: FeH        148: Cr         149: CrN        150: CrO        
;     151: Ca         152: CaH        153: CaOH       154: F          155: HF         
;     156: Li         157: LiCl       158: LiH        159: LiF        160: Cs         
;     161: CsCl       162: CsH        163: CsF        164: Rb         165: RbCl       
;     166: RbH        167: RbF        168: P          169: PH         170: PH3        
;     171: PO         172: P2         173: PS         174: PH2        175: P4O6       
;     176: H+         177: H_         178: Na+        179: K+         180: e_         
;     181: C+         182: He+        183: Ca+        184: Si+        185: Ca(a)      
;     186: Ca(b)      187: Ca(l)      188: Ti3O5(a)   189: Ti3O5(b)   190: Ti3O5(l)   
;     191: Ti(a)      192: Ti(b)      193: Ti(l)      194: TiO(a)     195: TiO(b)     
;     196: TiO(l)     197: V(cr)      198: V(l)       199: VO(s)      200: VO(l)      
;     201: V2O3(s)    202: V2O3(l)    203: V2O4(a)    204: V2O4(b)    205: V2O4(l)    
;     206: V2O5(s)    207: V2O5(l)    208: Cr(a)      209: Cr(b)      210: Cr(l)      
;     211: Na2S(a)    212: Na2S(b)    213: Na2S(l)    214: Na3AlF6(a) 215: Na3AlF6(b) 
;     216: Na3AlF6(l) 217: K2S(a)     218: K2S(b)     219: K2S(c)     220: K2S(l)     
;     221: KCl(l)     222: KCl(s)     223: Na(l)      224: Na(cr)     225: K(l)       
;     226: K(cr)      227: S(l)       228: S(cr1)     229: S(cr2)     230: NaCl(l)    
;     231: NaCl(s)    232: MgSiO3(a)  233: MgSiO3(b)  234: MgSiO3(c)  235: MgSiO3(l)  
;     236: Mg2SiO4(s) 237: Mg2SiO4(l) 238: MgAl2O4(s) 239: MgAl2O4(l) 240: Al2O3(a)   
;     241: Al2O3(l)   242: SiO2(hqz)  243: SiO2(l)    244: NaAlSi3O8( 245: Na2O(l)    
;     246: Na2O(a)    247: Na2O(c)    248: NaAlSiO4(  249: KAlSi2O6(  250: Fe(a)      
;     251: Fe(b)      252: Fe(c)      253: Fe(d)      254: Fe(l)      255: FeS(a)     
;     256: FeS(b)     257: FeS(c)     258: FeS(l)     259: H2O(l)     260: H2O(s)     
;     261: NH3(l)     262: NH3(s)     263: Cr2O3(a)   264: Cr2O3(b)   265: Cr2O3(c)   
;     266: Cr2O3(d)   267: Cr2O3(l)   268: LiAlO2(s)  269: LiAlO2(l)  270: LiF(cr)    
;     271: CsCl(a)    272: CsCl(b)    273: CsCl(l)    274: RbCl(cr)   275: RbCl(l)    
;     276: H3PO4(l)   277: H3PO4(cr)  

set_plot,'ps'
!p.font=0
loadct,0
device,filename=dirs+filenam+'.chem_eq_ARC.ps',/helvetica,/landscape,/color
ylabelposs=16
xminimum=1E-12;-8
xmaximum=1
bottom=n_elements(abundances(*,0))-1 

;select molecule in abundances-pressure array by index in list above -1
;e.x. H3PO4 is abundances(*,277-1)

plot,abundances(*,0),pressure*1E-6,/ylog,yrange=[40,1E-5],xrange=[xminimum,1],/xlog,/xstyle,/ystyle,thick=3.0,$
  ytitle=textoidl('Pressure (bar)'),xtitle=textoidl('Mixing Ratio'),$
  title='',charthick=4,xthick=4,ythick=4

oplot,abundances(*,1-1),pressure*1E-6,thick=3;,color=cgcolor('PUR5', !D.Table_Size-2)  ;H2
xyouts,abundances(bottom,0),ylabelposs,textoidl('H_2'),charthick=4;,color=cgcolor('Orange', !D.Table_Size-2) 

oplot,abundances(*,177-1),pressure*1E-6,color=cgcolor('TG6', !D.Table_Size-2),thick=3  ;H-
xyouts,abundances(bottom,177-1),ylabelposs,'H-',color=cgcolor('TG6', !D.Table_Size-2) 


if (abundances(bottom,34-1) gt xminimum) then oplot,abundances(*,33),pressure*1E-6,color=cgcolor('PUR5', !D.Table_Size-2),thick=3  ;H
if (abundances(bottom,34-1) gt xminimum) then xyouts,abundances(bottom,33),ylabelposs+10,'H',color=cgcolor('PUR5', !D.Table_Size-2) ,charthick=4

oplot,abundances(*,62-1),pressure*1E-6,color=cgcolor('Orange', !D.Table_Size-2),thick=3  ;He
xyouts,abundances(bottom,62-1),ylabelposs+10,'He',color=cgcolor('Orange', !D.Table_Size-2) 

oplot,abundances(*,9-1),pressure*1E-6,color=cgcolor('Red', !D.Table_Size-2),thick=3  ;H2O
xyouts,abundances(bottom,9-1),ylabelposs,textoidl('H_2O'),color=cgcolor('Red', !D.Table_Size-2) 

oplot,abundances(*,4-1),pressure*1E-6,color=cgcolor('Blue', !D.Table_Size-2),thick=3 ;CO
xyouts,abundances(bottom,4-1),ylabelposs+10,textoidl('CO'),color=cgcolor('Blue', !D.Table_Size-2) 

oplot,abundances(*,15-1),pressure*1E-6,color=cgcolor('dodger blue', !D.Table_Size-2),thick=3 ;CO2
xyouts,abundances(bottom,15-1),ylabelposs+10,textoidl('CO_2'),color=cgcolor('dodger blue', !D.Table_Size-2) 

oplot,abundances(*,12-1),pressure*1E-6,color=cgcolor('Green', !D.Table_Size-2),thick=3 ;CH4
xyouts,abundances(bottom,12-1),ylabelposs,textoidl('CH_4'),color=cgcolor('Green', !D.Table_Size-2)

oplot,abundances(*,68-1),pressure*1E-6,color=cgcolor('Purple', !D.Table_Size-2),thick=3  ;NH3
xyouts,abundances(bottom,68-1),ylabelposs+10,textoidl('NH_3'),color=cgcolor('Purple', !D.Table_Size-2)

oplot,abundances(*,108-1),pressure*1E-6,color=cgcolor('Gold', !D.Table_Size-2),thick=3  ;Na
xyouts,abundances(bottom,108-1),ylabelposs+10,'Na',color=cgcolor('Gold', !D.Table_Size-2) 

oplot,abundances(*,178-1),pressure*1E-6,color=cgcolor('Gold', !D.Table_Size-2),thick=3,linestyle=2 ;Na-
xyouts,abundances(bottom,178-1),ylabelposs,'Na+',color=cgcolor('Gold', !D.Table_Size-2) 
;oplot,abundances(*,180-1),pressure*1E-6,color=cgcolor('Gold', !D.Table_Size-2),thick=3,linestyle=2 ;Na-
;xyouts,abundances(bottom,180-1),ylabelposs,'3-',color=cgcolor('Gold', !D.Table_Size-2) 

oplot,abundances(*,114-1),pressure*1E-6,color=cgcolor('Sea Green', !D.Table_Size-2),thick=3  ;K
xyouts,abundances(bottom,114-1),ylabelposs+10,'K',color=cgcolor('Sea Green', !D.Table_Size-2) 

;oplot,abundances(*,130),pressure*1E-6 ;Li
;oplot,abundances(*,138),pressure*1E-6 ;Rb
;oplot,abundances(*,134),pressure*1E-6 ;Cs
oplot,abundances(*,128-1),pressure*1E-6,color=cgcolor('Violet Red', !D.Table_Size-2),thick=3  ;TiO
xyouts,abundances(bottom,128-1),ylabelposs+10,'TiO',color=cgcolor('Violet Red', !D.Table_Size-2) 

oplot,abundances(*,130-1),pressure*1E-6,color=cgcolor('Maroon', !D.Table_Size-2),thick=3  ;VO
xyouts,abundances(bottom,130-1),ylabelposs,'VO',color=cgcolor('Maroon', !D.Table_Size-2) 

oplot,abundances(*,147-1),pressure*1E-6,color=cgcolor('TG5', !D.Table_Size-2),thick=3  ;FeH
xyouts,abundances(bottom,147-1),ylabelposs,'FeH',color=cgcolor('TG5', !D.Table_Size-2) 

device,/close

;  spawn,'ps2eps -B -C -f '+filenam+'chem_eq_ARC.ps'
  spawn,'ps2pdf -dEPSCrop -dAutoRotatePages=/None -dAutoFilterColorImages=false -dColorImageFilter=/FlateEncode -dUseFlateCompression=true '+dirs+filenam+'.chem_eq_ARC.ps '+dirs+filenam+'.chem_eq_ARC.pdf'
  spawn,'rm '+dirs+filenam+'.chem_eq_ARC.ps'
;  spawn,'pdftk temporary.pdf cat 1L output '+filenam+'.pdf'
;  spawn,'open '+filenam+'.pdf'


end

