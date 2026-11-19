pro arc_TP_delta,ptfile,logg,Kir,gamma1,gamma2,beta,alpha,T,Tint,Press,delta,Pdelta
; parameterize T-P profile as Guillot A&A 520, A27
; and
; Retrieval as Line et al. (2013, ApJ, 775, 137, pg 6)
; 
; Assumes Two independant downward streams of visible radiation, and one infrared stream
;
; Input is:
;   1) NetCDF file of ATMO generated Pressure-Temperature profile
;   2) 5 parameters for T-P profile
;     Kir    = Planck mean thermal IR opacity (cm^2 g^-1)
;     gamma1 = ratios of the Planck mean opacities in the visible to thermal
;     gamma2 = ratios of the Planck mean opacities in the visible to thermal
;     beta   = albedo/redistribution factor
;     alpha  = partitions opt/ir flux
;   3) Internal temperature and Stellar Temperature
; Output is temperature, overwriting NetCDF file
; 
; DKS 8/16

  Common SystemParameters,rplanet,rstar,rorbit,stellarfile,Tstar,Cloudtype,Hazetype,Chemistry,Isothermal,Opacity_switch,Abund_fit_switch,Abund_fit_index,Abund_priorLow,Abund_priorUP,rain_switch,Element_fit_switch,Element_fit_index,Ele_priorLow,Ele_priorUP,Distance,nband
  Common cpcount,ChildProcessCount,dirs,allmodel,allmodelpoints,Emission,nlevel,pmin,pmax,allchemistry,allTP


  ;Common ARC_retrieval_data, x,y,err,waverr
  ;Common SystemParameters,rplanet,logg,rstar,rorbit,stellarfile,Tstar,Cloudtype,Hazetype,Chemistry,Isothermal,Opacity_switch,Abund_fit_switch,Abund_fit_index,Abund_priorLow,Abund_priorUP
  ;Common cpcount,ChildProcessCount,dirs,allmodel,allmodelpoints,Emission,nlevel,pmin,pmax

rstarcm=rstar*696000D5; cm
rorbitcm=rorbit*1.496E13 ; au->cm

;Open ATMO NCDF file of P-T profile
fid=ncdf_open(ptfile,/write)
ncdf_varget,fid,'pressure',pressure ;cgs Barye, bar is pressure*1D-6
ncdf_varget,fid,'temperature',temperature


; Line et al eq 16
Tau=Kir*Pressure/(10^logg)
; Line et al. eq 15
Tirr=beta*(rstarcm/(2*rorbitcm))^0.5*Tstar

; Equation 13 & 14 Line et al. (2013)
  ;The EXPINT function returns the value of the exponential integral En(x).
  EXPINT1 = EXPINT( 2, gamma1*tau , /DOUBLE)  ;, EPS=value, ITER=variable, ITMAX=value] )
  EXPINT2 = EXPINT( 2, gamma2*tau , /DOUBLE)  ;, EPS=value, ITER=variable, ITMAX=value] )
  E1=2./3+2./(3*gamma1)*(1+(gamma1*tau/2.-1)*exp(-gamma1*tau))+2*gamma1/3.*(1-tau^2./2.)*EXPINT1
  E2=2./3+2./(3*gamma2)*(1+(gamma2*tau/2.-1)*exp(-gamma2*tau))+2*gamma2/3.*(1-tau^2./2.)*EXPINT2
temperature=(3*Tint^4./4.*(2./3+tau)+3*Tirr^4./4*(1-alpha)*E1+3*Tirr^4./4*alpha*E2)^0.25

press=pressure*1D-6 ;bar
T=temperature       ;Kelvin
;Use Mollier paramerization to add flexability at top of atmosphere
;print,'passed params=',delta,pdelta
T=T*(1D0+delta/(1d0+press/Pdelta))
temperature=T



;Write Analytic Temperature to NCDF file (to feed ATMO)
NCDF_VARPUT,fid,'temperature',temperature
NCDF_CLOSE, fid 

end
