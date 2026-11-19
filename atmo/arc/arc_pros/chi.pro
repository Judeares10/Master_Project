pro chi,y,model,err,dof,chiout
;
; Pro Chi
; 
; PURPOSE:  output chi^2 statistics, and tell you how bad the fit is
;
; INPUT:
; y     = data array
; err   = data uncertainty
; model = model points
; dof   = degrees of freedom of the fit
;
; OUTPUT:
;
; chiout = chi^2 value 
;
; dks 2,10,2010
;

print,'==================================='
chi2=total((y-model)^2D0/err^2.D0)  ;calculate chi^2 from data and model
chi2v=chi2/dof                      ;calculate reduced chi^2

print,'Chi^2  = ',chi2
print,'Chi_v^2 = ',chi2v

p=chisqr_pdf(chi2,dof)                  ;Chi^2 table of probabilities cumulative distribution function for a chi-square distribution with Df degrees of freedom
pvalue=1.-p
confidence=chisqr_cvf(1.0D0-p,1)^0.5D0  ;convert p-value to sigma

print,'==================================='

print,'Good Chi_v^2 fit has values < ',1D0+sqrt(2D0/dof)

; check to see if the fit is "good" or "bad" ie good fit has: chi2v < 1 + (2/DOF)^1/2
; comment on quality
if ( chi2v le 1D0+sqrt(2D0/dof)) and (chi2v ge 1D0-sqrt(2D0/dof)) then print,'Chi_v^2 Does indicates a good fit!! :-)'
if ( chi2v ge 1D0+sqrt(2D0/dof)) then print,'Chi_v^2 indicates a bad fit, sorry'
if ( chi2v le 1D0-sqrt(2D0/dof)) then print,'Chi_v^2 indicates a bad fit, sorry, your errors are to small'

print,'P-value =',pvalue
print,'Sigma   =',confidence
print,'Prob to exceed Chi^2 value (randomly find worse fit) is :',100-p*100 ,' percent'
print,'Prob to find lower Chi^2 value (randomly find better fit) is :',p*100 ,' percent'
print,'So the model can be rejected at the ',confidence ,' sigma confidence level'

print,'==================================='

chiout=chi2

end
