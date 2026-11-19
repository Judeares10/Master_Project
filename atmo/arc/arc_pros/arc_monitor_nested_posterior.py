from __future__ import unicode_literals
import numpy as np
from os import path,mkdir
import os
from scipy.io import readsav
from matplotlib.pyplot import *
from matplotlib.backends.backend_pdf import PdfPages
from scipy.interpolate import interp1d
import sys
import glob
from matplotlib.patches import Rectangle
import corner as corner

rcParams_orig=rcParams
#----------------------------------------------------------------------------------------
# arc_monitor_nested_posterior.py
#
# Description:
#   Program Outputs corner posterior plot generated from a Nested sampling run in IDL.
#   IDL save file with nested sampling run is read in, and posterior T-P profile calculated.
#
# Usage:
#   Launch from command line as:
#   $ python arc_monitor_nested_plot_pt.py <directory_path_name>
#
# Outputs:
#   <directory_path_name>/arc_monitor_nested_PT_profile.pdf
#
# Last edited: DKS 7 July 2020
#-----------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------
# Read in path argument.  If absent use current directory
print('Number of arguments:', len(sys.argv), 'arguments.')
print('Argument List:', str(sys.argv))
if len(sys.argv) > 1:
    dirs=sys.argv[1]
else:
    dirs='./'
print('Working in Directory: ',dirs)
#-----------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------
# Read in IDL save file
list=glob.glob(dirs+"*.MULTINEST.sav")
number_of_images=len(list)         #Count number of images

sav = readsav(os.path.join(list[0]))
print('Reading in File: ',list[0])

#>>> print(sav.keys())
#dict_keys(['multinest_str', 'nestpars', 'lp', 'log_likelihood_func', 'prior_func', 'num_params', 'num', 'sample_num', 'expand', 'tol', 'tofit'])
#multinest_str is an IDL structure, rest are normal arrays
multinest_str = sav['multinest_str']
nestpars=sav['nestpars']
tofit = sav['tofit']
#Show the variables in the multinest_str structure
#>>> multinest_str.dtype.names  
#('SAMPLES', 'LOG_PROB', 'LOGZ', 'LOGZ_ERR', 'H', 'SAMPLES_MODEL', 'SAMPLES_ABUNDANCES', 'SAMPLES_TEMPERATURE', 'WS', 'PRESSURE', 'X', 'Y', 'ERR', 'WAVERR', 'MOLNAME')

pressure   =multinest_str['PRESSURE'][0]*1E-6       # barye -> bar 

temperature=multinest_str['SAMPLES_TEMPERATURE'][0] # for IDL structures, there is an extra array ontop of the variable you need to get rid of, so select [0]
prob=np.exp(multinest_str['LOG_PROB'][0])

#Get futher info
list=glob.glob(dirs+"*arc_retrieval.setup.sav")
number_of_images=len(list)         #Count number of images

sav = readsav(os.path.join(list[0]))
print('Reading in File: ',list[0])
UPPRIORSPARAMS = sav['UPPRIORSPARAMS']
LOWPRIORSPARAMS = sav['LOWPRIORSPARAMS']
Chemistry = sav['Chemistry']
Element_fit_switch = sav['Element_fit_switch'].astype(int)
ELEMENT_FIT_INDEX = sav['ELEMENT_FIT_INDEX'].astype(int)

#-----------------------------------------------------------------------------------------
#Define all the fit parameters possible for the different chemistries
parameternames=['Temperature (K)','ln(Z/Z$_{\odot}$)','log(g)','R$_{pl}$(Jup)','ln($\delta_{cloud}$)','ln($\delta_{haze}$)','Cloud Top (bar)','$\\alpha_{haze}$','K$_{ZZ}$ (cm$^2$ s)','fsed','$\sigma$_g','log($\kappa_{IR}$)','log($\gamma_{O/IR}$)','log($\gamma_{O/IR}$)','$\\beta$','$alpha$', 'Temp Int (K)','$\delta$','P$\delta$']

parameternames_main=parameternames

if Chemistry == 1:
    Molecule_names_ln    =  ['VMR ln(H$\tex_2$)','VMR ln(He)','VMR ln(H$\tex_2$O)','VMR ln(CO$\tex_2$)','VMR ln(CO)','VMR ln(CH$\tex_4$)','VMR ln(NH$\tex_3$)','VMR ln(Na)','VMR ln(K)','VMR ln(Li)','VMR ln(Rb)','VMR ln(Cs)','VMR ln(TiO)','VMR ln(VO)','VMR ln(FeH)','VMR ln(PH$\tex_3$)','VMR ln(H$\tex_2$S)','VMR ln(HCN)','VMR ln(C$\tex_2$H$\tex_2$)','VMR ln(SO$\tex_2$)','VMR ln(Fe)','VMR ln(H_)']
    Molecule_names_log   =  ['VMR log(H$\tex_2$)','VMR log(He)','VMR log(H$\tex_2$O)','VMR log(CO$\tex_2$)','VMR log(CO)','VMR log(CH$\tex_4$)','VMR log(NH$\tex_3$)','VMR log(Na)','VMR log(K)','VMR log(Li)','VMR log(Rb)','VMR log(Cs)','VMR log(TiO)','VMR log(VO)','VMR log(FeH)','VMR log(PH$\tex_3$)','VMR log(H$\tex_2$S)','VMR log(HCN)','VMR log(C$\tex_2$H$\tex_2$)','VMR ln(SO$\tex_2$)','VMR ln(Fe)','VMR log(H_)']
    parameternames=[parameternames,Molecule_names_ln]
    parameternameslog=[parameternames_main,Molecule_names_log]

if Chemistry == 3:
    Element_names      =  ['H','He','C','N','O','Na','K','Si','Ar','Ti','V','S','Cl','Mg','Al','Ca','Fe','Cr','Li','Cs','Rb','F','P'] #Spectrally active molecules
    Element_names_ln   =  ['ln(H/H$_{\odot}$)','ln(He/He$_{\odot}$)','ln(C/C$_{\odot}$)','ln(N/N$_{\odot}$)','ln(O/O$_{\odot}$)','ln(Na/Na$_{\odot}$)','ln(K/K$_{\odot}$)','ln(Si/Si$_{\odot}$)','ln(Ar/Ar$_{\odot}$)','ln(Ti/Ti$_{\odot}$)','ln(V/V$_{\odot}$)','ln(S/S$_{\odot}$)','ln(Cl/Cl$_{\odot}$)','ln(Mg/Mg$_{\odot}$)','ln(Al/Al$_{\odot}$)','ln(Ca/Ca$_{\odot}$)','ln(Fe/Fe$_{\odot}$)','ln(Cr/Cr$_{\odot}$)','ln(Li/Li$_{\odot}$)','ln(Cs/Cs$_{\odot}$)','ln(Rb/Rb$_{\odot}$)','ln(F/F$_{\odot}$)','ln(P/P$_{\odot}$)'] #Spectrally active molecules
    Element_names_log = ['log(H/H$_{\odot}$)','log(He/He$_{\odot}$)','log(C/C$_{\odot}$)','log(N/N$_{\odot}$)','log(O/O$_{\odot}$)','log(Na/Na$_{\odot}$)','log(K/K$_{\odot}$)','log(Si/Si$\tex_{\odot}$)','log(Ar/Ar$_{\odot}$)','log(Ti/Ti$_{\odot}$)','log(V/V$_{\odot}$)','log(S/S$_{v}$)','log(Cl/Cl$_{\odot}$)','log(Mg/Mg$\tex_{\odot}$)','log(Al/Al$_{\odot}$)','log(Ca/Ca$_{\odot}$)','log(Fe/Fe$_{\odot}$)','log(Cr/Cr$_{\odot}$)','log(Li/Li$\tex_{\odot}$)','log(Cs/Cs$_{\odot}$)','log(Rb/Rb$_{\odot}$)','log(F/F$_{\odot}$)','log(P/P$_{\odot}$)'] #Spectrally active molecules
    Element_names      = [Element_names[i] for i in ELEMENT_FIT_INDEX-1]  #select only list values which were fit
    Element_names_ln   = [Element_names_ln[i] for i in ELEMENT_FIT_INDEX-1]
    Element_names_log  = [Element_names_log[i] for i in ELEMENT_FIT_INDEX-1]
    parameternames= parameternames+  Element_names_ln       # parameternames=[parameternames,Element_names_ln]
    parameternameslog=parameternames_main+Element_names_log # parameternameslog=[parameternames_main,Element_names_log]
    parameternameslog[1]='log(Z/Z$_{\odot}$)'

p=[parameternames[i] for i in tofit]
parameternames=p
parameternameslog=[parameternameslog[i] for i in tofit]

#-----------------------------------------------------------------------------------------

#_ = ylim(1.35,1.42)
#_ = xlim(-1.5,1.75)


nfree=len(nestpars[0,:])
lenth=len(nestpars)


#Make Chains with Log10 instead of ln
freefit=(tofit >= 19).nonzero()
freefit=freefit[0]
LNESTPARS= np.zeros((lenth,nfree),dtype=float) #nestpars
LNESTPARS[:,:]=nestpars[:,:]
MdH=(tofit == 1) #true/false
MdHindex=(tofit == 1).nonzero() #index
if MdH[0] == True: LNESTPARS[:,MdHindex]=np.log10(np.exp(nestpars[:,MdHindex]))
for i in range(0,len(freefit)): LNESTPARS[:,freefit[i]]=np.log10(np.exp(nestpars[:,freefit[i]]))

#-----------------------------------------------------------------------------------------
# Generate Posterior
#
# Define Temperature posterior 1,2,3-sigma arrays
med=np.zeros(nfree)
sig1l=np.zeros(nfree)
sig1h=np.zeros(nfree)
sig2l=np.zeros(nfree)
sig2h=np.zeros(nfree)
sig3l=np.zeros(nfree)
sig3h=np.zeros(nfree)
sig4l=np.zeros(nfree)
sig4h=np.zeros(nfree)

# Posterior is generated by selecting a each parameter,
# sorting the dead+live sample values, sorting the associated probability for each sample with the same indicies,
# making the cummulative distribution of the probability, and selecting relevent 1,2,3 sigma values from the cummulative distribution.
nfree=len(LNESTPARS[0,:])
#Prange = [None] * nfree
for bin in range(0,nfree):              # Loop over each fit parameter
    #sort parameters
    srt=np.argsort(LNESTPARS[:,bin])  # sort temperatures at one pressure level, index
    parameter=LNESTPARS[srt,bin]      # sort temperatures at one pressure level
    parameter_prob=prob[srt]            # probability sorted with increasing temperature
    #weighted average parmaeter value
    w_ave=np.sum(parameter*parameter_prob) 
    #print(w_ave)
    value=0.5
    absolute_val_array = np.abs(np.cumsum(parameter_prob) - value)
    smallest_difference_index = absolute_val_array.argmin()
    #print(parameter[smallest_difference_index])
    medpar=parameter[smallest_difference_index]
    #Calculate 1,2,3 sigma intervals from cummulative distribution of probability which has sorted by the parameter value
    smallest_difference_index=np.abs(np.cumsum(parameter_prob) - 0.5E0+0.68268E0/2.).argmin()
    el=(parameter[smallest_difference_index])
    smallest_difference_index=np.abs(np.cumsum(parameter_prob) - 0.5E0-0.68268E0/2.).argmin()
    eh=(parameter[smallest_difference_index])
    smallest_difference_index=np.abs(np.cumsum(parameter_prob) - 0.5E0+0.95449E0/2.).argmin()
    e2l=(parameter[smallest_difference_index])
    smallest_difference_index=np.abs(np.cumsum(parameter_prob) - 0.5E0-0.95449E0/2.).argmin()
    e2h=(parameter[smallest_difference_index])
    smallest_difference_index=np.abs(np.cumsum(parameter_prob) - 0.5E0+0.9973E0/2.).argmin()
    e3l=(parameter[smallest_difference_index])
    smallest_difference_index=np.abs(np.cumsum(parameter_prob) - 0.5E0-0.9973E0/2.).argmin()
    e3h=(parameter[smallest_difference_index])
    smallest_difference_index=np.abs(np.cumsum(parameter_prob) - 0.5E0+0.999936E0/2.).argmin()
    e4l=(parameter[smallest_difference_index])
    smallest_difference_index=np.abs(np.cumsum(parameter_prob) - 0.5E0-0.999936E0/2.).argmin()
    e4h=(parameter[smallest_difference_index])
    #save values to arrays
    med[bin]=medpar
    sig1l[bin]=el
    sig1h[bin]=eh
    sig2l[bin]=e2l
    sig2h[bin]=e2h
    sig3l[bin]=e3l
    sig3h[bin]=e3h
    sig4l[bin]=e4l
    sig4h[bin]=e4h
    
#-----------------------------------------------------------------------------------------
#plot

# Calculate and plot C/O ratio re-constructed from posterior
a=(tofit >18).nonzero()  #fit elements indicies
a=a[0]
if Chemistry == 3:
    Cindex=(ELEMENT_FIT_INDEX == 3).nonzero()
    Oindex=(ELEMENT_FIT_INDEX == 5).nonzero()
    if Cindex[0] >= 0 :
        if Oindex[0] >= 0 :
            C=(10**(LNESTPARS[:,a[Cindex[0]]])*3.162E-04)  #C VMR
            O=(10**(LNESTPARS[:,a[Oindex[0]]])*5.754E-04)  #O VMR
            C_O=C/O
            srt=np.argsort(C_O[:,0])  # sort temperatures at one pressure level, index
            parameter=np.squeeze(C_O[srt,0])      # sort temperatures at one pressure level
            parameter_prob=np.squeeze(prob[srt])            # probability sorted with increasing temperature
            value=0.5
            absolute_val_array = np.abs(np.cumsum(parameter_prob) - value)
            smallest_difference_index = absolute_val_array.argmin()
            medpar=parameter[smallest_difference_index]
            #Calculate 1,2,3 sigma intervals from cummulative distribution of probability which has sorted by the parameter value
            smallest_difference_index=np.abs(np.cumsum(parameter_prob) - 0.5E0+0.68268E0/2.).argmin()
            el=(parameter[smallest_difference_index])
            smallest_difference_index=np.abs(np.cumsum(parameter_prob) - 0.5E0-0.68268E0/2.).argmin()
            eh=(parameter[smallest_difference_index])
            smallest_difference_index=np.abs(np.cumsum(parameter_prob) - 0.5E0+0.999936E0/2.).argmin()
            e4l=(parameter[smallest_difference_index])
            smallest_difference_index=np.abs(np.cumsum(parameter_prob) - 0.5E0-0.999936E0/2.).argmin()
            e4h=(parameter[smallest_difference_index])
            print(medpar,e4l,e4h)
            #
            rcParams=rcParams_orig
            rcParams['font.size']= 16  # default 10
            matplotlib.pyplot.hist(parameter,weights=parameter_prob,histtype='stepfilled',bins=30,color='#4682b4',range=[0,1.25],alpha=0.75)
            xlabel('C/O ratio')
            ylabel('Probability')
            _ = xlim(np.round(e4l,decimals=1),np.round(e4h,decimals=1))
            #axes = np.array(figure.axes)
            #ax = axes
            axvline(0.56, color="b",linestyle='dotted')
            axvline(medpar, color="r")
            axvline(el, color="r",linestyle='--')
            axvline(eh, color="r",linestyle='--')
            pp = PdfPages(dirs+'arc_monitor_nested_corner_COratio.pdf')
            savefig(pp,format='pdf',bbox_inches='tight',pad_inches=0.05)
            pp.close()
            clf()
#-----------------------------------------------------------------------------------------
#COrner plot
# print(rcParams)
rcParams['text.usetex']=True
rcParams['text.latex.unicode']=True
rcParams['axes.linewidth'] = 2 # set the value globally
rcParams['lines.linewidth'] = 2
rcParams['font.size']= 16  # default 10
rcParams['font.weight']= 900  ## normal=400 bold=700, lighter,
rcParams['axes.labelsize']= 'large'
rcParams['axes.labelweight']= 900
rcParams['ytick.labelsize']= 'large'
rcParams['xtick.labelsize']= 'large'
rcParams['xtick.major.width'] = 2
rcParams['ytick.major.width'] = 2

rcParams["patch.force_edgecolor"] = True
rcParams['patch.linewidth'] = 2            #histogram line thickness
rcParams['patch.edgecolor'] = '#4682b4'
rcParams['patch.facecolor'] = '#4682b4'
 
#---------------------------------------------------------
# Default corner plot to show 4-sigma range
# Put 4-sigma range in python list to pass to corner.corner
range=[]                  # define blank list
i=0                       # define Counter because Python sucks and I can't just use item
for item in sig4l:        # Python's for loops are a "for each" loop, need to do this over the numpy array, not the python list piece of crap
    range.append((sig4l[i],sig4h[i]))
    i +=1

range[0]=(-1,sig4h[0])

figure=corner.corner(LNESTPARS,bins=12,weights=prob,range=range,show_titles=True,color='#4682b4',labels=parameternameslog)

# Extract the axes
axes = np.array(figure.axes).reshape((nfree, nfree))

# Add mean and sigma bars to historams
i=0                       # define Counter because Python sucks and I can't just use item
for item in sig4l:        # Python's for loops are a "for each" loop, need to do this over the numpy array, not the python list piece of crap
    ax = axes[i, i]
    ax.axvline(med[i], color="r")
    ax.axvline(sig1l[i], color="r",linestyle='--')
    ax.axvline(sig1h[i], color="r",linestyle='--')
    i +=1


#rcParams["legend.loc"] = 'upper right'
#legend(fontsize='x-small')
#show()
pp = PdfPages(dirs+'arc_monitor_nested_corner.pdf')
savefig(pp,format='pdf',bbox_inches='tight',pad_inches=0.05)
pp.close()
clf()

