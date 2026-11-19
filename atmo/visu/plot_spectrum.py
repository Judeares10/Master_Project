from pylab import *
from netCDF4 import *
import numpy as np
import matplotlib.ticker
from scipy.interpolate import *


def plot_trans(nfig=1,fin='',tit='',label='',color='r',linewidth=2):

	spectrum = Dataset(fin)
	nu = spectrum.variables['nu'][:]
	rprs = spectrum.variables['transit_radius'][:]
	
	figure(nfig,figsize=(18,10))
	
	nu = 1./nu
	nu = nu/100.
	nu = nu*1E6
	
	#Convert to 'transit depth'
	#rprs = rprs**2.

	semilogx(nu[:],rprs[:],label=label,color=color,linewidth=linewidth)
	
	xlabel('Wavelength $\mu$m',fontsize=35)
	ylabel('R$_p$/R$_s$',fontsize=35)
	xlim(0.3,30)
	tick_params(axis='x',labelsize='30')
	tick_params(axis='y',labelsize='30')
	
	
	#Set xaxis ticker format
	ax = gca()
	majorLocator = FixedLocator([0.4,0.6,1.0,2,4,8,16,24])
	ax.xaxis.set_major_locator(majorLocator)
	ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
	
	title(tit,fontsize=18)
	
	legend(loc=0,fontsize=30).draw_frame(0)
	plt.savefig("plottr.png", dpi=200)

	show()
	
def plot_emiss(nfig=1,fin='',tit='',label='',color='r',linewidth=1.0):
	
	spectrum = Dataset(fin)
	nu = spectrum.variables['nu'][:]
	fnup = spectrum.variables['fnu'][:]
	fnus = spectrum.variables['fnu_star'][:]
	rp = spectrum.variables['R_planet_TOA'][:]
	rs = spectrum.variables['R_star'][:]
	
	figure(nfig,figsize=(18,10))

	#Convert to wavelength and SI
	nu = 1./nu
	nu = nu/100.
		
	#Convert to um
	nu = nu*1E6
	
	#Flux ratio
	fpfs = fnup/fnus
	
	fpfs = fpfs*(rp**2./rs**2.)

	semilogx(nu[:],fpfs[:],label=label,color=color,linewidth=linewidth)
	
	xlabel('Wavelength $\mu$m',fontsize=35)
	ylabel('F$_p$/F$_s$',fontsize=35)
	tick_params(axis='x',labelsize='30')
	tick_params(axis='y',labelsize='30')
	xlim(0.3,30)
	
	ylim(1E-6,0.0040)
	
	#Set xaxis ticker format
	ax = gca()
	majorLocator = FixedLocator([0.4,0.6,1.0,2,4,8,16,24])
	ax.xaxis.set_major_locator(majorLocator)
	ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())

	title(tit,fontsize=30)

	
	legend(loc=0,fontsize=30).draw_frame(0)
	plt.savefig("plotem.png", dpi=200)

	show()



