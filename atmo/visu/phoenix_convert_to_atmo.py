from numpy import *
from pylab import *
from scipy.interpolate import *
import netCDF4 as nc
import sys

def write_phoenix(fin=None,fout=None):

	# Check input and ouput files specified
	if fin==None or fout==None:
		print('Error: you must specify both fin and fout')
		sys.exit()

	# Read file and use first and second columns
	ff = loadtxt(fin,usecols=(0,1))
	print("Completed loading text")
	#Sort wavelengths
	ff=ff[ff[:,0].argsort()]
	print("Completed sort waevelgnths")
	#Convert flux from file units to erg/s/cm**2/A
	#ff[:,1] = 10**(ff[:,1]+(-8.))

	#Convert wavelength from Angstrom to cm
	ff[:,0] = ff[:,0]/1.0E8
	ff[:,1] = ff[:,1]/1E-8

	#Convert erg/s/cm --> erg/s/st
	ff[:,1] = ff[:,1]/(4.*pi)

	wavelength = ff[:,0]
	flux = ff[:,1]

	#Construct variables to write
	fint = interp1d(1./(wavelength[::-1]),flux[::-1]*wavelength[::-1]**2,fill_value=0.,bounds_error=False)

	nfreq = 500000 # Number of points
	# Maximum wavenumber
	numax = 50000.
	# Wavenumber spacing 
	dnu = numax/nfreq

	# Compute new wavenumber grid
	wn = linspace(0.5*dnu,numax-0.5*dnu,nfreq)
	# Interpolate flux onto new wavenumber grid
	fn = fint(wn)

	#Write to netCDF
	nout = nc.Dataset(fout,'w',format='NETCDF3_CLASSIC')

	nout.createDimension('nfreq',nfreq)
	nnu=nout.createVariable('nu','f8',('nfreq',))
	nflux=nout.createVariable('hnu','f8',('nfreq',))

	nnu.units  = 'cm-1'
	nflux.units  = 'erg s-1 cm-1 ster-1'

	nnu[:]  = wn[:]
	nflux[:]  = fn[:]

	nout.close()