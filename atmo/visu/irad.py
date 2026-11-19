
from constant          import *
from pylab             import *
from numpy             import *
from math              import *
from scipy.interpolate import *
from scipy.signal      import *
from scipy.optimize    import *
import matplotlib              as mpl
import matplotlib.font_manager as mpfm
import read_atmo as ra
import netCDF4 as nc

mpl.rc('text',usetex=True)    
mpl.rc('font',family='serif',size=12)

prop = mpfm.FontProperties(size=12)

def write_irad(fin='hd209458pr5.500000.txt',fout='kurucz_hd209.ncdf'):

    ff=loadtxt(fin)
    lam  = ff[:,0]
    flux = ff[:,1]

  #  lam  = loadtxt(flam)
  #  sflux = loadtxt(fflux,dtype='|S90')
  #  flux = np.zeros((sflux.size,8))

    # for i in range(sflux.size):
    #     for j in range(8):
    #         flux[i,j] = float(sflux[i][j*10:(j+1)*10])
            
    # lam = lam.ravel()
    # flux = flux.ravel()

    #convert in cgs

    lam = lam*1.0E-7
    flux = flux/1.0E-7

    fint = interp1d(1./(lam[::-1]),flux[::-1]*lam[::-1]**2,fill_value=0.,bounds_error=False)

    nfreq = 500000
    numax = 50000.
    dnu = numax/nfreq
    print(dnu)

    wn = linspace(0.5*dnu,numax-0.5*dnu,nfreq)
    fn = fint(wn)

    print(sum(fn)*dnu*4*pi)

    nout = nc.Dataset(fout,'w',format='NETCDF3_CLASSIC')

    nout.createDimension('nfreq',nfreq)
    nnu=nout.createVariable('nu','f8',('nfreq',))
    nflux=nout.createVariable('hnu','f8',('nfreq',))

    nnu.units  = 'cm-1'
    nflux.units  = 'erg s-1 cm-1 ster-1'
    
    nnu[:]  = wn[:]
    nflux[:]  = fn[:]

    nout.close()
