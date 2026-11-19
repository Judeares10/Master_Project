from pylab import *
from netCDF4 import *
from constant import *

from pylab             import *
from numpy             import *
from math              import *
from scipy.interpolate import *
from scipy.signal      import *
from scipy.optimize    import *
from scipy.special     import sph_harm
from scipy.integrate   import *
import matplotlib              as mpl
import matplotlib.font_manager as mpfm
import mpl_toolkits.mplot3d as m3d

mpl.rc('text',usetex=False)    
mpl.rc('font',family='serif',size=12)

prop = mpfm.FontProperties(size=12)

def plot_mixing(nfig=1,feq='chem_out.ncdf',fneq='chem_neq_out.ncdf'):

    fdeq = Dataset(feq)
    fdneq = Dataset(fneq)

    figure(nfig)
    clf()
    sp=subplot(111)

    imol=[0,33,61,3,8,11,67]
    molname=['H2','H','He','CO','H2O','CH4','NH3']
    color = ['b','y','r','g','black','m','c']

    for i in range(len(imol)):
        
        plot(fdneq.variables['abundances'][imol[i],:],fdneq.variables['pressure'][:]/1E6,color=color[i],label=molname[i])
        
        plot(fdeq.variables['abundances'][imol[i],:],fdeq.variables['pressure'][:]/1E6,color=color[i],linestyle='--')

    sp.set_yscale('log')
    sp.set_xscale('log')
    
    ylim((1.0E3,1E-7))
    xlim((1E-12,1))

    xlabel('Abundances')
    ylabel('Pressure [bar]')

    legend(loc=0,prop=prop).draw_frame(0)
    show()
