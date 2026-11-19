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

def plot_pt(nfig=1,file='t1000_g5_H2He.ncdf',tit='',clear=False,label='',cond = '',color = '',linewidth = 1.,linestyle='-',marker=''):
	
	
    # OPTIONAL PARAMETERS :
    
    # nfig : index of the figure, default 1
    # file : input ncdf P-T profile file, default 't1000_g5_H2He.ncdf'
    # tit : plot title, default ''
    # clear : erase previous curves if True, default False
    # label : curve label, default ''
    # color : curve color, python's default color choice if ='', default ''
    # cond : if ='H2O' or 'h2o', add H2O analytical condensation curve in red (from Lodders & Fegley 2001) ; if ='NH3' or 'nh3',
    #         add NH3 analytical condensation curve in magenta (same source) ; if ='both' or 'BOTH', add both condensation curves ; default '' (no condensation curve)
    # linewidth : line width of the curves, default 1
    # linestyle : default straight line '-'
    # marker : default none ('')


    ncdfFile = Dataset(file,'r')
    vars = ncdfFile.variables

    tt = vars['temperature'][:]
    pp = vars['pressure'][:]

    #conversion in bar
    pp = pp/1.0E6

    figure(nfig,figsize=(15,10))
    if clear:
        clf()
    
    if (cond == 'h2o' or cond == 'H2O' or cond == 'both' or cond == 'BOTH') :
    	    semilogy(10000./(38.84-3.83*numpy.log10(pp)-3.93*0.-0.20*0.*numpy.log10(pp)),pp,label = 'H2O condensation line',color='red',lw=linewidth)
    if (cond == 'nh3' or cond == 'NH3' or cond == 'both' or cond == 'BOTH') :
    	    semilogy(10000./(68.02-6.31*numpy.log10(pp)-6.19*0.),pp,label = 'NH3 condensation line',color='magenta',lw=linewidth)
   	    
    if (color=='') :
    	    semilogy(tt, pp,label=label,lw=linewidth,linestyle=linestyle,marker=marker)
    else :
    	    semilogy(tt, pp,label=label,color=color,lw=linewidth,linestyle=linestyle,marker=marker)
    	    
    title(tit,fontsize=12)
    ylim((pp.max(),pp.min()))
    xlabel('Temperature [K]')
    ylabel('Pressure [bar]')

    if not label=='':
        legend(prop=prop).draw_frame(0)
    show()



def plot_abundances(nfig=2,file='chem_out.ncdf',tit='',clear=True,linewidth = 1.,xlimits=[],atmo_imol=[1,4,9,12,34,62,68],linestyle='-'):
	
	
    # OPTIONAL PARAMETERS :
    
    # nfig : index of the figure, default 2
    # file : input ncdf abundances file, default 'chem_out.ncdf'
    # tit : plot title, default ''
    # clear : erase previous curves if True, default False
    # linewidth : line width of the curves, default 1
    # linestyle : line widthstyle of the curves, default '-'
    # atmo_imol : list of the molecules to plot


    ncdfFile = Dataset(file,'r')
    vars = ncdfFile.variables

    nmol = len(ncdfFile.dimensions['nmol'])
    nlevel = len(ncdfFile.dimensions['nlevel'])

    if atmo_imol==[]:
        atmo_imol = array(range(nmol))
    else:
        atmo_imol = array(atmo_imol)-1

    ab = vars['abundances'][:,:]
    pp = vars['pressure'][:]
    nm = vars['molname'][:,:]

    #conversion in bar
    pp = pp/1.0E6

    figure(nfig,figsize=(15,10))

    if clear:
        clf()
    sp=subplot(111)
    for i in range(atmo_imol.size):
        label = nm[atmo_imol[i],:].tostring()
        plot(ab[atmo_imol[i],:], pp,label=label,lw=linewidth)

    sp.set_yscale('log')
    sp.set_xscale('log')
    title(tit,fontsize=12)
    ylim((pp.max(),pp.min()))
    if (xlimits): # Check that list is not empty
      xlim(xlimits)
    xlabel('Abundances')
    ylabel('Pressure [bar]')
    legend(prop=prop).draw_frame(0)

    show()
 
def plot_spectrum(nfig = 3, file = 't700_spectrum.ncdf', tit = '',clear = False,label = '',coeff = 1.,color = '',temp = False,nbin = 1,linewidth=1.):
	
	
    # OPTIONAL PARAMETERS :
    
    # nfig : index of the figure, default 3
    # file : input ncdf spectrum file, default 't700_spectrum.ncdf'
    # tit : plot title, default ''
    # clear : erase previous curves if True, default False
    # label : curve label, default ''
    # color : curve color, python's default color choice if ='', default ''
    # temp : if True plot the brightness temperature instead of the spectral flux, default False
    # coeff : only when temp = False, plot the spectral flux scaled by coeff, default 1.
    # nbin : if > 1, calculate nbin-points arithmetic means on the frequencies and the fluxes or the temperatures to decrease the number of points to plot
    # linewidth : line width of the curves, default 1.


    ncdfFile = Dataset(file,'r')
    vars = ncdfFile.variables

    nu = vars['nu'][:]
    fnu = vars['fnu'][:]
    nfreq = size(nu)


    if (not temp) :
       
       if (nbin == 1) :
          #conversion
          lamb = 1./(nu * 1.e-4)
          flambda = fnu/(lamb**2)*1.e4 /1000.

          figure(nfig,figsize=(15,10))
          if clear:
              clf()
          if color == '': plot(lamb, flambda*coeff,label=label,lw=linewidth)  
          else :          plot(lamb, flambda*coeff,label=label,color=color,lw=linewidth)
          
       else :
  	  # arithmetic means
          nur = numpy.zeros((nfreq/nbin))
          fnur = numpy.zeros((nfreq/nbin))
          for i in range(1,nfreq/nbin):
       	       nur[i] = nu[nbin*i]
       	       fnur[i] = sum(fnu[nbin*i-nbin/2:nbin*i+nbin/2])/nbin
          nur[0] = nu[0]
          fnur[0] = sum(fnu[0:nbin/2 -1])/(nbin/2)
          
          #conversion
          lambr = 1./(nur*1.e-4)
          flambdar = fnur/(lambr**2)*1.e4 /1000.
          
          figure(nfig,figsize=(15,10))
          if clear:
              clf()
          if color == '': plot(lambr, flambdar*coeff,label=label,lw=linewidth)  
          else :          plot(lambr, flambdar*coeff,label=label,color=color,lw=linewidth)
          
       title(tit,fontsize=12)
    
       xlabel('Wavelength [$\mu$m]')
       #if coeff == 1. : ylabel('Spectral flux [W.m-2.$\mu$m-1]')
       #else :           ylabel('Spectral flux [arbitrary unit]')
       ylabel('Spectral flux [W.m-2.$\mu$m-1]')

       
    else :
    	    
       #conversion into brightness temperature
       tb = hplanck*c*nu/kb/numpy.log(1. + 2.*pi*hplanck*c2*nu**3/fnu)
       
       if (nbin == 1) :

          #conversion
          lamb = 1./(nu * 1.e-4)
          
          figure(nfig,figsize=(15,10))
          if clear:
             clf()
          if color == '': plot(lamb, tb,label=label,lw=linewidth)  
          else :          plot(lamb, tb,label=label,color=color,lw=linewidth)
          
       else :
          
          #arithmetic means
          nur = numpy.zeros((nfreq/nbin))
          tbr = numpy.zeros((nfreq/nbin))
          for i in range(1,nfreq/nbin):
       	       nur[i] = nu[nbin*i]
       	       tbr[i] = sum(tb[nbin*i-nbin/2:nbin*i+nbin/2])/nbin
          nur[0] = nu[0]
          tbr[0] = sum(tb[0:nbin/2 - 1])/(nbin/2)
          
          #conversion
          lambr = 1./(nur*1.e-4)
        
          figure(nfig,figsize=(15,10))
          if clear:
             clf()
          if color == '': plot(lambr, tbr,label=label,lw=linewidth)  
          else :          plot(lambr, tbr,label=label,color=color,lw=linewidth)
       
       title(tit,fontsize=12)
    
       xlabel('Wavelength [$\mu$m]')
       ylabel('Brightness temperature [K]')
       
 
    if not label=='':
        legend(prop=prop).draw_frame(0)
    show()



def plot_trans_spec(nfig = 4, file = 'transmision_spectrum.ncdf', tit = '',clear = False,label = '',color = '',xunit='wn',linewidth = 1.):


    # OPTIONAL PARAMETERS :
    
    # nfig : index of the figure, default 4
    # file : input ncdf effective radii file, default 't700_radii.ncdf'
    # tit : plot title, default ''
    # clear : erase previous curves if True, default False
    # label : curve label, default ''
    # color : curve color, python's default color choice if ='', default ''
    # nbin : if > 1, calculate nbin-points arithmetic means on the frequencies and the effective radii to decrease the number of points to plot
    # linewidth : line width of the curves, default 1.
    
    
    ncdfFile = Dataset(file,'r')
    vars = ncdfFile.variables

    nu = vars['nu'][:]
    Rp = vars['transit_radius'][:]
    
    figure(nfig,figsize=(15,10))
    title(tit)

    if (xunit == 'wn'):

      if clear:
        clf()
      if (color == ''):
        plot(nu,Rp,label=label,lw=linewidth)  
      else:
        plot(nu,Rp,label=label,color=color,lw=linewidth)
      
      xlabel('Wavenumber [cm$^{-1}$]')
      ylabel('Effective radius [$R_*$]')
    
    elif (xunit == 'wl'):
      print 
      lamb = 1./(nu*1e-4)
      figure(nfig,figsize=(15,10))
      if clear:
        clf()
      if (color == ''):
        semilogx(lamb,Rp,label=label,lw=linewidth)
      else:
        semilogx(lamb,Rp,label=label,color=color,lw=linewidth)
        
      xlabel('Wavelength [m]')
      ylabel('Effective radius [$R_*$]')
    else:
      print('Unit of x-axis not recognised')
      return

    if not label=='':
        legend(prop=prop).draw_frame(0)
    
    show()
