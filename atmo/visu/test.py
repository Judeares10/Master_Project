from pylab import *
from netCDF4 import *
from constant import *
import read_atmo as ra

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

def test1(nfig=1,fOH_irr='OH_irr.dat',fOH_rev='OH_rev.dat'):

    figure(nfig,figsize=(6,4))
    clf()
    sp = subplot(111)

    ff_rev= loadtxt(fOH_rev)
    ff_irr= loadtxt(fOH_irr)

    plot(ff_rev[:,0],ff_rev[:,1],label='Fully reversed network')
    plot(ff_irr[:,0],ff_irr[:,1],label='With irreversible reactions')

    sp.set_xscale('log')
    xlabel('Integration time [s]')
    ylabel('Abundance of OH')
    title('P = 0.1 bar, T = 1800 K',fontsize=12)
    legend(loc=0,prop=prop).draw_frame(0)
    
    subplots_adjust(bottom=0.15)


def test2(nfig=2,feq='chem_eq.ncdf',fneq='chem_neq.ncdf'):

    fdeq = Dataset(feq)
    fdneq = Dataset(fneq)

    figure(nfig,figsize=(6,5))
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

    pmax=np.max(fdeq.variables['pressure'][:]/1E6)
    pmin=np.min(fdeq.variables['pressure'][:]/1E6)

    
    ylim((pmax,pmin))
    xlim((1E-12,1))

    xlabel('Abundances')
    ylabel('Pressure [bar]')

    legend(loc=0,prop=prop).draw_frame(0)

    subplots_adjust(bottom=0.12)
    show()

def test2_uv(nfig=2,feq='chem_eq.ncdf',fneq='chem_neq.ncdf'):

    fdeq = Dataset(feq)
    fdneq = Dataset(fneq)

    figure(nfig,figsize=(6,5))
    clf()
    sp=subplot(111)

    imol=[0,33,61,3,8,14,60]
    molname=['H2','H','He','CO','H2O','CO2','N2']
    color = ['b','y','r','g','black','m','c']

    for i in range(len(imol)):
        
        plot(fdneq.variables['abundances'][imol[i],:],fdneq.variables['pressure'][:]/1E3,color=color[i],label=molname[i])
        
        plot(fdeq.variables['abundances'][imol[i],:],fdeq.variables['pressure'][:]/1E3,color=color[i],linestyle='--')

    sp.set_yscale('log')
    sp.set_xscale('log')

    pmax=np.max(fdeq.variables['pressure'][:]/1E3)
    pmin=np.min(fdeq.variables['pressure'][:]/1E3)

    
    ylim((pmax,pmin))
    ylim((1E3,1E-5))
    xlim((1E-10,1))

    xlabel('Abundances')
    ylabel('Pressure [mbar]')

    legend(loc=0,prop=prop).draw_frame(0)

    subplots_adjust(bottom=0.12)
    show()

def test3(nfig=3):

    ff_ck = loadtxt('fluxheat_ck.dat')
    ff_lbl1 = loadtxt('fluxheat_lbl_1.dat')
    ff_lbl2 = loadtxt('fluxheat_lbl_2.dat')
    ff_lbl3 = loadtxt('fluxheat_lbl_3.dat')

    figure(nfig,figsize=(6,4))
    clf()
    sp = subplot(121)

    
    plot(ff_ck[:,1]/1E6,ff_ck[:,0]/1.0E6,label='Correlated k')
    plot((ff_lbl1[:,1]+ff_lbl2[:,1]+ff_lbl3[:,1])/1E6,ff_lbl1[:,0]/1.0E6,label='Line by line')
    
    
    
    legend(loc=0,prop=prop).draw_frame(0)
    pmin = ff_ck[:,0].min()/1.0E6
    pmax = ff_ck[:,0].max()/1.0E6

    sp.set_yscale('log')
    ylim((pmax,pmin))
    xlabel('Flux [kW/m$^2$]')
    ylabel('Pressure [bar]')

    show()

    sp = subplot(122)
    plot(ff_ck[:,2]/10.,ff_ck[:,0]/1.0E6,label='Correlated k')
    plot((ff_lbl1[:,2]+ff_lbl2[:,2]+ff_lbl3[:,2])/10.,ff_lbl1[:,0]/1.0E6,label='Line by line')

    legend(loc=0,prop=prop).draw_frame(0)
    pmin = ff_ck[:,0].min()/1.0E6
    pmax = ff_ck[:,0].max()/1.0E6

    sp.set_yscale('log')
    ylim((pmax,pmin))
    xlabel('Heating rate [W/m$^3$]')

    
    show()
    subplots_adjust(bottom=0.12)
    

def test4(nfig=4):

    figure(nfig,figsize=(6,5))
    clf()
    sp=subplot(111)

    ra.plot_pt(nfig=nfig,file='pt_t1200_g5.5_k0.ncdf',linestyle='-.',label='Teff=1200K, logg=5.5, logkzz=0.0',color='b')
    ra.plot_pt(nfig=nfig,file='pt_t1200_g5.5_k4.ncdf',linestyle='--',label='Teff=1200K, logg=5.5, logkzz=4.0',color='b')
    ra.plot_pt(nfig=nfig,file='pt_t1200_g5.5_k8.ncdf',linestyle='-',label='Teff=1200K, logg=5.5, logkzz=8.0',color='b')

    ra.plot_pt(nfig=nfig,file='pt_t1400_g5.5_k0.ncdf',linestyle='-.',label='Teff=1400K, logg=5.5, logkzz=0.0',color='r')
    ra.plot_pt(nfig=nfig,file='pt_t1400_g5.5_k4.ncdf',linestyle='--',label='Teff=1400K, logg=5.5, logkzz=4.0',color='r')
    ra.plot_pt(nfig=nfig,file='pt_t1400_g5.5_k8.ncdf',linestyle='-',label='Teff=1400K, logg=5.5, logkzz=8.0',color='r')
    

    subplots_adjust(bottom=0.12)
    show()

def test4_spectrum(nfig=41,fin='spec1_t2000_g5.0_k8.ncdf',clear=False,color='b'):

    figure(nfig,figsize=(12,5))
    if clear==True:
        clf()
    sp=subplot(111)

    ff1 = Dataset(fin)

    nlam = 1000

    f1int=interp1d(1./(ff1.variables['nu'][::-1]*1E-4),ff1.variables['fnu'][::-1]*ff1.variables['nu'][::-1]**2/1000./1E4,bounds_error=False)


    lam2 = linspace(0.801,5.,nlam)
  
    plot(lam2,f1int(lam2),color=color,label=fin)
  

    xlim((0.4,5))
   

    legend(prop=prop).draw_frame(0)
    xlabel(r'Wavelength [$\mu$m]')
    ylabel(r'Spectral flux [W/m$^2$/$\mu$m]')

    show()


def plot_phoenix(nfig=42,fin='lte012-5.0-0.0.BT-NextGen.7'):

    figure(nfig,figsize=(12,5))
    clf()
    ff = loadtxt(fin,usecols=(0,1))
    ff=ff[ff[:,0].argsort()]
    
    nlam = 1000
    fint=interp1d(ff[:,0]/1.0E4,10**ff[:,1]/1.0E7)
    lam2 = linspace(0.801,5.,nlam)

    plot(lam2,fint(lam2),label=fin)
        
    xlabel('Wavelength $\lambda$ [$\mu$m]')
    ylabel('Spectral flux F$_\lambda$ [W m$^{-2}$ $\mu$m$^{-1}$]')
    xlim((0,5))
    legend(prop=prop).draw_frame(0)
    show()

    


def test5(nfig=5):

    figure(nfig,figsize=(6,5))
    clf()
    sp=subplot(111)

    ra.plot_pt(nfig=nfig,file='pt_t1200_g3.0_k0.ncdf',linestyle='--',label='Teff=1200K, logg=3.0, logkzz=0.0',color='b')
    ra.plot_pt(nfig=nfig,file='pt_t1200_g3.0_k8.ncdf',linestyle='-.',label='Teff=1200K, logg=3.0, logkzz=8.0',color='b')
    ra.plot_pt(nfig=nfig,file='pt_t1200_g3.0_k11.ncdf',linestyle='-',label='Teff=1200K, logg=3.0, logkzz=11.0',color='b')

    subplots_adjust(bottom=0.12)
    show()


def test5_spectrum(nfig=51):

    figure(nfig,figsize=(12,5))
    clf()
    sp=subplot(111)

    ff1 = Dataset('spec1_t1200_g3.0_k0.ncdf')
    ff2 = Dataset('spec2_t1200_g3.0_k0.ncdf')

    ff3 = Dataset('spec1_t1200_g3.0_k11.ncdf')
    ff4 = Dataset('spec2_t1200_g3.0_k11.ncdf')

    nlam = 500

    dnu = ff1.variables['nu'][1]-ff1.variables['nu'][0]

    ftot1=sum(ff1.variables['fnu'][:])*dnu+sum(ff2.variables['fnu'][:])*dnu

    ftot2=sum(ff3.variables['fnu'][:])*dnu+sum(ff4.variables['fnu'][:])*dnu

    print ftot1/1000.,ftot2/1000.

    f1int=interp1d(1./(ff1.variables['nu'][::-1]*1E-4),ff1.variables['fnu'][::-1]*ff1.variables['nu'][::-1]**2/1000./1E4,bounds_error=False)

    f2int=interp1d(1./(ff2.variables['nu'][::-1]*1E-4),ff2.variables['fnu'][::-1]*ff2.variables['nu'][::-1]**2/1000./1E4,bounds_error=False)

    f3int=interp1d(1./(ff3.variables['nu'][::-1]*1E-4),ff3.variables['fnu'][::-1]*ff3.variables['nu'][::-1]**2/1000./1E4,bounds_error=False)

    f4int=interp1d(1./(ff4.variables['nu'][::-1]*1E-4),ff4.variables['fnu'][::-1]*ff4.variables['nu'][::-1]**2/1000./1E4,bounds_error=False)

    lam1 = linspace(0.401,0.799,nlam)
    lam2 = linspace(0.801,5.,nlam)
  
    plot(lam2,f1int(lam2),color='b',label='Teff=1200K, logg=3.0, logkzz=0')
    plot(lam1,f2int(lam1),color='b')
    plot(lam2,f3int(lam2),color='r',label='Teff=1200K, logg=3.0, logkzz=11')
    plot(lam1,f4int(lam1),color='r')

    xlim((0.4,5))
    ylim((0,300000))

    legend(prop=prop).draw_frame(0)
    xlabel(r'Wavelength [$\mu$m]')
    ylabel(r'Spectral flux [W/m$^2$/$\mu$m]')

    show()

def nak_spectrum(nfig=51):

    figure(nfig,figsize=(12,5))
    clf()
    sp=subplot(111)

    ff1 = Dataset('spec1_t1000_g5.5_NaK.ncdf')
    ff2 = Dataset('spec2_t1000_g5.5_NaK.ncdf')

    nlam = 2000

    dnu = ff1.variables['nu'][1]-ff1.variables['nu'][0]
    ftot1=sum(ff1.variables['fnu'][:])*dnu+sum(ff2.variables['fnu'][:])*dnu

    
    print ftot1/1000.


    dnu = ff1.variables['nu'][1]-ff1.variables['nu'][0]

    f1int=interp1d(1./(ff1.variables['nu'][::-1]*1E-4),ff1.variables['fnu'][::-1]*ff1.variables['nu'][::-1]**2/1000./1E4,bounds_error=False)

    f2int=interp1d(1./(ff2.variables['nu'][::-1]*1E-4),ff2.variables['fnu'][::-1]*ff2.variables['nu'][::-1]**2/1000./1E4,bounds_error=False)

    lam1 = linspace(0.401,0.799,nlam)
    lam2 = linspace(0.801,5.,nlam)
    
    plot(lam2,f1int(lam2),color='b',label='Teff=1000K, logg=5.5, Na-K')
    plot(lam1,f2int(lam1),color='b')


    fin = 'lte010-5.5-0.0.BT-NextGen.7'
    ff = loadtxt(fin,usecols=(0,1))
    ff=ff[ff[:,0].argsort()]
    
    nlam = 1000
    fint=interp1d(ff[:,0]/1.0E4,10**ff[:,1]/1.0E7)
    lam2 = linspace(0.601,2.,nlam)
    
    plot(lam2,fint(lam2),label=fin) 


    fin = 'lte10-5.5-0.0.BrownDwarf-mod.spec.7'
    ff = loadtxt(fin,usecols=(0,1))
    ff=ff[ff[:,0].argsort()]
    
    nlam = 1000
    fint=interp1d(ff[:,0]/1.0E4,10**(ff[:,1]-25.9007901434))
    lam2 = linspace(0.601,2.,nlam)
    
    plot(lam2,fint(lam2),label=fin)


    fin = 'lte010-5.5-0.0a+0.0.BT-Cond.7'
    ff = loadtxt(fin,usecols=(0,1))
    ff=ff[ff[:,0].argsort()]
    
    nlam = 1000
    fint=interp1d(ff[:,0]/1.0E4,10**ff[:,1]/1.0E7)
    lam2 = linspace(0.601,2.,nlam)
    
    plot(lam2,fint(lam2),label=fin) 


    fin = 'lte10-5.5-0.0.AMES-cond.7'
    ff = loadtxt(fin,usecols=(0,1))
    ff=ff[ff[:,0].argsort()]
    
    nlam = 1000
    fint=interp1d(ff[:,0]/1.0E4,10**(ff[:,1]-7.))
    lam2 = linspace(0.601,2.,nlam)
    
    plot(lam2,fint(lam2),label=fin)
    



    xlim((0.55,1.4))
    ylim((1E-1,1E6))

    sp.set_yscale('log')
    legend(prop=prop,loc=4).draw_frame(0)
    xlabel(r'Wavelength [$\mu$m]')
    ylabel(r'Spectral flux [W/m$^2$/$\mu$m]')

    show()




def plot_color(nfig=100):

    ff8= loadtxt('color_g5.0_k8.dat')
    ff0= loadtxt('color_g5.0_k0.dat')

    figure(nfig,figsize=(6,4))
    clf()

    plot(ff8[:,1],ff8[:,0],label='logg=5.0, logkzz=8')
    plot(ff0[:,1],ff0[:,0],label='logg=5.0, logkzz=0')
    ylim(ff8[:,0].max(),ff8[:,0].min())
    legend(prop=prop).draw_frame(0)
    xlabel('J - Ks')
    ylabel('M$_J$')
