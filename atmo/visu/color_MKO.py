from constant          import *
from pylab             import *
from numpy             import *
from math              import *
from scipy.interpolate import *
from scipy.signal      import *
from scipy.optimize    import *
import matplotlib              as mpl
import matplotlib.font_manager as mpfm
import netCDF4 as nc

mpl.rc('text',usetex=True)
mpl.rc('font',family='serif',size=12)

prop = mpfm.FontProperties(size=12)

def calc_color(fin='spec1_t2000_g5.0_k8.ncdf', R = 0.09, d=10., dirfilter = '/Users/ptrembli/codes/atmo/visu/',verbose=True):  # t in K, radius R in Rsun

        dirmko     = dirfilter+'filter_MKO/'
        dirwise    = dirfilter+'response_wise/'
        dirhst     = dirfilter+'filter_HST/'
        dirspitzer = dirfilter+'filter_spitzer/'

        # filter files from http://irtfweb.ifa.hawaii.edu/~nsfcam/filters.html
        # and http://www.jach.hawaii.edu/UKIRT/instruments/wfcam/user_guide/description.html#Tab1.1 for WFCAM Y filter
        filj = loadtxt(dirmko+'nsfcam_jmk_trans.dat.txt')
        filk = loadtxt(dirmko+'nsfcam_kmk_trans.dat.txt')
        filks = loadtxt(dirmko+'nsfcam_ksmk_trans.dat.txt')
        filh = loadtxt(dirmko+'nsfcam_hmk_trans.dat.txt')
        fillp = loadtxt(dirmko+'nsfcam_lpmk_trans.dat.txt')
        fily = loadtxt(dirmko+'y-119.txt')

        filj_hst = loadtxt(dirhst+'HST_WFC3_IR.F125W.dat.txt')
        filh_hst = loadtxt(dirhst+'HST_WFC3_IR.F160W.dat.txt')

        #conversion from micrometer or nm to cm
        lambj =  filj[::-1,0]*1E-4
        lambk =  filk[::-1,0]*1E-4
        lambks =  filks[::-1,0]*1E-4
        lambh =  filh[::-1,0]*1E-4
        lamblp =  fillp[::-1,0]*1E-4
        lamby =  fily[::-1,0]*1E-7

        lambj_hst =  filj_hst[::-1,0]*1E-8
        lambh_hst =  filh_hst[::-1,0]*1E-8

        #transmission
        transj = np.maximum(filj[::-1,1],0.)/100.
        transk = np.maximum(filk[::-1,1],0.)/100.
        transks = np.maximum(filks[::-1,1],0.)/100.
        transh = np.maximum(filh[::-1,1],0.)/100.
        translp = np.maximum(fillp[::-1,1],0.)/100.
        transy = np.maximum(fily[::-1,1],0.)/100.

        transj_hst = np.maximum(filj_hst[::-1,1],0.)/filj_hst[::-1,1].max()
        transh_hst = np.maximum(filh_hst[::-1,1],0.)/filh_hst[::-1,1].max()

        #interpolation functions
        fintj = interp1d(lambj,transj,bounds_error = False, fill_value = 0.)
        fintk = interp1d(lambk,transk,bounds_error = False, fill_value = 0.)
        fintks = interp1d(lambks,transks,bounds_error = False, fill_value = 0.)
        finth = interp1d(lambh,transh,bounds_error = False, fill_value = 0.)
        fintlp = interp1d(lamblp,translp,bounds_error = False, fill_value = 0.)
        finty = interp1d(lamby,transy,bounds_error = False, fill_value = 0.)

        fintj_hst = interp1d(lambj_hst,transj_hst,bounds_error = False, fill_value = 0.)
        finth_hst = interp1d(lambh_hst,transh_hst,bounds_error = False, fill_value = 0.)

        fw1 = loadtxt(dirwise+'RSR-W1.txt')
        fw2 = loadtxt(dirwise+'RSR-W2.txt')
        fw3 = loadtxt(dirwise+'RSR-W3.txt')
        fw4 = loadtxt(dirwise+'RSR-W4.txt')

        fintw1 = interp1d(fw1[:,0]/1E4,abs(fw1[:,1])/sum(abs(fw1[:,1])),bounds_error = False, fill_value = 0.)
        fintw2 = interp1d(fw2[:,0]/1E4,abs(fw2[:,1])/sum(abs(fw2[:,1])),bounds_error = False, fill_value = 0.)
        fintw3 = interp1d(fw3[:,0]/1E4,abs(fw3[:,1])/sum(abs(fw3[:,1])),bounds_error = False, fill_value = 0.)
        fintw4 = interp1d(fw4[:,0]/1E4,abs(fw4[:,1])/sum(abs(fw4[:,1])),bounds_error = False, fill_value = 0.)


        #reference flux (Vega) in W/m**2/um
        #see http://irtfweb.ifa.hawaii.edu/IRrefdata/iwafdv.html and Hewett et al 2006, MNRAS 367, 454-468 (for the Y-band)
        Fveg_j = 3.01E-09
        Fveg_k = 4.00E-10
        Fveg_ks = 4.35E-10
        Fveg_h = 1.18E-09
        Fveg_lp = 5.31E-11
        Fveg_y = 5.71E-09

        #load spectrum
        spec = nc.Dataset(fin)
        lamb = np.zeros(spec.variables['nu'][:].size)
        lamb[:]= 1./spec.variables['nu'][:]            #wavelength in mc
        flamb  = spec.variables['fnu'][:]/lamb**2      #spectral flux in erg/s/cm**2/cm
        dlamb  = lamb[:-1]-lamb[1:]

        # mean spectral flux in bands in W/m**2/um
        flux_j = sum(fintj(lamb[:-1])*flamb[:-1]*lamb[:-1]*dlamb)/sum(fintj(lamb[:-1])*lamb[:-1]*dlamb)/1.e7
        flux_k = sum(fintk(lamb[:-1])*flamb[:-1]*lamb[:-1]*dlamb)/sum(fintk(lamb[:-1])*lamb[:-1]*dlamb)/1.e7
        flux_ks = sum(fintks(lamb[:-1])*flamb[:-1]*lamb[:-1]*dlamb)/sum(fintks(lamb[:-1])*lamb[:-1]*dlamb)/1.e7
        flux_h = sum(finth(lamb[:-1])*flamb[:-1]*lamb[:-1]*dlamb)/sum(finth(lamb[:-1])*lamb[:-1]*dlamb)/1.e7
        flux_lp = sum(fintlp(lamb[:-1])*flamb[:-1]*lamb[:-1]*dlamb)/sum(fintlp(lamb[:-1])*lamb[:-1]*dlamb)/1.e7
        flux_y = sum(finty(lamb[:-1])*flamb[:-1]*lamb[:-1]*dlamb)/sum(finty(lamb[:-1])*lamb[:-1]*dlamb)/1.e7


        flux_j_hst = sum(fintj_hst(lamb[:-1])*flamb[:-1]*lamb[:-1]*dlamb)/sum(fintj_hst(lamb[:-1])*lamb[:-1]*dlamb)/1.e7
        flux_h_hst = sum(finth_hst(lamb[:-1])*flamb[:-1]*lamb[:-1]*dlamb)/sum(finth_hst(lamb[:-1])*lamb[:-1]*dlamb)/1.e7

        # convert R to unit of 10pc
        R = R*696000/3.0857E14

        d = d/10.

        Nlamb = spec.variables['fnu'][:]/(hplanck*spec.variables['nu'][:])

        flux_w1 = sum(fintw1(lamb[:-1])*flamb[:-1]*lamb[:-1]*dlamb)/sum(fintw1(lamb[:-1])*lamb[:-1]*dlamb)/1.e7
        flux_w2 = sum(fintw2(lamb[:-1])*flamb[:-1]*lamb[:-1]*dlamb)/sum(fintw2(lamb[:-1])*lamb[:-1]*dlamb)/1.e7
        flux_w3 = sum(fintw3(lamb[:-1])*flamb[:-1]*lamb[:-1]*dlamb)/sum(fintw3(lamb[:-1])*lamb[:-1]*dlamb)/1.e7
        flux_w4 = sum(fintw4(lamb[:-1])*flamb[:-1]*lamb[:-1]*dlamb)/sum(fintw4(lamb[:-1])*lamb[:-1]*dlamb)/1.e7

	# compute colors
        jmag = -2.5*log10(flux_j/Fveg_j) - 5.*log10(R) + 5*log10(d)
        kmag = -2.5*log10(flux_k/Fveg_k) - 5.*log10(R) + 5*log10(d)
        ksmag = -2.5*log10(flux_ks/Fveg_ks) - 5.*log10(R) + 5*log10(d)
        hmag = -2.5*log10(flux_h/Fveg_h) - 5.*log10(R) + 5*log10(d)
        lpmag = -2.5*log10(flux_lp/Fveg_lp) - 5.*log10(R) + 5*log10(d)
        ymag = -2.5*log10(flux_y/Fveg_y) - 5.*log10(R) + 5*log10(d)

        jmag_hst = -2.5*log10(flux_j_hst/Fveg_j) - 5.*log10(R) + 5*log10(d)
        hmag_hst = -2.5*log10(flux_h_hst/Fveg_h) - 5.*log10(R) + 5*log10(d)


        w1mag =  - 2.5*log10(flux_w1/8.1787E-15/1E4)-5.*log10(R) + 5*log10(d)
        w2mag =  - 2.5*log10(flux_w2/2.4150E-15/1E4)-5.*log10(R) + 5*log10(d)
        w3mag =  - 2.5*log10(flux_w3/6.5151E-17/1E4)-5.*log10(R) + 5*log10(d)
        w4mag =  - 2.5*log10(flux_w4/5.0901E-18/1E4)-5.*log10(R) + 5*log10(d)

        if verbose:
                print('Y magnitude :', round(ymag,2))
                print('J magnitude :', round(jmag,2))
                print('H magnitude :', round(hmag,2))
                print('K magnitude :', round(kmag,2))
                print('Ks magnitude :', round(ksmag,2))
                print('Lp magnitude :', round(lpmag,2))
                print('H magnitude :', round(hmag,2))
                print('Y magnitude :', round(ymag,2))
                print('W1 mag: ',round(w1mag,2))
                print('W2 mag: ',round(w2mag,2))
                print('W3 mag: ',round(w3mag,2))
                print('W4 mag: ',round(w4mag,2))
                print('J mag HST : ',round(jmag_hst,2))
                print('H mag HST : ',round(hmag_hst,2))

                print('')
                print('Y - J =', ymag-jmag)
                print('J - H =', jmag-hmag)
                print('H - K =', hmag-kmag)
                print('J - K =', jmag-kmag)
                print('J - Ks =', jmag-ksmag)
                print('J - H (hst) =', jmag_hst-hmag_hst)
        return array((ymag,jmag,hmag,kmag,ksmag))

