from pylab import *
from netCDF4 import *
from constant import *
# Jude was here
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

def load_dat_file(filename):
    """
    Load a .dat file and return:
    data : numpy array
    columns : list of column names
    """
    with open(filename, "r") as f:
        lines = f.readlines()

    # find header line (first line not starting with '!')
    for i, line in enumerate(lines):
        if not line.startswith("!"):
            header_line = i
            break

    columns = lines[header_line].split()

    data = np.loadtxt(filename, skiprows=header_line + 1)

    return data, columns


def get_column(data, columns, name):
    """
    Extract a column by name.
    """
    idx = columns.index(name)
    return data[:, idx]


def plot_pt(t='Rename me',nfig=1,file='t1000_g5_H2He.ncdf',file2="",tit='',file3="",clear=False,label='',cond = '',color = '',color2 = '',linewidth = 2.,linestyle='-',marker='',LA1="",LP="",LA2=""):
	
	
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
    	    semilogy(tt, pp,label=LA1,lw=linewidth,linestyle=linestyle,marker=marker)
    else :
    	    semilogy(tt, pp,label=LA1,color=color,lw=linewidth,linestyle=linestyle,marker=marker)
    
    if not file2 ==  "" :
        
        data, columns = load_dat_file(file2)

        pressure = get_column(data, columns, "pressure[bar]")
        temperature = get_column(data, columns, "temperature[K]")
        plt.plot(temperature, pressure,label=LP)
    
    if not file3 ==  "" :
        ncdfFile = Dataset(file3,'r')
        vars = ncdfFile.variables

        tt = vars['temperature'][:]
        pp = vars['pressure'][:]

        #conversion in bar
        pp = pp/1.0E6

        if clear:
            clf()
    
        if (cond == 'h2o' or cond == 'H2O' or cond == 'both' or cond == 'BOTH') :
    	    semilogy(10000./(38.84-3.83*numpy.log10(pp)-3.93*0.-0.20*0.*numpy.log10(pp)),pp,label = 'H2O condensation line',color='red',lw=linewidth)
        if (cond == 'nh3' or cond == 'NH3' or cond == 'both' or cond == 'BOTH') :
    	    semilogy(10000./(68.02-6.31*numpy.log10(pp)-6.19*0.),pp,label = 'NH3 condensation line',color='magenta',lw=linewidth)
   	    
        if (color=='') :
    	    semilogy(tt, pp,label=LA2,lw=linewidth,linestyle=linestyle,marker=marker)
        else :
    	    semilogy(tt, pp,label=LA2,color=color,lw=linewidth,linestyle=linestyle,marker=marker)

        

    	    
    ylim((pp.max(),pp.min()))
    xlim(500,3000)
    xlabel('Temperature [K]',fontsize = 30)
    ylabel('Pressure [bar]',fontsize = 30)

    #custom_lines = [Line2D([0], [0]),Line2D([0], [0],color="orange")]
    #ax = plt.gca()
    #leg1 = legend(custom_lines, [l1, l2])
    #ax.add_artist(leg1)
    plt.tick_params(axis='both', which='major', labelsize=25)
    plt.legend(prop={'size': 30})
    plt.title(t,fontsize = 30)
    plt.savefig((t) +".png", dpi=200)