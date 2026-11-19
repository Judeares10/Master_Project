from pylab import *
from netCDF4 import *
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
import pickle
import matplotlib.patches as mpatches



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



def plot_pt(t='Rename me',nfig=1,ATMO='',PACT="",FREK="",VULC="",clear=False,label='',cond = '',color = '',color2 = '',linewidth = 2.,linestyle='-',marker='',LA1="",LP="",LA2=""):
	
	
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
    figure(nfig,figsize=(15,10))

    if not ATMO ==  "" :

        ncdfFile = Dataset(ATMO,'r')
        vars = ncdfFile.variables

        tt = vars['temperature'][:]
        pp = vars['pressure'][:]

        pp = pp/1.0E6


        plt.plot(tt, pp,label="ATMO")
    
    if not PACT ==  "" :
        
        data, columns = load_dat_file(PACT)

        pressure = get_column(data, columns, "pressure[bar]")
        temperature = get_column(data, columns, "temperature[K]")
        plt.plot(temperature, pressure,label="PACT")

    if not FREK ==  "" :
        ds = Dataset(FREK, "r")
        pressure = ds.variables["pressure"][:]
        temperature = ds.variables["temperature"][:]
        
        plt.plot(temperature,pressure/1000,label="FRECKLE")

    if not VULC == "" :
        with open(VULC, 'rb') as handle:
            data = pickle.load(handle)
        plt.plot(data['atm']['Tco'], data['atm']['pco']/1.e6, label="VULCAN")

    plt.gca().set_yscale('log')
    plt.ylim(1e3,1e-5) 
    plt.xlim(500,3000)
    plt.xlabel('Temperature [K]',fontsize = 30)
    plt.ylabel('Pressure [bar]',fontsize = 30)
    plt.tick_params(axis='both', which='major', labelsize=25)
    plt.legend(prop={'size': 30})
    plt.title(t,fontsize = 30)
    plt.savefig((t) +".png", dpi=200)


#plot_pt(t='HD209 PT Profile modelled by\n ATMO,FRECKLE,VULCAN And PACT',nfig=1,ATMO="/home/jh1639/Comparison/chem_mix_pho_kinetics(kzz=1e9).ncdf",PACT="hd209458b_out_15.dat",FREK="HD209F.nc",VULC="HD209.vul",clear=False,label='',cond = '',color = '',color2 = '',linewidth = 2.,linestyle='-',marker='')


def plot_abun(t='Rename me',nfig=2,ATMO="",PACT="",FREK="",VULC="",aclear=True,linewidth = 1.,xlimits=[],atmo_imol=[1,4,9,12,34,62,68],linestyle='-'):
    speciesstr = []
    line_styles = ['-', '--', '-.', ':', (0, (5, 2)), (0, (3, 1, 1, 1))]
    
    
    
    
    figure(nfig,figsize=(15,10))

    if not ATMO ==  "" :
        ncdfFile = Dataset(ATMO,'r')
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
        pp = pp/1.0E6
        for i in range(atmo_imol.size):
            fixforlabel = nm[atmo_imol[i],:].tobytes()
            fixforlabel = str(fixforlabel).split("'")
            label = fixforlabel[1]
            speciesstr.append(label.strip())
            plot(ab[atmo_imol[i],:],pp, label=label,linestyle=line_styles[i % len(line_styles)],color="tab:blue")
        ylim((pp.max(),pp.min()))


    if not PACT ==  "" :
        data, columns = load_dat_file(PACT)

        pressure = get_column(data, columns, "pressure[bar]")


        for spi in speciesstr:
            if spi in columns:
                pos = speciesstr.index(spi)
                abundance = get_column(data, columns, spi)
                plot(10**(abundance), pressure,linestyle=line_styles[pos % len(line_styles)],color="tab:red")
            else:
                print(f"{spi} not found in file")

    if not FREK ==  "" :
        ds = Dataset(FREK, "r")
        pressure = ds.variables["pressure"][:]
        for spi in speciesstr:
            print(spi)
            abundance = ds.variables[spi][:]
            pos = speciesstr.index(spi)
            plt.plot(abundance, pressure/1000,linestyle=line_styles[pos % len(line_styles)],color="tab:green")

    if not VULC ==  "" :

        with open(VULC, 'rb') as handle:
            data = pickle.load(handle)

        vulcan_spec = data['variable']['species']
        for spi in speciesstr:
            pos = speciesstr.index(spi)
            plt.plot(data['variable']['ymix'][:,vulcan_spec.index(spi)], data['atm']['pco']/1.e6,linestyle=line_styles[pos % len(line_styles)],color="tab:orange")



    ATMO = mpatches.Patch(color='tab:blue', label='ATMO')
    VULCAN = mpatches.Patch(color='tab:orange', label='VULCAN')
    FRECKLE = mpatches.Patch(color='tab:green', label='FRECKLE')
    PACT = mpatches.Patch(color='tab:red', label='PACT')

    ax = plt.gca()

    legend1 = ax.legend(handles=[ATMO, VULCAN, FRECKLE, PACT], loc="upper right")
    ax.add_artist(legend1)





    plt.ylim(1e3,1e-5) 
    plt.xlim(1e-11,1e1)
    xlabel('Abundances',fontsize = 25)
    ylabel("Pressure [bar]",fontsize = 25)
    plt.gca().set_yscale('log')
    plt.gca().set_xscale('log')
    leg = plt.legend()
    for line in leg.legendHandles:
        line.set_color('black')
    plt.title(t,fontsize = 25)
    plt.savefig((t) +".png", dpi=200)



        



#plot_abun(t='Abundance Profile',nfig=2,ATMO="chem_mix_pho_kinetics(kzz=1e9).ncdf",PACT="hd209458b_out_15.dat",FREK="HD209F.nc",VULC="HD209.vul",aclear=True,linewidth = 1.,xlimits=[],atmo_imol=[1,4])
