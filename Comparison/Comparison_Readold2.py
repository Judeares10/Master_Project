#Copied all python module imports from read_atmo.py as needed when plotting ATMO data
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
	  # Load .dat file and returns data as np array and name of each coloum
    # Used to read PACT output file
    
    # Opens and reads file 
    with open(filename, "r") as f:
        lines = f.readlines()

    # Serchs for name of coloum by searching for ! 
    for i, line in enumerate(lines):
        if not line.startswith("!"):
            header_line = i
            break

    columns = lines[header_line].split()

		# Puts data from file into np array but avoids including name of data
    data = np.loadtxt(filename, skiprows=header_line + 1)

    return data, columns


def get_column(data, columns, name):
    # Exracts data from coloumn, given the name of the coloumn
    idx = columns.index(name)
    return data[:, idx]



def plot_pt(t='Rename me',nfig=1,ATMO='',PACT="",FREK="",VULC="",clear=False,label='',cond = '',color = '',color2 = '',linewidth = 2.,linestyle='-',marker=''):
	
	
    # OPTIONAL PARAMETERS :
    
    # t: plot title, default ''
    # nfig : index of the figure, default 1
    # ATMO : input ncdf ATMO file, default ''
    # PACT : input ncdf PACT file, default ''
    # FREK : input ncdf FRECKLL file, default ''
    # VULC : input ncdf VULCAN file, default ''
    # clear : erase previous curves if True, default False
    # label : curve label, default ''
    # color : curve color, python's default color choice if ='', default ''
    # cond : if ='H2O' or 'h2o', add H2O analytical condensation curve in red (from Lodders & Fegley 2001) ; if ='NH3' or 'nh3',
    #         add NH3 analytical condensation curve in magenta (same source) ; if ='both' or 'BOTH', add both condensation curves ; default '' (no condensation curve)
    # linewidth : line width of the curves, default 1
    # linestyle : default straight line '-'
    # marker : default none ('')
    
    
    # Sets up figure and sets figure size
    figure(nfig,figsize=(15,10))

		# If loop check if ATMO file included
    if not ATMO ==  "" :

				# Imports data from ATMO file
        ncdfFile = Dataset(ATMO,'r')
        vars = ncdfFile.variables
				
				# Extracts pressure and temperature from dataset
        tt = vars['temperature'][:]
        pp = vars['pressure'][:]

				# Converts to bar
        pp = pp/1.0E6

				# Plots ATMO PT profile
        plt.plot(tt, pp,label="ATMO")
    
    # If loop check if PACT file included
    if not PACT ==  "" :
        
        # Uses load_data_file to import dataset from PACT file
        data, columns = load_dat_file(PACT)

				# Uses get_column to extract pressure and temperature from dataset
        pressure = get_column(data, columns, "pressure[bar]")
        temperature = get_column(data, columns, "temperature[K]")
				# Plots PACT PT profile
        plt.plot(temperature, pressure,label="PACT")

    # If loop check if FRECKLLfile included
    if not FREK ==  "" :
			
        # Imports FRECKLL data from file and extracts pressure and temp 
        ds = Dataset(FREK, "r")
        pressure = ds.variables["pressure"][:]
        temperature = ds.variables["temperature"][:]
        # Plots PT profile for VULCAN includeds /1000 as need to convert from mil Bar
        plt.plot(temperature,pressure/1000,label="FRECKLE")
    
    # If loop check if VULCAN file included
    if not VULC == "" :
        # Uses pickle module to extract data from VULCAN file
        with open(VULC, 'rb') as handle:
            data = pickle.load(handle)
        # Plots PT profile for vulcan
        plt.plot(data['atm']['Tco'], data['atm']['pco']/1.e6, label="VULCAN")     
     
     # Sets y axis to log scale and sets limits for plot to match ATMO plots
    plt.gca().set_yscale('log')
    plt.ylim(1e3,1e-5) 
    plt.xlim(500,3000)
    # Plots sets labels for graph
    plt.xlabel('Temperature [K]',fontsize = 30)
    plt.ylabel('Pressure [bar]',fontsize = 30)
    # Adjusts size of font for axis and adds legend
    plt.tick_params(axis='both', which='major', labelsize=25)
    plt.legend(prop={'size': 30})
    # Saves figure and sets title
    plt.title(t,fontsize = 30)
    plt.savefig((t) +".png", dpi=200)
     
     
     
def plot_abun(t='Rename me',nfig=2,ATMO="",PACT="",FREK="",VULC="",aclear=True,linewidth = 1.,xlimits=[],atmo_imol=[1,4,9,12,34,62,68],linestyle='-'):
    
    # Sets up line styles and list for species names for later use
    speciesstr = []
    line_styles = ['-', '--', '-.', ':', (0, (5, 2)), (0, (3, 1, 1, 1))]
    
    # OPTIONAL PARAMETERS :
    
    # t: plot title, default ''
    # nfig : index of the figure, default 2
    # ATMO : input ncdf ATMO file, default ''
    # PACT : input ncdf PACT file, default ''
    # FREK : input ncdf FRECKLL file, default ''
    # VULC : input ncdf VULCAN file, default ''
    # aclear : erase previous curves if True, default False
    # linewidth : line width of the curves, default 1
    # xlimits : limits for x axis, defuault ''
    # atmo_imol : list of the molecules to plot
    
    
    # Sets up figure and sets figure size
    figure(nfig,figsize=(15,10))

    if not ATMO ==  "" :
        # Extracts data from ATMO file 
        ncdfFile = Dataset(ATMO,'r')
        vars = ncdfFile.variables
				# Extracts names of species
        nmol = len(ncdfFile.dimensions['nmol'])
        nlevel = len(ncdfFile.dimensions['nlevel'])

				# Formats atmo_imol list to be used for indexing
        if atmo_imol==[]:
            atmo_imol = array(range(nmol))
        else:
            atmo_imol = array(atmo_imol)-1
				
				# Extracts abundances, pressure and temp from data
        ab = vars['abundances'][:,:]
        pp = vars['pressure'][:]
        nm = vars['molname'][:,:]
        
        # Converts to Bars
        pp = pp/1.0E6
        
        # For loop to inlcuded all species in atmo_imol
        for i in range(atmo_imol.size):
            # Fix for label see week 4
            fixforlabel = nm[atmo_imol[i],:].tobytes()
            fixforlabel = str(fixforlabel).split("'")
            label = fixforlabel[1]
            # Adds the name of species to list to be used by other models
            speciesstr.append(label.strip())
            # Plots PT for ATMO
            plot(ab[atmo_imol[i],:],pp, label=label,linestyle=line_styles[i % len(line_styles)],color="tab:blue")
        # Sets y axis limits based on ATMO limits
        #ylim((pp.max(),pp.min()))

		# Checks if PACT file included
    if not PACT ==  "" :
        # Uses load_dat_file to extract data from PACT file
        data, columns = load_dat_file(PACT)
				# Extracts pressure data using get_column
        pressure = get_column(data, columns, "pressure[bar]")

				# For loop to include all species
        for spi in speciesstr:
            # Chekcs if species is in model
            if spi in columns:
                
                # Finds index posion of specie
                pos = speciesstr.index(spi)
                # Gets abundances for specie
                abundance = get_column(data, columns, spi)
                # Plots abundaces, PACT outputs is the log of abundance
                plot(10**(abundance), pressure,linestyle=line_styles[pos % len(line_styles)],color="tab:red")
            else:
                print(f"{spi} not found in file")

    if not FREK ==  "" :
        # Extracts dataset from FRECKLL file
        ds = Dataset(FREK, "r")
        # Extracts pressure data using get_column
        pressure = ds.variables["pressure"][:]
        # For loop to include all species
        for spi in speciesstr:
            # Gets abundances for specie
            abundance = ds.variables[spi][:]
            pos = speciesstr.index(spi)
            # Plots abundances also has /1000 as converts from mil bar
            plt.plot(abundance, pressure/1000,linestyle=line_styles[pos % len(line_styles)],color="tab:green")

    
    if not VULC ==  "" :

				# Uses pickle to extract data from VULCAN file
        with open(VULC, 'rb') as handle:
            data = pickle.load(handle)
				# Extract species list from vulcan for indexing
        vulcan_spec = data['variable']['species']
        # For loop to include all species
        for spi in speciesstr:
            pos = speciesstr.index(spi)
            # Plots abundances.
            plt.plot(data['variable']['ymix'][:,vulcan_spec.index(spi)], data['atm']['pco']/1.e6,linestyle=line_styles[pos % len(line_styles)],color="tab:orange")


		# Creates artistic legend that differentiates models by colour
    ATMO = mpatches.Patch(color='tab:blue', label='ATMO')
    VULCAN = mpatches.Patch(color='tab:orange', label='VULCAN')
    FRECKLE = mpatches.Patch(color='tab:green', label='FRECKLE')
    PACT = mpatches.Patch(color='tab:red', label='PACT')
    ax = plt.gca()
    legend1 = ax.legend(handles=[ATMO, VULCAN, FRECKLE, PACT], loc="upper right")
    ax.add_artist(legend1)




		# Sets labels, limits, log scales, and plots legend for differentiating species based on line style
    plt.ylim(1e3,1e-5) 
    plt.xlim(1e-11,1e1)
    xlabel('Abundances',fontsize = 25)
    ylabel("Pressure [bar]",fontsize = 25)
    plt.gca().set_yscale('log')
    plt.gca().set_xscale('log')
    leg = plt.legend()
    for line in leg.legendHandles:
        line.set_color('black')
    # Sets title and save fig
    plt.title(t,fontsize = 25)
    plt.savefig((t) +".png", dpi=200)



        



#plot_abun(t='Abundance Profilefor HD209 modelled by\n ATMO,FRECKLE,VULCAN And PACT',nfig=2,ATMO="chem_mix_pho_kinetics(kzz=1e9).ncdf",PACT="hd209458b_out_15.dat",FREK="HD209F.nc",VULC="HD209.vul",aclear=True,linewidth = 1.,xlimits=[],atmo_imol=[1,4])
