# Packages
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib import cm
import re
from pathlib import Path

# Paths
file_paths = ['05_16_25_100mM_DA_Vesicles_10mM_MethylOrange_Trial5.csv']
output_paths = ['']

# General Parameters
exp_header_regex = r'.*_t(\d*)min'
num_conditions = 2
condition_names = ['Encapsulated', 'Not Encapsulated']
epsA464 = 21.6
show_figures = True
save_figures = False
img_dpi = 200 # want more for higher resolution. Slower when more.
fig_size_in_inches = (10, 10) # (width, height)

if save_figures: 
    assert len(file_paths)==len(output_paths), "\nTo save the figures, make sure an output folder path is given for each file path.\nAborted."
##########################################################################################

# Script
for fi, fp in enumerate(file_paths):
    data = dict()
    
    # Get the data
    with open(fp, 'r') as f:
        rows = f.readlines()
        experiment_headers = rows[0].strip('\n').split(',')
        header_inds, header_names = zip(*[(i, h) for i, h in enumerate(experiment_headers) if re.match(exp_header_regex, h) ])
        for row in rows[2:]:
            row = row.strip('\n').split(',')
            if row[0]:
                row_data = [float(row[i+1]) for i in header_inds if row[i]]
                data[float(row[0])] = row_data
            else: break
    
    # Plot absorption curves
    timesteps = np.unique(np.array([float(re.findall(exp_header_regex, h)[0]) for h in header_names]))
    timestep_unit = None

    fig, axs = plt.subplots(*(num_conditions, 1), 
                                   figsize=fig_size_in_inches, 
                                   sharex=False, 
                                   layout='tight')
    fig.suptitle(fp+'\n', fontsize=16)

    for i, t in enumerate(timesteps):
        for cond_num in range(num_conditions):
            axs[cond_num].plot(data.keys(), [data[wave][cond_num::2][i] for wave in data], color=cm.turbo(i/len(timesteps)), label=f'{t} min')

    # Set Plot Settings
    for cond_num in range(num_conditions):    
        axs[cond_num].set_title(condition_names[cond_num], fontsize=18)
        axs[cond_num].set_xlabel('Wavelength (nm)',fontsize=16)
        axs[cond_num].set_ylabel('Absorbance', fontsize=16)
        #axs[cond_num].title('10 mM Methyl Orange Encapsulated in Oleic Acid Vesicles Trial 2')
        axs[cond_num].legend(loc='best')
        axs[cond_num].tick_params(which='both', labelsize=14, pad=5)
        axs[cond_num].set_xlim(200,600)
    
    if show_figures:
        plt.show()
    if save_figures:
        plt.savefig(Path(output_paths[fi]) / 'absorbances.png', dpi=img_dpi)

    
    # Plot kinetics data
    def linear(x, a, b):
        return a*x + b
    
    key_A464 = list(data.keys())[list(map(round, data.keys())).index(464)] 
    
    fig, axs = plt.subplots(*(2, 1),
                            layout='constrained', 
                            figsize=fig_size_in_inches,
                            sharex=False)
    fig.suptitle(fp+'\n', fontsize=16)

    for cond_num in range(num_conditions):
        A464 = np.array(data[key_A464][cond_num::2])
        concMO = A464 / epsA464
        lnconcMO = np.log(concMO)
    
        lnMO = np.log(A464)
    
        popt, pcov = curve_fit(linear, timesteps, lnconcMO)
        #print(popt)
        fitconcMO = linear(timesteps, *popt)
        #print(fitconcMO)
        #print(concMO*1000)
        axs[cond_num].plot(timesteps, lnconcMO, color=cm.turbo(0.0), marker='o', ls='none', label='lnconcMO')
        axs[cond_num].plot(timesteps, fitconcMO, color=cm.turbo(0.0), label=f'y = {popt[0]: .4f}x + {popt[1]: .4f}')
        axs[cond_num].set_xlabel('Time (min)', fontsize=16)
        axs[cond_num].set_ylabel('ln[Methyl Orange]', fontsize=16)
        axs[cond_num].tick_params(which='both', labelsize=14, pad=5)
        #plt.title('Encapsulated Methyl Orange UV Decomposition')
        axs[cond_num].set_title(f"{condition_names[cond_num]}", fontsize=18)
        axs[cond_num].legend(loc='best')
        
    if show_figures:
        plt.show()
    if save_figures:
        plt.savefig(Path(output_paths[fi]) / 'kinetics.png', dpi=img_dpi)