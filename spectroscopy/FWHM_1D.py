# imports
import numpy as np
from scipy.optimize import curve_fit
import WrightTools as wt
import matplotlib.pyplot as plt
from matplotlib.transforms import TransformedBbox, Bbox
from scipy import signal as s
from scipy.stats import norm
from scipy.special import erfc, erf
from scipy.signal import find_peaks, peak_widths
import sys
import argparse
import os
import math
import easygui

# fitting functions
def gauss_curve(x, x0, a, sigma):
    return a * np.exp(-(x-x0)**2/(2*sigma**2))

def skew2(x, sigmag, mu, alpha, c, a):
    normpdf = (1 / (sigmag * np.sqrt(2 * math.pi))) * np.exp(-(np.power((x - mu), 2) / (2 * np.power(sigmag, 2))))
    normcdf = (0.5 * (1 + erf((alpha * ((x - mu) / sigmag)) / (np.sqrt(2)))))
    return 2 * a * normpdf * normcdf + c, np.max(normpdf)

def skew(x, mu, a, sigmag, alpha, c):
    return skew2(x, sigmag, mu, alpha, c, a)[0]

parser = argparse.ArgumentParser(prog='FWHM_1D', description='')
parser.add_argument('-d', '--default_directory')
parser.add_argument('-s', '--set_default') #need to store this in a sep file, not environment variable
args = parser.parse_args()

# file inputs and setup
try: 
    assert args.default_directory
    default_path = args.default_directory
except: #implement saving to other file and calling that directory, like a config
    default_path = os.path.abspath('.')

data = wt.open(easygui.diropenbox(default=os.path.normpath(default_path))+'\primary.wt5')

x = data.axes[0].full
y = data.channels[8].full
mean = np.average(x, weights=y)
sigma = np.sqrt(sum(y * (x - mean)**2) / sum(y))
test_amp = 3*max(y)
print(f'Estimated mean={mean}, sigma={sigma}, and test amplitude={test_amp} for fitting.')

fig = plt.figure()
gs = fig.add_gridspec(1,3, wspace=0)#fig.add_gridspec(2,3, wspace=0)
orig = fig.add_subplot(gs[0,0])
orig.set_title('Raw Data', wrap=True, fontsize='medium')
gfit_plot = fig.add_subplot(gs[0,1])
gfit_plot.set(yticks=[])
gfit_plot.set_title('Gaussian Fit', wrap=True, fontsize='medium')
sgfit_plot = fig.add_subplot(gs[0,2])
sgfit_plot.set(yticks=[])
sgfit_plot.set_title('Skewed Gaussian Fit', wrap=True, fontsize='medium')
'''out_text = fig.add_subplot(gs[1,:])
out_text.set(frame_on=False, xticks=[], yticks=[])'''
#data1.print_tree()
#wt.artists.quick1D(data, 'mono')

orig.plot(x, y)
orig.set_xlabel(f'{data.axes[0].natural_name} ({data.axes[0].units})')
orig.set_ylabel(f'{data.channels[8].natural_name}')
gfit_plot.plot(x, y)
orig.set_xlabel(f'{data.axes[0].natural_name} ({data.axes[0].units})')
orig.set_ylabel(f'{data.channels[8].natural_name}')
sgfit_plot.plot(x, y)
orig.set_xlabel(f'{data.axes[0].natural_name} ({data.axes[0].units})')
orig.set_ylabel(f'{data.channels[8].natural_name}')
plt.show(block=False)

# fitting, peak finding and FWHM
while True:
    try:
        #alpha = float(input("Enter guess for left(-) or right(+) skew in range [-1,1]. (dtype=float):  "))
        alpha = float(easygui.enterbox(msg="Enter guess for left(-) or right(+) skew in range [-1,1]. (dtype=float):  "))
        print(alpha)
        assert -1<=alpha<=1
        break
    except:
        pass
floor = 0.0001
gauss_fit = curve_fit(gauss_curve, x, y, p0=(mean, test_amp, sigma))[0]
skewed_fit = curve_fit(skew, x, y, p0=(mean, test_amp, sigma, alpha, floor))[0]


gauss_peak, _ = find_peaks([gauss_curve(i,*gauss_fit) for i in x])
gauss_width, _, gleft, gright = peak_widths([gauss_curve(i,*gauss_fit) for i in x], gauss_peak, rel_height=0.5)
sgauss_peak, _ = find_peaks([skew(i,*skewed_fit) for i in x])
sgauss_width, _, sgleft, sgright = peak_widths([skew(i,*skewed_fit) for i in x], sgauss_peak, rel_height=0.5)


gfit_plot.plot(x, [gauss_curve(i, *gauss_fit) for i in x], 'k-')

sgfit_plot.plot(x, [skew(i,*skewed_fit) for i in x], 'k-')
out_txt = f"The gaussian fit parameters are: \n\nmean={gauss_fit[0]}, \namplitude={gauss_fit[1]}, \nsigma={gauss_fit[2]}.\n\
FWHM={gauss_width[0]/len(x)*(max(x)-min(x))}\n\n\n\
The skewed gaussian fit parameters are: \n\nmean={skewed_fit[0]}, \namplitude={skewed_fit[1]}, \nsigma={skewed_fit[2]},\
\nalpha={skewed_fit[3]}, \nfloor={skewed_fit[4]}.\n\
FWHM={sgauss_width[0]/len(x)*(max(x)-min(x))}"
'''out_text.text(0.01, 0.99, 
              f'The gaussian fit parameters are: \n\nmean={gauss_fit[0]}, \namplitude={gauss_fit[1]}, \nsigma={gauss_fit[2]}.\n' \
                + f'FWHM={gauss_width[0]/len(x)*(max(x)-min(x))}\n\n\n' \
                    + f'The skewed gaussian fit parameters are: \n\nmean={skewed_fit[0]}, \namplitude={skewed_fit[1]}, \nsigma={skewed_fit[2]}, \nalpha={skewed_fit[3]}, \nfloor={skewed_fit[4]}.\n' \
                        + f'FWHM={sgauss_width[0]/len(x)*(max(x)-min(x))}',
                        va='top', ha='left', wrap=True, fontsize='small', clip_box=TransformedBbox(Bbox([[0.1, 0.1], [0.99, 0.99]]), out_text.transAxes))'''
plt.show(block=False)

easygui.msgbox(msg=out_txt)

