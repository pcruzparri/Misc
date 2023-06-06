# imports
import numpy as np
from scipy.optimize import curve_fit, root
import WrightTools as wt
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.widgets import TextBox
from matplotlib.transforms import TransformedBbox, Bbox
from scipy import signal as s
from scipy.stats import norm
from scipy.special import erfc, erf
from scipy.signal import find_peaks, peak_widths
import math
import sys


# fitting functions
def gauss_curve(x, x0, a, sigma):
    return a * np.exp(-(x-x0)**2/(2*sigma**2))

def skew2(x, sigmag, mu, alpha, c, a):
    normpdf = (1 / (sigmag * np.sqrt(2 * math.pi))) * np.exp(-(np.power((x - mu), 2) / (2 * np.power(sigmag, 2))))
    normcdf = (0.5 * (1 + erf((alpha * ((x - mu) / sigmag)) / (np.sqrt(2)))))
    return 2 * a * normpdf * normcdf + c, np.max(normpdf)

def skew(x, mu, a, sigmag, alpha, c):
    return skew2(x, sigmag, mu, alpha, c, a)[0]


# file inputs and setup
data = wt.open(str(sys.argv[1]))

fig = plt.figure()
gs = fig.add_gridspec(2,3, wspace=0)
orig = fig.add_subplot(gs[0,0])
orig.set_title('Raw Data', wrap=True, fontsize='medium')
gfit_plot = fig.add_subplot(gs[0,1])
gfit_plot.set(yticks=[])
gfit_plot.set_title('Gaussian Fit', wrap=True, fontsize='medium')
sgfit_plot = fig.add_subplot(gs[0,2])
sgfit_plot.set(yticks=[])
sgfit_plot.set_title('Skewed Gaussian Fit', wrap=True, fontsize='medium')
out_text = fig.add_subplot(gs[1,:])
out_text.set(frame_on=False, xticks=[], yticks=[])
#data1.print_tree()
#wt.artists.quick1D(data, 'mono')

x = data.axes[0].full
y = data.channels[8].full
mean = np.average(x, weights=y)
sigma = np.sqrt(sum(y * (x - mean)**2) / sum(y))
test_amp = 3*max(y)
print(f'Estimated mean={mean}, sigma={sigma}, and test amplitude={test_amp} for fitting.')
orig.plot(x, y)
gfit_plot.plot(x, y)
sgfit_plot.plot(x, y)
plt.show(block=False)

# fitting, peak finding and FWHM
while True:
    try:
        alpha = float(input("Enter guess for left(-) or right(+) skew in range [-1,1]. (dtype=float):  "))
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

out_text.text(0.01, 0.99, 
              f'The gaussian fit parameters are: \n\nmean={gauss_fit[0]}, \namplitude={gauss_fit[1]}, \nsigma={gauss_fit[2]}.\n' \
                + f'FWHM={gauss_width[0]/len(x)*(max(x)-min(x))}\n\n\n' \
                    + f'The skewed gaussian fit parameters are: \n\nmean={skewed_fit[0]}, \namplitude={skewed_fit[1]}, \nsigma={skewed_fit[2]}, \nalpha={skewed_fit[3]}, \nfloor={skewed_fit[4]}.\n' \
                        + f'FWHM={sgauss_width[0]/len(x)*(max(x)-min(x))}',
                        va='top', ha='left', wrap=True, fontsize='small', clip_box=TransformedBbox(Bbox([[0.1, 0.1], [0.99, 0.99]]), out_text.transAxes))
plt.show()

