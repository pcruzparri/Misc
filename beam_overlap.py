import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import rotate
from mpl_toolkits import mplot3d
import pyvista as pv

@np.vectorize
def makeGaussian(size, fwhm=3, center=None):
    """ Make a square gaussian kernel.

    size is the length of a side of the square
    fwhm is full-width-half-maximum, which
    can be thought of as an effective radius.
    """

    x = np.arange(0, size, 1, float)
    y = x[:, np.newaxis]

    if center is None:
        x0 = y0 = size // 2
    else:
        x0 = center[0]
        y0 = center[1]

    '''kernel = np.zeros((size, size))
    y, x = np.ogrid[-y0:y0 + 1, -x0:x0 + 1]
    mask = x ** 2 + y ** 2 <= fwhm ** 2
    kernel[mask] = 1
    return kernel'''

    return np.exp(-4*np.log(2) * ((x-x0)**2 + (y-y0)**2) / fwhm**2)


def clipped_arr(arr, thresh=0.5):
    return np.array(np.where(arr>thresh)).T


# Constants
print('Setting constants...')
thresh = 0.1
n = 61 # um
d = 201
abc_fwhm = np.array([10]*3) # um
abc_angles = np.array([(-1.5, 0, 0), (30, 0, 0), (0, 0, 0)]) #deg
thickness = 10

# Preparing Matrices
print('Preparing matrices...')
ls_beams = []
assert(len(abc_fwhm) == len(abc_angles))
for beam in range(len(abc_fwhm)):
    print(beam)
    arr2d = makeGaussian(n, fwhm=abc_fwhm[beam])
    arr3d = np.array([arr2d]*d)
    for rot_ind, rot in enumerate(abc_angles[beam]):
        arr3d = rotate(arr3d, rot, axes=(rot_ind, (rot_ind+1)%3), reshape=False, mode='nearest')
    ls_beams.append(arr3d)

# Calculating volume
print('Calculating volume...')
vol = np.sum(np.multiply(np.multiply(ls_beams[0], ls_beams[1]), ls_beams[2])) # um
vol2 = vol*(10**-12)
print(f'The calculated volume is {vol} $um^3$, or {vol2} $cm^3$')

vis = True
if vis:
    print('Visualizing...')
    # Creating plot
    plotter = pv.Plotter()
    plotter.add_mesh(clipped_arr(ls_beams[0], thresh=thresh), color='g')
    plotter.add_mesh(clipped_arr(ls_beams[1], thresh=thresh), color='r')
    plotter.add_mesh(clipped_arr(ls_beams[2], thresh=thresh), color='b')
    plotter.show()
    int_plotter = pv.PolyData(clipped_arr(np.multiply(np.multiply(ls_beams[0], ls_beams[1]), ls_beams[2]), thresh=0.1))
    #int_plotter.plot()
    int_volume = int_plotter.delaunay_3d(alpha=10)
    shell = int_volume.extract_geometry()
    print(shell.volume, 'um^3, or', shell.volume/(10**12), 'cm^3')
    shell.plot()
