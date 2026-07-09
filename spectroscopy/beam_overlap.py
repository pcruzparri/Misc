'''
Author: Peter L. Cruz Parrilla
'''
__all__ = ['interaction_volume']

import numpy as np
from scipy.ndimage import rotate
import warnings
from functools import reduce
from collections.abc import Sequence

@np.vectorize
def makeGaussian(size, fwhm=3, center=None, amp=1, norm=False):
    def gaussian():
        """Note: Copied without modification from StackOverflow user giessel.
        
        Make a square gaussian kernel.

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

        return np.exp(-4*np.log(2) * ((x-x0)**2 + (y-y0)**2) / fwhm**2)
    g = gaussian()
    if norm: return amp*g/np.sum(g)
    else: return g

def clipped_arr(arr, thresh=0.5):
    return np.array(np.where(arr>thresh)).T

def get_plane_overlap(arr1, arr2, thresh=0.5, clip=(0,20), return_points=False):
    '''Get overlap points of two 3D arrays along a plane.
    '''
    # Combine both arrays and get the counted unique elements.
    all = np.unique(np.concatenate((clipped_arr(arr1, thresh=thresh), (clipped_arr(arr2, thresh=thresh))), axis=0),
                   return_counts=True, axis=0)
    
    # Get items with counts greater than 1 as the overlap points.
    p1_arr = np.array([i for i in all[0][np.where(all[1]>1)] if i[0]==clip[0]])
    p2_arr = np.array([i for i in all[0][np.where(all[1]>1)] if i[0]==clip[1]])
    
    if return_points: return p1_arr, p2_arr #point locations in 3D matrix on two planes of the box along an axis
    else: return p1_arr.size #Due to symmetry, only get the size overlapping point list on one of the planes. 

def interaction_volume(angles, fwhm, size, thresh=0.5, vis=False, amp=[], norm=[], clip=(0, -1)):
    '''Calculate and visualize multiple beam interaction
    
    Parameters
    ----------

    angles: (N, 3) array 
        Contains the set of three rotation angles for N beams in degrees. Rotations are described in the [xy, yz, zx] planes.

    fwhm: (N, ) array
        Contains the FWHM for each beam in microns.
    
    size: (x or y, z) array
        Sets the size of the box with dimensions x*y*z (in microns) where the xy plane slices the beam and z is the default direction of propagation.
    
    thresh: float
        Sets the minimum value used to visualize the gaussian-like beams and the interaction volume. Only used for visualization.

    vis: bool
        Set whether or not to create pyvista point cloud plots. 
    
    amp: (N, ) float array
        Set the intensity of each beam such that the sum of a slice through the axis of propagation is equal to the beam's intensity.

    norm: (N, ) bool array
        Set whether the beam intensity should be normalized and weighed by amp. amp is not applied without norm=True for the beam.

    Return
    ----------

    (vol_um3, vol_cm3)
        Returns a tuple of with the interaction volume of the beams in cubic microns and in cubic centimeters.

    Example
    ----------
    angles = np.array([[-10, 0, 0], [0, 0, 0], [30, 0, 0]])
    fwhm = np.array([10]*len(angles))
    size = (50, 200)

    interaction_volume(angles, fwhm, size, thresh=0.1, vis=True)
    '''
    
    print('Checking parameters...')
    assert False not in list(map(lambda x: isinstance(x, Sequence) | isinstance(x, np.ndarray), 
                                 [angles, fwhm, size, amp, norm])), 'angles, fwhm, size, amp, and norm must be iterables.'
    assert max(fwhm)<0.75*size[0] and max(fwhm)<0.75*size[1], 'Consider smaller FWHM relative to smallest box dimension for more accurate results.'
    if amp or norm:
        assert len(fwhm) == len(angles) == len(amp) == len(norm), 'The first axis of angles, fwhm, amp, and norm must have the same size.' 
    else: 
        assert len(fwhm) == len(angles), 'The first axis of angles and fwhm must have the same size.'
        amp = np.ones(len(fwhm))
        norm = np.array([False]*len(fwhm))

    # Preparing Matrices
    print('Preparing matrices...')
    ls_beams = []
    for beam in range(len(fwhm)):
        print(beam)
        arr2d = makeGaussian(size[0], fwhm=fwhm[beam], amp=amp[beam], norm=norm[beam])
        arr3d = np.array([arr2d]*size[1])
        for rot_ind, rot in enumerate(angles[beam]):
            arr3d = rotate(arr3d, rot, axes=(rot_ind, (rot_ind+1)%3), reshape=False, mode='nearest')
        ls_beams.append(arr3d[clip[0]:clip[1]])
    
    beam_product = reduce(np.multiply, ls_beams)

    print('Checking inclusion of interaction volume in box')
    for beam in range(len(angles)):
        if get_plane_overlap(*(ls_beams[beam], ls_beams[(beam+1)%3]), thresh=thresh, clip=(0,size[1]-1)):
            warnings.warn(f'Beams {beam} and {(beam+1)%3} are overlapping at the edges. Consider extending the box size or changing the angles.')

    assert thresh<np.min([np.max(i) for i in ls_beams]), f'Thresh should be less than {np.min(np.max(ls_beams, axis=1))}'
    if vis:
        from pyvista import Plotter, PolyData, global_theme
        global_theme.color_cycler = 'default'
        print('Visualizing...')
        # Plotting beams
        plotter = Plotter(shape=(1, 2))
        plotter.subplot(0, 0)
        for beam in ls_beams:
            arr_pts = clipped_arr(beam, thresh=thresh)
            plotter.add_mesh(arr_pts,
                             scalars=beam[tuple(arr_pts.T)],
                             cmap='YlGnBu',
                             style='points_gaussian',
                             opacity='linear',
                             show_edges=False)
        plotter.show_grid()
        #Plotting interaction volume
        plotter.subplot(0, 1)
        int_pts = clipped_arr(beam_product, thresh=thresh**len(fwhm))
        plotter.add_mesh(int_pts,
                         scalars=beam_product[tuple(int_pts.T)],
                         cmap='YlGnBu',
                         point_size=10,
                         style='points_gaussian',
                         opacity='linear',
                         show_edges=False)
        plotter.show_grid()
        plotter.link_views()
        plotter.show()

    # Calculating volume
    print('Calculating volume...')
    vol_um3 = np.sum(beam_product) # um^3
    vol_cm3 = vol_um3*(10**-12) # cm^3
    return vol_um3, vol_cm3