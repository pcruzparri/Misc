'''
Author: Peter L. Cruz Parrilla
'''
__all__ = ['interaction_volume']

import numpy as np
from scipy.ndimage import rotate
import pyvista as pv
import warnings

@np.vectorize
def makeGaussian(size, fwhm=3, center=None):
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

def interaction_volume(angles, fwhm, size, thresh=0.5, vis=False):
    '''Calculate and visualize multiple beam interaction
    
    Parameters:
    angles: (N, 3) array 
        Contains the set of three rotation angles for N beams in degrees.

    fwhm: (N, ) array
        Contains the FWHM for each beam in microns.
    
    size: (xy, z) array
        Sets the size of the box with dimensions x*y*z (in microns) where the xy plane slices the beam and z is the default direction of propagation.
    
    thresh: float
        Sets the minimum value used to visualize the gaussian-like beams and the interaction volume. Only used for visualization.

    vis: bool
        Set whether or not to create pyvista point cloud plots.

    Return:
    (vol_um3, vol_cm3)
        Returns a tuple of with the interaction volume of the beams in cubic microns and in cubic centimeters.
    '''

    assert max(fwhm)<0.75*size[0] and max(fwhm)<0.75*size[1], 'Consider smaller FWHM relative to smallest box dimension for more accurate results.'

    # Preparing Matrices
    print('Preparing matrices...')
    ls_beams = []
    assert len(fwhm) == len(angles), 'The first axis of angles and fwhm must have the same size.'
    for beam in range(len(fwhm)):
        print(beam)
        arr2d = makeGaussian(size[0], fwhm=fwhm[beam])
        arr3d = np.array([arr2d]*size[1])
        for rot_ind, rot in enumerate(angles[beam]):
            arr3d = rotate(arr3d, rot, axes=(rot_ind, (rot_ind+1)%3), reshape=False, mode='nearest')
        ls_beams.append(arr3d)

    print('Checking inclusion of interaction volume in box')
    for beam in range(3):
        if get_plane_overlap(*(ls_beams[beam],ls_beams[(beam+1)%3]), thresh=thresh, clip=(0,size[1]-1)):
            warnings.warn(f'Beams {beam} and {(beam+1)%3} are overlapping at the edges. Consider extending the box size or changing the angles.')

    vis = True
    if vis:
        print('Visualizing...')
        # Plotting beams
        plotter = pv.Plotter()
        plotter.add_mesh(clipped_arr(ls_beams[0], thresh=thresh), color='g', point_size=10, opacity=0.7)
        plotter.add_mesh(clipped_arr(ls_beams[1], thresh=thresh), color='r', point_size=10, opacity=0.7)
        plotter.add_mesh(clipped_arr(ls_beams[2], thresh=thresh), color='b', point_size=10, opacity=0.7)
        plotter.show_grid()
        plotter.show()
        #Plotting interaction volume
        int_plotter = pv.PolyData(clipped_arr(np.multiply(np.multiply(ls_beams[0], ls_beams[1]), ls_beams[2]), thresh=thresh))
        int_plotter.plot()

    # Calculating volume
    print('Calculating volume...')
    vol_um3 = np.sum(np.multiply(np.multiply(ls_beams[0], ls_beams[1]), ls_beams[2])) # um^3
    vol_cm3 = vol_um3*(10**-12)
    return vol_um3, vol_cm3
