import pyvista as pv
import re
import numpy as np
from scipy.spatial.transform import Rotation as R

# Z-matrix reader

fpath = r"C:\Users\Peter\src\zmat.txt"



def from_angle(p1, p2, angle, p2_distance):
    r21 = (p1-p2)/np.linalg.norm(p1-p2)*p2_distance
    axis = np.cross(p1, p2)
    r = R.from_rotvec(axis/np.linalg.norm(axis)*angle/360*2*np.pi)
    return r.apply(r21)+p2


class Molecule:
    def __init__(self, struct_fpath, struct_ftype):
        if struct_ftype = 'zmat':
            self.structure = self._from_zmat(struct_fpath)

    def _from_zmat(self):
        with open(fpath) as f:
            data = []
            for row in f.readlines():
                data += [re.split('\s+', row.strip('\n| '))]

            zmat = {}
            zmat_const = dict((tuple(const) for const in data[data.index(['']) + 1:]))
            labels = ['atom', 'bond_atom', 'bond_length', 'angle_atom', 'angle', 'dihedral_atom', 'dihedral']

            for ind, atom in enumerate(data[:data.index([''])]):
                zmat[ind] = dict(zip(labels, atom))

            structure =








"""
plotter = pv.Plotter(shape=(1, 2))


p1 = np.array([1, 2, 3])
p2 = np.array([5, 3, 7])
p3 = from_angle(p1, p2, 180, 10)
points = [p1, p2, p3]
colors = ['r', 'g', 'b']



for ind, p in enumerate(points):
    sphere = pv.Sphere(center=p)
    plotter.add_mesh(sphere, color=colors[ind])
plotter.add_axes()

plotter.subplot(0, 1)
for ind, p in enumerate(points[:2]):
    sphere = pv.Sphere(center=p)
    plotter.add_mesh(sphere, color=colors[ind])
plotter.add_axes()

plotter.show()
"""



