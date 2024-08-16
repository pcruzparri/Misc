# imports
import numpy as np
import matplotlib.pyplot as plt
# import multiprocessing
import time


# useful constants
k = 1.380649e-23 # J/K
hbar = 1.054571817e-34 #Js
mu = 9.2740100783e-24 #J/T


# lambda functions
dissociation_time = lambda dyeR, molR, dyeD, molD: (dyeR + molR) ** 2 / (dyeD + molD)


# normal functions
def w_mixing(g_dye: float,
             g_mol: float,
             B0: float,
             m_dye: np.array,
             m_mol: np.array,
             A_dye: np.array,
             A_mol: np.array):
    """
    :param g_dye: float,
    g-factor for the dye
    :param g_mol: float,
    g-factor for the molecule
    :param B0: float,
    applied, static magnetic field
    :param m_dye: array
    nuclear spin state for NMR-active nuclei in the dye. Maps to A_dye.
    :param m_mol: array
    nuclear spin state for NMR-active nuclei in the dye. Maps to A_mol. Does not include the spin states of C_alpha.
    :param A_dye: array
    nuclear spin state for NMR-active nuclei in the dye. Maps to m_dye.
    :param A_mol: array
    nuclear spin state for NMR-active nuclei in the dye. Maps to m_mol. Does not include the hyperfine coupling constant of C_alpha.

    :return:
    """

    assert m_dye.shape[0] == A_dye.shape[0]
    assert m_mol.shape[0] == A_mol.shape[0]

    conf_dye = np.sum(np.multiply(m_dye, A_dye))
    conf_mol = np.sum(np.multiply(m_mol, A_dye))

    a = 0.5 * ((g_dye + g_mol) * mu * B0 / hbar + conf_dye - 0.5 * A0 - conf_mol)
    b = 0.5 * ((g_dye + g_mol) * mu * B0 / hbar + conf_dye + 0.5 * A0 - conf_mol)

    return a, b


