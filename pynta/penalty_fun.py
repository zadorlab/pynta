from typing import List, Tuple
from ase import Atoms
from ase.io.trajectory import Trajectory
import numpy as np
from scipy.optimize import minimize


def skew(x: np.ndarray) -> np.ndarray:
    assert len(x) == 3
    X = np.zeros((3, 3))
    X[0, 1] = -x[2]
    X[0, 2] = x[1]
    X[1, 2] = -x[0]
    X -= X.T
    return X


class AdsorbatePlacer:
    def __init__(self,
                 slab: Atoms,
                 adsorbate: Atoms,
                 bonds: List[Tuple[int, int]],
                 dists: List[float],
                 initial_height: float = 1.0,
                 weight: float = 1.0,
                 scale: float = 1.0,
                 trajectory: str = None) -> None:
        self.slab = slab
        self.adsorbate = adsorbate

        zmax_slab = np.max(self.slab.positions[:, 2])
        zmin_ads = np.min(self.adsorbate.positions[:, 2])
        self.adsorbate.positions[:, 2] += initial_height + zmax_slab - zmin_ads
        self.ads_ref = self.slab + self.adsorbate

        self.bonds = bonds
        self.dists = dists

        self.weight = weight
        self.scale = scale
        self.trajectory = trajectory
        if self.trajectory is not None:
            self.trajectory = Trajectory(self.trajectory, 'w', self.ads_ref)

    @property
    def nslab(self) -> int:
        ''' How many atoms are in slab

        Returns
        -------
        nslab : int
            number of atoms in slab

        '''
        return len(self.slab)

    @property
    def nadsorbate(self) -> int:
        ''' How many atoms are in adsorbate

        Returns
        -------
        nadsorbate : int
            number of atoms in adsorbate

        '''
        return len(self.adsorbate)

    @property
    def nx(self) -> int:
        ''' Return the product of array elements over a given axis

        Returns
        -------
        nx : int
            product of array of elements self.adsorbate.positions.shape

        '''
        return np.prod(self.adsorbate.positions.shape)

    def set_y(self, yin: np.ndarray) -> None:
        ''' Set y axis

        Parameters
        ----------
        yin : np.ndarray
            an old y which is used to set a new one

        '''
        ads = self.adsorbate.copy()
        center = yin[:3]
        axis = yin[3:]
        angle = np.linalg.norm(axis) * 180 / np.pi
        ads.positions += center - ads.positions.mean(0)
        ads.rotate(angle, v=axis, center=center)
        self.ads_ref.positions[self.nslab:] = ads.positions

    def penalty(self, y: np.ndarray) -> float:
        ''' Define a penalty function

        Parameters
        ----------
        y : np.ndarray
            an axis

        Returns
        -------
        pen : float
            a weighted penalty function

        '''
        self.set_y(y)
        pen = 0.
        for bond, dist in zip(self.bonds, self.dists):
            pen += (self.ads_ref.get_distance(*bond) - dist)**2
        return pen * self.weight

    def dpdx(self, y: np.ndarray) -> np.ndarray:
        ''' Return a flaten dpdx array

        Parameters
        ----------
        y : np.ndarray
            an axis

        Returns
        -------
        dpdx : np.ndarray
            a flaten dpdx array

        '''
        self.set_y(y)
        dpdx = np.zeros((self.nadsorbate, 3))
        for bond, dist in zip(self.bonds, self.dists):
            xij = self.ads_ref.get_distance(*bond, vector=True)
            dij = np.linalg.norm(xij)
            deriv = 2 * (dij - dist) * xij / dij
            if bond[0] >= self.nslab:
                dpdx[bond[0] - self.nslab] -= deriv
            if bond[1] >= self.nslab:
                dpdx[bond[1] - self.nslab] += deriv
        return dpdx.ravel()

    def penalty_jac(self, y: np.ndarray) -> np.ndarray:
        ''' Deal with jacobians

        Parameters
        ----------
        y : np.ndarray
            an axis

        Returns
        -------
        self.dpdx(y) @ self.dxdy(y) : np.ndarray
            a product of multiplied matrix

        '''
        return self.dpdx(y) @ self.dxdy(y)

    def x(self, y: np.ndarray) -> np.ndarray:
        ''' Return a flaten array of adsorbate positions

        Parameters
        ----------
        y : np.ndarray
            an axis

        Returns
        -------
        np.ndarray
            a flaten array of adsorbate positions

        '''
        self.set_y(y)
        return self.ads_ref.positions[self.nslab:].ravel().copy()

    def dxdy(self, y: np.ndarray) -> np.ndarray:
        ''' Get a dxdy array

        Parameters
        ----------
        y : np.ndarray
            an axis

        Returns
        -------
        dxdy : np.ndarray
            a reshaped (-1, 6) dxdy array

        '''
        self.set_y(y)
        dxdy = np.zeros((self.nadsorbate, 3, 6))
        dxdy[:, :, :3] = np.eye(3)[np.newaxis]

        v = y[3:]
        theta = np.linalg.norm(v)
        u = v / theta
        sint = np.sin(theta)
        cost = np.cos(theta)
        cosct = (1 - cost) / theta
        dX = self.adsorbate.positions.copy()
        dX -= dX.mean(0)
        for i, xi in enumerate(dX):
            dxdy[i, :, 3:] = (
                -sint * np.outer(xi, u)
                + cost * np.outer(np.cross(u, xi), u)
                - sint * (skew(xi) + np.outer(np.cross(u, xi), u)) / theta
                + sint * (u @ xi) * np.outer(u, u)
                + cosct * ((u @ xi) * np.eye(3)
                           + np.outer(u, xi)
                           - 2 * (u @ xi) * np.outer(u, u))
            )
        return dxdy.reshape((-1, 6))

    def energy(self, y: np.ndarray) -> float:
        ''' Return a sum of penalty function and total energy

        Parameters
        ----------
        y : np.ndarray
            an axis

        Returns
        -------
        energy : float
            a sum of penalty function and total energy

        '''
        return self.penalty(y) + self.ads_ref.get_potential_energy()

    def gradient(self, y: np.ndarray) -> np.ndarray:
        ''' Return a scaled gradient

        Parameters
        ----------
        y : np.ndarray
            an axis

        Returns
        -------
        np.ndarray
            a scaled gradient

        '''
        jac = self.penalty_jac(y)
        forces = self.ads_ref.get_forces()
        if self.trajectory is not None:
            self.trajectory.write()
        jac -= forces[self.nslab:].ravel() @ self.dxdy(y)
        return jac * self.scale

    def optimize(self) -> Atoms:
        ''' Minimize penalty function defined as:

        .. math::
            P = min\\left[E_{xTB} + \\sum_{i = 1}^{N}{\\left(r^{ts\\_guess}_{X_i M} - \\overline{{r^{min}_{X_{i} M}}}\\right)^2}\\right]

        Where:

        :math:`r^{ts\\_guess}_{X_i M}` -- distance between reacting atom and the nearest surface atom ina TS guess;

        :math:`\overline{{r^{min}_{X_{i} M}}}` -- target bond distance, i.e. an average distance betweenreacting atom and the nearest surface atom calculated for all symmetry distinct minima;

        :math:`N` -- number of atoms included in a penalty function definition;

        :math:`E_{xTB}` -- total energy of the system calculated using a robust semi-empirical xTB code

        Returns
        -------
        ads_ref : Atoms
            a new structure with improve location of TS guess for which the
            penalty function is in minimum

        '''
        y0 = np.zeros(6)
        y0[:3] = self.ads_ref.positions[self.nslab:].mean(0)
        y0[3:] = 0.001

        res1 = minimize(self.penalty, y0, jac=self.penalty_jac,
                        method='l-bfgs-b')
        y1 = res1['x']

        res2 = minimize(self.energy, y1, jac=self.gradient, method='bfgs',
                        options={'disp': True})

        # an alternative but usualy bfgs is better/more stable here
        # res2 = minimize(self.energy, y1, jac=self.gradient,
        #                 method='l-bfgs-b',
        #                 options={'iprint': 1})

        print(res2)

        y2 = res2['x']
        self.set_y(y2)
        return self.ads_ref.copy()
