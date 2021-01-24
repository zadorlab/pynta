from pathlib import PosixPath
from typing import Dict, Tuple

from ase.build import fcc111, fcc211, fcc100
from ase.build import bcc111, bcc110, hcp0001
from ase.build import diamond111, diamond100
from ase.optimize import BFGSLineSearch
from ase.io import write
from ase import Atoms

import os


class GetSlab:
    def __init__(
            self,
            socket_calculator: str,
            surface_type: str,
            symbol: str,
            a: str,
            repeats_surface: Tuple[int, int, int],
            vacuum: float,
            slab_name: str,
            pseudopotentials: Dict[str, str],
            pseudo_dir: str,
            balsam_exe_settings: Dict[str, float],
            calc_keywords: Dict[str, str],
            creation_dir: PosixPath) -> None:
        ''' A class for preparing and optimizing a user defined slab

        Parameters
        __________

        surface_type : str
            type of the surface. Available options are:
            fcc111, fcc211, fcc100
        symbol : str
            atomic symbol of the studied metal surface
            e.g. 'Cu'
        a : float
            a lattice constant
        repeats_surface : tuple
            specify reapeats in (x, y, z) direction,
            eg. (3, 3, 1)
        vacuum : float
            amount of vacuum in the z direction (Units: Angstrem)
        slab_name : str
            a user defined name of the slab
        pseudopotentials : dict(str='str')
            a dictionary with all pseudopotentials, as for Quantum Espresso
            keys() - symbol
            values() - file name of a pseudopotential
            e.g.

            >>> dict(Cu='Cu.pbe-spn-kjpaw_psl.1.0.0.UPF')
        pseudo_dir : str
            a path to the QE's pseudopotentials main directory
            e.g.
            ``'/home/mgierad/espresso/pseudo'``
        balsam_exe_settings : dict{str:int}
            a dictionary with balsam execute parameters (cores, nodes, etc.),
            e.g.

            >>> balsam_exe_settings = {'num_nodes': 1,
                                   'ranks_per_node': 48,
                                   'threads_per_rank': 1}

        calc_keywords : dict{str:str}
            a dictionary with parameters to run DFT package. Quantum Espresso
            is used as default, e.g.

            calc_keywords = {'kpts': (3, 3, 1), 'occupations': 'smearing',
                            'smearing':  'marzari-vanderbilt',
                            'degauss': 0.01, 'ecutwfc': 40, 'nosym': True,
                            'conv_thr': 1e-11, 'mixing_mode': 'local-TF'}
        creation_dir : posix
            a posix path to the main working directory

        '''
        self.socket_calculator = socket_calculator
        self.surface_type = surface_type
        self.symbol = symbol
        self.a = a
        self.repeats_surface = repeats_surface
        self.vacuum = vacuum
        self.slab_name = slab_name
        self.pseudopotentials = pseudopotentials
        self.pseudo_dir = pseudo_dir
        self.balsam_exe_settings = balsam_exe_settings
        self.calc_keywords = calc_keywords
        self.creation_dir = creation_dir

    def run_slab_opt(self) -> None:
        ''' Run slab optimization '''
        if self.surface_type == 'fcc111':
            GetSlab.opt_fcc111(self)
        elif self.surface_type == 'fcc211':
            GetSlab.opt_fcc211(self)
        elif self.surface_type == 'fcc100':
            GetSlab.opt_fcc100(self)
        elif self.surface_type == 'bcc111':
            GetSlab.opt_bcc111(self)
        elif self.surface_type == 'bcc110':
            GetSlab.opt_bcc110(self)
        elif self.surface_type == 'hcp0001':
            GetSlab.opt_hcp0001(self)
        elif self.surface_type == 'diamond111':
            GetSlab.opt_diamond111(self)
        elif self.surface_type == 'diamond100':
            GetSlab.opt_diamond100(self)
        else:
            print('{} not implemented. Avaiable parameters are:'.format(
                self.surface_type))
            print('fcc111, fcc211, fcc100, bcc111, bcc110, hcp0001, '
                  'diamond111, diamond100')

    def opt_fcc111(self) -> None:
        ''' Optimize fcc111 slab '''
        slab = fcc111(self.symbol, self.repeats_surface, self.a, self.vacuum)
        self.prepare_slab_opt(slab)

    def opt_fcc211(self) -> None:
        ''' Optimize fcc211 slab '''
        slab = fcc211(self.symbol, self.repeats_surface, self.a, self.vacuum)
        self.prepare_slab_opt(slab)

    def opt_fcc100(self) -> None:
        ''' Optimize fcc100 slab '''
        slab = fcc100(self.symbol, self.repeats_surface, self.a, self.vacuum)
        self.prepare_slab_opt(slab)

    def opt_bcc111(self) -> None:
        ''' Optimize bcc111 slab '''
        slab = bcc111(self.symbol, self.repeats_surface, self.a, self.vacuum)
        self.prepare_slab_opt(slab)

    def opt_bcc110(self) -> None:
        ''' Optimize bcc110 slab '''
        slab = bcc110(self.symbol, self.repeats_surface, self.a, self.vacuum)
        self.prepare_slab_opt(slab)

    def opt_hcp0001(
            self,
            c: int = None) -> None:
        ''' Optimize hcp0001 slab '''
        slab = hcp0001(self.symbol, self.repeats_surface, self.a, c,
                       self.vacuum)
        self.prepare_slab_opt(slab)

    def opt_diamond111(self) -> None:
        ''' Optimize diamond111 slab '''
        slab = diamond111(self.symbol, self.repeats_surface, self.a,
                          self.vacuum)
        self.prepare_slab_opt(slab)

    def opt_diamond100(self) -> None:
        ''' Optimize diamond100 slab '''
        slab = diamond100(self.symbol, self.repeats_surface, self.a,
                          self.vacuum)
        self.prepare_slab_opt(slab)

    def prepare_slab_opt(
            self,
            slab: Atoms) -> None:
        ''' Prepare slab optimization with Quantum Espresso '''

        balsamcalc_module = __import__('pynta.balsamcalc', fromlist=[
            self.socket_calculator])

        sock_calc = getattr(balsamcalc_module, self.socket_calculator)

        # n_kpts = IO().get_kpoints(self.repeats_surface)

        job_kwargs = self.balsam_exe_settings.copy()
        extra_calc_keywords = self.calc_keywords.copy()
        # add kpoints and distribute it among nodes = n_kpts
        # extra_calc_keywords['kpts'] = self.repeats_surface
        # extra_calc_keywords['job_args'] = '-nk {}'.format(n_kpts)
        # change how k-points are distrubuted among nodes
        # job_kwargs.update([('num_nodes', n_kpts)])

        slab.pbc = (True, True, False)

        if self.socket_calculator == 'EspressoBalsamSocketIO':
            slab.calc = sock_calc(
                workflow='QE_Socket',
                job_kwargs=job_kwargs,
                pseudopotentials=self.pseudopotentials,
                pseudo_dir=self.pseudo_dir,
                **extra_calc_keywords
            )
        else:
            slab.calc = sock_calc(
                workflow='QE_Socket',
                job_kwargs=job_kwargs,
                **extra_calc_keywords
            )

        fname = os.path.join(self.creation_dir, self.slab_name)

        opt = BFGSLineSearch(atoms=slab, trajectory=fname + '.traj')
        opt.run(fmax=0.01)
        slab.get_potential_energy()
        slab.get_forces()
        slab.calc.close()
        write(fname + '.xyz', slab)
