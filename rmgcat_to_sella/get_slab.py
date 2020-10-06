#!/usr/bin/env python3
from rmgcat_to_sella.balsamcalc import EspressoBalsamSocketIO

from ase.build import fcc111, fcc211, fcc100
from ase.build import bcc111, bcc110, hcp0001
from ase.build import diamond111, diamond100
from ase.optimize import BFGSLineSearch
from ase.io import write


class GetSlab:
    def __init__(
            self,
            surface_type,
            symbol,
            a,
            repeats,
            vacuum,
            slab_name,
            pseudopotentials,
            pseudo_dir,
            balsam_exe_settings,
            calc_keywords,
            creation_dir
    ):
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
        repeats : tuple
            specify reapeats in (x, y, z) direction,
            eg. (3, 3, 1)
        vacuum : float
            amount of vacuum in the z direction (Units: Angstrem)
        slab_name : str
            a user defined name of the slab
        slab_extension : str
            usually, it should be 'xyz'
        pseudopotentials : dict(str='str')
            a dictionary with all pseudopotentials, as for Quantum Espresso
            keys() - symbol
            values() - file name of a pseudopotential
            e.g. dict(Cu='Cu.pbe-spn-kjpaw_psl.1.0.0.UPF')
        pseudo_dir : str
            a path to the QE's pseudopotentials main directory
            e.g.
            '/home/mgierad/espresso/pseudo'

        '''
        self.surface_type = surface_type
        self.symbol = symbol
        self.a = a
        self.repeats = repeats
        self.vacuum = vacuum
        self.slab_name = slab_name
        self.pseudopotentials = pseudopotentials
        self.pseudo_dir = pseudo_dir
        self.balsam_exe_settings = balsam_exe_settings
        self.calc_keywords = calc_keywords
        self.creation_dir = creation_dir

    def run_slab_opt(self):
        ''' Run slab optimization '''
        if self.surface_type == 'fcc111':
            GetSlab.opt_fcc111(self)
        elif self.surface_type == 'fcc211':
            GetSlab.opt_fcc211(self)
        elif self.surface_type == 'fcc100':
            GetSlab.opt_fcc100(self)
        else:
            print('{} not implemented. Avaiable parameters are:'.format(
                self.surface_type))
            print('fcc111, fcc100, fcc211')

    def opt_fcc111(self):
        ''' Optimize fcc111 slab '''
        slab = fcc111(self.symbol, self.repeats_surface, self.a, self.vacuum)
        self.prepare_slab_opt(slab)

    def opt_fcc211(self):
        ''' Optimize fcc211 slab '''
        slab = fcc211(self.symbol, self.repeats_surface, self.a, self.vacuum)
        self.prepare_slab_opt(slab)

    def opt_fcc100(self):
        ''' Optimize fcc100 slab '''
        slab = fcc100(self.symbol, self.repeats_surface, self.a, self.vacuum)
        self.prepare_slab_opt(slab)

    def opt_bcc111(self):
        ''' Optimize bcc111 slab '''
        slab = bcc111(self.symbol, self.repeats_surface, self.a, self.vacuum)
        self.prepare_slab_opt(slab)

    def opt_bcc110(self):
        ''' Optimize bcc110 slab '''
        slab = bcc110(self.symbol, self.repeats_surface, self.a, self.vacuum)
        self.prepare_slab_opt(slab)

    def opt_hcp0001(self, c=None):
        ''' Optimize hcp0001 slab '''
        slab = hcp0001(self.symbol, self.repeats_surface, self.a, c,
                       self.vacuum)
        self.prepare_slab_opt(slab)

    def opt_diamond111(self):
        ''' Optimize diamond111 slab '''
        slab = diamond111(self.symbol, self.repeats_surface, self.a,
                          self.vacuum)
        self.prepare_slab_opt(slab)

    def opt_diamond100(self):
        ''' Optimize diamond100 slab '''
        slab = diamond100(self.symbol, self.repeats_surface, self.a,
                          self.vacuum)
        self.prepare_slab_opt(slab)

    def prepare_slab_opt(self, slab):
        ''' Prepare slab optimization with Quantum Espresso '''

        job_kwargs = self.balsam_exe_settings.copy()
        # job_kwargs.update([('user_workdir',cwd)])
        QE_keywords_slab = self.calc_keywords.copy()
        # QE_keywords.update([('kpts',self.repeats)])
        # Not sure of intended behavior, but an example to show
        # you can change keys as necessary here

        slab.calc = EspressoBalsamSocketIO(
            workflow='QE_Socket',
            job_kwargs=job_kwargs,
            pseudopotentials=self.pseudopotentials,
            pseudo_dir=self.pseudo_dir,
            **QE_keywords_slab
        )

        label = self.slab_name

        opt = BFGSLineSearch(atoms=slab, trajectory=label + '.traj')
        opt.run(fmax=0.01)
        slab.get_forces()
        slab.calc.close()
        write(self.slab_name + '.xyz', slab)
