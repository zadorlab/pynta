import os
import shutil

from ase.build import fcc111, fcc211, fcc100
from ase.io import write
from ase.calculators.espresso import Espresso
from ase.calculators.socketio import SocketIOCalculator

from sella import Sella

# from ase.optimize import LBFGS
# from gpaw import GPAW, PW
# os.environ['GPAW_SETUP_PATH']
# avaiable options for slab build are: fcc100, fcc110, bcc100, bcc110, bcc111, fcc111, hcp0001, hcp10m10, diamond100, diamond111


class GetSlab:
    def __init__(self, surface_type, symbol, a, repeats, vacuum, slab_name,
                 pseudopotentials, pseudo_dir):
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
        # slab = fcc111(self.symbol, self.repeats, self.a,
        #               self.vacuum, orthogonal=False, periodic=True)
        slab = fcc111(self.symbol, self.repeats, self.a,
                      self.vacuum)
        GetSlab.prepare_slab_opt(self, slab)

    def opt_fcc211(self):
        ''' Optimize fcc211 slab '''
        slab = fcc211(self.symbol, self.repeats, self.a,
                      self.vacuum)
        GetSlab.prepare_slab_opt(self, slab, )

    def opt_fcc100(self):
        ''' Optimize fcc100 slab '''
        slab = fcc100(self.symbol, self.repeats, self.a,
                      self.vacuum)
        GetSlab.prepare_slab_opt(self, slab)

    def prepare_slab_opt(self, slab):
        ''' Prepare slab optimization with Quantum Espresso '''
        # setting up calculator
        # calc = GPAW(xc='PBE', mode = 'pw', kpts=(4, 4, 4))
        # slab.set_calculator(calc)
        # dyn = LBFGS(slab, trajectory = 'slab_Cu.traj')
        # dyn.run(fmax=0.01)

        unixsocket = self.slab_name
        socketpath = f'/tmp/ipi_{unixsocket}'
        if os.path.exists(socketpath):
            os.remove(socketpath)
        if os.path.exists(unixsocket):
            shutil.rmtree(unixsocket)
        os.makedirs(unixsocket)

        label = os.path.join(unixsocket, self.slab_name)

        espresso = Espresso(command='/home/ehermes/local/bin/mpirun -np 48 /home/ehermes/local/bin/pw.x -inp PREFIX.pwi --ipi {{unixsocket}}:UNIX > PREFIX.pwo'
                            .format(unixsocket=unixsocket),
                            label=label,
                            pseudopotentials=self.pseudopotentials,
                            pseudo_dir=self.pseudo_dir,
                            kpts=(3, 3, 1),
                            occupations='smearing',
                            smearing='marzari-vanderbilt',
                            degauss=0.01,  # Rydberg
                            ecutwfc=40,  # Rydberg
                            nosym=True,  # Allow symmetry breaking during optimization
                            conv_thr=1e-11,
                            mixing_mode='local-TF',
                            )
        with SocketIOCalculator(espresso, unixsocket=unixsocket) as calc:
            slab.calc = calc
            opt = Sella(slab, order=0, delta0=1e-2, trajectory=label + '.traj')
            opt.run(fmax=0.01)
        ener = slab.get_potential_energy()
        force = slab.get_forces()
        write(self.slab_name + '.xyz', slab)

        # To run quickly on my Mac
        # espresso = Espresso(command='mpirun -np 8 /Users/mgierad/00_SANDIA_WORK/03_codes/build/q-e-qe-6.4.1/bin/pw.x -inp PREFIX.pwi --ipi {{unixsocket}}:UNIX > PREFIX.pwo'
        #                 .format(unixsocket=unixsocket),
        #                 label=label,
        #                 pseudopotentials=self.pseudopotentials,
        #                 pseudo_dir=self.pseudo_dir,
        #                 kpts=(3, 3, 1),
        #                 occupations='smearing',
        #                 smearing='marzari-vanderbilt',
        #                 degauss=0.01,  # Rydberg
        #                 ecutwfc=40,  # Rydberg
        #                 nosym=True,  # Allow symmetry breaking during optimization
        #                 conv_thr=1e-11,
        #                 mixing_mode='local-TF',
        #                 )

        # with SocketIOCalculator(espresso, unixsocket=unixsocket) as calc:
        #     slab.calc = calc
        #     opt = Sella(slab, order=0, delta0=1e-2, trajectory=label + '.traj')
        #     opt.run(fmax=0.01)
        # ener = slab.get_potential_energy()
        # force = slab.get_forces()
        # write(self.slab_name + '.xyz', slab)
