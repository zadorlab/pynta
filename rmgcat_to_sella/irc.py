import os

from ase.io import read, write
# from pathlib import Path
from rmgcat_to_sella.ts import TS


class IRC():
    def __init__(
            self,
            facetpath,
            slab,
            repeats,
            ts_dir,
            yamlfile,
            pseudopotentials,
            pseudo_dir,
            executable,
            balsam_exe_settings,
            calc_keywords
            ):
        ''' Initializing

        Parameters:
        ___________
        facetpath : str
            a path to the workflow's main dir
            e.g. 'Cu_111'
        slab : str
            a '.xyz' file name with the optimized slab
            e.g.
            'Cu_111_slab_opt.xyz'
        repeats: tuple
            specify reapeats in (x, y, z) direction,
            eg. (3, 3, 1)
        ts_dir : str
            a path to directory with TSs
            e.g. 'TS_estimate'
        yamlfile : str
            a name of the .yaml file with reaction list
        pseudopotentials : dict(str: str)
            a dictionary with QE pseudopotentials for all species.
            e.g.
            dict(Cu='Cu.pbe-spn-kjpaw_psl.1.0.0.UPF',
                H='H.pbe-kjpaw_psl.1.0.0.UPF',
                O='O.pbe-n-kjpaw_psl.1.0.0.UPF',
                C='C.pbe-n-kjpaw_psl.1.0.0.UPF',
                )
        pseudo_dir : str
            a path to the QE's pseudopotentials main directory
            e.g.
            '/home/mgierad/espresso/pseudo'

        '''
        self.facetpath = facetpath
        self.slab = slab
        self.repeats = repeats
        self.ts_dir = ts_dir
        self.yamlfile = yamlfile
        self.pseudopotentials = pseudopotentials
        self.pseudo_dir = pseudo_dir
        self.executable = executable
        self.balsam_exe_settings = balsam_exe_settings
        self.calc_keywords = calc_keywords

    def set_up_irc(self, pytemplate_f, pytemplate_r):
        ''' Set up IRC calculations

        Parameters
        __________
        pytemplate_f, pytemplate_r : python scripts
            python scripts templates for irc calculations
        '''

        ts_path = os.path.join(self.facetpath, self.ts_dir)
        ts_uq_dir = ts_path + '_unique'
        ts = TS(self.facetpath, self.slab, self.ts_dir,
                self.yamlfile, self.repeats)
        rxn = ts.get_rxn_name()
        unique_ts_index = ts.check_symm(ts_uq_dir)

        # load templates for irc_f and irc_r calculations
        with open(pytemplate_f, 'r') as f:
            template_f = f.read()
        with open(pytemplate_r, 'r') as r:
            template_r = r.read()

        for i, prefix in enumerate(unique_ts_index):
            prefix = prefix[1:]
            irc_dir = os.path.join(self.facetpath, 'IRC', prefix)
            irc_py_dir = os.path.join(self.facetpath, 'IRC')
            ts_file_name = os.path.join(prefix + '_' + rxn + '_ts')
            ts_file_name_xyz = os.path.join(prefix, ts_file_name + '.xyz')
            os.makedirs(irc_dir, exist_ok=True)
            src_ts_xyz_path = os.path.join(
                ts_uq_dir, prefix, prefix + '_' + rxn + '_ts_final.xyz')
            dest_ts_path = os.path.join(irc_dir, ts_file_name)
            try:
                write(dest_ts_path + '.xyz', read(src_ts_xyz_path))
                write(dest_ts_path + '.png', read(src_ts_xyz_path))
                IRC.create_job_files_irc(self, irc_py_dir, ts_file_name,
                                         template_f, '_irc_f.py', prefix, rxn,
                                         ts_file_name_xyz)
                IRC.create_job_files_irc(self, irc_py_dir, ts_file_name,
                                         template_r, '_irc_r.py', prefix, rxn,
                                         ts_file_name_xyz)
            except FileNotFoundError:
                # skip the file because calculation did not finished
                print('Calculations for {} probably did not finish'.format(
                    src_ts_xyz_path))
                print('Skipping...')
                pass
                # raise

    def create_job_files_irc(self, irc_py_dir, ts_file_name, template,
                             which_irc, prefix, rxn, ts_file_name_xyz):
        ''' Create python scripts to submit jobs

        Parameters:
        __________
        irc_py_dir : str
            a path to directory where .py scripts are to be saved
            e.g. 'Cu_111/IRC/'
        ts_file_name : str
            a file name of the TS (without 'xyz.')
            e.g. '01_CO2+H_CHO2_ts'
        template : python file
            a template for irc_f or irc_r calculations
        which_irc : str
            what to add at the end of .py file, before '.py'?
            e.g.
            '_irc_f.py'
        prefix : str
            a prefix for the given geometry
            e.g. '00'
        rxn : str
            name of the reaction
            e.g. 'O+H_OH'
        ts_file_name.xyz : str
            a path to ts
            e.g. '00/00_CO2+H_CHO2_ts.xyz'

        '''
        job_name = os.path.join(irc_py_dir, ts_file_name[:-3] + which_irc)
        with open(job_name, 'w') as f:
            f.write(template.format(
                prefix=prefix, rxn=rxn, TS_xyz=ts_file_name_xyz,
                pseudopotentials=self.pseudopotentials,
                pseudo_dir=self.pseudo_dir,
                executable=self.executable,
                balsam_exe_settings=self.balsam_exe_settings,
                calc_keywords=self.calc_keywords
            ))
        f.close()

    def prepare_opt_irc(self, struc_path, irc, traj, pytemplate_irc_opt):
        ''' Preapare files for IRC optimization

        Parameters
        __________
        struc_path :
            path to all sets of IRC calculations, eg. IRC/00, IRC/01 ...
        irc : str
            specify between irc_f and irc_r
        traj : ase trajectory file
            e.g *irc_f.traj from previous irc calculation
        pytemplate_irc_opt : python script
            template file for IRC optimization job

        '''
        irc_opt_pth = os.path.join(struc_path, irc + '_opt')
        os.makedirs(irc_opt_pth, exist_ok=True)
        trajPath = os.path.join(struc_path, traj)
        initXYZ = os.path.join(irc_opt_pth, traj[:-5])
        write(initXYZ + '.xyz', read(trajPath))
        write(initXYZ + '_initial.png', read(trajPath))
        geom = os.path.join(traj[:-4] + 'xyz')
        IRC.create_job_files(self, pytemplate_irc_opt, irc_opt_pth, traj, geom)

    def create_job_files(self, pytemplate_irc_opt, irc_opt_pth, traj, geom):
        ''' Create slurm files for IRC optimization

        Parameters
        __________
        pytemplate_irc_opt : python script
            template file for IRC optimization job
        irc_opt_pth : str
            directory where irc optimization will we placed,
            e.g Cu_111/IRC/00/irc_r_opt
        traj : ase trajectory file
            e.g *irc_f.traj from previous irc calculation
        geom : str
            a name of a .xyz file with the coordinates of the TS

        '''
        with open(pytemplate_irc_opt, 'r') as f:
            pytemplate_irc_opt = f.read()
            ts = TS(self.facetpath, self.slab, self.ts_dir,
                    self.yamlfile, self.repeats)
            rxn = ts.get_rxn_name()
            prefix = traj[:2]
            fname = os.path.join(irc_opt_pth, prefix + '_' +
                                 rxn + '_' + os.path.split(irc_opt_pth)[1]
                                 + '.py')
            with open(fname, 'w') as f:
                f.write(pytemplate_irc_opt.format(
                    geom=geom, rxn=rxn, prefix=prefix,
                    pseudopotentials=self.pseudopotentials,
                    pseudo_dir=self.pseudo_dir,
                    executable=self.executable,
                    balsam_exe_settings=self.balsam_exe_settings,
                    calc_keywords=self.calc_keywords))
            f.close()
        f.close()

    def opt_after_IRC(self, irc_dir, pytemplate_irc_opt):
        ''' Set up optimization to minimas after IRC calculation

        Parameters
        __________
        facetpath : str
            a path to main directory
            e.g. Cu_111/
        irc_dir : str
            a path to main IRC directory
        pytemplate_irc_opt : python script
            slurm template for IRC optimization

        The function checks if *traj files containing IRC trajectories exists
        and have some content. If so, a minimum optimization will be set up.
        Otherwise, the given geometry is skipped.

        '''
        irc_path = os.path.join(self.facetpath, irc_dir)
        for struc in sorted(os.listdir(irc_path), key=str):
            struc_path = os.path.join(irc_path, struc)
            if os.path.isdir(struc_path):
                for traj in os.listdir(struc_path):
                    if traj.endswith('irc_f.traj'):
                        trajPath = os.path.join(struc_path, traj)
                        if os.stat(trajPath).st_size == 0:
                            pass
                        else:
                            try:
                                IRC.prepare_opt_irc(
                                    self, struc_path, 'irc_f', traj,
                                    pytemplate_irc_opt)
                            except FileNotFoundError:
                                # Error handling to be developed
                                raise
                    elif traj.endswith('irc_r.traj'):
                        trajPath = os.path.join(struc_path, traj)
                        if os.stat(trajPath).st_size == 0:
                            pass
                        else:
                            try:
                                IRC.prepare_opt_irc(
                                    self, struc_path, 'irc_r', traj,
                                    pytemplate_irc_opt)
                            except FileNotFoundError:
                                # Error handling to be developed
                                raise
