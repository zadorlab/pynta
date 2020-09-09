import os

from pathlib import Path

from ase.io import read, write
from rmgcat_to_sella.ts import TS
from rmgcat_to_sella.io import IO


class IRC():
    def __init__(
            self,
            facetpath,
            slab,
            repeats,
            ts_estimate_dir,
            yamlfile,
            pseudopotentials,
            pseudo_dir,
            balsam_exe_settings,
            calc_keywords,
            creation_dir
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
        ts_estimate_dir : str
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
        self.ts_estimate_dir = ts_estimate_dir
        self.yamlfile = yamlfile
        self.pseudopotentials = pseudopotentials
        self.pseudo_dir = pseudo_dir
        self.balsam_exe_settings = balsam_exe_settings
        self.calc_keywords = calc_keywords
        self.creation_dir = creation_dir
        self.io = IO()

    def set_up_irc(
            self,
            rxn,
            pytemplate_f,
            pytemplate_r):
        ''' Set up IRC calculations

        Parameters:
        __________
        rxn : dict(yaml[str:str])
            a dictionary with info about the paricular reaction. This can be
            view as a splitted many reaction .yaml file into a single reaction
            .yaml file
        pytemplate_f, pytemplate_r : python scripts
            python scripts templates for irc calculations

        '''
        rxn_name = self.io.get_rxn_name(rxn)
        ts_path = os.path.join(self.facetpath, rxn_name, self.ts_estimate_dir)
        ts_uq_dir = ts_path + '_unique'
        ts = TS(self.facetpath, self.slab, self.ts_estimate_dir,
                self.yamlfile, self.repeats, self.creation_dir)
        unique_ts_index = ts.check_symm(ts_uq_dir)

        # load templates for irc_f and irc_r calculations
        with open(pytemplate_f, 'r') as f:
            template_f = f.read()
        with open(pytemplate_r, 'r') as r:
            template_r = r.read()

        for i, prefix in enumerate(unique_ts_index):
            prefix = prefix[1:]
            irc_dir_prefix = os.path.join(
                self.facetpath, rxn_name, 'IRC', prefix)
            irc_dir, _ = os.path.split(irc_dir_prefix)
            ts_file_name = os.path.join(prefix + '_' + rxn_name + '_ts')
            ts_file_name_xyz = os.path.join(prefix, ts_file_name + '.xyz')
            os.makedirs(irc_dir_prefix, exist_ok=True)
            src_ts_xyz_path = os.path.join(
                ts_uq_dir, prefix, prefix + '_' + rxn_name + '_ts_final.xyz')
            dest_ts_path = os.path.join(irc_dir_prefix, ts_file_name)
            try:
                write(dest_ts_path + '.xyz', read(src_ts_xyz_path))
                write(dest_ts_path + '.png', read(src_ts_xyz_path))
                self.create_job_files_irc(irc_dir,
                                          ts_file_name,
                                          template_f,
                                          '_irc_f.py',
                                          prefix,
                                          rxn_name,
                                          ts_file_name_xyz)
                self.create_job_files_irc(irc_dir,
                                          ts_file_name,
                                          template_r,
                                          '_irc_r.py',
                                          prefix,
                                          rxn_name,
                                          ts_file_name_xyz)
            except FileNotFoundError:
                # skip the file because calculation did not finished
                print('Calculations for {} probably did not finish'.format(
                    src_ts_xyz_path))
                print('Skipping...')
                pass

    def create_job_files_irc(
            self,
            irc_dir,
            ts_file_name,
            template,
            which_irc,
            prefix,
            rxn_name,
            ts_file_name_xyz):
        ''' Create python scripts to submit jobs

        Parameters:
        __________
        irc_dir : str
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
        rxn_name : str
            name of the reaction
            e.g. 'O+H_OH'
        ts_file_name.xyz : str
            a path to ts
            e.g. '00/00_CO2+H_CHO2_ts.xyz'

        '''
        job_name = os.path.join(irc_dir, ts_file_name[:-3] + which_irc)
        with open(job_name, 'w') as f:
            f.write(template.format(
                prefix=prefix, rxn=rxn_name, ts_xyz=ts_file_name_xyz,
                pseudopotentials=self.pseudopotentials,
                pseudo_dir=self.pseudo_dir,
                balsam_exe_settings=self.balsam_exe_settings,
                calc_keywords=self.calc_keywords
            ))

    def opt_after_IRC(
            self,
            pytemplate_irc_opt):
        ''' Create opt IRC calculations for reactions

        Parameters:
        ___________

        pytemplate_irc_opt : python script
            template file for IRC optimization job

        '''

        reactions = self.io.open_yaml_file(self.yamlfile)
        for rxn in reactions:
            try:
                self.check_irc_finished(rxn, pytemplate_irc_opt)
            except FileNotFoundError:
                irc_error = irc_path.split('/')[1]
                print('IRC calculations for {} '
                      'did not finished.'.format(irc_error))
                print('Skipping...')
                pass

    def check_irc_finished(
            self,
            rxn,
            pytemplate_irc_opt):
        ''' Preset opt irc calculations

        Parameter:
        __________
        rxn : dict(yaml[str:str])
            a dictionary with info about the paricular reaction. This can be
            view as a splitted many reaction .yaml file into a single reaction
            .yaml file
        pytemplate_irc_opt : python script
            template file for IRC optimization job

        '''
        rxn_name = self.io.get_rxn_name(rxn)
        irc_path = os.path.join(self.facetpath, rxn_name, 'IRC')
        for irc_name in ['irc_f', 'irc_r']:
            irc_traj_paths = Path(irc_path).glob(
                '**/*{}.traj'.format(irc_name))
            for irc_traj_path in irc_traj_paths:
                print(irc_traj_path)
                # if file is empty, there is somethig wrong with calculation
                # skip at this moment
                if os.stat(irc_traj_path).st_size == 0:
                    pass
                else:
                    self.prepare_opt_irc(irc_traj_path, irc_name,
                                         pytemplate_irc_opt, rxn_name)

    def prepare_opt_irc(
            self,
            irc_traj_path,
            irc_name,
            pytemplate_irc_opt,
            rxn_name):
        ''' Preapare files for IRC optimization

        Parameters
        __________
        irc_traj_path :
            get path to IRC trajectory
        irc_name : str
            specify 'irc_f' or 'irc_r' depending which type of calculations to
            set up
        pytemplate_irc_opt : python script
            template file for IRC optimization job
        rxn_name : str
            The name of the reaction in the following format:
            'OH_H+O'

        '''
        # get path to set of IRC calculations, eg. ../IRC/00, ../IRC/01
        struc_path, irc_traj = os.path.split(irc_traj_path)
        irc_opt_path = os.path.join(struc_path, irc_name + '_opt')
        os.makedirs(irc_opt_path, exist_ok=True)

        init_xyz = os.path.join(irc_traj_path, irc_traj[:-5])
        write(init_xyz + '.xyz', read(irc_traj_path))
        write(init_xyz + '_initial.png', read(irc_traj_path))

        xyz_geom_file = os.path.jpin(irc_traj_path[:-4] + 'xyz')
        self.create_job_files(rxn_name,
                              pytemplate_irc_opt,
                              irc_opt_path,
                              irc_traj_path,
                              xyz_geom_file)

    def create_job_files(
            self,
            rxn_name,
            pytemplate_irc_opt,
            irc_opt_path,
            irc_traj_path,
            xyz_geom_file):
        ''' Create slurm files for IRC optimization

        Parameters
        __________
        rxn_name : str
            The name of the reaction in the following format:
            'OH_H+O'
        pytemplate_irc_opt : python script
            template file for IRC optimization job
        irc_opt_path : str
            directory where irc optimization will we placed,
            e.g Cu_111/IRC/00/irc_r_opt
        irc_traj_path : ase trajectory file
            e.g *irc_f.traj from previous irc calculation
        xyz_geom_file : str
            a name of a .xyz file with the coordinates of the TS

        '''
        with open(pytemplate_irc_opt, 'r') as f:
            pytemplate_irc_opt = f.read()
            prefix = irc_traj_path[:2]
            fname = os.path.join(irc_opt_path, prefix + '_'
                                 + rxn_name + '_' +
                                 os.path.split(irc_opt_path)[1]
                                 + '.py')
            with open(fname, 'w') as f:
                f.write(pytemplate_irc_opt.format(
                    geom=xyz_geom_file, rxn=rxn_name, prefix=prefix,
                    pseudopotentials=self.pseudopotentials,
                    pseudo_dir=self.pseudo_dir,
                    balsam_exe_settings=self.balsam_exe_settings,
                    calc_keywords=self.calc_keywords))
