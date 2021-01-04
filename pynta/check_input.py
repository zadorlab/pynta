import os
import sys


class InputChecker():
    ''' A class to check for input files '''

    def __init__(
            self,
            yamlfile: str,
            inputR2S: str,
            run_me_py: str,
            run_me_sh: str,
            working_dir: str) -> None:
        ''' Initialize

        Parameters
        ----------

        yamlfile : str
            a name of the .yaml file with reaction list
        slab : str
            a :literal:`*.xyz` file name with the optimized slab
            e.g.
            ``'Cu_100_slab_opt.xyz'``
        inputR2S : python file
            an input file with paramters to the workflow
        run_me : str
            the workflow execution script
            e.g.
            ``'run_me.py'``
        working_dir : str
            a main dir where all files are stored

        '''
        self.yamlfile = yamlfile
        self.inputR2S = inputR2S
        self.run_me_py = run_me_py
        self.run_me_sh = run_me_sh
        self.working_dir = working_dir

    def check_all(self) -> None:
        ''' Print info about checking input files '''
        print('Checking input...')
        print('---')
        if not self.is_input_file():
            print('---')
            print('Error')
            print('Make sure all input files are '
                  'in your working directory')
            print('---')
            # exit if at least one error
            sys.exit()
        else:
            print('Passed!')
            print('---')

    def is_input_file(self) -> bool:
        ''' Check if there are input files in the working directory '''
        # create check list
        check_list = []

        # check inputs
        check_yaml = self.check_yaml()
        check_list.append(check_yaml)
        check_inputR2S = self.check_inputR2S()
        check_list.append(check_inputR2S)
        check_run_me_py = self.check_run_me_py()
        check_list.append(check_run_me_py)
        check_run_me_sh = self.check_run_me_sh()
        check_list.append(check_run_me_sh)

        # There is an error if at least one element of check_list is False
        if all(check_list):
            return True
        else:
            return False

    def check_yaml(self) -> bool:
        ''' Check for :literal:`*.yaml` file with list of reactions '''
        if not os.path.isfile(self.yamlfile):
            print('!    .yaml file ({}) is not in your working '
                  'directory: \n{}'.format(self.yamlfile, self.working_dir))
            return False
        else:
            return True

    def check_slab(self) -> bool:
        ''' Check for slab :literal:`*.xyz` file with slab geometry '''
        if not os.path.isfile(self.slab):
            print('!    .xyz file ({}) with optimized slab is not in your '
                  'working directory: '
                  '\n{}'.format(self.slab, self.working_dir))
            return False
        else:
            return True

    def check_inputR2S(self) -> bool:
        ''' Check for inputR2S file '''
        if not os.path.isfile(self.inputR2S):
            print('!    inputR2S.py file ({}) is not in your current working '
                  'directory: \n{}'.format(self.inputR2S, self.working_dir))
            return False
        else:
            return True

    def check_run_me_py(self) -> bool:
        ''' Check for ``'run_me.py'`` file '''
        if not os.path.isfile(self.run_me_py):
            print('!    run_me.py file ({}) is not in your current working '
                  'directory: \n{}'.format(self.run_me_py, self.working_dir))
            return False
        else:
            return True

    def check_run_me_sh(self) -> bool:
        ''' Check for ``'run_me.sh'`` file '''
        if not os.path.isfile(self.run_me_sh):
            print('!    run_me.sh file ({}) is not in your current working '
                  'directory: \n{}'.format(self.run_me_sh, self.working_dir))
            return False
        else:
            return True
