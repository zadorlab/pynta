import os


class InputChecker():
    ''' A class to check for input files '''

    def __init__(self,
                 yamlfile,
                 slab,
                 inputR2S,
                 run_me):
        ''' Initialize

        Parameters:
        ___________

        yamlfile : str
            a name of the .yaml file with reaction list
        slab : str
            a '.xyz' file name with the optimized slab
            e.g.
            'Cu_100_slab_opt.xyz'
        inputR2S : python file
            an input file with paramters to the workflow
        run_me : python script
            the workflow execution script

        '''
        self.yamlfile = yamlfile
        self.slab = slab
        self.inputR2S = inputR2S
        self.run_me = run_me
        # Get the path to the working directory
        self.working_dir = os.getcwd()

    def check_all(self):
        ''' Print info about checking input files '''
        print('Checking input...')
        print('---')
        if not self.is_input_file():
            print('---')
            print('Error')
            print('Make sure all input files are '
                  'in your working directory')
        else:
            print('Passed')

    def is_input_file(self):
        ''' Check if there are input files in the working directory '''
        # create check list
        check_list = []

        # check inputs
        check_yaml = self.check_yaml()
        check_list.append(check_yaml)
        check_slab = self.check_slab()
        check_list.append(check_slab)
        check_inputR2S = self.check_inputR2S()
        check_list.append(check_inputR2S)
        check_run_me = self.check_run_me()
        check_list.append(check_run_me)
        
        # There is an error if at least one element of check_list is False
        if any(check_list):
            return False
        else:
            return True

    def check_yaml(self):
        if not os.path.isfile(self.yamlfile):
            print('!    .yaml file is not in '
                  'your working directory: \n{}'.format(self.working_dir))
            return False
        else:
            return True

    def check_slab(self):
        if not os.path.isfile(self.slab):
            print('!    .xyz file with optimized slab is not in '
                  'your working directory: \n{}'.format(self.working_dir))
            return False
        else:
            return True

    def check_inputR2S(self):
        if not os.path.isfile(self.inputR2S):
            print('!    inputR2S.py file is not in your current'
                  ' working directory: \n{}'.format(self.working_dir))
            return False
        else:
            return True

    def check_run_me(self):
        if not os.path.isfile(self.run_me):
            print('!    run_me.py file is not in your current'
                  ' working directory: \n{}'.format(self.working_dir))
            return False
        else:
            return True
