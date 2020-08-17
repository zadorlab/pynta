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

    def check(self):
        ''' Print info about checking input files '''
        print('Checking input...')
        print('---')
        if not self.is_input_file():
            print('---')
            print('Error')
            print('Make sure all input files are '
                  'in your working directory')
            exit()
        else:
            print('Passed')

    def is_input_file(self):
        ''' Check if there are input files in the working directory '''
        # Get the path to the working directory
        working_dir = os.getcwd()

        # check inputs
        if not os.path.isfile(self.yamlfile):
            print('!    .yaml file is not in '
                  'your working directory: \n{}'.format(working_dir))
            return False
        elif not os.path.isfile(self.slab):
            print('!    .xyz file with optimized slab is not in '
                  'your working directory: \n{}'.format(working_dir))
            return False
        elif not os.path.isfile(self.inputR2S):
            print('!    inputR2S.py file is not in '
                  'your current working directory: \n{}'.format(working_dir))
            return False
        elif not os.path.isfile(self.run_me):
            print('!    run_me.py file is not in '
                  'your current working directory: \n{}'.format(working_dir))
            return False
        else:
            return True
