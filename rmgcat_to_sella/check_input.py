import yaml
import os


class YamlSplitter():
    ''' A class to handle input preparations - Currently only .yaml file '''

    def __init__(self,
                 big_yaml_file):
        '''
        Parameters:
        ___________

        big_yaml_file : str
            path to the yaml file with all reactions

        '''
        self.big_yaml_file = big_yaml_file

    def split(self):
        ''' Split big_yaml_file that contains many reactions into
            single reaction yaml files'''

        # load all-reaction .yaml file
        all_rxns = self.open_big_yaml_file()
        # load .yaml template
        template = self.open_template_yaml_file()

        # go through all reactions and creates single reaction .yaml files
        for rxn in all_rxns:
            index = rxn['index']
            reaction = rxn['reaction']
            reaction_family = rxn['reaction_family']
            reactant = self.fix_format(rxn['reactant'].split('\n'))
            product = self.fix_format(rxn['product'].split('\n'))

            # name of the .yaml file to be created
            new_yaml = 'reaction_{}.yaml'.format(str(index).zfill(2))

            with open(new_yaml, 'w') as f:
                f.write(template.format(
                    index=index,
                    reaction=reaction,
                    reaction_family=reaction_family,
                    reactant=reactant,
                    product=product))

    def fix_format(self,
                   reactant):
        ''' Fix indentation of the reactant and product

        Parameters:
        ___________

        reactants : list(str)
            reactant info in a form of a list, where each line is an element
            of that list. No indentation
            e.g.
            [
            'multiplicity -187',
            '1 *1 O u0 p0 c0 {2,S} {4,S}[]',
            '2 *2 H u0 p0 c0 {1,S}',
            '3 *3 X u0 p0 c0 ',
            '4    X u0 p0 c0 {1,S}', ''
            ]


        Returns:
        ________

        fixed_format : str
            a string with fixed indentations
            e.g.
                multiplicity -187
                1 *1 O u0 p0 c0 {2,S} {4,S}[]
                2 *2 H u0 p0 c0 {1,S}
                3 *3 X u0 p0 c0
                4    X u0 p0 c0 {1,S}
        '''
        fixed_format = []

        # go through each line of the reactant and add 4 spaces
        # at the beginning
        for line in reactant:
            # oh! the multiplicity line has to be indentate differently
            if line.startswith('multiplicity'):
                fixed_format.append('   ' + line)
            else:
                fixed_format.append('        ' + line)

        # join every element of the list into a single string using new line
        # character as a separator
        separator = '\n '
        fixed_format = separator.join(fixed_format)
        return fixed_format

    def open_big_yaml_file(self):
        ''' Load the .yaml file with all reactions

        Returns:
        ________

        all_rxns : list[dict[str: str]]
            loaded yaml file

        '''
        with open(self.big_yaml_file, 'r') as f:
            yamltxt = f.read()
        all_rxns = yaml.safe_load(yamltxt)
        return all_rxns

    def open_template_yaml_file(self):
        ''' Load the .yaml file serving as a template

        Returns:
        ________

        template : list[dict[str: str]]
            loaded template yaml file that is used to generate one-reaction
            .yanl files

        '''
        path = os.path.abspath(__file__)
        dir_path = os.path.dirname(path)
        yaml_template = os.path.join(
            dir_path, 'yaml_template', 'yaml_template.py')

        with open(yaml_template, 'r') as f:
            template = f.read()
        return template
