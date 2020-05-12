import unittest
# import sys
import os
# import rmgcat_to_sella.ts
# import rmgcat_to_sella.irc
# import rmgcat_to_sella.main
# import rmgcat_to_sella.scan_bond


class Test(unittest.TestCase):

    def test_jobtemplate(self):
        path = os.path.abspath(__file__)
        dir_path = os.path.dirname(path)
        path_template = os.path.join(dir_path, 'jobtemplate/')
        self.assertTrue(os.path.exists(path_template))
        for i in os.listdir(path_template):
            jpath = os.path.join(path_template, i)
            self.assertTrue(os.path.exists(jpath))

    def test_pytemplate(self):
        path = os.path.abspath(__file__)
        dir_path = os.path.dirname(path)
        path_template = os.path.join(dir_path, 'pytemplate/')
        self.assertTrue(os.path.exists(path_template))
        for i in os.listdir(path_template):
            ppath = os.path.join(path_template, i)
            self.assertTrue(os.path.exists(ppath))

    # def test_ts_module(self):
    #     # print('Testing if all packages script can be loaded')
    #     package_scripts = ['ts', 'irc', 'scan_bond']
    #     for script in package_scripts:
    #         script = os.path.join('rmgcat_to_sella.' + script)
    #         if script not in sys.modules:
    #             return False
    #             print('OK')


if __name__ == '__main__':
    unittest.main()
