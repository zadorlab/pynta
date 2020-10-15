import unittest
import os


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


if __name__ == '__main__':
    unittest.main()
