#!/usr/bin/env python3

from setuptools import setup
from setuptools.config import read_configuration
from os import path

with open('README.md', 'r') as f:
    long_description = f.read()
with open('requirements.txt', 'r') as f:
    requirements = f.read().split('\n')

current_dir = path.dirname(__file__)
conf_dict = read_configuration(path.join(current_dir, 'setup.cfg'))

setup(install_requires=requirements,
      include_package_data=True,

      )

# setup(name='pynta',
#       version='1.0.1',
#       summary='A workflow for automated reaction path exploration on metallic surfaces.',
#       homepage='https://zadorlab.github.io/pynta/',
#       author='Maciej Gierada, Eric D. Hermes, David H. Bross',
#       author_email='maciek.gierada@gmai.com ehermes@sandia.gov',
#       long_description=long_description,
#       long_description_type='text/markdown',
#       license='GPL-3.0',
#       include_package_data=True,
#       packages=['pynta',
#                 'pynta.pytemplate',
#                 'pynta.jobtemplate',
#                 'pynta.excatkit',
#                 'pynta.license'],
#       python_requires='>=3.6, <4',
#       install_requires=requirements,
#       classifiers=[
#           'Development Status :: 4 - Beta',
#           'Intended Audience :: Developers',
#           'Intended Audience :: Science/Research',
#           'Topic :: Scientific/Engineering :: Chemistry',
#           'Topic :: Scientific/Engineering :: Mathematics',
#           'Topic :: Scientific/Engineering :: Physics',
#           'Natural Language :: English',
#           'Operating System :: OS Independent',
#           'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
#           'Programming Language :: Python :: 3.6',
#           'Programming Language :: Python :: 3.7',
#           'Programming Language :: Python :: 3.8',
#       ],
#       )
