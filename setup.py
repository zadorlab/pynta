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
