#!/usr/bin/env python3

from setuptools import setup, find_packages

with open('README.md', 'r') as f:
    long_description = f.read()


setup(name='rmgcat_to_sella',
      version='0.0.1',
      summary='Workflow for automatic search of 1st order saddle points (transition states)',
      author='Maciej Gierada, Eric Hermes',
      author_email='mgierad@sandia.gov. ehermes@sandia.gov',
      long_description=long_description,
      long_description_type='text/markdown',
      # packages=find_packages(),
      packages=['rmgcat_to_sella', 'rmgcat_to_sella.pytemplate','rmgcat_to_sella.jobtemplate'],
      python_requires='>=3.6',
      install_requires=['numpy', 'ase', 'catkit', 'spglib', 'matplotlib<3.2', 'networkx<2.4'],
      )
