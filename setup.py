#!/usr/bin/env python3

from setuptools import setup, find_packages

with open('README.md', 'r') as f:
    long_description = f.read()


setup(name='pynta',
      version='1.0.0',
      summary='A workflow for automated reaction path exploration on '
      'metallic surfaces.',
      author='Maciej Gierada, Eric D. Hermes, David Bross',
      author_email='mgierad@sandia.gov. ehermes@sandia.gov',
      long_description=long_description,
      long_description_type='text/markdown',
      #   packages=find_packages(),
      packages=['pynta',
                'pynta.pytemplate',
                'pynta.jobtemplate',
                'pynta.excatkit'],
      python_requires='>=3.6',
      install_requires=['numpy', 'ase', 'spglib',
                        'matplotlib<3.2', 'networkx<2.4'],
      )
