#!/usr/bin/env python3

from setuptools import setup, find_packages

with open('README.md', 'r') as f:
    long_description = f.read()


setup(name='rmgcat_to_sella',
      version='0.0.1',
      author='Eric Hermes',
      author_email='ehermes@sandia.gov',
      long_description=long_description,
      long_description_type='text/markdown',
      packages=find_packages(),
      python_requires='>=3.6',
      install_requires=['numpy', 'ase', 'catkit', 'spglib'],
      )
