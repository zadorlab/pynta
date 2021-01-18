#!/usr/bin/env python3

from setuptools import setup

with open('README.md', 'r') as f:
    long_description = f.read()
with open('requirements.txt', 'r') as f:
    requirements = f.read().split('\n')

setup(name='pynta',
      version='1.0.0',
      summary='A workflow for automated reaction path exploration on '
      'metallic surfaces.',
      author='Maciej Gierada, Eric D. Hermes, David H. Bross',
      author_email='maciek.gierada@gmai.com ehermes@sandia.gov',
      long_description=long_description,
      long_description_type='text/markdown',
      license='GPL-3.0',
      include_package_data=True,
      packages=['pynta',
                'pynta.pytemplate',
                'pynta.jobtemplate',
                'pynta.excatkit',
                'pynta.license'],
      python_requires='>=3.6, <4',
      install_requires=requirements,
      classifiers=[
          'Development Status :: 4 - Beta',
          'Intended Audience :: Developers',
          'Intended Audience :: Science/Research',
          'Topic :: Scientific/Engineering :: Chemistry',
          'Topic :: Scientific/Engineering :: Mathematics',
          'Topic :: Scientific/Engineering :: Physics'
          'Natural Language :: English',
          'Operating System :: OS Independent'
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7',
      ],
      )
