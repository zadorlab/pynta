![Build status](https://github.com/zadorlab/pynta/workflows/Build%20status/badge.svg)
![GitHub tag (latest by date)](https://img.shields.io/github/v/tag/zadorlab/pynta?label=version)
![GitHub last commit](https://img.shields.io/github/last-commit/zadorlab/pynta?label=last%20modified)
![GitHub commits since tagged version](https://img.shields.io/github/commits-since/zadorlab/pynta/v1.1.0?color=9cf)
![GitHub](https://img.shields.io/badge/License-GPLv3-orange)
![GitHub top language](https://img.shields.io/github/languages/top/zadorlab/pynta?color=brightgreen)
[![Total alerts](https://img.shields.io/lgtm/alerts/g/zadorlab/pynta.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/zadorlab/pynta/alerts/)
![GitHub all releases](https://img.shields.io/github/downloads/zadorlab/pynta/total?color=yellow)

# Pynta

![image](https://github.com/zadorlab/pynta/assets/40675474/3a55408a-ab09-446e-8efb-a1fbf73880c1)

``pynta`` is an automated workflow code to enable the calculation of thermochemistry
and rate coefficients for reactions involving metallic surfaces.

The ``pynta`` code is designed to automatically characterize chemical reactions
relevant to heterogeneous catalysis. In particular, it spawns and processes a
large number of ab initio quantum chemistry calculations to study gas-phase
reactions on crystal facets. It is designed to run both on local clusters and
on petascale and upcoming exascale machines.

Pynta generates initial guesses for all gas phase species and adsorbates on the surface accounting
for the symmetry of the surface. These guesses are then optimized and frequencies are calculated for
the unique resulting geometries. Based on the provided reaction semi-atom mappings the
adsorbate are placed on the surface and perturbed in many unique ways to generate
transition state guesses, the best guesses are optimized and vibrational frequencies
and IRCs are generated for unique successfully optimized transition states. This
information can then be processed to generate thermochemistry and rate coefficients
for the adsorbates and reactions that can be used in microkinetic models.

Pynta is designed to work with the workflow code
[Fireworks](https://materialsproject.github.io/fireworks/), which provides users
with significant flexiblity in running Pynta.




# Pynta Documentation

Instructions on how to install and run Pynta are available in the documentation [here](https://zadorlab.github.io/pynta/).

# Acknowledgments

Please cite our paper if you are using Pynta:

Matthew S. Johnson, Maciej Gierada, Eric D. Hermes, David H. Bross, Khachik Sargsyan, Habib N. Najm, and Judit Zádor
Journal of Chemical Information and Modeling 2023 63 (16), 5153-5168 https://doi.org/10.1021/acs.jcim.3c00948 

```
@article{Johnson2023,
   author = {Johnson, Matthew S. and Gierada, Maciej and Hermes, Eric D. and Bross, David H. and Sargsyan, Khachik and Najm, Habib N. and Zádor, Judit},
   title = {Pynta─An Automated Workflow for Calculation of Surface and Gas–Surface Kinetics},
   journal = {Journal of Chemical Information and Modeling},
   volume = {63},
   number = {16},
   pages = {5153-5168},
   note = {doi: 10.1021/acs.jcim.3c00948},
   ISSN = {1549-9596},
   DOI = {10.1021/acs.jcim.3c00948},
   url = {https://doi.org/10.1021/acs.jcim.3c00948},
   year = {2023},
   type = {Journal Article}
}
```

# Copyright Notice

For five (5) years from 1/19/2021 the United States Government is granted for
itself and others acting on its behalf a paid-up, nonexclusive, irrevocable
worldwide license in this data to reproduce, prepare derivative works, and
perform publicly and display publicly, by or on behalf of the Government. There
is provision for the possible extension of the term of this license. Subsequent
to that period or any extension granted, the United States Government is
granted for itself and others acting on its behalf a paid-up, nonexclusive,
irrevocable worldwide license in this data to reproduce, prepare derivative
works, distribute copies to the public, perform publicly and display publicly,
and to permit others to do so. The specific term of the license can be
identified by inquiry made to National Technology and Engineering Solutions of
Sandia, LLC or DOE.

NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT OF
ENERGY, NOR NATIONAL TECHNOLOGY AND ENGINEERING SOLUTIONS OF SANDIA, LLC, NOR
ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY
LEGAL RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS
USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

Any licensee of this software has the obligation and responsibility to abide by
the applicable export control laws, regulations, and general prohibitions
relating to the export of technical data. Failure to obtain an export control
license or other authority from the Government may result in criminal liability
under U.S. laws.

# License Notice

Copyright 2021 National Technology & Engineering Solutions of Sandia,
LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the
U.S. Government retains certain rights in this software.
