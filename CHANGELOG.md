# Changelog

## 0.0.1

- Develpmet version of Pynta supporting SLURM queue system

## 1.0.0

- Initial Pynta release

  - Added [NWChem](https://www.nwchem-sw.org) and [Quantum Espresso](https://www.quantum-espresso.org) calculators.

  - Added support for monodentate adsorption.

  - Added support for diatomic reactions (two framgents of a molecule or two molecules/atoms react with each other). e.g

    - CH --> C + H
    - H<sub>2</sub>O --> H + OH
    - HCOOCH<sub>3</sub> --> HCO + CH<sub>3</sub>O

  - Supported facets:

    - FCC111
    - FCC100
    - FCC211
    - BCC111
    - BCC100
    - HCP0001
    - diamond111
    - diamond100

  - Added support for all metalic surfaces, provided that there is a pseudopotential implemented for a desired metal in either NWChem of Quantum Espresso and this metal is supported by ASE.

  - Added support for Balsam DB. This allows a user to define only one job managing all Pynta subtasks, which can be later on submitted to any queue system on HPC machinses and will be seen as a one job.

## 1.0.1

- Minor bug fixes

## 1.1.0

- added an automated logic to specify which adsorbate atom connects to the surface, rather that specifying it by hand
- fixed atoms numbering during 02 step input generation that led to specifying wrong atomic indicies for a penalty function calculation for some molecules
- added a new method `Restart().remove_empty_pickle_files()` that removes all empty `*.pckl` files from unfinished calculations
- added a new method `Restart().remove_error_files()` that removes all temp files from unfinished calculations
- added a new method
- rewritten part where various balsam calculators are imported
- added Dockerfile
- new graphics
- updated docs
- setup.cfg added and setup.py refactored - getting ready for PyPI released
- fixed requirements.txt
- other bug fixes and improvements
