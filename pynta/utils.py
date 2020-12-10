#!/usr/bin/env python
# by default, each atom is bonded to the surface through the atom
# at index 0. This can be overwrite by using this dict
edge_cases_bonded_dict = dict(CO3=1, CH3O=1, CH3OH=1)

# by default, the first topology structure is used to generate
# adsorbates. This can be modified using this dict
edge_cases_topology_dict = dict(COOH=1, HCOOH=1, CH3O2=1, HCOOCH3=17)
