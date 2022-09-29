from ase.utils.structure_comparator import SymmetryEquivalenceCheck
from ase.io import write, read
from copy import deepcopy

def get_unique_sym(geoms):
    ''' Check for the symmetry equivalent structures in the given files

    Parameters
    ___________
    geoms: list of paths to .xyz or .traj files to compare

    Returns
    ________
    idx_list : list(str)
        a list with prefixes of all symmetrically distinct sites

    '''
    comparator = SymmetryEquivalenceCheck()

    good_adsorbates_atom_obj_list = []
    geos_out = []

    for geom in geoms:
        adsorbate_atom_obj = read(geom)
        adsorbate_atom_obj.pbc = True
        comparision = comparator.compare(
            adsorbate_atom_obj, good_adsorbates_atom_obj_list)

        if comparision is False:
            good_adsorbates_atom_obj_list.append(adsorbate_atom_obj)
            geos_out.append(geom)

    return geos_out

def get_unique_sym_indices(geoms):
    ''' Check for the symmetry equivalent structures in the given files

    Parameters
    ___________
    geoms: list of paths to .xyz or .traj files to compare

    '''
    comparator = SymmetryEquivalenceCheck()

    good_adsorbates_atom_obj_list = []
    geos_out = []

    for geom in geoms:
        adsorbate_atom_obj = read(geom)
        adsorbate_atom_obj.pbc = True
        comparision = comparator.compare(
            adsorbate_atom_obj, good_adsorbates_atom_obj_list)

        if comparision is False:
            good_adsorbates_atom_obj_list.append(adsorbate_atom_obj)
            geos_out.append(geom)

    indices = [geoms.index(g) for g in geos_out]

    return indices

def get_unique_sym_structs(geoms):
    ''' Check for the symmetry equivalent structures in the given files

    Parameters
    ___________
    geoms: list of Atoms objects to compare


    '''
    comparator = SymmetryEquivalenceCheck()

    geoms_copy = deepcopy(geoms)

    good_adsorbates_atom_obj_list = []
    geos_out = []

    for i,geom in enumerate(geoms_copy):
        adsorbate_atom_obj = geom
        adsorbate_atom_obj.pbc = True
        comparision = comparator.compare(
            adsorbate_atom_obj, good_adsorbates_atom_obj_list)

        if comparision is False:
            good_adsorbates_atom_obj_list.append(adsorbate_atom_obj)
            geos_out.append(geoms[i])

    return geos_out

def get_unique_sym_struct_indices(geoms):
    ''' Check for the symmetry equivalent structures in the given files

    Parameters
    ___________
    geoms: list of Atoms objects to compare


    '''
    comparator = SymmetryEquivalenceCheck()

    geoms_copy = deepcopy(geoms)

    good_adsorbates_atom_obj_list = []
    geos_out = []

    for i,geom in enumerate(geoms_copy):
        adsorbate_atom_obj = geom
        adsorbate_atom_obj.pbc = True
        comparision = comparator.compare(
            adsorbate_atom_obj, good_adsorbates_atom_obj_list)

        if comparision is False:
            good_adsorbates_atom_obj_list.append(adsorbate_atom_obj)
            geos_out.append(geoms[i])

    indices = [geoms_copy.index(g) for g in geos_out]

    return indices

def filter_nonunique_TS_guess_indices(geoms,Es):
    ''' Check for the symmetry equivalent structures in the given files

    Parameters
    ___________
    geoms: list of paths to .xyz or .traj files to compare

    '''
    comparator = SymmetryEquivalenceCheck()

    good_adsorbates_atom_obj_list = []
    geos_out = []
    Esout = []

    for j,geom in enumerate(geoms):
        adsorbate_atom_obj = read(geom)
        adsorbate_atom_obj.pbc = True
        for i,good_adsorbate in enumerate(good_adsorbates_atom_obj_list):
            comparison = comparator.compare(adsorbate_atom_obj,good_adsorbate)
            if comparison and Es[j] < Esout[i]:
                geos_out[i] = geom
                good_adsorbates_atom_obj_list[i] = adsorbate_atom_obj
                Esout[i] = Es[j]
                break
            elif comparison:
                break
        else:
            good_adsorbates_atom_obj_list.append(adsorbate_atom_obj)
            geos_out.append(geom)
            Esout.append(Es[j])

    return geos_out,Esout
