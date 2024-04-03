import shutil
import os
import ase
import sys
from ase.utils.structure_comparator import SymmetryEquivalenceCheck
from ase.io import write, read
import ase.constraints
from copy import deepcopy
from importlib import import_module
import numpy as np
import copy


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
        adsorbate_atom_obj.set_constraint() # Reset constraints for comparison
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
        adsorbate_atom_obj.set_constraint() # Reset constraints before comparison
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
        adsorbate_atom_obj.set_constraint() # Reset constraints before comparison
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
    indices = []

    for i,geom in enumerate(geoms_copy):
        adsorbate_atom_obj = geom
        adsorbate_atom_obj.pbc = True
        adsorbate_atom_obj.set_constraint() # Reset constraints before comparison
        comparision = comparator.compare(
            adsorbate_atom_obj, good_adsorbates_atom_obj_list)

        if comparision is False:
            good_adsorbates_atom_obj_list.append(adsorbate_atom_obj)
            indices.append(i)

    return indices

def get_unique_sym_struct_index_clusters(geoms):
    ''' Check for the symmetry equivalent structures in the given files

    Parameters
    ___________
    geoms: list of Atoms objects to compare


    '''
    comparator = SymmetryEquivalenceCheck()

    geoms_copy = deepcopy(geoms)

    good_adsorbates_atom_obj_list = []
    indices = []

    for i,geom in enumerate(geoms_copy):
        adsorbate_atom_obj = geom
        adsorbate_atom_obj.pbc = True
        adsorbate_atom_obj.set_constraint() # Reset constraints before comparison
        comparison = None
        for j,adlist in enumerate(good_adsorbates_atom_obj_list):
            comparison = comparator.compare(adsorbate_atom_obj, [adlist[0]])
            ind = j
            if comparison:
                break
        else:
            comparison = False

        if comparison is False:
            good_adsorbates_atom_obj_list.append([adsorbate_atom_obj])
            indices.append([i])
        else:
            indices[ind].append(i)

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
        adsorbate_atom_obj.set_constraint() # Reset constraints before comparison
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

def get_fmax(at):
    return np.max(np.abs([np.linalg.norm(at.get_forces()[i,:]) for i in range(at.get_forces().shape[0])]))

def name_to_ase_software(software_name):
    """
    go from software_name to the associated
    ASE calculator constructor
    """
    if software_name == "XTB":
        module = import_module("xtb.ase.calculator")
        return getattr(module, software_name)
    else:
        module = import_module("ase.calculators."+software_name.lower())
        return getattr(module, software_name)

def name_to_ase_opt(opt_name):
    """
    go from the optimizer name to the
    ASE optimizer
    """
    module = import_module("ase.optimize")
    return getattr(module, opt_name)

def clean_pynta_path(path,save_initial_guess=True):
    assert save_initial_guess

    for p in os.listdir(path):
        if p[:2] == "TS": #delete TSs
            shutil.rmtree(os.path.join(path,p))
        elif p == "Adsorbates":
            for ad in os.listdir(os.path.join(path,p)):
                if ad == ".DS_Store":
                    os.remove(os.path.join(path,p,ad))
                    continue
                for ind in os.listdir(os.path.join(path,p,ad)):
                    if ind == "info.json":
                        continue
                    elif ind.isdigit():
                        for file in os.listdir(os.path.join(path,p,ad,ind)):
                            if not "_init.xyz" in file:
                                pa = os.path.join(path,p,ad,ind,file)
                                if os.path.isdir(pa):
                                    shutil.rmtree(pa)
                                else:
                                    os.remove(pa)
                    else:
                        os.remove(os.path.join(path,p,ad,ind))

def construct_constraint(d):
    """
    construct a constrain from a dictionary that is the input to the constraint
    constructor plus an additional "type" key that indices the name of the constraint
    returns the constraint
    """
    constraint_dict = copy.deepcopy(d)
    constructor = getattr(ase.constraints,constraint_dict["type"])
    del constraint_dict["type"]
    return constructor(**constraint_dict)

def copyDataAndSave(origin, destination, file):
    '''
		Function to copy a file from a origin subdirectory to a destination subdirectory.
        If the file already exists, a "_BK" string is added to the name that it corresponds
        to a backup copy. This function needs another function 'check_copy', which exchanges
        the files, so that the most current file always does not have the _BK string,
        in this way the files are not affected when rerun the workflow.

        Parameters
        ___________
        origin: str
        destination: str
        file : str

    '''
    src = os.path.join(origin, file)
    dst = os.path.join(destination, file)

    count = 1

    while os.path.exists(dst):
        base, ext = os.path.splitext(file)
        new_name = f'{base}_BK{count}{ext}'
        dst = os.path.join(destination, new_name)
        count += 1

    shutil.copy(src, dst)

    check_copy(file, destination)

def check_copy(file2cpy, subdir):
    files = os.listdir(subdir)

    file_maxval = None
    number_maxval = -1

    for file in files:
        if file.startswith(file2cpy) and "_BK" in file:
            num = int(file.split("_BK")[1].split(".")[0])
            if num > number_maxval:
                number_maxval = num
                file_maxval = file

    if file_maxval:
        file_origin = os.path.join(subdir, file2cpy)
        file_bk = os.path.join(subdir, file_maxval)

        os.rename(file_origin, file_origin + "_tmp")
        os.rename(file_bk, file_origin)
        os.rename(file_origin + "_tmp", file_bk)