import shutil
import os
import ase
from ase import Atoms 
from molecule.quantity import ScalarQuantity
from ase.utils.structure_comparator import SymmetryEquivalenceCheck
from ase.io import write, read
import ase.constraints
from ase.geometry import get_distances
from copy import deepcopy
from importlib import import_module
import numpy as np
import copy
import json

def sites_match(site1,site2,slab,tol=0.5):
    """determine if two sites match

    Args:
        site1 (dict): site dictionary for first site
        site2 (dict): site dictionary for second site
        slab (ase.Atoms): the ase.Atoms object for the slab the sites are on
        tol (float, optional): Distance tolerance defaults to 0.5 Angstroms.

    Returns:
        _type_: _description_
    """
    v,dist = get_distances([site1["position"]], [site2["position"]], cell=slab.cell, pbc=slab.pbc)
    if dist > tol:
        return False
    elif site1["site"] != site2["site"]:
        return False
    elif site1["morphology"] != site2["morphology"]:
        return False
    else:
        return True
    
def get_occupied_sites(struct,sites,nslab,allowed_site_dict=dict(),site_bond_cutoff=2.5,
                       site_bond_disruption_cutoff=0.5,cutoff_corrections={"ontop": 0.3}):
    """determine what sites are occupied by what atoms

    Args:
        struct (ase.Atoms): Atoms object for the structure
        sites (list): list of site dictionaries
        nslab (int): number of atoms in the surface
        allowed_site_dict (dict, optional): dictionary mapping atom index to a list of allowed (site,morphology) for that atom
        site_bond_cutoff (float, optional): _description_. Defaults to 2.5.

    Raises:
        ValueError: _description_

    Returns:
        _type_: _description_
    """
    occ_sites = []
    for i in range(nslab,len(struct)):
        pos = struct.positions[i]
        if i in allowed_site_dict.keys():
            allowed_sites = allowed_site_dict[i]
        else:
            allowed_sites = None
        siteout = None
        mindist = None
        for site in sites:
            if allowed_sites and (site["site"],site["morphology"]) not in allowed_sites: #skip disallowed sites
                continue
            v,dist = get_distances([site["position"]], [pos], cell=struct.cell, pbc=struct.pbc)
            if siteout is None:
                siteout = site
                mindist = dist
                n = v
            else:
                if dist < mindist:
                    mindist = dist
                    siteout = site
                    n = v
        
        if mindist is None:
            #print(i)
            #view(struct)
            raise ValueError
        
        if (siteout["site"] not in cutoff_corrections.keys() and mindist < site_bond_cutoff) or (siteout["site"] in cutoff_corrections.keys() and mindist < site_bond_cutoff + cutoff_corrections[siteout["site"]]):
            mindn = None
            for j in range(nslab,len(struct)): #check for site bond disruption by other adsorbed atoms
                if i == j:
                    continue
                else:
                    AB,ABdist = get_distances([struct.positions[j]], [pos], cell=struct.cell, pbc=struct.pbc)
                    x = np.dot(AB[0,0,:],n[0,0,:])
                    if x < 0: #this means our target atom is closer than the other atom, which is okay (for this atom)
                        continue
                    dn = np.linalg.norm(AB[0,0,:] - x/mindist**2*n[0,0,:])
                    if mindn is None:
                        mindn = dn
                    elif dn < mindn:
                        mindn = dn
                    
                
            if mindn is None or mindn > site_bond_disruption_cutoff:
                s = deepcopy(siteout)
                s["bonding_index"] = i
                s["bond_length"] = mindist
                occ_sites.append(s)
    
    return occ_sites

class SiteOccupationException(Exception):
    pass

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
        adsorbate_atom_obj.set_constraint()
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
        adsorbate_atom_obj.set_constraint()
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
        adsorbate_atom_obj.set_constraint()
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
        adsorbate_atom_obj.set_constraint()
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
        adsorbate_atom_obj.set_constraint()
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
        adsorbate_atom_obj.set_constraint()
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
    if software_name == "TBLite" or software_name == "XTB":
        module = import_module("tblite.ase")
        return getattr(module, "TBLite")
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
    if save_initial_guess:
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
    else:
         for p in os.listdir(path):
            if p[:2] == "TS" or p == "Adsorbates": #delete TSs
                shutil.rmtree(os.path.join(path,p))
    
    if os.path.exists(os.path.join(path,"reaction_library")):
        shutil.rmtree(os.path.join(path,"reaction_library"))
    if os.path.exists(os.path.join(path,"thermo_library.py")):
        os.remove(os.path.join(path,"thermo_library.py"))

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

def to_dict(self):
    out_dict = dict()
    attrs = [attr for attr in dir(self) if not attr.startswith("_")]
    for attr in attrs:
        val = getattr(self,attr)
    
        if not isinstance(val,list) and not isinstance(val,np.ndarray) and not isinstance(val,dict) and callable(val):
            continue
        
        try:
            json.dumps(val)
            out_dict[attr] = val 
            assert not isinstance(val,np.ndarray)
        except:
            if isinstance(val, ScalarQuantity):
                out_dict[attr] = {
                "class": val.__class__.__name__,
                "value": val.value,
                "units": val.units,
                "uncertainty": val.uncertainty,
                "uncertainty_type": val.uncertainty_type,
            }
            elif isinstance(val, Atoms):
                new_dict = val.todict()
                new_dict["class"] = val.__class__.__name__
                out_dict[attr] = new_dict
            else:
                out_dict[attr] = to_dict
                
    return out_dict