from molecule.molecule import Molecule,Atom,Bond
from molecule.exceptions import AtomTypeError
from ase.io import read, write
from ase.geometry import get_distances
from ase.visualize import view
from acat.adsorption_sites import SlabAdsorptionSites
from pynta.utils import get_unique_sym, get_occupied_sites, sites_match
from pynta.mol import *
from pynta.geometricanalysis import *
from pynta.tasks import *
from fireworks import LaunchPad, Workflow
from fireworks.queue.queue_launcher import rapidfire as rapidfirequeue
from fireworks.features.multi_launcher import launch_multiprocess
from fireworks.utilities.fw_serializers import load_object_from_file
from fireworks.core.rocket_launcher import rapidfire
from fireworks.core.fworker import FWorker
from copy import deepcopy
import numpy as np
import json
import shutil
import os
import itertools 
import logging

def get_unstable_pairs(pairsdir,adsorbate_dir,sites,site_adjacency,nslab,max_dist=3.0,show=False):
    out_pairs = []
    coadname_dict = {"O=[Pt]": 1, "N#[Pt]": 1, "O[Pt]": 2, "[Pt]": 1}
    allowed_structure_site_structure_map = generate_allowed_structure_site_structure_map(adsorbate_dir,sites,site_adjacency,nslab,max_dist=max_dist,cut_coads_off_num=None)
    if show:
        config_show = []
    for pair in os.listdir(pairsdir):
        if not "_" in pair or pair[0] == ".":
            continue
        coadname = pair.split("_")[1]
        for num in os.listdir(os.path.join(pairsdir,pair)):
            p = os.path.join(pairsdir,pair,num) 
            if num.isdigit() and os.path.isdir(p):
                init = read(os.path.join(p,"init.xyz"))
                with open(os.path.join(p,"info.json"),'r') as f:
                    initinfo = json.load(f)
                try:
                    m = Molecule().from_adjacency_list(initinfo["adjlist"])
                    #m.update(sort_atoms=False)
                except AtomTypeError:
                    m = Molecule().from_adjacency_list(initinfo["adjlist"],raise_atomtype_exception=False,
                            raise_charge_exception=False,check_consistency=False)
                    for a in m.atoms:
                        if a.charge != 0:
                            a.charge = 0
                            v = 0
                            for a2,edge in a.edges.items():
                                v += edge.order*2
                            if not a.is_hydrogen():
                                a.lone_pairs = int(np.round((8-v)/2))
#                     try:
#                         m.update(sort_atoms=False)
#                     except Exception as e:
#                         print(m.to_adjacency_list())
#                         raise e
                
                g = extract_pair_graph(init,sites,site_adjacency,nslab,max_dist=max_dist,
                                       cut_coads_off_num=coadname_dict[coadname],allowed_structure_site_structure_map=allowed_structure_site_structure_map)
                    
                
                #g.update(sort_atoms=False)
                outpath = os.path.join(p,"out.xyz")
                if not os.path.exists(outpath):
                    if show:
                        pass
                        
                    out_pairs.append(g.to_group())
                else:
                    final = read(outpath)
                    try:
                        gout = extract_pair_graph(final,sites,site_adjacency,nslab,max_dist=max_dist,allowed_structure_site_structure_map=allowed_structure_site_structure_map)
                        if len(gout.atoms) != len(g.atoms):
                            out_pairs.append(g.to_group())
                            if show:
                                config_show.append(init)
                                config_show.append(final)
                            continue
                        #gout.update(sort_atoms=False)
                    except FindingPathError:
                        continue
                    
                    if not gout.is_isomorphic(g,save_order=True,strict=False):
                        out_pairs.append(g.to_group())
                        if show:
                            config_show.append(init)
                            config_show.append(final)
    if show:
        view(config_show)
    return out_pairs

def copy_stable_pairs(pairsdir,sites,site_adjacency,nslab,max_dist=3.0):
    out_pairs = []
    count = 0
    coadname_dict = {"O=[Pt]": 1, "N#[Pt]": 1, "O[Pt]": 2, "[Pt]": 1}
    for pair in os.listdir(pairsdir):
        coadname = pair.split("_")[1]
        for num in os.listdir(os.path.join(pairsdir,pair)):
            p = os.path.join(pairsdir,pair,num) 
            if num.isdigit() and os.path.isdir(p):
                init = read(os.path.join(p,"init.xyz"))
                with open(os.path.join(p,"info.json"),'r') as f:
                    initinfo = json.load(f)
                try:
                    m = Molecule().from_adjacency_list(initinfo["adjlist"])
                    #m.update(sort_atoms=False)
                except AtomTypeError:
                    m = Molecule().from_adjacency_list(initinfo["adjlist"],raise_atomtype_exception=False,
                            raise_charge_exception=False,check_consistency=False)
                    for a in m.atoms:
                        if a.charge != 0:
                            a.charge = 0
                            v = 0
                            for a2,edge in a.edges.items():
                                v += edge.order*2
                            if not a.is_hydrogen():
                                a.lone_pairs = int(np.round((8-v)/2))
#                     try:
#                         m.update(sort_atoms=False)
#                     except Exception as e:
#                         print(m.to_adjacency_list())
#                         raise e
                
                g = extract_pair_graph(init,sites,site_adjacency,nslab,max_dist=max_dist,
                                       cut_coads_off_num=coadname_dict[coadname])
                    
                
                #g.update(sort_atoms=False)
                outpath = os.path.join(p,"out.xyz")
                if not os.path.exists(outpath):
                    out_pairs.append(g.to_group())
                else:
                    final = read(outpath)
                    try:
                        gout = extract_pair_graph(final,sites,site_adjacency,nslab,max_dist=max_dist)
                        if len(gout.atoms) != len(g.atoms):
                            out_pairs.append(g.to_group())
                            continue
                        #gout.update(sort_atoms=False)
                    except FindingPathError:
                        continue
                    
                    if not gout.is_isomorphic(g,save_order=True,strict=False):
                        out_pairs.append(g.to_group())
                    else: #matches
                        os.makedirs(os.path.join(os.path.split(pairsdir)[0],"sampling",str(count)))
                        shutil.copyfile(os.path.join(p,"out.xyz"),
                                        os.path.join(os.path.split(pairsdir)[0],"sampling",str(count),"out.xyz")
                                       )
                        shutil.copyfile(os.path.join(p,"info.json"),
                                        os.path.join(os.path.split(pairsdir)[0],"sampling",str(count),"info.json")
                                       )
                        shutil.copyfile(os.path.join(p,"init.xyz"),
                                        os.path.join(os.path.split(pairsdir)[0],"sampling",str(count),"init.xyz")
                                       )
                        count += 1


    return None
