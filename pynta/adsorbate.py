from pynta.mol import *
from pysidt.sidt import read_nodes
from pynta.sidt import SurfaceBondLengthSIDT
from pynta.geometricanalysis import generate_site_molecule
import logging 
import os 
import pynta

def generate_adsorbate_guesses(mol,ads,slab,mol_to_atoms_map,metal,
                               single_site_bond_params_lists,single_sites_lists,double_site_bond_params_lists,double_sites_lists,
                               Eharmtol,Eharmfiltertol,Ntsmin,slab_sites,site_adjacency):
    mol_surf_inds = [mol.atoms.index(a) for a in mol.get_adatoms()]
    atom_surf_inds = [mol_to_atoms_map[i] for i in mol_surf_inds]
    nslab = len(slab)
    slabmol,neighbor_sites,ninds = generate_site_molecule(slab, slab_sites, slab_sites, site_adjacency, max_dist=None)
    admols = []
    surface_site_indices = []
    if len(atom_surf_inds) == 1:
        site_bond_params_lists = deepcopy(single_site_bond_params_lists)
        sites_lists = single_sites_lists
        for i,site_bond_params_list in enumerate(site_bond_params_lists):
            site_bond_params_list[0]["ind"] = atom_surf_inds[0]+len(slab)
            
            admol = slabmol.copy(deep=True)
            m = mol.copy(deep=True)
            m_site = m.get_surface_sites()[0]
            bd = list(m_site.bonds.values())[0]
            order = bd.order 
            m_adatom = list(m_site.bonds.keys())[0]
            admol_site_index = [j for j,s in enumerate(neighbor_sites) if sites_match(s,sites_lists[i][0],slab)][0]
            admol_site = admol.atoms[admol_site_index]
            m.remove_atom(m_site)
            admol = admol.merge(m)
            admol.add_bond(Bond(m_adatom,admol_site,order=order))
            admol.update_multiplicity()
            admol.update_atomtypes()
            admol.update_connectivity_values()
            admol.identify_ring_membership()
            admols.append(admol)
            surface_site_indices.append([admol_site_index])
            
            #add up pulling potential
            for ind in range(len(ads)):
                if ind in atom_surf_inds:
                    continue
                pos = deepcopy(site_bond_params_list[0]["site_pos"])
                pos[2] += 8.5
                site_bond_params_list.append({"site_pos": pos,"ind": ind+len(slab), "k": 0.1, "deq": 0.0})

    elif len(atom_surf_inds) == 2:
        site_bond_params_lists = deepcopy(double_site_bond_params_lists)
        sites_lists = double_sites_lists
        for i,site_bond_params_list in enumerate(site_bond_params_lists):
            site_bond_params_list[0]["ind"] = atom_surf_inds[0]+len(slab)
            site_bond_params_list[1]["ind"] = atom_surf_inds[1]+len(slab)
            
            admol = slabmol.copy()
            m = mol.copy(deep=True)
            m_sites = m.get_surface_sites()
            bds = list(m_site.bonds.values())
            orders = [bd.order for bd in bds]
            m_adatoms = [list(m_site.bonds.keys())[0] for msite in m_sites]
            admol_site_indices = [[j for j,s in enumerate(neighbor_sites) if sites_match(s,sites_lists[i][0],slab)][0],[j for j,s in enumerate(neighbor_sites) if sites_match(s,sites_lists[i][1],slab)][0]]
            admol_sites = [admol.atoms[admol_site_index] for admol_site_index in admol_site_indices]
            for m_site in m_sites:
                m.remove_atom(m_site)
            admol = admol.merge(m)
            for i,m_adatom in m_adatoms:
                admol.add_bond(Bond(m_adatom,admol_sites[i],order=orders[i]))
            admol.update_multiplicity()
            admol.update_atomtypes()
            admol.update_connectivity_values()
            admol.identify_ring_membership()
            admols.append(admol)
            surface_site_indices.append(admol_site_indices)
    else:
        raise ValueError("Only monodentate and bidentate guesses currently allowed. The infrastructure can support tridenate and higher, but the filtering process may be very expensive above bidentate.")


    mol_fixed_bond_pairs = [[mol.atoms.index(bd.atom1),mol.atoms.index(bd.atom2)] for bd in mol.get_all_edges() if (not bd.atom1.is_surface_site()) and (not bd.atom2.is_surface_site())]
    atom_fixed_bond_pairs = [[mol_to_atoms_map[pair[0]]+len(slab),mol_to_atoms_map[pair[1]]+len(slab)]for pair in mol_fixed_bond_pairs]
    constraint_list = [{"type": "FixBondLength", "a1": pair[0], "a2": pair[1]} for pair in atom_fixed_bond_pairs]+["freeze slab"]
    nodes_file = os.path.join(os.path.split(pynta.models.__file__)[0],"adsorbate_bond_length_tree_nodes.json")
    nodes = read_nodes(nodes_file)
    sidt_model = SurfaceBondLengthSIDT(nodes=nodes)
    
    geos = []
    for i,sites_list in enumerate(sites_lists):
        geo,h1,h2 = place_adsorbate(ads,slab,atom_surf_inds,sites_list,admols[i],surface_site_indices[i],metal,sidt_model)
        if h1:
            site_bond_params_lists[i][0]["site_pos"][2] += h1
        if h2:
            site_bond_params_lists[i][1]["site_pos"][2] += h2
        geos.append(geo)

    print("initial geometries")
    print(len(geos))
    geos_out = []
    Eharms = []
    site_bond_params_lists_out = []
    for i,geo in enumerate(geos):
        #freeze bonds for messier first opt
        geo_out,Eharm,Fharm = run_harmonically_forced_xtb(geo,[],site_bond_params_lists[i],len(slab),
                                molecule_to_atom_maps=mol_to_atoms_map,ase_to_mol_num=None,
                                method="GFN1-xTB",constraints=constraint_list)
        if geo_out:
            geo_out.calc = None
            geos_out.append(geo_out)
            Eharms.append(Eharm)
            site_bond_params_lists_out.append(site_bond_params_lists[i])

    print("optimized geometries")
    print(len(geos_out))
    inds = get_unique_sym_struct_indices(geos_out)

    print("after symmetry")
    print(len(inds))

    geos_out = [geos_out[ind] for ind in inds]
    Eharms = [Eharms[ind] for ind in inds]

    if len(atom_surf_inds) == 1: #should be small, don't bother filtering
        xyzsout = geos_out
        site_bond_params_lists_final = [site_bond_params_lists_out[ind] for ind in inds]
        return xyzsout
    else:
        Einds = np.argsort(np.array(Eharms))
        Emin = np.min(np.array(Eharms))
        xyzsout = []
        site_bond_params_lists_final = []
        for Eind in Einds:
            if Eharms[Eind]/Emin < Eharmtol: #include all TSs with energy close to Emin
                xyzsout.append(geos_out[Eind])
                site_bond_params_lists_final.append(site_bond_params_lists_out[Eind])
            elif Eharms[Eind]/Emin > Eharmfiltertol: #if the energy is much larger than Emin skip it
                continue
            elif len(xyzsout) < Ntsmin: #if the energy isn't similar, but isn't much larger include the smallest until Ntsmin is reached
                xyzsout.append(geos_out[Eind])
                site_bond_params_lists_final.append(site_bond_params_lists_out[Eind])

    return xyzsout

def place_adsorbate(ads,slab,atom_surf_inds,sites,admol,admol_site_indices,metal,sidt_model):
    if len(atom_surf_inds) == 1:
        geo = slab.copy()
        h = estimate_surface_bond_length(admol,admol_site_indices[0],sidt_model,metal)
        add_adsorbate_to_site(geo, ads, atom_surf_inds[0], sites[0], height=h)
        return geo,h,None
    elif len(atom_surf_inds) == 2:
        geo = slab.copy()
        h1 = estimate_surface_bond_length(admol,admol_site_indices[0],sidt_model,metal)
        h2 = estimate_surface_bond_length(admol,admol_site_indices[1],sidt_model,metal)
        ori = get_mic(sites[0]['position'], sites[1]['position'], geo.cell)
        add_adsorbate_to_site(geo, deepcopy(ads), atom_surf_inds[0], sites[0], height=h1, orientation=ori)
        if np.isnan(geo.positions).any(): #if nans just ignore orientation and let it optimize
            geo = slab.copy()
            add_adsorbate_to_site(geo, deepcopy(ads), atom_surf_inds[0], sites[0], height=h1, orientation=None)
        return geo,h1,h2
    else:
        raise ValueError

def estimate_surface_bond_length(admol,admol_site_index,sidt_model,metal):
    m = admol.copy(deep=True)
    site_atom = m.atoms[admol_site_index]
    adatom = [a for a in site_atom.edges.keys() if not a.is_surface_site()][0]
    site_atom.label = "*"
    adatom.label = "*"
    L,tr = sidt_model.evaluate(m,metal,trace=True)
    if tr == "Root":
        logging.warning("Estimation used SIDT Root group")
    return L 