from ase.io import read
from ase.io.trajectory import Trajectory
from ase.visualize import view
import json
import os
import numpy as np
import matplotlib.pyplot as plt
from ase.vibrations import VibrationsData
from ase.thermochemistry import HarmonicThermo, IdealGasThermo
from molecule.kinetics import SurfaceArrhenius
import ase.units
from molecule.molecule import Molecule
from acat.adsorption_sites import SlabAdsorptionSites
from molecule.thermo import Wilhoit
from typing import List, Union, Tuple, Dict
from ase import db
from ase.db import connect

eV_to_Jmol = 9.648328e4

def plot_eharm(path,Eharmtol=3.0,Eharmfiltertol=30.0):
    eharms = []
    guess_dirs = os.listdir(path)
    for guess in guess_dirs:
        try:
            with open(os.path.join(path,guess,"harm.json")) as f:
                out = json.load(f)
                eharms.append(out["harmonic energy"])
        except:
            pass
    eharms = sorted(eharms)
    eharmmin = np.min(eharms)

    plt.plot(eharms)
    plt.xlabel("TS Guess in Order of Increasing Energy")
    plt.ylabel("Harmonic Energy [eV]")
    plt.yscale("log")
    plt.plot(range(len(eharms)),np.ones(len(eharms))*eharmmin*Eharmtol)
    plt.plot(range(len(eharms)),np.ones(len(eharms))*eharmmin*Eharmfiltertol)

def get_opt_dirs(path):
    opt_dirs = []
    guess_dirs = os.listdir(path)
    for guess in guess_dirs:
        p = os.path.join(path,guess,"opt.xyz")
        if os.path.exists(p):
            opt_dirs.append(p)
    return opt_dirs

def get_opt_traj_dirs(path):
    opt_dirs = []
    guess_dirs = os.listdir(path)
    for guess in guess_dirs:
        p = os.path.join(path,guess,"opt.xyz.traj")
        if os.path.exists(p):
            opt_dirs.append(p)
    return opt_dirs

def get_freq_dirs(path):
    freq_dirs = []
    guess_dirs = os.listdir(path)
    for guess in guess_dirs:
        p = os.path.join(path,guess,"vib.0.traj")
        if os.path.exists(p):
            freq_dirs.append(p)
    return freq_dirs

def get_irc_dirs(path):
    freq_dirs = []
    guess_dirs = os.listdir(path)
    for guess in guess_dirs:
        p = os.path.join(path,guess,"irc_forward.traj")
        if os.path.exists(p):
            freq_dirs.append(p)
        p = os.path.join(path,guess,"irc_reverse.traj")
        if os.path.exists(p):
            freq_dirs.append(p)
    return freq_dirs

def get_energies(path,atom_corrections=None):
    Es = dict()
    thermos = dict()
    fs = dict()
    guess_dirs = os.listdir(path)
    slab = read(os.path.join(os.path.split(path)[0],"slab.xyz"))
    Eslab = slab.get_potential_energy()
    with open(os.path.join(path,"info.json"),'r') as f:
        info = json.load(f)
    
    m = Molecule().from_adjacency_list(info["reactants"])
    if atom_corrections:
        AEC = 0.0
        for atom in m.atoms:
            s = str(atom.element)
            if not "X" in s:
                AEC += atom_corrections[s]
    else: #if no corrections don't add a correction
        AEC = 0.0
    
    for guess in guess_dirs:
        p = os.path.join(path,guess,"vib.json_vib.json")
        if os.path.exists(p):
            sp = read(os.path.join(path,guess,"opt.xyz"))
            freqdir = os.path.join(path,guess,"vib.json_vib.json")
            if not os.path.exists(freqdir):
                Es[guess] = None
                thermos[guess] = None
                continue
            vibdata = get_vibdata(os.path.join(path,guess,"opt.xyz"),freqdir,len(slab))
            freqs = vibdata.get_frequencies()
            fs[guess] = freqs
            freqs = [np.complex128(x) for x in freqs]
            ind = np.argmax(np.abs(np.imag(np.array(freqs))))

            ZPE = vibdata.get_zero_point_energy()
            E = sp.get_potential_energy()
            Es[guess] = E+ZPE-Eslab-AEC
            thermos[guess] = (HarmonicThermo(np.array([x for x in np.real(vibdata.get_energies()) if x > 0.0]),potentialenergy=E-Eslab-AEC))

    return Es,thermos,fs

def get_kinetics(path,metal,facet):
    info = json.load(open(os.path.join(path,"info.json"),'r'))
    slab = read(os.path.join(os.path.split(path)[0],"slab.xyz"))
    site_density = get_site_density(slab,metal,facet)
    rnames = info["species_names"]
    pnames = info["reverse_names"]
    if info["forward"]:
        rnum = len(Molecule().from_adjacency_list(info["reactants"]).split())
        pnum = len(Molecule().from_adjacency_list(info["products"]).split())
    else:
        rnum = len(Molecule().from_adjacency_list(info["products"]).split())
        pnum = len(Molecule().from_adjacency_list(info["reactants"]).split())
    rE,pE,rthermos,pthermos = get_reactant_products_energy(path,rnames,pnames)
    Es,thermos,_ = get_energies(path)

    fEs = {k: E-rE if E else None for k,E in Es.items()}
    rEs = {k: E-pE if E else None for k,E in Es.items()}
    fks = dict()
    if rthermos:
        for i,fE in fEs.items():
            if thermos[i]:
                arr = fit_rate_coefficient(rthermos,thermos[i],fE,rnum,s0=site_density)
                fks[i] = arr
            else:
                fks[i] = None
    rks = dict()
    if pthermos:
        for i,rE in rEs.items():
            if thermos[i]:
                arr = fit_rate_coefficient(pthermos,thermos[i],rE,pnum,s0=site_density)
                rks[i] = arr
            else:
                rks[i] = None
    return fEs,rEs,fks,rks,[get_nasa_for_species(th) for th in rthermos],[get_nasa_for_species(th) for th in pthermos]

def fit_rate_coefficient(thermoRs,thermoTS,dE,rnum,s0,Ts=None):
    if Ts is None:
        Ts = np.linspace(298.0,2500.0,50)
    kB = ase.units.kB
    h = 6.582119569e-16 # eV * s
    c0 = 2.4282*10**22*1000.0 #molecules/m^3
    Lunits = None
    ks = []

    for T in Ts:
        GTS = thermoTS.get_helmholtz_energy(T,verbose=False)
        GR = 0.0
        factor = s0
        Lunits = -2

        for therm in thermoRs:
            if isinstance(therm,HarmonicThermo):
                GR += therm.get_helmholtz_energy(T,verbose=False)
                factor /= s0
                Lunits += 2
            elif isinstance(therm,IdealGasThermo):
                GR += therm.get_gibbs_energy(T,1.0e5,verbose=False)
                factor /= c0
                Lunits += 3
        for i in range(rnum-len(thermoRs)): # account for sites
            factor /= s0
            Lunits += 2

        k = kB*T/h * np.exp(-(GTS-GR)/(kB*T)) * factor
        ks.append(k)

    molecules = rnum-1
    if Lunits == 0 and molecules==0:
        arr = SurfaceArrhenius().fit_to_data(Ts,np.array(ks),"s^-1")
    elif molecules == 1:
        arr = SurfaceArrhenius().fit_to_data(Ts,np.array(ks),"m^{Lunits}/(molecule*s)".format(Lunits=Lunits))
    else:
        arr = SurfaceArrhenius().fit_to_data(Ts,np.array(ks),"m^{Lunits}/(molecules^{mols}*s)".format(Lunits=Lunits,mols=molecules))

    return arr

def get_adsorbate_energies(ad_path,atom_corrections=None,include_zpe=True):
    """
    get ZPE corrected adsorbate energies
    """
    dirs = os.listdir(ad_path)
    slab = read(os.path.join(os.path.split(os.path.split(ad_path)[0])[0],"slab.xyz"))
    Eslab = slab.get_potential_energy()
    with open(os.path.join(ad_path,"info.json")) as f:
        info = json.load(f)

    m = Molecule().from_adjacency_list(info["adjlist"])
    if atom_corrections:
        AEC = 0.0
        for atom in m.atoms:
            s = str(atom.element)
            if not "X" in s:
                AEC += atom_corrections[s]
    else: #if no corrections don't add a correction
        AEC = 0.0
    
    gasphase = len(info["gratom_to_molecule_surface_atom_map"]) == 0
    spin = (Molecule().from_adjacency_list(info["adjlist"]).multiplicity - 1.0)/2.0
    Es = dict()
    thermos = dict()
    fs = dict()
    for d in dirs:
        if d == "info.json":
            continue
        optdir = os.path.join(ad_path,d,d+".xyz")
        freqdir = os.path.join(ad_path,d,"vib.json_vib.json")
        if not (os.path.exists(optdir) and os.path.exists(freqdir)):
            continue
        sp = read(optdir)
        E = sp.get_potential_energy()
        if gasphase:
            vibdata = get_vibdata(optdir,freqdir,0)
        else:
            vibdata = get_vibdata(optdir,freqdir,len(slab))

        freqs = vibdata.get_frequencies().tolist()
        fs[d] = freqs

        ZPE = vibdata.get_zero_point_energy()

        if not gasphase:

            if include_zpe:
                Es[d] = E + ZPE - Eslab - AEC
            else:
                Es[d] = E - Eslab - AEC
            thermos[d] = HarmonicThermo(np.array([x for x in np.real(vibdata.get_energies()) if x > 0.0]),potentialenergy=E-Eslab-AEC)
        else:
            if len(sp) == 1:
                geometry = 'monatomic'
            elif any([ I < 1.0e-3 for I in sp.get_moments_of_inertia()]):
                geometry = "linear"
            else:
                geometry = "nonlinear"

            if include_zpe:
                Es[d] = E + ZPE - AEC
            else:
                Es[d] = E - AEC

            thermos[d] = IdealGasThermo(np.real(vibdata.get_energies()),geometry,
                                        potentialenergy=E - AEC,atoms=sp,symmetrynumber=1,
                                        natoms=len(sp),spin=spin)

    return Es,thermos,fs



def get_reactant_products_energy(ts_path,reactants,products):
    path = os.path.split(ts_path)[0]
    ads_path = os.path.join(path,"Adsorbates")
    rthermos = []
    pthermos = []
    rE = 0.0

    for r in reactants:
        dr,rthermo,_ = get_adsorbate_energies(os.path.join(ads_path,r))
        Emin = np.inf
        ind = -1
        for key,val in dr.items():
            if val and val < Emin:
                ind = key
                Emin = val
        if ind == -1:
            rthermos = []
            break
        rE += dr[ind]
        
        rthermos.append(rthermo[ind])

    pE = 0.0
    for p in products:
        dp,pthermo,_ = get_adsorbate_energies(os.path.join(ads_path,p))
        Emin = np.inf
        ind = -1
        for key,val in dp.items():
            if val and val < Emin:
                ind = key
                Emin = val
        if ind == -1:
            pthermos = []
            break
        pE += dp[ind]
        pthermos.append(pthermo[ind])
    return rE,pE,rthermos,pthermos

def get_vibdata(dopt,dvib,nslab):
    xyz = read(dopt)
    inds = range(nslab,len(xyz))
    with open(dvib,'r') as f:
        out = json.load(f)
    sh = (len(inds),3,len(inds),3)
    H = np.asarray(out["hessian"])
    if H.shape != sh:
        H = H.reshape(sh)
    vib = VibrationsData(xyz,H,inds)
    return vib

def get_animated_mode(dopt,dvib,nslab,n=0):

    vib = get_vibdata(dopt,dvib,nslab)
    atoms = vib.iter_animated_mode(n)

    tr = Trajectory("temp_vib_vib.traj",mode='a')
    for a in atoms:
        tr.write(a)
    tr2 = Trajectory("temp_vib_vib.traj")
    view(tr2)
    os.remove("temp_vib_vib.traj")

def get_site_density(slab,metal,facet):
    cas = SlabAdsorptionSites(slab, facet,allow_6fold=False,composition_effect=False,
                        label_sites=True,
                        surrogate_metal=metal)
    S = len(cas.get_sites())
    cell = slab.cell
    n = np.cross(cell[0],cell[1])
    A = np.linalg.norm(n)
    return S/A * 10**20 #molecules/m^2

def get_cp(th,T,dT=0.01):
    if isinstance(th,HarmonicThermo):
        # The following is the approximate derivative using the finite difference method to obtain cp
        return ((th.get_helmholtz_energy(T+dT,verbose=False) +
                 (T+dT)*th.get_entropy(T+dT,verbose=False)) -
                (th.get_helmholtz_energy(T-dT,verbose=False) +
                 (T-dT)*th.get_entropy(T-dT,verbose=False)))/(2*dT)
    elif isinstance(th,IdealGasThermo):
        return ((th.get_enthalpy(T+dT,verbose=False) +
                 (T+dT)*th.get_entropy(T+dT,pressure=1.0e5,verbose=False)) -
                (th.get_enthalpy(T-dT,verbose=False) +
                 (T-dT)*th.get_entropy(T-dT,pressure=1.0e5,verbose=False)))/(2*dT)

def get_nasa_for_species(th,dT=0.01):
    
    if isinstance(th,HarmonicThermo):
        G298 = th.get_helmholtz_energy(298.0,verbose=False) * eV_to_Jmol
        S298 = th.get_entropy(298.0,verbose=False) * eV_to_Jmol
    elif isinstance(th,IdealGasThermo):
        G298 = th.get_gibbs_energy(298.0,pressure=1.0e5,verbose=False) * eV_to_Jmol
        S298 = th.get_entropy(298.0,pressure=1.0e5,verbose=False) * eV_to_Jmol
    else:
        raise ValueError
    H298 = G298 + 298.0*S298
    Cp0 = get_cp(th,1.0,dT=0.1) * eV_to_Jmol
    CpInf = get_cp(th,1e7,dT=dT) * eV_to_Jmol
    Cps = []
    Ts = np.arange(10.0, 3001.0, 10.0, np.float64)
    for T in Ts:
        Cps.append(get_cp(th,T,dT=dT)*eV_to_Jmol)

    wh = Wilhoit().fit_to_data(Tdata=Ts,Cpdata=np.array(Cps),Cp0=Cp0,CpInf=CpInf,H298=H298,S298=S298)

    nasa = wh.to_nasa(Tmin=10.0, Tmax=3000.0, Tint=500.0)
    return nasa

def get_gibbs_energy_reaction(rnasas,pnasas,T,dT=0.01): 
    dG = 0.0
    for th in rnasas:
        dG -= th.get_free_energy(T)
    for th in pnasas:
        dG += th.get_free_energy(T)
    return dG

def get_entropy_reaction(rnasas,pnasas,T,dT=0.01): 
    dG = 0.0
    for th in rnasas:
        dG -= th.get_entropy(T)
    for th in pnasas:
        dG += th.get_entropy(T)
    return dG

def get_enthalpy_reaction(rnasas,pnasas,T,dT=0.01): 
    dG = 0.0
    for th in rnasas:
        dG -= th.get_enthalpy(T)
    for th in pnasas:
        dG += th.get_enthalpy(T)
    return dG

def write_min_en_species_db(
        ad_path: str,
        slab_path: str,
) -> db:
    '''
    Generates an ase database of the minimum energy configurations for all
    Pynta optimized adsorbates and gas-phase species.

    Parameters:
        ad_path: Absolute path to Pynta generated Adsorbates directory.
        slab_path: Absolute path to Pynta optimized slab. ('your_path/slab.xyz')
    Returns:
        db: ase database
    '''
    os.chdir(ad_path)
    ad_dirs = [dir for dir in filter(os.path.isdir, os.listdir(ad_path))]
    db = connect('min_E_ads.db')
    for ad in ad_dirs:
        path_to_ad = os.path.join(ad_path, ad)
        slab = read(slab_path)

        try:
            with open(os.path.join(path_to_ad, "info.json")) as f:
                info = json.load(f)
        except:
            print(f"No information for {path_to_ad}")
            continue

        adj_list = info['adjlist']
        gasphase = len(info["gratom_to_molecule_surface_atom_map"]) == 0
        spin = (Molecule().from_adjacency_list(info["adjlist"]).multiplicity - 1.0) / 2.0
        name = info["name"]
        os.chdir(path_to_ad)
        ad_configs = [dir for dir in filter(os.path.isdir, os.listdir(path_to_ad))]

        # Setting up empty dictionaries to collect information for each config
        ad_dict = {}
        DFT_en_dict = {}
        freq_dict = {}
        atoms_dict = {}
        for config in ad_configs:
            optdir = os.path.join(path_to_ad, config, config + ".xyz")
            freqdir = os.path.join(path_to_ad, config, "vib.json_vib.json")
            atoms = read(optdir)
            atoms_dict[config] = atoms
            energy = atoms.get_potential_energy()
            if gasphase:
                vibdata = get_vibdata(optdir, freqdir, 0)
            else:
                vibdata = get_vibdata(optdir, freqdir, len(slab))
            freqs = vibdata.get_frequencies().tolist()

            # Handling of imaginary frequencies.
            for index, freq in enumerate(freqs):
                if freq.imag != 0:
                    # print(f" for {ad_path}/{d} \n we have imaginary freq \n setting to 12cm^-1")
                    freqs[index] = 12
                else:
                    freqs[index] = freq.real
            zpe = vibdata.get_zero_point_energy()
            freq_dict[config] = freqs
            DFT_en_dict[config] = energy
            ad_dict[config] = energy + zpe
        min_E_config = min(ad_dict, key=ad_dict.get)
        atoms = atoms_dict[min_E_config]
        data = {}
        data['frequencies'] = freq_dict[min_E_config]
        db.write(atoms, name=name, adj_list=adj_list, spin=spin, gasphase=gasphase, data=data)


class thermo:
    """Base class for calculating thermodynamic properties from Pynta optimized
    adsorbates."""
    def __init__(self,
                 atoms,
                 DFT_energy,
                 ZPE_energy,
                 name,
                 adj_list,
                 spin,
                 gasphase,
                 data,
                 dft_en_slab,
                 cat_element,
                 sites_per_cell,
                 slab_len,
                dft_en_C_ref = 0,
                dft_en_O_ref = 0,
                dft_en_H_ref = 0,
                dft_en_N_ref = 0):

        # start by defining some physical constants
        self.R = 8.3144621E-3  # ideal Gas constant in kJ/mol-K
        self.kB = 1.38065e-23  # Boltzmann constant in J/K
        self.h = 6.62607e-34  # Planck constant in J*s
        self.c = 2.99792458e8  # speed of light in m/s
        self.amu = 1.6605e-27  # atomic mass unit in kg
        self.Avogadro = 6.0221E23  # mole^-1
        self.GHz_to_Hz = 1.0E9  # convert rotational constants from GHz to Hz
        self.invcm_to_invm = 1.0E2  # convert cm^-1 to m^-1, for frequencies
        self.P_ref = 1.0E5  # reference pressure, 1 bar = 1E5 Pascal
        self.hartree_to_kcalpermole = 627.5095  # convert hartree/molecule to kcal/mol
        self.hartree_to_kJpermole = 2627.25677  # convert hartree/molecule to kJ/mol
        self.eV_to_kJpermole = 96.485  # convert eV/molecule to kJ/mol


        self.T_switch = 1000.0  # K, switching temperature in NASA polynomial. Default. can overwrite.

        # System specific stuff
        #self.name.count(
        '''THE NUMBER OF SITES IS USED FOR TRANSLATION ENTROPY, NOT SURE RIGHT NOW HOW TO HANDLE THIS WITH BIDENTATE'''
        area = np.norm(np.cross(atoms.cell[0], atoms.cell[1]))
        self.unit_cell_area = (area * 10**-20) / sites_per_cell  # m2 - using surface area per binding site

        self.cutoff_frequency = 1.0  # cm^-1
        self.twoD_gas = False
        self.dHfatct = {'CH4': -66.540, 'H2O': -238.903, 'H2': 0,
                        'NH3': -38.563}  # heats of formation (0K) in kJ/mol from the ATcT database for the reference species, version 1.202
        self.dHrxnatct = {'H2-2H': 432.0680600, 'O2-2O': 493.6871,
                          'N2-2N': 941.165}  # Heats of the dissociation reactions in the gas phase from the ATcT database, version 1.202
        """
        These are the ZPE-correct DFT energies of your gas phase references
        """
        self.Eref = {'CH4': dft_en_C_ref, 'H2O': dft_en_O_ref, 'H2': dft_en_H_ref, 'NH3': dft_en_N_ref}  # eV
        self.Eslab = dft_en_slab  # DFT energy of the relaxed bare slab in eV

        #NEED TO DOUBLE CHECK THAT THIS THERMOCHEM MODEL IS APPROPRIATE FOR BIDENTATE.
        if not gasphase:
            self.sites = name.count('[Pt]')
        else:
            self.sites = 0
        formula = atoms.get_chemical_formula(mode='all')

        # THE FOLLOWING DOES NOT WORK FOR BIMETALLIC OR OXIDE SURFACES
        self.composition = {'H': formula.count('H'), 'C': formula.count('C'),
                            'N': formula.count('N'), 'O': formula.count('O'),
                            cat_element: self.sites}
        self.N_ad_atoms = len([atom for atom in atoms if atom.symbol != cat_element])
        print(f"number of adsorbate atoms = {self.N_ad_atoms}")
        self.adsorbate_mass = sum(atoms.get_masses()[slab_len:])
        self.adsorbate_mass_units = 'amu'

        # Currently we don't calculate gas phase species for all intermediates
        self.energy_gas = 0
        self.DFT_energy = DFT_energy
        self.ZPE_energy = ZPE_energy
        self.DFT_energy_gas_units = 'eV'
        self.DFT_energy_units = 'eV'
        self.ZPE_energy_units = 'eV'

        new_list = adj_list.replace("['", "").replace("']", "").replace(" '", "").strip().split("',")
        self.adjacency_list = adj_list
        self.formatted_adj_list = '\n'.join(new_list)

        self.spin = spin
        self.data = data
        self.name = name

        self.frequencies_units = 'cm-1'

        # Determining whether translational entropy will be included
        if self.data['frequencies'][1] < self.cutoff_frequency:
            self.twoD_gas = True

        



