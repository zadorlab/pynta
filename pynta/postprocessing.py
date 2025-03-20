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
from pynta.geometricanalysis import get_adsorbate_energies, get_vibdata

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
        return ((th.get_helmholtz_energy(T+dT,verbose=False) + (T+dT)*th.get_entropy(T+dT,verbose=False)) - (th.get_helmholtz_energy(T-dT,verbose=False) + (T-dT)*th.get_entropy(T-dT,verbose=False)))/(2*dT)
    elif isinstance(th,IdealGasThermo):
        return ((th.get_enthalpy(T+dT,verbose=False) + (T+dT)*th.get_entropy(T+dT,pressure=1.0e5,verbose=False)) - (th.get_enthalpy(T-dT,verbose=False) + (T-dT)*th.get_entropy(T-dT,pressure=1.0e5,verbose=False)))/(2*dT)

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
    # Removing old db
    if os.path.exists(os.path.join(os.getcwd(), 'min_E_ads.db')):
        os.system('rm min_E_ads.db')
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
        ad_dict = dict()
        DFT_en_dict = dict()
        freq_dict = dict()
        atoms_dict = dict()
        for config in ad_configs:
            optdir = os.path.join(path_to_ad, config, config + ".xyz")
            freqdir = os.path.join(path_to_ad, config, "vib.json_vib.json")
            if not os.path.exists(freqdir):
                continue

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
                    print(f" for {ad_path}/{ad}/{config} \n we have imaginary freq \n setting to 12cm^-1")
                    freqs[index] = 12
                    if freq.imag > 150:
                        print(f"Imag. freq = {freq.imag}; consider running finer optimization")
                else:
                    freqs[index] = freq.real
            zpe = vibdata.get_zero_point_energy()
            freq_dict[config] = freqs
            DFT_en_dict[config] = energy
            ad_dict[config] = energy + zpe
        if len(ad_dict) == 0:
            continue
        min_E_config = min(ad_dict, key=ad_dict.get)
        atoms = atoms_dict[min_E_config]
        data = {}
        data['frequencies'] = freq_dict[min_E_config]
        db.write(atoms, name=name, zpe=zpe, adj_list=adj_list, spin=spin,
                 gasphase=gasphase, data=data, num_sites=len(info["gratom_to_molecule_surface_atom_map"]))
        return

class Thermo:
    """Base class for calculating thermodynamic properties from Pynta optimized
    adsorbates.

    parameters:
        atoms: ase atoms object
        DFT_energy: potential_energy eV with no ZPE correction
        ZPE_energy: ZPE energy eV
        num_sites: 1 for monodentate, 2 for bidentate, etc
        name: adsorbate name generated from Pynta
        adj_list: adjency list of adsorbate
        spin: multiplicity of adsorbate
        freqs: vibrational frequncies cm^-1
        dft_en_slab: energy of relaxed slab eV
        cat_element: Atomic symbol for catalyst
        slab_len: number of atoms in slab
        author: Name to be used within RMG thermochemistry library
        institution: Name of research instituion to be used in RMG library
        index: used for labeling species in RMG library
        sites_per_cell: We set this to 1 by default, assuming that the 2D lattice gas
            can freely diffuse over the entire area of our unit cell.
        dft_en_C_ref: ZPE corrected gas phase reference for carbon
        dft_en_O_ref: ZPE corrected gas phase reference for oxygen
        dft_en_H_ref: ZPE corrected gas phase reference for hydrogen
        dft_en_N_ref: ZPE corrected gas phase reference for nitrogen
        cutoff_freq: Cutoff used to determine whether to use 2D lattice gas

    Much of this code was copied/adapted from the Goldsmith group at Brown
    University. If you use this code, please cite the following work by Blondal et al.:
    https://doi.org/10.1021/acscatal.2c03378

    """
    def __init__(self,
                 atoms,
                 DFT_energy,
                 ZPE_energy,
                 num_sites,
                 name,
                 adj_list,
                 spin,
                 freqs,
                 dft_en_slab,
                 cat_element,
                 slab_len,
                 author,
                 institution,
                 index,
                 sites_per_cell=1,
                 dft_en_C_ref=0,
                 dft_en_O_ref=0,
                 dft_en_H_ref=0,
                 dft_en_N_ref=0,
                 cutoff_freq=100):

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
        area = np.linalg.norm(np.cross(atoms.cell[0], atoms.cell[1]))
        self.unit_cell_area_per_site = (area * 10 ** -20) / sites_per_cell  # m2 - using surface area per binding site
        self.cutoff_frequency = cutoff_freq  # cm^-1
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

        self.sites = num_sites
        if self.sites == 0:
            self.gasphase = True
        ad_atoms = atoms[slab_len:]
        formula = ad_atoms.get_chemical_formula(mode='all')

        # THE FOLLOWING DOES NOT WORK FOR BIMETALLIC OR OXIDE SURFACES
        # THERE WILL BE MULTIPLE CATALYST ELEMENTS
        self.composition = {'H': formula.count('H'), 'C': formula.count('C'),
                            'N': formula.count('N'), 'O': formula.count('O'),
                            cat_element: self.sites}

        self.N_ad_atoms = len(ad_atoms)
        print(f"number of adsorbate atoms = {self.N_ad_atoms}")
        self.adsorbate_mass = sum(ad_atoms.get_masses())
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
        self.freqs = freqs
        self.name = name

        self.frequencies_units = 'cm-1'

        if self.freqs[1] < self.cutoff_frequency:
            self.twoD_gas = True

        temperature = [298.15]  # NOTE 298.15 must be first for the NASA polynomial routine to work!
        T_low = 300.0
        T_high = 2000.0
        dT = 10.0  # temperature increment
        self.temperature = np.append(temperature, np.arange(T_low, T_high + dT, dT))

        self.cat_element = cat_element

        self.author = author
        self.institution = institution
        self.index = index

    def compute_heat_of_formation(self):

        # Leaving this here in case we want to start calculating gas-phase intermediates
        '''
        # Atomic hydrogen
        if self.name == '[Pt]':
            self.energy_gas = (self.DFT_energy_gas + self.ZPE_energy_gas + self.dHrxnatct[
                'H2-2H'] / self.eV_to_kJpermole) / 2
        # Atomic oxygen
        elif self.name == 'O=[Pt]':
            self.energy_gas = (self.DFT_energy_gas + self.ZPE_energy_gas + self.dHrxnatct[
                'O2-2O'] / self.eV_to_kJpermole) / 2
        #Atomic nitrogen
        elif self.name == 'N#[Pt]':
            self.energy_gas = (self.DFT_energy_gas + self.ZPE_energy_gas + self.dHrxnatct[
                'N2-2N'] / self.eV_to_kJpermole) / 2
        else:
            self.energy_gas = self.DFT_energy_gas + self.ZPE_energy_gas'''
        self.energy = self.DFT_energy + self.ZPE_energy

        self.dHrxndftgas = (self.energy_gas - self.composition['C'] * self.Eref['CH4']
                            - self.composition['O'] * self.Eref['H2O']
                            - self.composition['N'] * self.Eref['NH3']
                            - (self.composition['H'] / 2 - 2 * self.composition['C'] - self.composition[
                    'O'] - 3 / 2 * self.composition['N']) * self.Eref['H2'])
        self.dHfgas = (self.composition['C'] * self.dHfatct['CH4']
                       + self.composition['O'] * self.dHfatct['H2O']
                       + self.composition['N'] * self.dHfatct['NH3']
                       + (self.composition['H'] / 2 - 2 * self.composition['C'] - self.composition[
                    'O'] - 3 / 2 * self.composition['N']) * self.dHfatct['H2']
                       + self.dHrxndftgas * self.eV_to_kJpermole)

        self.dHads = self.energy - self.energy_gas - self.Eslab
        self.heat_of_formation_0K = self.dHfgas + self.dHads * self.eV_to_kJpermole
        self.heat_of_formation_0K_units = 'kJ/mol'
        print(f"heat of formation adsorbate= {self.heat_of_formation_0K:.4} kJ/mol")
        print(f"heat of formation precursor= {self.dHfgas:.4} kJ/mol")
        print(f"heat of reaction precursor= {self.dHrxndftgas:.4} eV")
        print(f"DFT binding energy= {self.dHads:.4} eV")
        print(f"molecular mass= {self.adsorbate_mass} {self.adsorbate_mass_units}")

        return

    def _get_translation_thermo(self):
        # unpack the constants (not essential, but makes it easier to read)
        R = self.R
        kB = self.kB
        h = self.h
        amu = self.amu
        P_ref = self.P_ref
        m = self.adsorbate_mass
        pi = np.pi
        area = self.unit_cell_area_per_site
        sites = self.sites

        # initialize the arrays for the partition function, entropy, enthalpy,
        # and heat capacity.
        Q_trans = np.ones(len(self.temperature))
        S_trans = np.zeros(len(self.temperature))
        dH_trans = np.zeros(len(self.temperature))
        Cp_trans = np.zeros(len(self.temperature))

        if self.twoD_gas:
            print("switching to 2D-gas for 2 lowest modes for %s" % self.name)
            # cycle through each temperature
            for (i, T) in enumerate(self.temperature):
                # partition function is: (2*pi*mass*kB*T/h**2)^(2/2) * area
                if (1 == 0):  # 3D gas, really here just for inspiration
                    V = kB * T / P_ref
                    Q_trans[i] = (2 * pi * m * amu * kB * T / h ** 2) ** (1.5) * V
                    S_trans[i] = R * (2.5 + np.log(Q_trans[i]))  #
                    Cp_trans[i] = R * 2.5  # NOTE: Cp = Cv + R
                    dH_trans[i] = R * 2.5 * T
                else:  # surface
                    if (1 == 0):  # Campbell + Arnadottir
                        V = kB * T / P_ref
                        Q_trans[i] = (2 * pi * m * amu * kB * T / h ** 2) ** (1.0) * V ** 0.66667
                        S_trans[i] = R * (2.0 + np.log(Q_trans[i]))
                        Cp_trans[i] = R * 1.66667  # NOTE: Cp = Cv + 2/3R
                        dH_trans[i] = R * 1.66667 * T

                    else:  # area is not a function of temperature (This is what we use for our calculations)
                        Q_trans[i] = (2 * pi * m * amu * kB * T / h ** 2) * area * sites
                        S_trans[i] = R * (2.0 + np.log(Q_trans[i]))
                        Cp_trans[i] = R * 1.0  # NOTE: Cp = Cv
                        dH_trans[i] = R * 1.0 * T

                        # add the results to the thermo object
        self.Q_trans = Q_trans
        self.S_trans = S_trans
        self.dH_trans = dH_trans
        self.Cp_trans = Cp_trans

        return

    def _get_vibrational_thermo(self):
        units = 1.0
        units *= self.h * self.c / self.kB * self.invcm_to_invm  # K * cm
        amu = self.amu
        kB = self.kB
        h = self.h
        P_ref = self.P_ref
        mass = float(self.adsorbate_mass)

        # initialize the arrays for the partition function, entropy, enthalpy,
        # and heat capacity.
        Q_vib = np.ones(len(self.temperature))
        S_vib = np.zeros(len(self.temperature))
        dH_vib = np.zeros(len(self.temperature))
        Cv_vib = np.zeros(len(self.temperature))

        for (t, temp) in enumerate(self.temperature):
            for (n, nu) in enumerate(self.freqs):
                if self.twoD_gas == True and n <= 1:  # skip the first two if we do 2D gas
                    # do nothing!
                    Q_vib[t] *= 1.0
                    S_vib[t] += 0.0
                    dH_vib[t] += 0.0
                    Cv_vib[t] += 0.0
                else:
                    # where did these equations come from
                    x = nu * units / temp  # cm^-1 * K cm / K = dimensionless
                    # I think this is to collect the ZPE
                    Q_vib[t] *= 1.0 / (1.0 - np.exp(- x))
                    # This one I found in thermochem source
                    S_vib[t] += -np.log(1.0 - np.exp(- x)) + x * np.exp(- x) / (1.0 - np.exp(- x))
                    # This one I also found
                    dH_vib[t] += x * np.exp(- x) / (1.0 - np.exp(- x))
                    # Where is this from?
                    Cv_vib[t] += x ** 2.0 * np.exp(- x) / (1.0 - np.exp(- x)) ** 2.0
            S_vib[t] *= self.R
            dH_vib[t] *= self.R * temp
            Cv_vib[t] *= self.R

        # add the results to the thermo object
        self.Q_vib = Q_vib
        self.S_vib = S_vib
        self.dH_vib = dH_vib
        self.Cv_vib = Cv_vib  # NOTE: the correction from Cv to Cp is handled in the translation partition function.
        # if the molecule is tightly bound and thus the 2D-gas is not used,
        # then we assume that Cp=Cv for the adsorbate.

        return

    def _fit_NASA(self):
        R = self.R
        heat_capacity = self.Cp
        reference_enthalpy = self.H[0]
        reference_entropy = self.S[0]
        T_switch = self.T_switch
        temperature = self.temperature
        i_switch = -1
        for i in range(len(temperature)):
            if temperature[i] == T_switch:
                i_switch = i
        if i_switch == -1:
            print("We have a problem! Cannot find switching temperature")

        # start by creating the independent variable matrix for the low-temperature fit
        YT = np.array(
            [np.ones(len(temperature[:i_switch + 1])), temperature[:i_switch + 1],
             temperature[:i_switch + 1] ** 2.0,
             temperature[:i_switch + 1] ** 3.0, temperature[:i_switch + 1] ** 4.0],
            dtype=np.float64)  # this is transpose of our Y
        Y = YT.transpose()  # this is the desired Y

        b = heat_capacity[:i_switch + 1] / R
        a_low = np.linalg.lstsq(Y, b, rcond=None)[0]

        T_ref = 298.15
        # now determine the enthalpy coefficient for the low-T region
        subtract = a_low[0] + a_low[1] / 2.0 * T_ref + a_low[2] / 3.0 * T_ref ** 2.0 + a_low[
            3] / 4.0 * T_ref ** 3.0 + \
                   a_low[4] / 5.0 * T_ref ** 4.0
        a_low = np.append(a_low, reference_enthalpy / R - subtract * T_ref)
        # now determine the entropy coefficient for the low-T region
        subtract = a_low[0] * np.log(T_ref) + a_low[1] * T_ref + a_low[2] / 2.0 * T_ref ** 2.0 + a_low[
            3] / 3.0 * T_ref ** 3.0 + a_low[4] / 4.0 * T_ref ** 4.0
        a_low = np.append(a_low, reference_entropy / R - subtract)

        #
        # NOW SWITCH TO HIGH-TEMPERATURE REGIME!
        #
        T_ref = T_switch
        # compute the heat capacity, enthalpy, and entropy at the switching point
        Cp_switch = a_low[0] + a_low[1] * T_ref + a_low[2] * T_ref ** 2.0 + a_low[3] * T_ref ** 3.0 + a_low[
            4] * T_ref ** 4.0
        H_switch = a_low[0] * T_ref + a_low[1] / 2.0 * T_ref ** 2.0 + a_low[2] / 3.0 * T_ref ** 3.0 + a_low[
            3] / 4.0 * T_ref ** 4.0 + a_low[4] / 5.0 * T_ref ** 5.0 + a_low[5]
        S_switch = a_low[0] * np.log(T_ref) + a_low[1] * T_ref + a_low[2] / 2.0 * T_ref ** 2.0 + a_low[
            3] / 3.0 * T_ref ** 3.0 + a_low[4] / 4.0 * T_ref ** 4.0 + a_low[6]

        # now repeat the process for the high-temperature regime
        a_high = [0.0]
        YT = np.array([temperature[i_switch:], temperature[i_switch:] ** 2.0, temperature[i_switch:] ** 3.0,
                       temperature[i_switch:] ** 4.0], dtype=np.float64)  # this is transpose of our Y
        Y = YT.transpose()  # this is the desired Y

        b = heat_capacity[i_switch:] / R - Cp_switch
        a_high = np.append(a_high, np.linalg.lstsq(Y, b, rcond=None)[0])
        a_high[0] = Cp_switch - (
                a_high[0] + a_high[1] * T_switch + a_high[2] * T_switch ** 2.0 + a_high[3] * T_switch ** 3.0 +
                a_high[
                    4] * T_switch ** 4.0)

        a_high = np.append(a_high, H_switch - (
                a_high[0] + a_high[1] / 2.0 * T_ref + a_high[2] / 3.0 * T_ref ** 2.0 + a_high[
            3] / 4.0 * T_ref ** 3.0 +
                a_high[4] / 5.0 * T_ref ** 4.0) * T_ref)
        a_high = np.append(a_high, S_switch - (
                a_high[0] * np.log(T_ref) + a_high[1] * T_ref + a_high[2] / 2.0 * T_ref ** 2.0 + a_high[
            3] / 3.0 * T_ref ** 3.0 + a_high[4] / 4.0 * T_ref ** 4.0))

        # Check to see if there is a discontinuity
        if (1 == 0):
            print("\ncheck for discontinuities:")
            cp_low_Tswitch = a_low[0] + a_low[1] * T_switch + a_low[2] * T_switch ** 2.0 + a_low[
                3] * T_switch ** 3.0 + \
                             a_low[4] * T_switch ** 4.0
            cp_high_Tswitch = a_high[0] + a_high[1] * T_switch + a_high[2] * T_switch ** 2.0 + a_high[
                3] * T_switch ** 3.0 + \
                              a_high[4] * T_switch ** 4.0
            H_low_Tswitch = a_low[0] * T_switch + a_low[1] / 2.0 * T_switch ** 2.0 + a_low[
                2] / 3.0 * T_switch ** 3.0 + \
                            a_low[3] / 4.0 * T_switch ** 4.0 + a_low[4] / 5.0 * T_switch ** 5.0 + a_low[5]
            H_high_Tswitch = a_high[0] * T_switch + a_high[1] / 2.0 * T_switch ** 2.0 + a_high[
                2] / 3.0 * T_switch ** 3.0 + \
                             a_high[3] / 4.0 * T_switch ** 4.0 + a_high[4] / 5.0 * T_switch ** 5.0 + a_high[5]
            S_low_Tswitch = a_low[0] * np.log(T_switch) + a_low[1] * T_switch + a_low[2] / 2.0 * T_switch ** 2.0 + \
                            a_low[
                                3] / 3.0 * T_switch ** 3.0 + a_low[4] / 4.0 * T_switch ** 4.0 + a_low[6]
            S_high_Tswitch = a_high[0] * np.log(T_switch) + a_high[1] * T_switch + a_high[
                2] / 2.0 * T_switch ** 2.0 + \
                             a_high[3] / 3.0 * T_switch ** 3.0 + a_high[4] / 4.0 * T_switch ** 4.0 + a_high[6]

            print("discontinuity at T_switch for Cp/R is %.4F" % (cp_low_Tswitch - cp_high_Tswitch))
            print("discontinuity at T_switch for H/R is %.4F" % (H_low_Tswitch - H_high_Tswitch))
            print("discontinuity at T_switch for S/R is %.4F" % (S_low_Tswitch - S_high_Tswitch))

            # line = '\n\t !cut and paste this value into the cti file!\n'
        line = '\tthermo = (\n'
        line += "\t\tNASA( [%.1F, %.1F], [%.8E, %.8E,\n \t\t %.8E, %.8E, %.8E,\n \t\t %.8E, %.8E]), \n" % (
            300.0, 1000.0, a_low[0], a_low[1], a_low[2], a_low[3], a_low[4], a_low[5], a_low[6])
        line += "\t\tNASA( [%.1F, %.1F], [%.8E, %.8E,\n \t\t %.8E, %.8E, %.8E,\n \t\t %.8E, %.8E]), \n" % (
            1000.0, max(temperature), a_high[0], a_high[1], a_high[2], a_high[3], a_high[4], a_high[5], a_high[6])
        line += "\t\t ),\n"

        rmg_line = '\tthermo = NASA(\n'
        rmg_line += '\t\tpolynomials = [\n'
        rmg_line += '\t\tNASAPolynomial(coeffs=[\n'
        rmg_line += "\t\t %.8E, %.8E, %.8E, %.8E, \n \t\t %.8E, %.8E, %.8E], Tmin=(%.1F,'K'), Tmax=(%.1F, 'K')), \n" % (
            a_low[0], a_low[1], a_low[2], a_low[3], a_low[4], a_low[5], a_low[6], 298.0, 1000.0)
        rmg_line += '\t\tNASAPolynomial(coeffs=[\n'
        rmg_line += "\t\t %.8E, %.8E, %.8E, %.8E, \n \t\t %.8E, %.8E, %.8E], Tmin=(%.1F,'K'), Tmax=(%.1F, 'K')), \n" % (
            a_high[0], a_high[1], a_high[2], a_high[3], a_high[4], a_high[5], a_high[6], 1000.0, max(temperature))
        rmg_line += "\t\t ],\n"
        rmg_line += "\t\t Tmin = (%.1F, 'K'),\n" % (298.0)
        rmg_line += "\t\t Tmax = (%.1F, 'K'),\n" % (max(temperature))
        rmg_line += "\t),\n"
        rmg_line += f'\tmetal = "{self.cat_element}",\n'
        rmg_line += '\tfacet = "111",\n'

        self.thermo_lines = rmg_line
        # molecule.thermo_lines_for_rmg =
        self.a_low = a_low
        self.a_high = a_high

        return

    def _format_RMG_output(self):
        line = '\n'
        line += 'entry(\n    index = %s,\n' % (self.index)
        line += f'    label = "{self.name}",\n'
        line += '    molecule = \n"""\n%s\n""",\n' % (self.formatted_adj_list)
        line += self.thermo_lines
        line += f'    longDesc = u"""Calculated by {self.author} at {self.institution} using statistical mechanics using Pynta Thermo in postprocessing. \n'
        line += "                   These methods were based on approach by Blondal et. al in https://doi.org/10.1021/acscatal.2c03378. If you use this database in your work, please cite the publication mentioned above.\n"
        if self.twoD_gas:
            line += '\n            The two lowest frequencies, %.1F and %.1F %s, where replaced by the 2D gas model.' % (
                self.freqs[0], self.freqs[1], self.frequencies_units.replace("'", ""))
        line += '""",\n)\n'

        self.species_lines = line

        return
