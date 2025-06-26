from ase.io import read, write
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
from pynta.mol import split_adsorbed_structures
from pynta.geometricanalysis import *
from pynta.utils import to_dict
from pynta.coveragedependence import mol_to_atoms, process_calculation
from pysidt.sidt import *
import logging 

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
        return (th.get_helmholtz_energy(T+dT,verbose=False) - th.get_helmholtz_energy(T-dT,verbose=False))/(2*dT)
    elif isinstance(th,IdealGasThermo):
        return (th.get_enthalpy(T+dT,verbose=False) - th.get_enthalpy(T-dT,verbose=False))/(2*dT)

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
):
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
           
class GasConfiguration:
    """Base class for calculating statmech of gas phase species
    """
    def __init__(self,
                 atoms,
                 mol,
                 vibdata,
                 name,
                 c_ref=0,
                 o_ref=0,
                 h_ref=0,
                 n_ref=0,
                 valid=True):
        """
        Args:
            atoms (_type_): ase.Atoms object for configuration
            mol (_type_): Molecule object corresponding to configuration
            vibdata (_type_): ase VibrationData object corresponding to configuration
            name (_type_): name of configuration
            c_ref (int, optional): ZPE corrected reference energy for C (CH4). Defaults to 0.
            o_ref (int, optional): ZPE corrected reference energy for O (H2O). Defaults to 0.
            h_ref (int, optional): ZPE corrected reference energy for H (H2). Defaults to 0.
            n_ref (int, optional): ZPE corrected reference energy for N (NH3). Defaults to 0.
            valid (bool, optional): whether the configuration is valid for thermochemistry/kinetics. Defaults to True.
        """
        
        # start by defining some physical constants
        self.R = 8.3144621  # ideal Gas constant in J/mol-K
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
        
        self.Eref = {'CH4': c_ref, 'H2O': o_ref, 'H2': h_ref, 'NH3': n_ref}  # eV
        self.dHfatct = {'CH4': -66.540, 'H2O': -238.903, 'H2': 0,
                        'NH3': -38.563}  # heats of formation (0K) in kJ/mol from the ATcT database for the reference species, version 1.202
        self.dHrxnatct = {'H2-2H': 432.0680600, 'O2-2O': 493.6871,
                          'N2-2N': 941.165}  # Heats of the dissociation reactions in the gas phase from the ATcT database, version 1.202
        
        self.atoms = atoms 
        self.mol = mol
        self.vibdata = vibdata 
        self.name = name
        self.valid = valid
        temperature = [298.15]  # NOTE 298.15 must be first for the NASA polynomial routine to work!
        self.T_low = 300.0
        self.T_switch = 1000.0
        self.T_high = 2000.0
        dT = 10.0  # temperature increment
        self.temperature = np.append(temperature, np.arange(self.T_low, self.T_high + dT, dT))

        formula = atoms.get_chemical_formula(mode='all')
        self.composition = {'H': formula.count('H'), 'C': formula.count('C'),
                            'N': formula.count('N'), 'O': formula.count('O')}

        self.energy_gas = 0 #cancels out later in calculations

        self.DFT_energy = self.atoms.get_potential_energy()
        self.ZPE_energy = self.vibdata.get_zero_point_energy()
        self.DFT_energy_gas_units = 'eV'
        self.DFT_energy_units = 'eV'
        self.ZPE_energy_units = 'eV'
        
    def compute_heat_of_formation(self):
        self.energy = self.DFT_energy + self.ZPE_energy

        self.dHrxndftgas = (self.energy - self.composition['C'] * self.Eref['CH4']
                            - self.composition['O'] * self.Eref['H2O']
                            - self.composition['N'] * self.Eref['NH3']
                            - (self.composition['H'] / 2 - 2 * self.composition['C'] - self.composition[
                    'O'] - 3 / 2 * self.composition['N']) * self.Eref['H2']) #eV
        self.dHfgas = (self.composition['C'] * self.dHfatct['CH4']
                       + self.composition['O'] * self.dHfatct['H2O']
                       + self.composition['N'] * self.dHfatct['NH3']
                       + (self.composition['H'] / 2 - 2 * self.composition['C'] - self.composition[
                    'O'] - 3 / 2 * self.composition['N']) * self.dHfatct['H2']
                       + self.dHrxndftgas * self.eV_to_kJpermole) * 1000.0 #J/mol

        self.heat_of_formation_0K = self.dHfgas #J/mol
    
    def to_dict(self):
        return to_dict(self)
    
    def fit_NASA(self):
        wh = Wilhoit().fit_to_data(Tdata=self.temperature,Cpdata=np.array(self.Cp),Cp0=self.Cp0,CpInf=self.Cpinf,H298=self.H[0],S298=self.S[0])

        self.nasa = wh.to_nasa(Tmin=298.15, Tmax=self.T_high, Tint=self.T_switch)
        
        self.thermo_lines = "thermo="+repr(self.nasa)+","
        
    def format_RMG_output(self):
        xyz = ""
        for i,a in enumerate(self.atoms):
            xyz += a.symbol+" "+str(self.atoms.positions[i][0])+" "+str(self.atoms.positions[i][1])+" "+str(self.atoms.positions[i][1])+"\n"
        
        line = '\n'
        line += 'entry(\n    index = {index},\n'
        line += f'    label = "{self.name}",\n'
        line += '    molecule = """\n%s""",\n' % (self.mol.to_adjacency_list())
        line += self.thermo_lines
        line += f'    \nlongDesc = u"""Calculated using Pynta https://doi.org/10.1021/acs.jcim.3c00948. \n'
        line += "                   Thermochemistry computed using approach from Blondal et. al in https://doi.org/10.1021/acscatal.2c03378. \nIf you use this library in your work, please cite the publications mentioned above.\n"
        line += "Hf298: {} [kcal/mol]\n".format(self.H[0]/4184.0)
        line += "Sf298: {} [cal/(mol-K)]\n".format(self.S[0]/4.184)
        line += xyz
        line += '""",\n)\n'

        self.rmg_species_text = line
    
    def create_RMG_header(self,
                          lib_name,
                          lib_short_desc,
                          lib_long_desc,
                          metal,facet):
        line = f'''#!/usr/bin/env python
# encoding: utf-8

name = "{lib_name}"
shortDesc = u"{lib_short_desc}"
longDesc = u"""{lib_long_desc}"""
#

entry(
    index = 1,
    label = "vacant",
    molecule =
"""
1 X  u0 p0 c0
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[
             0.000000000E+00,   0.000000000E+00,   0.000000000E+00,   0.000000000E+00,
             0.000000000E+00,   0.000000000E+00,   0.000000000E+00], Tmin=(298.0,'K'), Tmax=(1000.0, 'K')),
            NASAPolynomial(coeffs=[
             0.000000000E+00,   0.000000000E+00,   0.000000000E+00,   0.000000000E+00,
             0.000000000E+00,   0.000000000E+00,   0.000000000E+00], Tmin=(1000.0,'K'), Tmax=(3000.0, 'K')),
        ],
        Tmin = (298.0, 'K'),
        Tmax = (3000.0, 'K'),
    ),
    metal = "{metal}",
    facet = "{facet}",
)
        '''
        return line
        
    def run(self):
        self.compute_heat_of_formation()
        
        if len(self.atoms) == 1:
            geometry = 'monatomic'
        elif any([ I < 1.0e-3 for I in self.atoms.get_moments_of_inertia()]):
            geometry = "linear"
        else:
            geometry = "nonlinear"

        self.thermo = IdealGasThermo(np.real(self.vibdata.get_energies()),geometry,
                                        potentialenergy=(self.heat_of_formation_0K/1000.0)/self.eV_to_kJpermole-self.ZPE_energy,atoms=self.atoms,symmetrynumber=1,
                                        natoms=len(self.atoms),spin=(self.mol.multiplicity - 1)/2)
        self.G = []
        self.H = []
        self.S = []
        self.Cp = []
        for T in self.temperature:
            G = self.thermo.get_gibbs_energy(T,1.0e5,verbose=False)
            H = self.thermo.get_enthalpy(T,verbose=False)
            S = self.thermo.get_entropy(T,1.0e5,verbose=False)
            Cp = get_cp(self.thermo,T,dT=0.01)
            
            self.G.append(G*self.eV_to_kJpermole*1000.0)
            self.H.append(H*self.eV_to_kJpermole*1000.0)
            self.S.append(S*self.eV_to_kJpermole*1000.0)
            self.Cp.append(Cp*self.eV_to_kJpermole*1000.0)

        self.Cp0 = get_cp(self.thermo,10.0,dT=0.01)
        self.Cpinf = get_cp(self.thermo,1.0e7,dT=0.01)
        
        self.entropy_of_formation_298K = self.S[0]
        self.heat_of_formation_298K = self.H[0]
        
        self.fit_NASA()
        self.format_RMG_output()
        
class SurfaceConfiguration:
    """Base class for calculating statmech and thermodynamic properties for Pynta optimized
    adsorbates and TSs.

    parameters:
        atoms: ase atoms object of configuration
        slab: ase atoms object of slab
        admol, Molecule object representing the configuration on the slab
        vibdata, ase VibrationData object for the configuration
        name: adsorbate name generated from Pynta
        metal: metal of slab
        facet: facet of slab
        is_TS: if the configuration is a TS
        sites_per_cell: We set this to 1 by default, assuming that the 2D lattice gas
            can freely diffuse over the entire area of our unit cell.
        c_ref: ZPE corrected gas phase reference for carbon (CH4) [eV]
        o_ref: ZPE corrected gas phase reference for oxygen (H2O) [eV]
        h_ref: ZPE corrected gas phase reference for hydrogen(H2) [eV]
        n_ref: ZPE corrected gas phase reference for nitrogen (NH3) [eV]
        cutoff_freq: Cutoff used to determine whether to use 2D lattice gas
        valid: if the configuration is valid for thermochemistry and kinetics calculations
        valid_info: information dictionary about the validity
        mol: Molecule object representing the configuration without slab information (unnecessary if admol is specified)

    Much of this code was copied/adapted from the Goldsmith group at Brown
    University. If you use this code, please cite the following work by Blondal et al.:
    https://doi.org/10.1021/acscatal.2c03378

    """
    def __init__(self,
                 atoms,
                 slab,
                 admol,
                 vibdata,
                 name,
                 metal,
                 facet,
                 is_TS=False,
                 sites_per_cell=1,
                 c_ref=0.0,
                 o_ref=0.0,
                 h_ref=0.0,
                 n_ref=0.0,
                 cutoff_freq=100.0,
                 valid=True,
                 valid_info=None,
                 mol=None):

        # start by defining some physical constants
        self.R = 8.3144621  # ideal Gas constant in J/mol-K
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

        self.atoms = atoms
        self.slab = slab
        self.admol = admol
        if mol:
            self.mol = mol 
        else:
            self.mol = split_adsorbed_structures(admol)[0]
        self.vibdata = vibdata 
        self.metal = metal
        self.facet = facet 
        self.nslab = len(self.slab)
        self.is_TS = is_TS
        self.valid = valid
        self.valid_info = valid_info
        self.Eref = {'CH4': c_ref, 'H2O': o_ref, 'H2': h_ref, 'NH3': n_ref}  # eV
        self.Eslab = self.slab.get_potential_energy()  # DFT energy of the relaxed bare slab in eV

        self.site_num = len(self.mol.get_surface_sites())
        if self.site_num == 0:
            self.gasphase = True
        ad_atoms = atoms[self.nslab:]
        formula = ad_atoms.get_chemical_formula(mode='all')

        self.composition = {'H': formula.count('H'), 'C': formula.count('C'),
                            'N': formula.count('N'), 'O': formula.count('O')}

        self.N_ad_atoms = len(ad_atoms)
        self.adsorbate_mass = sum(ad_atoms.get_masses())
        self.adsorbate_mass_units = 'amu'

        self.energy_gas = 0 #cancels out in calculations

        self.DFT_energy = self.atoms.get_potential_energy()
        self.ZPE_energy = self.vibdata.get_zero_point_energy()
        self.DFT_energy_gas_units = 'eV'
        self.DFT_energy_units = 'eV'
        self.ZPE_energy_units = 'eV'

        self.adjacency_list = self.mol.to_adjacency_list()
        self.formatted_adj_list = self.adjacency_list

        self.freqs = [np.complex128(x).real for x in self.vibdata.get_frequencies()]
        self.name = name

        self.frequencies_units = 'cm-1'

        if not self.is_TS and self.freqs[1].real < self.cutoff_frequency:
            self.twoD_gas = True

        temperature = [298.15]  # NOTE 298.15 must be first for the NASA polynomial routine to work!
        self.T_low = 300.0
        self.T_high = 2000.0
        self.dT = 10.0  # temperature increment
        self.temperature = np.append(temperature, np.arange(self.T_low, self.T_high + self.dT, self.dT))

        self.cat_element = metal
              
    def to_dict(self):
        return to_dict(self)
        
    def compute_heat_of_formation(self):
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
        self.heat_of_formation_0K = (self.dHfgas + self.dHads * self.eV_to_kJpermole) * 1000.0
        self.heat_of_formation_0K_units = 'J/mol'

    def get_translation_thermo(self):
        # unpack the constants (not essential, but makes it easier to read)
        R = self.R
        kB = self.kB
        h = self.h
        amu = self.amu
        P_ref = self.P_ref
        m = self.adsorbate_mass
        pi = np.pi
        area = self.unit_cell_area_per_site
        sites = self.site_num

        # initialize the arrays for the partition function, entropy, enthalpy,
        # and heat capacity.
        Q_trans = np.ones(len(self.temperature)+2)
        S_trans = np.zeros(len(self.temperature)+2)
        dH_trans = np.zeros(len(self.temperature)+2)
        Cp_trans = np.zeros(len(self.temperature)+2)

        if self.twoD_gas:
            # cycle through each temperature
            for (i, T) in enumerate([10.0]+self.temperature.tolist()+[1.0e7]):
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
        self.Cp0_trans = Cp_trans[0] #all J-mol-K
        self.Cpinf_trans = Cp_trans[-1]
        self.Q_trans = Q_trans[1:-1]
        self.S_trans = S_trans[1:-1]
        self.dH_trans = dH_trans[1:-1]
        self.Cp_trans = Cp_trans[1:-1]

    def get_vibrational_thermo(self):
        units = self.h * self.c / self.kB * self.invcm_to_invm  # K * cm
        amu = self.amu
        kB = self.kB
        h = self.h
        P_ref = self.P_ref
        mass = float(self.adsorbate_mass)

        # initialize the arrays for the partition function, entropy, enthalpy,
        # and heat capacity.
        Q_vib = np.ones(len(self.temperature)+2)
        S_vib = np.zeros(len(self.temperature)+2)
        dH_vib = np.zeros(len(self.temperature)+2)
        Cv_vib = np.zeros(len(self.temperature)+2)

        for (t, temp) in enumerate([10.0]+self.temperature.tolist()+[1.0e7]):
            for (n, nu) in enumerate(self.freqs):
                if (self.twoD_gas == True and n <= 1) or nu == 0:  # skip the first two if we do 2D gas and skip if frequency was imaginary (no real component)
                    continue
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

        self.Cv0_vib = Cv_vib[0] #all J-mol-K
        self.Cvinf_vib = Cv_vib[-1]
        # add the results to the thermo object
        self.Q_vib = Q_vib[1:-1]
        self.S_vib = S_vib[1:-1]
        self.dH_vib = dH_vib[1:-1]
        self.Cv_vib = Cv_vib[1:-1]  # NOTE: the correction from Cv to Cp is handled in the translation partition function.
        # if the molecule is tightly bound and thus the 2D-gas is not used,
        # then we assume that Cp=Cv for the adsorbate.

        return

    def fit_NASA(self):
        wh = Wilhoit(comment='metal = {metal},facet = {facet}'.format(metal=self.metal,facet=self.facet)).fit_to_data(Tdata=self.temperature,Cpdata=self.Cp,Cp0=self.Cp0,CpInf=self.Cpinf,H298=self.H[0],S298=self.S[0])

        self.nasa = wh.to_nasa(Tmin=298.15, Tmax=self.T_high, Tint=self.T_switch)
        
        self.thermo_lines = "thermo="+repr(self.nasa)+","

    def format_RMG_output(self):
        xyz = str(len(self.atoms)) + "\n"
        for s,(x,y,z) in zip(self.atoms.symbols,self.atoms.positions):
            xyz += s + " " + str(x) + " " + str(y) + " " + str(z) + "\n"
        line = '\n'
        line += 'entry(\n    index = {index},\n'
        line += f'    label = "{self.name}",\n'
        line += '    molecule = \n"""\n%s""",\n' % (self.formatted_adj_list)
        line += self.thermo_lines
        line += f'    \nlongDesc = u"""Calculated using Pynta https://doi.org/10.1021/acs.jcim.3c00948. \n'
        line += "                   Thermochemistry computed using approach from Blondal et. al in https://doi.org/10.1021/acscatal.2c03378. \nIf you use this library in your work, please cite the publications mentioned above.\n"
        line += "Hf298: {} [kcal/mol]\n".format(self.H[0]/4184.0)
        line += "Sf298: {} [cal/(mol-K)]\n".format(self.S[0]/4.184)
        line += xyz
        if self.twoD_gas:
            line += '\n            The two lowest frequencies, %.1F and %.1F %s, where replaced by the 2D gas model.' % (
                self.freqs[0], self.freqs[1], self.frequencies_units.replace("'", ""))
        line += '""",\n)\n'

        self.rmg_species_text = line

        return

    def create_RMG_header(self,
                          lib_name,
                          lib_short_desc,
                          lib_long_desc):
        line = f'''#!/usr/bin/env python
# encoding: utf-8

name = "{lib_name}"
shortDesc = u"{lib_short_desc}"
longDesc = u"""{lib_long_desc}"""
#

entry(
    index = 1,
    label = "vacant",
    molecule =
"""
1 X  u0 p0 c0
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[
             0.000000000E+00,   0.000000000E+00,   0.000000000E+00,   0.000000000E+00,
             0.000000000E+00,   0.000000000E+00,   0.000000000E+00], Tmin=(298.0,'K'), Tmax=(1000.0, 'K')),
            NASAPolynomial(coeffs=[
             0.000000000E+00,   0.000000000E+00,   0.000000000E+00,   0.000000000E+00,
             0.000000000E+00,   0.000000000E+00,   0.000000000E+00], Tmin=(1000.0,'K'), Tmax=(3000.0, 'K')),
        ],
        Tmin = (298.0, 'K'),
        Tmax = (3000.0, 'K'),
    ),
    metal = "{self.cat_element}",
    facet = "{self.facet}",
)
        '''
        return line

    def calculate_thermochemistry(self):
        # call the subroutine for the translational and vibrational partition functions
        self.get_translation_thermo()
        self.get_vibrational_thermo()
        temperature = self.temperature

        # Atom corrections from Thermochemistry chapter by Ruscic and Bross
        # Page 75 of DOI = https://doi.org/10.1016/B978-0-444-64087-1.00001-2
        h_correction = 4.234  # kJ/mol. enthalpy_H(298) - enthalpy_H(0)
        c_correction = 1.051  # kJ/mol. enthalpy_C(298) - enthalpy_C(0)
        n_correction = 4.335  # kJ/mol. enthalpy_N(298) - enthalpy_N(0)
        o_correction = 4.340  # kJ/mol. enthalpy_O(298) - enthalpy_O(0)

        # Still not certain why these are needed. Need to read book
        self.heat_of_formation_correction = 0.0
        self.heat_of_formation_correction += self.composition['H'] * h_correction * 1000.0 #J/mol
        self.heat_of_formation_correction += self.composition['C'] * c_correction * 1000.0
        self.heat_of_formation_correction += self.composition['N'] * n_correction * 1000.0
        self.heat_of_formation_correction += self.composition['O'] * o_correction * 1000.0

        # note that the partition function is the product of the individual terms,
        # whereas the thermodynamic properties are additive
        self.Q = self.Q_trans * self.Q_vib
        self.S = self.S_trans + self.S_vib
        self.dH = self.dH_trans + self.dH_vib
        self.Cp = self.Cp_trans + self.Cv_vib  # see comments in each section regarding Cp vs Cv
        self.Cp0 = self.Cp0_trans + self.Cv0_vib
        self.Cpinf = self.Cpinf_trans + self.Cvinf_vib
        self.heat_of_formation_298K = self.heat_of_formation_0K + self.dH[
            0] - self.heat_of_formation_correction
        self.H = self.heat_of_formation_298K + self.dH - self.dH[0]
        self.G = self.H - self.temperature*self.S
        self.entropy_of_formation_298K = self.S[0]
        # now that we've computed the thermo properties, go ahead and fit them to a NASA polynomial
        self.fit_NASA()
        self.format_RMG_output()
        return

    def run(self):
        self.compute_heat_of_formation()
        self.calculate_thermochemistry()
        
class Kinetics:
    """
    Class for computing the kinetic parameters of reactions
    """
    def __init__(self,reactants,transition_state,reactant_mol,product_mol,site_density,
                 metal,facet,products=None,family_comment=None,valid=True):
        """
        Args:
            reactants: list of GasConfiguration/SurfaceConfiguration reactants for the reaction
            transition_state: SurfaceConfiguration object corresponding to the transition state
            reactant_mol: Molecule object representing the tagged reactants (not slab resolved)
            product_mol: Molecule object representing the tagged products (not slab resolved)
            site_density: site density of the slab in molecules/m^2
            metal: metal of the slab
            facet: facet of the slab
            products: list of GasConfiguration/SurfaceConfiguration product for the reaction. Defaults to None.
            family_comment: string containing information about the reaction, transition state and reactants/products. Defaults to None.
            valid: if the class is valid for calculating kinetics
        """
        self.metal = metal
        self.facet = facet
        temperature = [298.15]
        self.T_low = 300.0
        self.T_high = 2000.0
        self.dT = 10.0  # temperature increment
        self.temperature = np.append(temperature, np.arange(self.T_low, self.T_high + self.dT, self.dT))
        self.reactants = reactants
        self.transition_state = transition_state
        self.site_density = site_density
        self.products = products
        self.reactant_mol = reactant_mol
        self.product_mol = product_mol
        self.valid = valid
        self.valid_info = self.transition_state.valid_info
        self.rnum = len(reactant_mol.split())
        if self.products:
            self.pnum = len(product_mol.split())
        else:
            self.pnum = None
        self.reactant_names = [x.name for x in self.reactants]
        if self.products:
            self.product_names = [x.name for x in self.products]
        else:
            self.product_names = None
        self.reaction_str = " + ".join(self.reactant_names+["vacantX" for i in range(self.rnum-len(self.reactant_names))])+" <=> "+" + ".join(self.product_names+["vacantX" for i in range(self.pnum-len(self.product_names))])
        self.family_comment = family_comment

    def to_dict(self):
        return to_dict(self)
    
    def calculate_barrier(self):
        bar_f = self.transition_state.energy - self.transition_state.Eslab
        bar_r = bar_f
        for r in self.reactants:
            if isinstance(r,SurfaceConfiguration):
                bar_f -= r.energy - r.Eslab
            elif isinstance(r,GasConfiguration):
                bar_f -= r.energy
        for p in self.products:
            if isinstance(p,SurfaceConfiguration):
                bar_r -= p.energy - p.Eslab
            elif isinstance(p,GasConfiguration):
                bar_r -= p.energy

        self.barrier_f = bar_f
        self.barrier_r = bar_r
        
    def calculate_kinetic_parameters(self):
        kB = ase.units.kB
        R = 8.314 #J/(mol K)
        h = 6.582119569e-16 # eV * s
        c0 = 2.4282*10**22*1000.0 #molecules/m^3

        kfs = []
        krevs = []
        for i,T in enumerate(self.temperature):
            factor_f = self.site_density
            Lunits_f = -2
            factor_r = self.site_density
            Lunits_r = -2
            
            GTS = self.transition_state.G[i]
            
            GR = 0.0
            for r in self.reactants:
                GR += r.G[i]
                if isinstance(r,SurfaceConfiguration):
                    factor_f /= self.site_density
                    Lunits_f += 2
                elif isinstance(r,GasConfiguration):
                    factor_f /= c0
                    Lunits_f += 3

            for i in range(self.rnum-len(self.reactants)):
                factor_f /= self.site_density
                Lunits_f += 2
            
            kf = kB*T/h * np.exp(-(GTS-GR)/(R*T)) * factor_f
            kfs.append(kf)
            
            if self.products:
                GP = 0.0
                for p in self.products:
                    GP += p.G[i]
                    if isinstance(p,SurfaceConfiguration):
                        factor_r /= self.site_density
                        Lunits_r += 2
                    elif isinstance(p,GasConfiguration):
                        factor_r /= c0
                        Lunits_r += 3
                for i in range(self.pnum-len(self.products)):
                    factor_r /= self.site_density
                    Lunits_r += 2
                
                krev = kB*T/h * np.exp(-(GTS-GP)/(R*T)) * factor_r
                krevs.append(krev)
                
        if Lunits_f == 0 and self.rnum - 1 == 0:
            self.arr_f = SurfaceArrhenius().fit_to_data(self.temperature,np.array(kfs),"s^-1")
        elif self.rnum - 1 == 1:
            self.arr_f = SurfaceArrhenius().fit_to_data(self.temperature,np.array(kfs),"m^{Lunits}/(molecule*s)".format(Lunits=Lunits_f))
        else:
            self.arr_f = SurfaceArrhenius().fit_to_data(self.temperature,np.array(kfs),"m^{Lunits}/(molecules^{mols}*s)".format(Lunits=Lunits_f,mols=self.rnum-1))
        
        self.kfs = kfs
        
        if self.products:
            if Lunits_r == 0 and self.pnum - 1 == 0:
                self.arr_r = SurfaceArrhenius().fit_to_data(self.temperature,np.array(krevs),"s^-1")
            elif self.pnum - 1 == 1:
                self.arr_r = SurfaceArrhenius().fit_to_data(self.temperature,np.array(krevs),"m^{Lunits}/(molecule*s)".format(Lunits=Lunits_r))
            else:
                self.arr_r = SurfaceArrhenius().fit_to_data(self.temperature,np.array(krevs),"m^{Lunits}/(molecules^{mols}*s)".format(Lunits=Lunits_r,mols=self.pnum-1))
        
        self.krevs = krevs 
            
    def generate_kinetics_entry(self):
        xyz = str(len(self.transition_state.atoms)) + "\n"
        for s,(x,y,z) in zip(self.transition_state.atoms.symbols,self.transition_state.atoms.positions):
            xyz += s + " " + str(x) + " " + str(y) + " " + str(z) + "\n"
        txt = '''entry(
    index = {index},
    label = "{rlabel}",
    kinetics = {kinetics},
    shortDesc = u"{short_desc}",
    longDesc = u"""{long_desc}""",
    metal = "{metal}",
    facet = "{facet}",
)'''
        txt = txt.format(rlabel=self.reaction_str,
                        kinetics=repr(self.arr_f),
                        short_desc="Computed using Pynta",
                        long_desc="Computed using Pynta\n" + self.family_comment + "\n" + xyz,
                        metal=self.metal,
                        facet=self.facet[3:],index="{index}")
        self.rmg_kinetics_text = txt

    def run(self):
        self.calculate_barrier()
        self.calculate_kinetic_parameters()
        self.generate_kinetics_entry()
        
    def get_gibbs_energy_reaction(self,T):
        dG = 0.0
        for th in self.reactants:
            dG -= th.nasa.get_free_energy(T)
        for th in self.products:
            dG += th.nasa.get_free_energy(T)
        return dG
        
    def get_entropy_reaction(self,T):
        dS = 0.0
        for th in self.reactants:
            dS -= th.nasa.get_entropy(T)
        for th in self.products:
            dS += th.nasa.get_entropy(T)
        return dS

    def get_enthalpy_reaction(self,T):
        dH = 0.0
        for th in self.reactants:
            dH -= th.nasa.get_enthalpy(T)
        for th in self.products:
            dH += th.nasa.get_enthalpy(T)
        return dH
    
    def create_RMG_header(self,
                          lib_name,
                          lib_short_desc,
                          lib_long_desc):
        line = f'''#!/usr/bin/env python
# encoding: utf-8

name = "{lib_name}"
shortDesc = u"{lib_short_desc}"
longDesc = u"""{lib_long_desc}"""
'''
        return line
    
def get_species(path,adsorbates_path,metal,facet,slab,sites,site_adjacency,nslab,c_ref=0.0,o_ref=0.0,h_ref=0.0,
               n_ref=0.0):
    """
    Generate a dictionary of GasConfiguration/SurfaceConfiguration objects corresponding to a Pynta 
    species calculation
    Args:
        path: directory of the corresponding adsorbate/species in the Pynta run
        adsorbates_path: Adsorbates directory of the Pynta run
        metal: metal of the slab
        facet: facet of the slab
        slab: ase.Atoms object for the slab
        sites: list of sites correspondingn to the slab
        site_adjacency: site_adjacency corresponding to the slab
        nslab: number of atoms in the slab
        c_ref: ZPE corrected gas phase reference for carbon (CH4) [eV]
        o_ref: ZPE corrected gas phase reference for oxygen (H2O) [eV]
        h_ref: ZPE corrected gas phase reference for hydrogen(H2) [eV]
        n_ref: ZPE corrected gas phase reference for nitrogen (NH3) [eV]

    Returns:
        dictionary mapping string index to GasConfiguration/SurfaceConfiguration objects
    """
    with open(os.path.join(path,"info.json"),'r') as f:
        info = json.load(f)
        name = info["name"]
        mol = Molecule().from_adjacency_list(info["adjlist"])
        keep_binding_vdW_bonds = False
        keep_vdW_surface_bonds = False
        for bd in mol.get_all_edges():
            if bd.order == 0:
                if bd.atom1.is_surface_site() or bd.atom2.is_surface_site():
                    keep_binding_vdW_bonds = True
                    m = mol.copy(deep=True)
                    b = m.get_bond(m.atoms[mol.atoms.index(bd.atom1)],m.atoms[mol.atoms.index(bd.atom2)])
                    m.remove_bond(b)
                    out = m.split()
                    if len(out) == 1: #vdW bond is not only thing connecting adsorbate to surface
                        keep_vdW_surface_bonds = True

    
    gas_phase = len(mol.get_surface_sites()) == 0
    species_dict = dict()
    for ind in os.listdir(path):
        if os.path.exists(os.path.join(path,ind,"vib.json_vib.json")):
            dopt = os.path.join(path,ind,ind+".xyz")
            dvib = os.path.join(path,ind,"vib.json_vib.json")
            atoms = read(dopt)
            if not gas_phase:
                vibdata = get_vibdata(dopt,dvib,nslab)
                try: 
                    admol,neighbor_sites,ninds = generate_adsorbate_2D(atoms, sites, site_adjacency, nslab, max_dist=np.inf, 
                          keep_binding_vdW_bonds=keep_binding_vdW_bonds, keep_vdW_surface_bonds=keep_vdW_surface_bonds)
                    split_structs = split_adsorbed_structures(admol)
                    if len(split_structs) != 1 or not split_structs[0].is_isomorphic(mol,save_order=True):
                        valid = False
                    else:
                        valid = True
                except (SiteOccupationException,TooManyElectronsException,FailedFixBondsException) as e:
                    valid = False
            else:
                valid = True
                vibdata = get_vibdata(dopt,dvib,0)
            
            if not gas_phase:
                spc = SurfaceConfiguration(atoms,slab,admol,vibdata,name,metal,facet,is_TS=False,sites_per_cell=1,
                  c_ref=c_ref,o_ref=o_ref,h_ref=h_ref,n_ref=n_ref,valid=valid)
                spc.run()
            else:
                spc = GasConfiguration(atoms,mol,vibdata,name,
                  c_ref=c_ref,o_ref=o_ref,h_ref=h_ref,n_ref=n_ref,valid=valid)
                spc.run()
                
            species_dict[ind] = spc
            
    return species_dict 
                
def get_TS(path,adsorbates_path,metal,facet,slab,sites,site_adjacency,nslab,c_ref=0.0,o_ref=0.0,h_ref=0.0,n_ref=0.0):
    """
    Generate a dictionary of GasConfiguration/SurfaceConfiguration objects corresponding to a Pynta 
    species calculation
    Args:
        path: directory of the corresponding TS in the Pynta run
        adsorbates_path: Adsorbates directory of the Pynta run
        metal: metal of the slab
        facet: facet of the slab
        slab: ase.Atoms object for the slab
        sites: list of sites correspondingn to the slab
        site_adjacency: site_adjacency corresponding to the slab
        nslab: number of atoms in the slab
        c_ref: ZPE corrected gas phase reference for carbon (CH4) [eV]
        o_ref: ZPE corrected gas phase reference for oxygen (H2O) [eV]
        h_ref: ZPE corrected gas phase reference for hydrogen(H2) [eV]
        n_ref: ZPE corrected gas phase reference for nitrogen (NH3) [eV]

    Returns:
        dictionary mapping string index to TS SurfaceConfiguration objects
    """
    allowed_structure_site_structures = generate_allowed_structure_site_structures(adsorbates_path,sites,
                                                                                   site_adjacency,nslab,max_dist=np.inf)
    valid_dict,valid_info = validate_TS_configs(path,sites,site_adjacency,nslab,irc_concern_len=8)
    
    with open(os.path.join(path,"info.json"),'r') as f:
        info = json.load(f)
        
    reactants = Molecule().from_adjacency_list(info["reactants"])
    products = Molecule().from_adjacency_list(info["products"])
    
    broken_bonds,formed_bonds = get_broken_formed_bonds(reactants,products)
    target_TS = reactants.copy(deep=True)
    
    for labels in list(broken_bonds)+list(formed_bonds):
        label1,label2 = list(labels)
        a1 = target_TS.get_labeled_atoms(label1)[0]
        a2 = target_TS.get_labeled_atoms(label2)[0]
        if target_TS.has_bond(a1,a2):
            bd = target_TS.get_bond(a1,a2)
            bd.set_order_str('R')
        else:
            bd = Bond(a1,a2,order='R')
            target_TS.add_bond(bd)
    
    target_TS.clear_labeled_atoms()
    
    keep_binding_vdW_bonds_in_reactants=False
    keep_vdW_surface_bonds_in_reactants=False
    mol = reactants
    for bd in mol.get_all_edges():
        if bd.order == 0:
            if bd.atom1.is_surface_site() or bd.atom2.is_surface_site():
                keep_binding_vdW_bonds_in_reactants = True
                m = mol.copy(deep=True)
                b = m.get_bond(m.atoms[mol.atoms.index(bd.atom1)],m.atoms[mol.atoms.index(bd.atom2)])
                m.remove_bond(b)
                out = m.split()
                if len(out) == 1: #vdW bond is not only thing connecting adsorbate to surface
                    keep_vdW_surface_bonds_in_reactants = True
    keep_binding_vdW_bonds_in_products=False
    keep_vdW_surface_bonds_in_products=False
    mol = products
    for bd in mol.get_all_edges():
        if bd.order == 0:
            if bd.atom1.is_surface_site() or bd.atom2.is_surface_site():
                keep_binding_vdW_bonds_in_products = True
                m = mol.copy(deep=True)
                b = m.get_bond(m.atoms[mol.atoms.index(bd.atom1)],m.atoms[mol.atoms.index(bd.atom2)])
                m.remove_bond(b)
                out = m.split()
                if len(out) == 1: #vdW bond is not only thing connecting adsorbate to surface
                    keep_vdW_surface_bonds_in_products = True
        
    keep_binding_vdW_bonds = keep_binding_vdW_bonds_in_reactants and keep_binding_vdW_bonds_in_products
    keep_vdW_surface_bonds = keep_vdW_surface_bonds_in_reactants and keep_vdW_surface_bonds_in_products

    ts_dict = dict()
    for k in valid_dict.keys():
        dopt = os.path.join(path,k,"opt.xyz")
        dvib = os.path.join(path,k,"vib.json_vib.json")
        vibdata = vibdata = get_vibdata(dopt,dvib,nslab)
        atoms = read(dopt)
        valid = valid_dict[k]
        vinfo = valid_info[k]
        try:
            admol,neighbor_sites,ninds = generate_TS_2D(atoms, os.path.join(path,"info.json"),  metal, facet, sites, site_adjacency, nslab, 
                                                    imag_freq_path=os.path.join(path,k,"vib.0.traj"),
                     max_dist=np.inf, allowed_structure_site_structures=allowed_structure_site_structures,
                     keep_binding_vdW_bonds=keep_binding_vdW_bonds,keep_vdW_surface_bonds=keep_vdW_surface_bonds)
            
            spc = SurfaceConfiguration(atoms,slab,admol,vibdata,os.path.split(path)[1],metal,facet,is_TS=True,sites_per_cell=1,
                      c_ref=c_ref,o_ref=o_ref,h_ref=h_ref,n_ref=n_ref,valid=valid,valid_info=vinfo,)
            
            spc.run()
            
            ts_dict[k] = spc
        except (SiteOccupationException,TooManyElectronsException,ValueError):
            spc = SurfaceConfiguration(atoms,slab,None,vibdata,os.path.split(path)[1],metal,facet,is_TS=True,sites_per_cell=1,
                      c_ref=c_ref,o_ref=o_ref,h_ref=h_ref,n_ref=n_ref,valid=valid,valid_info=vinfo,mol=target_TS)
            
            spc.run()
            
            ts_dict[k] = spc
            
    return ts_dict

def get_kinetics(path,adsorbates_path,metal,facet,slab,sites,site_adjacency,nslab,site_density,config_dict=None,
                c_ref=0.0,o_ref=0.0,h_ref=0.0,n_ref=0.0):
    """
    Generate a dictionary of GasConfiguration/SurfaceConfiguration objects corresponding to a Pynta 
    species calculation
    Args:
        path: directory of the corresponding TS in the Pynta run
        adsorbates_path: Adsorbates directory of the Pynta run
        metal: metal of the slab
        facet: facet of the slab
        slab: ase.Atoms object for the slab
        sites: list of sites correspondingn to the slab
        site_adjacency: site_adjacency corresponding to the slab
        nslab: number of atoms in the slab
        site_density: the site density of the slab in molecules/m^2
        config_dict: dictionary mapping names to lowest energy GasConfiguration/SurfaceConfiguration objects, will generate needed ones automatically if not provided
        c_ref: ZPE corrected gas phase reference for carbon (CH4) [eV]
        o_ref: ZPE corrected gas phase reference for oxygen (H2O) [eV]
        h_ref: ZPE corrected gas phase reference for hydrogen(H2) [eV]
        n_ref: ZPE corrected gas phase reference for nitrogen (NH3) [eV]

    Returns:
        dictionary mapping string index to TS SurfaceConfiguration objects
    """
    ts_dict = get_TS(path,adsorbates_path,metal,facet,slab,sites,site_adjacency,nslab,
                c_ref=c_ref,o_ref=o_ref,h_ref=h_ref,n_ref=n_ref)
    
    with open(os.path.join(path,"info.json"),'r') as f:
        info = json.load(f)
    
    if config_dict is None:
        config_dict = dict()
        for spc_name in np.unique(np.array(info["species_names"]+info["reverse_names"])):
            spcs = get_species(os.path.join(adsorbates_path,spc_name),adsorbates_path,metal,facet,slab,sites,site_adjacency,nslab,c_ref=c_ref,o_ref=o_ref,h_ref=h_ref,n_ref=n_ref)
            minspc = spcs[min({k:v for k,v in spcs.items() if v.valid},key=lambda x: spcs[x].energy)]
            config_dict[spc_name] = minspc


    if info["forward"]:
        reactant_mol = Molecule().from_adjacency_list(info["reactants"])
        product_mol = Molecule().from_adjacency_list(info["products"])
    else:
        product_mol = Molecule().from_adjacency_list(info["reactants"])
        reactant_mol = Molecule().from_adjacency_list(info["products"])
    
    reactants = [config_dict[name] for name in info["species_names"]]
    products = [config_dict[name] for name in info["reverse_names"]]

    rstr = reactant_mol.to_adjacency_list()
    reactant_mol.clear_labeled_atoms()
    pstr = product_mol.to_adjacency_list()
    product_mol.clear_labeled_atoms()

    family_comment = ""
    if "family_name" in info.keys():
        family_comment += info["family_name"]+"\n"

    family_comment += "reactants:\n" + rstr
    family_comment += "products:\n" + pstr

    kin_dict = dict()
    for k,ts in ts_dict.items():
        kin = Kinetics(reactants,ts,reactant_mol,product_mol,
                    site_density,metal=metal,facet=facet,products=products,family_comment=family_comment,valid=ts.valid)
        kin.run()
        kin_dict[k] = kin

    return kin_dict

def get_reference_energies(adsorbates_path,nslab,check_finished=False):
    """
    Extract references energies for reference species for thermochemistry calculations
    Args:
        adsorbates_path: path to the Adsorbates directory in the pynta run
        nslab: size of the slab

    Returns:
        Reference energies: c_ref,o_ref,h_ref,n_ref,
         and finished_atoms a list of the atoms that appropriate references are provided for. 
        
    """
    spc_names = [x for x in os.listdir(adsorbates_path) if os.path.isdir(os.path.join(adsorbates_path,x))]
    finished_names = []
    
    #get reference species
    if "C" in spc_names and (not check_finished or os.path.exists(os.path.join(adsorbates_path,"C","complete.sgnl"))):
        c_energies = []
        for ind in os.listdir(os.path.join(adsorbates_path,"C")):
            if not os.path.isdir(os.path.join(adsorbates_path,"C",ind)):
                continue
            cdopt = os.path.join(adsorbates_path,"C",ind,ind+".xyz")
            cdvib = os.path.join(adsorbates_path,"C",ind,"vib.json_vib.json")
            c_energy = read(cdopt).get_potential_energy() + get_vibdata(cdopt,cdvib,nslab,gas_phase=True).get_zero_point_energy()
            c_energies.append(c_energy)
        c_ref = min(c_energies)
        finished_names.append("C")
    else:
        logging.warning("No C reference (no CH4 calculation), thermochemistry with C will be wrong")
        c_ref = 0.0

    if "O" in spc_names and (not check_finished or os.path.exists(os.path.join(adsorbates_path,"O","complete.sgnl"))):
        o_energies = []
        for ind in os.listdir(os.path.join(adsorbates_path,"O")):
            if not os.path.isdir(os.path.join(adsorbates_path,"O",ind)):
                continue
            odopt = os.path.join(adsorbates_path,"O",ind,ind+".xyz")
            odvib = os.path.join(adsorbates_path,"O",ind,"vib.json_vib.json")
            o_energy = read(odopt).get_potential_energy() + get_vibdata(odopt,odvib,nslab,gas_phase=True).get_zero_point_energy()
            o_energies.append(o_energy)
        o_ref = min(o_energies)
        finished_names.append("O")
    else:
        logging.warning("No O reference (no H2O calculation), thermochemistry with O will be wrong")
        o_ref = 0.0

    if "[H][H]" in spc_names and (not check_finished or os.path.exists(os.path.join(adsorbates_path,"[H][H]","complete.sgnl"))):
        h_energies = []
        for ind in os.listdir(os.path.join(adsorbates_path,"[H][H]")):
            if not os.path.isdir(os.path.join(adsorbates_path,"[H][H]",ind)):
                continue
            hdopt = os.path.join(adsorbates_path,"[H][H]",ind,ind+".xyz")
            hdvib = os.path.join(adsorbates_path,"[H][H]",ind,"vib.json_vib.json")
            h_energy = read(hdopt).get_potential_energy() + get_vibdata(hdopt,hdvib,nslab,gas_phase=True).get_zero_point_energy()
            h_energies.append(h_energy)
        h_ref = min(h_energies)
        finished_names.append("[H][H]")
    else:
        logging.warning("No H reference (no H2 calculation), thermochemistry with H will be wrong")
        h_ref = 0.0

    if "N" in spc_names and (not check_finished or os.path.exists(os.path.join(adsorbates_path,"N","complete.sgnl"))):
        n_energies = []
        for ind in os.listdir(os.path.join(adsorbates_path,"N")):
            if not os.path.isdir(os.path.join(adsorbates_path,"N",ind)):
                continue
            ndopt = os.path.join(adsorbates_path,"N",ind,ind+".xyz")
            ndvib = os.path.join(adsorbates_path,"N",ind,"vib.json_vib.json")
            n_energy = read(ndopt).get_potential_energy() + get_vibdata(ndopt,ndvib,nslab,gas_phase=True).get_zero_point_energy()
            n_energies.append(n_energy)
        n_ref = min(n_energies)
        finished_names.append("N")
    else:
        logging.warning("No N reference (no NH3 calculation), thermochemistry with N will be wrong")
        n_ref = 0.0
    
    if "[H][H]" not in finished_names:
        finished_atoms = []
    else:
        finished_atoms = ["H"] + [x for x in finished_names if x != "[H][H]"]
    
    return c_ref,o_ref,h_ref,n_ref,finished_atoms
    
def postprocess(path,metal,facet,sites,site_adjacency,slab_path=None,check_finished=False):
    """
    Postprocess Pynta run into GasConfiguration/SurfaceConfiguration/Kinetics objects
    Args:
        path: directory of the pynta run
        metal: metal of the slab
        facet: facet of the slab
        sites: list of slab sites
        site_adjacency: slab adjacency
        slab_path: path to the slab file. Defaults to None.
        check_finished: If True the function will only postprocess directories with a complete.sgnl file in them
            During a Pynta run this allows us to avoid postprocessing unfinished calculations. Defaults to False.

    Returns:
        spc_dict,kin_dict dictionaries mapping names of species and transition states to lowest energy valid GasConfiguration/SurfaceConfiguration and 
            lowest barrier valid Kinetics objects. 
    """
    if slab_path:
        slab = read(slab_path)
    else:
        slab = read(os.path.join(path,"slab.xyz"))
        
    nslab = len(slab)

    site_density = get_site_density(slab,metal,facet)
    
    spc_names = [x for x in os.listdir(os.path.join(path,"Adsorbates")) if os.path.isdir(os.path.join(path,"Adsorbates",x))]

    c_ref,o_ref,h_ref,n_ref,finished_atoms = get_reference_energies(os.path.join(path,"Adsorbates"),nslab,check_finished=check_finished)

    spc_dict = dict()
    spc_dict_thermo = dict()
    for name in spc_names:
        if check_finished and not os.path.exists(os.path.join(path,"Adsorbates",name,"complete.sgnl")):
            continue # not ready to process this species yet
        spcs = get_species(os.path.join(path,"Adsorbates",name),os.path.join(path,"Adsorbates"),metal,facet,slab,sites,site_adjacency,
                           nslab,c_ref=c_ref,o_ref=o_ref,h_ref=h_ref,n_ref=n_ref)
        
        reference_missing = False
        if len(spcs) > 0:
            for elm in list(spcs.values())[0].mol.get_element_count().keys():
                if elm not in finished_atoms and elm != 'X': #we don't have appropriate references for this species
                    reference_missing = True 
                
        min_energy = np.inf
        mink = None
        for k,spc in spcs.items():
            if spc.valid and spc.energy < min_energy: #ZPE corrected energies
                min_energy = spc.energy
                mink = k
        if mink:
            spc = spcs[mink]
            spc_dict[name] = spc
            if not reference_missing: #if missing references don't do thermochemistry
                spc_dict_thermo[name] = spc
    
    ts_dict = dict()
    for ts in os.listdir(path):
        
        if ts[:2] != "TS":
            continue
        
        if check_finished and not os.path.exists(os.path.join(path,ts,"complete.sgnl")):
            continue # not ready to process this species yet
            
        try:
            kinetics = get_kinetics(os.path.join(path,ts),os.path.join(path,"Adsorbates"),metal,facet,slab,sites,site_adjacency,nslab,site_density,config_dict=spc_dict,
                                   c_ref=c_ref,o_ref=o_ref,h_ref=h_ref,n_ref=n_ref)
        except KeyError: #species required isn't in spc_dict yet
            continue
            
        min_barrier = np.inf
        mink = None
        for k,kinetic in kinetics.items():
            if kinetic.valid and kinetic.barrier_f < min_barrier: #ZPE corrected barriers
                min_barrier = kinetic.barrier_f
                mink = k
        if mink:
            kin = kinetics[mink]
            ts_dict[ts] = kin

    return spc_dict,ts_dict,spc_dict_thermo

def write_rmg_libraries(path,spc_dict,spc_dict_thermo,ts_dict,metal,facet):
    """
    Writes RMG thermo and kinetics libraries in the pynta directory
    Args:
        path: directory of the pynta run
        spc_dict: dictionary mapping names to lowest energy GasConfiguration/SurfaceConfiguration objects
        ts_dict: dictionary mapping TS names (TS1...etc) to Kinetics objects corresponding to the lowest energy valid TS/reactants/products
    """
    index = 0
    thermo_text = ""
    spc_dictionary_txt = "vacantX\n1 X u0 p0 c0\n\n"
    for name,spc in spc_dict_thermo.items():
        if thermo_text == "":
            if isinstance(spc,SurfaceConfiguration):
                thermo_text += spc.create_RMG_header("thermo_library","","")
            else:
                thermo_text += spc.create_RMG_header("thermo_library","","",metal,facet)
        thermo_text += spc.rmg_species_text.replace("{index}",str(index))
        index += 1
        thermo_text += "\n"
    
    for name,spc in spc_dict.items():
        spc_dictionary_txt += name + "\n"
        spc_dictionary_txt += spc.mol.to_adjacency_list()
        spc_dictionary_txt += "\n"

    with open(os.path.join(path,"thermo_library.py"),'w') as f:
        f.write(thermo_text)

    index = 0
    reaction_text = ""
    for ts,kinetics in ts_dict.items():
        if reaction_text == "":
            reaction_text += kinetics.create_RMG_header("reaction_library",lib_short_desc="",lib_long_desc="")
        reaction_text += kinetics.rmg_kinetics_text.replace("{index}",str(index)) 
        index += 1
        reaction_text += "\n"

    if os.path.exists(os.path.join(path,"reaction_library")):
        shutil.rmtree(os.path.join(path,"reaction_library"))

    os.makedirs(os.path.join(path,"reaction_library"))

    
    with open(os.path.join(path,"reaction_library","reactions.py"),'w') as f:
        f.write(reaction_text)

    with open(os.path.join(path,"reaction_library","dictionary.txt"),'w') as f:
        f.write(spc_dictionary_txt)
        
def get_energy_correction_configuration(Ncoad_energy_dict,ts_dict,config_name,coad_name,iter,Ncoad,reactant_names=None):
    """Compute the energy corrections (difference in energy between the isolated configuration and non-isolated configuraitons)
        at each coverage for adsorbate/TS with name config_name
        
        Note the TS energy correction can depend on the direction of reaction so it needs the reactant_names
    """
    if config_name not in Ncoad_energy_dict[coad_name][0].keys(): #gas phase
        return 0.0
    elif config_name not in ts_dict.keys() and config_name != coad_name: #adsorbate that is not the co-adsorbate
        return Ncoad_energy_dict[coad_name][iter][config_name][Ncoad]-Ncoad_energy_dict[coad_name][iter][coad_name][Ncoad-1]
    elif config_name not in ts_dict.keys() and config_name == coad_name: #co-adsorbate
        try:
            return Ncoad_energy_dict[coad_name][iter][config_name][Ncoad-1]/Ncoad
        except Exception as e:
            print((config_name,coad_name,iter,Ncoad))
            raise e
    else: #TS
        assert reactant_names is not None
        Ncoad_reactants = reactant_names.count(coad_name)
        if Ncoad-Ncoad_reactants <= 0:
            return 0.0
        else:
            return Ncoad_energy_dict[coad_name][iter][config_name][Ncoad-Ncoad_reactants]-Ncoad_energy_dict[coad_name][iter][coad_name][Ncoad-1]

def get_barrier_correction(Ncoad_energy_dict,ts_dict,ts_name,coad_name,iter,Ncoad,reactant_names):
    ts_correction = get_energy_correction_configuration(Ncoad_energy_dict,ts_dict,ts_name,coad_name,iter,Ncoad,reactant_names=reactant_names)
    reactant_correction = 0.0
    for rname in reactant_names:
        reactant_correction += get_energy_correction_configuration(Ncoad_energy_dict,ts_dict,rname,coad_name,iter,Ncoad)

    return ts_correction-reactant_correction

def get_reaction_energy_correction(Ncoad_energy_dict,ts_dict,ts_name,coad_name,iter,Ncoad,reactant_names,product_names):
    reactant_correction = 0.0
    for rname in reactant_names:
        reactant_correction += get_energy_correction_configuration(Ncoad_energy_dict,ts_dict,rname,coad_name,iter,Ncoad)
    product_correction = 0.0
    for pname in product_names:
        product_correction += get_energy_correction_configuration(Ncoad_energy_dict,ts_dict,pname,coad_name,iter,Ncoad)

    return product_correction - reactant_correction

def extract_covdep_data(path,pynta_path,ts_dict,metal,facet,sites,site_adjacency,nslab,coad_names):
    admol_name_path_dict = {k: os.path.join(pynta_path,k,v,"opt.xyz") for k,v in ts_dict.items()}
    admol_name_structure_dict = dict()
    ads = [x for x in os.listdir(os.path.join(pynta_path,"Adsorbates")) if x[0] != "."]
    allowed_structure_site_structures = generate_allowed_structure_site_structures(os.path.join(pynta_path,"Adsorbates"),sites,site_adjacency,nslab,max_dist=np.inf)

    for ts in ts_dict.keys():
        info_path = os.path.join(pynta_path,ts,"info.json")
        atoms = read(admol_name_path_dict[ts])
        st,_,_ = generate_TS_2D(atoms, info_path,  metal, facet, sites, site_adjacency, nslab,
                max_dist=np.inf, allowed_structure_site_structures=allowed_structure_site_structures)
        admol_name_structure_dict[ts] = st
        with open(info_path,"r") as f:
            info = json.load(f)
            for name in info["species_names"]+info["reverse_names"]:
                if name not in ads:
                    ads.append(name)

    for ad in ads:
        p = os.path.join(pynta_path,"Adsorbates",ad)
        with open(os.path.join(p,"info.json")) as f:
            info = json.load(f)
        mol = Molecule().from_adjacency_list(info["adjlist"])
        if mol.contains_surface_site():
            keep_binding_vdW_bonds=False 
            keep_vdW_surface_bonds=False
            for bd in mol.get_all_edges():
                if bd.order == 0:
                    if bd.atom1.is_surface_site() or bd.atom2.is_surface_site():
                        keep_binding_vdW_bonds = True
                        m = mol.copy(deep=True)
                        b = m.get_bond(m.atoms[mol.atoms.index(bd.atom1)],m.atoms[mol.atoms.index(bd.atom2)])
                        m.remove_bond(b)
                        out = m.split()
                        if len(out) == 1: #vdW bond is not only thing connecting adsorbate to surface
                            keep_vdW_surface_bonds = True
            ad_xyz = get_best_adsorbate_xyz(p,sites,site_adjacency,nslab,allowed_structure_site_structures,False,False)
            if ad_xyz is None:
                continue
            admol_name_path_dict[ad] = ad_xyz 
            atoms = read(ad_xyz)
            st,_,_ = generate_adsorbate_2D(atoms, sites, site_adjacency, nslab, max_dist=np.inf, allowed_structure_site_structures=allowed_structure_site_structures)
            admol_name_structure_dict[ad] = st

    ad_energy_dict = get_lowest_adsorbate_energies(os.path.join(pynta_path,"Adsorbates"))
    
    coadmol_E_dict = dict()
    coadmol_stability_dict = dict()
    for coad_name in coad_names:
        coad_path = os.path.join(pynta_path,"Adsorbates",coad_name)
        coadmol_E_dict[coad_name] = dict()
        coadmol_stability_dict[coad_name] = dict()
        Es = get_adsorbate_energies(coad_path)[0]
        for p in os.listdir(coad_path):
            if p == "info.json" or (p not in Es.keys()):
                continue
            admol_init,neighbor_sites_init,ninds_init = generate_adsorbate_2D(read(os.path.join(coad_path,p,p+"_init.xyz")),sites,site_adjacency,nslab,max_dist=np.inf,allowed_structure_site_structures=allowed_structure_site_structures)
            admol,neighbor_sites,ninds = generate_adsorbate_2D(read(os.path.join(coad_path,p,p+".xyz")),sites,site_adjacency,nslab,max_dist=np.inf,allowed_structure_site_structures=allowed_structure_site_structures)
            out_struct = split_adsorbed_structures(admol,clear_site_info=False)[0]
            out_struct_init = split_adsorbed_structures(admol_init,clear_site_info=False)[0]
            coadmol_E_dict[coad_name][out_struct] = Es[p] 
            if admol_init.is_isomorphic(admol,save_order=True):
                coadmol_stability_dict[coad_name][out_struct_init] = True
            else:
                coadmol_stability_dict[coad_name][out_struct_init] = False

    i = 0
    Ncoad_energy_dict = dict()
    Ncoad_config_dict = dict()
    tree_dict = dict()
    iter_path = os.path.join(path,"Iterations",str(i))
    
    #check if old format files
    files_are_old_format = len([file for file in os.listdir(iter_path) if file.startswith("Ncoad_energy_")][0].split("_")) == 3
    
    if files_are_old_format:
        coad_name = coad_names[0]
        Ncoad_energy_dict[coad_name] = dict()
        Ncoad_config_dict[coad_name] = dict()
        while os.path.exists(iter_path):
            nodes = read_nodes(os.path.join(iter_path,"regressor.json"))
            root = [n for n in nodes.values() if n.parent is None][0]
            tree = MultiEvalSubgraphIsomorphicDecisionTreeRegressor([adsorbate_interaction_decomposition],
                                                        nodes=nodes)
            tree_dict[i] = tree
            Ncoad_energy_dict[coad_name][i] = dict()
            Ncoad_config_dict[coad_name][i] = dict()
            files = os.listdir(iter_path)
            for file in files:
                if file.startswith("Ncoad_energy_"):
                    name = file.split("_")[-1].split(".")[0]
                    with open(os.path.join(iter_path,file),'r') as f:
                        out = {int(k): v for k,v in json.load(f).items()}
                    Ncoad_energy_dict[coad_name][i][name] = out
                elif file.startswith("Ncoad_config_"):
                    name = file.split("_")[-1].split(".")[0]
                    with open(os.path.join(iter_path,file),'r') as f:
                        out = {int(k): v for k,v in json.load(f).items()}
                    for k in out.keys():
                        if isinstance(out[k],list):
                            out[k] = [Molecule().from_adjacency_list(x,check_consistency=False) for x in out[k]]
                        else:
                            out[k] = [Molecule().from_adjacency_list(out[k],check_consistency=False)]
                    Ncoad_config_dict[coad_name][i][name] = out
                    
            i += 1
            iter_path = os.path.join(path,"Iterations",str(i))
    else:
        for coad_name in coad_names:
            Ncoad_energy_dict[coad_name] = dict()
            Ncoad_config_dict[coad_name] = dict()
            i = 0
            iter_path = os.path.join(path,"Iterations",str(i))
            while os.path.exists(iter_path):
                Ncoad_energy_dict[coad_name][i] = dict()
                Ncoad_config_dict[coad_name][i] = dict()
                i += 1
                iter_path = os.path.join(path,"Iterations",str(i))
        
        i = 0
        iter_path = os.path.join(path,"Iterations",str(i))
        while os.path.exists(iter_path):
            nodes = read_nodes(os.path.join(iter_path,"regressor.json"))
            root = [n for n in nodes.values() if n.parent is None][0]
            tree = MultiEvalSubgraphIsomorphicDecisionTreeRegressor([adsorbate_interaction_decomposition],
                                                        nodes=nodes)
            tree_dict[i] = tree
            Ncoad_energy_dict[coad_name][i] = dict()
            Ncoad_config_dict[coad_name][i] = dict()
            files = os.listdir(iter_path)
            for file in files:
                if file.startswith("Ncoad_energy_"):
                    spfile = file.split("_")
                    coad_name = spfile[-1].split(".")[0]
                    name = spfile[-2]
                    with open(os.path.join(iter_path,file),'r') as f:
                        out = {int(k): v for k,v in json.load(f).items()}
                    Ncoad_energy_dict[coad_name][i][name] = out
                elif file.startswith("Ncoad_config_"):
                    spfile = file.split("_")
                    coad_name = spfile[-1].split(".")[0]
                    name = spfile[-2]
                    with open(os.path.join(iter_path,file),'r') as f:
                        out = {int(k): v for k,v in json.load(f).items()}
                    for k in out.keys():
                        if isinstance(out[k],list):
                            out[k] = [Molecule().from_adjacency_list(x,check_consistency=False) for x in out[k]]
                        else:
                            out[k] = [Molecule().from_adjacency_list(out[k],check_consistency=False)]
                    Ncoad_config_dict[coad_name][i][name] = out
                    
            i += 1
            iter_path = os.path.join(path,"Iterations",str(i))

    ts_info = dict()
    for ts_name in ts_dict.keys():
        with open(os.path.join(pynta_path,ts_name,"info.json"),'r') as f:
            ts_info[ts_name] = json.load(f)

    max_coad_indexes = {coad_name: {i:max(Ncoad_energy_dict[coad_name][i].keys()) for i in sorted(Ncoad_energy_dict[coad_name].keys())} for coad_name in coad_names}
    
    return Ncoad_energy_dict,Ncoad_config_dict,tree_dict,admol_name_structure_dict,admol_name_path_dict,ts_info,max_coad_indexes,ad_energy_dict,coadmol_E_dict

def analyze_covdep_lowest_energy(Ncoad_config_dict,iter_configs,config_name,coad_name,metal,slab,sites,admol_name_structure_dict,admol_name_path_dict,tree_dict):
    configs_3D = []
    configs_2D = []
    sidt_Es = []
    sidt_traces = []
    partial_admol = admol_name_structure_dict[config_name]
    admol_path = admol_name_path_dict[config_name]
    partial_atoms = read(admol_path)
    for k,v in Ncoad_config_dict[coad_name][iter_configs][config_name].items():
        for x in v:
            atoms = mol_to_atoms(x,slab,sites,metal,partial_atoms=partial_atoms,partial_admol=partial_admol)
            configs_3D.append(atoms)
            Einteraction,std,tr = tree_dict[iter_configs].evaluate(x,trace=True, estimate_uncertainty=True)
            sidt_traces.append(tr)
            configs_2D.append(x)
            sidt_Es.append(Einteraction)
            
    return configs_3D,configs_2D,sidt_Es,sidt_traces

def analyze_covdep_sample_data(config_name,coad_name,Ncoad_energy_dict,path,pynta_path,
                               slab,metal,facet,sites,site_adjacency,ad_energy_dict,ts_dict,coadmol_E_dict,reactant_names=None):
    
    to_eV = 1.0/(96.48530749925793*1000.0) #converts from J/mol to eV
    configs_3D = []
    config_Es = []
    config_E_correction = []
    config_xyzs = []
    config_mols = []
    for k in range(len(Ncoad_energy_dict[coad_name])):
        if os.path.exists(os.path.join(path,"Iterations",str(k),"Samples")):
            for i in os.listdir(os.path.join(path,"Iterations",str(k),"Samples")):
                infopath = os.path.join(path,"Iterations",str(k),"Samples",i,"info.json")
                with open(infopath,"r") as f:
                    info = json.load(f)
                if os.path.split(os.path.split(os.path.split(info["xyz"])[0])[0])[1] == config_name:
                    d = os.path.join(path,"Iterations",str(k),"Samples",i)
                    datum_E,datums_stability = process_calculation(d,ad_energy_dict,slab,metal,facet,sites,site_adjacency,pynta_path,coadmol_E_dict[coad_name],max_dist=3.0,rxn_alignment_min=0.7,
                    coad_disruption_tol=1.1,out_file_name="out",init_file_name="init",vib_file_name="vib_vib",is_ad=config_name not in ts_dict.keys())
                    xyz = os.path.join(path,"Iterations",str(k),"Samples",i,"out.xyz")
                    if os.path.exists(xyz):
                        config_xyzs.append(xyz)
                        configs_3D.append(read(xyz))
                        if datum_E is not None:
                            config_Es.append(datum_E.value*to_eV)
                            config_mols.append(datum_E.mol)
                            Ncoad = len(split_adsorbed_structures(datum_E.mol)) - 1
                            if config_name not in ts_dict.keys() and config_name != coad_name: #adsorbate that is not the co-adsorbate
                                Ecorr = (datum_E.value-Ncoad_energy_dict[coad_name][len(Ncoad_energy_dict)-1][coad_name][Ncoad-1])*to_eV
                            elif config_name not in ts_dict.keys() and config_name == coad_name: #co-adsorbate
                                Ecorr = Ncoad_energy_dict[coad_name][k][config_name][Ncoad-1]/Ncoad*to_eV
                            else: #TS
                                assert reactant_names is not None
                                Ncoad_reactants = reactant_names.count(coad_name)
                                if Ncoad-Ncoad_reactants <= 0:
                                    Ecorr = 0.0
                                else:
                                    Ecorr = (datum_E.values-Ncoad_energy_dict[coad_name][len(Ncoad_energy_dict[coad_name])-1][coad_name][Ncoad-1]/Ncoad*(Ncoad-Ncoad_reactants))*to_eV

                            config_E_correction.append(Ecorr)
                        else:
                            config_Es.append(None)
                            
    return configs_3D,config_Es,config_E_correction,config_xyzs,config_mols