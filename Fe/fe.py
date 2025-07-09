from pynta.main import Pynta

pyn = Pynta(
    path="/projects/westgroup/lekia.p/pynta/Fe",
    launchpad_path="/projects/westgroup/lekia.p/pynta/my_launchpad.yaml",
    fworker_path="/projects/westgroup/lekia.p/pynta/my_fworker.yaml",
    queue_adapter_path="/projects/westgroup/lekia.p/pynta/my_qadapter.yaml",
    rxns_file="/projects/westgroup/lekia.p/pynta/Fe/reaction.yaml",
    # rxns_file="/projects/westgroup/lekia.p/pynta/fischer_tropsch_rxns.yaml",
    surface_type="bcc110",
    metal="Fe",
    socket=False,
    queue=True,
    njobs_queue=900,
    a=2.831,
    repeats=(3, 3, 4),
    label="Nitridation_on_Fe110",
    software='Espresso',
    software_kwargs={
        'kpts': (5, 5, 1),
        'tprnfor': True,
        'occupations': 'smearing',
        'smearing': 'methfessel-paxton',
        'input_dft': 'PBE',
        'vdw_corr': 'grimme-d3', # newly added
        'diagonalization': 'david', # newly added
        'diago_david_ndim': 4, # newly added
        'nspin': 2,
        'starting_magnetization(1)': 0.7,
        'degauss': 0.01,
        'ecutwfc': 40,
        'ecutrho': 400.0, # newly added
        'nosym': True,
        'conv_thr': 1e-6,
        'mixing_mode': 'local-TF',
        'mixing_beta': 0.2,
        'mixing_ndim': 10,
        'electron_maxstep': 200, # newly added
        'pseudo_dir': "/home/lekia.p/SSSP_1.3.0_PBE_efficiency",
        "pseudopotentials": {
            "Fe": 'Fe.pbe-spn-kjpaw_psl.0.2.1.UPF',
            "H":  'H.pbe-rrkjus_psl.1.0.0.UPF',
            "O":  'O.pbe-n-kjpaw_psl.0.1.UPF',
            "N":  'N.pbe-n-radius_5.UPF',
            "C": 'C.pbe-n-kjpaw_psl.1.0.0.UPF',
            "Pt": 'pt_pbe_v1.4.uspp.F.UPF',
        },
        "command": '/projects/westgroup/espresso/qe-7.0/bin/pw.x < PREFIX.pwi > PREFIX.pwo',
        # "command": 'pw.x < PREFIX.pwi > PREFIX.pwo'
    },
    software_kwargs_gas={
        'kpts': 'gamma',
        'tprnfor': True,
        'occupations': 'smearing',
        'smearing': 'methfessel-paxton',
        'input_dft': 'PBE',
        'degauss': 0.005,
        'ecutwfc': 40,
        'nosym': True,
        'conv_thr': 1e-8,
        'mixing_mode': 'local-TF',
        'mixing_beta': 0.2,
        'mixing_ndim': 10,
        'electron_maxstep': 200,
        'pseudo_dir': "/home/lekia.p/SSSP_1.3.0_PBE_efficiency",
        "pseudopotentials": {
            "Fe": 'Fe.pbe-spn-kjpaw_psl.0.2.1.UPF',
            "H":  'H.pbe-rrkjus_psl.1.0.0.UPF',
            "O":  'O.pbe-n-kjpaw_psl.0.1.UPF',
            "N":  'N.pbe-n-radius_5.UPF',
            "C": 'C.pbe-n-kjpaw_psl.1.0.0.UPF',
            "Pt": 'pt_pbe_v1.4.uspp.F.UPF',
        },
        "command": '/projects/westgroup/espresso/qe-7.0/bin/pw.x < PREFIX.pwi > PREFIX.pwo',
        # "command": 'pw.x < PREFIX.pwi > PREFIX.pwo'
    },
    TS_opt_software_kwargs={"conv_thr": 1e-11},
    lattice_opt_software_kwargs={
        'kpts': (25,25,25),
        'ecutwfc': 30,
        'ecutrho': 400.0, # newly added
        'degauss': 0.02,
        'mixing_mode': 'plain',
        'mixing_beta': 0.05,            # Very conservative mixing
        'conv_thr': 1e-4,               # Relaxed convergence
    },
    harm_f_software="TBLite",
    harm_f_software_kwargs={"method": "GFN2-xTB"},
    # slab_path="/projects/westgroup/lekia.p/pynta/Fe/fe_slab.xyz",
    surrogate_metal="Pt",
)

# pyn.generate_slab()

pyn.execute(
    calculate_adsorbates=True,
    calculate_transition_states=False,
    launch=True
)
