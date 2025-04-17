from pynta.main import Pynta

pyn = Pynta(path="/home/mjohns9/jobs/pynta/Ru0001_fischer_tropsch_3_17_25",launchpad_path="/home/mjohns9/local_launchpad.yaml",
		fworker_path="/home/mjohns9/my_fworker.yaml",
		queue_adapter_path="/home/mjohns9/local_qadapter.yaml",
		rxns_file="/home/mjohns9/jobs/pynta/Ru0001_fischer_tropsch_3_17_25/fischer_tropsch_rxns.yaml",
                surface_type="hcp0001",metal="Ru",socket=False,queue=True,njobs_queue=900,a=2.72989644,c=4.30731137,
		repeats=(3,3,4),label="Ru0001_fischer_tropsch_3_17_25_6",
software_kwargs={'kpts': (4, 4, 1), 'tprnfor': True, 'occupations': 'smearing',
                            'smearing':  'gauss', 'input_dft': 'BEEF-vdW',
                            'degauss': 0.01, 'ecutwfc': 40, 'nosym': True,
                            'conv_thr': 1e-6, 'mixing_mode': 'local-TF',
                            "pseudopotentials": {"Ru": 'Ru.pbe-spn-kjpaw_psl.1.0.0.UPF',"H": 'H.pbe-kjpaw_psl.1.0.0.UPF',"O": 'O.pbe-n-kjpaw_psl.1.0.0.UPF',"C": 'C.pbe-n-kjpaw_psl.1.0.0.UPF',"N": 
'N.pbe-n-kjpaw_psl.1.0.0.UPF',
                            },
                            "command": 'srun /opt/custom/espresso/6.6_nostress/bin/pw.x < PREFIX.pwi > PREFIX.pwo'},
                software_kwargs_gas={'kpts': 'gamma', 'tprnfor': True, 'occupations': 'smearing',
                            'smearing':  'gauss', 'input_dft': 'BEEF-vdW',
                            'degauss': 0.005, 'ecutwfc': 40, 'nosym': True,
                            'conv_thr': 1e-6, 'mixing_mode': 'local-TF',
                            'mixing_beta': 0.2, 'mixing_ndim': 10,
                            "pseudopotentials": {"Ru": 'Ru.pbe-spn-kjpaw_psl.1.0.0.UPF',"H": 'H.pbe-kjpaw_psl.1.0.0.UPF',"O": 'O.pbe-n-kjpaw_psl.1.0.0.UPF',"C": 'C.pbe-n-kjpaw_psl.1.0.0.UPF',"N": 
'N.pbe-n-kjpaw_psl.1.0.0.UPF',
                            },
                            "command": 'mpirun -np 1 /home/mjohns9/qe-7.1/bin/pw.x < PREFIX.pwi > PREFIX.pwo'},
               TS_opt_software_kwargs={"conv_thr":1e-11}, 
               lattice_opt_software_kwargs={'kpts': (25,25,25), 'ecutwfc': 40, 'degauss':0.02, 'mixing_mode': 'plain'},
               slab_path = "/home/mjohns9/jobs/pynta/Ru0001_fischer_tropsch_3_17_25/slab.xyz",
		surrogate_metal="Pt",
)

pyn.execute(generate_initial_ad_guesses=False,calculate_adsorbates=False,calculate_transition_states=True,launch=False)

