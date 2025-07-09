import os

# Disable MPI for ASE
os.environ['ASE_MPI'] = '0'

# Clear MPI-related environment variables
if 'OMPI_MCA_plm_slurm_args' in os.environ:
    del os.environ['OMPI_MCA_plm_slurm_args']
if 'I_MPI_HYDRA_BOOTSTRAP' in os.environ:
    del os.environ['I_MPI_HYDRA_BOOTSTRAP']
if 'I_MPI_HYDRA_BOOTSTRAP_EXEC_EXTRA_ARGS' in os.environ:
    del os.environ['I_MPI_HYDRA_BOOTSTRAP_EXEC_EXTRA_ARGS']

# Force serial execution
os.environ['OMPI_COMM_WORLD_SIZE'] = '1'
os.environ['OMPI_COMM_WORLD_RANK'] = '0'

from pynta.main import Pynta
from fairchem.core.calculate import pretrained_mlip
import torch
if torch.cuda.is_available():
    device = "cuda"
    device_name = f"CUDA (GPU: {torch.cuda.get_device_name()})"
    print(f"  ✅ CUDA available, using GPU")
else:
    device = "cpu"
    device_name = "CPU"
    print(f"  ⚠️  CUDA not available, using CPU")
predictor = pretrained_mlip.get_predict_unit("uma-s-1", device=device)

print("predictor:", predictor)
print("type:", type(predictor))


pyn = Pynta(
    path="/projects/westgroup/lekia.p/pynta/Fe",
    launchpad_path="/projects/westgroup/lekia.p/pynta/my_launchpad.yaml",
    fworker_path="/projects/westgroup/lekia.p/pynta/my_fworker.yaml",
    queue_adapter_path="/projects/westgroup/lekia.p/pynta/my_qadapter.yaml",
    rxns_file="/projects/westgroup/lekia.p/pynta/Fe/reaction.yaml",
    surface_type="bcc110",
    metal="Fe",
    socket=False,
    queue=True,
    njobs_queue=10,
    a=2.831,
    pbc=(True,True,True),
    repeats=(3, 3, 3),
    label="Nitridation_on_Fe110",
    software='FairChemCalculator',
    software_kwargs={
        # QE internal optimization settings
        'predict_unit': predictor,
        'task_name': 'oc20'
    },
    
    # Gas phase calculations
    software_kwargs_gas={
        'predict_unit': predictor,
        'task_name': 'omol'
    },
    
    # Transition state optimization settings
    TS_opt_software_kwargs={
        "forc_conv_thr": 5e-4,
        "etot_conv_thr": 1e-5, 
    },
    
    # Lattice optimization settings
    lattice_opt_software_kwargs={
        'predict_unit': predictor,
        'task_name': 'oc20'
    },
    
    # Use TBLite for frequencies (fast)
    # harm_f_software="TBLite",
    # harm_f_software_kwargs={ 
    #     "method": "GFN1-xTB",
    # },
    slab_path="/projects/westgroup/lekia.p/pynta/Fe/slab.xyz",
    surrogate_metal="Fe",
)

pyn.execute(
    calculate_adsorbates=True,
    calculate_transition_states=True,
    launch=True
)