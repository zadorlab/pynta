#import sys
#sys.path.insert(0,'/home/tdprice/pynta_vasp/pynta/')
from pynta.main import Pynta
import os
import ast

# Reading in vasp stuff

with open('vasp_slab_params.txt') as f:
    slab_params = f.read()
with open('vasp_bulk_params.txt') as f:
    bulk_params = f.read()
with open('vasp_gas_params.txt') as f:
    gas_params = f.read()
vasp_slab_params = ast.literal_eval(slab_params)
vasp_gas_params = ast.literal_eval(gas_params)
vasp_bulk_params = ast.literal_eval(bulk_params)


#os.environ['VASP_PP_PATH']= '/home/sours/programs/vasp_PP' # '/home/ark245/programs/pseudopotentials/pseudo54
#vasp_exe = 'vasp_std'
#os.environ['VASP_COMMAND']='NTASKS=`echo $SLURM_TASKS_PER_NODE|tr \'(\' \' \'|awk \'{print $1}\'`; NNODES=`scontrol show hostnames $SLURM_JOB_NODELIST|wc -l`; NCPU=`echo " $NTASKS * $NNODES " | bc`; echo "num_cpu=" $NCPU; $(which mpirun) --map-by core --display-map --report-bindings --mca btl_openib_allow_ib true --mca btl_openib_if_include mlx5_0:1 --mca btl_openib_warn_nonexistent_if 0 --mca btl_openib_warn_no_device_params_found 0 --mca pml ob1 --mca btl openib,self,vader --mca mpi_cuda_support 0 -np $NCPU %s | tee -a op.vasp' % vasp_exe
#vasp_pp_path = os.environ['VASP_PP_PATH']
#vasp_command = os.environ['VASP_COMMAND']

pyn = Pynta(path=os.getcwd(),launchpad_path=os.getcwd(),
                fworker_path=os.getcwd(),
                queue_adapter_path=os.getcwd(),
                rxns_file=os.getcwd(),
                surface_type="fcc111",software="Vasp", vacuum=10, metal="Ag",socket=False,queue=True,njobs_queue=7,
                repeats=(3,3,4),label="tdprice",num_jobs=4,frozen_layers=2,
                software_kwargs=vasp_slab_params,
                software_kwargs_gas=vasp_gas_params,
                lattice_opt_software_kwargs=bulk_params)
pyn.generate_slab()
#pyn.execute(generate_initial_ad_guesses=True, calculate_adsorbates=True, calculate_transition_states=True, launch=True)