# Map the FW-task on polaris single-queue allocation

""" Polaris Mapping """

import os
from fireworks.core.fworker import FWorker

def createCommand(node, software):
    binary = os.environ.get("EXE")
    if binary is None:
        binary = '/soft/applications/quantum_espresso/7.3.1-nvhpc23.1-libxc610/bin/pw.x'

    if software == 'Espresso':
        #command = 'mpiexec --hosts {} -n 4 --ppn 4 --depth=8 --cpu-bind depth --env OMP_NUM_THREADS=8 --env CUDA_VISIBLE_DEVICES=0,1,2,3 {} -nk 4 -in PREFIX.pwi > PREFIX.pwo'.format(node, binary)
        command = 'mpiexec --hosts {} -n 4 -ppn 4 --depth=4 --cpu-bind --cpu-bind=list:0,1:2,3  -env OMP_NUM_THREADS=1 /soft/tools/mpi_wrapper_utils/gpu_tile_compact.sh  {} -nk 4 -in PREFIX.pwi > PREFIX.pwo'.format(node, affinity, binary)
    elif software == 'PWDFT':
        #affinity = '/lus/eagle/projects/catalysis_aesp/raymundohe/bin/affinity.sh'
        affinity = ''
        #command = 'mpiexec --hosts {} -n 4 --ppn 4 --mem-bin list:0:1:2:3 --cpu-bind list:0-7:8-15:16-23:24-31 --env OMP_NUM_THREADS=1 {}  {} PREFIX.nwxi > PREFIX.nwxo'.format(node, affinity, binary)
        #command = 'mpiexec -n 4 -ppn 4 --depth=2 --cpu-bind=list:0,1:2,3  -env OMP_NUM_THREADS=1 {} {} PREFIX.nwxi > PREFIX.nwxo'.format(node, affinity, binary)
        #command = 'mpiexec --hosts {} -n 4 -ppn 4 --depth=4 --cpu-bind --cpu-bind=list:0,1:2,3  -env OMP_NUM_THREADS=1 /soft/tools/mpi_wrapper_utils/gpu_tile_compact.sh  {} {} PREFIX.nwxi > PREFIX.nwxo'.format(node, affinity, binary)
        command = 'mpiexec --hosts {} -n 12 --ppn 12 --cpu-bind list:0-7:8-15:16-23:24-31:32-39:40-47:52-59:60-67:68-75:76-83:84-91:92-99 --mem-bind list:0:0:0:0:0:0:1:1:1:1:1:1 --env OMP_NUM_THREADS=1  {} {} PREFIX.nwxi > PREFIX.nwxo'.format(node, affinity, binary)
        #command = 'mpiexec -n 4 -ppn 4 --depth=4 --cpu-bind --cpu-bind=list:0-7:8-15:16-23:24-31:32-39:40-47:52-59:60-67:68-75:76-83:84-91:92-99  -env OMP_NUM_THREADS=1 {} {} PREFIX.nwxi > PREFIX.nwxo'.format(node, affinity, binary)

    return command


def createFWorkers(num_jobs):

    nodes_list = []
    nodefile = os.environ.get("PBS_NODEFILE")

    if nodefile is None:
        for i in range(num_jobs):
            nodes_list.append(f'localhost_{i}')
    else:
        with open(nodefile) as f:
            for line in f:
                nodes_list.append(line.strip('\n'))

    #nodes_list.pop(0) # delete the first node for dispatch

    fworkers = []

    for i in range(num_jobs):
        host = nodes_list[i]
        env = {'host': host}

        fworkers.append(FWorker(env=env))

    return fworkers

