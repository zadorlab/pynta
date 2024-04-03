# Map the FW-task on polaris single-queue allocation

""" Polaris Mapping """

import os
from fireworks.core.fworker import FWorker

def createCommand(node, software):
    binary = os.environ.get("EXE")
    if binary is None:
        binary = '/lus/eagle/projects/catalysis_aesp/raymundohe/soft/qe-7.1/build_cuda_nogpuawarempi_openmp/bin/pw.x'

    if software == 'Espresso':
        command = 'mpiexec --hosts {} -n 4 --ppn 4 --depth=8 --cpu-bind depth --env OMP_NUM_THREADS=8 --env CUDA_VISIBLE_DEVICES=0,1,2,3 {} -in PREFIX.pwi > PREFIX.pwo'.format(node, binary)
    elif software == 'PWDFT':
        command = 'mpiexec --hosts {} -n 4 --ppn 4 --cpu-bind depth --env OMP_NUM_THREADS=1 --env CUDA_VISIBLE_DEVICES=0,1,2,3 {} PREFIX.nwxi > PREFIX.nwxo'.format(node, binary)
#add NWChem?
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

    fworkers = []

    for i in range(num_jobs):
        host = nodes_list[i]
        print("host", host)
        env = {'host': host}

        fworkers.append(FWorker(env=env))

    return fworkers

