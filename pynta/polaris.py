# Map the FW-task on polaris single-queue allocation

""" Polaris Mapping """

import os


def getNodes():
    nodes_list = []
    nodefile = os.environ.get("PBS_NODEFILE")

    if nodefile is None:
        nodes_list.append('localhost')
        return nodes_list

    with open(nodefile) as f:
        for line in f:
            nodes_list.append(line.strip('\n'))

    return nodes_list


def getBin():
    binary = os.environ.get('QE_EXE')
    if binary is None:
        binary = '/lus/eagle/projects/catalysis_aesp/abagusetty/qe-7.1/build_cuda_nompiaware/bin/pw.x'

    return binary


def readMapTxt(fname):
    maplist = []
    try:
        with open(fname, "r") as file:
            for line in file:
                num = int(line.strip())
                maplist.append(num)
    except FileNotFoundError:
        maplist = []

    return maplist


def writeMapTxt(fname, maplist):
    with open(fname, "w") as file:
        for num in maplist:
            file.write(str(num) + "\n")


def compare_returnmax(list1, list2):
    if list1 and list2:
        last1 = list1[-1]
        last2 = list2[-1]

        if last1 == last2:
            return last1
        else:
            return max(last1, last2)
    else:
        return None


class MapTaskToNodes():
    count = 0  # Class variable to count the number of instances created.
    listCount = []
    dir_root = os.getcwd()
    fname = f'{dir_root}/map.txt'
    nodelist = []

    flag = True
    flag_init = False

    def init_data(self):
        if MapTaskToNodes.flag_init == False:
            nodelist = getNodes()
            nodelist.pop(0)
            MapTaskToNodes.nnodes = len(nodelist)
            MapTaskToNodes.nodelist = nodelist

            MapTaskToNodes.flag_init = True

    def _get_index(self):
        return self.icount % MapTaskToNodes.nnodes

    def getCommand(self):
        affinity = ['list:24-31:56-63', 'list:16-23:48-55',
                    'list:8-15:40-47', 'list:0-7:32-39']

        inode =  self._get_index()

        binary = getBin()

        command = 'mpiexec --hosts {} -n 4 --ppn 4 -d 8 --cpu-bind depth --env OMP_NUM_THREADS=8 --env CUDA_VISIBLE_DEVICES=0,1,2,3  {} -nk 4 -in PREFIX.pwi > PREFIX.pwo'.format(MapTaskToNodes.nodelist[inode], binary)

        return command


    def __init__(self, _method=None):
        self.method = _method
        self.icount = MapTaskToNodes.count
        self.maplist = []


        self.init_data()

        if MapTaskToNodes.count == 0 and MapTaskToNodes.flag == True:
            if os.path.exists(MapTaskToNodes.fname):
                os.remove(MapTaskToNodes.fname)
                MapTaskToNodes.flag = False

        if self.method != 'MDMin':
            MapTaskToNodes.count += 1

        if MapTaskToNodes.count > 2:
            self.maplist = readMapTxt(MapTaskToNodes.fname)

        writeMapTxt(MapTaskToNodes.fname, MapTaskToNodes.listCount)

        if compare_returnmax(self.maplist, MapTaskToNodes.listCount) == None:
            MapTaskToNodes.listCount.append(MapTaskToNodes.count)
        else:
            MapTaskToNodes.listCount.append(self.icount)
