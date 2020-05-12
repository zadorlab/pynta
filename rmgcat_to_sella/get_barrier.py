import os
from ase.io import read
import re

# barrier is calculate in referece to the forward or reverse trajectory as
# after the irc calculations. Note, that some forward or reverse
# strucutres are not converge to the local minima - further optimization
# required.
def getBarrierAfterIRC(path, distCutoff):
    confList = []
    for conf in sorted(os.listdir(path), key=str):
        confPath = os.path.join(path, conf)
        if os.path.isdir(confPath):
            confList.append(conf)
    # print(confList)

    reactants = []  # forward
    reactmp = []
    products = []  # reverse
    prodtmp = []
    TS_abs_ener = []

    for conf in sorted(os.listdir(path), key=str):
        confPath = os.path.join(path, conf)
        if os.path.isdir(confPath):
            for f in os.listdir(confPath):
                if f.endswith('irc_f.traj'):
                    ftrajPath = os.path.join(confPath, f)
                    ftrajFile = read(ftrajPath)[36:]
                    if ftrajFile.get_distance(0, 1) < distCutoff:
                        reactmp.append(f)
                    elif ftrajFile.get_distance(0, 1) > distCutoff:
                        prodtmp.append(f)
                elif f.endswith('irc_r.traj'):
                    ftrajPath = os.path.join(confPath, f)
                    ftrajFile = read(ftrajPath)[36:]
                    # for el in reactmp:
                    #     if f_prefix in el:
                    #         print(f)
                    if ftrajFile.get_distance(0, 1) < distCutoff:
                        reactmp.append(f)
                    elif ftrajFile.get_distance(0, 1) > distCutoff:
                        prodtmp.append(f)
    # print(reactmp)
    # print(prodtmp)

    # This is a messy part. It works but far from beeing clean
    tmp_r = []
    good_index_react = []
    good_reactants = []
    for num, r in enumerate(reactmp):
        r_prefix = r[:3]
        if r_prefix not in tmp_r:
            tmp_r.append(r_prefix)
            good_index_react.append(num)
    for inx in good_index_react:
        good_reactants.append(reactmp[inx])

    set_diff_react = set(reactmp) - set(good_reactants)
    diff_react = list(set_diff_react)
    prodtmp = diff_react + prodtmp

    tmp_p = []
    good_index_prod = []
    good_products = []
    for num, p in enumerate(prodtmp):
        p_prefix = p[:3]
        if p_prefix not in tmp_p:
            tmp_p.append(p_prefix)
            good_index_prod.append(num)
    for inx in good_index_prod:
        good_products.append(prodtmp[inx])
    good_products = sorted(good_products)

    set_diff_prod = set(prodtmp) - set(good_products)
    diff_prod = list(set_diff_prod)
    good_reactants = sorted(diff_prod + good_reactants)
    # print(good_reactants)
    # print(good_products)
    #

    for geom in good_reactants:
        ffile = os.path.join(path, geom[:-5] + '.out')
        with open(ffile, 'r') as f:
            enerVal = None
            reactant_ener = []
            for line in f:
                if 'IRC' in line:
                    line = line.split()
                    # get the energy of irc_f structure (last matching
                    # line)
                    enerVal = float(line[3])
                    reactant_ener.append(line[3])
            if enerVal is not None:
                reactants.append(enerVal)
        f.close()
        print(reactant_ener)
        # TS_abs_ener.append(float(reactant_ener[0]))

    for geom in good_products:
        ffile = os.path.join(path, geom[:-5] + '.out')
        with open(ffile, 'r') as f:
            enerVal = None
            for line in f:
                if 'IRC' in line:
                    line = line.split()
                    # get the energy of irc_f structure (last matching
                    # line)
                    enerVal = float(line[3])
            if enerVal is not None:
                products.append(enerVal)
        f.close()

    react_ener = []
    for react, prod in zip(reactants, products):
        react_ener.append("{:8.3f}".format((prod - react) * 23.06035 * 4.184))

    TS_ener = []
    for react, ts in zip(reactants, TS_abs_ener):
        TS_ener.append("{:8.3f}".format((ts - react) * 23.06035 * 4.184))

    print('CO* + O* --> CO2* + *' + '\t' + 'CO* + O* --> TS')
    for rener, ts in zip(react_ener, TS_ener):
        print(rener + ' kJ/mol' + '\t' + '\t' + ts + ' kJ/mol')


def getBarrier(path, distCutoff):
    confList = []
    for conf in sorted(os.listdir(path), key=str):
        confPath = os.path.join(path, conf)
        if os.path.isdir(confPath):
            confList.append(conf)
    # print(confList)

    reactants = []  # forward
    reactmp = []
    products = []  # reverse
    prodtmp = []
    TS_abs_ener = []

    for conf in sorted(os.listdir(path), key=str):
        confPath = os.path.join(path, conf)
        if os.path.isdir(confPath):
            for f in os.listdir(confPath):
                if f.endswith('irc_f.traj'):
                    ftrajPath = os.path.join(confPath, f)
                    ftrajFile = read(ftrajPath)[36:]
                    if ftrajFile.get_distance(0, 1) < distCutoff:
                        reactmp.append(f)
                    elif ftrajFile.get_distance(0, 1) > distCutoff:
                        prodtmp.append(f)
                elif f.endswith('irc_r.traj'):
                    ftrajPath = os.path.join(confPath, f)
                    ftrajFile = read(ftrajPath)[36:]
                    # for el in reactmp:
                    #     if f_prefix in el:
                    #         print(f)
                    if ftrajFile.get_distance(0, 1) < distCutoff:
                        reactmp.append(f)
                    elif ftrajFile.get_distance(0, 1) > distCutoff:
                        prodtmp.append(f)
    # print(reactmp)
    # print(prodtmp)

    # This is messy part. It works but far from beeing clean
    tmp_r = []
    good_index_react = []
    good_reactants = []
    for num, r in enumerate(reactmp):
        r_prefix = r[:3]
        if r_prefix not in tmp_r:
            tmp_r.append(r_prefix)
            good_index_react.append(num)
    for inx in good_index_react:
        good_reactants.append(reactmp[inx])

    set_diff_react = set(reactmp) - set(good_reactants)
    diff_react = list(set_diff_react)
    prodtmp = diff_react + prodtmp

    tmp_p = []
    good_index_prod = []
    good_products = []
    for num, p in enumerate(prodtmp):
        p_prefix = p[:3]
        if p_prefix not in tmp_p:
            tmp_p.append(p_prefix)
            good_index_prod.append(num)
    for inx in good_index_prod:
        good_products.append(prodtmp[inx])
    good_products = sorted(good_products)

    set_diff_prod = set(prodtmp) - set(good_products)
    diff_prod = list(set_diff_prod)
    good_reactants = sorted(diff_prod + good_reactants)
    # print(good_reactants)
    # print(good_products)
    ###
    # calculating how many elemet
    '''
    find match between '_ _' --> will give rxn name
    '''
    rxn = re.search('_(.+?)_', good_reactants[0])
    if rxn:
        rxn = rxn.group(1)
    ####

    underscore_len = good_reactants[0].count('_') - 1
    howManyRemove = len(rxn) + len(good_reactants[:2]) + underscore_len


    for geom in good_reactants:
        whichIRC = os.path.join(geom[howManyRemove:][:-5] + '_opt')
        prefix = geom[:2]
        ffile = os.path.join(path, prefix, whichIRC, geom[:-7] + '_relax.out')
        with open(ffile, 'r') as f:
            enerVal = None
            # reactant_ener = []
            for line in f:
                if 'Sella' in line:
                    line = line.split()
                    # get the energy of irc_f structure (last matching
                    # line)
                    enerVal = float(line[3])
                    # reactant_ener.append(line[3])
            if enerVal is not None:
                reactants.append(enerVal)
        f.close()
        tsfile = os.path.join(path, geom[:-5] + '.out')
        with open(tsfile, 'r') as ts:
            reactant_ener = []
            for line in ts:
                if 'IRC' in line:
                    line = line.split()
                    reactant_ener.append(line[3])
        ts.close()
        TS_abs_ener.append(float(reactant_ener[0]))

    for geom in good_products:
        whichIRC = os.path.join(geom[howManyRemove:][:-5] + '_opt')
        prefix = geom[:2]
        ffile = os.path.join(path, prefix, whichIRC, geom[:-7] + '_relax.out')
        with open(ffile, 'r') as f:
            enerVal = None
            for line in f:
                if 'Sella' in line:
                    line = line.split()
                    # get the energy of irc_f structure (last matching
                    # line)
                    enerVal = float(line[3])
            if enerVal is not None:
                products.append(enerVal)
        f.close()

    react_ener = []
    for react, prod in zip(reactants, products):
        react_ener.append("{:8.3f}".format((prod - react) * 23.06035 * 4.184))

    TS_ener = []
    for react, ts in zip(reactants, TS_abs_ener):
        TS_ener.append("{:8.3f}".format((ts - react) * 23.06035 * 4.184))

    print('dE_react' + '\t' + 'dE_TS')
    for rener, ts in zip(react_ener, TS_ener):
        print(rener + ' kJ/mol' + '\t' + '\t' + ts + ' kJ/mol')

