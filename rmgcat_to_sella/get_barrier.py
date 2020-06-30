import os
from ase.io import read
import re
from pathlib import Path


def find_finished_prefixes(path):
    ''' Get prefixes of sucesfully finished IRC calculations '''

    # Find all TSs for which IRC calculation (forward and reverse trajectory)
    # and subsequent minimization finished succesfully. Return prefixes.

    finished_dict = {}
    finished_list = []
    how_many_ircs = []
    all_final_png = Path(path).glob('**/irc*/*final.png')

    for final_png in all_final_png:
        # split path and return the last column with *final.png file name
        final_png = str(final_png).split('/')[-1]
        finished_list.append(final_png)

    for finished in finished_list:
        count = finished_list.count(finished)
        # generate a list of touples *png and how many of them
        how_many_ircs.append((finished, count))

    for finished, count in how_many_ircs:
        finished_dict[finished] = int(count)

    # return only those values for witch there are two entries,
    # i.e. both irc files, irc_f and irc_r (the same png filename)
    # additionally print only first two characters of keys (prefixes)

    irc_pairs = {k[:2]: v for (k, v) in finished_dict.items() if v == 2}
    prefix_list = list(irc_pairs.keys())

    return sorted(prefix_list)


def check_adsAt_type(adsAt0, adsAt1):
    ''' check if adsAt0 and adsAt1 are integers '''
    for adsAt in [adsAt0, adsAt1]:
        try:
            assert isinstance(adsAt, int)
        except AssertionError as err:
            print('{} has to be integer. It is {} now'.
                  format(adsAt, type(adsAt)))
            print(err)
            raise err


def getBarrier(facetpath, repeats, path, distCutoff, adsAt0, adsAt1):
    ''' Calculate barrier heights '''

    check_adsAt_type(adsAt0, adsAt1)

    slab = read(facetpath) * repeats
    nslab = len(slab)

    # confList = []
    # for conf in sorted(os.listdir(path), key=str):
    #     confPath = os.path.join(path, conf)
    #     if os.path.isdir(confPath):
    #         confList.append(conf)
    conf_list = find_finished_prefixes(path)

    reactants = []  # forward
    reactmp = []
    products = []  # reverse
    prodtmp = []
    TS_abs_ener = []

    # for conf in sorted(os.listdir(path), key=str):
    for conf in conf_list:
        confPath = os.path.join(path, conf)
        if os.path.isdir(confPath):
            for f in os.listdir(confPath):
                if f.endswith('irc_f.traj'):
                    ftrajPath = os.path.join(confPath, f)
                    ftrajFile = read(ftrajPath)[nslab:]
                    if ftrajFile.get_distance(adsAt0, adsAt1) < distCutoff:
                        # print(conf)
                        # print(ftrajFile.get_distance(adsAt0, adsAt1))
                        reactmp.append(f)
                    else:
                        # print(conf)
                        # print(ftrajFile.get_distance(adsAt0, adsAt1))
                        prodtmp.append(f)
                elif f.endswith('irc_r.traj'):
                    rtrajPath = os.path.join(confPath, f)
                    rtrajFile = read(rtrajPath)[nslab:]
                    if rtrajFile.get_distance(adsAt0, adsAt1) < distCutoff:
                        reactmp.append(f)
                    else:
                        prodtmp.append(f)
    # This is messy part. It works but far from beeing clean
    tmp_r = []
    good_index_react = []
    good_reactants = []
    for num, r in enumerate(reactmp):
        r_prefix = r[:2]
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
        p_prefix = p[:2]
        if p_prefix not in tmp_p:
            tmp_p.append(p_prefix)
            good_index_prod.append(num)
    for inx in good_index_prod:
        good_products.append(prodtmp[inx])
    good_products = sorted(good_products)

    set_diff_prod = set(prodtmp) - set(good_products)
    diff_prod = list(set_diff_prod)
    good_reactants = sorted(diff_prod + good_reactants)
    print(good_reactants)
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
        # print(ffile)
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
        # print(react)
        TS_ener.append("{:8.3f}".format((ts - react) * 23.06035 * 4.184))

    print('\t' + 'dE_react' + '\t' + 'dE_TS')
    for conf, rener, ts in zip(conf_list, react_ener, TS_ener):
        print(conf + '\t' + rener + ' kJ/mol' + '\t' + '\t' + ts + ' kJ/mol')


##########

# from rmgcat_to_sella.adjacency_to_3d import rmgcat_to_gratoms

def calc_TS_ener_relative_to_min(minima_path, ts_path, reactants_list):
    ''' Calculate reaction energy relatively to the most stable reactants '''
    lowest_reactant_ener = get_lowest_reactant_ener(
        minima_path, reactants_list)
    tss_ener = get_ts_ener(ts_path)
    # print(ts_ener)

    activation_barriers = []
    for ts_ener in tss_ener:
        barrier = (ts_ener - lowest_reactant_ener) * 23.06035 * 4.184
        activation_barriers.append(round(barrier, 3))
    print(activation_barriers)


def get_ts_ener(ts_path):
    ''' Get energy of all TSs'''
    ts_ener_dict = {}
    tss = get_ts_out_files(ts_path)
    for ts in tss:
        with open(ts, 'r') as f:
            data = f.readlines()
            enerLine = data[-1]
            enerVal = enerLine.split()
            ts_ener_dict[ts] = float(enerVal[3])
    return list(ts_ener_dict.values())


def get_lowest_reactant_ener(minima_path, reactants_list):
    ''' Get lowest energy of the most stable reactant. '''
    reactant_ener_dict = {}
    reactants = get_reactant_out_files(minima_path, reactants_list)
    for reactant in reactants:
        with open(reactant, 'r') as f:
            data = f.readlines()
            enerLine = data[-1]
            enerVal = enerLine.split()
            reactant_ener_dict[reactant] = float(enerVal[3])
    lowest_reactant_ener = min(reactant_ener_dict.values())
    return lowest_reactant_ener


def get_ts_out_files(ts_path):
    ''' Get TS .out files. '''
    tss = []
    ts_file_list = Path(ts_path).glob('*out')
    for ts_out_file in ts_file_list:
        tss.append(str(ts_out_file))
    return sorted(tss)


def get_reactant_out_files(minima_path, reactants_list):
    ''' Get .out files for each reactants. '''
    reactants = []
    for reactant in reactants_list:
        reactant = reactant + '_'
        outfile = '{}*out'.format(reactant)
        reactant_out_list = Path(minima_path).glob(outfile)
        for reactant_out_file in reactant_out_list:
            reactants.append(str(reactant_out_file))
        return reactants
