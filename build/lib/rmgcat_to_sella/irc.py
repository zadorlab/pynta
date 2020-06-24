import os

from ase.io import read, write
from rmgcat_to_sella.ts import checkSymm
from pathlib import Path


def set_up_irc(facetpath, TSdir, pytemplate_f, pytemplate_r):
    ts_path = os.path.join(facetpath, TSdir)
    rxn_path = os.path.join(ts_path, '00')
    rxn_path_list = Path(rxn_path).glob('*final.xyz')
    for rxn in rxn_path_list:
        rxn = str(rxn)
        rxn = os.path.split(rxn)[1] # only *final.xyz part of path to be extracted
        rxn = rxn.split('_')[1]
    unique_ts_index = checkSymm(ts_path)

    for i, prefix in enumerate(unique_ts_index):
        prefix = prefix[1:]
        ircDir = os.path.join(facetpath, 'IRC', prefix)
        ircpy = os.path.join(facetpath, 'IRC')
        ts_file = os.path.join(prefix + '_' + rxn + '_ts')
        TS_xyz = os.path.join(prefix, ts_file + '.xyz')
        os.makedirs(ircDir, exist_ok=True)
        ts_src_TSxyz_path = os.path.join(ts_path, prefix, prefix + '_' + rxn + '_final.xyz')
        ts_dest_path = os.path.join(ircDir, ts_file)
        try:
            write(ts_dest_path + '.xyz', read(ts_src_TSxyz_path))
            write(ts_dest_path + '.png', read(ts_src_TSxyz_path))
        except FileNotFoundError:
            pass
        '''create run job scripts'''
        with open(pytemplate_f, 'r') as f:
            template = f.read()
            job_name_f = os.path.join(ircpy, ts_file[:-3] + '_irc_f.py')
            with open(job_name_f, 'w') as f:
                f.write(template.format(
                    prefix=prefix, rxn=rxn, TS_xyz=TS_xyz))
                f.close()
        f.close()
        with open(pytemplate_r, 'r') as r:
            template = r.read()
            job_name_r = os.path.join(ircpy, ts_file[:-3] + '_irc_r.py')
            with open(job_name_r, 'w') as r:
                r.write(template.format(
                    prefix=prefix, rxn=rxn, TS_xyz=TS_xyz))
                r.close()
        r.close()
        

        # if os.path.isdir(uniqueTSdir):
        #     shutil.rmtree(uniqueTSdir)
        #     os.makedirs(uniqueTSdir, exist_ok=True)
        # else:
        #     os.makedirs(uniqueTSdir, exist_ok=True)
    # get ts geoms and create png

    # for conf in os.listdir(ts_path):
    #     conf_path = os.path.join(ts_path, conf)
    #     if os.path.isdir(conf_path):
    #         for TSxyz in os.listdir(conf_path):
    #             if TSxyz.endswith('_final.xyz'):
    #                 print(TSxyz)
                    # rxn = TSxyz.split('_')[1]
                    # prefix = TSxyz[:2]
                    # fname = os.path.join(prefix + '_' + rxn)
                    # ircDir = os.path.join(facetpath, 'IRC', prefix)
                    # ircpy = os.path.join(facetpath, 'IRC')
                    # os.makedirs(ircDir, exist_ok=True)
                    # ts_src_TSxyz_path = os.path.join(conf_path, TSxyz)
                    # ts_dest_path = os.path.join(ircDir, fname)
                    # write(ts_dest_path + '.xyz', read(ts_src_TSxyz_path))
                    # write(ts_dest_path + '.png', read(ts_src_TSxyz_path))
                    # TS_xyz = os.path.join(prefix, fname + '.xyz')
                    # # create run job scripts
                    # with open(pytemplate_f, 'r') as f:
                    #     template = f.read()
                    #     job_name_f = os.path.join(ircpy, fname + '_irc_f.py')
                    #     with open(job_name_f, 'w') as f:
                    #         f.write(template.format(
                    #             prefix=prefix, rxn=rxn, TS_xyz=TS_xyz))
                    #         f.close()
                    # f.close()
                    # with open(pytemplate_r, 'r') as r:
                    #     template = r.read()
                    #     job_name_r = os.path.join(ircpy, fname + '_irc_r.py')
                    #     with open(job_name_r, 'w') as r:
                    #         r.write(template.format(
                    #             prefix=prefix, rxn=rxn, TS_xyz=TS_xyz))
                    #         r.close()
                    # r.close()


def set_up_irc_choosen(path, list_struc, pytemplate_f, pytemplate_r):
    ts_path_list = []
    for struc in list_struc:
        struc_path = os.path.join(path, struc)
        ts_path_list.append(struc_path)

    for tsDir in ts_path_list:
        for traj in os.listdir(tsDir):
            if traj.endswith('traj'):
                trajPath = os.path.join(tsDir, traj)
                xyzPath = os.path.join(tsDir, traj[:-5] + '_final.xyz')
                write(xyzPath, read(trajPath))
            if traj.endswith('final.xyz'):
                fname = os.path.join(traj[:-10] + '_ts')
                rxn = fname[4:][:4]
                prefix = traj[:3]
                facetpath = os.path.split(path)
                ircDir = os.path.join(facetpath[0], 'IRC', prefix)
                ircpy = os.path.join(facetpath[0], 'IRC')
                os.makedirs(ircDir, exist_ok=True)
                ts_src_TSxyz_path = os.path.join(tsDir, traj)
                ts_dest_path = os.path.join(ircDir, fname)
                write(ts_dest_path + '.xyz', read(ts_src_TSxyz_path))
                write(ts_dest_path + '.png', read(ts_src_TSxyz_path))
                TS_xyz = os.path.join(prefix, fname + '.xyz')
                # create run job scripts
                with open(pytemplate_f, 'r') as f:
                    template = f.read()
                    job_name_f = os.path.join(ircpy, fname + '_irc_f.py')
                    with open(job_name_f, 'w') as f:
                        f.write(template.format(
                            prefix=prefix, rxn=rxn, TS_xyz=TS_xyz))
                        f.close()
                f.close()
                with open(pytemplate_r, 'r') as r:
                    template = r.read()
                    job_name_r = os.path.join(ircpy, fname + '_irc_r.py')
                    with open(job_name_r, 'w') as r:
                        r.write(template.format(
                            prefix=prefix, rxn=rxn, TS_xyz=TS_xyz))
                        r.close()
                r.close()


def set_up_irc_ts_estimate(facetpath, TSdir, pytemplate_f, pytemplate_r):
    # get ts geoms and create png
    ts_path = os.path.join(facetpath, TSdir + '_unique/')
    for conf in os.listdir(ts_path):
        conf_path = os.path.join(ts_path, conf)
        if os.path.isdir(conf_path):
            for TStraj in os.listdir(conf_path):
                if TStraj.endswith('.traj'):
                    srcTStraj = os.path.join(conf_path, TStraj)
                    destTStraj = os.path.join(
                        conf_path, TStraj[:-5] + '_final.xyz')
                    write(destTStraj, read(srcTStraj))
                if TStraj.endswith('.xyz'):
                    fname = os.path.join(TStraj[:8] + '_ts')
                    # fname_irc = TSxyz[:6]
                    rxn = fname[4:][:4]
                    prefix = TStraj[:3]
                    ircDir = os.path.join(facetpath, 'IRC_' + TSdir, prefix)
                    ircpy = os.path.join(facetpath, 'IRC_' + TSdir)
                    os.makedirs(ircDir, exist_ok=True)
                    ts_src_TSxyz_path = os.path.join(conf_path, TStraj)
                    ts_dest_path = os.path.join(ircDir, fname)
                    write(ts_dest_path + '.xyz', read(ts_src_TSxyz_path))
                    write(ts_dest_path + '.png', read(ts_src_TSxyz_path))
                    TS_xyz = os.path.join(prefix, fname + '.xyz')
                    # create run job scripts
                    with open(pytemplate_f, 'r') as f:
                        template = f.read()
                        job_name_f = os.path.join(ircpy, fname + '_irc_f.py')
                        with open(job_name_f, 'w') as f:
                            f.write(template.format(
                                prefix=prefix, rxn=rxn, TS_xyz=TS_xyz))
                            f.close()
                    f.close()
                    with open(pytemplate_r, 'r') as r:
                        template = r.read()
                        job_name_r = os.path.join(ircpy, fname + '_irc_r.py')
                        with open(job_name_r, 'w') as r:
                            r.write(template.format(
                                prefix=prefix, rxn=rxn, TS_xyz=TS_xyz))
                            r.close()
                    r.close()


def prepare_opt(strucPath, irc, traj, pytemplate):
    optDir = os.path.join(strucPath, irc + '_opt')
    os.makedirs(optDir, exist_ok=True)
    trajPath = os.path.join(strucPath, traj)
    initXYZ = os.path.join(optDir, traj[:-5])
    write(initXYZ + '.xyz', read(trajPath))
    write(initXYZ + '_initial.png', read(trajPath))
    geom = os.path.join(traj[:-4] + 'xyz')
    create_job_files(pytemplate, optDir, traj, geom)


def create_job_files(pytemplate, optDir, traj, geom):
    with open(pytemplate, 'r') as f:
        pytemplate = f.read()
        fname = os.path.join(optDir, traj[:-5] + '_opt.py')
        prefix = traj[:2]
        # rxn = traj[:-11][3:]
        rxn = traj.split('_')[1]
        with open(fname, 'w') as f:
            f.write(pytemplate.format(geom=geom, rxn=rxn, prefix=prefix))
        f.close()
    f.close()


def optAfterIRC(path, pytemplate):
    ''' Function to set up optimization to minimas after IRC calculation
        
        path - path to main IRC directory, e.g. Cu_111/IRC
        pytemplate - slurm template for this calculations

        The function checks if *traj files containing IRC trajectories exists and have some content. If so, a minimum optimization will be set up. Otherwise, the given geometry is skipped.
    '''
    for struc in sorted(os.listdir(path), key=str):
        strucPath = os.path.join(path, struc)
        if os.path.isdir(strucPath):
            for traj in os.listdir(strucPath):
                if traj.endswith('irc_f.traj'):
                    trajPath = os.path.join(strucPath, traj)
                    if os.stat(trajPath).st_size == 0:
                        pass
                    else:
                        try:
                            prepare_opt(strucPath, 'irc_f', traj, pytemplate)
                        except FileNotFoundError:
                            pass
                elif traj.endswith('irc_r.traj'):
                    trajPath = os.path.join(strucPath, traj)
                    if os.stat(trajPath).st_size == 0:
                        pass
                    else:
                        try:
                            prepare_opt(strucPath, 'irc_r', traj, pytemplate)
                        except FileNotFoundError:
                            pass
            # Error handling to be developed