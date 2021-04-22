from ase.io import Trajectory, write
import os

path = '/Volumes/Samsung_T3/Google_Drive/00_DYPLOM/04_POST_DOC/00_pynta_code/00_results/000_menten/00_production_runs/01/Cu_111/CH_C+H/TS_estimate_unique/00/00_Cu_111_CH_C+H.traj'


def convert_traj_to_xyzs(path):
    os.chdir(os.path.dirname(path))

    traj = Trajectory(path)

    with open('opt_steps.xyz', 'a+') as main:
        for i, step in enumerate(traj):
            i = str(i).zfill(2)
            f_name = 'step_{}.xyz'.format(i)
            write(f_name, step)
            with open(f_name, 'r') as inner:
                inner_file = inner.read()
            main.write(inner_file)
