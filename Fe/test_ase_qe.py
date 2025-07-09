from ase import Atoms
from ase.calculators.espresso import Espresso

atoms = Atoms('Fe8', 
              positions=[
              (0.0000000000,2.0018193000, 8.0000000000),
              (2.8310000000,2.0018193000,8.0000000000),
              (1.4155000000,4.0036386000,8.0000000000),
              (4.2465000000,4.0036386000,8.0000000000),
              (0.0000000000,0.0000000000,10.0018193000),
              (2.8310000000,0.0000000000,10.0018193000),
              (1.4155000000,2.0018193000,10.0018193000),
              (4.2465000000,2.0018193000,10.0018193000)
              ],  # use your final coords
              cell=[
              (5.66200000000000, 0.00000000000000, 0.00000000000000),
              (2.83100000000000, 4.00363859507823, 0.00000000000000),
              (0.00000000000000, 0.00000000000000, 18.00181929753911)
              ], pbc=[True,True,False])

input_data = {
        'calculation':'scf',
        'tprnfor':True,
        'tstress':True,
        'verbosity':'medium',
        'nspin':2,
        'starting_magnetization(1)':0.6,
        'ecutwfc':50,
        'ecutrho':600.0,
        'occupations':'smearing',
        'smearing':'marzari-vanderbilt',
        'degauss':0.03,
        'conv_thr':1e-5,
        'mixing_mode':'plain',
        'mixing_beta':0.05,
        'diagonalization':'david'
}
calc = Espresso(
    pseudo_dir="/home/lekia.p/SSSP_1.3.0_PBE_efficiency",
    pseudopotentials={"Fe":"Fe.pbe-spn-kjpaw_psl.0.2.1.UPF"},
    input_data=input_data,
    kpts=(3,3,1),
    input_dft="PBE",
    command="pw.x -in espresso.pwi > espresso.pwo 2>&1",
)

atoms.set_calculator(calc)

try:
    f = atoms.get_forces()
except Exception:
    # now espresso.pwo contains both stdout and stderr
    with open("espresso.pwo") as out:
        print(out.read()[-2000:])   # last 2000 chars of the log
    raise
