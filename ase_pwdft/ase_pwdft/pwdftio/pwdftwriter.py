import os
import numpy as np
from copy import deepcopy
try:
    from ase.dft.kpoints import KPoints, kpts2kpts
except ImportError:
    # Fallback for older ASE versions
    KPoints = None
    kpts2kpts = None

_special_kws = ['center', 'autosym', 'autoz', 'theory', 'basis', 'xc', 'task',
                'set', 'symmetry', 'label', 'geompar', 'basispar', 'kpts',
                'bandpath', 'restart_kw']

_system_type = {1: 'polymer', 2: 'surface', 3: 'crystal'}


def _get_geom(atoms, **params):
    geom_header = ['geometry units angstrom']
    for geomkw in ['center', 'autosym', 'autoz']:
        geom_header.append(geomkw if params.get(geomkw) else 'no' + geomkw)
    if 'geompar' in params:
        geom_header.append(params['geompar'])
    geom = [' '.join(geom_header)]

    outpos = atoms.get_positions()
    pbc = atoms.pbc
    if np.any(pbc):
        scpos = atoms.get_scaled_positions()
        for i, pbci in enumerate(pbc):
            if pbci:
                outpos[:, i] = scpos[:, i]
        npbc = pbc.sum()
        cellpars = atoms.cell.cellpar()
        
        if npbc == 2:  # Check if npbc is 2
            geom.append('  system {}'.format(_system_type[npbc]))
            geom.append('    lat_a {:20.16e}'.format(cellpars[0]))
            geom.append('    lat_b {:20.16e}'.format(cellpars[1]))
            geom.append('    lat_c {:20.16e}'.format(cellpars[2]))
            geom.append('    alpha {:20.16e}'.format(cellpars[3]))
            geom.append('    beta {:20.16e}'.format(cellpars[4]))
            geom.append('    gamma {:20.16e}'.format(cellpars[5]))
        
        elif npbc == 3:  # Check if npbc is 3
            geom.append('  system {} units angstrom'.format(_system_type[3]))  # if it is surface, follow the same input as crystal
            geom.append('    lattice_vectors')
            for row in atoms.cell:
                geom.append('      {:20.16e} {:20.16e} {:20.16e}'.format(*row))
        
        else:
            if pbc[0]:
                geom.append('    lat_a {:20.16e}'.format(cellpars[0]))
            if pbc[1]:
                geom.append('    lat_b {:20.16e}'.format(cellpars[1]))
            if pbc[2]:
                geom.append('    lat_c {:20.16e}'.format(cellpars[2]))
            if pbc[1] and pbc[2]:
                geom.append('    alpha {:20.16e}'.format(cellpars[3]))
            if pbc[0] and pbc[2]:
                geom.append('    beta {:20.16e}'.format(cellpars[4]))
            if pbc[1] and pbc[0]:
                geom.append('    gamma {:20.16e}'.format(cellpars[5]))
        
        geom.append('  end')

    for i, atom in enumerate(atoms):
        geom.append('  {:<2} {:20.16e} {:20.16e} {:20.16e}'
                    ''.format(atom.symbol, *outpos[i]))
    symm = params.get('symmetry')
    if symm is not None:
        geom.append('  symmetry {}'.format(symm))
    geom.append('end')
    return geom


def _get_basis(theory, **params):
    return []


"""
    if 'basis' not in params:
        if theory in ['pspw', 'band', 'paw']:
            return []
    basis_in = params.get('basis', '3-21G')
    if 'basispar' in params:
        header = 'basis {} noprint'.format(params['basispar'])
    else:
        header = 'basis noprint'
    basis_out = [header]
    if isinstance(basis_in, str):
        basis_out.append('   * library {}'.format(basis_in))
    else:
        for symbol, ibasis in basis_in.items():
            basis_out.append('{:>4} library {}'.format(symbol, ibasis))
    basis_out.append('end')
    return basis_out
"""

_special_keypairs = [('nwpw', 'simulation_cell'),
                     ('nwpw', 'carr-parinello'),
                     ('nwpw', 'brillouin_zone'),
                     ('tddft', 'grad'),
                     ]


def _format_brillouin_zone(array, name=None):
    out = ['  brillouin_zone']
    if name is not None:
        out += ['    zone_name {}'.format(name)]
    template = '    kvector' + ' {:20.16e}' * array.shape[1]
    for row in array:
        out.append(template.format(*row))
    out.append('  end')
    return out


def _get_bandpath(bp):
    if bp is None:
        return []
    out = ['nwpw']
    out += _format_brillouin_zone(bp.kpts, name=bp.path)
    out += ['  zone_structure_name {}'.format(bp.path),
            'end',
            'task band structure']
    return out


def _format_line(key, val):
    if val is None:
        return key
    if isinstance(val, bool):
        return '{} .{}.'.format(key, str(val).lower())
    else:
        return ' '.join([key, str(val)])


def _fix_pseudopotentials_block(params):
    """
    Ensure pseudopotentials are always written as a block in nwpw.
    Modifies params in-place if needed.
    """
    nwpw = params.get('nwpw', {})
    psp = nwpw.get('pseudopotentials', None)
    if psp is not None:
        # If already a block dict, do nothing
        if isinstance(psp, dict) and '__block__' in psp:
            return
        block = []
        # If dict, convert to block
        if isinstance(psp, dict):
            for symbol, fname in psp.items():
                block.append(f"  {symbol} library {fname}")
        elif isinstance(psp, str):
            lines = psp.strip().split('\n')
            for line in lines:
                toks = line.strip().split()
                if len(toks) == 2:
                    block.append(f"  {toks[0]} library {toks[1]}")
        elif isinstance(psp, list):
            for entry in psp:
                if isinstance(entry, (list, tuple)) and len(entry) == 2:
                    block.append(f"  {entry[0]} library {entry[1]}")
        # Replace with block dict
        nwpw['pseudopotentials'] = {'__block__': block}
        params['nwpw'] = nwpw


def _format_block(key, val, twod_hcurve, lmbfgs, nindent=0):
    print(f"[DEBUG] _format_block called: key={key}, val={val}")
    prefix = '  ' * nindent
    prefix2 = '  ' * (nindent + 1)
    if val is None:
        return [prefix + key]

    # Special handling for pseudopotentials block
    if key == 'pseudopotentials':
        # Accept both __block__ and plain dict
        block = None
        if isinstance(val, dict) and '__block__' in val:
            block = val['__block__']
        elif isinstance(val, dict):
            block = [f"  {symbol} library {fname}" for symbol, fname in val.items()]
        if block is not None:
            out = [prefix + key]
            # Emit each line in the block as a separate line
            for line in block:
                out.append(prefix2 + line.strip())
            out.append(prefix + 'end')
            return out

    if not isinstance(val, dict):
        return [prefix + _format_line(key, val)]

    out = [prefix + key]
    if key == 'nwpw' and twod_hcurve is True:
        out.append(prefix2 + '2d-hcurve')
    if key == 'nwpw' and lmbfgs is True:
        out.append(prefix2 + 'lmbfgs')
    for subkey, subval in val.items():
        if (key, subkey) in _special_keypairs:
            if (key, subkey) == ('nwpw', 'brillouin_zone'):
                out += _format_brillouin_zone(subval)
            else:
                out += _format_block(subkey, subval, twod_hcurve, lmbfgs, nindent=nindent + 1)
        elif subkey == 'pseudopotentials' and isinstance(subval, dict) and '__block__' in subval:
            out += _format_block('pseudopotentials', subval, twod_hcurve, lmbfgs, nindent=nindent + 1)
        else:
            if isinstance(subval, dict):
                subval = ' '.join([_format_line(a, b)
                                   for a, b in subval.items()])
            out.append(prefix2 + ' '.join([_format_line(subkey, subval)]))
    out.append(prefix + 'end')
    return out


def _get_other(twod_hcurve, lmbfgs, **params):
    out = []

    for kw, block in params.items():
        if kw in _special_kws:
            continue
        out += _format_block(kw, block, twod_hcurve, lmbfgs)
    return out


def _get_set(**params):
    return ['set ' + _format_line(key, val) for key, val in params.items()]


_gto_theories = ['tce', 'ccsd', 'mp2', 'tddft', 'scf', 'dft']
_pw_theories = ['band', 'pspw', 'paw']
_all_theories = _gto_theories + _pw_theories


def _get_theory(**params):
    # Default: user-provided theory
    theory = params.get('theory')
    if theory is not None:
        return theory

    # Check if the user passed a theory to xc
    xc = params.get('xc')
    if xc in _all_theories:
        return xc

    # Check for input blocks that correspond to a particular level of
    # theory. Correlated theories (e.g. CCSD) are checked first.
    for kw in _gto_theories:
        if kw in params:
            return kw

    # If the user passed an 'nwpw' block, then they want a plane-wave
    # calculation, but what kind? If they request k-points, then
    # they want 'band', otherwise assume 'pspw' (if the user wants
    # to use 'paw', they will have to ask for it specifically).
    nwpw = params.get('nwpw')
    if nwpw is not None:
        if 'monkhorst-pack' in nwpw or 'brillouin_zone' in nwpw:
            return 'band'
        return 'pspw'

    # When all else fails, default to dft.
    return 'dft'


_xc_conv = dict(lda='slater pw91lda',
                pbe='xpbe96 cpbe96',
                revpbe='revpbe cpbe96',
                rpbe='rpbe cpbe96',
                pw91='xperdew91 perdew91',
                )


def _update_mult(magmom_tot, **params):
    theory = params['theory']
    if magmom_tot == 0:
        magmom_mult = 1
    else:
        magmom_mult = np.sign(magmom_tot) * (abs(magmom_tot) + 1)
    if 'scf' in params:
        for kw in ['nopen', 'singlet', 'doublet', 'triplet', 'quartet',
                   'quintet', 'sextet', 'septet', 'octet']:
            if kw in params['scf']:
                break
        else:
            params['scf']['nopen'] = magmom_tot
    elif theory in ['scf', 'mp2', 'ccsd', 'tce']:
        params['scf'] = dict(nopen=magmom_tot)

    if 'dft' in params:
        if 'mult' not in params['dft']:
            params['dft']['mult'] = magmom_mult
    elif theory in ['dft', 'tddft']:
        params['dft'] = dict(mult=magmom_mult)

    if 'nwpw' in params:
        if 'mult' not in params['nwpw']:
            params['nwpw']['mult'] = magmom_mult
    elif theory in ['pspw', 'band', 'paw']:
        params['nwpw'] = dict(mult=magmom_mult)

    return params


def _get_simulation_cell(atoms, **params):
    """Automatically generate simulation_cell from ASE cell if not provided"""
    nwpw = params.get('nwpw', dict())
    
    # If simulation_cell is already provided, don't override
    if 'simulation_cell' in nwpw:
        return params
    
    # Get cell dimensions
    cell = atoms.cell
    cell_diag = cell.diagonal()
    
    # Check if cell is cubic (all diagonal elements are the same)
    if np.allclose(cell_diag[0], cell_diag[1]) and np.allclose(cell_diag[1], cell_diag[2]):
        # Use SC format for cubic cells
        nwpw['simulation_cell'] = {
            'SC': cell_diag[0]  # Single value for cubic cell
        }
    else:
        # Use lattice_vectors format for non-cubic cells
        nwpw['simulation_cell'] = {
            'lattice_vectors': [
                [cell[0, 0], cell[0, 1], cell[0, 2]],
                [cell[1, 0], cell[1, 1], cell[1, 2]],
                [cell[2, 0], cell[2, 1], cell[2, 2]]
            ]
        }
    
    params['nwpw'] = nwpw
    return params


def _get_kpts(atoms, **params):
    """Converts top-level 'kpts' argument to native keywords"""
    kpts = params.get('kpts')
    if kpts is None:
        return params

    nwpw = params.get('nwpw', dict())

    if 'monkhorst-pack' in nwpw or 'brillouin_zone' in nwpw:
        raise ValueError("Redundant k-points specified!")

    if KPoints is not None and isinstance(kpts, KPoints):
        nwpw['brillouin_zone'] = kpts.kpts
    elif isinstance(kpts, dict):
        if kpts.get('gamma', False) or 'size' not in kpts:
            if kpts2kpts is not None:
                nwpw['brillouin_zone'] = kpts2kpts(kpts, atoms).kpts
            else:
                raise ValueError("kpts2kpts not available - please install ASE with k-points support")
        else:
            nwpw['monkhorst-pack'] = ' '.join(map(str, kpts['size']))
    elif isinstance(kpts, np.ndarray):
        nwpw['brillouin_zone'] = kpts
    else:
        nwpw['monkhorst-pack'] = ' '.join(map(str, kpts))

    params['nwpw'] = nwpw
    return params


def write_pwdft_in(fd, atoms, properties=None, echo=True, twod_hcurve=True, lmbfgs=True, **params):
    """Writes PWDFT input file.

    Parameters
    ----------
    fd
        file descriptor
    atoms
        atomic configuration
    properties
        list of properties to compute; by default only the
        calculation of the energy is requested
    echo
        if True include the `echo` keyword at the top of the file,
        which causes the content of the input file to be included
        in the output file
    params
        dict of instructions blocks to be included
    """
    params = deepcopy(params)

    if properties is None:
        properties = ['energy']

    if 'stress' in properties:
        if 'set' not in params:
            params['set'] = dict()
        params['set']['includestress'] = True

    task = params.get('task')
    if task is None:
        if 'stress' in properties or 'forces' in properties:
            task = 'gradient'
        else:
            task = 'energy'

    # params = _get_simulation_cell(atoms, **params)  # Disabled for now
    params = _get_kpts(atoms, **params)

    theory = _get_theory(**params)
    params['theory'] = theory
    xc = params.get('xc')
    if 'xc' in params:
        xc = _xc_conv.get(params['xc'].lower(), params['xc'])
        if theory in ['dft', 'tddft']:
            if 'dft' not in params:
                params['dft'] = dict()
            params['dft']['xc'] = xc
        elif theory in ['pspw', 'band', 'paw']:
            if 'nwpw' not in params:
                params['nwpw'] = dict()
            params['nwpw']['xc'] = xc

    magmom_tot = int(atoms.get_initial_magnetic_moments().sum())
    params = _update_mult(magmom_tot, **params)

    # --- Fix pseudopotentials block ---
    _fix_pseudopotentials_block(params)

    label = params.get('label', 'pwdft')
    perm = os.path.abspath(params.pop('perm', label))
    scratch = os.path.abspath(params.pop('scratch', label))
    restart_kw = params.get('restart_kw', 'start')
    if restart_kw not in ('start', 'restart'):
        raise ValueError("Unrecognised restart keyword: {}!"
                         .format(restart_kw))
    short_label = label.rsplit('/', 1)[-1]
    if echo:
        out = ['echo']
    else:
        out = []

    perm = './perm'
    scratch = './perm'

    out.extend(['title "{}"'.format(short_label),
                'permanent_dir {}'.format(perm),
                'scratch_dir {}'.format(scratch),
                '{} {}'.format(restart_kw, short_label),
                '\n'.join(_get_geom(atoms, **params)),
                '\n'.join(_get_other(twod_hcurve, lmbfgs, **params)),
                '\n'.join(_get_set(**params.get('set', dict()))),
                'task {} {}'.format(theory, task),
                '\n'.join(_get_bandpath(params.get('bandpath', None)))])

    # Remove any accidental editor/terminal artifacts
    out = [line for line in out if not (line.strip().startswith('1G:q') or line.strip().startswith(':q'))]

    fd.write('\n\n'.join(out))
