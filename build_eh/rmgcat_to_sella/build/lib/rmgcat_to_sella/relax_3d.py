import os

def create_relax_jobs(facetpath, pytemplate, shtemplate=None):
    with open(pytemplate, 'r') as f:
        pytemplate = f.read()

    if shtemplate is not None:
        with open(shtemplate, 'r') as f:
            shtemplate = f.read()

    for species in os.listdir(facetpath):
        speciespath = os.path.join(facetpath, species)
        if not os.path.isdir(speciespath):
            continue
        for structure in os.listdir(speciespath):
            if structure.endswith('xyz'):
                prefix = os.path.splitext(structure)[0]
                fname = os.path.join(facetpath, f'{species}_{prefix}_relax.py')
                with open(fname, 'w') as f:
                    f.write(pytemplate.format(adsorbate=species,
                                              prefix=prefix))
                if shtemplate is None:
                    continue

                fname = os.path.join(facetpath, f'{species}_{prefix}_relax.sh')
                with open(fname, 'w') as f:
                    f.write(shtemplate.format(adsorbate=species,
                                              prefix=prefix))
