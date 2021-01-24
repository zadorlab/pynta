import os
import sys
import sphinx_rtd_theme
from datetime import date
# sys.path.insert(0, os.path.abspath(
#     '/Users/mgierad/00_SANDIA_WORK/05_rmgcat_to_stella/pynta/'))
# sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath('../'))

# -- Project information -----------------------------------------------------

project = 'Pynta'
copyright = 'Copyright 2021 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.'
author = 'Maciej Gierada, Eric D. Hermes, David H. Bross'

# The full version, including alpha/beta/rc tags
# release = 'Jan 7, 2020'
today = date.today()
release = today.strftime('%B %d, %Y')


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.napoleon',
              'sphinx.ext.viewcode',
              'sphinx.ext.todo',
              'sphinx.ext.mathjax']
# 'm2r'

# extensions = ['sphinx.ext.autodoc', 'sphinx.ext.napoleon']
pygments_style = 'sphinx'

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

napoleon_google_docstring = False
napoleon_use_param = True
napoleon_use_ivar = True
napoleon_include_init_with_doc = True
todo_include_todos = True

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# html_theme = 'alabaster'
html_theme = 'sphinx_rtd_theme'

html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
# html_style = 'ase.css'
# html_favicon = 'static/ase.ico'
html_static_path = ['_static']
html_last_updated_fmt = '%a, %d %b %Y %H:%M:%S'

autodoc_mock_imports = [
    'balsam', 'inputR2S.optimize_slab', 'inputR2S.surface_types_and_repeats',
    'inputR2S.a', 'inputR2S.vacuum', 'inputR2S.pseudo_dir',
    'inputR2S.metal_atom', 'inputR2S.pseudopotentials', 'inputR2S.yamlfile',
    'inputR2S.scfactor', 'inputR2S.scaled1', 'inputR2S.scaled2',
    'inputR2S.executable', 'inputR2S.node_packing_count',
    'inputR2S.balsam_exe_settings', 'inputR2S.calc_keywords',
    'inputR2S.creation_dir']

# # -- Options for LaTeX output ---------------------------------------------

# latex_elements = {
#     # Additional stuff for the LaTeX preamble.
#     'preamble': '\\usepackage{typo3}'
# }
