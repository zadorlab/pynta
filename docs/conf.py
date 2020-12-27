import os
import sys
import sphinx_rtd_theme
# sys.path.insert(0, os.path.abspath(
#     '/Users/mgierad/00_SANDIA_WORK/05_rmgcat_to_stella/pynta/'))
# sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath('../'))

# -- Project information -----------------------------------------------------

project = 'Pynta'
copyright = '2020, Maciej Gierada'
author = 'Maciej Gierada'

# The full version, including alpha/beta/rc tags
release = 'Jan 7, 2020'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.napoleon',
              'sphinx.ext.viewcode',
              'm2r']

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

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# html_theme = 'alabaster'
html_theme = 'sphinx_rtd_theme'

html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
html_style = 'ase.css'
# html_favicon = 'static/ase.ico'
html_static_path = ['_static']
html_last_updated_fmt = '%a, %d %b %Y %H:%M:%S'

autodoc_mock_imports = ['balsam', 'inputR2S']
