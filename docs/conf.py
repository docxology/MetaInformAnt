"""
METAINFORMANT Documentation Configuration

Sphinx configuration for generating comprehensive API documentation
with auto-documentation from Google-style docstrings.
"""

import os
import sys
from pathlib import Path

# Add the src directory to the Python path for imports
sys.path.insert(0, os.path.abspath('../src'))

# -- Project information -----------------------------------------------------
project = 'METAINFORMANT'
copyright = '2024, METAINFORMANT Contributors'
author = 'METAINFORMANT Contributors'

# Get version from pyproject.toml
def get_version():
    """Extract version from pyproject.toml."""
    try:
        import tomllib
        with open('../pyproject.toml', 'rb') as f:
            data = tomllib.load(f)
            return data['project']['version']
    except Exception:
        return '0.2.0'

version = get_version()
release = version

# -- General configuration ---------------------------------------------------
extensions = [
    # Core Sphinx extensions
    'sphinx.ext.autodoc',        # Auto-documentation from docstrings
    'sphinx.ext.autosummary',    # Generate autosummary tables
    'sphinx.ext.viewcode',       # Add source code links
    'sphinx.ext.napoleon',       # Google/NumPy style docstrings
    'sphinx.ext.intersphinx',    # Cross-references to other projects
    'sphinx.ext.todo',           # TODO directives
    'sphinx.ext.coverage',       # Coverage reporting

    # Type hints support
    'sphinx_autodoc_typehints',

    # Markdown support
    'myst_parser',

    # Jupyter notebook support (optional)
    # 'nbsphinx',
]

# Source file parsing
source_suffix = {
    '.rst': None,
    '.md': None,
}

# Master document
master_doc = 'index'

# Language and locale
language = 'en'
locale_dirs = ['locale/']
gettext_compact = False

# -- Auto-documentation settings ---------------------------------------------
autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'show-inheritance': True,
    'member-order': 'bysource',
    'special-members': '__init__',
}

# Type hints settings
autodoc_typehints = 'description'
autodoc_typehints_description_target = 'documented'

# Mock imports for optional dependencies during documentation build
autodoc_mock_imports = [
    'Bio',          # Biopython
    'sklearn',      # scikit-learn
    'networkx',     # NetworkX
    'plotly',       # Plotly
    'bokeh',        # Bokeh
    'scanpy',       # ScanPy
    'anndata',      # AnnData
    'umap',         # UMAP
    'scipy',        # SciPy
    'torch',        # PyTorch
    'tensorflow',   # TensorFlow
    'pysam',        # PySAM
    'ncbi_datasets', # NCBI Datasets
]

# -- Napoleon (Google/NumPy docstring) settings ----------------------------
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_use_keyword = True
napoleon_custom_sections = None

# -- MyST (Markdown) settings -----------------------------------------------
myst_enable_extensions = [
    'colon_fence',
    'deflist',
    'dollarmath',
    'fieldlist',
    'html_admonition',
    'html_image',
    'replacements',
    'smartquotes',
    'strikethrough',
    'substitution',
    'tasklist',
]

myst_heading_anchors = 3
myst_linkify_fuzzy_links = False

# -- Intersphinx settings ---------------------------------------------------
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'pandas': ('https://pandas.pydata.org/docs/', None),
    'matplotlib': ('https://matplotlib.org/stable/', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/', None),
    'sklearn': ('https://scikit-learn.org/stable/', None),
    'networkx': ('https://networkx.org/documentation/stable/', None),
}

# -- HTML output ------------------------------------------------------------
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_css_files = []

# Theme options
html_theme_options = {
    'canonical_url': '',
    'analytics_id': '',
    'display_version': True,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    'vcs_pageview_mode': '',
    'style_nav_header_background': '#2980B9',
    'collapse_navigation': True,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False,
}

# Logo and favicon
# html_logo = '_static/logo.png'
# html_favicon = '_static/favicon.ico'

# Custom CSS
html_css_files = [
    'css/custom.css',
]

# -- Autosummary settings ---------------------------------------------------
autosummary_generate = True
autosummary_imported_members = True

# -- Coverage settings ------------------------------------------------------
coverage_modules = []
coverage_ignore_modules = []
coverage_ignore_functions = []
coverage_ignore_classes = []

# -- Todo settings ----------------------------------------------------------
todo_include_todos = True

# -- Linkcheck settings -----------------------------------------------------
linkcheck_retries = 3
linkcheck_timeout = 30
linkcheck_workers = 5
linkcheck_ignore = [
    r'https://github\.com/.*#.*',  # GitHub anchors
    r'https://pypi\.org/.*',       # PyPI links
]

# -- Options for HTML help output -------------------------------------------
htmlhelp_basename = 'metainformantdoc'

# -- Options for LaTeX output ------------------------------------------------
latex_elements = {
    'papersize': 'letterpaper',
    'pointsize': '10pt',
    'preamble': '',
    'figure_align': 'htbp',
}

latex_documents = [
    (master_doc, 'metainformant.tex', 'METAINFORMANT Documentation',
     'METAINFORMANT Contributors', 'manual'),
]

# -- Options for manual page output ------------------------------------------
man_pages = [
    (master_doc, 'metainformant', 'METAINFORMANT Documentation',
     [author], 1)
]

# -- Options for Texinfo output ----------------------------------------------
texinfo_documents = [
    (master_doc, 'metainformant', 'METAINFORMANT Documentation',
     author, 'metainformant', 'Comprehensive bioinformatics toolkit',
     'Miscellaneous'),
]

# -- Extension configuration --------------------------------------------------

# Sphinx AutoAPI (alternative to autodoc for API docs)
# autoapi_dirs = ['../src/metainformant']
# autoapi_generate_api_docs = True
# autoapi_add_toctree_entry = True

# -- Custom functions ---------------------------------------------------------

def skip_member(app, what, name, obj, skip, options):
    """Skip certain members from documentation."""
    # Skip private methods that start with single underscore
    if name.startswith('_') and not name.startswith('__'):
        return True

    # Skip certain inherited methods
    if name in ['__weakref__', '__dict__', '__module__', '__doc__']:
        return True

    return skip

def setup(app):
    """Setup function for Sphinx extensions."""
    app.connect('autodoc-skip-member', skip_member)
