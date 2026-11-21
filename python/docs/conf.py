# SPDX-License-Identifier: MIT
"""Sphinx configuration for the FVQual Python documentation.

Mirrors the openswmm.engine documentation setup — pydata theme, MyST, design
grids, intersphinx — so the two sites read the same way and can link to each
other without a reader noticing a seam.

**Status.** The compiled ``fvqual`` package does not exist yet; the C API it
will bind is in ``include/fvqual/capi/``. Nothing here uses ``autodoc`` against
a module that is not built, because a documentation site whose build depends on
a package that has not been written is a site that is red from the first day
and stays that way. When the extension lands, add ``autodoc`` and the API pages
alongside these.
"""

project = "FVQual"
author = "Caleb Buahin"
copyright = "2026, Caleb Buahin"

extensions = [
    "myst_parser",
    "sphinx.ext.intersphinx",
    "sphinx.ext.extlinks",
    "sphinx_design",
]

templates_path = ["_templates"]
exclude_patterns = ["_build"]

language = "en"
pygments_style = "sphinx"

# Cross-links to the sibling projects, so the C++ API, the SWMM engine, and
# these pages form one set of documentation rather than three islands.
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
}

extlinks = {
    "cxx": ("https://hydrocouple.github.io/FVQual/%s", "C++ API: %s"),
    "hydrocouple": ("https://hydrocouple.github.io/HydroCouple/%s", "HydroCouple: %s"),
    "sdk": ("https://hydrocouple.github.io/HydroCoupleSDK/%s", "HydroCoupleSDK: %s"),
    "openswmm": ("https://hydrocouple.github.io/openswmm.engine/%s", "openswmm: %s"),
}

html_theme = "pydata_sphinx_theme"
html_static_path = ["_static"]
html_css_files = ['hydrocouple-palette.css']

html_theme_options = {
    "github_url": "https://github.com/hydrocouple/FVQual",
    "icon_links": [
        {
            "name": "C++ API (Doxygen)",
            "url": "https://hydrocouple.github.io/FVQual/",
            "icon": "fa-solid fa-code",
        },
        {
            "name": "HydroCouple",
            "url": "https://github.com/HydroCouple/HydroCouple",
            "icon": "fa-brands fa-github",
        },
    ],
    "navigation_with_keys": False,
}

myst_enable_extensions = ["colon_fence", "deflist"]
