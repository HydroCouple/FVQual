# SPDX-License-Identifier: MIT
"""Sphinx configuration for the FVQual Python documentation.

Mirrors the openswmm.engine documentation setup — pydata theme, MyST, design
grids, intersphinx — so the two sites read the same way and can link to each
other without a reader noticing a seam.

The Python reference is generated from the ctypes bindings in python/fvqual.
"""

project = "FVQual"
author = "Caleb Buahin"
copyright = "2026, Caleb Buahin"

extensions = [
    "sphinx.ext.viewcode",
    "autoapi.extension",
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

html_theme_options = {
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

# Parse the source tree so optional engines and services need not run during a docs build.
autoapi_dirs = ['../fvqual']
autoapi_options = ["members", "undoc-members", "show-module-summary"]
autoapi_add_toctree_entry = False


# -- HydroCouple documentation theme ---------------------------------------
import os  # noqa: E402
import sys  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "_hydrocouple"))
import hydrocouple_theme  # noqa: E402

hydrocouple_theme.configure(
    globals(),
    manual_id="FVQual.python",
    project_label="FVQual · Python",
    repository="https://github.com/HydroCouple/FVQual",
    hero={
        "kicker": "FVQual · Python",
        "title": "Drive FVQual from a script",
        "lede": "The fvqual package opens, steps and inspects FVQual's three-dimensional "
                "hydrodynamic and water-quality simulations through its C API, using "
                "standard-library ctypes.",
        "actions": [
            {"label": "Read the design", "url": "design", "style": "solid"},
            {"label": "API reference", "url": "api", "style": "ghost"},
            {"label": "C++ reference", "url": "https://www.hydrocouple.org/FVQual/", "style": "ghost"},
        ],
    },
)
