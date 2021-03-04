#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

import os

from setuptools import find_packages, setup

here = os.path.abspath(os.path.dirname(__file__))
README = open(os.path.join(here, "README.rst")).read()
CHANGES = open(os.path.join(here, "CHANGES.rst")).read()
REQUIRES_PYTHON = ">=3.6.0"

about = {}
with open(os.path.join(here, "raven", "__version__.py"), "r") as f:
    exec(f.read(), about)

reqs = [line.strip() for line in open("requirements.txt")]
dev_reqs = [line.strip() for line in open("requirements_dev.txt")]

docs_reqs = [
    "sphinx>=1.7",
    "sphinx-autoapi",
    "nbsphinx",
    "sphinx_rtd_theme"
]

classifiers = [
    "Development Status :: 3 - Alpha",
    "Intended Audience :: Developers",
    "Intended Audience :: Science/Research",
    "Operating System :: MacOS :: MacOS X",
    "Operating System :: POSIX",
    "Programming Language :: Python",
    "Natural Language :: English",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.6",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: 3.8",
    "Topic :: Scientific/Engineering :: Atmospheric Science",
    "License :: OSI Approved :: MIT License",
]

setup(
    name="raven",
    version=about["__version__"],
    description="Raven offers processes related to hydrological modeling, and in particular, the Raven  "
    "hydrological modeling framework.",
    long_description=README + "\n\n" + CHANGES,
    long_description_content_type="text/x-rst",
    author=about["__author__"],
    author_email=about["__email__"],
    url="https://github.com/Ouranosinc/raven",
    python_requires=REQUIRES_PYTHON,
    classifiers=classifiers,
    license="MIT license",
    keywords="wps pywps birdhouse raven hydrology gis",
    packages=find_packages(),
    include_package_data=True,
    install_requires=reqs,
    extras_require={
        "dev": dev_reqs,  # pip install ".[dev]"
        "docs": docs_reqs,
    },
    entry_points={
        "console_scripts": [
            "raven-wps=raven.cli:cli",
        ]
    },
)
