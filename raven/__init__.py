# -*- coding: utf-8 -*-

"""Top-level package for Raven."""

import os
from .wsgi import application
from pathlib import Path
import warnings

__author__ = """David Huard"""
__email__ = 'huard.david@ouranos.ca'
__version__ = '0.8.0'

if 'DO_NOT_CHECK_EXECUTABLE_EXISTENCE' not in os.environ:
    raven_exec = Path(__file__).parent.parent / 'bin' / 'raven'
    if not raven_exec.exists():
        raise IOError("The raven executable is not installed.")

    ostrich_exec = Path(__file__).parent.parent / 'bin' / 'ostrich'
    if not ostrich_exec.exists():
        raise IOError("The ostrich executable is not installed.")

raven_simg = Path(__file__).parent.parent / 'bin' / 'hydro-raven-latest.simg'
if not raven_simg.exists():
    warnings.warn("The Raven Singularity image has not been downloaded. Execute \n"
                  "$ singularity pull shub://132.217.141.54/hydro/raven:latest \n"
                  "and store the image in raven/bin/")
