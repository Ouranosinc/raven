# -*- coding: utf-8 -*-

"""Top-level package for Raven."""

from .wsgi import application
from pathlib import Path

__author__ = """David Huard"""
__email__ = 'huard.david@ouranos.ca'
__version__ = '0.1.0'

raven_exec = Path(__file__).parent.parent / 'bin' / 'raven'
if not raven_exec.exists:
    raise IOError("The raven executable is not installed.")

ostrich_exec = Path(__file__).parent.parent / 'bin' / 'ostrich'
if not ostrich_exec.exists:
    raise IOError("The ostrich executable is not installed.")

raven_simg = Path(__file__).parent.parent / 'bin' / 'hydro-raven-latest.simg'
if not raven_simg.exists:
    raise IOError("The Raven Singularity image has not been downloaded. Execute \n"
                  "$ singularity pull shub://132.217.141.54/hydro/raven:latest \n"
                  "and store the image in raven/bin/")
