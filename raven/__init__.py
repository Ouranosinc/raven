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
