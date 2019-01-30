# -*- coding: utf-8 -*-

"""Top-level package for Raven."""

from .wsgi import application
from pathlib import Path

__author__ = """David Huard"""
__email__ = 'huard.david@ouranos.ca'
__version__ = '0.1.0'

raven_exec = str(Path(__file__).parent.parent / 'bin' / 'raven')
ostrich_exec = str(Path(__file__).parent.parent / 'bin' / 'raven')
