from .gr4j_cemaneige import simulation as gr4j
import os
from .basic import Raven, GR4JCemaneige, MOHYSE, HMETS
from .rv import RV

_dir = os.path.abspath(os.path.dirname(__file__))

raven_templates = {'raven-gr4j-cemaneige': os.path.join(_dir, 'raven-gr4j-cemaneige'),
                   'raven-mohyse': os.path.join(_dir, 'raven-mohyse'),
                   'raven-hmets': os.path.join(_dir, 'raven-hmets')}
