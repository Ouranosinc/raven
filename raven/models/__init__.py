from .gr4j_cemaneige import simulation as gr4j
import os

_dir = os.path.abspath(os.path.dirname(__file__))

raven_templates = {'raven-gr4j-cemaneige': os.path.join(_dir, 'raven-gr4j-cemaneige'),
                   'raven-hmets': os.path.join(_dir, 'raven-hmets')}
