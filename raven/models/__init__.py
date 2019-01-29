from .gr4j_cemaneige import simulation as gr4j
import os
from .basic import Raven, GR4JCN, MOHYSE, HMETS, HBVEC
from .basic_ostrich import Ostrich, GR4JCN_OST
from .rv import RV, RVI

_dir = os.path.abspath(os.path.dirname(__file__))

raven_templates = {'raven-gr4j-cemaneige': os.path.join(_dir, 'raven-gr4j-cemaneige'),
                   'raven-mohyse': os.path.join(_dir, 'raven-mohyse'),
                   'raven-hmets': os.path.join(_dir, 'raven-hmets'),
                   'raven-hbv-ec': os.path.join(_dir, 'raven-hbv-ec')}

ostrich_templates = {'ostrich-gr4j-cemaneige': os.path.join(_dir, 'ostrich-gr4j-cemaneige')}


def get_model(name):
    """Return the corresponding Raven emulated model instance.

    Parameters
    ----------
    name : str
      Model name.

    Returns
    -------
    Raven model instance
    """
    from raven.models import basic
    model_cls = getattr(basic, name, None)

    if model_cls is None:
        raise ValueError("Model {} is not recognized.".format(model_cls))

    return model_cls()
