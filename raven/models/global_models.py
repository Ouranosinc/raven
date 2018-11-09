from .basic import Raven
from pathlib import Path


class GR4JCemaneige(Raven):
    templates = tuple((Path(__file__) / 'raven-gr4j-cemaneige').glob("*.rv?"))

