"""Web Service Gateway Interface for PyWPS processes."""

import os
from pathlib import Path

from pywps.app.Service import Service

from .processes import processes


def create_app(cfgfiles=None):
    """Create PyWPS application."""
    config_files = [Path(__file__).parent.joinpath("default.cfg")]
    if cfgfiles:
        config_files.extend(cfgfiles)
    if "PYWPS_CFG" in os.environ:
        config_files.append(os.environ["PYWPS_CFG"])
    service = Service(processes=processes, cfgfiles=config_files)
    return service


application = create_app()
