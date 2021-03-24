Raven : Hydrological modeling and analytics
===========================================

.. image:: https://readthedocs.org/projects/pavics-raven/badge/?version=latest
    :target: https://pavics-raven.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://travis-ci.com/Ouranosinc/raven.svg?branch=master
   :target: https://travis-ci.com/Ouranosinc/raven
   :alt: Travis Build

.. image:: https://img.shields.io/github/license/Ouranosinc/raven.svg
    :target: https://github.com/Ouranosinc/raven/blob/master/LICENSE.txt
    :alt: GitHub license

.. image:: https://badges.gitter.im/bird-house/birdhouse.svg
    :target: https://gitter.im/bird-house/birdhouse?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge
    :alt: Join the chat at https://gitter.im/bird-house/birdhouse

.. image:: https://zenodo.org/badge/135511617.svg
   :target: https://zenodo.org/badge/latestdoi/135511617
   :alt: DOI

Raven (the bird)
  *Raven offers processes related to hydrological modeling, and in particular, the Raven hydrological modeling framework.*

Raven is an open source server project offering hydrological modeling and analysis capabilities through the Web Processing Service (WPS) standard. Raven processes can be embedded in a graphical user interface or accessed directly from a programming environment. From Python, birdy_ WPSClient provides a user-friendly python interface to Raven's WPS processes.

Raven was made to help scientists run hydrological modeling experiments with climate change projections. It includes four lumped daily hydrological models (GR4J-CN, HBV-EC, HMETS, MOHYSE) that can be run in multi-model experiments. Meteorological input variables as well as streamflow and storage outputs use the netCDF format. Raven bundles model calibration processes, time series analysis (with xarray_), hydrological indicators and frequency analysis (using xclim_). On top of this, a database of pre-calibrated model parameters over North America is available to perform model regionalization, allowing simulations in watersheds with no streamflow observations. The properties of custom watersheds can be extracted from a Digital Elevation Model and a land-use database.

Raven can be compiled and installed, or simply deployed using docker. A hosted version is available at  https://pavics.ouranos.ca/twitcher/ows/proxy/raven.


Documentation
-------------

Learn more about Raven in its official documentation at
https://pavics-raven.readthedocs.io.

Submit bug reports, questions and feature requests at
https://github.com/Ouranosinc/raven/issues

Contributing
------------

You can find information about contributing in our `Developer Guide`_.

Please use bumpversion_ to release a new version.

License
-------

Free software: MIT license

Credits
-------

This project was funded by the CANARIE_ research software program.

Hydrological models are based on the `Raven`_ modeling framework.

This package was created with Cookiecutter_ and the `bird-house/cookiecutter-birdhouse`_ project template.

.. _`birdy`: https://birdy.readthedocs.io
.. _`xarray`: http://xarray.pydata.org
.. _`xclim`: https://xclim.readthedocs.io
.. _`Raven`: http://raven.uwaterloo.ca
.. _`CANARIE`: https://www.canarie.ca
.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`bird-house/cookiecutter-birdhouse`: https://github.com/bird-house/cookiecutter-birdhouse
.. _`Developer Guide`: https://pavics-raven.readthedocs.io/en/latest/dev_guide.html
.. _bumpversion: https://pavics-raven.readthedocs.io/en/latest/dev_guide.html#bump-a-new-version
