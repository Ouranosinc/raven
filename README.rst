===========================================
Raven : Hydrological modeling and analytics
===========================================

+----------------------------+----------------------------------------+
| Versions                   | |pypi| |conda| |platforms|             |
+----------------------------+----------------------------------------+
| Documentation and Support  | |docs| |gitter| |versions|             |
+----------------------------+----------------------------------------+
| Open Source                | |license|                              |
+----------------------------+----------------------------------------+
| Coding Standards           | |black| |ruff| |pre-commit|            |
+----------------------------+----------------------------------------+
| Development Status         | |status| |build| |coveralls|           |
+----------------------------+----------------------------------------+

Raven (the bird)
  *A WPS service that offers processes related to hydrological modelling.*

Raven is an open source server project offering data collection and preparation, as well as geoprocessing and catchment delineation through the Web Processing Service (WPS) standard. Raven processes can be embedded in a graphical user interface or accessed directly from a programming environment. From Python, birdy_ WPSClient provides a user-friendly python interface to Raven's WPS processes for geospatial processing.

The properties of custom watersheds can be extracted from a Digital Elevation Model and a land-use database.

Raven can be compiled and installed, or simply deployed using docker. A hosted version is available at https://pavics.ouranos.ca/twitcher/ows/proxy/raven.

Documentation
-------------

Learn more about Raven in its official documentation at https://pavics-raven.readthedocs.io.

Submit bug reports, questions and feature requests at https://github.com/Ouranosinc/raven/issues

Contributing
------------

You can find information about contributing in our `Developer Guide`_.

Please use `bump-my-version`_ to release a new version.

License
-------

* Free software: MIT license
* Documentation: https://raven.readthedocs.io.

Credits
-------

This project was funded by the `CANARIE`_ research software program.

Hydrological models are based on the `Raven`_ modeling framework.

This package was created with `Cookiecutter`_ and the `bird-house/cookiecutter-birdhouse`_ project template.

.. _`birdy`: https://birdy.readthedocs.io
.. _`xarray`: http://xarray.pydata.org
.. _`xclim`: https://xclim.readthedocs.io
.. _`Raven`: http://raven.uwaterloo.ca
.. _`CANARIE`: https://www.canarie.ca
.. _`Cookiecutter`: https://github.com/audreyr/cookiecutter
.. _`bird-house/cookiecutter-birdhouse`: https://github.com/bird-house/cookiecutter-birdhouse
.. _`Developer Guide`: https://pavics-raven.readthedocs.io/en/latest/dev_guide.html
.. _`bump-my-version`: https://pavics-raven.readthedocs.io/en/latest/dev_guide.html#bump-a-new-version

.. |black| image:: https://img.shields.io/badge/code%20style-black-000000.svg
        :target: https://github.com/psf/black
        :alt: Python Black

.. |build| image:: https://github.com/Ouranosinc/raven/actions/workflows/main.yml/badge.svg
    :target: https://github.com/Ouranosinc/raven/actions/workflows/main.yml
    :alt: Build Status

.. |conda| image:: https://img.shields.io/conda/vn/conda-forge/raven-wps.svg
    :target: https://anaconda.org/conda-forge/raven-wps
    :alt: Anaconda Version

.. |coveralls| image:: https://coveralls.io/repos/github/Ouranosinc/raven/badge.svg
    :target: https://coveralls.io/github/Ouranosinc/raven
    :alt: Coveralls

.. |docs| image:: https://readthedocs.org/projects/pavics-raven/badge/?version=latest
    :target: https://pavics-raven.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. |fossa| image:: https://app.fossa.com/api/projects/git%2Bgithub.com%2FOuranosinc%2Fraven.svg?type=shield
    :target: https://app.fossa.com/projects/git%2Bgithub.com%2FOuranosinc%2Fraven?ref=badge_shield
    :alt: FOSSA

.. |gitter| image:: https://badges.gitter.im/bird-house/birdhouse.svg
    :target: https://gitter.im/bird-house/birdhouse?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge
    :alt: Join the chat at https://gitter.im/bird-house/birdhouse

.. |license| image:: https://img.shields.io/github/license/Ouranosinc/raven.svg
    :target: https://github.com/Ouranosinc/raven/blob/main/LICENSE
    :alt: License

.. |platforms| image:: https://anaconda.org/conda-forge/raven-wps/badges/platforms.svg
    :target: https://anaconda.org/conda-forge/raven-wps
    :alt: Supported Python Versions

.. |pre-commit| image:: https://results.pre-commit.ci/badge/github/Ouranosinc/raven/main.svg
    :target: https://results.pre-commit.ci/latest/github/Ouranosinc/raven/main
    :alt: pre-commit.ci status

.. |pypi| image:: https://img.shields.io/pypi/v/birdhouse-raven.svg
    :target: https://pypi.python.org/pypi/birdhouse-raven
    :alt: PyPI

.. |ruff| image:: https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json
    :target: https://github.com/astral-sh/ruff
    :alt: Ruff

.. |status| image:: https://www.repostatus.org/badges/latest/active.svg
    :target: https://www.repostatus.org/#active
    :alt: Project Status: Active - The project has reached a stable, usable state and is being actively developed.

.. |versions| image:: https://img.shields.io/pypi/pyversions/birdhouse-raven.svg
    :target: https://pypi.python.org/pypi/birdhouse-raven
    :alt: Supported Python Versions

.. |zenodo| image:: https://zenodo.org/badge/135511617.svg
    :target: https://zenodo.org/badge/latestdoi/135511617
    :alt: DOI
