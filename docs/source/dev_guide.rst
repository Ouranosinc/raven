.. _devguide:

Developer Guide
===============

.. contents::
    :local:
    :depth: 1

.. WARNING:: To create new processes look at examples in Emu_.


Re-create a fresh environment
-----------------------------

.. code-block:: console

  $ make stop  # in case you previously did 'make start'
  $ conda deactivate  # exit the current 'raven' conda env so we can destroy it
  $ conda env remove -n raven  # destroy the current conda env to recreate one from scratch
  $ conda env create -f environment.yml
  $ conda activate raven
  $ make develop  # install raven-wps and additional dev tools


Building the docs
-----------------

First install dependencies for the documentation:

.. code-block:: console

  $ make develop

Run the Sphinx docs generator:

.. code-block:: console

  $ make docs

.. _testing:

Running tests
-------------

Run tests using pytest_.

First activate the ``raven`` Conda environment and install ``pytest``.

.. code-block:: console

   $ source activate raven
   $ pip install -r requirements_dev.txt  # if not already installed
   OR
   $ make develop

Run quick tests (skip slow and online):

.. code-block:: console

    $ pytest -m 'not slow and not online'"

Run all tests:

.. code-block:: console

    $ pytest

Check pep8:

.. code-block:: console

    $ flake8

Run tests the lazy way
----------------------

Do the same as above using the ``Makefile``.

.. code-block:: console

    $ make test
    $ make test-all
    $ make lint


Running notebooks tests
-----------------------

Assuming that the ``raven`` conda env has already been created and is up-to-date and
raven-wps has been installed with ``make develop``.

.. code-block:: console

    # start local raven-wps server to test against
    $ make start  # remember to make stop once done

    # to test all notebooks
    $ make test-notebooks

    OR

    # to test a single notebook (note the .run at the end of the notebook path)
    $ make docs/source/notebooks/Assess_probabilistic_flood_risk.ipynb.run

The notebooks may also require other WPS services (Finch and Flyingpigeon).  By
default these are from the production server but we can point the notebooks to
local servers if needed for development purposes:

.. code-block:: console

    # to test all notebooks
    $ make FLYINGPIGEON_WPS_URL=http://localhost:8093 FINCH_WPS_URL=http://localhost:5000 test-notebooks

    OR

    # to test a single notebook (note the .run at the end of the notebook path)
    $ make FLYINGPIGEON_WPS_URL=http://localhost:8093 FINCH_WPS_URL=http://localhost:5000 docs/source/notebooks/Assess_probabilistic_flood_risk.ipynb.run

If instead we want to run the notebooks against the production raven-wps server
or any other raven-wps servers:

.. code-block:: console

    # to test all notebooks
    $ make WPS_URL=https://pavics.ouranos.ca/twitcher/ows/proxy/raven/wps test-notebooks

    OR

    # to test juste 1 notebook (note the .run at the end of the notebook path)
    $ make WPS_URL=https://pavics.ouranos.ca/twitcher/ows/proxy/raven/wps docs/source/notebooks/Assess_probabilistic_flood_risk.ipynb.run

We can also override all three of the server variables (WPS_URL, FINCH_WPS_URL,
FLYINGPIGEON_WPS_URL) and be able to choose and pick any servers anywhere we
want.


Starting local Jupyter server to edit/develop notebooks
-------------------------------------------------------

Asume the ``raven`` conda env has already been created and is up-to-date and
raven-wps has been installed with ``make develop``.

.. code-block:: console

    # start local raven-wps server to test against
    $ make start  # remember to make stop once done

    # to start local jupyter notebook server listing all current notebooks
    $ make notebook  # Control-C to terminate once done

    # Can also override all 3 WPS_URL, FINCH_WPS_URL and FLYINGPIGEON_WPS_URL here as well,
    # just like 'make test-notebooks' to be able to pick and choose any servers anywhere we want.

    # By overridding these variables at 'make notebook' step we will not need to
    # override them one by one in each notebook because each notebook also look
    # for those variables as environment variables.


Bulk refresh all notebooks output
---------------------------------

This automated refresh only works for notebooks that passed ``make
test-notebooks`` above.  For those that failed, manually starting a local
Jupyter server and refresh them manually.

Asume the ``raven`` conda env has already been created and is up-to-date and
raven-wps has been installed with ``make develop``.

.. code-block:: console

    # start local raven-wps server to test against
    $ make start  # remember to make stop once done

    # to refresh all notebooks
    $ make refresh-notebooks

    OR

    # to refresh juste 1 notebook (note the .refresh at the end of the notebook path)
    $ make docs/source/notebooks/Assess_probabilistic_flood_risk.ipynb.refresh

    # Can also override all 3 WPS_URL, FINCH_WPS_URL and FLYINGPIGEON_WPS_URL here as well,
    # just like 'make test-notebooks' to be able to pick and choose any servers anywhere we want.


Prepare a release
-----------------

Update the Conda specification file to build identical environments_ on a specific OS.

.. note:: You should run this on your target OS, in our case Linux.

.. code-block:: console

  $ conda env create -f environment.yml
  $ source activate raven
  $ make clean
  $ make install
  $ conda list -n raven --explicit > spec-file.txt

.. _environments: https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#building-identical-conda-environments


Bump a new version
------------------

Make a new version of Raven in the following steps:

* Make sure everything is commit to GitHub.
* Update ``CHANGES.rst`` with the next version.
* Dry Run: ``bumpversion --dry-run --verbose --new-version 0.8.1 patch``
* Do it: ``bumpversion --new-version 0.8.1 patch``
* ... or: ``bumpversion --new-version 0.9.0 minor``
* Push it: ``git push``
* Push tag: ``git push --tags``

See the bumpversion_ documentation for details.

.. _bumpversion: https://pypi.org/project/bumpversion/
.. _pytest: https://docs.pytest.org/en/latest/
.. _Emu: https://github.com/bird-house/emu
