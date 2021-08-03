.. _installation:

Installation
============

.. contents::
    :local:
    :depth: 1

Install from Conda-Forge (suggested)
------------------------------------

.. warning::
  The Raven-WPS conda package is not yet available on conda-forge. This documentation will be updated when it is.

Create an Anaconda environment named `ravenwps-env`:

.. code-block:: console

   $ conda env create -n ravenwps-env python=3.7
   $ source activate ravenwps-env

This should now prepend the environment to your shell commands (ie: `(ravenwps-env) $`).
Now install directly from `conda-forge`:

.. code-block:: console

   (ravenwps-env) $ conda install -c conda-forge raven-wps

Install from GitHub
-------------------

Check out code from the Raven GitHub repo and start the installation:

.. code-block:: console

   $ git clone https://github.com/Ouranosinc/raven.git
   $ cd raven

Environment Setup with Anaconda (macOS/Linux)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create Conda environment named `raven`:

.. code-block:: console

   $ conda env create -n ravenwps-env -f environment.yml
   # or alternatively,
   $ make conda_env

The environment can then be activated with:

.. code-block::

   $ source activate ravenwps-env

This should now prepend the environment to your shell commands (ie: `(ravenwps-env) $`).

Environment Setup using System Libraries (Linux)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. warning::
    This approach is not formally supported, but is presently working for the time being.
    It is up to the user to install the `raven` model and `ostrich` model optimization binaries.
    Those can be downloaded from source via the following links:
        - RAVEN: http://raven.uwaterloo.ca/Downloads.html
        - OSTRICH: https://github.com/usbr/ostrich/

.. note::
   While the following shows how to install `raven` for an Deb-based Linux, if the OS-equivalent dependencies
   are available to Python, `raven` should be able to run on any modern operating system (macOS/Windows/*nix).

First we need to install several system libraries that RavenWPS and RavenPy depend upon and make a virtual environment:

.. code-block:: console

   $ sudo apt-get install libhdf5-dev netcdf-bin libnetcdf-dev libgdal-dev libproj-dev libgeos-dev libspatialindex-dev python3-dev
   $ pip3 install virtualenv
   $ virtualenv ravenwps-env
   $ . ravenwps-env/bin/activate

We then need to install the python RavenPy library from sources:

.. code-block:: console

   (ravenwps-env) $ git clone https://github.com/CSHS-CWRA/RavenPy/
   (ravenwps-env) $ pip install RavenPy/.[gis]
   (ravenwps-env) $ pip install RavenPy/. --verbose --install-option="--with-binaries"

Installing and Launching RavenWPS
---------------------------------

Now we can install the raven-wps service:

.. code-block:: console

  (ravenwps-env) $ pip install -e .
  # or alternatively,
  (ravenwps-env) $ make install

For development you can use this command:

.. code-block:: console

  (ravenwps-env) $ pip install -e .[dev]
  # or alternatively,
  (ravenwps-env) $ make develop

Then clone the Raven Test Data repo somewhere on your disk:

.. code-block:: console

   (ravenwps-env) $ git clone https://github.com/Ouranosinc/raven-testdata.git

You can then run the test suite by doing:

.. code-block:: console

   (ravenwps-env) $ export RAVENPY_TESTDATA_PATH=/path/to/raven-testdata
   (ravenwps-env) $ pytest

Start Raven PyWPS service
-------------------------

After successful installation you can start the service using the ``raven`` command-line.

.. code-block:: console

   (ravenwps-env) $ raven-wps --help # show help
   (ravenwps-env) $ raven-wps start  # start service with default configuration
   # or alternatively,
   (ravenwps-env) $ raven-wps start --daemon # start service as daemon
   loading configuration
   forked process id: 42

The deployed WPS service is by default available on:

http://localhost:9099/wps?service=WPS&version=1.0.0&request=GetCapabilities.

You can find which process uses a given port using the following command (here for port 5000):

.. code-block:: console

   $ netstat -nlp | grep :5000


Check the log files for errors:

.. code-block:: console

   $ tail -f  pywps.log

... or do it the lazy way
+++++++++++++++++++++++++

You can also use the ``Makefile`` to start and stop the service:

.. code-block:: console

  (ravenwps-env) $ make start
  (ravenwps-env) $ make status
  (ravenwps-env) $ tail -f pywps.log
  (ravenwps-env) $ make stop

..
    Run Raven as Docker container
    -----------------------------

    You can also run Raven as a Docker container, see the :ref:`Tutorial <tutorial>`.

You can also run Raven as a Docker container.

.. warning::

  TODO: Describe Docker container support.

Use Ansible to deploy Raven on your System
------------------------------------------

Use the `Ansible playbook`_ for PyWPS to deploy Raven on your system.

.. _Ansible playbook: http://ansible-wps-playbook.readthedocs.io/en/latest/index.html
