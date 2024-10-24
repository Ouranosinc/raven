.. _installation:

Installation
============

.. contents::
    :local:
    :depth: 1

Install from Conda-Forge (suggested)
------------------------------------

Create an Anaconda environment named `ravenwps-env`:

.. code-block:: bash

    conda env create -n ravenwps-env python=3.9
    source activate ravenwps-env

This should now prepend the environment to your shell commands (ie: `(ravenwps-env) $`).
Now install directly from `conda-forge`:

.. code-block:: bash

    (ravenwps-env) conda install -c conda-forge raven-wps

Install from GitHub
-------------------

Check out code from the Raven GitHub repo and start the installation:

.. code-block:: bash

    git clone https://github.com/Ouranosinc/raven.git
    cd raven

Environment Setup with Anaconda (macOS/Linux)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create Conda environment named `raven`:

.. code-block:: bash

    conda env create -n ravenwps-env -f environment.yml

The environment can then be activated with:

.. code-block:: bash

    source activate ravenwps-env

This should now prepend the environment to your shell commands (ie: `(ravenwps-env) $`).

Installing and Launching RavenWPS
---------------------------------

Now we can install the raven-wps service:

    .. code-block:: bash

        pip install -e .

    Or, alternatively:

    .. code-block:: bash

        make install

For development you can use this command:

    .. code-block:: bash

    pip install -e .[dev]

    Or, alternatively:

    .. code-block:: bash

        make develop

Start Raven PyWPS service
~~~~~~~~~~~~~~~~~~~~~~~~~

After successful installation you can start the service using the ``raven`` command-line:

    .. code-block:: bash

        (ravenwps-env) $ raven-wps --help # show help
        (ravenwps-env) $ raven-wps start  # start service with default configuration

    Or, alternatively:

    .. code-block:: bash

        (ravenwps-env) $ raven-wps start --daemon # start service as daemon
        loading configuration
        forked process id: 42

The deployed WPS service is by default available on:

http://localhost:9099/wps?service=WPS&version=1.0.0&request=GetCapabilities.

You can find which process uses a given port using the following command (here for port 5000):

.. code-block:: bash

    netstat -nlp | grep :5000

Check the log files for errors:

.. code-block:: bash

    tail -f pywps.log

... or do it the lazy way
~~~~~~~~~~~~~~~~~~~~~~~~~

You can also use the ``Makefile`` to start and stop the service:

.. code-block:: bash

    (ravenwps-env) make start
    (ravenwps-env) make status
    (ravenwps-env) tail -f pywps.log
    (ravenwps-env) make stop

..
    Run Raven as Docker container
    -----------------------------

    You can also run Raven as a Docker container, see the :ref:`Tutorial <tutorial>`.

You can also run Raven as a Docker container.

.. ::

  TODO: Describe Docker container support.
