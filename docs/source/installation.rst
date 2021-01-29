.. _installation:

Installation
============

.. contents::
    :local:
    :depth: 1

Install from Conda
------------------

.. warning::

   TODO: Prepare Conda package.

Install from GitHub
-------------------

Check out code from the Raven GitHub repo and start the installation:

.. code-block:: console

   $ git clone https://github.com/Ouranosinc/raven.git
   $ cd raven

Create Conda environment named `raven`:

.. code-block:: console

   $ conda env create -f environment.yml
   $ source activate raven

Install Raven app:

.. code-block:: console

  $ pip install -e .
  OR
  make install

For development you can use this command:

.. code-block:: console

  $ pip install -e .[dev]
  OR
  $ make develop

Then clone the Raven Test Data repo somewhere on your disk:

.. code-block:: console

   (ravenpy-env) $ git clone git@github.com:Ouranosinc/raven-testdata.git

You can then run the test suite by doing:

.. code-block:: console

   (ravenpy-env) $ RAVENPY_TESTDATA_PATH=/path/to/raven-testdata pytest

Start Raven PyWPS service
-------------------------

After successful installation you can start the service using the ``raven`` command-line.

.. code-block:: console

   $ raven-wps --help # show help
   $ raven-wps start  # start service with default configuration

   OR

   $ raven-wps start --daemon # start service as daemon
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

  $ make start
  $ make status
  $ tail -f pywps.log
  $ make stop

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
