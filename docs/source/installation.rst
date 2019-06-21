.. _installation:

Installation
============

Install from Anaconda
---------------------

.. todo::

   Prepare Conda package.

Install from GitHub
-------------------

Check out code from the Raven GitHub repo and start the installation:

.. code-block:: sh

   $ git clone https://github.com/huard/raven.git
   $ cd raven
   $ conda env create -f environment.yml
   $ source activate raven
   $ python setup.py develop

... or do it the lazy way
+++++++++++++++++++++++++

The previous installation instructions assume you have Anaconda installed.
We provide also a ``Makefile`` to run this installation without additional steps:

.. code-block:: sh

   $ git clone https://github.com/huard/raven.git
   $ cd raven
   $ make clean    # cleans up a previous Conda environment
   $ make install  # installs Conda if necessary and runs the above installation steps

Start Raven PyWPS service
-------------------------

After successful installation you can start the service using the ``raven`` command-line.

.. code-block:: sh

   $ raven --help # show help
   $ raven start  # start service with default configuration

   OR

   $ raven start --daemon # start service as daemon
   loading configuration
   forked process id: 42

The deployed WPS service is by default available on:

http://localhost:9099/wps?service=WPS&version=1.0.0&request=GetCapabilities.

.. NOTE:: Remember the process ID (PID) so you can stop the service with ``kill PID``.

Check the log files for errors:

.. code-block:: sh

   $ tail -f  pywps.log

... or do it the lazy way
+++++++++++++++++++++++++

You can also use the ``Makefile`` to start and stop the service:

.. code-block:: sh

  $ make start
  $ make status
  $ tail -f pywps.log
  $ make stop

..
    Run Raven as Docker container
    -----------------------------

    You can also run Raven as a Docker container, see the :ref:`Tutorial <tutorial>`.


Use Ansible to deploy Raven on your System
------------------------------------------

Use the `Ansible playbook`_ for PyWPS to deploy Raven on your system.
Follow the `example`_ for Raven given in the playbook.

Building the docs
-----------------

First install dependencies for the documentation::

  $ make bootstrap_dev
  $ make docs


.. _Ansible playbook: http://ansible-wps-playbook.readthedocs.io/en/latest/index.html
.. _example: http://ansible-wps-playbook.readthedocs.io/en/latest/tutorial.html
