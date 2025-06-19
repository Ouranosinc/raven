Raven : Hydrological modeling and analytics
===========================================

Raven is a server providing access to hydrological modeling and analysis through the Web Processing Service (WPS) standard. It was made to help hydrologists work with climate data and automate some of the tedious work required to calibrate and run hydrological models. Its modeling engine is the Raven_ hydrological modeling framework, which can emulate a variety of lumped and distributed models.

Raven is part of the birdhouse_, a community building a suite of WPS servers supporting climate services and sustainable development goals. The idea is that instead of downloading large volumes of data locally and then analyzing it, part of the analysis can be done remotely by servers, close to the source data. Instead of users sending a plain download request, users send a request for pre-processed data. Work with low-added value can be delegated to a server, and the real research is performed on reduced datasets on local machines.

In this model, scientists need to interact closely with a server to submit requests and poll the server for its response once the job is complete. Because the boilerplate code and formats used to communicate with a server can detract from the science, we've built a generic WPS client interface (see birdy_) that hides the WPS protocol behind a native looking python interface. Remote WPS processes can be called just like a python function, returning an asynchronous response whose progress can be easily monitored.

User documentation
------------------

.. toctree::
   :maxdepth: 1

   getting_started
   installation
   configuration
   notebooks/index
   dev_guide
   processes
   authors
   changelog

Credits
-------
This project was funded by the CANARIE_ research software program.

Hydrological models are based on the `Raven`_ modeling framework.

Indices and tables
------------------
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. _`birdhouse`: https://birdhouse.readthedocs.io
.. _`birdy`: https://birdy.readthedocs.io
.. _`Raven`: https://raven.uwaterloo.ca
.. _`CANARIE`: https://www.canarie.ca
