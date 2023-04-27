Getting started
===============

* Install `birdy`_ with `pip install birdhouse-birdy`
* Connect to a Raven WPS server. You may deploy your own version or use the server hosted at Ouranos_.

.. code-block:: python

    from birdy import WPSClient

    url = "https://pavics.ouranos.ca/twitcher/ows/proxy/raven"
    wps = WPSClient(url)

The :class:`wps` object behaves as a module, holding functions that, instead of being executed on your machine, call a remote process on the server. See the notebook tutorials for examples.

If you don't want to install anything and just try it, go to https://pavics.ouranos.ca/jupyter and login with the `public` account and the `public` password. Note that your notebooks will be publicly visible and modifiable, so don't leave anything valuable there. Also, from time to time we'll reset the public folders, so make sure you keep a local copy of your work.

.. _`birdy`: https://birdy.readthedocs.io
.. _`Ouranos`: https://ouranos.ca
