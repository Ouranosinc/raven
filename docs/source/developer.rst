Developer Documentation
=======================


Running Raven from Python
-------------------------

The RAVEN hydrological model is designed to be run from the terminal; it's not available as a routine in a library. To use it within Python means that we actually need to create configuration files on disk and then launch the RAVEN model in that directory. There is a considerable amount of boilerplate code involved that we've hidden in the :class:`Raven` class and subclasses for each emulated model.


Here is an example of how the GR4J emulator would be called using input test data.

.. code-block:: python

   from raven.models import GR4JCemaneige
   from tests.common import TESTDATA

   ts = TESTDATA['raven-gr4j-cemaneige-nc-ts']
   gr4j = GR4JCemaneige(workdir='/tmp/test')
   params = gr4j.RVP.params(0.529, -3.396, 407.29, 1.072, 16.9, 0.947)
   gr4j([ts,],  rvp={'params':params}, rvi={'start_date':dt.datetime(2000, 1, 1), 'end_date':dt.datetime(2002, 1, 1)}, rvh={'area':4100})
   gr4j.diagnostics['DIAG_RMSE']
   >>> 39.701
