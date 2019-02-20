from raven.processes import RavenProcess
from pywps import LiteralInput
from raven.models import get_model, HMETS, GR4JCN, HBVEC, MOHYSE
from . import wpsio as wio

hmets = LiteralInput('hmets', 'Comma separated list of HMETS parameters',
                     abstract='Parameters: ' + ', '.join(HMETS.params._fields),
                     data_type='string',
                     min_occurs=0)

gr4jcn = LiteralInput('gr4jcn', 'Comma separated list of GR4JCN parameters',
                      abstract='Parameters: ' + ', '.join(GR4JCN.params._fields),
                      data_type='string',
                      min_occurs=0)


class RavenMultiModelProcess(RavenProcess):
    identifier = 'raven-multi-model'
    abstract = 'Multi model simulation'
    title = ''
    version = ''
    tuple_inputs = {'hmets': HMETS.params,
                    'gr4jcn': GR4JCN.params,
                    'hbvec': HBVEC.params,
                    'mohyse': MOHYSE.params}

    inputs = [wio.ts, hmets, gr4jcn, wio.start_date, wio.end_date, wio.duration, wio.run_name,
              wio.name, wio.area, wio.latitude, wio.longitude, wio.elevation]

    def __init__(self):

        super(RavenProcess, self).__init__(
            self._handler,
            identifier=self.identifier,
            title=self.title,
            version=self.version,
            abstract=self.abstract,
            inputs=self.inputs,
            outputs=self.outputs,
            status_supported=True,
            store_supported=True
        )

    def _handler(self, request, response):

        # Extract the individual model parameters
        params = {}
        for key, val in self.tuple_inputs.items():
            if key in request.inputs:
                params[key] = request.inputs.pop(key)

        # Launch the individual models
        resp = {}
        for key, val in params.items():
            req = request.copy()
            req['params'] = val
            self.model_cls = get_model(key)
            resp[key] = RavenProcess._handler(self, req, response)

        # Aggregate the results
