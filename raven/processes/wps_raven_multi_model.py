from raven.processes import RavenProcess
from pywps import LiteralInput
from raven.models import HMETS, GR4JCN, HBVEC, MOHYSE, RavenMultiModel
from . import wpsio as wio

hmets = LiteralInput('hmets', 'Comma separated list of HMETS parameters',
                     abstract='Parameters: ' + ', '.join(HMETS.params._fields),
                     data_type='string',
                     min_occurs=0)

gr4jcn = LiteralInput('gr4jcn', 'Comma separated list of GR4JCN parameters',
                      abstract='Parameters: ' + ', '.join(GR4JCN.params._fields),
                      data_type='string',
                      min_occurs=0)

# This won't work until we merge the hrus with the parameters
mohyse = LiteralInput('mohyse', 'Comma separated list of MOHYSE parameters',
                      abstract='Parameters: ' + ', '.join(MOHYSE.params._fields),
                      data_type='string',
                      min_occurs=0)

hbvec = LiteralInput('hbvec', 'Comma separated list of HBV-EC parameters',
                     abstract='Parameters: ' + ', '.join(HBVEC.params._fields),
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

    inputs = [wio.ts, hmets, gr4jcn, hbvec, wio.start_date, wio.end_date, wio.duration, wio.run_name,
              wio.name, wio.area, wio.latitude, wio.longitude, wio.elevation]
    model_cls = RavenMultiModel

    def model(self, request):
        models = list(set(request.inputs.keys()).intersection(self.tuple_inputs.keys()))
        return self.model_cls(models=models, workdir=self.workdir)
