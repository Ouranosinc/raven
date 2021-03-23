from raven.processes import RavenProcess
from ravenpy.models import HMETS, GR4JCN, HBVEC, MOHYSE, RavenMultiModel
from . import wpsio as wio


class RavenMultiModelProcess(RavenProcess):
    identifier = 'raven-multi-model'
    abstract = 'Multi model simulation'
    title = ''
    version = ''
    tuple_inputs = {'hmets': HMETS.params,
                    'gr4jcn': GR4JCN.params,
                    'hbvec': HBVEC.params,
                    'mohyse': MOHYSE.params}

    inputs = [wio.ts, wio.hmets, wio.gr4jcn, wio.hbvec, wio.start_date, wio.end_date, wio.duration, wio.run_name,
              wio.name, wio.area, wio.latitude, wio.longitude, wio.elevation]
    model_cls = RavenMultiModel

    def model(self, request):
        solution = self.get_config(request, ids=("rvc",))
        models = list(set(request.inputs.keys()).intersection(self.tuple_inputs.keys()))
        model = self.model_cls(models=models, workdir=self.workdir)
        return model
