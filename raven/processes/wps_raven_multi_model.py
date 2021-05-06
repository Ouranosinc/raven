from ravenpy.models import GR4JCN, HBVEC, HMETS, MOHYSE, RavenMultiModel

from raven.processes import RavenProcess

from . import wpsio as wio


class RavenMultiModelProcess(RavenProcess):
    identifier = "raven-multi-model"
    abstract = "Multi model simulation"
    title = ""
    version = ""
    tuple_inputs = {
        "hmets": HMETS.Params,
        "gr4jcn": GR4JCN.Params,
        "hbvec": HBVEC.Params,
        "mohyse": MOHYSE.Params,
    }

    inputs = [
        wio.ts,
        wio.hmets,
        wio.gr4jcn,
        wio.hbvec,
        wio.start_date,
        wio.end_date,
        wio.duration,
        wio.run_name,
        wio.area,
        wio.latitude,
        wio.longitude,
        wio.elevation,
    ]
    model_cls = RavenMultiModel

    def model(self, request):
        models = list(set(request.inputs.keys()).intersection(self.tuple_inputs.keys()))
        model = self.model_cls(models=models, workdir=self.workdir)
        return model
