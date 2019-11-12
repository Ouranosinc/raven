from . base import Raven
from . emulators import GR4JCN, HBVEC, HMETS, MOHYSE, get_model
from . rv import RV, RVI


class RavenMultiModel(Raven):
    identifier = 'raven-multi-model'

    rvt = RV(pr=None, prsn=None, tasmin=None, tasmax=None, evspsbl=None, water_volume_transport_in_river_channel=None)
    rvi = RVI()
    rvh = RV(name=None, area=None, elevation=None, latitude=None, longitude=None)

    def __init__(self, models, workdir=None):
        """Create multi-model raven instance.

        Parameters
        ----------
        models : sequence
          Model identifiers ('gr4jcn', 'hmets', 'mohyse', 'hbvec').
        """
        import tempfile

        self._names = models
        self._models = []

        workdir = workdir or tempfile.mkdtemp()
        Raven.__init__(self, workdir)

        for name in models:
            m = get_model(name)(workdir)
            m.model_dir = m.name
            self._models.append(m)

    def _rename_run_name(self, run_name=None):
        rns = set([m.rvi.run_name for m in self._models])
        if (run_name is not None) or (len(rns) < len(self._models)):
            for m in self._models:
                rn = run_name or m.rvi.run_name
                m.rvi.run_name = rn + '_' + m.identifier

    def assign(self, key, value):
        """Assign key to all models, unless it's model parameters."""
        if key in self._names:
            m = self._models[self._names.index(key)]
            m.assign('params', value)
        else:
            for m in self._models:
                m.assign(key, value)

    @property
    def rvs(self):
        out = []
        for m in self._models:
            out.extend(m.rvs)
        return out

    def run(self, ts, overwrite=False, **kwds):
        """Run model.

        Parameters
        ----------
        kwds : dict
          model_name : array
            Parameter array.
        """
        if overwrite:
            self.setup(overwrite)

        self._rename_run_name(kwds.pop('run_name', None))

        p = {}
        for m in self._models:
            p[m.identifier] = kwds.pop(m.identifier, None)

        procs = []
        for m in self._models:
            # Add params to kwds if passed in run.
            kw = kwds.copy()
            if p[m.identifier]:
                kw['params'] = p[m.identifier]

            procs.extend(m.run(ts, **kw))

        return procs
