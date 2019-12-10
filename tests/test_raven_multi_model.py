from .common import TESTDATA
from raven.models import RavenMultiModel
import datetime as dt
import zipfile


class TestRavenMultiModel:

    def test_simple(self):
        ts = TESTDATA['raven-hmets-nc-ts']
        model = RavenMultiModel(models=['gr4jcn', 'hmets'])
        gr4jcn = (0.529, -3.396, 407.29, 1.072, 16.9, 0.947)
        hmets = (9.5019, 0.2774, 6.3942, 0.6884, 1.2875, 5.4134, 2.3641, 0.0973, 0.0464, 0.1998, 0.0222, -1.0919,
                 2.6851, 0.3740, 1.0000, 0.4739, 0.0114, 0.0243, 0.0069, 310.7211, 916.1947)

        model(ts,
              start_date=dt.datetime(2000, 1, 1),
              end_date=dt.datetime(2002, 1, 1),
              area=4250.6,
              elevation=843.0,
              latitude=54.4848,
              longitude=-123.3659,
              gr4jcn=gr4jcn,
              hmets=hmets,
              )

        assert len(model.q_sim) == 2
        z = zipfile.ZipFile(model.outputs['rv_config'])
        assert len(z.filelist) == 10
