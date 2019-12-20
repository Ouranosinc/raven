from collections import namedtuple
from pathlib import Path

from raven.models import Raven, Ostrich
from .rv import RV, RVT, RVI, Ost


class GR4JCN(Raven):
    """GR4J + Cemaneige"""
    identifier = 'gr4jcn'
    templates = tuple((Path(__file__).parent / 'raven-gr4j-cemaneige').glob("*.rv?"))

    params = namedtuple('GR4JParams', ('GR4J_X1', 'GR4J_X2', 'GR4J_X3', 'GR4J_X4', 'CEMANEIGE_X1', 'CEMANEIGE_X2'))

    rvp = RV(params=params(None, None, None, None, None, None))
    rvt = RVT(pr=None, prsn=None, tasmin=None, tasmax=None, evspsbl=None, water_volume_transport_in_river_channel=None)
    rvi = RVI()
    rvh = RV(name=None, area=None, elevation=None, latitude=None, longitude=None)
    rvd = RV(one_minus_CEMANEIGE_X2=None, GR4J_X1_hlf=None)

    def derived_parameters(self):
        self.rvd.GR4J_X1_hlf = self.rvp.params.GR4J_X1 * 1000. / 2.
        self.rvd.one_minus_CEMANEIGE_X2 = 1.0 - self.rvp.params.CEMANEIGE_X2


class GR4JCN_OST(Ostrich, GR4JCN):
    _p = Path(__file__).parent / 'ostrich-gr4j-cemaneige'
    templates = tuple(_p.glob("model/*.rv?")) + tuple(_p.glob('*.t??'))
    rvi = RVI(supress_output=True)
    txt = Ost(algorithm='DSS',
              max_iterations=50,
              lowerBounds=GR4JCN.params(None, None, None, None, None, None),
              upperBounds=GR4JCN.params(None, None, None, None, None, None),
              )

    def derived_parameters(self):
        """Derived parameters are computed by Ostrich."""
        pass


class MOHYSE(Raven):
    identifier = 'mohyse'
    templates = tuple((Path(__file__).parent / 'raven-mohyse').glob("*.rv?"))

    params = namedtuple('MOHYSEParams', ', '.join(['par_x{:02}'.format(i) for i in range(1, 9)]))
    hrus = namedtuple('MOHYSEHRU', ('par_x09', 'par_x10'))

    rvp = RV(params=params(*((None,) * 8)))
    rvh = RV(name=None, area=None, elevation=None, latitude=None, longitude=None, hrus=hrus(None, None))
    rvt = RVT(pr=None, prsn=None, tasmin=None, tasmax=None, evspsbl=None, water_volume_transport_in_river_channel=None)
    rvi = RVI()
    rvd = RV(par_rezi_x10=None)

    def derived_parameters(self):
        self.rvd['par_rezi_x10'] = 1.0 / self.rvh.hrus.par_x10


class MOHYSE_OST(Ostrich, MOHYSE):
    _p = Path(__file__).parent / 'ostrich-mohyse'
    templates = tuple(_p.glob("model/*.rv?")) + tuple(_p.glob('*.t??'))

    rvi = RVI(supress_output=True)
    txt = Ost(algorithm='DSS',
              max_iterations=50,
              lowerBounds=MOHYSE.params(None, None, None, None, None, None, None, None),
              upperBounds=MOHYSE.params(None, None, None, None, None, None, None, None),
              hruslowerBounds=MOHYSE.hrus(None, None),
              hrusupperBounds=MOHYSE.hrus(None, None)
              )

    def derived_parameters(self):
        """  Derived parameters are computed by Ostrich.  """
        pass


class HMETS(GR4JCN):
    identifier = 'hmets'
    templates = tuple((Path(__file__).parent / 'raven-hmets').glob("*.rv?"))

    params = namedtuple('HMETSParams', ('GAMMA_SHAPE', 'GAMMA_SCALE', 'GAMMA_SHAPE2', 'GAMMA_SCALE2',
                                        'MIN_MELT_FACTOR', 'MAX_MELT_FACTOR', 'DD_MELT_TEMP', 'DD_AGGRADATION',
                                        'SNOW_SWI_MIN', 'SNOW_SWI_MAX', 'SWI_REDUCT_COEFF', 'DD_REFREEZE_TEMP',
                                        'REFREEZE_FACTOR', 'REFREEZE_EXP', 'PET_CORRECTION',
                                        'HMETS_RUNOFF_COEFF', 'PERC_COEFF', 'BASEFLOW_COEFF_1',
                                        'BASEFLOW_COEFF_2', 'TOPSOIL', 'PHREATIC'))

    rvp = RV(params=params(*((None,) * len(params._fields))))
    rvt = RVT(pr=None, prsn=None, tasmin=None, tasmax=None, evspsbl=None, water_volume_transport_in_river_channel=None)
    rvi = RVI()
    rvd = RV(TOPSOIL_m=None, PHREATIC_m=None, SUM_MELT_FACTOR=None, SUM_SNOW_SWI=None, TOPSOIL_hlf=None,
             PHREATIC_hlf=None)

    def derived_parameters(self):
        self.rvd['TOPSOIL_hlf'] = self.rvp.params.TOPSOIL * 0.5
        self.rvd['PHREATIC_hlf'] = self.rvp.params.PHREATIC * 0.5
        self.rvd['TOPSOIL_m'] = self.rvp.params.TOPSOIL / 1000.
        self.rvd['PHREATIC_m'] = self.rvp.params.PHREATIC / 1000.
        self.rvd['SUM_MELT_FACTOR'] = self.rvp.params.MIN_MELT_FACTOR + self.rvp.params.MAX_MELT_FACTOR
        self.rvd['SUM_SNOW_SWI'] = self.rvp.params.SNOW_SWI_MIN + self.rvp.params.SNOW_SWI_MAX


class HMETS_OST(Ostrich, HMETS):
    _p = Path(__file__).parent / 'ostrich-hmets'
    templates = tuple(_p.glob("model/*.rv?")) + tuple(_p.glob('*.t??'))

    rvi = RVI(supress_output=True)
    txt = Ost(algorithm='DSS',
              max_iterations=50,
              lowerBounds=HMETS.params(None, None, None, None, None, None, None, None, None, None, None,
                                       None, None, None, None, None, None, None, None, None, None),
              upperBounds=HMETS.params(None, None, None, None, None, None, None, None, None, None, None,
                                       None, None, None, None, None, None, None, None, None, None),
              )

    def derived_parameters(self):
        """Derived parameters are computed by Ostrich."""
        pass

    def ost2raven(self, ops):
        """Return a list of parameter names calibrated by Ostrich that match Raven's parameters.

        Parameters
        ----------
        ops: dict
          Optimal parameter set returned by Ostrich.

        Returns
        -------
        HMETSParams named tuple
          Parameters expected by Raven.
        """
        names = ['par_x{:02}'.format(i) for i in range(1, 22)]
        names[5] = 'par_sum_x05_x06'
        names[9] = 'par_sum_x09_x10'

        out = [ops[n] for n in names]
        out[19] *= 1000
        out[20] *= 1000
        return self.params(*out)


class HBVEC(GR4JCN):
    identifier = 'hbvec'
    templates = tuple((Path(__file__).parent / 'raven-hbv-ec').glob("*.rv?"))

    params = namedtuple('HBVECParams', ('par_x{:02}'.format(i) for i in range(1, 22)))
    mae = namedtuple('MeanAverageEvap', ('x{:02}'.format(i) for i in range(1, 13)))
    mat = namedtuple('MeanAverageTemp', ('x{:02}'.format(i) for i in range(1, 13)))

    rvp = RV(params=params(*((None,) * len(params._fields))))
    rvd = RV(one_plus_par_x15=None, par_x11_half=None, mae=mae, mat=mat)
    rvt = RVT(pr=None, prsn=None, tasmin=None, tasmax=None, evspsbl=None,
              water_volume_transport_in_river_channel=None)
    rvh = RV(name=None, area=None, elevation=None, latitude=None, longitude=None)

    def derived_parameters(self):
        import xarray as xr

        self.rvd['one_plus_par_x15'] = self.rvp.params.par_x15 + 1.0
        self.rvd['par_x11_half'] = self.rvp.params.par_x11 / 2.0

        tasmax = xr.open_dataset(self.rvt.tasmax)[self.rvt.tasmax_var]
        tasmin = xr.open_dataset(self.rvt.tasmin)[self.rvt.tasmin_var]
        evap = xr.open_dataset(self.rvt.evspsbl)[self.rvt.evspsbl_var]

        tas = (tasmax + tasmin) / 2.
        self.rvd.mat = self.mat(*tas.groupby('time.month').mean().values)
        self.rvd.mae = self.mae(*evap.groupby('time.month').mean().values)


class HBVEC_OST(Ostrich, HBVEC):
    _p = Path(__file__).parent / 'ostrich-hbv-ec'
    templates = tuple(_p.glob("model/*.rv?")) + tuple(_p.glob('*.t??'))

    rvi = RVI(supress_output=True)
    low = HBVEC.params
    high = HBVEC.params
    txt = Ost(algorithm='DSS',
              max_iterations=50,
              lowerBounds=low(None, None, None, None, None, None, None, None, None, None, None,
                              None, None, None, None, None, None, None, None, None, None),
              upperBounds=high(None, None, None, None, None, None, None, None, None, None, None,
                               None, None, None, None, None, None, None, None, None, None),
              )

    def derived_parameters(self):
        import xarray as xr

        tasmax = xr.open_dataset(self.rvt.tasmax)[self.rvt.tasmax_var]
        tasmin = xr.open_dataset(self.rvt.tasmin)[self.rvt.tasmin_var]
        evap = xr.open_dataset(self.rvt.evspsbl)[self.rvt.evspsbl_var]

        tas = (tasmax + tasmin) / 2.
        self.rvd.mat = self.mat(*tas.groupby('time.month').mean().values)
        self.rvd.mae = self.mae(*evap.groupby('time.month').mean().values)


def get_model(name):
    """Return the corresponding Raven emulated model instance.

    Parameters
    ----------
    name : str
      Model class name or model identifier.

    Returns
    -------
    Raven model instance
    """
    from raven.models import emulators

    model_cls = getattr(emulators, name, None)

    if model_cls is None:
        for m in [GR4JCN, MOHYSE, HMETS, HBVEC]:
            if m.identifier == name:
                model_cls = m

    if model_cls is None:
        raise ValueError("Model {} is not recognized.".format(name))

    return model_cls
