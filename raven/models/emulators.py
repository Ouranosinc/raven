from raven.models import Raven, Ostrich
from pathlib import Path
from collections import namedtuple
from .rv import RV, RVI, Ost


class GR4JCN(Raven):
    """GR4J + Cemaneige"""
    templates = tuple((Path(__file__).parent / 'raven-gr4j-cemaneige').glob("*.rv?"))

    params = namedtuple('GR4JParams', ('GR4J_X1', 'GR4J_X2', 'GR4J_X3', 'GR4J_X4', 'CEMANEIGE_X1', 'CEMANEIGE_X2'))

    rvp = RV(params=params(None, None, None, None, None, None))
    rvt = RV(pr=None, prsn=None, tasmin=None, tasmax=None, evspsbl=None, water_volume_transport_in_river_channel=None)
    rvi = RVI()
    rvh = RV(name=None, area=None, elevation=None, latitude=None, longitude=None)
    rvd = RV(one_minus_CEMANEIGE_X2=None, GR4J_X1_hlf=None)

    def derived_parameters(self):
        self.rvd.GR4J_X1_hlf = self.rvp.params.GR4J_X1 * 1000. / 2.
        self.rvd.one_minus_CEMANEIGE_X2 = 1.0 - self.rvp.params.CEMANEIGE_X2


class GR4JCN_OST(Ostrich, GR4JCN):
    _p = Path(__file__).parent / 'ostrich-gr4j-cemaneige'
    templates = tuple(_p.glob("model/*.rv?")) + tuple(_p.glob('*.t??'))

    low = GR4JCN.params
    high = GR4JCN.params
    txt = RV(max_iterations=None,
             low=low(None, None, None, None, None, None),
             high=high(None, None, None, None, None, None),
             )

    def derived_parameters(self):
        """Derived parameters are computed by Ostrich."""
        pass


class MOHYSE(Raven):
    templates = tuple((Path(__file__).parent / 'raven-mohyse').glob("*.rv?"))

    params = namedtuple('MOHYSEParams', ', '.join(['par_x{:02}'.format(i) for i in range(1, 9)]))
    hrus = namedtuple('MOHYSEHRU', ('par_x09', 'par_x10'))

    rvp = RV(params=params(*((None, ) * 8)))
    rvh = RV(name=None, area=None, elevation=None, latitude=None, longitude=None, hrus=hrus(None, None))
    rvt = RV(pr=None, prsn=None, tasmin=None, tasmax=None, evspsbl=None, water_volume_transport_in_river_channel=None)
    rvi = RVI()
    rvd = RV(par_rezi_x10=None)

    def derived_parameters(self):
        self.rvd['par_rezi_x10'] = 1.0 / self.rvh.hrus.par_x10


class HMETS(GR4JCN):
    templates = tuple((Path(__file__).parent / 'raven-hmets').glob("*.rv?"))

    params = namedtuple('HMETSParams', ('GAMMA_SHAPE', 'GAMMA_SCALE', 'GAMMA_SHAPE2', 'GAMMA_SCALE2',
                                        'MIN_MELT_FACTOR', 'MAX_MELT_FACTOR', 'DD_MELT_TEMP', 'DD_AGGRADATION',
                                        'SNOW_SWI_MIN', 'SNOW_SWI_MAX', 'SWI_REDUCT_COEFF', 'DD_REFREEZE_TEMP',
                                        'REFREEZE_FACTOR', 'REFREEZE_EXP', 'PET_CORRECTION',
                                        'HMETS_RUNOFF_COEFF', 'PERC_COEFF', 'BASEFLOW_COEFF_1',
                                        'BASEFLOW_COEFF_2', 'TOPSOIL', 'PHREATIC'))

    rvp = RV(params=params(*((None,) * len(params._fields))))
    rvt = RV(pr=None, prsn=None, tasmin=None, tasmax=None, evspsbl=None, water_volume_transport_in_river_channel=None)
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


class HBVEC(GR4JCN):
    templates = tuple((Path(__file__).parent / 'raven-hbv-ec').glob("*.rv?"))

    params = namedtuple('HBVECParams', ('par_x{:02}'.format(i) for i in range(1, 22)))
    mae = namedtuple('MeanAverageEvap', ('x{:02}'.format(i) for i in range(1, 13)))
    mat = namedtuple('MeanAverageTemp', ('x{:02}'.format(i) for i in range(1, 13)))

    rvp = RV(params=params(*((None,) * len(params._fields))))
    rvd = RV(one_plus_par_x15=None, par_x11_half=None, mae=mae, mat=mat)
    rvt = RV(pr=None, prsn=None, tasmin=None, tasmax=None, evspsbl=None,
             water_volume_transport_in_river_channel=None)
    rvh = RV(name=None, area=None, elevation=None, latitude=None, longitude=None)

    def derived_parameters(self):
        import xarray as xr

        self.rvd['one_plus_par_x15'] = self.rvp.params.par_x15 + 1.0
        self.rvd['par_x11_half'] = self.rvp.params.par_x11 / 2.0

        tasmax = xr.open_dataset(self.rvt.tasmax)[self.rvd.tasmax_var]
        tasmin = xr.open_dataset(self.rvt.tasmin)[self.rvd.tasmin_var]
        evap = xr.open_dataset(self.rvt.evspsbl)[self.rvd.evspsbl_var]

        tas = (tasmax + tasmin) / 2.
        self.rvd.mat = self.mat(*tas.groupby('time.month').mean().values)
        self.rvd.mae = self.mae(*evap.groupby('time.month').mean().values)
