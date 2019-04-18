from xclim.streamflow import Stats, generic
from .base_xclim import make_xclim_indicator_process


stats = Stats(identifier='ts_stats',
              long_name='{freq} {op} of {indexer} daily flow ',
              description="{freq} {op} of {indexer} daily flow",
              compute=generic.select_resample_op)

TSStatsProcess = make_xclim_indicator_process('TSStats', stats)
