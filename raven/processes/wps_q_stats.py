from xclim.streamflow import Stats, generic, base_flow_index
from .base_xclim import make_xclim_indicator_process


stats = Stats(identifier='ts_stats',
              long_name='{freq} {op} of {indexer} daily flow ',
              description="{freq} {op} of {indexer} daily flow",
              compute=generic.select_resample_op)

freq = Stats(identifier='freq_analysis',
             long_name='N-year return period {mode} {indexer} {window}-day flow',
             description="Streamflow frequency analysis for the {mode} {indexer} {window}-day flow "
                         "estimated using the {dist} distribution.",
             compute=generic.frequency_analysis)

TSStatsProcess = make_xclim_indicator_process('TSStats', stats)

FreqAnalysisProcess = make_xclim_indicator_process('FreqAnalysis', freq)

BaseFlowIndexProcess = make_xclim_indicator_process('BaseFlowIndex', base_flow_index)
