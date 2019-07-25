@echo off
@TITLE SAVE BEST SOLUTION

echo saving input files for the best solution found...

IF NOT EXIST best mkdir best

REM copy model\function_out.txt best\function_out.txt
copy model\output\Diagnostics.csv best\Diagnostics.csv 
copy model\raven-mohyse-salmon.rvp best\raven-mohyse-salmon.rvp
copy model\raven-mohyse-salmon.rvh best\raven-mohyse-salmon.rvh
copy model\output\Hydrographs.csv best\Hydrographs.csv
