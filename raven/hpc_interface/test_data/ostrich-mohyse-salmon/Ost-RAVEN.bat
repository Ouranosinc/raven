@echo off
copy .\raven-mohyse-salmon.rvp model\raven-mohyse-salmon.rvp
copy .\raven-mohyse-salmon.rvh model\raven-mohyse-salmon.rvh
cd model

REM run the raven executable, which creates the diagnostics file
Raven_rev153_MacOS.exe raven-mohyse-salmon -o output\

cd ..
