#!/bin/bash

set -e

cp ./*.rvp model/
cp ./*.rvh model/

cd model

/opt/Raven_rev157/raven_rev.exe `basename *.rvp .rvp` -o output/
#opt/en_rev153_MacOS.exe raven-mohyse-salmon -o output/

cd ..

exit 0

