#!/bin/bash

set -e

echo "saving input files for the best solution found ..."

if [ ! -e best ] ; then
    mkdir best
fi

cp model/raven-mohyse-salmon.rvp  best/raven-mohyse-salmon.rvp
cp model/raven-mohyse-salmon.rvh  best/raven-mohyse-salmon.rvh
cp model/output/Diagnostics.csv best/Diagnostics.csv
cp model/output/Hydrographs.csv best/Hydrographs.csv

exit 0
