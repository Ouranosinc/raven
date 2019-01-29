#!/bin/bash

set -e

echo "saving input files for the best solution found ..."

if [ ! -e best ] ; then
    mkdir best
fi

cp model/raven-gr4j-cemaneige.rvp  best/raven-gr4j-cemaneige.rvp
cp model/raven-gr4j-cemaneige.rvc  best/raven-gr4j-cemaneige.rvc
cp model/output/Diagnostics.csv best/. 
cp model/output/Hydrographs.*   best/.

exit 0

