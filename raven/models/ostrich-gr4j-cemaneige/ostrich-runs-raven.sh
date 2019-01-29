#!/bin/bash

set -e

cp ./raven-gr4j-cemaneige.rvp model/raven-gr4j-cemaneige.rvp
cp ./raven-gr4j-cemaneige.rvc model/raven-gr4j-cemaneige.rvc

cd model

./raven raven-gr4j-cemaneige -o output/

cd ..

exit 0

