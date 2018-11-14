#!/bin/bash

\rm -f ubuntu.img
singularity create ubuntu.img
sudo singularity bootstrap ubuntu.img ubuntu.def
# Expecting already logged to the singularity registry
export SREGISTRY_CLIENT=registry
sregistry push --name hydro/raven --tag latest ubuntu.img


