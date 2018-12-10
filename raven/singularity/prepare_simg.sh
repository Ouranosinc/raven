#!/bin/bash

\rm -f raven.img
singularity create raven.img
sudo singularity bootstrap raven.img ubuntu.def
# Expecting already logged to the singularity registry
export SREGISTRY_CLIENT=registry
sregistry push --name hydro/raven --tag latest raven.img


