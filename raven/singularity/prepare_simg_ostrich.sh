#!/bin/bash

\rm -f ostrich.img
singularity create ostrich.img
sudo singularity bootstrap ostrich.img ubuntu_ostrich.def
# Expecting already logged to the singularity registry
export SREGISTRY_CLIENT=registry
sregistry push --name hydro/ostrich --tag latest ostrich.img


