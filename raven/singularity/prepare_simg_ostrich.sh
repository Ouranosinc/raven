#!/bin/bash

\rm -f ubuntu.img
singularity create ubuntu.img
sudo singularity bootstrap ubuntu.img ubuntu_ostrich.def
# Expecting already logged to the singularity registry
export SREGISTRY_CLIENT=registry
sregistry push --name hydro/ostrich --tag latest ubuntu.img


