#!/bin/bash

if [ -d "src" ]; then
   \rm -rf src/
fi
mkdir src
cd src
wget http://raven.uwaterloo.ca/files/v2.8.1/Raven_2.8.1_Source.zip
unzip Raven_2.8.1_Source.zip
cd ..
