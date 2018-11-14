#!/bin/bash

if [ -d "src" ]; then
   \rm -rf src/
fi
mkdir src
cd src
wget http://www.eng.buffalo.edu/~lsmatott/Ostrich/Ostrich_v17.12.19.zip 
unzip Ostrich_v17.12.19.zip
cd ..

