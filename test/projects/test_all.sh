#!/bin/bash

echo "Passed argument: $*"

for i in g eight_BK semileptonic_smeared
do
    cd $i
    bash test.sh $*
    cd ..
    echo
done