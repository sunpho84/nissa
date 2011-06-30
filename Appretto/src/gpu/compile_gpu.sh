#!/bin/bash

comp="nvcc gpu.cu -I /usr/lib64/openmpi/1.4-gcc/include -L /usr/lib64/openmpi/1.4-gcc/lib/ -lmpi"

old_sha=$(cat gpu.sha 2>/dev/null)
new_sha=$($comp -E|sha1sum|tee gpu.sha)

if [ "$old_sha" != "$new_sha" ]
then
    echo Recompiling gpu.o
    $comp -c -o gpu.o
else
    echo Do not need to recompile gpu.o
fi

