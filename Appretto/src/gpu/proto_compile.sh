#!/bin/bash

<<<<<<< .mine
comp="nvcc gpu.cu -L/usr/lib/openmpi/lib -lmpi -I/usr/lib/openmpi/include -I/home/francesco/Prace/nissa/Appretto/src/gpu/ -ccbin /home/francesco/Prace/nissa/Appretto/gcc"
=======
comp="nvcc gpu.cu -I /usr/lib64/openmpi/1.4-gcc/include -L /usr/lib64/openmpi/1.4-gcc/lib/ -lmpi -arch sm_13"
>>>>>>> .r335

old_sha=$(cat gpu.sha 2>/dev/null)
new_sha=$($comp -E|sha1sum|tee gpu.sha)

if [ "$old_sha" != "$new_sha" ] || [ ! -f gpu.o ]
then
    rm -f gpu.o
    echo Recompiling gpu.o
    $comp -c -o gpu.o
else
    echo Do not need to recompile gpu.o
fi

