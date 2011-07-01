#!/bin/bash

bash compile_gpu.sh

for i in choose copy
do
  nvcc -o test/$i test/$i.c gpu.o -L/usr/lib/openmpi/lib -lmpi -I/usr/lib/openmpi/include
done