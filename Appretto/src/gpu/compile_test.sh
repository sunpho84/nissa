#!/bin/bash

bash compile_gpu.sh

for i in choose copy
do
  nvcc -o test/$i test/$i.c gpu.o -L /usr/lib64/openmpi/1.4-gcc/lib/ -lmpi -I /usr/lib64/openmpi/1.4-gcc/include
done