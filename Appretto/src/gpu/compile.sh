#!/bin/bash

nvcc -c -o gpu.o gpu.cu -I /usr/lib64/openmpi/1.4-gcc/include -L /usr/lib64/openmpi/1.4-gcc/lib/ -lmpi

nvcc -o test test.c gpu.o -L /usr/lib64/openmpi/1.4-gcc/lib/ -lmpi -I /usr/lib64/openmpi/1.4-gcc/include