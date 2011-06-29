#!/bin/bash

nvcc -c -o gpu.o gpu.cu

nvcc -o test test.c gpu.o