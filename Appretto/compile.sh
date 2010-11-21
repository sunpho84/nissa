#!/bin/bash

mpicxx -o $1 $1.cpp -I../src -Wall \
    -llemon -L/home/prace/Prace/Programs/lemon/lib/ -I/home/prace/Prace/Programs/lemon/include/ \
    -llime  -L/home/prace/Prace/Programs/lime/lib/  -I/home/prace/Prace/Programs/lime/include/