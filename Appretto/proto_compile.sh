#!/bin/bash

#copy this file as "compile.sh" and put the correct paths
#launch it from the path with executable

if [ "$1" == "" ]
then
    echo "Error: use "$0" program"
    exit
fi

mpicxx -o $1 $1.cpp -I../src -Wall \
    -llemon -L/home/prace/Prace/Programs/lemon/lib/ -I/home/prace/Prace/Programs/lemon/include/ \
    -llime  -L/home/prace/Prace/Programs/lime/lib/  -I/home/prace/Prace/Programs/lime/include/