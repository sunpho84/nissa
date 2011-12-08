#!/bin/bash

if [ -f $1.c ]
then
    a=$1.c
else
    if [ -f $1.cpp ]
    then
	a=$1.cpp
    else
	echo "Error, neither file $1.c and $1.cpp present!"
	exit
    fi
fi

g++ -o $1 $a -I /Users/francesco/Prace/nissa/Analysis_programs/src $suff -Wall `rootlib` -I $(dirname $1) -O2
