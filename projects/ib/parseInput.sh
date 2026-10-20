#!/bin/bash

if [ -z $2 ]
then
    echo "Use: $0 input output"
    exit
fi

awk -f ${0/sh/awk} $1 | dot -Tpng > $2
