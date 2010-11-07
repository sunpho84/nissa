#!/bin/bash

#for now the configuration are not readed

if [ ! $1 ]
then
    echo "Use: "$0" script"
    exit
fi

if [ ! -f $1 ]
then
    echo $1 script not found!
    exit 1
fi

source /home/prace/prace_conf.sh

qsub -q verylong -o jobout -e joberr -l nodes=4:ppn=4  $1
