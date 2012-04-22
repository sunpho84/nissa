#!/bin/bash

#get the subversion
SVN_VERS=$(svnversion -n 2>/dev/null)

if [ ! -f params.sh ]
then
    echo -e "COMP=\nLEMON_PATH=\nINC_PATH=\"-I$PWD/src\"\nLIB_PATH=\nCOMP_FLAG=\"-std=c99 -O2 -Wall\"\n" > params.sh
    echo "Error, file 'params.sh' not found. Created: fill it and relaunch."
    exit
fi

source params.sh
if [ "$COMP" == "" ];then echo "File 'params.sh' not filled!";exit;fi

if [ "$1" == "" ]
then
    echo "Error: use "$0" program_to_compile"
    exit
fi

comp()
{
    $COMP $3 -o $1 $2 $COMP_FLAG $INC_PATH -lm -llemon -L$LEMON_PATH/lib/ -I$LEMON_PATH/include/ -Isrc -D SVN_VERS=\"$SVN_VERS\"
}

comp $1 $1.cpp
