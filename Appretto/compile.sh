#!/bin/bash

if [ ! -f paths.sh ]
then
    echo -e "COMP=\nLEMON_PATH=\nLIME_PATH=\nINC_PATH=\"-I$PWD/src\"\nLIB_PATH=\nCOMP_FLAG=\"-std=c99 -O5 -Wall\"\n" > paths.sh
    echo "Error, file 'paths.sh' not found. Created: fill it and relaunch."
    exit
fi

source paths.sh
if [ "$COMP" == "" ];then echo "File 'paths.sh' not filled!";exit;fi

if [ "$1" == "" ]
then
    echo "Error: use "$0" program_to_compile"
    exit
fi

comp()
{
    $COMP $3 -o $1 $2 $COMP_FLAG $INC_PATH -llemon -L$LEMON_PATH/lib/ -I$LEMON_PATH/include/ -llime -L$LIME_PATH/lib/ -I$LIME_PATH/include/
}

recomp_appretto()
{
    if [ -f src/appretto_checksum ];then old_checksum=$(cat src/appretto_checksum);else old_checksum="0";fi
    checksum=$(comp - src/appretto.c "-c -E"|sha1sum|awk '{print $1}')
    if [ $old_checksum != $checksum ]
    then
	echo "Recompiling Appretto library..."
	comp src/appretto src/appretto.c -c
        
	if [ $? == 0 ]
	then
	    echo $checksum > src/appretto_checksum
	else
	    echo "Error during Appretto library compilation"
	    exit
	fi
    else
	echo "Appretto library not recompiled"
    fi
}

recomp_appretto
comp $1 $1.c "src/appretto"
