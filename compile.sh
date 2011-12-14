#!/bin/bash

#get the subversion
SVN_VERS=$(svnversion -n 2>/dev/null)

if [ ! -f params.sh ]
then
    echo -e "COMP=\nLEMON_PATH=\nLIME_PATH=\nINC_PATH=\"-I$PWD/src\"\nLIB_PATH=\nCOMP_FLAG=\"-std=c99 -O2 -Wall\"\n" > params.sh
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
    $COMP $3 -o $1 $2 $COMP_FLAG $INC_PATH -lm -llemon -L$LEMON_PATH/lib/ -I$LEMON_PATH/include/ -llime -L$LIME_PATH/lib/ -I$LIME_PATH/include/ -D SVN_VERS=\"$SVN_VERS\"
}

precomp()
{
    $COMP $3 -o $1 $2 $COMP_FLAG $INC_PATH -L$LEMON_PATH/lib/ -I$LEMON_PATH/include/ -L$LIME_PATH/lib/ -I$LIME_PATH/include/ -D SVN_VERS=\"$SVN_VERS\"
}

recomp_nissa()
{
    if [ -f src/nissa_checksum ];then old_checksum=$(cat src/nissa_checksum);else old_checksum="0";fi
    checksum=$(precomp - src/nissa.c "-c -E"|tee /tmp/nissa.out|sha1sum|awk '{print $1}')
    
    if [ "$old_checksum" != "$checksum" ] || [ ! -f src/nissa ]
    then
	echo "Recompiling nissa library..."
	rm -f src/nissa
	precomp src/nissa src/nissa.c -c
        
	if [ $? == 0 ] && [ -f src/nissa ]
	then
	    echo $checksum > src/nissa_checksum
	else
	    echo "Error during nissa library compilation"
	    exit
	fi
    else
	echo "nissa library not recompiled"
    fi
}

recomp_nissa
comp $1 $1.c "src/nissa"

#in case of error compile showing the output
if [ $? != 0 ]
then
    comp $1 $1.c "-E src/nissa"
fi
