#!/bin/bash

if [ "$configured" == "" ]
then
    
    echo ------loading base parameters-------
    
    #set path for ahmidas
    base_ahmidas_pa=$base_programs/ahmidas_pa
    base_ahmidas_fe=$base_programs/ahmidas_fe
    
    #set path for tmLQCD
    base_tmLQCD_pa=$base_programs/tmLQCD_pa
    base_tmLQCD_fe=$base_programs/tmLQCD_fe
    
    #set path for scripts
    base_scripts=$base_nissa/Scripts
    
    #set path for input prototypes
    base_protos=$base_nissa/Protos
    
    #if the job have been sent to a queque, cd to the folder
    #nb: all the job should be called on the folder containing the 
    #analysis for the current configuration
    if [ "$PBS_O_WORKDIR" != "" ]
    then
	base_ahmidas=$base_ahmidas_pa
	base_tmLQCD=$base_tmLQCD_pa
	cd $PBS_O_WORKDIR
	echo "We are in a job queque"
    else
	base_ahmidas=$base_ahmidas_fe
	base_tmLQCD=$base_tmLQCD_fe
	echo "Not in a job queque"
    fi
    
    echo

    echo -e Base Programs path:\\t $base_programs
    echo -e Ahmidas location:\\t $base_ahmidas
    echo -e tmLQCD location:\\t $base_tmLQCD
    echo -e Main script path:\\t $base_scripts
    echo -e Program prototypes:\\t $base_protos

    configured=1
fi