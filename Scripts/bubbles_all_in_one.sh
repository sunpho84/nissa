#!/bin/bash

#set the base folder
source ~/nissa_conf.sh

echo
echo " ##################################################################################"
echo " ##################################################################################"
echo " ####                                                                          ####"
echo " ####                   Disconnected diagrams calculation                      ####"
echo " ####                                                                          ####"
echo " ##################################################################################"
echo " ##################################################################################"
echo
echo "Started at : " $(date)
echo
echo "Working dir: $PWD"
echo

#take initial time
tic=$(date +%s)

#reset the time log
if [ -f "$base_conf/time_log" ]
then
    mv -f $base_conf/time_log $base_conf/time_log_$tic
fi

echo "Bubbles calculation started at : "$(date) >> $base_conf/time_log

#count the number of mu
nmu=${#list_mu[@]}

#Generate the wall4 source
(
    echo "L "$L
    echo "T "$T
    echo "Ncontr 16"
    echo $(seq 0 15)
    echo "NChromoContr 16"
    echo 0 5
    echo "GaugeConf Conf"
    echo "Seed "$source_seed
    echo "StartingSource "$starting_source
    echo "EndingSource "$ending_source
    echo "NoiseType "$source_noise
    echo "Kappa "$kappa
    echo "NMass "$nmu
    echo ${mu[@]}
    echo "Residue "$inversion_precision
    echo "StoppingCriterion standard"
    echo "MinimalResidue "$minimal_precision
    echo "NIterMax 10000"
    echo "Output bubble"
)  > input

$MPI_TM_PREF $base_nissa/Appretto/projects/bubbles/bubbles_standalone input
	  
echo "Bubbles calculation ended at : "$(date) | tee >> $base_conf/time_log
