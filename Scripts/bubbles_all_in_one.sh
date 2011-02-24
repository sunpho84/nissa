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

#If not present, generate the source position (for nucleon)
if [ ! -f $base_conf/source_pos ]
then
    (
        echo -n $(($(bash $base_scripts/casual.sh)%$T))
        for((i=1;i<4;i++));do echo -n \ $(($(bash $base_scripts/casual.sh)%$L));done
    ) > $base_conf/source_pos
fi

#Generate the input
(
    echo "L "$L
    echo "T "$T
    echo "Ncontr 16"
    echo $(seq 0 15)
    echo "NChromoContr 2"
    echo 0 5
    echo "GaugeConf Conf"
    echo "Seed "$source_seed
    echo "SourcePosition "$(cat source_pos)
    echo "StartingSource "$starting_source
    echo "EndingSource "$ending_source
    echo "NoiseType "$source_noise
    echo "Kappa "$kappa
    echo "NMass "$nmu
    echo ${list_mu[@]}
    echo "Residue "$inversion_precision
    echo "StoppingCriterion standard"
    echo "MinimalResidue "$minimal_precision
    echo "NIterMax "$niter_max
    echo "Output bubble"
)  > input

$MPI_TM_PREF $base_nissa/Appretto/projects/bubbles/bubbles_standalone input
	  
echo "Bubbles calculation ended at : "$(date) | tee >> $base_conf/time_log
