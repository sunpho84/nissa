#!/bin/bash

#base path are already set

#now we are in the folder of the analized conf
#it's time to load info about the analysis
echo
echo ------analysis parameters-------
echo -e musea:\\t $musea
echo -e kappa:\\t $kappa
echo -e beta:\\t $beta
echo -e path:\\t $base_analysis

#parse the source
nsource=${#list_source_name[@]}
tn0=$nsource
tn1=${#list_source_type[@]}
tn2=${#list_source_pars[@]}
tn3=${#list_source_seed[@]}
if [ "$tn1" -ne "$tn0" ] || [ "$tn2" -ne "$tn0" ] || [ "$tn3" -ne "$tn0" ]
then
    echo "Source name, type and pars do not match"
    echo "Check configuration file!"
    exit
fi

#setup parallelization
echo
echo ------mpi parameters-------
if [ "$WHERE" == JUGENE ]
then
    
    echo -e Number of cores: $LOADL_BG_SIZE
    
    MPI_TM_PREF="mpirun -np "$((4*$LOADL_BG_SIZE))" -mode VN"
    MPI_AH_PREF="mpirun -np "$LOADL_BG_SIZE" -mode SMP"

    USE_MPI=1

elif [ "$WHERE" == ROMA3 ]
then
    
    echo -e Number of nodes: $MPI_nn
    echo -e Number of cores: $MPI_nc
    echo -e Total instance:  $MPI_np
    
    if [ "$PBS_O_WORKDIR" != "" ]
    then
	MPI_TM_PREF="mpiexec -np "$MPI_np
	MPI_AH_PREF="mpiexec -np "$MPI_np
	USE_MPI=1
    else
	MPI_TM_PREF=""
	MPI_AH_PREF=""
	USE_MPI=0
    fi
else	
    echo "Erorr: unknown location $WHERE"
fi

#load the info about the conf, basically nothing
echo
echo ------configuration parameters-------
base_conf=$base_analysis/$confno
echo -e Conf\#:\\t $confno
echo -e Path: \\t $base_conf

