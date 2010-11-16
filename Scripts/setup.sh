#!/bin/bash

#base path are already set

#now we are in the folder of the analized conf
#it's time to load info about the analysis
echo ------loading analysis parameters-------
source ../analysis_parameters.sh

echo -e musea:\\t $musea
echo -e kappa:\\t $kappa
echo -e beta:\\t $beta
echo -e path:\\t $base_analysis

#setup parallelization
echo ------loading mpi parameters-------

if [ "$WHERE" == JUGENE ]
then
    
    echo -e Number of cores: $LOADL_BG_SIZE
    
    MPI_TM_PREF="mpirun -np "$((4*$LOADL_BG_SIZE))" -mode VN"
    MPI_AH_PREF="mpirun -np "$LOADL_BG_SIZE" -mode SMP"
elif [ "$WHERE" == ROMA3 ]
then
    
    echo -e Number of nodes: $MPI_nn
    echo -e Number of cores: $MPI_nc
    echo -e Total instance:  $MPI_np
    
    if [ "$PBS_O_WORKDIR" != "" ]
    then
	MPI_TM_PREF="mpiexec -np "$MPI_np
	MPI_AH_PREF=""
    else
	MPI_TM_PREF="mpiexec -np 2"
	MPI_AH_PREF=""
    fi
else	
    echo "Erorr: unknown location $WHERE"
fi

#load the info about the conf, basically nothing
echo ------loading configuration parameters-------
source gaugeconf_parameters.sh

base_conf=$base_analysis/$confno
echo -e Conf\#:\\t $confno
echo -e Path: \\t $base_conf
echo
