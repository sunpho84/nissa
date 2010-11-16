#!/bin/bash

echo ------loading mpi parameters-------
#well actually these will be useless
echo Use MPI: $use_mpi

if [ "$use_mpi" != 0 ]
then
    echo -e Number of nodes: $MPI_nn
    echo -e Number of cores: $MPI_nc
    echo -e Total instance:  $MPI_np

    if [ "$WHERE" == JUGENE ]
    then
#	if [ "$PBS_O_WORKDIR" != "" ]
#	then
	    MPI_TM_PREF="mpirun -np "$((4*$LOADL_BG_SIZE))" -mode VN"
	    MPI_AH_PREF="mpirun -np "$LOADL_BG_SIZE" -mode SMP"
#	    MPI_AH_PREF="mpirun -np 32 -mode SMP"
#	else
#	    MPI_TM_PREF="llrun -np 128"
#	    MPI_AH_PREF=""
#	fi
    elif [ "$WHERE" == ROMA3 ]
    then
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
fi

echo
