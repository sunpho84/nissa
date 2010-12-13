#!/bin/bash

echo 
echo " ############## Set up ################# "
echo

######################## 1 - Setup parallelization #######################

echo
if [ "$WHERE" == JUGENE ]
then
    
    echo -e Number of cores: $LOADL_BG_SIZE
    
    MPI_TM_PREF="mpirun -np "$((4*$LOADL_BG_SIZE))" -mode VN -env \"DCMF_EAGER=500000000\" -mapfile TXYZ"
    MPI_AH_PREF="mpirun -np "$((4*$LOADL_BG_SIZE))" -mode VN -env \"DCMF_EAGER=500000000\" -mapfile TXYZ"

    echo "Number of nodes and cores: "$LOADL_BG_SIZE", "$((4*$LOADL_BG_SIZE))
    
    USE_MPI=1

elif [ "$WHERE" == AURORA ]
then
    
    ENV=$(echo $AURORA_ENV | sed 's/:/ -x /g')

    echo -e Number of cores: $AURORA_NP
    
    MPI_TM_PREF="mpirun $ENV $AURORA_BTL -np $AURORA_NP -hostfile $AURORA_HOSTFILE -rf $AURORA_RANKFILE"
    MPI_AH_PREF=$MPI_TM_PREF

    USE_MPI=1

elif [ "$WHERE" == ROMA3 ]
then
    
    if [ "$PBS_O_WORKDIR" != "" ]
    then
	echo  "Total instance: 16"
    
	MPI_TM_PREF="mpirun -np 16"
	MPI_AH_PREF="mpirun -np 16"
	USE_MPI=1
    else
	MPI_TM_PREF=""
	MPI_AH_PREF=""
	USE_MPI=0
    fi

elif [ "$WHERE" == MACSILV ]
then
    echo -e Total instance:  2
    
    MPI_TM_PREF="/usr/local/bin/mpiexec -np 2"
    MPI_AH_PREF="/usr/local/bin/mpiexec -np 2"
    USE_MPI=1
    
elif [ "$WHERE" == PCFRA ]
then
    MPI_TM_PREF=""
    MPI_AH_PREF=""
    USE_MPI=0

    echo "Poor Samsung, don't stress it!"
    
else	
    echo "Erorr: unknown location $WHERE"
fi

########################### 2 - Path settings ######################

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

if [ "$USE_MPI" == "1" ]
then
    base_ahmidas=$base_ahmidas_pa
    base_tmLQCD=$base_tmLQCD_pa
    echo "We'll use MPI"
else
    base_ahmidas=$base_ahmidas_fe
    base_tmLQCD=$base_tmLQCD_fe
    echo "We'll not use MPI"
fi

########################### 3 - Information ######################

echo
echo -e Base Programs path:\\t $base_programs
echo -e Ahmidas location:\\t $base_ahmidas
echo -e tmLQCD location:\\t $base_tmLQCD
echo -e Main script path:\\t $base_scripts
echo -e Program prototypes:\\t $base_protos
