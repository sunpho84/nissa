# @ job_name         = test_ll
# @ error            = $(job_name).$(jobid).err
# @ output           = $(job_name).$(jobid).out
# @ environment      = COPY_ALL
# @ wall_clock_limit = 00:30:00
# @ job_type         = BLUEGENE
# @ bg_size          = 64
# @ notification     = never
# @ notify_user      = fr.sanfilippo@gmail.com
# @ account_no       = INFN_RM123
# @ queue

#!/bin/bash

echo "Starting"
date

#prog and arguments
EXE=/gpfs/scratch/userexternal/fsanfili/programs/nissa/work_in_prog/bgq/test_spi
ARGS=""

#launch
touch run

inst=0

while [ -f run ]
do
    sleep 1
    if [ -f start ]
    then
        rm -fr start
	date
        echo "===========================starting============================"
	runjob --envs OMP_NUM_THREADS=64 --ranks-per-node 1 --np 64 $ARGS --exe $EXE --timeout 60 > out_$inst 2>&1
	inst=$(($inst+1))
    else
        echo "not starting"
    fi
done

