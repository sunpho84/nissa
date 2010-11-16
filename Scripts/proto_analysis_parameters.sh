#!/bin/bash

L=24
T=48

list_source_name=(Wall6  Wall7  Wall8  Wall9  Wall10 Wall11 Wall12 Wall13 Wall14 Wall15 Wall16 Wall17 Wall18 Wall19)
list_source_type=(Wall4  Wall4  Wall4  Wall4  Wall4  Wall4  Wall4  Wall4  Wall4  Wall4  Wall4  Wall4  Wall4  Wall4 )
list_source_pars=(0      0      0      0      0      0      0      0      0      0      0      0      0      0     )
list_source_nois=(4      4      4      4      4      4      4      4      4      4      4      4      4      4     )
list_source_seed=(1      1      1      1      1      1      1      1      1      1      1      1      1      1     )
list_source_prec=(1.e-6  1.e-7  1.e-8  1.e-9  1.e-10 1.e-11 1.e-12 1.e-13 1.e-14 1.e-15 1.e-16 1.e-17 1.e-18 1.e-19)

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


list_theta=(0.0)
list_mu=(0.0040 0.0064 0.0085 0.0100 0.0150 0.0220 0.0270 0.0320 0.2500 0.3200 0.3900 0.4600)
#list_mu=(0.0040 0.0177 0.2123)
kappa=0.160856
beta=3.90
musea=0.0040
setup_name=3.90

list_itheta_spec="0"
list_imu_spec="0"
list_f_spec="1"

base_analysis=/work/hch02/hch02g/Analysis/3.90/24/0.0064

use_mpi=1

MPI_nn=16
MPI_nc=32
MPI_np=$(( $MPI_nn * $MPI_nc ))
