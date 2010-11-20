# @ shell = /bin/bash
# @ job_name = job_strunz
# @ error  = $(job_name).$(jobid).err
# @ output = $(job_name).$(jobid).out
# @ environment = COPY_ALL
# @ notification = always
# @ notify_user = Francesco.Sanfilippo@roma1.infn.it
# @ wall_clock_limit = 01:00:00
# @ job_type         = BLUEGENE
# @ bg_size          = 512
# @ bg_connection    = TORUS
# @ queue
#

#!/bin/bash

#setup everythig
source ~/nissa_conf.sh

#set base path for the analysis
base_analysis=[path_to_the_analysis]

L=[L]
T=[T]

#Path to the configurations (which have to have named as 'conf.xxxx'
source_confs=[path_to_confs]

#Different source to be generated, in different columns
list_source_name=( [name] )
list_source_type=( [type] ) # Point12  Wall4  Wall1
list_source_pars=( [posi] ) # t_x_y_z  t      t
list_source_nois=( [type] ) # possible: -1 1 2 4
list_source_seed=( [seed] ) # seed
list_source_prec=( [resd] ) # residual

#this is the list of theta and sea for which we will do first inversion
list_theta=( [theta1] [theta2] )
list_mu=( [mu1] [mu2] )
beta=[beta]
kappa=[kappac]
musea=[musea]

#List of two points functions
two_points_correlations=(A0A0 A0P5 AKAK P5A0 P5P5 P5V0 S0S0 S0P5 TKAK TKTK V0A0 V0P5 V0V0 VKAK VKVK)
two_points_theta1=(0) #in the case you want more theta for the spectator

#This is for the three points
three_points_correlations=(A0A0 A0P5 AKAK P5A0 P5P5 P5V0 S0S0 S0P5 TKAK TKTK V0A0 V0P5 V0V0 VKAK VKVK)
list_itheta_spec=(0 1)
list_imu_spec=(0 1)
list_f_spec=(0 1)

#Useless
#MPI_nn=
#MPI_nc=
#MPI_np=$(( $MPI_nn * $MPI_nc ))

#perform all the setting up
source $base_scripts/setup.sh

#now cd to the analysis
cd $base_analysis

#search the gauge configuration to analyse,
#exit if everything calculated
source $base_scripts/select_new_conf.sh

#go to the conf
base_conf=$base_analysis/$conf
mkdir -vp $base_conf
cd $base_conf

#adapt here in the case you do not want to have seed and offset different for each configs                                                                                                                      
if [ ! -f global_seeds.sh ]
then
    echo additive_seed=$(bash $base_scripts/casual.sh) > global_seeds.sh
    echo additive_time_offset=$(($(bash $base_scripts/casual.sh)%$T)) >> global_seeds.sh
fi

source global_seeds.sh

#list of the program to run
#comment unwanted things
source $base_scripts/two_points.sh
source $base_scripts/three_points.sh

cd $base_analysis
llsubmit analysis.sh