# @ shell = /bin/bash
# @ job_name = job
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

#setup everything
source ~/nissa_conf.sh

############ - 1 - Paths and names #########

#set base path for the analysis
base_analysis=[path_to_the_analysis]

#file with the list of configurations
conf_list_file=[conf_list_file]

#name of the analysis (eg: "2pts", or "Bk")
analysis_name=[name]

#Path to the configurations (which have to have named as 'conf.xxxx'
source_confs=[path_to_confs]

#precision for the I/O (32 or 64 bits)
IO_prec=[32 or 64]

########### - 2 - Physical information ########

#volume and parallelization (xyz)
L=[L]
T=[T]
NProc=([X] [Y] [Z])

#Different source to be generated, in different columns
list_source_name=( [name] )
list_source_type=( [type] ) # Point12  Wall4  Wall1
list_source_pars=( [posi] ) # t_x_y_z  t      t
list_source_nois=( [type] ) # possible: -1 1 2 4
list_source_seed=( [seed] ) # seed
list_source_prec=( [resd] ) # residual

#list of theta and sea for which we will do first inversion
list_theta=( [theta1] [theta2] )
list_mu=( [mu1] [mu2] )
beta=[beta]
kappa=[kappac]
musea=[musea]

#List of two points functions
two_points_correlations=(A0A0 A0P5 A0S0 AKAK BKBK P5A0 P5P5 P5S0 P5V0 S0A0 S0S0 S0V0 S0P5 TKTK V0A0 V0P5 V0S0 V0V0 VKAK VKVK)
two_points_theta1=(0) #in the case you want more theta for the spectator in the 2pts contractions

#List of info for the three points
three_points_correlations=(A0P5 AKP5 BKP5 P5P5 S0P5 TKP5 V0P5 VKP5)
list_seq_theta=(0)    #specify the theta values for the sequential propagator
list_itheta_spec=(0)  #specify the theta to use for the spectator
list_imu_spec=(0)     #specify the mu to use for the spectator
list_r_spec=(0)       #specify the r to use for the spectator

############## - 3 Setups ################

nsource=${#list_source_type[@]}

#now cd to the analysis
cd $base_analysis

#search the gauge configuration to analyse, and setup it
#exit if everything calculated
source $base_scripts/select_new_conf.sh

#if not present, generate the seed for the conf
if [ ! -f $base_conf/global_seeds.sh ]
then
    echo additive_seed=$(bash $base_scripts/casual.sh) > $base_conf/global_seeds.sh
    echo additive_time_offset=$(($(bash $base_scripts/casual.sh)%$T)) >> $base_conf/global_seeds.sh
fi

source $base_conf/global_seeds.sh

#list of the program to run
#comment unwanted things
source $base_scripts/two_points.sh
source $base_scripts/three_points.sh

#prohibite future re-run of the same analysis on this conf
touch $base_conf/analysis_"$analysis_name"_completed

cd $base_analysis
llsubmit analysis.sh
