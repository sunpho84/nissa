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
set echo

echo
echo "====> JOB STARTED AT : " $(date)
echo

#setup everything
source ~/nissa_conf.sh

############ - 1 - Paths and names #########

#set base path for the analysis
base_analysis=[path_to_the_analysis]

#file with the list of configurations
conf_list_file=[conf_list_file]

#name of the analysis (eg: "2pts", or "Bk")
analysis_name=Bk

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
list_source_name=( Wall   )
list_source_type=( Wall4  ) # Point12  Wall4  Wall1
list_source_pars=( 0      ) # t_x_y_z  t      t
list_source_nois=( 4      ) # possible: -1 1 2 4
list_source_seed=( 0      ) # seed
list_source_prec=( [resi] ) # residual

#list of theta and sea for which we will do first inversion
list_theta=( 0 ) #multiply by pi/L (periodic = 2)
list_mu=( [mu1] [mu2] )
nmu_low=[nmu]
beta=[beta]
kappa=[kappac]
musea=[musea]

list_ops_eight=(0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15)
two_points_correlations=(P5P5 A0P5 P5A0)

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
source $base_scripts/Bk.sh

#prohibite future re-run of the same analysis on this conf
touch $base_conf/analysis_"$analysis_name"_completed

cd $base_analysis

llsubmit analysis.sh

echo
echo "====> JOB ENDED AT : " $(date)
echo

