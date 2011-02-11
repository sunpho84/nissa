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
analysis_name=[name]

#Path to the configurations (which have to have named as 'conf.xxxx'
source_confs=[path_to_confs]

########### - 2 - Physical information ########

#volume and parallelization (xyz)
L=[L]
T=[T]

#number of source, precision of the inversion
starting_source=[starting] #number of the starting sources (included)
ending_source=[ending] #number of the ending sources (excluded)
source_noise=[-1 1 2 4]
inversion_precision=[inverter_precision]
niter_max=[niter]
minimal_precision=[minimal_precision]

#mu
list_mu=( [mu1] [mu2] )
kappa=[kappa]

##################### end of tunable parameters ###################

#now cd to the analysis
cd $base_analysis

#search the gauge configuration to analyse, and setup it
#exit if everything calculated
source $base_scripts/select_new_conf.sh

#if not present, generate the seed for the conf
if [ ! -f $base_conf/global_seeds.sh ]
then
    echo source_seed=$(bash $base_scripts/casual.sh) > $base_conf/global_seeds.sh
fi

source $base_conf/global_seeds.sh

#list of the program to run
#comment unwanted things
source $base_scripts/bubbles_all_in_one.sh

#prohibite future re-run of the same analysis on this conf
touch $base_conf/analysis_"$analysis_name"_completed

cd $base_analysis

#llsubmit analysis.sh

echo
echo "====> JOB ENDED AT : " $(date)
echo

