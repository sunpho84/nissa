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

#volume
L=[L]
T=[T]

#noise type residual for the inverter and maximal number of iterations
noise_type=4 #-1 1 2 4
stopping_residue=[resd]
num_max_iter=[num]
minimal_residue=[min]

#list of masses
kappa=[kappac]
list_mu=( [mu1] [mu2] )

#List of two points functions
two_points_correlations=(P5P5 A0P5 P5A0)

############## - 3 Setups ################

#now cd to the analysis
cd $base_analysis

#search the gauge configuration to analyse, and setup it
#exit if everything calculated
source $base_scripts/select_new_conf.sh

#if not present, generate the seed for the conf
if [ ! -f $base_conf/global_seeds.sh ]
then
    echo seed=$(bash $base_scripts/casual.sh) > $base_conf/global_seeds.sh
    echo twall=$(($(bash $base_scripts/casual.sh)%$T)) >> $base_conf/global_seeds.sh
fi

source $base_conf/global_seeds.sh

################# main script #################

#take the list of the micro-correlations

cd $base_nissa/Data/Correlations_content

list_2micro=$(cat ${two_points_correlations[@]}|awk '{print $1,$2}')

cd $base_conf

#prepare the input file

echo "\
L $L
T $T
GaugeConfPath Conf
Kappa $kappa
Seed $seed
TLeftWall $twall
TRightWall "$(( ( $twall + $T/2 ) % $T ))"
NoiseType $noise_type
NMass ${#list_mu[@]}
${list_mu[@]}
Residue $stopping_residue
StoppingCriterion standard
MinimalResidue $minimal_residue
NiterMax $num_max_iter

OutFileOtto otto

OutFileTwoPoints two_points
NContrTwoPoints "$(echo "$list_2micro"|wc -l)"
${list_2micro[@]}
" > input

#shell the program
$MPI_TM_PREF $base_nissa/Appretto/projects/eight_BK/eight_BK_all_in_one input

######################### finalization ###########################

#prohibite future re-run of the same analysis on this conf
touch $base_conf/analysis_"$analysis_name"_completed

cd $base_analysis

llsubmit analysis.sh

echo
echo "====> JOB ENDED AT : " $(date)
echo

