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

#noise type residual for the inverter and maximal number of iterations
noise_type=4 #-1 1 2 4
stopping_residue=[resd]
minimal_residue=[resd2]
num_max_iter=[num]

#list of theta and sea for which we will do first inversion
kappa=[kappac]
list_mu=( [mu1] [mu2] )
list_theta=( [theta1] [theta2] ) #multiply by pi/L (periodic = 2)

#List of two points functions (and with chromo insertion)
two_points_correlations=(A0A0 A0P5 A0S0 AKAK BKBK P5A0 P5P5 P5S0 P5V0 S0A0 S0S0 S0V0 S0P5 TKTK V0A0 V0P5 V0S0 V0V0 VKAK VKVK)
two_points_ch_correlations=(S0P5 P5P5)

#List of info for the three points
three_points_correlations=(A0P5 AKP5 BKP5 P5P5 S0P5 TKP5 V0P5 VKP5)
three_points_ch_correlations=(S0P5 P5P5)
itheta_mu_r_spec=(0 0 0) #specify the index of theta, mu and r to use for the spectator

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
list_3micro=$(cat ${three_points_correlations[@]}|awk '{print $1,$2}')

list_2ch_micro=$(cat ${two_points_ch_correlations[@]}|awk '{print $1,$2}')
list_3ch_micro=$(cat ${three_points_ch_correlations[@]}|awk '{print $1,$2}')

cd $base_conf

#prepare the input file

echo "\
L $L
T $T
GaugeConfPath Conf
Kappa $kappa
Seed $seed
Twall $twall
NoiseType $noise_type
Nmass ${#list_mu[@]}
${list_mu[@]}
NTheta ${#list_theta[@]}
${list_theta[@]}
Residue $stopping_residue
StoppingCriterion standard
MinimalResidue $minimal_residue
NiterMax $num_max_iter

NContrTwoPoints "$(echo "$list_2micro"|wc -l)"
${list_2micro[@]}
NChromoContrTwoPoints "$(echo "$list_2ch_micro"|wc -l)"
${list_2ch_micro[@]}
OutfileTwoPoints two_points

NSpec 1
iThetaMassR ${itheta_mu_r_spec[@]}

NContrThreePoints "$(echo "$list_3micro"|wc -l)"
${list_3micro[@]}
NChromoContrThreePoints "$(echo "$list_3ch_micro"|wc -l)"
${list_3ch_micro[@]}
OutfileThreePoints three_points" > input

#shell the program
$MPI_TM_PREF $base_nissa/Appretto/projects/semileptonic/semileptonic input

######################### finalization ###########################

#prohibite future re-run of the same analysis on this conf
touch $base_conf/analysis_"$analysis_name"_completed

cd $base_analysis

llsubmit analysis.sh

echo
echo "====> JOB ENDED AT : " $(date)
echo

