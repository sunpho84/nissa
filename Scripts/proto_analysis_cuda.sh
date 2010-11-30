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

#number of source, precision of the inversion
nsources=[number of sources]
source_noise=[-1 1 2 4]
inversion_precision=[inverter_precision]

#mu
list_mu=( [mu1] [mu2] )
beta=[beta]
kappa=[kappac]

#perform all the setting up
source $base_scripts/setup.sh

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
source $base_scripts/bubbles.sh

#prohibite future re-run of the same analysis on this conf
touch $base_conf/analysis_"$analysis_name"_completed

cd $base_analysis
#llsubmit analysis.sh
