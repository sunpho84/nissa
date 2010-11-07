#!/bin/bash

#set the base folder
source ~/nissa_conf.sh

rm -vf instruct_*

#This reads all the user needed correlation functions
corr_list=$(cat correlations_needed)

#Path where to find list of micro-correlations needed
path_instruct=$(for corr in $corr_list;do echo $base_nissa/Data/Correlations_content/$corr;done)

cat $path_instruct|gawk '{print $1,$2}'|sort|uniq > micro_correlations

cat micro_correlations|while read a b
do
    echo "    <operator>" $a"</operator>"
    echo "    <operator>" $b"</operator>"
done > contract_two_lines_middle.xml