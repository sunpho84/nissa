#!/bin/bash

L=SED_L
T=SED_T

#Different source to be generated, in different columns
list_source_name=( [name] )
list_source_type=( [type] ) # Point12  Wall4  Wall1
list_source_pars=( [posi] ) # t_x_y_z  t      t
list_source_nois=( [type] ) # possible: -1 1 2 4
list_source_seed=( [seed] ) # seed
list_source_prec=( [resd] ) # residual

#some flag
use_external_additive_seed=1
use_external_additive_time_offset=1

#this is the list of theta and sea for which we will do first inversion
list_theta=( [theta1] [theta2] )
list_mu=( [mu1] [mu2] )
beta=[beta]
kappa=[kappac]
musea=[musea]

#List of two points functions
two_points_correlations=(A0A0 A0P5 AKAK P5A0 P5P5 P5V0 S0S0 S0P5 TKAK TKTK V0A0 V0P5 V0V0 VKAK VKVK)

#This is for the three points
list_itheta_spec="0 1"
list_imu_spec="0 1"
list_f_spec="0 1"

base_analysis=SED_Base_analysis

MPI_nn=
MPI_nc=
MPI_np=$(( $MPI_nn * $MPI_nc ))

#perform all the setting up
source $base_scripts/setup.sh

#list of the program to run
#comment unwanted things
source $base_scripts/two_points.sh
source $base_scripts/three_points.sh
