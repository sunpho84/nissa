#!/bin/bash

#base path are already set

#now we are in the folder of the analized conf
#it's time to load info about the analysis
source $base_scripts/Setup/analysis.sh

#setup parallelization
source $base_scripts/Setup/mpi.sh

#load the info about the conf, basically nothing
source $base_scripts/Setup/gaugeconf.sh


