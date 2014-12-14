#!/bin/bash
rm -fr nissa_config_* mapping_*
runjob -n 1 --envs OMP_NUM_THREADS=1 : /fermi/home/userexternal/fsanfili/programs/nissa/build/tools/find_best_vir_partitioning/find input
mv mapping_* mapping
echo "verbosity_lv 1" >> nissa_config_*
mv nissa_config_* nissa_config