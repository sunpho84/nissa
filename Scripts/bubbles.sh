#!/bin/bash

#set the base folder
source ~/nissa_conf.sh

echo
echo " ##################################################################################"
echo " ##################################################################################"
echo " ####                                                                          ####"
echo " ####                   Disconnected diagrams calculation                      ####"
echo " ####                                                                          ####"
echo " ##################################################################################"
echo " ##################################################################################"
echo
echo "Working dir: $PWD"
echo

#reset the time log
rm -f $base_conf/time_log

#take initial time
tic=$(date +%s)

#count the number of mu
nmu=${#list_mu[@]}

echo "Number of sources: "$nsource
echo

for((is=0;is<nsource;is++))
do

    echo "######################## FIRST STEP: Source generation ##########################"
    echo
    echo "Generation of "$is" source"
    echo
    
    source_seed=$base_seed+$is
    source_name=$(printf %02d $is)

    targ_dir=Sources/$source_name

    if [ -d "$targ_dir" ]
    then
        echo "Source "$source_name" already existing"
        echo
    else
        
        mkdir -pv $targ_dir
        cd $targ_dir
	
        #Generate the wall4 source
        (
	    echo "L "$L
	    echo "T "$T
	    echo "Seed "$source_seed
	    echo "TakeSlice 0"
	    echo "TWall -1"
	    echo "NoiseType "$source_nois
	    echo "Filename source"
	    echo "Precision "$IO_prec
	)  > input

        $MPI_TM_PREF $base_nissa/Appretto/tools/generate_stochastic_source/generate_stochastic_source input
        
        if [ ! "$?" ]
        then
	    echo "Source generated with success!"
#            rm input
        fi

	cd -

    fi

    #take time
    tac=$tic
    tic=$(date +%s)
    echo "Time to generate source: "$(($tic-$tac)) >> $base_conf/time_log

    echo "######################## SECOND STEP: Inversion of doublet ############################"

    for mu in ${list_mu[@]}
    do

      echo
      echo "Starting inversion for mass "$mu
      echo
      
      base_inv=$base_conf/Props/$source_name\_$mu/
      mkdir -pv $base_inv

      echo "Inverting: "$source_name" "$theta
      
      if [ ! -f $base_inv/completed ]
      then
        
          #prepare input for inverter
	  two_kappamu=$(echo $kappa|awk '{printf("%.12f\n",2*$1*'$mu')}')
	  cat $base_protos/cg_cuda.input|sed '
            s|SED_Conf|'$conf'|
            s|SED_NL|'$L'|;
            s|SED_NT|'$T'|;
            s|SED_ThetaX|0|;
            s|SED_ThetaY|0|;
            s|SED_ThetaZ|0|;
            s|SED_NrXProcs|'${NProc[0]}'|;
            s|SED_NrYProcs|'${NProc[1]}'|;
            s|SED_NrZProcs|'${NProc[2]}'|;
            s|SED_Beta|'$beta'|;
            s|SED_Kappa|'$kappa'|;
            s|SED_2kappamu|'$two_kappamu'|;
            s|SED_GaugeConfigInputFile|Conf|;
            s|SED_IndexEnd|3|;
            s|SED_SolverPrecision|'$inversion_prec'|
            ' > $base_inv/inverter.input

	  cd $base_inv
	  
          #link configuration and source
	  ln -vfs $base_conf/Conf Conf.$conf
	  for i in $(seq -f%02.0f 00 $last_prop_index)
	  do
            ln -vsf $base_conf/Sources/$source_name/source.$i $base_inv/source.$conf.00.$i
	  done
          
	  $MPI_TM_PREF $base_tmLQCD/invert -f inverter.input
	  
          #take time
	  tac=$tic
	  tic=$(date +%s)
	  echo "Time to invert the source: "$(($tic-$tac)) >> $base_conf/time_log
        
	  echo
	  echo "Inversion finished"
	  echo
	  
      else
	  echo "Already inverted"
	  echo
      fi 

      cd $base_conf
      
    done
    
done