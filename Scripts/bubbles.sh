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
echo "Started at : " $(date)
echo
echo "Working dir: $PWD"
echo

#take initial time
tic=$(date +%s)

#reset the time log
if [ -f "$base_conf/time_log" ]
then
    mv -f $base_conf/time_log $base_conf/time_log_$tic
fi

echo "Bubbles calculation started at : "$(date) >> $base_conf/time_log

#count the number of mu
nmu=${#list_mu[@]}

echo "Number of sources: "$nsources
echo

for((is=0;is<nsources;is++))
do

    echo "######################## FIRST STEP: Source generation ##########################"
    echo
    echo "Generation of "$is" source"
    echo
    
    source_seed=$base_seed+$is
    source_name=$(printf %02d $is)

    source_dir=$base_conf/Sources/$source_name

    if [ -d "$source_dir" ]
    then
        echo "Source "$source_name" already existing"
        echo
    else
        
        mkdir -pv $source_dir
        cd $source_dir
	
        #Generate the wall4 source
        (
	    echo "L "$L
	    echo "T "$T
	    echo "Seed "$source_seed
	    echo "TakeSlice 0"
	    echo "TWall -1"
	    echo "NoiseType "$source_noise
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
            s|SED_SolverPrecision|'$inversion_precision'|
            ' > $base_inv/inverter.input

	  cd $base_inv
	  
          #link configuration and source
	  ln -vfs $base_conf/Conf Conf.$conf
	  for i in $(seq -f%02.0f 00 03)
	  do
            ln -vsf $base_conf/Sources/$source_name/source.$i $base_inv/source.$conf.00.$i
	  done
          
	  $MPI_TM_PREF $base_tmLQCD/invert -f inverter.input
	  
	  for i in $(seq -f%02.0f 00 03)
	  do
              mv $base_inv/source.$conf.00.$i.inverted $base_inv/prop.$i
	  done

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

    echo "######################## THIRD STEP: Three point function calculation ############################"
    echo
    
    base_bubbles=$base_conf/bubbles/$source_name\_$is
    if [ ! -f $base_bubbles/completed ]
    then
        mkdir -vp $base_bubbles
    fi
    
    cd $base_bubbles
    
    echo "L "$L  > $base_bubbles/input
    echo "T "$T >> $base_bubbles/input
    echo "TWall 0" >> $base_bubbles/input
    echo "Source "$source_dir"/source" >> $base_bubbles/input
    echo "NContr 16" >> $base_bubbles/input
    seq 0 15 >> $base_bubbles/input
    echo "Output bubbles" >> $base_bubbles/input
    echo "NProp "$nmu >> $base_bubbles/input
    for mu in ${list_mu[@]}
    do
	base_inv=$base_conf/Props/$source_name\_$mu/
      	r=$( echo $mu | awk '{print ($1>0)}' )
        echo $base_inv/prop "m="$mu 0 $r >> $base_bubbles/input
    done
    
    $MPI_TM_PREF $base_nissa/Appretto/projects/bubbles/bubbles input
    
done

echo "Bubbles calculation ended at : "$(date) >> $base_conf/time_log

echo
echo "Bubbles calculation ended at : " $(date)
echo