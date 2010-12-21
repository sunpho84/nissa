#!/bin/bash

#set the base folder
source ~/nissa_conf.sh

echo
echo " ###############################################################"
echo " ###############################################################"
echo " ####                                                       ####"
echo " ####                   Bk calculation                      ####"
echo " ####                                                       ####"
echo " ###############################################################"
echo " ###############################################################"
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

seedLR=( 0 89723)

echo "Bk calculation started at : "$(date) >> $base_conf/time_log

#count the number of mu
nmu=${#list_mu[@]}

echo "Number of sources: "$nsource
echo

TH=$(( $T / 2 ))

for((is=0;is<nsource;is++))
do
  
  #parse source parameters                                                                                                 
  source_type=${list_source_type[$is]}
  source_nois=${list_source_nois[$is]}
  source_pars=($(echo ${list_source_pars[$is]}|sed 's|_| |g'))
  source_prec=${list_source_prec[$is]}
  
  is_tsource=$(( $((${source_pars[0]}+$additive_time_offset)) % $TH ))
  is_source_seed=$(( $base_seed + $additive_seed ))
  is_source_name=${list_source_name[$is]}

  for LR in 0 1
  do
    
    echo "######################## FIRST STEP: Source generation ##########################"
    echo
    echo "Generation of "$is" source"
    echo
    
    
    source_name=$is_source_name"_"$LR

    source_dir=$base_conf/Sources/$source_name

    if [ -d "$source_dir" ]
    then
        echo "Source "$source_name" already existing"
        echo
    else
        
        mkdir -pv $source_dir
        cd $source_dir
	
	if [ "$LR" -eq 1 ]
	then
	    tsource=$(( $is_tsource + $TH ))
	    source_seed=$(( $is_source_seed + ${seedLR[$LR]} ))
	else
	    tsource=$is_tsource
	    source_seed=$is_source_seed
	fi
	
        #Generate the wall4 source
        (
	    echo "L "$L
	    echo "T "$T
	    echo "Seed "$source_seed
	    echo "TakeSlice 1"
	    echo "TWall "$tsource
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

    echo
    
    base_inv=$base_conf/Props/$source_name/
    mkdir -pv $base_inv
    
    echo "Inverting: "$source_name" "$theta
    
    if [ ! -f $base_inv/completed ]
    then
	
        #create additional masses list: first mass is included as base
	for mu in ${list_mu[@]};do echo $mu;done|awk 'NR>1{printf("%.12f\n",2*$1*'$kappa')}' > $base_inv/extra_masses.input

        #prepare input for inverter
	mu=${list_mu[0]}
	two_kappamu=$(echo $kappa|awk '{printf("%.12f\n",2*$1*'$mu')}')
	cat $base_protos/cgmms.input|sed '
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
            s|SED_CgmmsNoExtraMasses|'$(( ${#list_mu[@]} - 1 ))'|;
            s|SED_GaugeConfigInputFile|Conf|;
            s|SED_IndexEnd|3|;
            s|SED_SolverPrecision|'$source_prec'|
            ' > $base_inv/inverter.input
	
        #write the instructions of what to write
	(
	    echo 1
	    echo 1
	    echo $IO_prec
	) > $base_inv/cg_mms_instructions
	
	cd $base_inv
	
        #link configuration and source
	ln -vfs $base_conf/Conf Conf.$conf
	for i in $(seq -f%02.0f 00 03)
	do
	  ln -vsf $base_conf/Sources/$source_name/source.$i $base_inv/source.$conf.00.$i
	done
          
	$MPI_TM_PREF $base_tmLQCD/invert -f inverter.input
	
        #Move each mass in appropriate folder
	for im in $(seq -f%02.0f 00 $(( $nmu - 1 )) )
        do
	  
	  mu=${list_mu[10#$im]}
	  
	  for ics in $(seq -f%02.0f 00 3)
          do
	    
	    for r in 0 1
            do
	      
	      mkdir -vp $base_inv/$mu/$r
	      
	      orig=$base_inv/source.$conf.00.$ics.cgmms.$im.inverted.$r
	      dest=$base_inv/$mu/$r/prop.$ics.fuf
	      
	      if [ ! -f $orig ]
		  then
		  echo "Could not find: "$orig
		  exit
	      else
		  mv -v $orig $dest
	      fi
	      
	    done
	    
	  done
	  
	done
	
        #take time
	tac=$tic
	tic=$(date +%s)
	echo "Time to invert the source: "$LR $(($tic-$tac)) >> $base_conf/time_log
        
	echo
	echo "Inversion finished"
	echo
	
    else
	echo "Already inverted"
	echo
    fi 
    
    cd $base_conf
    
  done
  
  name_list=(First Second Third Fourth)
  nmu_list=($nmu_low $nmu $nmu_low $nmu)
  LR_list=(0 0 1 1)

  echo "######################## THIRD STEP: Three point function calculation ############################"
  echo
  
  base_Bk=$base_conf/Bk/$is_source_name
  if [ ! -f $base_Bk/completed ]
  then
      mkdir -vp $base_Bk
  fi
  
  cd $base_Bk
  
  (
      cd $base_nissa/Data/Correlations_content/
      cat ${two_points_correlations[@]}
      cd $OLDPWD
  ) | awk '{print $1,$2}' > $base_Bk/micro_correlations
  nmicro=$(wc $base_Bk/micro_correlations | awk '{print $1}')

  echo "L "$L  > $base_Bk/input
  echo "T "$T >> $base_Bk/input
  echo "TWall "$is_tsource >> $base_Bk/input
  for((ilist=0;ilist<4;ilist++))
  do
    echo "NProp"${name_list[$ilist]}"List "$((${nmu_list[$ilist]} * 2)) >> $base_Bk/input
    for((imu=0;imu<${nmu_list[$ilist]};imu++))
    do
      mu=${list_mu[$imu]}
      for((r=0;r<2;r++))
      do
	echo " "$base_conf/Props/$is_source_name\_${LR_list[$ilist]}/$mu/$r/prop fuf $mu 0.00 0 $r >> $base_Bk/input
      done
    done
  done
  echo "NContr_eight "${#list_ops_eight[@]} >> $base_Bk/input
  echo ${list_ops_eight[@]} >> $base_Bk/input
  echo "NContr_2points "$nmicro >> $base_Bk/input
  cat $base_Bk/micro_correlations >> $base_Bk/input
  echo "Output_eight eights" >> $base_Bk/input
  echo "Output_2points two_points" >> $base_Bk/input
  
  $MPI_TM_PREF $base_nissa/Appretto/projects/eight_BK/eight_BK input
  
done

echo "Bk calculation ended at : "$(date) >> $base_conf/time_log
