#!/bin/bash

#set the base folder
source ~/nissa_conf.sh

echo
echo " ##################################################################################"
echo " ##################################################################################"
echo " ####                                                                          ####"
echo " ####                    Two-point functions calculation                       ####"
echo " ####                                                                          ####"
echo " ##################################################################################"
echo " ##################################################################################"
echo
echo "Started at : ", $(date)
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

echo "Two-point functions calculation started at : "$(date) >> $base_conf/time_log

#count the number of mu and theta
nmu=${#list_mu[@]}
ntheta=${#list_theta[@]}

echo "Number of sources: "$nsource
echo

for((is=0;is<nsource;is++))
do

    echo "######################## FIRST STEP: Source generation ##########################"
    echo
    echo "Initial Source number "$is" generation"
    echo

    #parse source parameters
    source_type=${list_source_type[$is]}
    source_nois=${list_source_nois[$is]}
    source_name=${list_source_name[$is]}
    source_pars=($(echo ${list_source_pars[$is]}|sed 's|_| |g'))
    source_seed=$((${list_source_seed[$is]}+$additive_seed))
    source_prec=${list_source_prec[$is]}

    tsource=$(( $((${source_pars[0]}+$additive_time_offset)) % $T ))

    if [ $source_type == "Point12" ]
    then
        last_prop_index=11
    elif [ $source_type == "Wall4" ]
    then
        last_prop_index=3
    elif [ $source_type == "Wall1" ]
    then
        last_prop_index=0
    else
	echo "Unknown source type"$source_type
	exit
    fi
    
    echo "Generating source named: "$source_name" of type: "$source_type" with pars: "${source_pars[@]}" and seed: "$source_seed

    targ_dir=Sources/$source_name

    if [ -d "$targ_dir" ]
    then
        echo "Source "$source_name" already existing"
        echo
    else
        
        mkdir -pv $targ_dir
        cd $targ_dir
	
	if [ $source_type == Point12 ]
        then
            
            #Generate the point12 source
            sed '
                 s|SED_NL|'$L'|;
                 s|SED_NT|'$T'|;
                 s|SED_Pos|'${source_pars[0]}'\ '${source_pars[1]}'\ '${source_pars[2]}'\ '${source_pars[3]}'|
                ' $base_protos/generate_point_source_input.xml > generate_point_source_input.xml
            
	    #invoke the program
	    $MPI_AH_PREF $base_ahmidas/example/generate_point_source
        elif [ $source_type == Wall4 ]
        then
            
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
            
	elif [ $source_type == Wall1 ]
        then
            
            #Generate the wall1 source
            sed '
                 s|SED_NL|'$L'|;
                 s|SED_NT|'$T'|;
                 s|SED_NoiseType|'$source_nois'|;
                 s|SED_Wall_T_Pos|'$tsource'|;
                 s|SED_Seed|'$source_seed'|' $base_protos/generate_ultrasthoc_wall_source_input.xml > generate_stochastic_source_input.xml
            $MPI_AH_PREF $base_ahmidas/example/generate_stochastic_source
            
        else
            echo "Unknown source type: "$source_type
	    exit
        fi

        if [ ! "$?" ]
        then
	    echo "Source generated with success!"
#            rm generate_*source*.xml
        fi

	cd -

    fi

    #take time
    tac=$tic
    tic=$(date +%s)
    echo "Time to generate source: "$(($tic-$tac)) >> $base_conf/time_log

    echo "######################## SECOND STEP: Inversion of doublet ############################"

    for theta in ${list_theta[@]}
    do

      echo
      echo "First inversion"
      echo 
      echo "In this workflow we invert only for all theta "
      echo
      
      base_inv=$base_conf/Props/$source_name/$theta/
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
            s|SED_ThetaX|'$theta'|;
            s|SED_ThetaY|'$theta'|;
            s|SED_ThetaZ|'$theta'|;
            s|SED_NrXProcs|'${NProc[0]}'|;
            s|SED_NrYProcs|'${NProc[1]}'|;
            s|SED_NrZProcs|'${NProc[2]}'|;
            s|SED_Beta|'$beta'|;
            s|SED_Kappa|'$kappa'|;
            s|SED_2kappamu|'$two_kappamu'|;
            s|SED_CgmmsNoExtraMasses|'$(( ${#list_mu[@]} - 1 ))'|;
            s|SED_GaugeConfigInputFile|Conf|;
            s|SED_IndexEnd|'$last_prop_index'|;
            s|SED_SolverPrecision|'$source_prec'|
            ' > $base_inv/inverter.input

	  #write the instructions of what to write
	  (
	      echo 1
	      echo 1
	      echo $IO_prec
	  ) > $base_inv/cg_mms_instructions

	  OLD=$PWD
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
        
          #clean input and log file
#         rm -vf extra_masses.input inverter.input output.para
        
          #clean link to conf and source
#         rm -vf source.$conf.00.?? Conf.$conf
        
          #clean useless first mass up output
#         rm -vf source.$conf.00.??.inverted
#         rm -vf .source*
	  
	  echo
	  echo "Inversion finished"
	  echo
	  
          #Move each mass in appropriate folder
	  for im in $(seq -f%02.0f 00 $(( ${#list_mu[@]} - 1 )) )
          do
	    
	    mu=${list_mu[10#$im]}
	    
            for ics in $(seq -f%02.0f 00 $last_prop_index)
            do

	      for r1 in 0 1
	      do
		
		mkdir -vp $base_inv/$mu/$r1

		orig=$base_inv/source.$conf.00.$ics.cgmms.$im.inverted.$r1
		dest=$base_inv/$mu/$r1/prop.$ics
		
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
	  
	  touch completed
	  
	  cd $OLD
	  
      else
	  echo "Already inverted"
	  echo
      fi 
      
      echo

    done # theta
      
    echo "######################## THIRD STEP: Contractions ###########################"
    echo

    if [ $source_type == "Point12" ]
    then
        vol_fact=1
        prog_contr=applications/contract_meson_2pts
    elif [ $source_type == "Wall4" ]
    then
	prog_contr="Appretto/projects/meson_2pts/meson_2pts input"
    elif [ $source_type == "Wall1" ]
    then
        vol_fact=$(( $L * $L * $L ))
        prog_contr=applications/contract_meson_ultrastochastic_2pts
    fi

    base_2pts=$base_conf/2pts/$source_name

    if [ ! -f $base_2pts/completed ]
    then

        mkdir -vp $base_2pts
        cd $base_2pts
        	
	(
	    cd $base_nissa/Data/Correlations_content/
	    cat ${two_points_correlations[@]}
	    cd $OLDPWD
	) | awk '{print $1,$2}' > $base_2pts/micro_correlations
	nmicro=$(wc $base_2pts/micro_correlations | awk '{print $1}')
	
        nprop1=$(( 2 * $nmu ))
        nprop2=$(( $nprop1 * $ntheta ))
        echo "Nprop1: "$nprop1
        echo "Nprop2: "$nprop2

	echo "L "$L  > $base_2pts/input
	echo "T "$T >> $base_2pts/input
	echo "TWall "$tsource >> $base_2pts/input
	echo "NPropFirstList "$nprop1 >> $base_2pts/input
	for itheta1 in ${two_points_theta1[@]}
	do
	  theta1=${list_theta[$itheta1]}
	  for((imu1=0;imu1<nmu;imu1++))
	  do
	    mu1=${list_mu[$imu1]}
	    for((r1=0;r1<2;r1++))
	    do
	      echo " "$base_conf/Props/$source_name/$theta1/$mu1/$r1/prop $mu1 $theta1 0 $r1 >> $base_2pts/input
	    done
	  done
	done
	echo "NPropSecondList "$nprop2 >> $base_2pts/input
	for((itheta2=0;itheta2<ntheta;itheta2++))
	do
	  theta2=${list_theta[$itheta2]}
	  for((imu2=0;imu2<nmu;imu2++))
	  do
	    mu2=${list_mu[$imu2]}
	    for((r2=0;r2<2;r2++))
	    do
	      echo " "$base_conf/Props/$source_name/$theta2/$mu2/$r2/prop $mu2 $theta2 0 $r2 >> $base_2pts/input
	    done
	  done
	done

	echo "Ncontr "$nmicro >> $base_2pts/input
	cat $base_2pts/micro_correlations >> $base_2pts/input
	echo "Output "$base_2pts"/two_points_contractions" >> $base_2pts/input

	echo
	echo "Launching program: "$prog_contr
	echo
        $MPI_TM_PREF $base_nissa/$prog_contr $base_2pts/input
	
        #take time
	tac=$tic
	tic=$(date +%s)
	echo "Time to make contractions: "$(($tic-$tac)) >> $base_conf/time_log
	
	touch completed
	
	cd $base_conf

    else
	
	echo "Contractions already performed"

    fi

done

echo "Two-point functions calculation ended at : "$(date) >> $base_conf/time_log

echo
echo "Two-point functions calculation ended at : " $(date)
echo

