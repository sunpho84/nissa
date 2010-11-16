#!/bin/bash

#set the base folder
source ~/nissa_conf.sh

echo
echo "Two points function calculation"
echo "-------------------------------"
echo "Working dir: $PWD"
echo

#take initial time
tic=$(date +%s)
rm -f $base_conf/time_log

nmu=${#list_mu[@]}

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
    source_seed=${list_source_seed[$is]}
    source_prec=${list_source_prec[$is]}
    tsource=${source_pars[0]}

    if [ $source_type == "Point12" ]
    then
        last_prop_index=11
        reco=reconstruct_doublet
    elif [ $source_type == "Wall4" ]
    then
        last_prop_index=3
        reco=reconstruct_stochastic_doublet
    elif [ $source_type == "Wall1" ]
    then
        last_prop_index=0
        reco=reconstruct_ultrastochastic_doublet
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
	    echo $MPI_AH_PREF $base_ahmidas/example/generate_point_source
	    $MPI_AH_PREF $base_ahmidas/example/generate_point_source
        elif [ $source_type == Wall4 ]
        then
            
            #Generate the wall4 source
            sed '
                 s|SED_NL|'$L'|;
                 s|SED_NT|'$T'|;
                 s|SED_NoiseType|'$source_nois'|;
                 s|SED_Wall_T_Pos|'${source_pars[0]}'|;
                 s|SED_Seed|'$source_seed'|' $base_protos/generate_sthoc_wall_source_input.xml > generate_stochastic_source_input.xml
	    echo $MPI_AH_PREF $base_ahmidas/example/generate_stochastic_source
            $MPI_AH_PREF $base_ahmidas/example/generate_stochastic_source
            
	elif [ $source_type == Wall1 ]
        then
            
            #Generate the wall1 source
            sed '
                 s|SED_NL|'$L'|;
                 s|SED_NT|'$T'|;
                 s|SED_NoiseType|'$source_nois'|;
                 s|SED_Wall_T_Pos|'${source_pars[0]}'|;
                 s|SED_Seed|'$source_seed'|' $base_protos/generate_ultrasthoc_wall_source_input.xml > generate_stochastic_source_input.xml
	    echo $MPI_AH_PREF $base_ahmidas/example/generate_stochastic_source
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
      
      targ=$base_conf/Props/$source_name/$theta/
      
      echo "Inverting: "$source_name" "$theta
      
      if [ ! -f $targ/completed ]
      then
        #prepare the list of folders where to put final data
	  for mu in ${list_mu[@]}
	  do
            mkdir -pv $targ/$mu/
	  done
        
          #create additional masses list: first mass is included as base
	  for mu in ${list_mu[@]};do echo $mu;done|awk 'NR>1{printf("%.12f\n",2*$1*'$kappa')}' > $targ/extra_masses.input
	  
          #prepare input for inverter
	  mu=${list_mu[0]}
	  two_kappamu=$(echo $kappa|awk '{printf("%.12f\n",2*$1*'$mu')}')
	  cat $base_protos/cgmms.input|sed '
            s|SED_Confno|'$confno'|
            s|SED_NL|'$L'|;
            s|SED_NT|'$T'|;
            s|SED_ThetaX|'$theta'|;
            s|SED_ThetaY|'$theta'|;
            s|SED_ThetaZ|'$theta'|;
            s|SED_Beta|'$beta'|;
            s|SED_Kappa|'$kappa'|;
            s|SED_2kappamu|'$two_kappamu'|;
            s|SED_CgmmsNoExtraMasses|'$(( ${#list_mu[@]} - 1 ))'|;
            s|SED_GaugeConfigInputFile|Conf|;
            s|SED_IndexEnd|'$last_prop_index'|;
            s|SED_SolverPrecision|'$source_prec'|
            ' > $targ/inverter.input
	  
	  OLD=$PWD
	  cd $targ
	  
          #link configuration and source
	  ln -vfs $base_conf/Conf Conf.$confno
	  for i in $(seq -f%02.0f 00 $last_prop_index)
	  do
            ln -vsf $base_conf/Sources/$source_name/source.$i source.$confno.00.$i
	  done
        
	  $MPI_TM_PREF $base_tmLQCD/invert -f inverter.input

          #take time
	  tac=$tic
	  tic=$(date +%s)
	  echo "Time to invert the source: "$(($tic-$tac)) >> $base_conf/time_log
        
          #clean input and log file
#         rm -vf extra_masses.input inverter.input output.para
        
          #clean link to conf and source
#         rm -vf source.$confno.00.?? Conf.$confno
        
          #clean useless first mass up output
#         rm -vf source.$confno.00.??.inverted
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
	      
	      orig=source.$confno.00.$ics.cgmms.$im.inverted
	      dest=$mu/prop.$ics
	      
	      if [ ! -f $orig ]
	      then
		  echo "Could not find: "$orig
		  exit
	      else
		  mv -v $orig $dest
	      fi		    
	      
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
        prog_contr=applications/contract_multi_lines
    elif [ $source_type == "Wall4" ]
    then
        vol_fact=$(( $L * $L * $L ))
        prog_contr=applications/contract_multi_stochastic_lines
    elif [ $source_type == "Wall1" ]
    then
        vol_fact=$(( $L * $L * $L ))
        prog_contr=applications/contract_multi_ultrastochastic_lines
    fi
    
    targ=$base_conf/2pts/$source_name/$theta/
    
    if [ ! -f $targ/correlators.dat ]
    then
        mkdir -vp $targ
        cd $targ
        
        (
            echo $L $T
            echo $kappa
            echo $theta $theta $theta 1
            echo $tsource
            echo $vol_fact
            echo $base_conf/Conf
            echo 4
            echo 5 5
            echo -1 5
            echo 5 -1
            echo -1 -1
            echo $nmu
            for((imu1=0;imu1<nmu;imu1++))
            do
                mu=${list_mu[$imu1]}
                echo $mu $base_conf/Props/$source_name/$theta/$mu/prop.
            done
        ) > input
	
        echo $MPI_AH_PREF $base_ahmidas/$prog_contr input
        $MPI_AH_PREF $base_ahmidas/$prog_contr input
  
        #take time
	tac=$tic
	tic=$(date +%s)
	echo "Time to make contractions: "$(($tic-$tac)) >> $base_conf/time_log

	cd -

    fi

done