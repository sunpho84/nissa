#!/bin/bash

#set the base folder
source ~/nissa_conf.sh

echo
echo "Two points function calculation"
echo "-------------------------------"
echo "Working dir: $PWD"
echo

#reset the time log
rm -f $base_conf/time_log

#take initial time
tic=$(date +%s)

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
            sed '
                 s|SED_NL|'$L'|;
                 s|SED_NT|'$T'|;
                 s|SED_NoiseType|'$source_nois'|;
                 s|SED_Wall_T_Pos|'$tsource'|;
                 s|SED_Seed|'$source_seed'|' $base_protos/generate_sthoc_wall_source_input.xml > generate_stochastic_source_input.xml
            $MPI_AH_PREF $base_ahmidas/example/generate_stochastic_source
            
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
        prog_contr=applications/contract_meson_2pts
    elif [ $source_type == "Wall4" ]
    then
        vol_fact=$(( $L * $L * $L ))
        prog_contr=applications/contract_meson_stochastic_2pts
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
	
	ncombo=$(( 4 * $ntheta * $(( $ntheta + 1 )) / 2 * $nmu * $(( $nmu + 1 )) / 2 ))
	echo "Ncombo: "$ncombo
	
        nprop=$(( $ntheta * $nmu ))
        echo "Nprop: "$nprop

	echo $L $T >> $base_2pts/input
        echo $kappa >> $base_2pts/input
        echo $tsource >> $base_2pts/input
        echo $vol_fact >> $base_2pts/input
        echo $base_conf/Conf >> $base_2pts/input
        echo $nmicro >> $base_2pts/input
	cat $base_2pts/micro_correlations >> $base_2pts/input
        echo $nprop >> $base_2pts/input
	for((itheta=0;itheta<ntheta;itheta++))
	do
	    
	    theta=${list_theta[$itheta]}
	    
	    for((imu1=0;imu1<nmu;imu1++))
	    do
                mu=${list_mu[$imu1]}
                echo $base_conf/Props/$source_name/$theta/$mu/prop. $mu $theta  >> $base_2pts/input
	    done
	    
	done
	
	echo $ncombo >> $base_2pts/input
	
	list_targ_combo=""
	for itheta1 in $two_points_theta1
	do
	    theta1=${list_theta[$itheta1]}
	    for((itheta2=0;itheta2<ntheta;itheta2++))
	    do
		theta2=${list_theta[$itheta2]}
		for((if1=0;if1<2;if1++))
		do
		    for((if2=0;if2<2;if2++))
		    do
			for((imu1=0;imu1<nmu;imu1++))
			do
			    mu1=${list_mu[$imu1]}
			    for((imu2=0;imu2<nmu;imu2++))
			    do
				mu2=${list_mu[$imu2]}
				
				iprop1=$(( $itheta1 * $nmu + $imu1 ))
				iprop2=$(( $itheta2 * $nmu + $imu2 ))
				
				targ_combo=$base_2pts/$theta1/$mu1/$if1/$theta2/$mu2/$if2/
				mkdir -pv $targ_combo
				
				echo $iprop1 $if1 $iprop2 $if2 $targ_combo >> $base_2pts/input
				
				list_targ_combo=$list_targ_combo$targ_combo" "
			    done
			done
		    done
		done
	    done
	done
	
        $MPI_AH_PREF $base_ahmidas/$prog_contr $base_2pts/input
	
        #let's put together all the micro-correlations needed for the macro correlations
        #put all the needed volume factor, and translate the corr. to the origin
	for targ_combo in $list_targ_combo
	do

	    echo $targ_combo
	    cd $targ_combo

            for contr in ${two_points_correlations[@]}
            do
                awk '
                    {a=a$3" ";b=b$1"_"$2" "}
                 END{print a;system("paste "b)}' $base_nissa/Data/Correlations_content/$contr|awk '
               NR==1{n=NF;for(i=1;i<=n;i++)c[i]=$i/n}
                NR>1{x=0;y=0;for(i=1;i<=n;i++){x+=$(2*i-1)*c[i];y+=$(2*i)*c[i]};printf("%.12g\t%.12g\n",x,y)}' > $targ_combo/$contr
            done
            #rm -fv *_*

	    cd $OLDPWD

	done

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