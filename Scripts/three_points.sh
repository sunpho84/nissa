#!/bin/bash

#set the base folder
source ~/nissa_conf.sh

echo
echo "Three points function calculation"
echo "-------------------------------"
echo "Working dir: $PWD"
echo

nmu=${#list_mu[@]}

for((is=0;is<nsource;is++))
do

    echo "######################## FIRST STEP: Sequential source generation ##########################"
    echo
    echo "Initial Source number "$is" generation"
    echo

    #parse source parameters
    source_type=${list_source_type[$is]}
    source_name=${list_source_name[$is]}
    source_pars=($(echo ${list_source_pars[$is]}|sed 's|_| |g'))
    source_seed=$((${list_source_seed[$is]}+$additive_seed))
    source_prec=${list_source_prec[$is]}

    tsource=$(( $((${source_pars[0]}+$additive_time_offset)) % $T ))

    if [ $source_type == "Point12" ]
    then
        last_prop_index=11
        finalize=finalize_meson_sequential_propagator
        sequ=generate_meson_sequential_source
    elif [ $source_type == "Wall4" ]
    then
        last_prop_index=3
        finalize=finalize_meson_sequential_stochastic_propagator
        sequ=generate_meson_sequential_stochastic_source
    elif [ $source_type == "Wall1" ]
    then
        last_prop_index=0
        finalize=finalize_meson_sequential_ultrastochastic_propagator
        sequ=generate_meson_sequential_ultrastochastic_source
    fi

    #loops over the spectator
    for itheta_spec in $list_itheta_spec
    do

      theta_spec=${list_theta[$itheta_spec]}

      for imu_spec in $list_imu_spec
      do

	mu_spec=${list_mu[$imu_spec]}
      
	for f_spec in $list_f_spec
	do
	
	  echo "Generating sequential source from the original named: "$source_name" of type: "$source_type" with pars: "${source_pars[@]}" and seed: "$source_seed", for the "$theta_spec" theta, "$mu_spec" mass, "$f_spec" flavour"
	
	  source_dir=$base_conf/Props/$source_name/$theta_spec/$mu_spec/$f_spec
	  targ_dir=$base_conf/SeqSources/$source_name/$theta_spec/$mu_spec/$f_spec
	  
	  if [  -d "$targ_dir" ]
	  then
	      echo "Sequential Source "$targ_dir" already existing"
	      echo
	  else
	      
	      mkdir -pv $targ_dir
	      cd $targ_dir
	      
	      #link the S0 propagator
	      for ics in $(seq -f%02.0f 00 $last_prop_index)
	      do
		ln -sfv $source_dir/prop.$ics
	      done
	    
              #Generate the sequential p5 source input file
	      TH=$(( $T / 2 ))

	      if [ $tsource -lt $TH ]
	      then
		  TSlice=$(( $tsource + $TH ))
	      else
		  TSlice=$(( $tsource - $TH ))
	      fi
	      sed '
                s|SED_NL|'$L'|;
                s|SED_NT|'$T'|;
                s|SED_Take_Slice|1|;
                s|SED_Chosen_Slice|'$TSlice'|;
                s|SED_S0_Flav|'$f_spec'|
                s|SED_IndexEnd|'$last_prop_index'|;
                s|SED_ThetaX|'$theta_spec'|;
                s|SED_ThetaY|'$theta_spec'|;
                s|SED_ThetaZ|'$theta_spec'|;
                ' $base_protos/generate_sequential_source_input.xml > generate_sequential_source_input.xml
	    
   	      #invoke the program
	      $MPI_AH_PREF $base_ahmidas/applications/$sequ generate_sequential_source_input.xml
	      
	      if [ "$?" ]
	      then
		  
		  echo "Source generated with success!"
  		  #rm generate_*source*.xml
		  
		  for ics in $(seq -f%02.0f 00 $last_prop_index)
		  do
		    #rm -fv prop.$ics
		    mv -v prop.$ics.seq source.$ics
		  done
		  
	      fi
	      
	      cd -
	      
	  fi

	  #ok now the source have been created
	  
	  #choose the f1 flavor to have the charged combo
	  if [ $f_spec == 0 ]
	  then
	      f1=1
	  else
	      f1=0
	  fi
	  
	  echo "######################## SECOND STEP: Inversion of doublet ############################"
	  echo
	  echo "Second inversion"
	  echo 
	
	  #this is the list of theta to use on the sequential line
	  for theta1 in ${list_theta[@]}
	  do
	    
	    targ=$base_conf/SeqProps/$source_name/$theta_spec/$mu_spec/$f_spec/$theta1/
	    
	    echo "Inverting: "$source_name" "$theta1
	    
	    if [ ! -f $targ/completed ]
	    then
                #prepare the list of folders where to put final data
		for mu1 in ${list_mu[@]}
		do
		    mkdir -pv $targ/$mu1/
		done
		
                #create additional masses list: first mass is included as base
		for mu1 in ${list_mu[@]};do echo $mu1;done|awk 'NR>1{printf("%.12f\n",2*$1*'$kappa')}' > $targ/extra_masses.input

                #prepare input for inverter
		mu1=${list_mu[0]} #controlla
		two_kappamu=$(echo $kappa|awk '{printf("%.12f\n",2*$1*'$mu1')}')
		cat $base_protos/cgmms.input|sed '
                    s|SED_Confno|'$confno'|;
                    s|SED_NL|'$L'|;
                    s|SED_NT|'$T'|;
                    s|SED_ThetaX|'$theta1'|;
                    s|SED_ThetaY|'$theta1'|;
                    s|SED_ThetaZ|'$theta1'|;
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
		    ln -vsf $base_conf/SeqSources/$source_name/$theta_spec/$mu_spec/$f_spec/source.$i source.$confno.00.$i
		done
		
		$MPI_TM_PREF $base_tmLQCD/invert -f inverter.input
		
                #clean input and log file
		#rm -vf extra_masses.input inverter.input output.para
		
                #clean link to conf and source
		#rm -vf source.$confno.00.?? Conf.$confno
		
                #clean useless first mass up output
		#rm -vf source.$confno.00.??.inverted
		#rm -vf .source*
		
		echo
		echo "Inversion finished"
		echo
        
                #Move each mass in appropriate folder and finalize the sequential propagator
		for imu1 in $(seq -f%02.0f 00 $(( ${#list_mu[@]} - 1 )) )
		do
		    
		    mu1=${list_mu[$imu1]}
		    
		    for ics in $(seq -f%02.0f 00 $last_prop_index)
		    do
			
			orig=source.$confno.00.$ics.cgmms.$imu1.inverted
			dest=$mu1/prop.$ics
			
			if [ ! -f $orig ]
			then
			    echo "Could not find: "$orig
#			    exit
			else
			    mv -v $orig $dest
			fi		    
			
		    done
		
		    targ=$base_conf/SeqProps/$source_name/$theta_spec/$mu_spec/$f_spec/$theta1/$mu1/
		    
		    echo "Finalizing $source_name $theta1 $mu1 doublet"
		    
		    cd $targ
		    
                    #link the configuration
		    ln -vfs $base_conf/Conf Conf0
		    
                    cat $base_protos/finalize_meson_sequential_propagator_input.xml|sed '
                         s|SED_NL|'$L'|;
                         s|SED_NT|'$T'|;
                         s|SED_ThetaX|'$theta1'|;
                         s|SED_ThetaY|'$theta1'|;
                         s|SED_ThetaZ|'$theta1'|;
                         s|SED_Beta|'$beta'|;
                         s|SED_Kappa|'$kappa'|;
                         s|SED_Mu|'$mu1'|;
                         s|SED_S0_Flav|'$f_spec'|;
                         s|SED_Last_prop_index|'$last_prop_index'|;
                         ' > $targ/finalize_meson_sequential_propagator_input.xml
		    
		    $MPI_AH_PREF $base_ahmidas/applications/$finalize finalize_meson_sequential_propagator_input.xml
		    
                    #move the two flavour to appropriate folder
		    mkdir -pv $targ/$f1
		    for ics in $(seq -f%02.0f 00 $last_prop_index)
		    do
			mv -v $targ/prop.$ics.final $targ/$f1/prop.$ics
		    done
		    #rm -fv $targ/prop.*
		    
                    #remove the input file and the link to conf
		    #rm -vf Conf0 reconstruct_doublet_input.xml
		    
		    cd -
		    
		    echo "Doublet reconstructed"
		    echo

		done #loop over output mass

		touch completed
		
		cd $OLD
		
	    else
		echo "Already inverted"
		echo
	    fi
	    
	  done #loop over "sequential" theta
	  

	  echo "######################## THIRD STEP: Three point function calculation ############################"
	  echo
	
	  for((imu1=0;imu1<nmu;imu1++))
	  do
	    
	    mu1=${list_mu[$imu1]}
	    
	    for((imu2=imu1;imu2<nmu;imu2++))
	    do
	    
	      f2=$f_spec #this is the charged
	      theta2=$theta1 #to put -

	      mu2=${list_mu[$imu2]}
	    
	      echo "Contracting combination: "$mu_spec $mu1 $mu2 $f_spec $f1 $f2

	      source1=$base_conf/SeqProps/$source_name/$theta_spec/$mu_spec/$f_spec/$theta1/$mu1/$f1/
	      source2=$base_conf/Props/$source_name/$theta2/$mu2/$f2
	    
	      if [ ! -d $source1 ] || [ ! -d $source2 ]
	      then
		  echo "Could not find "$source1" or "$source2
		  exit
	      fi
	      
	      targ=$base_conf/3pts/$source_name/$theta_spec/$mu_spec/$f_spec/$theta1/$mu1/$f1/$theta2/$mu2/$f2/
	      
	      mkdir -pv $targ
	      
	      (
		  (
		      cat $base_protos/contract_two_lines_head.xml
		      cat $base_analysis/contract_two_lines_middle.xml
		      cat $base_protos/contract_two_lines_tail.xml
		      )|sed 's|SED_NL|'$L'|;
                             s|SED_NT|'$T'|;
                             s|SED_Line_a|'$source1'|;
                             s|SED_IndexEnd|'$last_prop_index'|;
                             s|SED_Line_b|'$source2'|;'
	      ) > $targ/contract_two_lines_input.xml
	      
	      cd $targ
	      
	      if [ -f P5P5 ]
	      then
		  echo "Already calculated"
	      else
		  
		  if [ $source_type == "Point12" ]
		  then
		      vol_fact=1
		      $MPI_AH_PREF $base_ahmidas/applications/contract_two_lines contract_two_lines_input.xml
		  elif [ $source_type == "Wall4" ]
		  then
		      vol_fact=$(( $L * $L * $L ))
		      $MPI_AH_PREF $base_ahmidas/applications/contract_two_stochastic_lines contract_two_lines_input.xml
		  elif [ $source_type == "Wall1" ]
		  then
		      vol_fact=$(( $L * $L * $L ))
		      $MPI_AH_PREF $base_ahmidas/applications/contract_two_ultrastochastic_lines contract_two_lines_input.xml
		  fi
		  
                  #rm -vf contract_two_lines_input.xml
		  
		  split -l $(( $T + 1 )) -d correlators.dat micro
		  
		  lim=($(awk '{print $1"_"$2}' $base_analysis/micro_correlations))
		  n=$(( ${#lim[@]} -1 ))
		  lin=($(seq -w 0 $n))
		  
		  for i in $(seq 0 $n)
		  do
		    echo ${lim[$i]}
		    mv micro${lin[$i]} name_micro_${lim[$i]}
		  done
		
 		  #let's put together all the micro-correlations needed for the macro correlations
                  #put all the needed volume factor, and translate the corr. to the origin                             
		  for op in $(cat $base_analysis/correlations_needed)
		  do
		    list_coef=""
		    list_file=""
		    awk '{a=a$3" ";b=b"name_micro_"$1"_"$2" "}
                      END{print a;system("paste "b)}' $base_nissa/Data/Correlations_content/$op|awk '
                    BEGIN{norm=1.0/'$vol_fact';T='$T';t=(T-'$tsource')%T}
                    NR==1{n=NF;for(i=1;i<=n;i++)c[i]=$i/n}
                     NR>2{for(i=1;i<=n;i++)
                             {x[t]+=$(3*i-1)*c[i];y[t]+=$(3*i)*c[i]}
                              t=(t+1)%T}
                      END{for(t=0;t<T;t++)printf("%.10e\t%.10e\n",x[t]*norm,y[t]*norm)}' > $op
		  done
		  rm -fv name_micro*
	      fi
	      
	      cd -
	      
	    done #loop over mu2
	  
	  done #loop over mu1

          #########################Checking P5P5#######################
	  
	  echo "Now checking two points functions"
	  echo
	  
	  for((imu1=imu_spec;imu1<nmu;imu1++))
	  do
	      
	      mu1=${list_mu[$imu1]}
	      theta1=$theta_spec
	      
	      targ=$base_conf/3pts/2pts_check/$source_name/$theta_spec/$mu_spec/$f_spec/$mu1/$f1
	      mkdir -pv $targ
	      
	      if [ -f $targ/P5P5 ]
	      then
		  echo "Two Points funcion already checked"
	      else
		  
		  cd $targ
		  
		  source1=$base_conf/SeqProps/$source_name/$theta_spec/$mu_spec/$f_spec/$theta1/$mu1/$f1
		  source2=$base_conf/Sources/$source_name/
		  
		  echo "Checking: "$source1
		  
		  sed '                                                                                              
                      s|SED_NL|'$L'|;
                      s|SED_NT|'$T'|;
                      s|SED_S0_Flav|'$f_spec'|
                      s|SED_Line_a|'$source1'|;
                      s|SED_IndexEnd|'$last_prop_index'|;
                      s|SED_Line_b|'$source2'|;' $base_protos/check_meson_3pts_input.xml > $targ/check_meson_3pts_input.xml
		  
		  if [ $source_type == "Point12" ]
		  then
		      vol_fact=1
		      $MPI_AH_PREF $base_ahmidas/example/check_meson_3pts
		  elif [ $source_type == "Wall4" ]
		  then
		      vol_fact=$(( $L * $L * $L ))
		      $MPI_AH_PREF $base_ahmidas/example/check_meson_stochastic_3pts
		  elif [ $source_type == "Wall1" ]
		  then
		      vol_fact=$(( $L * $L * $L ))
		      $MPI_AH_PREF $base_ahmidas/example/check_meson_ultrastochastic_3pts
		  fi
		  
		  awk 'NR==2+'$tsource'{print $1,$2/'$vol_fact',$3/'$vol_fact'}' $targ/correlators.dat > $targ/P5P5
		  
                    #rm -vf contract_two_lines_input.xml
		  
	      fi
	      
	  done
	  
	done #spectator flavor
	
      done #spectator mu
      
    done #spectator theta

done #source name