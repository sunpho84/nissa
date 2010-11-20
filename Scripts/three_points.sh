#!/bin/bash

#set the base folder
source ~/nissa_conf.sh

echo
echo "Three points function calculation"
echo "-------------------------------"
echo "Working dir: $PWD"
echo

#take initial time
tic=$(date +%s)

#count the number of mu and theta
nmu=${#list_mu[@]}
ntheta=${#list_theta[@]}

echo "Number of sources: "$nsource
echo

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
	vol_fact=1
        last_prop_index=11
        generate=generate_meson_sequential_propagator
        finalize=finalize_meson_sequential_propagator
        prog_contr=contract_meson_3pts
    elif [ $source_type == "Wall4" ]
    then
        vol_fact=$(( $L * $L * $L ))
        last_prop_index=3
        generate=generate_meson_sequential_stochastic_propagator
        finalize=finalize_meson_sequential_stochastic_propagator
        prog_contr=contract_meson_stochastic_3pts
    elif [ $source_type == "Wall1" ]
    then
        vol_fact=$(( $L * $L * $L ))
        last_prop_index=0
        generate=generate_meson_sequential_ultrastochastic_propagator
        finalize=finalize_meson_sequential_ultrastochastic_propagator
        prog_contr=contract_meson_ultrasthocastic_3pts
    fi

    #loops over the spectator
    for itheta_spec in ${list_itheta_spec[@]}
    do

      theta_spec=${list_theta[$itheta_spec]}

      for imu_spec in ${list_imu_spec[@]}
      do

	mu_spec=${list_mu[$imu_spec]}
      
	for f_spec in ${list_f_spec[@]}
	do
	
	  echo "Generating sequential source from the original named: "$source_name" of type: "$source_type" with pars: "${source_pars[@]}" and seed: "$source_seed", for the "$theta_spec" theta, "$mu_spec" mass, "$f_spec" flavour"
	
	  source_dir=$base_conf/Props/$source_name/$theta_spec/$mu_spec
	  targ_dir=$base_conf/SeqSources/$source_name/$theta_spec/$mu_spec/$f_spec
	  
	  if [ -d "$targ_dir" ]
	  then
	      echo "Sequential Source "$targ_dir" already existing"
	      echo
	  else
	      
	      mkdir -pv $targ_dir
	      cd $targ_dir
	      
	      #link the DD propagator
	      for ics in $(seq -f%02.0f 00 $last_prop_index)
	      do
		ln -sfv $source_dir/prop.$ics
	      done
	    
	      #link the configuration
	      ln -sfv $base_conf/Conf Conf0
	    
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
                s|SED_S0_Flav|'$f_spec'|;
                s|SED_Beta|'$beta'|;
                s|SED_Kappa|'$kappa'|;
                s|SED_Mu|'$mu_spec'|;
                s|SED_IndexEnd|'$last_prop_index'|;
                s|SED_ThetaX|'$theta_spec'|;
                s|SED_ThetaY|'$theta_spec'|;
                s|SED_ThetaZ|'$theta_spec'|;
                ' $base_protos/generate_sequential_source_input.xml > generate_sequential_source_input.xml
	    
   	      #invoke the program
	      $MPI_AH_PREF $base_ahmidas/applications/$generate generate_sequential_source_input.xml
	      
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
                    s|SED_Conf|'$conf'|;
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
		ln -vfs $base_conf/Conf Conf.$conf
		for i in $(seq -f%02.0f 00 $last_prop_index)
		  do
		    ln -vsf $base_conf/SeqSources/$source_name/$theta_spec/$mu_spec/$f_spec/source.$i source.$conf.00.$i
		done
		
		$MPI_TM_PREF $base_tmLQCD/invert -f inverter.input
		
                #clean input and log file
		#rm -vf extra_masses.input inverter.input output.para
		
                #clean link to conf and source
		#rm -vf source.$conf.00.?? Conf.$conf
		
                #clean useless first mass up output
		#rm -vf source.$conf.00.??.inverted
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
			
			orig=source.$conf.00.$ics.cgmms.$imu1.inverted
			dest=$mu1/prop.$ics
			
			if [ ! -f $orig ]
			then
			    echo "Could not find: "$orig
#			    exit
			else
			    mv -v $orig $dest
			fi		    
			
		    done
		
		    
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
	
	  base_3pts=$base_conf/3pts/$source_name/$itheta_spec/$imu_spec/$mu_spec
	  
	  if [ ! -f $base_3pts/completed ]
	  then
	      
              mkdir -vp $base_3pts
              cd $base_3pts
              
	      (
		  cd $base_nissa/Data/Correlations_content/
		  cat ${three_points_correlations[@]}
		  cd $OLDPWD
	      ) | awk '{print $1,$2}' > $base_3pts/micro_correlations 
	      nmicro=$(wc $base_3pts/micro_correlations | awk '{print $1}')
	      
              nprop=$(( $ntheta * $nmu ))
              echo "Nprop: "$nprop

	      ncombo=$(( $nprop * $nprop ))
	      echo "Ncombo: "$ncombo

	      rm -vf input
	      echo $L $T >> $base_3pts/input
              echo $kappa >> $base_3pts/input
	      echo $f_spec >> $base_3pts/input
              echo $tsource >> $base_3pts/input
              echo $vol_fact >> $base_3pts/input
              echo $base_conf/Conf >> $base_3pts/input
              echo $nmicro >> $base_3pts/input
	      cat $base_3pts/micro_correlations >> $base_3pts/input
              echo $nprop >> $base_3pts/input
	      for((itheta0=0;itheta0<ntheta;itheta0++))
	      do
		  theta0=${list_theta[$itheta0]}
		  for((imu0=0;imu0<nmu;imu0++))
		  do
                      mu0=${list_mu[$imu0]}
                      echo $base_conf/Props/$source_name/$theta0/$mu0/prop. $mu0 $theta0 >> $base_3pts/input
		  done
	      done
              echo $nprop >> $base_3pts/input
	      for((itheta1=0;itheta1<ntheta;itheta1++))
	      do
		  theta1=${list_theta[$itheta1]}
		  for((imu1=0;imu1<nmu;imu1++))
		  do
                      mu1=${list_mu[$imu1]}
                      echo $base_conf/SeqProps/$source_name/$theta_spec/$mu_spec/$f_spec/$theta1/$mu1/prop. $mu1 $theta1 >> $base_3pts/input
		  done
	      done
	      
	      echo $ncombo >> $base_3pts/input
	      
	      list_targ_combo=""
	      for((itheta0=0;itheta0<ntheta;itheta0++))
	      do
		  theta0=${list_theta[$itheta0]}
		  for((itheta1=0;itheta1<ntheta;itheta1++))
		  do
		      theta1=${list_theta[$itheta1]}
		      for((imu0=0;imu0<nmu;imu0++))
		      do
			  mu0=${list_mu[$imu0]}
			  for((imu1=0;imu1<nmu;imu1++))
			  do
			      mu1=${list_mu[$imu1]}
			      
			      iprop0=$(( $itheta0 * $nmu + $imu0 ))
			      iprop1=$(( $itheta1 * $nmu + $imu1 ))
			      
			      targ_combo=$base_3pts/$theta0\_$mu0\_$theta1\_$mu1\_
			      
			      echo $iprop0 $iprop1 $targ_combo >> $base_3pts/input
			      
			      list_targ_combo=$list_targ_combo$targ_combo" "
			  done
		      done
		  done
	      done
              echo $base_conf/Sources/$source_name/source. $base_3pts/2pts_check_P5P5 >> $base_3pts/input

              $MPI_AH_PREF $base_ahmidas/applications/$prog_contr $base_3pts/input
	      
              #let's put together all the micro-correlations needed for the macro correlations
              #put all the needed volume factor, and translate the corr. to the origin
	      for targ_combo in $list_targ_combo
	      do
		  
		  echo $targ_combo
		  
		  for contr in ${two_points_correlations[@]}
		  do
                      awk '
                           {a=a$3" ";b="'$targ_combo'"b$1"_"$2" "}
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
	  
	done #spectator flavor
	
      done #spectator mu
      
    done #spectator theta

done #source