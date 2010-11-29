#!/bin/bash

#set the base folder
source ~/nissa_conf.sh

echo
echo " ##################################################################################"
echo " ##################################################################################"
echo " ####                                                                          ####"
echo " ####                    Three points functions calculation                    ####"
echo " ####                                                                          ####"
echo " ##################################################################################"
echo " ##################################################################################"
echo

#take initial time
tic=$(date +%s)

#count the number of mu and theta
nmu=${#list_mu[@]}
ntheta=${#list_theta[@]}

#create the list of theta for the sequential, whcih is reverted in sign
for((itheta=0;itheta<ntheta;itheta++))
do
    list_seq_theta[$itheta]=$(echo ${list_theta[$itheta]}|awk '{c=substr($1,0,1);if(int(c)==c)ini=1;else ini=2;if($1>0)s="-";print ini,s substr($1,ini)}'
done

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

    #loops over the spectator
    for itheta_spec in ${list_itheta_spec[@]}
    do

      theta_spec=${list_theta[$itheta_spec]}

      for imu_spec in ${list_imu_spec[@]}
      do

	mu_spec=${list_mu[$imu_spec]}
      
	for r_spec in ${list_r_spec[@]}
	do
	
	  #choose the rS1 flavor to have the charged combo
	  if [ $r_spec == 0 ]
	  then
	      rS1=1
	  else
	      rS1=0
	  fi

	  echo "Generating sequential source from the original named: "$source_name" of type: "$source_type" with pars: "${source_pars[@]}" and seed: "$source_seed", for the "$theta_spec" theta, "$mu_spec" mass, "$r_spec" flavour"
	
	  source_dir=$base_conf/Props/$source_name/$theta_spec/$mu_spec/$r_spec
	  targ_dir=$base_conf/SeqSources/$source_name/$theta_spec/$mu_spec/$r_spec
	  
	  if [ -d "$targ_dir" ]
	  then
	      echo "Sequential Source "$targ_dir" already existing"
	      echo
	  else
	      
	      mkdir -pv $targ_dir
	      
              #Generate the sequential p5 source input file
	      TH=$(( $T / 2 ))

	      if [ $tsource -lt $TH ]
	      then
		  TSlice=$(( $tsource + $TH ))
	      else
		  TSlice=$(( $tsource - $TH ))
	      fi

   	      #invoke the program
	      echo $MPI_AH_PREF $base_nissa/Appretto/tools/generate_meson_seqsource/generate_meson_seqsource $L $T $TSlice $source_dir/prop $targ_dir/source $IO_prec
	      $MPI_AH_PREF $base_nissa/Appretto/tools/generate_meson_seqsource/generate_meson_seqsource $L $T $TSlice $source_dir/prop $targ_dir/source $IO_prec
	      
	  fi

          #take time 
	  tac=$tic
	  tic=$(date +%s)
	  echo "Time to generate the sequential source: "$(($tic-$tac)) >> $base_conf/time_log

	  echo "######################## SECOND STEP: Inversion of doublet ############################"
	  echo
	  echo "Second inversion"
	  echo 	  

	  #this is the list of theta to use on the sequential line
	  for theta_seq in ${list_seq_theta[@]}
	  do
	      
	      base_inv=$base_conf/SeqProps/$source_name/$theta_spec/$mu_spec/$r_spec/$theta_seq/
	      mkdir -pv $base_inv
	      
	      echo "Inverting: "$source_name" "$theta_seq
	      
	      if [ ! -f $base_inv/completed ]
	      then
	  	  
                  #create additional masses list: first mass is included as base
		  for mu1 in ${list_mu[@]};do echo $mu1;done|awk 'NR>1{printf("%.12f\n",2*$1*'$kappa')}' > $base_inv/extra_masses.input
		  
                  #prepare input for inverter
		  mu1=${list_mu[0]} #controlla
		  two_kappamu=$(echo $kappa|awk '{printf("%.12f\n",2*$1*'$mu1')}')
		  cat $base_protos/cgmms.input|sed '
                    s|SED_Conf|'$conf'|;
                    s|SED_NL|'$L'|;
                    s|SED_NT|'$T'|;
                    s|SED_ThetaX|'$theta_seq'|;
                    s|SED_ThetaY|'$theta_seq'|;
                    s|SED_ThetaZ|'$theta_seq'|;
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
		      echo $r_spec   #if r_spec=0, don't save the up
		      echo $rS1      #and as rS1=!r_spec, save the down
		      echo $IO_prec  #and vice versa
		  ) > $base_inv/cg_mms_instructions

		  OLD=$PWD
		  cd $base_inv
		  
                  #link configuration and source
		  ln -vfs $base_conf/Conf Conf.$conf
		  for i in $(seq -f%02.0f 00 $last_prop_index)
		  do
		      ln -vsf $base_conf/SeqSources/$source_name/$theta_spec/$mu_spec/$r_spec/source.$i source.$conf.00.$i
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
		  for((imu1=0;imu1<${#list_mu[@]};imu1++))
		  do
		      
		      mu1=${list_mu[$imu1]}
		      mkdir -pv $base_inv/$mu1
		      
		      for ics in $(seq -f%02.0f 00 $last_prop_index)
		      do
			  
			file=$(echo $imu1|awk '{printf("'$base_inv/source.$conf.00.$ics.cgmms."%02d.inverted."$rS1"\n"'",$1)}')
			dest=$base_inv/$mu1/prop.$ics
			  
			if [ ! -f $file ]
			then
			    echo "Could not find: "$file
	                    exit
			else
			    mv -v $file $dest
			fi
			
		      done
		      
		      
		  done #loop over output mass
		  
		  touch completed
		  
		  cd $OLD
		  
	      else
		  echo "Already inverted"
		  echo
	      fi
	      
              #take time
	      tac=$tic
	      tic=$(date +%s)
	      echo "Time to invert sequential source: "$(($tic-$tac)) >> $base_conf/time_log
	      
	  done #loop over "sequential" theta_seq
	  
	  echo "######################## THIRD STEP: Three point function calculation ############################"
	  echo
	
	  if [ $source_type == "Point12" ]
	  then
              vol_fact=1
              prog_contr=applications/contract_meson_2pts
	  elif [ $source_type == "Wall4" ]
	  then
              prog_contr="Appretto/projects/meson_2pts/meson_2pts"
	  elif [ $source_type == "Wall1" ]
	  then
              vol_fact=$(( $L * $L * $L ))
              prog_contr=applications/contract_meson_ultrastochastic_2pts
	  fi

	  base_3pts=$base_conf/3pts/$source_name\_$theta_spec\_$mu_spec\_$r_spec
	  
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

	      echo "L "$L  > $base_3pts/input
              echo "T "$T >> $base_3pts/input
              echo "TWall "$tsource >> $base_3pts/input
	      echo "NPropFirstList "$nprop >> $base_3pts/input #list of quark propagators to be reverted, which close the S1
              for theta1 in ${list_theta[@]}
              do
		  for mu1 in ${list_mu[@]}
		  do #to make the charged, take the same r of the spectator
		    echo " "$base_conf/Props/$source_name/$theta1/$mu1/$r_spec/prop $mu1 $theta1 0 $r_spec >> $base_3pts/input
		  done
	      done
	      echo "NPropSecondList "$nprop >> $base_3pts/input #list of sequential propagators
              for theta2 in ${list_seq_theta[@]}
              do
		  for mu2 in ${list_mu[@]}
		  do
		    echo " "$base_conf/SeqProps/$source_name/$theta_spec/$mu_spec/$r_spec/$theta2/$mu2/prop $mu2 $theta2 2 $rS1 >> $base_3pts/input
		  done
              done
	      echo "Ncontr "$nmicro >> $base_3pts/input
	      cat $base_3pts/micro_correlations >> $base_3pts/input
	      echo "Output "$base_3pts"/three_points_contractions" >> $base_3pts/input

              echo
              echo "Launching program: "$prog_contr
              echo
              $MPI_TM_PREF $base_nissa/$prog_contr $base_3pts/input
	      
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