#!/bin/bash

#set the base folder
source ~/nissa_conf.sh

echo
echo "Two points function calculation"
echo "-------------------------------"
echo "Working dir: $PWD"
echo

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
            s|SED_SolverPrecision|'$prec'|
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
	    
	    mu=${list_mu[$im]}
	    
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
	    
	    targ=$base_conf/Props/$source_name/$theta/$mu/

	    echo "Reconstructing $source_name $theta $mu doublet"
	    
	    cd $targ
            
            #link the configuration
            ln -vfs $base_conf/Conf Conf0
            
            cat $base_protos/reconstruct_doublet_input.xml|sed '
                 s|SED_NL|'$L'|;
                 s|SED_NT|'$T'|;
                 s|SED_ThetaX|'$theta'|;
                 s|SED_ThetaY|'$theta'|;
                 s|SED_ThetaZ|'$theta'|;
                 s|SED_Beta|'$beta'|;
                 s|SED_Kappa|'$kappa'|;
                 s|SED_Mu|'$mu'|;
                 s|SED_Flavored_Source|0|;
                 s|SED_Source_Flavor|0|;
                 s|SED_Last_prop_index|'$last_prop_index'|;
                 ' > $targ/reconstruct_doublet_input.xml
		
            $MPI_AH_PREF $base_ahmidas/example/$reco reconstruct_doublet_input.xml
            
            #move the two flavour to appropriate folder
            for f in 0 1
            do
	      mkdir -pv $targ/$f
	      for ics in $(seq -f%02.0f 00 $last_prop_index)
                do
		mv -v $targ/prop.$ics.$f $targ/$f/prop.$ics
	      done
            done
            
            #remove the input file and the link to conf
            #rm -vf Conf0 reconstruct_doublet_input.xml
            
            cd -
            
            echo "Doublet reconstructed"
            echo
	    
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
      
    for theta1 in ${list_theta[@]}
    do

      if [ "$Standing2pts" == "1" ]
      then
	  list_theta2=${list_theta[0]}
	  echo "Contraction for only standing 2pts"
      else
	  list_theta2=${list_theta[@]}
	  echo "Contraction for all theta combinations of 2pts"
      fi
      
      for theta2 in $list_theta2
      do

	for((imu1=0;imu1<nmu;imu1++))
	do
	
	  for((imu2=imu1;imu2<nmu;imu2++))
	  do
	    
	    mu1=${list_mu[$imu1]}
	    mu2=${list_mu[$imu2]}
	    
	    for s1 in 0 1
	    do
	    
	      for s2 in 0 1
	      do
	      
		echo "Contracting combination:" $mu1 $mu2 $s1 $s2 $theta1 $theta2
	      
		source1=$base_conf/Props/$source_name/$theta1/$mu1/$s1
		source2=$base_conf/Props/$source_name/$theta2/$mu2/$s2
	      
		if [ ! -d $source1 ] || [ ! -d $source2 ]
		then
		    echo "Could not find "$source1" or "$source2
		    exit
		fi
	      
		targ=$base_conf/2pts/$source_name/$theta1/$mu1/$s1/$theta2/$mu2/$s2/
		
		mkdir -pv $targ
		
		(
		    (
			cat $base_protos/contract_two_lines_head.xml
			cat $base_analysis/contract_two_lines_middle.xml
			cat $base_protos/contract_two_lines_tail.xml
			)|sed '
                           s|SED_NL|'$L'|;
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
		    rm -fv micro*
		  
		    #let's put together all the micro-correlations needed for the macro correlations
		    #put all the needed volume factor, and translate the corr. to the origin
		    for op in $(cat $base_analysis/correlations_needed)
		    do
		      list_coef=""
		      list_file=""
		      awk '
                         {a=a$3" ";b=b"name_micro_"$1"_"$2" "}
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
		
	      done
	    
	    done
	    
	  done
	  
	done
      
      done
      
    done
    
done