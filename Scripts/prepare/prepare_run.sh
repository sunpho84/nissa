#!/bin/bash

source ~/nissa_conf.sh
echo

if [ ! -f run_instructions.sh ]
then
    (
	echo "L="
	echo "T="
	echo "source_confs= #ABSOLUTE path to the configs"
	echo "conf_list=\" \""
	echo "use_external_additive_seed=1 #unless you are doing some test, live it like it"
	echo "use_external_additive_time_offset=1 #idem"
    ) > run_instructions.sh

    echo "Adapt the file 'run_instructions.sh' and rerun"
    echo
    exit
else
    source run_instructions.sh
fi

for tconf in $conf_list
do
    conf=$(printf %.4d $tconf)

    mkdir -pv $conf
    cd $conf
    
    ln -sfv $source_confs/conf.$tconf Conf
    (
	if [ "$use_external_additive_seed" == 1 ]
	then
	    additive_seed=$(($RANDOM%100000))
	else
	    additive_seed=0
	fi
	   
	if [ "$use_external_additive_time_offset" == 1 ]
	then
	    additive_time_offset=$(($RANDOM%$T))
	else
	    additive_time_offset=0
	fi

	echo "#!/bin/bash"
	echo "confno="$conf
	echo "additive_seed="$additive_seed
	echo "additive_time_offset="$additive_time_offset
    ) > gauge_parameters.sh
    cd -
done

sed 's|SED_L|'$L'|;s|SED_T|'$T'|;s|SED_Base_analysis|'$PWD'|' $base_scripts/proto/proto_analysis_parameters.sh > analysis_parameters.sh

echo
echo "Now adapt the 'analysis_parameters.sh' and delete 'run_instructions.sh'"
echo