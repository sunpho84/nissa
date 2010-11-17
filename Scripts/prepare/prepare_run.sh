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

    if [ "$WHERE" == "JUGENE" ]
    then
	cp $base_scripts/proto/work_jugene_header.sh work_header.sh
    else
	touch work_header.sh
    fi

    echo "Adapt the file 'run_instructions.sh' and 'work_header.sh' (if needed) and rerun"
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

	cat ../work_header.sh

	echo "#!/bin/bash"
	echo
	echo "confno="$conf
	echo "additive_seed="$additive_seed
	echo "additive_time_offset="$additive_time_offset
	echo
	echo "source ../analysis.sh"
	echo
    ) > work.sh
    cd -
done

sed 's|SED_L|'$L'|;s|SED_T|'$T'|;s|SED_Base_analysis|'$PWD'|' $base_scripts/proto/proto_analysis.sh > analysis.sh

echo
echo "Now adapt the 'analysis.sh' and delete 'run_instructions.sh'"
echo