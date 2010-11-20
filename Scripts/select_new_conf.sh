if [ ! -f "$conf_list_file" ]
then
    echo "Configuration list '"$conf_list_file"' not present!!!"
    exit
fi
pwd
conf_list=$(cat $conf_list_file)
for conf_test in $conf_list
do
    if [ ! -f "$conf_test/analysis_"$analysis_name"_completed" ] && [ ! -f "$conf_test/analysis_"$analysis_name"_running" ]
    then
        conf=$conf_test
        echo "Will work on conf: "$conf
        break
    fi
done
if [ "$conf" == "" ]
then
    echo "All configurations calculated!"
    exit
fi

#set the directory name and create it
base_conf=$base_analysis/$conf
mkdir -pv $base_conf

#take posses of the filder
touch "$conf_test/analysis_"$analysis_name"_running"
trap "rm -f "$conf"/analysis_"$analysis_name"_running" EXIT

#link the conf
cd $base_conf
ln -sf $source_confs/conf.$conf Conf
