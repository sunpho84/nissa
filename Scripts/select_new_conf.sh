if [ ! -f "$conf_list_file" ]
then
    echo "Configuration list '"$conf_list_file"' not present!!!"
    exit
fi
pwd
conf_list=$(cat $conf_list_file)
for conf_test in $conf_list
do
    if [ ! -f "$conf_test/analysis_"$analysis_name"_completed" ] $$ [ ! -f "$conf_test/analysis_"$analysis_name"_running" ]
    then
        conf=$conf_test
        echo "Will work on conf: "$conf
	touch "$conf_test/analysis_"$analysis_name"_running"
	trap "rm -f "$conf_test"/analysis_"$analysis_name"_running" EXIT
        break
    fi
done
if [ "$conf" == "" ]
then
    echo "All configurations calculated!"
    exit
fi
