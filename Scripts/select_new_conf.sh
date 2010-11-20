if [ ! -f conf_list ]
then
    echo "Configuration list 'conf_list' not present!!!"
    exit
fi
conf_list=$(cat conf_list)
for conf_test in $conf_list
do
    if [ ! -d $conf_test ]
    then
        conf=$conf_test
        echo "Will work on "$conf
        break
    fi
done
if [ "$conf" == "" ]
then
    echo "All configuration calculated!"
    exit
fi
