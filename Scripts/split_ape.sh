#!/bin/bash

if [ ! "$2" ]
then
    echo "Use: "$0" file combo"
    exit
fi

sed 's|(||;s|)||g;s|,| |g' $1 |awk '
$1=="Kcrit"{cart=sprintf("%.4f_%.4f/'$2'",$15,$6);
     system("mkdir -p "cart)}
st{print $0 > cart"/"oss;st--}
$1!=$1+0 && NF==1{oss=$1;st=32}'
