#!/bin/bash

if [ "$1" == "" ]
then
    echo Use $0 file
    exit
fi

HFILE=$(echo $1|sed 's|cpp|h|')
MACRO=_$(echo $(basename $HFILE)|tr [:lower:] [:upper:]|sed 's|\.|_|')

echo "#ifndef $MACRO"
echo "#define $MACRO"
if [ -f "$HFILE" ]
then
    grep include $HFILE
fi

grep "(" $1|grep -v "//"|grep -v "basetype"|awk '{a=substr($0,0,1)}a!=" " && a!="\t" && a!= "{" && NF>=2'|grep -v define|grep -v pragma|grep -v main|sort|uniq|awk -F { '{print $1";"}' 
echo "#endif"