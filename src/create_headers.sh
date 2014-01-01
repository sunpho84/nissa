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

egrep ' +[[:alnum:]_]+ [[:alnum:]_]+\([[:alnum:]_, &*]+\)' $1|grep -v "//"|sort|uniq|awk -F { '{print $1";"}' 
echo "#endif"
