#!/bin/bash

if [ "$1" == "" ]
then
    echo Use $0 file
    exit
fi

grep "(" $1|grep -v "//"|grep -v "basetype"|awk '{a=substr($0,0,1)}a!=" " && a!="\t" && a!= "{" && NF>=2'|grep -v define|grep -v pragma|grep -v main|sort|uniq|awk -F { '{print $1";"}' 
