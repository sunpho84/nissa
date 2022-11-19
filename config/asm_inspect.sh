#!/bin/bash

if [ ! -z "$1" ]
then
    
    cat - > $1
    rm -fr $1.?*
    awk 'function escape_pattern(pat){safe=pat;gsub(/[][^$.*?+{}\\()|]/, "\\\\&", safe);return safe}
    	{b=($1=="#Bookmark")}
	(b && $0 ~ /BEGIN/){a=$NF" "a}
	{split(a,ar," ");for(i in ar) print $0 > "'$1'."ar[i]}
	(b && $0 ~ /END/){gsub(escape_pattern($NF)" ","",a)}' $1
else
    echo "Use $0 file"
    exit
fi
