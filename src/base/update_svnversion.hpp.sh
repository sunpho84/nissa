#!/bin/bash

if [ "$1" == "" ]
then
    echo "Error, use $0 dir"
    exit 1
fi

#take old version
old_svnversion=$(if [ -f svnversion.hpp ];then sed 's|\"||g' svnversion.hpp|awk '{print $NF}';else echo 0;fi)

#take new version
cd $1
new_svnversion=$(svnversion 2>&1)
cd $OLDPWD

#if svn unavailable or not svn downloaded reset to old
if \
    [ "$new_version" == "" ] || \
    [ "$new_svnversion" == exported ] || \
    [ "$new_svnversion" == "Unversioned directory" ]
then
    new_svnversion="$old_svnversion"
fi

#compare and decide if to update
echo \#define SVN_VERSION \""$new_svnversion"\" > test
if [ ! -f svnversion.hpp ]
then
    mv test svnversion.hpp
else
    diff test svnversion.hpp > /dev/null
    if [ "$?" != "0" ]
    then
	mv test svnversion.hpp
    else
	rm test
    fi
fi
