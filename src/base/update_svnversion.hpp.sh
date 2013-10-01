#!/bin/bash

#take old version
old_svnversion=$(if [ -f svnversion.hpp ];then sed 's|\"||g' svnversion.hpp|awk '{print $NF}';else echo 0;fi)

#take new version
new_svnversion=$(svnversion)
if [ "$new_svnversion" == exported ] || [ "$new_svnversion" == "Unversioned directory" ]
then
    new_svnversion="$old_svnversion"
fi

#compare and decide if to update
echo \#define SVN_VERSION \""$new_svnversion"\"

