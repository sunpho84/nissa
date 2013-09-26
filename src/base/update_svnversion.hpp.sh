#!/bin/bash

#take new version
new_svnversion=$(svnversion)

#take old version
old_svnversion=$(if [ -f svnversion.hpp ];then sed 's|\"||g' svnversion.hpp|awk '{print $NF}';else echo 0;fi)

#compare and decide if to update
if [ "$new_svnversion" != exported ] && [ "$new_svnversion" != "$old_svnversion" ]
then
    echo \#define SVN_VERSION \""$new_svnversion"\" > svnversion.hpp
fi 
