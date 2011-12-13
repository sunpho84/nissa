cat `find . -name "*.c"`|grep "("|grep -v "//"|awk '{a=substr($0,0,1)}a!=" " && a!="\t" && a!= "{" && NF>=2'|grep -v define|sort|uniq|awk -F { '{print $1";"}' > headers.h
