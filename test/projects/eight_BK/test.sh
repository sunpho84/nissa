#!/bin/bash
echo "Starting test: smeared_BK_all_in_one"
rm -frv out? > log
$* ../../../build_bgq_emu/bin/smeared_BK_all_in_one input >> log
for AB in A B
do
    err=$(
	for i in otto_source$AB\_00_00 otto_source$AB\_00_08 otto_source$AB\_08_00 otto_source$AB\_08_08 twopoints_source$AB\_00_00 twopoints_source$AB\_00_08 twopoints_source$AB\_08_00 twopoints_source$AB\_08_08
	do
	    echo $AB $i $(paste out$AB/$i comparison_out/$i|grep -v "#"|awk 'NF==4{d+=($1-$3)**2+($2-$4)**2;n++}END{print d/n}')
	done|tee -a log|awk '{a+=$3}END{print a}'
	)
    echo "Error for source $AB: $err"
done
echo "Test finished"