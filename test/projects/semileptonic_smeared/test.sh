#!/bin/bash
echo "Starting test: semileptonic_smeared"
rm -frv out_source? > log
$* ../../../build_pc/bin/semileptonic_smeared input >> log

err=$(
    for i in 2pts_00_00 2pts_00_30 2pts_30_00 2pts_30_30 3pts_sp0_00_00 3pts_sp0_00_30 3pts_sp0_30_00 3pts_sp0_30_30 3pts_sp1_00_00 3pts_sp1_00_30 3pts_sp1_30_00 3pts_sp1_30_30
    do
	echo B $i $(paste out_sourceB/$i comparison_out_sourceB/$i|grep -v "#"|awk 'NF==4{d+=($1-$3)**2+($2-$4)**2;n++}END{print sqrt(d/n)}')
    done |tee -a log|awk '{a+=$3}END{print a}'
    )
echo "Error for source B: $err"

echo "Test finished"