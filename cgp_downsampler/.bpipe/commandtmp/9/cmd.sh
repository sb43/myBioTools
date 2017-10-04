#!/bin/bash
set -exu
   set -o pipefail
(mkdir -p ~/sb43_scratch116/rt_user/downsampling_564469/scripts/cgp_downsampler/downsampledBams &&             DS_FACTOR=$( awk -v bases=3137454505 '($0!~/#/) {sum=sum+$9} END { print (17+1)/(sum/bases)}' ../../downsample_input/CAL-27_E16.A20sH12qE11.bam.bas) && 					 SAMPLE=$(awk '(NR==2) {print $2}' ../../downsample_input/CAL-27_E16.A20sH12qE11.bam.bas) && 					 /software/CGP/canpipe/test/bin/canpipe_test bamdownsamplerandom seed=3137454505 p=$DS_FACTOR I=../../downsample_input/CAL-27_E16.A20sH12qE11.bam level=1 O=~/sb43_scratch116/rt_user/downsampling_564469/scripts/cgp_downsampler/downsampledBams/$SAMPLE.bam &&  					 date >CAL-27_E16.A20sH12qE11.do_downsampling.track) > .bpipe/commandtmp/9/cmd.out
