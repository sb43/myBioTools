#!/bin/bash
set -exu
   set -o pipefail
(SAMPLE=$(awk '(NR==2) {print $2}' ../../downsample_input/CAL-27.bam.bas) &&          mkdir -p ~/sb43_scratch116/rt_user/downsampling_564469/scripts/cgp_downsampler/downsampledBams/$SAMPLE && 				 /software/CGP/canpipe/test/bin/canpipe_test bwa_mem.pl -b '-T 30 -Y' -mt 6 -o ~/sb43_scratch116/rt_user/downsampling_564469/scripts/cgp_downsampler/downsampledBams/$SAMPLE -r /lustre/scratch112/sanger/cgppipe/canpipe/live/ref/human/GRCh37d5/genome.fa -t 6  -s $SAMPLE ~/sb43_scratch116/rt_user/downsampling_564469/scripts/cgp_downsampler/downsampledBams/$SAMPLE.bam &&  				 rm ~/sb43_scratch116/rt_user/downsampling_564469/scripts/cgp_downsampler/downsampledBams/$SAMPLE.bam &&  				 date >CAL-27.bwa_mem.track) > .bpipe/commandtmp/16/cmd.out
