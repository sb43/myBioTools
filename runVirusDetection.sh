#!/bin/bash
set -ue
##################
##Author: sb43   #
##Date:9/02/2016 #
##################
in_bam=${1:-}
out_dir=${2:-}
num_threads=${3:-}

usage ()
{
  echo 'Usage : runVirusDetection.sh <BAM_FILE> <OUTDIR> <number_of_cpu>'
  exit
}

if [[ -z "$in_bam" || -z "$out_dir" || -z "$num_threads" ]]
then
  usage
fi

mkdir -p $out_dir/tmp

echo "Feteching unmapped reads......"
#Unmapped read whose mate is mapped
samtools view -u  -f 4 -F264  $in_bam > $out_dir/tmp/tmp1.bam &
#mapped read whose mate is unmapped
samtools view -u  -f 8 -F260  $in_bam > $out_dir/tmp/tmp2.bam &
#Both reads unmapped
samtools view -u  -f 12 -F256  $in_bam > $out_dir/tmp/tmp3.bam &
# wait till all background processes finish
wait
# merge reads 
echo "Merging unmapped reads......"
samtools merge -u - $out_dir/tmp/tmp[123].bam | samtools sort -n - $out_dir/tmp/unmapped
#create fastq
echo "Creating fastq input for GOTTCHA script......"
bedtools bamtofastq -i $out_dir/tmp/unmapped.bam -fq $out_dir/tmp/unmapped.fastq 

echo " Running GOTTCHA script....."

/software/CGP/external-apps/GOTTCHA/bin/gottcha.pl \
--threads $num_threads \
--mode all \
--minQ 10 \
--outdir $out_dir \
--input $out_dir/tmp/unmapped.fastq \
--database /lustre/scratch112/sanger/cgppipe/canpipe/test/ref/human/38/gottcha_db/database/GOTTCHA_VIRUSES_c5900_k24_u30_xHUMAN3x.strain

echo "Completed virus detection using GOTTCHA....."
echo "Check result file unmapped.gottcha.tsv .... "
