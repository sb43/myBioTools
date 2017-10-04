
my($bam_file)=$ARGV[0];

my $data=`samtools view -H  $bam_file`;

my @lines=split("\n",$data);
my $count=0;
my $rg_hash={'SM'=>'NA','PL'=>'NA', 'DS'=>'NA'};
foreach my $line(@lines) {
 if($line=~/^\@RG/) {
	last if $count > 0;
  my(@record)=(split /\t/, $line);
  foreach my $rg (@record){
    my($tag,$val)=split(/:/,$rg);
    $rg_hash->{$tag}=$val;
  }
  $count++;
 }
}

print "#SM:$rg_hash->{'SM'}\tPL:$rg_hash->{'PL'}\tDS:$rg_hash->{'DS'}\tBAM:$bam_file\n";


