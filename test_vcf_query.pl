use strict;
use Vcf;
my $vcf_file = $ARGV[0];

print "using vcf $vcf_file\n";
my $vcf = Vcf->new(file => $vcf_file);
      $vcf->parse_header();
      $vcf->recalc_ac_an(0);
foreach my $chr_location(1..22) {
  print "$chr_location\n";
  $vcf->open(region => $chr_location);
  while (my $x=$vcf->next_data_array()) {
     #next if $$x[6] ne 'PASS';
     my $location_key="$$x[0]:$$x[1]:$$x[3]:$$x[4]";
     print "$location_key\n";
   }
}
