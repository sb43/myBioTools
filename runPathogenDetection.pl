#!/software/perl-5.16.3/bin/perl 
##########LICENCE############################################################
# Copyright (c) 2016 Genome Research Ltd.
# 
# Author: Cancer Genome Project cgpit@sanger.ac.uk
# 
# This file is part of runPathogenDetection.
# 
# cgpVAF is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the Free
# Software Foundation; either version 3 of the License, or (at your option) any
# later version.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
# details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
##########LICENCE##############################################################
BEGIN {
  use Cwd qw(abs_path);
  use File::Basename;
  $SIG{__WARN__} = sub {warn $_[0] unless(( $_[0] =~ m/^Subroutine Tabix.* redefined/) || ($_[0] =~ m/^Use of uninitialized value \$buf/))};
};

use strict;
use English qw( -no_match_vars );
use Pod::Usage qw(pod2usage);
use Const::Fast qw(const);
use Getopt::Long;
use Try::Tiny qw(try catch finally);
use File::Path qw(mkpath);
use Capture::Tiny qw(:all);

const my $gottcha_path => qw(/software/CGP/external-apps/GOTTCHA/bin);
#awk command to create FASTQ file
const my $awk_cmd => '{print "@"$1"\n"$10"\n+\n"$11}';
const my $BACTERIA => '/lustre/scratch112/sanger/cgppipe/canpipe/test/ref/human/38/gottcha_db/database/GOTTCHA_BACTERIA_c4937_k24_u30_xHUMAN3x.strain';
const my $VIRUSES => '/lustre/scratch112/sanger/cgppipe/canpipe/test/ref/human/38/gottcha_db/database/GOTTCHA_VIRUSES_c5900_k24_u30_xHUMAN3x.strain';

try {
  my ($options) = option_builder();
  my $cmd=undef;	

  if(defined $options->{'i'} && lc($options->{'i'}) eq 'y' ) {
		print "Feteching all unmapped reads and their mates ......\n";
		#unmapped read whose mate is mapped
		$cmd="samtools view -f 4 -F264 $options->{'bam'} | awk \'$awk_cmd\'  >$options->{'o'}/tmp/unmapped.fastq";
		_run_cmd($cmd);
		#mapped read whose mate is unmapped
		$cmd="samtools view -f 8 -F260 $options->{'bam'} | awk \'$awk_cmd\'  >>$options->{'o'}/tmp/unmapped.fastq";
		_run_cmd($cmd);
		#Both reads unmapped
		my $cmd="samtools view -f 12 -F256 $options->{'bam'} | awk \'$awk_cmd\' >>$options->{'o'}/tmp/unmapped.fastq";
		_run_cmd($cmd);
  }else {
		print "Feteching only unmapped reads ......\n";
		$cmd="samtools view -f 4 $options->{'bam'} | awk \'$awk_cmd\'  >$options->{'o'}/tmp/unmapped.fastq";
		_run_cmd($cmd);
	} 
  print "Running detection for : $options->{'t'} .....\n";
	my $gottcha_db = $VIRUSES;
  if($options->{'t'} eq 'BACTERIA'){
	 $gottcha_db=$BACTERIA;
  }
  if(defined $options->{'d'}){
		print "Using user defined fasta file as detection database";
		$gottcha_db=$options->{'d'}; 
	}
  
  $cmd = "$gottcha_path/gottcha.pl --threads $options->{'n'} --mode summary --minQ 10 --outdir $options->{'o'} ".
  " --input $options->{'o'}/tmp/unmapped.fastq".
  " --noPlasmidHit ".
  " --database $gottcha_db";
  # run command
  _run_cmd($cmd) if(defined $cmd);
  
  print "\nCompleted virus detection using GOTTCHA.....\nCheck result file unmapped.gottcha.tsv .... \n";
}
catch {
print "$_";
};


sub option_builder {
  my ($factory) = @_;
  my %opts;
  &GetOptions (
          'h|help'    => \$opts{'h'},
          'bam|sampleBam=s' => \$opts{'bam'},
          'i|includeMate=s' => \$opts{'i'},
          't|analysisType=s' => \$opts{'t'},
          'd|userDb=s' => \$opts{'d'},
          'n|numCpu=i' => \$opts{'n'},
          'o|outdir=s'  => \$opts{'o'},
  );

  pod2usage(-verbose => 1) if(defined $opts{'h'});
  if(!defined $opts{'t'}) {
    $opts{'t'}='VIRUSES';
  }
  if(!defined $opts{'n'}) {
    $opts{'n'}=1;
  }
	if(uc($opts{'t'}) ne 'VIRUSES' && uc($opts{'t'}) ne 'BACTERIA') {
		print "\nAnalysis type doesn't match with VIRUS or BACTERIA\n";
	}
	if(defined $opts{'o'}) {
		mkpath("$opts{'o'}/tmp");
	}
  pod2usage(q{'-bam' bam file must be specified.}) unless(defined $opts{'bam'}) ;
  pod2usage(q{'-t' analysis type[VIRUSES,BACTERIA:default VIRUSES].}) unless(defined $opts{'t'});
  pod2usage(q{'-n' number of threads }) unless(defined $opts{'n'});
  pod2usage(q{'-o' output location must be specified}) unless(defined $opts{'o'});
  return \%opts;
}


=head2 _run_cmd
runs external command
Inputs
=over 2
=item cmd - command to run
=back
=cut

sub _run_cmd {
	my($cmd)=@_;
	my ($out,$stderr,$exit)=capture{system($cmd)};
	if($exit) {
		print "Failed to run <<<<<<< \n $cmd  <<<<<< \n with status <<<<<< \n OUT:\n $out  :ERR:\n $stderr EXIT:\n $exit \n <<<<<<< \n";
	}
	return $out;
}


__END__

=head1 NAME

runPathogenDetection.pl - run virus detection a set of unmapped reads extracted from bam file and mapped against GOTTCHA signatures 

=head1 SYNOPSIS

runPathogenDetection.pl -bam -o [-i -t -d -n -h ]

Required Options (bam and outdir must be defined):

  --sampleBam         (-bam) sample bam file 
  --outdir            (-o) outdir [ Path to output directory ]
Optional:

  --includeMate       (-i) include mate sequence of unmapped read [Y,N:default N], an option to include/exclude mapped mate of an unmapped read, takes more time to fetch the reads 
  --databaseType      (-t) databaseType [VIRUSES,BACTERIA:default VIRUSES]
  --userDb            (-d) The path of signature database. The database can be
                           in FASTA format or BWA index (5 files).
  --numCpu            (-n) fasta reference genome file 
  --help              (-h) This message 
	Example:
     perl runPathogenDetection.pl -bam test.bam -o testdir

=cut

