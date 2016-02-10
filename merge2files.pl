#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use File::Spec;
use Getopt::Long;
use Pod::Usage qw(pod2usage);
use Data::Dumper;

my ($options) = option_builder();


if(!-s $options->{'f1'} && !-s $options->{'f2'}) {
	warn "Unable to find files\n";
	exit;
} 

open(my $fh1, $options->{'f1'}) || warn "print $!";
open(my $fh2, $options->{'f2'}) || warn "print $!";




my ($f1,$d1,$s1) = fileparse($options->{'f1'},qr/\.[^.]*/);
my ($f2,$d2,$s2) = fileparse($options->{'f2'},qr/\.[^.]*/);

my $col_f1=$options->{'c1'};
my $col_f2=$options->{'c2'};
$options->{'hdr1'}=0 unless($options->{'hdr1'}); 
$options->{'hdr2'}=0 unless($options->{'hdr2'}); 
chomp $col_f1;
chomp $col_f2;


# create file handle for filed to write [ common columns and different columns
open(my $fh_common, '>', $options->{'o'}.'/'.$f1.'_'.$f2.'common'.$s1); 

# create hash of both the files....

my ($file1_hash)=get_file_hash($fh1,$col_f1,$options->{'hdr1'},$f1);
my ($file2_hash)=get_file_hash($fh2,$col_f2,$options->{'hdr2'},$f2);

#combine files
combine_hash($file1_hash,$file2_hash,$fh_common,$options->{'hdr1'});

# sub to combine file data...
sub combine_hash {
	my($hash_f1,$hash_f2,$fh_combined,$header)=@_;
	my $count=0;
	# Write header lines to  file...
	if($header){
		print $fh_common $hash_f1->{'pre_header_lines'};
		print $fh_common $hash_f1->{'header'}."\t".$hash_f2->{'header'}."\n";
	}
	#print Dumper $hash_f1;
  # write common rows
	foreach my $key(keys %$hash_f1) {
		next if $key eq 'header';
		next if $key eq 'pre_header_lines';
		if ($hash_f2->{$key}) {
		  print $fh_combined $hash_f1->{$key}."\t".$hash_f2->{$key}."\n";
		  $count++;
		}
	}	
	print "Found $count matching lines\n"; 
}

# create hash from file
sub get_file_hash {
	my ($fh,$columns,$header_row,$f)=@_;
	my $file_hash;
	my $count=0;	
	while (my $line=<$fh>) {	
		$count++;
		next if $line=~/^\n/;
		$line=~s/\n//g;
		my ($col_key) = get_col_key($line,$columns) if $count >= $header_row;
		if($count <= $header_row ) {
			if($count==$header_row) {
			 $file_hash->{'header'}=$line;
				print "Columns to compare from file: : $col_key\n" if defined $col_key; 
			}
			else {
				$file_hash->{'pre_header_lines'}.="$line\n";
			}
			next;
		} 
		$file_hash->{$col_key}="$line";
	}

return $file_hash


}


# get common key from both files
sub get_col_key {
	my ($line,$columns)=@_;
	my $col_key;
	
                my @col=split(',',$columns);
                foreach my $col_num(@col) {
                        my($field)=(split "\t", $line)[$col_num - 1];
                        $col_key.='_'.$field;
                }	
return $col_key;	
}



1;


#########################################

sub option_builder {
	my ($factory) = @_;

	my %opts;

	&GetOptions (
					'h|help'    => \$opts{'h'},
					'f1|file1=s' => \$opts{'f1'},
					'f2|file2=s'  => \$opts{'f2'},
					'c1|col1=s'  => \$opts{'c1'},
					'c2|col2=s'  => \$opts{'c2'},
					'o|output=s'  => \$opts{'o'},
					'hdr1|header1=s'  => \$opts{'hdr1'},
					'hdr2|header2=s'  => \$opts{'hdr2'},
					'v|version'  => \$opts{'v'},
	);

	pod2usage(q{'-f1' file1 must be specified.}) unless(defined $opts{'f1'}) ;
	pod2usage(q{'-f2' file2 must be specified.}) unless(defined $opts{'f2'}) ;
	return \%opts;
}

__END__
=head1 NAME

merge2files.pl - creates a single merged file using user defined matching columns from both files    

=head1 SYNOPSIS

merge2files.pl [-h] -f1 -f2 -c1 -c2 [-hdr1 -hdr2]

  Required Options:
    --help          (-h)  This message and format of input file
     One or more of the following:
    --file1  (-f1) file1 
    --file2  (-f2) file2 
    --col1  (-c1) comma separated column numbers fom file1  
    --col2  (-c2) comma separated columns numbers fom file2 
    --header1  (-hdr1)  header row number first file 
    --header2  (-hdr2)  header row number second file
    --output  (-o) output folder 

=cut


