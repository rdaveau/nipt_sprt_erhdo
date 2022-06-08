#!/usr/bin/perl
use strict;
use warnings FATAL => qw[ numeric uninitialized ];
use Getopt::Long qw[ GetOptions ];

sub usage {
	my $args = shift;
	print "Usage: perl $0\n";
	print "--input-dir\n";
	print "--input-file\n";
	print "--output-file\n";
	if(defined($args)){
		print "ERROR: $args\n";
		exit 1;
	}
	exit 0;
}
my %opts = ();
GetOptions(\%opts, 'input-dir=s', 'input-file=s', 'output-file=s');

# USAGE: perl sprt_summary_rxli-2.pl \
# 	--input-dir ~/projects/moabi/nipt/data/HM00008/sprt_summary/v2 \
# 	--input-file  ~/projects/moabi/nipt/data/HM00008/HM00008-1.tsv \
# 	--output-file ~/projects/moabi/nipt/data/HM00008/sprt_summary/v2/HM00008_sprt.tsv

my $i = 0;
my %h = ();
my($d, @d, @p);
my @sprt = qw[ #SNP4 #FWD #REV #MUX %FWD %REV PC2 NC2 PC3 NC3 CS2 CS3 BS ];

open IN, $opts{'input-file'} or die $!;

chomp(my $header = readline(IN));
my @header = split /\t/, $header;

while(<IN>){
	$i++;
	chomp;
	@_ = split /\t/, $_;
	$h{$i}->{$_} = shift @_ foreach @header;
}
close IN;

for $i (1 .. $i){
	@p = ();
	@d = ('sprt');
	unshift @d, $h{$i}->{$_} foreach qw[ fam_id design ];
	push @d, $h{$i}->{'analysis_id'};
	push @p, $_ . $h{$i}->{$_} foreach qw[ NDP PDP HMZ HTZ SNP DIST ALPHA BETA ];
	push @d, join('_', @p), "r-out";
	$d = join('/', @d);

	open IN, "sed 1d $opts{'input-dir'}/${d}/$h{$i}->{fam_id}_FF.TABLE.tsv |" or die $!;
	chomp($a = readline(IN));
	close IN;
	@p = split /\t/, $a;
	$h{$i}->{'FF_' . $_} = shift @p foreach qw[ PAF MAF ZFXY ];

	open IN, "sed 1d $opts{'input-dir'}/${d}/$h{$i}->{fam_id}_COV.TABLE.tsv |" or die $!;
	chomp($a = readline(IN));
	close IN;
	@p = split /\t/, $a;
	$h{$i}->{'DP_' . $_} = shift @p foreach qw[ ZFX ZFY DPNI ];

	open IN, "grep $h{$i}->{target_gene} $opts{'input-dir'}/${d}/$h{$i}->{fam_id}_SPRT.TABLE.tsv | cut -f 2- |" or die $!;
	chomp($a = readline(IN));
	close IN;
	@p = split /\t/, $a;
	$h{$i}->{$_} = shift @p foreach @sprt;
}

unshift @sprt, $_ foreach qw[ DP_DPNI DP_ZFY DP_ZFX FF_ZFXY FF_MAF FF_PAF ];
push @header, @sprt;
open OUT, ">$opts{'output-file'}" or die $!;
print OUT join("\t", @header), "\n";
foreach $i (sort{ $a <=> $b } keys %h){
	@p = ();
	push @p, $h{$i}->{$_} foreach @header;
	print OUT join("\t", @p), "\n";
}
close OUT;