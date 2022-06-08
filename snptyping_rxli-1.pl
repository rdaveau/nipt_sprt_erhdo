#!/usr/bin/perl
use strict;
use warnings FATAL => qw[ numeric uninitialized ];
use Getopt::Long qw[ GetOptions ];
use List::Util qw[ sum min max ];
use File::Basename;

sub usage {
	my $args = shift;
	print "Usage: perl $0\n";
	print "--gender sexe of fetus and proband [used for RXLI]\n";
	print "--proband is proband (un)affected child or relative?\n";
	print "--Mo read pile from mother nuclear dna\n";
	print "--Fa read pile from father nuclear dna [optional]\n";
	print "--CI read pile from proband nuclear dna\n";
	print "--DPNI read pile from mother plasmatic dna\n";
	print "--DPN read pile from fetal nuclear dna\n";
	print "--HMZ minimal allele ratio to call a homozygote snv\n";
	print "--HTZ minimal allele ratio to call a heterozygote snv\n";
	print "--NDP minimal sequencing depth in nuclear dna\n";
	print "--PDP minimal sequencing depth in plasmatic dna\n";
	print "--out output file\n";
	if(defined($args)){
		print "ERROR: $args\n";
		exit 1;
	}
	exit 0;
}
my %opts = (HMZ => 0.85, HTZ => 0.15, NDP => 8, PDP => 15); # Default thresholds for HMZ/HTZ ratios and (N)uclear/(P)lasmatic depth
GetOptions(\%opts, 'gender:s', 'proband=s', 'Mo=s', 'Fa:s', 'CI=s', 'DPNI=s', 'DPN=s', 'HMZ:f', 'HTZ:f', 'NDP:i', 'PDP:i', 'out=s');

my @piles = qw[ Mo CI DPN DPNI ];
foreach (@piles){
	usage("Missing option: <$_>") unless defined $opts{$_};
	-f $opts{$_} or usage("File not found: $opts{$_}");
}
usage("Missing option: <gender>") unless defined $opts{gender};

my @proband = qw[ AC UC AR ];
usage("Unknown proband: <$opts{proband}>") unless $opts{proband} ~~ @proband; # Check if proband is properly set

my @gender = qw[ FF FM MF MM ];
usage("Unknown gender: <$opts{gender}>") unless $opts{gender} ~~ @gender; # Check if gender is properly set

my %snpmap = ();
my $Fa = (($opts{gender} eq 'MM') || ($opts{gender} eq 'MF' && $opts{proband} eq 'AR')) ? 0 : 1;
if($Fa){ # Deprecated
	usage("Missing option: <Fa>") unless defined $opts{Fa};
	-f $opts{Fa} or usage("File not found: $opts{Fa}");
	unshift @piles, 'Fa';
	%snpmap = (
		FF => {
			AC => {
				A => {
					AA => { type => 'AA', rhap => 'A', q0 => 1, q1 => 3 },
					AB => { type => 'DD', rhap => 'B', q0 => 2, q1 => 1 }
				},
				B => {
					BB => { type => 'BB', rhap => 'B', q0 => 1, q1 => 3 },
					AB => { type => 'CC', rhap => 'A', q0 => 2, q1 => 1 }
				}
			},
			UC => {
				A => {
					AA => { type => 'LL', rhap => 'B', q0 => 2, q1 => 1 },
					AB => { type => 'II', rhap => 'A', q0 => 1, q1 => 3 }
				},
				B => {
					BB => { type => 'KK', rhap => 'A', q0 => 2, q1 => 1 },
					AB => { type => 'JJ', rhap => 'B', q0 => 1, q1 => 3 }
				}
			},
			AR => {
				A => {
					AA => { type => 'AA', rhap => 'A', q0 => 1, q1 => 3 },
					BB => { type => 'QQ', rhap => 'B', q0 => 2, q1 => 1 }
				},
				B => {
					AA => { type => 'RR', rhap => 'A', q0 => 2, q1 => 1 },
					BB => { type => 'BB', rhap => 'B', q0 => 1, q1 => 3 }
				}
			}
		},
		FM => {
			AC => {
				A => {
					A => { type => 'EE', rhap => 'A', q0 => 1, q1 => 3 },
					B => { type => 'HH', rhap => 'B', q0 => 2, q1 => 1 }
				},
				B => {
					A => { type => 'GG', rhap => 'A', q0 => 2, q1 => 1 },
					B => { type => 'FF', rhap => 'B', q0 => 1, q1 => 3 }
				}
			},
			UC => {
				A => {
					A => { type => 'PP', rhap => 'B', q0 => 2, q1 => 1 },
					B => { type => 'MM', rhap => 'A', q0 => 1, q1 => 3 }
				},
				B => {
					A => { type => 'NN', rhap => 'B', q0 => 1, q1 => 3 },
					B => { type => 'OO', rhap => 'A', q0 => 2, q1 => 1 }
				}
			},
			AR => {
				A => {
					A => { type => 'EE', rhap => 'A', q0 => 1, q1 => 3 },
					B => { type => 'HH', rhap => 'B', q0 => 2, q1 => 1 }
				},
				B => {
					A => { type => 'GG', rhap => 'A', q0 => 2, q1 => 1 },
					B => { type => 'FF', rhap => 'B', q0 => 1, q1 => 3 }
				}
			}
		},
		MF => {
			AC => {
				A => {
					AA => { type => 'A', rhap => 'A', q0 => 2, q1 => 3 },
					AB => { type => 'D', rhap => 'B', q0 => 2, q1 => 3 }
				},
				B => {
					BB => { type => 'B', rhap => 'B', q0 => 2, q1 => 3 },
					AB => { type => 'C', rhap => 'A', q0 => 2, q1 => 3 }
				}
			},
			UC => {
				A => {
					AA => { type => 'L', rhap => 'B', q0 => 2, q1 => 3 },
					AB => { type => 'I', rhap => 'A', q0 => 2, q1 => 3 }
				},
				B => {
					BB => { type => 'K', rhap => 'A', q0 => 2, q1 => 3 },
					AB => { type => 'J', rhap => 'B', q0 => 2, q1 => 3 }
				}
			}
		}
	);
}else{
	%snpmap = (
		MM => {
			AC => {
				A => { type => 'E/G', rhap => 'A', q0 => 2, q1 => 3 },
				B => { type => 'F/H', rhap => 'B', q0 => 2, q1 => 3 }
			},
			UC => {
				A => { type => 'M/O', rhap => 'B', q0 => 2, q1 => 3 },
				B => { type => 'N/P', rhap => 'A', q0 => 2, q1 => 3 }
			},
			AR => {
				A => { type => 'E/G', rhap => 'A', q0 => 2, q1 => 3 },
				B => { type => 'F/H', rhap => 'B', q0 => 2, q1 => 3 }
			}
		},
		MF => {
			AR => {
				AA => { type => 'A/R', rhap => 'A', q0 => 2, q1 => 3 },
				BB => { type => 'B/Q', rhap => 'B', q0 => 2, q1 => 3 }
			}
		}
	);
}

my %DP = (DPNI => $opts{PDP});
$DP{$piles[$_]} = $opts{NDP} for 0 .. $#piles - 1;

my($r, $k);
my(@pile, @call, @n, @g, @ab);
my @desc = qw[ GT TRC RRC ARC BRC ];
my(%pile, %call, %r);

sub genotype {

	my $dp = shift; # Minimal DP threshold
	@pile = split /\t/, shift; # Pile readline
	$pile{$_} = shift @pile foreach qw[ POS REF TRC A T C G N Ins Del ];

	return() if $pile{TRC} < $dp # Skip low covered sites
	|| $pile{N} != 0 # Skip N's
	|| $pile{Ins} =~ /:/ # Skip insertions
	|| $pile{Del} =~ /:/; # Skip deletions

	$r = 0;
	@$_ = () foreach (\@n, \@g);

	foreach (qw[ A T C G ]){
		$pile{$_ . 'RR'} = $pile{$_} / $pile{TRC}; # Allele ratio
		$r += $pile{$_ . 'RR'}; # Cumulative allele ratio
		push @n, $pile{$_ . 'RR'} <= $opts{HTZ} ? 0 : (($pile{$_ . 'RR'} >= .5 - $opts{HTZ}) & ($pile{$_ . 'RR'} <= .5 + $opts{HTZ})) ? 1 : ($pile{$_ . 'RR'} >= $opts{HMZ}) ? 2 : 3; # Allele copy number (n=3 for ambiguous allele ratio)
		push @g, ($_) x $n[-1];
	}
	return() if $r != 1 || $#g != 1; # Skip remaining indels and ambiguous genotype

	@ab = grep(!/$pile{REF}/, @g);
	($pile{GT}, $pile{ALT}) = ($#ab < 0 ? 'AA' : grep(/1/, @n) ? 'AB' : 'BB', shift @ab);

	if(exists $call{$pile{POS}}->{ALT}){ # Used for properly defining ARC and BRC in DPNI[AA] for X-linked SNP4 genotyping
		@call = keys %{$call{$pile{POS}}->{ALT}};
		$pile{ALT} = shift @call;
	}

	$pile{BRC} = $pile{TRC} - ($pile{RRC} = $pile{$pile{REF}}) - ($pile{ARC} = defined $pile{ALT} ? $pile{$pile{ALT}} : 0);

	# Estimate bakground read ratio/enrichment
	$pile{BRE} = 0;
	if($pile{BRC} != 0 && $pile{ARC} != 0){
		$pile{BRR} = $pile{BRC} / $pile{TRC};
		$pile{BRE} = $pile{BRR} / sum($pile{ARC} / $pile{TRC}, $pile{BRR});
	}
	next if $pile{BRE} >= 0.3; # Skip multi-allelic/ambiguous sites

	return(\%pile);
}

open MAF, ">$opts{out}_MAF.tsv" or die "Can't open file: $opts{out}_MAF.tsv\n";
print MAF join("\t", 'CHR', 'POS', 'REF', 'ALT', 'GT', 'TRC', 'RRC', 'ARC', 'BRC'), "\n";

open PAF, ">$opts{out}_PAF.tsv" or die "Can't open file: $opts{out}_MAF.tsv\n";
print PAF join("\t", 'CHR', 'POS', 'REF', 'ALT', 'Fa', 'TRC', 'RRC', 'ARC', 'BRC'), "\n";

foreach my $chr (1 .. 22){

	print "Processing contig: $chr ...\n";

	%call = ();
	foreach my $sample (qw[ Mo Fa ]){
		open IN, "sed 1d $opts{$sample} | cut -f 2- | grep -e '^${chr}[[:space:]]' | cut -f 2- |" or die "Can't open file: $opts{$sample}\n";
		while(<IN>){

			$r = &genotype($opts{NDP}, $_);
			next if !defined($r);
			next if $$r{GT} eq 'AB';
			$call{$$r{POS}}->{REF} = $$r{REF};
			$call{$$r{POS}}->{ALT}->{$$r{ALT}}++ if defined $$r{ALT};
			$call{$$r{POS}}->{$sample}->{$_} = $$r{$_} foreach @desc;
		}
		close IN;
	}

	foreach $k (keys %call){
		delete $call{$k} and next if((!exists $call{$k}->{Mo}) || (!exists $call{$k}->{Fa}));
		delete $call{$k} if($call{$k}->{Mo}->{GT} eq $call{$k}->{Fa}->{GT});
	}

	open IN, "sed 1d $opts{DPNI} | cut -f 2- | grep -e '^${chr}[[:space:]]' | cut -f 2- |" or die "Can't open file: $opts{DPNI}\n";
	while(<IN>){

		$r = &genotype($opts{PDP}, $_);
		next if !defined($r);

		if(exists $call{$$r{POS}}->{REF}){
			$call{$$r{POS}}->{DPNI}->{$_} = $$r{$_} foreach @desc;
		}

		next if !$$r{ARC};
		$$r{ALT} = (defined $$r{ALT}) ? ($$r{ALT}) : ('NA');

		@call = ($chr);
		push @call, $$r{$_} foreach qw[ POS REF ALT GT TRC RRC ARC BRC ];
		print MAF join("\t", @call), "\n";
	}
	close IN;

	foreach $k (sort{ $a <=> $b } keys %call){
		delete $call{$k} and next if !exists $call{$k}->{DPNI};
		@call = keys %{$call{$k}->{ALT}};
		delete $call{$k} and next if $#call != 0; # Skip multi-allelic sites
		$call{$k}->{ALT} = $call[0];
		@call = ($chr, $k);
		push @call, $call{$k}->{$_} foreach qw[ REF ALT ];
		push @call, $call{$k}->{Mo}->{GT} eq 'BB' ? 'R' : 'A'; # Paternal-specific allele (Ref|Alt)
		push @call, $call{$k}->{DPNI}->{$_} foreach qw[ TRC RRC ARC BRC ];
		print PAF join("\t", @call), "\n";
	}

}
close MAF;
close PAF;

# X-linked SNP4 genotyping
%call = ();

print "Processing contig: X ...\n";

foreach my $sample (@piles){

	print "\tProcessing pile: $sample ...\n";

	open IN, "sed 1d $opts{$sample} | cut -f 2- | grep -e '^X[[:space:]]' | cut -f 2- |" or die "Can't open file: $opts{$sample}\n";
	while(<IN>){

		$r = &genotype($DP{$sample}, $_);
		next if !defined($r);

		$call{$$r{POS}}->{REF} = $$r{REF};
		$call{$$r{POS}}->{ALT}->{$$r{ALT}}++ if defined $$r{ALT};
		$call{$$r{POS}}->{$sample}->{$_} = $$r{$_} foreach @desc;
	}
	close IN;
}

foreach $k (keys %call){
	foreach (@piles){
		delete $call{$k} and last unless exists $call{$k}->{$_};
	}
	next unless exists $call{$k};
	delete $call{$k} and next if $call{$k}->{Mo}->{GT} ne $call{$k}->{DPNI}->{GT};
	delete $call{$k} and next if $call{$k}->{Mo}->{GT} ne 'AB';
	if(exists $call{$k}->{ALT}){
		@call = keys %{$call{$k}->{ALT}};
		delete $call{$k} and next if $#call != 0; # Skip multi-allelic sites
		$call{$k}->{ALT} = $call[0];
	}else{
		$call{$k}->{ALT} = 'NA';
	}
}

my @HMZ = ();

push @HMZ, 'DPN' if $opts{gender} =~ /^M/; # Male fetus
push @HMZ, 'CI' if $opts{gender} =~ /M$/; # Male proband
push @HMZ, 'Fa' if $Fa;

if($#HMZ >= 0){
	foreach (@HMZ){
		foreach $k (keys %call){
			delete $call{$k} and next if $call{$k}->{$_}->{GT} eq 'AB'; # Skip HTZ SNP
			$call{$k}->{$_}->{GT} = substr($call{$k}->{$_}->{GT}, 0, 1); # Translate AA/BB into A/B
		}
	}
}

if($Fa){
	foreach $k (keys %call){
		$call{$k}->{$_} = $snpmap{$opts{gender}}->{$opts{proband}}->{$call{$k}->{Fa}->{GT}}->{$call{$k}->{CI}->{GT}}->{$_} foreach qw[ type rhap q0 q1 ];
	}
}else{
	foreach $k (keys %call){
		$call{$k}->{$_} = $snpmap{$opts{gender}}->{$opts{proband}}->{$call{$k}->{CI}->{GT}}->{$_} foreach qw[ type rhap q0 q1 ];
	}
}

my @header = qw[ CHR POS REF ALT SNP.TYPE Q0.TYPE Q1.TYPE RISK.HAP ];
push @header, "${_}.GT" foreach @piles;
foreach $a (@piles){
	foreach $b (qw[ TRC RRC ARC BRC ]){
		push @header, "${a}.${b}";
	}
}
open OUT, ">$opts{out}_chrX.tsv" or die "Can't open file: $opts{out}_chrX.tsv\n";
open ERR, ">$opts{out}_undef.tsv" or die "Can't open file: $opts{out}_undef.tsv\n";

print OUT join("\t", @header), "\n";
print ERR join("\t", qw[ CHR POS REF ALT ], @piles), "\n";
foreach $k (sort { $a <=> $b } keys %call){
	@call = ('X', $k, $call{$k}->{REF}, $call{$k}->{ALT});
	if(!defined $call{$k}->{type}){
		push @call, $call{$k}->{$_}->{GT} foreach @piles;
		print ERR join("\t", @call), "\n";
	}else{
		push @call, $call{$k}->{$_} foreach qw[ type q0 q1 rhap ];
		push @call, $call{$k}->{$_}->{GT} foreach @piles;
		foreach $a (@piles){
			foreach $b (qw[ TRC RRC ARC BRC ]){
				push @call, $call{$k}->{$a}->{$b};
			}
		}
		print OUT join("\t", @call), "\n";
	}
}
close OUT;
close ERR;
my $root_dir = dirname($opts{out});
`touch ${root_dir}/done`;