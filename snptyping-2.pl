#!/usr/bin/perl
use strict;
use warnings FATAL => qw[ numeric uninitialized ];
use Getopt::Long qw[ GetOptions ];
use List::Util qw[ sum min max ];
use File::Basename;

sub usage {
	my $args = shift;
	print "Usage: perl $0\n";
	print "--genmode genetic model of disease hereditary\n";
	print "--proband is proband (un)affected child or relative?\n";
	print "--Mo read pile from mother nuclear dna\n";
	print "--Fa read pile from father nuclear dna\n";
	print "--CI read pile from proband nuclear dna\n";
	print "--DPNI read pile from mother plasmatic dna\n";
	print "--HMZ minimal allele ratio to call a homozygote snv\n";
	print "--HTZ minimal allele ratio to call a heterozygote snv\n";
	print "--NDP minimal sequencing depth in nuclear dna\n";
	print "--PDP minimal sequencing depth in plasmatic dna\n";
	print "--out output file prefix\n";
	if(defined($args)){
		print "ERROR: $args\n";
		exit 1;
	}
	exit 0;
}
my %opts = (HMZ => 0.85, HTZ => 0.15, NDP => 8, PDP => 15); # Default thresholds for HMZ/HTZ ratios and (N)uclear/(P)lasmatic depth
GetOptions(\%opts, 'genmode=s', 'proband=s', 'Mo=s', 'Fa=s', 'CI=s', 'DPNI=s', 'HMZ:f', 'HTZ:f', 'NDP:i', 'PDP:i', 'out=s');

foreach (qw[ Mo Fa CI DPNI ]){
	usage("Missing option: <$_>") unless defined $opts{$_};
	-f $opts{$_} or usage("File not found: $opts{$_}");
}
my @proband = qw[ AC UC AR UR ];
usage("Unknown proband: <$opts{proband}>") unless $opts{proband} ~~ @proband; # Check if proband is properly set
my @genmode = qw[ ADMI ADPI ARI ];
usage("Unknown genmode: <$opts{genmode}>") unless $opts{genmode} ~~ @genmode; # Check if genmode is properly set

my @piles = qw[ Mo Fa CI DPNI ];
my @header = qw[ CHR POS REF ALT SNP.TYPE SNP.SUBTYPE RISK.HAP Fa.ALLELE Fa.HAP ];
push @header, "${_}.GT" foreach @piles;
foreach $a (@piles){
	foreach $b (qw[ TRC RRC ARC BRC ]){
		push @header, "${a}.${b}";
	}
}

my($r, $k);
my(@call, @n, @g, @ab);
my(@pile, %pile, %call);
my @chr = 1..22;
my %snpmap = (
	AC => {
		AA => {
			AA => {
				AA => { type => '2', subtype => 'A', rhap => undef, snp3 => undef, hap => undef },
				AB => { type => 'W1', subtype => 'E', rhap => undef, snp3 => undef, hap => undef },
				BB => { type => 'W1', subtype => 'F', rhap => undef, snp3 => undef, hap => undef }
			},
			AB => {
				AA => { type => '3', subtype => 'A', rhap => 'A', snp3 => 'B', hap => '4' },
				AB => { type => '3', subtype => 'B', rhap => 'B', snp3 => 'B', hap => '3' },
				BB => { type => 'W1', subtype => 'I', rhap => undef, snp3 => undef, hap => undef }
			},
			BB => {
				AA => { type => 'W1', subtype => 'A', rhap => undef, snp3 => undef, hap => undef },
				AB => { type => '1', subtype => 'A', rhap => undef, snp3 => undef, hap => undef },
				BB => { type => 'W1', subtype => 'B', rhap => undef, snp3 => undef, hap => undef }
			}
		},
		AB => {
			AA => {
				AA => { type => '4a', subtype => 'A', rhap => 'A', snp3 => undef, hap => undef },
				AB => { type => '4b', subtype => 'B', rhap => 'B', snp3 => undef, hap => undef },
				BB => { type => 'W1', subtype => 'L', rhap => undef, snp3 => undef, hap => undef }
			},
			AB => {
				AA => { type => 'W2', subtype => 'A', rhap => undef, snp3 => undef, hap => undef },
				AB => { type => 'W2', subtype => 'B', rhap => undef, snp3 => undef, hap => undef },
				BB => { type => 'W2', subtype => 'C', rhap => undef, snp3 => undef, hap => undef }
			},
			BB => {
				AA => { type => 'W1', subtype => 'K', rhap => undef, snp3 => undef, hap => undef },
				AB => { type => '4b', subtype => 'A', rhap => 'A', snp3 => undef, hap => undef },
				BB => { type => '4a', subtype => 'B', rhap => 'B', snp3 => undef, hap => undef }
			}
		},
		BB => {
			AA => {
				AA => { type => 'W1', subtype => 'C', rhap => undef, snp3 => undef, hap => undef },
				AB => { type => '1', subtype => 'B', rhap => undef, snp3 => undef, hap => undef },
				BB => { type => 'W1', subtype => 'D', rhap => undef, snp3 => undef, hap => undef }
			},
			AB => {
				AA => { type => 'W1', subtype => 'J', rhap => undef, snp3 => undef, hap => undef },
				AB => { type => '3', subtype => 'C', rhap => 'A', snp3 => 'A', hap => '3' },
				BB => { type => '3', subtype => 'D', rhap => 'B', snp3 => 'A', hap => '4' }
			},
			BB => {
				AA => { type => 'W1', subtype => 'G', rhap => undef, snp3 => undef, hap => undef },
				AB => { type => 'W1', subtype => 'H', rhap => undef, snp3 => undef, hap => undef },
				BB => { type => '2', subtype => 'B', rhap => undef, snp3 => undef, hap => undef }
			}
		}
	},
	UC => {
		AA => {
			AA => {
				AA => { type => '2', subtype => 'A', rhap => undef, snp3 => undef, hap => undef },
				AB => { type => 'W1', subtype => 'E', rhap => undef, snp3 => undef, hap => undef },
				BB => { type => 'W1', subtype => 'F', rhap => undef, snp3 => undef, hap => undef }
			},
			AB => {
				AA => { type => '3', subtype => 'F', rhap => 'B', snp3 => 'B', hap => '3' },
				AB => { type => '3', subtype => 'E', rhap => 'A', snp3 => 'B', hap => '4' },
				BB => { type => 'W1', subtype => 'I', rhap => undef, snp3 => undef, hap => undef }
			},
			BB => {
				AA => { type => 'W1', subtype => 'A', rhap => undef, snp3 => undef, hap => undef },
				AB => { type => '1', subtype => 'A', rhap => undef, snp3 => undef, hap => undef },
				BB => { type => 'W1', subtype => 'B', rhap => undef, snp3 => undef, hap => undef }
			}
		},
		AB => {
			AA => {
				AA => { type => '4b', subtype => 'D', rhap => 'B', snp3 => undef, hap => undef },
				AB => { type => '4a', subtype => 'C', rhap => 'A', snp3 => undef, hap => undef },
				BB => { type => 'W1', subtype => 'L', rhap => undef, snp3 => undef, hap => undef }
			},
			AB => {
				AA => { type => 'W2', subtype => 'A', rhap => undef, snp3 => undef, hap => undef },
				AB => { type => 'W2', subtype => 'B', rhap => undef, snp3 => undef, hap => undef },
				BB => { type => 'W2', subtype => 'C', rhap => undef, snp3 => undef, hap => undef }
			},
			BB => {
				AA => { type => 'W1', subtype => 'K', rhap => undef, snp3 => undef, hap => undef },
				AB => { type => '4a', subtype => 'D', rhap => 'B', snp3 => undef, hap => undef },
				BB => { type => '4b', subtype => 'C', rhap => 'A', snp3 => undef, hap => undef }
			}
		},
		BB => {
			AA => {
				AA => { type => 'W1', subtype => 'C', rhap => undef, snp3 => undef, hap => undef },
				AB => { type => '1', subtype => 'B', rhap => undef, snp3 => undef, hap => undef },
				BB => { type => 'W1', subtype => 'D', rhap => undef, snp3 => undef, hap => undef }
			},
			AB => {
				AA => { type => 'W1', subtype => 'J', rhap => undef, snp3 => undef, hap => undef },
				AB => { type => '3', subtype => 'H', rhap => 'B', snp3 => 'A', hap => '4' },
				BB => { type => '3', subtype => 'G', rhap => 'A', snp3 => 'A', hap => '3' }
			},
			BB => {
				AA => { type => 'W1', subtype => 'G', rhap => undef, snp3 => undef, hap => undef },
				AB => { type => 'W1', subtype => 'H', rhap => undef, snp3 => undef, hap => undef },
				BB => { type => '2', subtype => 'B', rhap => undef, snp3 => undef, hap => undef }
			}
		}
	},
	AR => {
		AA => {
			AA => {
				AA => { type => '2', subtype => 'A', rhap => undef, snp3 => undef, hap => undef },
				AB => { type => '2', subtype => 'C', rhap => undef, snp3 => undef, hap => undef },
				BB => { type => '2', subtype => 'D', rhap => undef, snp3 => undef, hap => undef }
			},
			AB => {
				AA => { type => '3', subtype => 'A', rhap => 'A', snp3 => 'B', hap => '4' },
				AB => { type => 'W2', subtype => 'D', rhap => undef, snp3 => undef, hap => undef },
				BB => { type => '3', subtype => 'J', rhap => 'B', snp3 => 'B', hap => '3' }
			},
			BB => {
				AA => { type => '1', subtype => 'C', rhap => undef, snp3 => undef, hap => undef },
				AB => { type => '1', subtype => 'A', rhap => undef, snp3 => undef, hap => undef },
				BB => { type => '1', subtype => 'E', rhap => undef, snp3 => undef, hap => undef }
			}
		},
		AB => {
			AA => {
				AA => { type => '4a', subtype => 'A', rhap => 'A', snp3 => undef, hap => undef },
				AB => { type => 'W2', subtype => 'F', rhap => undef, snp3 => undef, hap => undef },
				BB => { type => '4b', subtype => 'F', rhap => 'B', snp3 => undef, hap => undef }
			},
			AB => {
				AA => { type => 'W2', subtype => 'A', rhap => undef, snp3 => undef, hap => undef },
				AB => { type => 'W2', subtype => 'B', rhap => undef, snp3 => undef, hap => undef },
				BB => { type => 'W2', subtype => 'C', rhap => undef, snp3 => undef, hap => undef }
			},
			BB => {
				AA => { type => '4b', subtype => 'E', rhap => 'A', snp3 => undef, hap => undef },
				AB => { type => 'W2', subtype => 'G', rhap => undef, snp3 => undef, hap => undef },
				BB => { type => '4a', subtype => 'B', rhap => 'B', snp3 => undef, hap => undef }
			}
		},
		BB => {
			AA => {
				AA => { type => '1', subtype => 'D', rhap => undef, snp3 => undef, hap => undef },
				AB => { type => '1', subtype => 'B', rhap => undef, snp3 => undef, hap => undef },
				BB => { type => '1', subtype => 'F', rhap => undef, snp3 => undef, hap => undef }
			},
			AB => {
				AA => { type => '3', subtype => 'I', rhap => 'A', snp3 => 'A', hap => '3' },
				AB => { type => 'W2', subtype => 'E', rhap => undef, snp3 => undef, hap => undef },
				BB => { type => '3', subtype => 'D', rhap => 'B', snp3 => 'A', hap => '4' }
			},
			BB => {
				AA => { type => '2', subtype => 'E', rhap => undef, snp3 => undef, hap => undef },
				AB => { type => '2', subtype => 'F', rhap => undef, snp3 => undef, hap => undef },
				BB => { type => '2', subtype => 'B', rhap => undef, snp3 => undef, hap => undef }
			}
		}
	},
	UR => {
		AA => {
			AA => {
				AA => { type => '2', subtype => 'A', rhap => undef, snp3 => undef, hap => undef },
				AB => { type => '2', subtype => 'C', rhap => undef, snp3 => undef, hap => undef },
				BB => { type => '2', subtype => 'D', rhap => undef, snp3 => undef, hap => undef }
			},
			AB => {
				AA => { type => '3', subtype => 'E', rhap => 'B', snp3 => 'B', hap => '3' },
				AB => { type => 'W2', subtype => 'D', rhap => undef, snp3 => undef, hap => undef },
				BB => { type => '3', subtype => 'L', rhap => 'A', snp3 => 'B', hap => '4' }
			},
			BB => {
				AA => { type => '1', subtype => 'C', rhap => undef, snp3 => undef, hap => undef },
				AB => { type => '1', subtype => 'A', rhap => undef, snp3 => undef, hap => undef },
				BB => { type => '1', subtype => 'E', rhap => undef, snp3 => undef, hap => undef }
			}
		},
		AB => {
			AA => {
				AA => { type => '4b', subtype => 'D', rhap => 'B', snp3 => undef, hap => undef },
				AB => { type => 'W2', subtype => 'F', rhap => undef, snp3 => undef, hap => undef },
				BB => { type => '4a', subtype => 'F', rhap => 'A', snp3 => undef, hap => undef }
			},
			AB => {
				AA => { type => 'W2', subtype => 'A', rhap => undef, snp3 => undef, hap => undef },
				AB => { type => 'W2', subtype => 'B', rhap => undef, snp3 => undef, hap => undef },
				BB => { type => 'W2', subtype => 'C', rhap => undef, snp3 => undef, hap => undef }
			},
			BB => {
				AA => { type => '4a', subtype => 'E', rhap => 'B', snp3 => undef, hap => undef },
				AB => { type => 'W2', subtype => 'G', rhap => undef, snp3 => undef, hap => undef },
				BB => { type => '4b', subtype => 'C', rhap => 'A', snp3 => undef, hap => undef }
			}
		},
		BB => {
			AA => {
				AA => { type => '1', subtype => 'D', rhap => undef, snp3 => undef, hap => undef },
				AB => { type => '1', subtype => 'B', rhap => undef, snp3 => undef, hap => undef },
				BB => { type => '1', subtype => 'F', rhap => undef, snp3 => undef, hap => undef }
			},
			AB => {
				AA => { type => '3', subtype => 'K', rhap => 'B', snp3 => 'A', hap => '4' },
				AB => { type => 'W2', subtype => 'E', rhap => undef, snp3 => undef, hap => undef },
				BB => { type => '3', subtype => 'H', rhap => 'A', snp3 => 'A', hap => '3' }
			},
			BB => {
				AA => { type => '2', subtype => 'E', rhap => undef, snp3 => undef, hap => undef },
				AB => { type => '2', subtype => 'F', rhap => undef, snp3 => undef, hap => undef },
				BB => { type => '2', subtype => 'B', rhap => undef, snp3 => undef, hap => undef }
			}
		}
	}
);

open OUT, ">$opts{out}_snpmap.tsv" or die "Can't open file: $opts{out}_snpmap.tsv\n";
print OUT join("\t", 'SNP.TYPE', 'SNP.SUBTYPE', 'Mo.GT', 'Fa.GT', 'CI.GT'), "\n";
foreach $a (sort keys %{$snpmap{$opts{proband}}}){
	foreach $b (sort keys %{$snpmap{$opts{proband}}->{$a}}){
		foreach (sort keys %{$snpmap{$opts{proband}}->{$a}->{$b}}){
			print OUT join("\t",
				$snpmap{$opts{proband}}->{$a}->{$b}->{$_}->{type},
				$snpmap{$opts{proband}}->{$a}->{$b}->{$_}->{subtype},
				$a, $b, $_
			), "\n";
		}
	}
}
close OUT;

my %DP = (DPNI => $opts{PDP});
$DP{$piles[$_]} = $opts{NDP} for 0 .. $#piles - 1;

foreach my $chr (@chr){

	open OUT, ">$opts{out}_chr${chr}.tsv" or die "Can't open file: $opts{out}_chr${chr}.tsv\n";
	print OUT join("\t", @header), "\n";

	print "Processing contig: $chr ...\n";

	undef %call;

	foreach my $sample (@piles){

		print "\tProcessing pile: $sample ...\n";

		open IN, "sed 1d $opts{$sample} | cut -f 2- | grep -e '^${chr}[[:space:]]' | cut -f 2- |" or die "Can't open file: $opts{$sample}\n";
		while(<IN>){

			@pile = split /\t/, $_;
			$pile{$_} = shift @pile foreach qw[ POS REF TRC A T C G N Ins Del ];

			next if $pile{TRC} < $DP{$sample} # Skip low covered sites
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
			next if $r != 1 || $#g != 1; # Skip remaining indels and ambiguous genotype

			@ab = grep(!/$pile{REF}/, @g);
			($pile{GT}, $pile{ALT}) = ($#ab < 0 ? 'AA' : grep(/1/, @n) ? 'AB' : 'BB', shift @ab);

			if(exists $call{$pile{POS}}->{ALT}){ # Used for properly defining ARC and BRC in DPNI[AA]
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

			$call{$pile{POS}}->{REF} = $pile{REF};
			$call{$pile{POS}}->{ALT}->{$pile{ALT}}++ if defined $pile{ALT};
			$call{$pile{POS}}->{$sample}->{$_} = $pile{$_} foreach qw[ GT TRC RRC ARC BRC ];
		}
		close IN;
	}

	foreach $k (sort { $a <=> $b } keys %call){
		foreach (@piles){
			delete $call{$k} and last unless exists $call{$k}->{$_};
		}
		next unless exists $call{$k};
		delete $call{$k} and next if $call{$k}->{Mo}->{GT} ne $call{$k}->{DPNI}->{GT};
		if(exists $call{$k}->{ALT}){
			@call = keys %{$call{$k}->{ALT}};
			delete $call{$k} and next if $#call != 0; # Skip multi-allelic sites
		}
		@call = ($chr, $k, $call{$k}->{REF}, (exists $call{$k}->{ALT}) ? ($call[0]) : ('NA'));
		push @call, $snpmap{$opts{proband}}->{$call{$k}->{Mo}->{GT}}->{$call{$k}->{Fa}->{GT}}->{$call{$k}->{CI}->{GT}}->{$_} foreach qw[ type subtype ];
		push @call, (defined $snpmap{$opts{proband}}->{$call{$k}->{Mo}->{GT}}->{$call{$k}->{Fa}->{GT}}->{$call{$k}->{CI}->{GT}}->{$_}) ?
			($snpmap{$opts{proband}}->{$call{$k}->{Mo}->{GT}}->{$call{$k}->{Fa}->{GT}}->{$call{$k}->{CI}->{GT}}->{$_}) : ('NA') foreach qw[ rhap snp3 hap ];
		push @call, $call{$k}->{$_}->{GT} foreach @piles;
		foreach $a (@piles){
			foreach $b (qw[ TRC RRC ARC BRC ]){
				push @call, $call{$k}->{$a}->{$b};
			}
		}
		print OUT join("\t", @call), "\n";
	}
	close OUT;
}
my $root_dir = dirname($opts{out});
`touch ${root_dir}/done`;