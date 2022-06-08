#!/usr/bin/perl
use strict;
use warnings FATAL => qw[ numeric uninitialized ];
use Getopt::Long qw[ GetOptions ];
use File::Basename;

sub usage {
	my $args = shift;
	print "Usage: perl $0\n";
	print "--input-file\n"; # Input batch descriptor (tab-delimited flat file)
	print "--output-dir\n"; # Output analysis directory
	print "--output-filename\n"; # Output file basename (w/o file extension)
	print "--read-length\n"; # Minimal read length after trimming
	print "--site\n"; # BWA-MEM read group PU tag
	print "--adapters\n"; # Illumina adapter sequences to trim (fasta)
	print "--genome\n"; # Human reference genome assembly (genomeFile.txt)
	print "--fasta\n"; # Human reference genome sequences (fasta)
	print "--fastq\n"; # Path to input directory containing R1/R2.fastq files
	print "--indels\n"; # Gold standard indels used for GATK bam post-processing (vcf)
	print "--snps\n"; # Gold standard snps used for GATK bam post-processing (vcf)
	print "--bed\n"; # Path to input directory containing capture bed design
	if(defined($args)){
		print "ERROR: $args\n";
		exit 1;
	}
	exit 0;
}
my %opts = ();
GetOptions(\%opts, 'input-file=s', 'output-dir=s', 'output-filename=s', 'read-length=i', 'site=s', 'adapters=s', 'genome=s', 'fasta=s', 'fastq=s', 'indels=s', 'snps=s', 'bed=s');

my($i, $j, $k, $fastq);
my $n = 1;
my @items = qw[ fastq gender_fetus gender_proband proband CI DPNI NDP PDP HMZ HTZ SNP DIST ALPHA BETA ];
my @gender = ();
my @all_params = ();
my @new_params = ();
my @all_fastq = ();
my @new_fastq = ();
my %fastq_ids = ();
my @readline = ();
my %readline = ();

chomp($k = `wc -l < $opts{'input-file'}`);
open IN, "$opts{'input-file'}" or die $!;

chomp(my $header = readline(IN));
my @header = split /\t/, $header;

$opts{'output-basename'} = basename($opts{'output-filename'});

open TSV, ">$opts{'output-filename'}.tsv" or die $!;
print TSV join("\t", qw[ design run_id fam_id matmut gender_fetus gender_proband proband Mo Fa CI DPNI NDP PDP HMZ HTZ SNP DIST ALPHA BETA target_gene analysis_id ]), "\n";

open JSON, ">$opts{'output-filename'}.json" or die $!;
print JSON qq[{
"description": "",
"author": "Romain Daveau",
"path": {
"OUT": "$opts{'output-dir'}"
},
"info": {
"RUN_DESCRIPTOR": "$opts{'output-basename'}",
"READ_LENGTH": "$opts{'read-length'}",
"PANEL_NAME": "NIPT_SPRT_RHDO",
"SEQUENCER": "ILLUMINA",
"TARGET_TYPE": "CAPTURE",
"HOSPIT": "$opts{site}",
"ADAPTERS": "$opts{adapters}",
"FASTA_FILE": "$opts{fasta}",
"GENOME_FILE": "$opts{genome}",
"KNOWN_INDELS": "$opts{indels}",
"KNOWN_SNPS": "$opts{snps}"
},
"fam_ids": {];

while(<IN>){
	$n++;
	chomp;
	@readline = split /\t/, $_;
	$readline{$header[$_]} = shift @readline for 0 .. $#header;

	print JSON qq!
"$readline{fam_id}": {
"design": "$readline{design}", "bed": "$opts{bed}/$readline{design}.bed",
"matmut": "$readline{matmut}", "target_gene": "$readline{target_gene}",
"sprt": [
!;

	foreach (@items){
		@{$readline{$_}} = split /;/, delete $readline{$_};
	}

	%fastq_ids = ();
	foreach $fastq (@{$readline{fastq}}){
		@_ = split /=/, $fastq;
		$fastq_ids{$_[0]} = $_[1];
	}

	@all_fastq = ();
	@all_params = ();

	for $i (0 .. $#{$readline{proband}}){

		@readline = ();
		push @readline, $readline{$_} foreach qw[ design run_id fam_id matmut ];
		push @readline, ${$readline{$_}}[$i] foreach qw[ gender_fetus gender_proband proband ];
		push @readline, $fastq_ids{$_} foreach qw[ Mo Fa ];
		push @readline, $fastq_ids{${$readline{$_}}[$i]} foreach qw[ CI DPNI ];
		push @readline, ${$readline{$_}}[$i] foreach qw[ NDP PDP HMZ HTZ SNP DIST ALPHA BETA ];
		print TSV join("\t", @readline, $readline{target_gene}, $i), "\n";

		@gender = ();
		push @gender, ${$readline{$_}}[$i] =~ 'Y' ? 'M' : 'F' foreach qw[ gender_fetus gender_proband ];
		@new_params = join(": ", "\"gender\"", "\"" . join('', @gender) . "\"");

		for $j (3 .. $#items){
			push @new_params, join(": ", "\"$items[$j]\"", "\"${$readline{$items[$j]}}[$i]\"");
		}
		push @all_params, join(" ", "{", join(", ", @new_params), "}");
	}
	print JSON join(",\n", @all_params);
	print JSON "\n],\n\"samples\": {\n";

	foreach $fastq (keys %fastq_ids){
		@new_fastq = ();
		for $i (1 .. 2){
			push @new_fastq, join(": ", "\"R$i\"", "\"$opts{fastq}/$readline{design}/$readline{run_id}/$fastq_ids{$fastq}_R${i}_001.fastq.gz\"");
		}
		push @all_fastq, join("\n", "\"$fastq\": {", join(",\n", @new_fastq), "}");
	}

	print JSON join(",\n", @all_fastq);
	print JSON "\n}";
	print JSON $n == $k ? "\n}" : "\n},";
}
close IN;
close TSV;

print JSON qq!
},
"trimfaq": { "OPTIONS": "threads=4 hdist=1 qtrim=rl trimq=20 tpe tbo" },
"bwa_mem": { "OPTIONS": "-t 8 -M" },
"dedup": { "OPTIONS": "VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true ASSUME_SORT_ORDER=coordinate" },
"gatkRTC": { "OPTIONS": "-nt 4" },
"gatkBQSR": { "OPTIONS": "-nct 4" },
"gatkPR": { "OPTIONS": "-nct 4" },
"pileup": { "OPTIONS": "-A -B -Q 0 -d 1000000" }
}
!;

close JSON;