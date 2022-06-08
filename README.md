## Purpose

## Prerequisites
For a given sample the entry point is a R1/R2 pair of fastq.gz files.
A demultiplexing step may therefore be required, for example:
`bcl2fastq --ignore-missing-positions --ignore-missing-controls --ignore-missing-filter --ignore-missing-bcls -R "$input_dir" --sample-sheet "SampleSheet.csv" -o "$output_dir" --no-lane-splitting -r 2 -p 6 -w 2`
**NB**: do not forget the argument `--no-lane-splitting`

## Input data tree
Input data should be organized by capture design and sequencing run as follows:
`${capture_design_id}/${sequencing_run_id}/*_R[1/2]_001.fastq.gz`
So for a given capture design (parent directory) and a given sequencing run (sub-directory) we find the all set of R1/R2 pairs (fastq.gz) of each sample associated to a given family eg. maternal, paternal, proband nuclear dna and maternal plasmatic dna.

## Capture design
Each design must be associated with two bed files:
- A standard design.bed file describing the genomic intervals sequenced outside the ZFX/ZFY capture probes.
- A design_ZFXY.bed file restricted to all specific ZFX/ZFY probes.

Each file must contain 4 columns with the identifier of the targeted gene in column 4.
For off-gene probes, please use the mandatory keyword "Ethnie".
For ZFX/ZFY specific probes, please use the mandatory keywords "ZFX" or "ZFY".
All bed design files must be placed in the same parent directory.

## Input file describing a given analysis
It is possible to conduct the analysis of several families in parallel.
To do this you need to create a tabulated flat file describing all the entry points and parameters associated with a series of analyses. This file is very important because it describes all the inputs/outputs and parameters of the analysis. You should be very careful when editing this file.

### X-linked inheritance (rxli)
The 20 columns of the file are in this order:
- **design**: design identifier. Must be both the name of the bed files and the parent input directory containing the fastq.gz files arranged by sequencing run subdirectories (see Input data tree section).
- **run_id**: identifier of the sequencing run. Name of the subdirectory containing the fastq.gz files.
- **fam_id**: family identifier. Defines the name of the output directory of the analysis associated with the family. This identifier must be **unique** (see Output data tree section).
- **matmut**: known maternal mutation.
- **fastq**: dictionary of type key=value defining input fastq.gz files. Each key identifies the sample, for example **Mo** for mother, **Fa** for father, etc. Note that **Mo** and **Fa** are mandatory names; the naming of the other keys is free (for proband, plasmatic dna, etc.). Each value represents the base name of the file eg. to the value x corresponds a couple of files x_R1_001.fastq.gz and x_R2_001.fastq.gz. Finally the couples key=value are separated by a ";" eg. Mo=a;Fa=b; etc.
- **\[specific\] gender_fetus**: sex of the fetus ie. XX (female) or XY (male).
- **\[specific\] gender_proband**: sex of the proband ie. XX or XY.
- **proband**: type of proband ie. AC (affected child), UC (unaffected child), AR (affected close relative) or UR (unaffected close relative).
- **CI**: proband identifier (key defined in column **fastq**).
- **\[optional\] DPN**: fetal nuclear dna identifier (key defined in column **fastq**). If available and with the use of the **v2_rxli.Snakefile** only (see Pipeline versions section).
- **DPNI**: plasmatic maternal dna identifier (key defined in column **fastq**).
- **NDP**: minimum sequencing depth in nuclear dna samples to call a snp (default to 8X).
- **PDP**: minimum sequencing depth in plasmatic dna samples to call a snp (default to 15X).
- **HMZ**: threshold used to call a homozygous snp (default to 0.85 ie. \[0.85-1\]).
- **HTZ**: threshold used to call a heterozygous snp (default to 0.15 ie. \[0.35-0.65\]).
- **SNP**: minimum number of snp to conduct a sprt analysis (default to 10).
- **DIST**: minimum inter-snp distance within a given sprt block (default to 50).
- **ALPHA**: alpha risk of sprt ie. type I error.
- **BETA**: beta risk of sprt ie. type II error.
- **target_gene**: identifier of the pathological gene as described in the bed file.

### Autosomal inheritance
The 20 columns of the file are in this order:
- **design**, **run_id**, **fam_id**, **matmut**: same as rxli.
- **\[specific\] patmut**: known paternal mutation.
- **fastq**: same as rxli.
- **\[specific\] genmode**: mode of genetic inheritance ie. ADMI (autosomal dominant maternal inheritance), ADPI (autosomal dominant paternal inheritance) or ARI (autosomal recessive inheritance).
- **proband**, **CI**: same as rxli.
- **\[optional\] DPN**: same as rxli. If available and with the use of the **v2.Snakefile** only (see Pipeline versions section).
- **DPNI**, **NDP**, **PDP**, **HMZ**, **HTZ**, **SNP**, **DIST**, **ALPHA**, **BETA**, **target_gene**: same as rxli.

## Input files to launch an analysis
2 files are required to run a given analysis. These 2 files (one .tsv and one .json) are automatically generated by a case specific perl script ie. according to the mode of inheritance and the availability of fetal nuclear dna (DPN sample).

### X-linked inheritance (rxli)
In case the DPN sample is available:
```
perl batch2json_rxli-1.pl \
  --input-file $input-file \           # path to input batch descriptor (tsv)
  --output-dir $output_dir \           # path to output analysis directory
  --output-filename $output_filename \ # path to output filename (no extension)
  --read-length $read-length \         # minimal read length after trimming
  --site $site \                       # bwa mem read group PU tag
  --adapters $adapters \               # Illumina adapter sequences to trim (fasta)
  --genome $genome \                   # reference genome assembly (genomeFile.txt)
  --fasta $fasta \                     # reference genome sequences (fasta)
  --fastq $fastq \                     # directory containing R1/R2.fastq.gz files
  --indels $indels \                   # gold standard indels used for GATK (vcf)
  --snps $snps \                       # gold standard snps used for GATK (vcf)
  --bed $bed                           # directory containing bed design files
```
If the DPN sample is **not** available, use the `batch2json_rxli-2.pl` script instead.

### Autosomal inheritance
Same as above, use `batch2json-1.pl` if the DPN sample is available, `batch2json-2.pl` otherwise. 


