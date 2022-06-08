workdir: config["path"]["OUT"]
output_path = config["path"]["OUT"]
mate_ids = ["R1", "R2"]
sample_ids = dict()
fam_ids = []
sprt_out = []

for fam_id in config["fam_ids"].keys():
	fam_ids.append(fam_id)
	sample_ids[fam_id] = [sample for sample in config["fam_ids"][fam_id]["samples"].keys()]
	for (RANK, sprt) in enumerate(config["fam_ids"][fam_id]["sprt"]):
		sprt_out.append(output_path + "/{}/{}/sprt/{RANK}/NDP{NDP}_PDP{PDP}_HMZ{HMZ}_HTZ{HTZ}_SNP{SNP}_DIST{DIST}_ALPHA{ALPHA}_BETA{BETA}/r-out/done".format(config["fam_ids"][fam_id]["design"], fam_id, RANK=RANK, **sprt))

pile2tab = []
for fam_id in fam_ids:
	for sample_id in sample_ids[fam_id]:
		pile2tab.append(output_path + "/{}/{}/samples/{}/{}_final.pile".format(config["fam_ids"][fam_id]["design"], fam_id, sample_id, sample_id))

rule all:
	input:
		sprt = sprt_out,
		tab = output_path + "/" + config["info"]["RUN_DESCRIPTOR"] + ".tsv"
	output:
		output_path + "/" + config["info"]["RUN_DESCRIPTOR"] + "_sprt.tsv"
	shell:
		"sprt_summary_rxli-3.pl \
		--input-dir " + output_path + " \
		--input-file {input.tab} \
		--output-file {output}"

rule getfaq:
	input:
		R1 = lambda wildcards: config["fam_ids"][wildcards.fam_id]["samples"][wildcards.sample]["R1"],
		R2 = lambda wildcards: config["fam_ids"][wildcards.fam_id]["samples"][wildcards.sample]["R2"]
	output:
		R1 = "{prefix}/{design}/{fam_id}/samples/{sample}/{sample}.R1.fastq.gz",
		R2 = "{prefix}/{design}/{fam_id}/samples/{sample}/{sample}.R2.fastq.gz"
	shell:
		'ln -s {input.R1} {output.R1} \
		;ln -s {input.R2} {output.R2} \
		2>&1 | tee -a {log}'

rule trimfaq:
	input:
		R1 = "{prefix}/{design}/{fam_id}/samples/{sample}/{sample}.R1.fastq.gz",
		R2 = "{prefix}/{design}/{fam_id}/samples/{sample}/{sample}.R2.fastq.gz"
	output:
		R1 = "{prefix}/{design}/{fam_id}/samples/{sample}/{sample}.R1_cleaned.fastq.gz",
		R2 = "{prefix}/{design}/{fam_id}/samples/{sample}/{sample}.R2_cleaned.fastq.gz"
	log:
		"{prefix}/{design}/{fam_id}/samples/{sample}/logs/{sample}_trimfaq.log"
	params:
		trimfaq_options = config["trimfaq"]["OPTIONS"],
		reference = config["info"]["ADAPTERS"],
		read_length = config["info"]["READ_LENGTH"]
	shell:
		'bbduk.sh -Xmx8g -Xms8g \
		in1={input.R1} in2={input.R2} \
		out1={output.R1} out2={output.R2} \
		ref={params.reference} minlen={params.read_length} {params.trimfaq_options} \
		2>&1 | tee -a {log}'

rule bwa_mem:
	input:
		R1 = "{prefix}/{design}/{fam_id}/samples/{sample}/{sample}.R1_cleaned.fastq.gz",
		R2 = "{prefix}/{design}/{fam_id}/samples/{sample}/{sample}.R2_cleaned.fastq.gz"
	output:
		temp("{prefix}/{design}/{fam_id}/samples/{sample}/{sample}_bwa_mem.sam")
	log:
		"{prefix}/{design}/{fam_id}/samples/{sample}/logs/{sample}_bwa_mem.log"
	params:
		bwa_mem_options = config["bwa_mem"]["OPTIONS"],
		reference = config["info"]["FASTA_FILE"],
		RG_ID = config["info"]["PANEL_NAME"],
		RG_PL = config["info"]["SEQUENCER"],
		RG_LB = config["info"]["TARGET_TYPE"],
		RG_PU = config["info"]["HOSPIT"]
	shell:
		'docker_wrapper -c \"bwa mem {params.bwa_mem_options} \
		-R @RG\\tID:{params.RG_ID}\\tPL:{params.RG_PL}\\tLB:{params.RG_LB}\\tPU:{params.RG_PU}\\tSM:{wildcards.sample} \
		{params.reference} {input.R1} {input.R2}\" -o {output} \
		2>&1 | tee -a {log}'

rule sam2bam:
	input:
		"{filename}.sam"
	output:
		temp("{filename}.bam")
	shell:
		'docker_wrapper -c \"samtools view -hbS -o {output} {input} \" 2>&1'

rule bamsort:
	input:
		"{filename}.bam"
	output:
		temp("{filename}_sorted.bam")
	shell:
		'docker_wrapper -c \"samtools sort -O bam -o {output} {input} \" 2>&1'

rule bam2bai:
	input:
		"{filename}.bam"
	output:
		temp("{filename}.bam.bai")
	shell:
		'docker_wrapper -c \"samtools index {input} {output} \" 2>&1'

rule dedup:
	input:
		"{prefix}/{design}/{fam_id}/samples/{sample}/{sample}_{ext}.bam"
	output:
		BAM = temp("{prefix}/{design}/{fam_id}/samples/{sample}/{sample}_{ext}_dedup.bam"),
		TXT = temp("{prefix}/{design}/{fam_id}/samples/{sample}/{sample}_{ext}_dedup.txt")
	log:
		"{prefix}/{design}/{fam_id}/samples/{sample}/logs/{sample}_dedup.log"
	params:
		dedup_options = config["dedup"]["OPTIONS"]
	shell:
		'java -Xmx7500m -Xms7500m -XX:ParallelGCThreads=1 -XX:+AggressiveHeap -jar /usr/share/java/picard.jar MarkDuplicates \
		{params.dedup_options} INPUT={input} OUTPUT={output.BAM} METRICS_FILE={output.TXT} \
		2>&1 | tee -a {log}'

rule gatkRTC:
	input:
		BAM = "{prefix}/{design}/{fam_id}/samples/{sample}/{sample}_{ext}.bam",
		BAI = "{prefix}/{design}/{fam_id}/samples/{sample}/{sample}_{ext}.bam.bai"
	output:
		temp("{prefix}/{design}/{fam_id}/samples/{sample}/{sample}_{ext}.intervals")
	log:
		"{prefix}/{design}/{fam_id}/samples/{sample}/logs/{sample}_gatkRTC.log"
	params:
		gatkRTC_options = config["gatkRTC"]["OPTIONS"],
		reference = config["info"]["FASTA_FILE"],
		refindels = config["info"]["KNOWN_INDELS"],
		intervals = lambda wildcards: config["fam_ids"][wildcards.fam_id]["bed"]
	shell:
		'java -Xmx31000m -Xms31000m -XX:ParallelGCThreads=4 -XX:+AggressiveHeap -jar /usr/share/java/GenomeAnalysisTK.jar -T RealignerTargetCreator \
		{params.gatkRTC_options} -R {params.reference} -L {params.intervals} --known {params.refindels} -I {input.BAM} -o {output} \
		2>&1 | tee -a {log}'

rule gatkIR:
	input:
		BAM = "{prefix}/{design}/{fam_id}/samples/{sample}/{sample}_{ext}.bam",
		BAI = "{prefix}/{design}/{fam_id}/samples/{sample}/{sample}_{ext}.bam.bai",
		INT = "{prefix}/{design}/{fam_id}/samples/{sample}/{sample}_{ext}.intervals"
	output:
		BAM = temp("{prefix}/{design}/{fam_id}/samples/{sample}/{sample}_{ext}_remap.bam"),
		BAI = temp("{prefix}/{design}/{fam_id}/samples/{sample}/{sample}_{ext}_remap.bai")
	log:
		"{prefix}/{design}/{fam_id}/samples/{sample}/logs/{sample}_gatkIR.log"
	params:
		reference = config["info"]["FASTA_FILE"],
		refindels = config["info"]["KNOWN_INDELS"],
		intervals = lambda wildcards: config["fam_ids"][wildcards.fam_id]["bed"]
	shell:
		'java -Xmx11000m -Xms11000m -XX:ParallelGCThreads=1 -XX:+AggressiveHeap -jar /usr/share/java/GenomeAnalysisTK.jar -T IndelRealigner \
		-R {params.reference} -L {params.intervals} -known {params.refindels} -targetIntervals {input.INT} -I {input.BAM} -o {output.BAM} \
		2>&1 | tee -a {log}'

rule gatkBQSR:
	input:
		BAM = "{prefix}/{design}/{fam_id}/samples/{sample}/{sample}_{ext}.bam",
		BAI = "{prefix}/{design}/{fam_id}/samples/{sample}/{sample}_{ext}.bam.bai"
	output:
		"{prefix}/{design}/{fam_id}/samples/{sample}/{sample}_{ext}_qtable.txt"
	log:
		"{prefix}/{design}/{fam_id}/samples/{sample}/logs/{sample}_gatkBQSR.log"
	params:
		gatkBQSR_options = config["gatkBQSR"]["OPTIONS"],
		reference = config["info"]["FASTA_FILE"],
		intervals = lambda wildcards: config["fam_ids"][wildcards.fam_id]["bed"],
		refindels = config["info"]["KNOWN_INDELS"],
		refsnps = config["info"]["KNOWN_SNPS"]
	shell:
		'java -Xmx31000m -Xms31000m -XX:ParallelGCThreads=4 -XX:+AggressiveHeap -jar /usr/share/java/GenomeAnalysisTK.jar -T BaseRecalibrator \
		{params.gatkBQSR_options} -R {params.reference} -L {params.intervals} -knownSites {params.refsnps} -knownSites {params.refindels} -I {input.BAM} -o {output} \
		2>&1 | tee -a {log}'

rule gatkPR:
	input:
		BAM = "{prefix}/{design}/{fam_id}/samples/{sample}/{sample}_bwa_mem_sorted_dedup_remap.bam",
		BAI = "{prefix}/{design}/{fam_id}/samples/{sample}/{sample}_bwa_mem_sorted_dedup_remap.bam.bai",
		TXT = "{prefix}/{design}/{fam_id}/samples/{sample}/{sample}_bwa_mem_sorted_dedup_remap_qtable.txt"
	output:
		BAM = "{prefix}/{design}/{fam_id}/samples/{sample}/{sample}_final.bam",
		BAI = "{prefix}/{design}/{fam_id}/samples/{sample}/{sample}_final.bai"
	log:
		"{prefix}/{design}/{fam_id}/samples/{sample}/logs/{sample}_gatkPR.log"
	params:
		gatkPR_options = config["gatkPR"]["OPTIONS"],
		reference = config["info"]["FASTA_FILE"],
		intervals = lambda wildcards: config["fam_ids"][wildcards.fam_id]["bed"]
	shell:
		'java -Xmx11000m -Xms11000m -XX:ParallelGCThreads=4 -XX:+AggressiveHeap -jar /usr/share/java/GenomeAnalysisTK.jar -T PrintReads \
		{params.gatkPR_options} -R {params.reference} -L {params.intervals} -BQSR {input.TXT} -I {input.BAM} -o {output.BAM} \
		2>&1 | tee -a {log}'

rule pileup:
	input:
		BAM = "{prefix}/{design}/{fam_id}/samples/{sample}/{sample}_final.bam",
		BAI = "{prefix}/{design}/{fam_id}/samples/{sample}/{sample}_final.bai"
	output:
		temp("{prefix}/{design}/{fam_id}/samples/{sample}/{sample}_final.pileup")
	log:
		"{prefix}/{design}/{fam_id}/samples/{sample}/logs/{sample}_pileup.log"
	params:
		pileup_options = config["pileup"]["OPTIONS"],
		reference = config["info"]["FASTA_FILE"],
		bedtarget = lambda wildcards: config["fam_ids"][wildcards.fam_id]["bed"]
	shell:
		'docker_wrapper -c \'samtools mpileup {params.pileup_options} \
		-f {params.reference} -l {params.bedtarget} {input.BAM} -o {output}\' \
		2>&1 | tee -a {log}'

rule ZFXY:
	input:
		"{prefix}/{design}/{fam_id}/samples/{sample}/{sample}_final.bam"
	output:
		"{prefix}/{design}/{fam_id}/samples/{sample}/{sample}_ZFXY.tsv"
	log:
		"{prefix}/{design}/{fam_id}/samples/{sample}/logs/{sample}_ZFXY.log"
	params:
		bedtarget = lambda wildcards: str(config["fam_ids"][wildcards.fam_id]["bed"]).replace(".bed", "_ZFXY.bed")
	shell:
		'bedtools coverage -a {params.bedtarget} -b {input} -d 1> {output} 2>&1 | tee -a {log}'

rule pile2tab:
	input:
		"{prefix}/{design}/{fam_id}/samples/{sample}/{sample}_final.pileup"
	output:
		"{prefix}/{design}/{fam_id}/samples/{sample}/{sample}_final.pile"
	log:
		"{prefix}/{design}/{fam_id}/samples/{sample}/logs/{sample}_pile2tab.log"
	shell:
		'check_variants.py <(echo -e \"{wildcards.sample},{input}\")> {output} \
		2>&1 | tee -a {log}'

rule sprt_in:
	input:
		Mo = lambda wildcards: "{}/{}/{}/samples/Mo/Mo_final.pile".format(wildcards.prefix, wildcards.design, wildcards.fam_id),
		Fa = lambda wildcards: "{}/{}/{}/samples/Fa/Fa_final.pile".format(wildcards.prefix, wildcards.design, wildcards.fam_id),
		CI = lambda wildcards: "{}/{}/{}/samples/{}/{}_final.pile".format(wildcards.prefix, wildcards.design, wildcards.fam_id, config["fam_ids"][wildcards.fam_id]["sprt"][int(wildcards.RANK)]["CI"], config["fam_ids"][wildcards.fam_id]["sprt"][int(wildcards.RANK)]["CI"]),
		DPNI = lambda wildcards: "{}/{}/{}/samples/{}/{}_final.pile".format(wildcards.prefix, wildcards.design, wildcards.fam_id, config["fam_ids"][wildcards.fam_id]["sprt"][int(wildcards.RANK)]["DPNI"], config["fam_ids"][wildcards.fam_id]["sprt"][int(wildcards.RANK)]["DPNI"])
	params:
		gender = lambda wildcards: config["fam_ids"][wildcards.fam_id]["sprt"][int(wildcards.RANK)]["gender"],
		proband = lambda wildcards: config["fam_ids"][wildcards.fam_id]["sprt"][int(wildcards.RANK)]["proband"]
	output:
		"{prefix}/{design}/{fam_id}/sprt/{RANK}/NDP{NDP}_PDP{PDP}_HMZ{HMZ}_HTZ{HTZ}_SNP{SNP}_DIST{DIST}_ALPHA{ALPHA}_BETA{BETA}/perl-out/done"
	log:
		"{prefix}/{design}/{fam_id}/sprt/{RANK}/NDP{NDP}_PDP{PDP}_HMZ{HMZ}_HTZ{HTZ}_SNP{SNP}_DIST{DIST}_ALPHA{ALPHA}_BETA{BETA}/logs/{fam_id}_sprt_in.log"
	shell:
		'snptyping_rxli-2.pl \
		--Mo {input.Mo} --Fa {input.Fa} --CI {input.CI} --DPNI {input.DPNI} --proband {params.proband} \
		--gender {params.gender} --NPD {wildcards.NDP} --PDP {wildcards.PDP} --HMZ {wildcards.HMZ} --HTZ {wildcards.HTZ} \
		--out {wildcards.prefix}/{wildcards.design}/{wildcards.fam_id}/sprt/{wildcards.RANK}/NDP{wildcards.NDP}_PDP{wildcards.PDP}_HMZ{wildcards.HMZ}_HTZ{wildcards.HTZ}_SNP{wildcards.SNP}_DIST{wildcards.DIST}_ALPHA{wildcards.ALPHA}_BETA{wildcards.BETA}/perl-out/{wildcards.fam_id} \
		2>&1 | tee -a {log}'

rule sprt_out:
	input:
		DPNI = lambda wildcards: "{}/{}/{}/samples/{}/{}_ZFXY.tsv".format(wildcards.prefix, wildcards.design, wildcards.fam_id, config["fam_ids"][wildcards.fam_id]["sprt"][int(wildcards.RANK)]["DPNI"], config["fam_ids"][wildcards.fam_id]["sprt"][int(wildcards.RANK)]["DPNI"]),
		SPRT = "{prefix}/{design}/{fam_id}/sprt/{RANK}/NDP{NDP}_PDP{PDP}_HMZ{HMZ}_HTZ{HTZ}_SNP{SNP}_DIST{DIST}_ALPHA{ALPHA}_BETA{BETA}/perl-out/done"
	params:
		matmut = lambda wildcards: config["fam_ids"][wildcards.fam_id]["matmut"],
		bedtarget = lambda wildcards: config["fam_ids"][wildcards.fam_id]["bed"],
		gender = lambda wildcards: config["fam_ids"][wildcards.fam_id]["sprt"][int(wildcards.RANK)]["gender"],
		proband = lambda wildcards: config["fam_ids"][wildcards.fam_id]["sprt"][int(wildcards.RANK)]["proband"]
	output:
		"{prefix}/{design}/{fam_id}/sprt/{RANK}/NDP{NDP}_PDP{PDP}_HMZ{HMZ}_HTZ{HTZ}_SNP{SNP}_DIST{DIST}_ALPHA{ALPHA}_BETA{BETA}/r-out/done"
	log:
		"{prefix}/{design}/{fam_id}/sprt/{RANK}/NDP{NDP}_PDP{PDP}_HMZ{HMZ}_HTZ{HTZ}_SNP{SNP}_DIST{DIST}_ALPHA{ALPHA}_BETA{BETA}/logs/{fam_id}_sprt_out.log"
	shell:
		'snptyping_rxli-3.R \
		--ZFXY={input.DPNI} \
		--pref={wildcards.fam_id} \
		--idir={wildcards.prefix}/{wildcards.design}/{wildcards.fam_id}/sprt/{wildcards.RANK}/NDP{wildcards.NDP}_PDP{wildcards.PDP}_HMZ{wildcards.HMZ}_HTZ{wildcards.HTZ}_SNP{wildcards.SNP}_DIST{wildcards.DIST}_ALPHA{wildcards.ALPHA}_BETA{wildcards.BETA}/perl-out \
		--odir={wildcards.prefix}/{wildcards.design}/{wildcards.fam_id}/sprt/{wildcards.RANK}/NDP{wildcards.NDP}_PDP{wildcards.PDP}_HMZ{wildcards.HMZ}_HTZ{wildcards.HTZ}_SNP{wildcards.SNP}_DIST{wildcards.DIST}_ALPHA{wildcards.ALPHA}_BETA{wildcards.BETA}/r-out \
		--bed={params.bedtarget} --matmut={params.matmut} --proband={params.proband} \
		--gender={params.gender} --HTZ={wildcards.HTZ} \
		--min.NDP={wildcards.NDP} --min.PDP={wildcards.PDP} \
		--min.SNP={wildcards.SNP} --min.dist={wildcards.DIST} \
		--alpha={wildcards.ALPHA} --beta={wildcards.BETA} \
		2>&1 | tee -a {log}'