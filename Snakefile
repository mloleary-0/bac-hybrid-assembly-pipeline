shell.executable("/bin/bash")

# on our server, the best way to run this will to navigate to the directory with 'Snakefile', and run:
# activate conda asm-pipeline-env
# 'snakemake -j 14'
#
# where -j 14 specificies that many cores to use at once as a maximum for the pipeline.  Our server is 16 cores, so this is a reasonable number.  Most steps that benefit from multithreading are set to use the provided number of cores, and won't automatically grab available cores.

# "rule medaka" may need to be modified to update the basecalling model depending on what's available.  r941_min_high_g330 means R9.4.1 flow cell, min = minion, high = high accuracy model (hac), g330 = guppy 3.3.0 or newer.# Available: r103_fast_g507, r103_fast_snp_g507, r103_fast_variant_g507, r103_hac_g507, r103_hac_snp_g507, r103_hac_variant_g507, r103_min_high_g345, r103_min_high_g360, r103_prom_high_g360, r103_prom_snp_g3210, r103_prom_variant_g3210, r103_sup_g507, r103_sup_snp_g507, r103_sup_variant_g507, r10_min_high_g303, r10_min_high_g340, r941_min_fast_g303, r941_min_fast_g507, r941_min_fast_snp_g507, r941_min_fast_variant_g507, r941_min_hac_g507, r941_min_hac_snp_g507, r941_min_hac_variant_g507, r941_min_high_g303, r941_min_high_g330, r941_min_high_g340_rle, r941_min_high_g344, r941_min_high_g351, r941_min_high_g360, r941_min_sup_g507, r941_min_sup_snp_g507, r941_min_sup_variant_g507, r941_prom_fast_g303, r941_prom_fast_g507, r941_prom_fast_snp_g507, r941_prom_fast_variant_g507, r941_prom_hac_g507, r941_prom_hac_snp_g507, r941_prom_hac_variant_g507, r941_prom_high_g303, r941_prom_high_g330, r941_prom_high_g344, r941_prom_high_g360, r941_prom_high_g4011, r941_prom_snp_g303, r941_prom_snp_g322, r941_prom_snp_g360, r941_prom_sup_g507, r941_prom_sup_snp_g507, r941_prom_sup_variant_g507, r941_prom_variant_g303, r941_prom_variant_g322, r941_prom_variant_g360

# For example the model named r941_min_fast_g303 should be used with data from MinION (or GridION) R9.4.1 flowcells using the fast Guppy basecaller version 3.0.3. By contrast the model r941_prom_hac_g303 should be used with PromethION data and the high accuracy basecaller (termed "hac" in Guppy configuration files). Where a version of Guppy has been used without an exactly corresponding medaka model, the medaka model with the highest version equal to or less than the guppy version should be selected.

# If you need to add rules, FYI snakemake will automatically make directories when an output file goes into a directory that doesn't exist, but will not make directories when a directory itself is an output target (without a specific mkdir call).

# view dag output (workflow) in dot, in graphviz package (apt install graphviz): snakemake --dag | dot -Tpdf > dag.pdf

#mamba create -n medaka-test -c conda-forge -c bioconda medaka

# mamba create -n asm-pipeline-env -c conda-forge -c bioconda medaka fastqc multiqc flye circlator checkm-genome prokka filtlong

# be aware that the channel order seems to matter a lot?  conda forge must be before bioconda or medaka causes the install to fail. may be resolveable by 'conda config --set channel_priority flexible' or 'strict'


#---- Notes ----#
# The logic for this pipeline is to replace pilon with the RR Wick's Polypolish, and probably POLCA.  Instead of trying to be clever and auto-detect when polishing fails to make changes, I will run Polypolish a flat three times after Medaka (and probably POLCA).  Polypolish apparently will not introduce errors and is less likely to introduce errors in repeat regions.
#
#
# View the workflow with: snakemake --rulegraph | dot -Tpdf > dag.pdf
# note that this won't show many steps without at least dummy input files in the appropriate directories.
#
#---- Malleable variables ----#
# Set these with --config
# flyemethod = [nano-raw|nano-hq] (default: nano-raw)
# guppymodel = (default: r941_min_high_g330)
# readfmt = "fastq.gz"


#------------------------------#

#---- Fixed variables ----#
# changes here may need changes elsewhere

INPUT_DIR = "long_reads"
output_dir = "assembly"
flye_dir = os.path.join(output_dir, "02-flye_assembly")
IDS, = glob_wildcards(INPUT_DIR + "/{sample}.fastq.gz")
circlator_out_dir = os.path.join(output_dir, "04-rotated")

#--------------------------#



#---- Pipeline ----#

configfile: "config/pipeline_config.yaml"

rule all:
    input:
        "01-QC_inputs/fastqc/multiqc_report.html", # multiqc
        expand("01-QC_inputs/read_info_histograms/{sample}.txt", sample = IDS),
        "07-checkM/checkM.results.tab",
        expand("08-prokka_annotation/{sample}.polished/{sample}.gbk", sample = IDS),
        expand("09-assembly_QC/long_read_mapping/{sample}/{sample}.coverage.tab", sample = IDS),
        expand("09-assembly_QC/long_read_mapping/{sample}/{sample}.idxstats.tab", sample = IDS), 
        expand("09-assembly_QC/long_read_mapping/{sample}/{sample}.reads.sorted.bam", sample = IDS),
        expand("09-assembly_QC/long_read_mapping/{sample}/{sample}.unmapped.fastq.gz", sample = IDS),
        expand("09-assembly_QC/short_read_mapping/{sample}/{sample}.unmapped.fastq.gz", sample = IDS),
        expand("09-assembly_QC/short_read_mapping/{sample}/{sample}.idxstats.tab", sample = IDS),
        expand("09-assembly_QC/short_read_mapping/{sample}/{sample}.coverage.tab", sample = IDS)
        


#        expand("02-flye_assembly/{sample}/assembly.fasta", sample = IDS), # flye assembly
#        expand("01-QC_inputs/fastqc/{sample}_fastqc.html", sample = IDS), # fastqc

rule QC_long_reads:
    input:
        fastq = "long_reads/{sample}.fastq.gz",
    output:
        "01-QC_inputs/fastqc/{sample}_fastqc.html",
#        directory("01-QC_inputs/fastqc/")
    threads: 
        workflow.cores
    log: 
        "logs/QC_long_reads.{sample}.log"
    conda:
        "config/environment.yml"
    shell:
        "fastqc --threads {threads} --nano --outdir 01-QC_inputs/fastqc {input.fastq} 2> {log}"

rule QC_short_reads:
    input:
        fastq = "short_reads/{sample}.fastq.gz",
    output:
        "01-QC_inputs/fastqc/{sample}_fastqc.html",
#        directory("01-QC_inputs/fastqc/")
    log: 
        "logs/QC_short_reads.{sample}.log"
    threads: 
        workflow.cores
    conda:
        "config/environment.yml"
    shell:
        "fastqc --threads {threads} --outdir 01-QC_inputs/fastqc {input.fastq} 2> {log}"


rule QC_multiqc:
    input: expand(["01-QC_inputs/fastqc/{sample}_fastqc.html"], sample = IDS)
    output: "01-QC_inputs/fastqc/multiqc_report.html"
    log: 
        expand("logs/QC_multiqc.{sample}.log", sample = IDS)
    conda:
        "config/environment.yml"
    params:
        multiqc_err_log = "logs/multiqc.log"
    shell: "multiqc -o 01-QC_inputs/fastqc 01-QC_inputs/fastqc 2> {params.multiqc_err_log}"
    


rule read_info_hists:
    input:
        fastq = "long_reads/{sample}.fastq.gz"
    output:
        file = "01-QC_inputs/read_info_histograms/{sample}.txt"
    log: 
        "logs/read_info_hists.{sample}.log"
    conda:
        "config/environment.yml"
    shell:
        "bin/read_info_histograms.sh {input.fastq} 2> {log} 1> {output.file}" 

# Note: I modified this version of read_info_histograms.  It requires both histograms.py and filtlong, but assumes they are in basically the same directory or in a specific structure.  Unfortunately, read_info_histograms nor histograms.py come with the conda installation of filtlong as far as I can tell.  The modification I made was to point at a conda install of filtlong, though '$(which filtlong)'.  Also, read_info_histograms will take .gzipped or not files, and apparently will simultaniously look at illumina reads if you add them as inputs 2 and 3.
   
   
rule flye_assembly:  # Note: change '--nano-raw' to '-nano-hq' if using guppy 5+
    input:
        fastq = "long_reads/{sample}.fastq.gz"
    output:
        "02-flye_assembly/{sample}/assembly.fasta"
    threads: 
        workflow.cores
    log: 
        "logs/flye_assembly.{sample}.log"
    conda:
        "config/environment.yml"
    params:
        method = expand("{param}", param = config["flyemethod"])
    shell:
        "flye --{params.method} {input.fastq} --out-dir 02-flye_assembly/{wildcards.sample} --threads {threads} 2> {log}" 
# add --resume option for flye?

rule medaka:
    input: 
        assembly = "02-flye_assembly/{sample}/assembly.fasta",
        fastq = "long_reads/{sample}.fastq.gz"
    output: "03-medaka_polish/{sample}/consensus.fasta"
    threads:
        workflow.cores
    log: 
        "logs/medaka.{sample}.log"
    conda:
        "config/environment.yml"
    params:
        model = config["guppymodel"],
        outdir = "03-medaka_polish/{sample}"
    shell:
        "medaka_consensus -i {input.fastq} -d {input.assembly} -o {params.outdir} -t {threads} -m {params.model} 2> {log}"

# see this for pilon to repeat this 6 times or until a condition is met: https://stackoverflow.com/questions/59535549/execute-snakemake-rule-repeatedly-until-certain-conditions-are-met

# trying to make this one iterative for 3x, check this SO question.  https://stackoverflow.com/questions/56274065/snakemake-using-a-rule-in-a-loop
# however I don't know that I need to, according to the FAQ for polypolish it doesn't usually make changes in subsequent rounds.

rule pre_rotate_bwa:
    input:
        draft = rules.medaka.output,
        SR1 = "short_reads/{sample}.R1.fastq.gz",
        SR2 = "short_reads/{sample}.R2.fastq.gz"
    output:
        SR1align = "04-polypolish/{sample}/alignments_1.sam",
        SR2align = "04-polypolish/{sample}/alignments_2.sam"
    threads:
        workflow.cores
    log: 
        SR1 = "logs/pre_rotate_bwa.{sample}.R1.log",
        SR2 = "logs/pre_rotate_bwa.{sample}.R2.log"
    run:
        shell("bwa index {input.draft}"),
        shell("bwa mem -t {threads} -a {input.draft} {input.SR1} 2> {log.SR1} 1> {output.SR1align}"),
        shell("bwa mem -t {threads} -a {input.draft} {input.SR2} 2> {log.SR2} 1> {output.SR2align}")

        
rule pre_rotate_polypolish_filter:
    input: 
        SR1 = rules.pre_rotate_bwa.output.SR1align,
        SR2 = rules.pre_rotate_bwa.output.SR2align
    output: 
        SR1filt = "04-polypolish/{sample}/filtered_1.sam",
        SR2filt = "04-polypolish/{sample}/filtered_2.sam"
    threads:
        workflow.cores
    log: 
        "logs/pre_rotate_polypolish_filter.{sample}.log"
    conda:
        "config/environment.yml"
    params:
    shell:
        "polypolish_insert_filter.py --in1 {input.SR1} --in2 {input.SR2} --out1 {output.SR1filt} --out2 {output.SR2filt} 2> {log}"

rule pre_rotate_polypolish:
    input: 
        draft = rules.medaka.output,
        SR1 = rules.pre_rotate_polypolish_filter.output.SR1filt,
        SR2 = rules.pre_rotate_polypolish_filter.output.SR2filt
    output: "04-polypolish/{sample}/consensus.fasta"
    threads:
        workflow.cores
    log: 
        "logs/pre_rotate_polypolish.{sample}.log"
    conda:
        "config/environment.yml"
    params:
    shell:
        "polypolish {input.draft} {input.SR1} {input.SR2} 2> {log} 1> {output}"



rule circlator:
    input: 
#        "02-flye_assembly/{sample}/assembly.fasta"# placeholder until I get medaka running
#        "03-medaka_polish/{sample}/consensus.fasta"
#        "03-pilon_polish/{sample}/pilon_final_result/{sample}.pilon.polished.fasta"
        rules.pre_rotate_polypolish.output
    output: 
        file = "05-rotated/{sample}.assembly.fasta"
    log: 
        "logs/circlator.{sample}.log"
    params:
        prefix = "05-rotated/{sample}.assembly"
    conda:
        "config/environment.yml"
    shell: 
        "circlator fixstart {input} {params.prefix} 2> {log}"



rule post_rotate_bwa:
    input:
        draft = rules.circlator.output.file,
        SR1 = "short_reads/{sample}.R1.fastq.gz",
        SR2 = "short_reads/{sample}.R2.fastq.gz"
    output:
        SR1align = "06-rotated_polypolish/{sample}/alignments_1.sam",
        SR2align = "06-rotated_polypolish/{sample}/alignments_2.sam"
    threads:
        workflow.cores
    log: 
        SR1 = "logs/post_rotate_bwa.{sample}.R1.log",
        SR2 = "logs/post_rotate_bwa.{sample}.R1.log"
    run:
        shell("bwa index {input.draft}"),
        shell("bwa mem -t {threads} -a {input.draft} {input.SR1} 2> {log} 1> {output.SR1align}"),
        shell("bwa mem -t {threads} -a {input.draft} {input.SR2} 2> {log} 1> {output.SR2align}")

        
rule post_rotate_polypolish_filter:
    input: 
        SR1 = rules.post_rotate_bwa.output.SR1align,
        SR2 = rules.post_rotate_bwa.output.SR2align
    output: 
        SR1filt = "06-rotated_polypolish/{sample}/filtered_1.sam",
        SR2filt = "06-rotated_polypolish/{sample}/filtered_2.sam"
    threads:
        workflow.cores
    log: 
        "logs/post_rotate_polypolish_filter.{sample}.log"
    conda:
        "config/environment.yml"
    params:
    shell:
        "polypolish_insert_filter.py --in1 {input.SR1} --in2 {input.SR2} --out1 {output.SR1filt} --out2 {output.SR2filt} 2> {log}"

rule post_rotate_polypolish:
    input: 
        draft = rules.circlator.output.file,
        SR1 = rules.post_rotate_polypolish_filter.output.SR1filt,
        SR2 = rules.post_rotate_polypolish_filter.output.SR2filt
    output: "06-rotated_polypolish/{sample}/{sample}.assembly.consensus.fasta"
    threads:
        workflow.cores
    log: 
        "logs/post_rotate_polypolish.{sample}.log"
    conda:
        "config/environment.yml"
    params:
        consolidated_dir = "06-rotated_polypolish/consolidated"
    shell:
        "polypolish {input.draft} {input.SR1} {input.SR2} 2> {log} 1> {output}"


rule consolidate_genomes:
    input: rules.post_rotate_polypolish.output
    output: "07-checkM/genomes/{sample}.assembly.consensus.fasta"
    log: 
        "logs/consolidate_genomes.{sample}.log"
    run:
        shell("cp {input} {output} 2> {log}"),
        shell("sed -i 's/_polypolish_polypolish//g' {output} 2>> {log}")

rule checkM:
    input: 
        expand(["07-checkM/genomes/{sample}.assembly.consensus.fasta"], sample = IDS) # Critically, this is how I get snakemake to key off the results of circlator, while only running checkm once.  Had to hard-code the target directory but that's ok.
    output: # had to hard-code these values as well, doing the directory as an output 
#        dir = directory("05-checkM"),
        file = "07-checkM/checkM.results.tab"
    threads: 
        workflow.cores
    log:
        expand("logs/checkM.{sample}.log", sample = IDS)
    conda:
        "config/environment.yml"
    params:
        genomes = "07-checkM/genomes",
        checkm_output = "07-checkM",
        name = "checkM.results.tab",
        checkm_err_log = "logs/checkM.log"
    shell: 
        "checkm lineage_wf -t {threads} --file {params.checkm_output}/{params.name} --tab_table -x .fasta {params.genomes} {params.checkm_output} 2> {params.checkm_err_log}"


rule prokka:
    input:
        polished = rules.consolidate_genomes.output # Medaka output, placeholder until medaka works
    output:
        polished = "08-prokka_annotation/{sample}.polished/{sample}.gbk" # medaka
    params:
        polished_outdir = "08-prokka_annotation/{sample}.polished/"
    threads:
        workflow.cores
    log: 
        "logs/prokka.{sample}.log"
    conda:
        "config/environment.yml"
    shell:
        "prokka --outdir {params.polished_outdir} --force --cpus {threads} --locustag {wildcards.sample} --prefix {wildcards.sample} {input.polished} 2> {log}"


# add rule for mapping reads against the assembly for QC:

    
rule map_long_reads_to_assembly:
    input:
        assembly = rules.consolidate_genomes.output,
        reads = "long_reads/{sample}.fastq.gz"
    output:
        samfile = "09-assembly_QC/long_read_mapping/{sample}/{sample}.aln.sam"
    threads:
        workflow.cores
    log: 
        "logs/map_long_reads_to_assembly.{sample}.log"
    run:
        shell("bwa index {input.assembly} 2> {log}"),
        shell("minimap2 -ax map-ont -t {threads} {input.assembly} {input.reads} 2>> {log} 1> {output.samfile}"),


rule samtools_stats_long_reads:
    input: 
        samfile = "09-assembly_QC/long_read_mapping/{sample}/{sample}.aln.sam",
    output:
        bamfile = "09-assembly_QC/long_read_mapping/{sample}/{sample}.reads.sorted.bam",
        idxstats = "09-assembly_QC/long_read_mapping/{sample}/{sample}.idxstats.tab", 
        coverage = "09-assembly_QC/long_read_mapping/{sample}/{sample}.coverage.tab",
        
    threads:
        workflow.cores
    log: 
        "logs/samtools_stats_long_reads.{sample}.log"
    run:
        shell("samtools sort -@ {threads} -o {output.bamfile} -T reads.tmp {input.samfile} 2> {log}"),
        shell("samtools index -@ {threads} {output.bamfile} 2>> {log}"),
        shell("samtools idxstats {output.bamfile} 2>> {log} 1> {output.idxstats}"),
        shell("samtools coverage {output.bamfile} 2>> {log} 1> {output.coverage}")

# idxstats prints: contig contig_length #mapped #unmapped


rule samtools_get_unmapped_long_reads:
    input: 
        bamfile = "09-assembly_QC/long_read_mapping/{sample}/{sample}.reads.sorted.bam"
    output:
        bamfile = "09-assembly_QC/long_read_mapping/{sample}/{sample}.unmapped.bam",
        sorted_bamfile = "09-assembly_QC/long_read_mapping/{sample}/{sample}.unmapped.sorted.bam",
        fastq = "09-assembly_QC/long_read_mapping/{sample}/{sample}.unmapped.fastq.gz"
    threads:
        workflow.cores
    log: 
        "logs/samtools_get_unmapped_long_reads.{sample}.log"
    run:
        shell("samtools view -h -f 4 {input.bamfile} 2> {log} 1> {output.bamfile}"),
        shell("samtools sort -n -o {output.sorted_bamfile} {output.bamfile} 2>> {log}"),
        shell("samtools fastq {output.bamfile} | gzip 2>> {log} 1> {output.fastq}")


rule map_short_reads_to_assembly:
    input:
        assembly = rules.consolidate_genomes.output,
        SR1 = "short_reads/{sample}.R1.fastq.gz",
        SR2 = "short_reads/{sample}.R2.fastq.gz"
    output:
        samfile = "09-assembly_QC/short_read_mapping/{sample}/{sample}.aln.sam"
    threads:
        workflow.cores
    log: 
        "logs/map_short_reads_to_assembly.{sample}.log"
    run:
        shell("bwa index {input.assembly} 2> {log}"),
        shell("bwa mem -t {threads} -a {input.assembly} {input.SR1} 2>> {log} 1> {output.samfile}")

rule samtools_stats_short_reads:
    input: 
        samfile = rules.map_short_reads_to_assembly.output
    output:
        bamfile = "09-assembly_QC/short_read_mapping/{sample}/{sample}.reads.sorted.bam",
        idxstats = "09-assembly_QC/short_read_mapping/{sample}/{sample}.idxstats.tab", 
        coverage = "09-assembly_QC/short_read_mapping/{sample}/{sample}.coverage.tab",
        
    threads:
        workflow.cores
    log: 
        "logs/samtools_stats_short_reads.{sample}.log"
    run:
        shell("samtools sort -@ {threads} -o {output.bamfile} -T reads.tmp {input.samfile} 2> {log}"),
        shell("samtools index -@ {threads} {output.bamfile} 2>> {log}"),
        shell("samtools idxstats {output.bamfile} 2>> {log} 1> {output.idxstats}"),
        shell("samtools coverage {output.bamfile} 2>> {log} 1> {output.coverage}")

# idxstats prints: contig contig_length #mapped #unmapped


rule samtools_get_unmapped_short_reads:
    input: 
        bamfile = "09-assembly_QC/short_read_mapping/{sample}/{sample}.reads.sorted.bam"
    output:
        bamfile = "09-assembly_QC/short_read_mapping/{sample}/{sample}.unmapped.bam",
        sorted_bamfile = "09-assembly_QC/short_read_mapping/{sample}/{sample}.unmapped.sorted.bam",
        fastq = "09-assembly_QC/short_read_mapping/{sample}/{sample}.unmapped.fastq.gz"
    threads:
        workflow.cores
    log: 
        "logs/samtools_get_unmapped_short_reads.{sample}.log"
    run:
        shell("samtools view -h -f 4 {input.bamfile} 2> {log} 1> {output.bamfile}"),
        shell("samtools sort -n -o {output.sorted_bamfile} {output.bamfile} 2>> {log}"),
        shell("samtools fastq {output.bamfile} | gzip 2>> {log} 1> {output.fastq}")






















# I can't add the ideel code here, because it requires a diamond database that would need to be independently set up, which may be too much of an ask.  Also can only have one snakemake pipeline active at once.







# IDS, = glob_wildcards("input_reads/{sample}.fastq*")
