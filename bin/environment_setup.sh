conda create -n asm-fastqc -c bioconda fastqc -y

conda create -n asm-multiqc -c bioconda multiqc -y

conda create -n asm-filtlong -c bioconda filtlong -y

conda create -n asm-flye -c bioconda flye -y

conda create -n asm-medaka -c conda-forge -c bioconda medaka -y

# may need to add conda-forge before bioconda for medaka

conda create -n asm-bwa -c bioconda bwa -y

conda create -n asm-polypolish -c bioconda polypolish -y

conda create -n asm-circlator -c bioconda circlator -y

conda create -n asm-checkm -c bioconda checkm-genome -y

conda create -n asm-prokka -c bioconda prokka -y

conda create -n asm-minimap2 -c bioconda minimap2 -y
# split out minimap and bwa envs

conda create -n asm-samtools -c bioconda samtools -y

conda create -n snakemake -c conda-forge -c bioconda -c anaconda snakemake graphviz -y

