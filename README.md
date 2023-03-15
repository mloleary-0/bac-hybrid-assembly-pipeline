# bac-hybrid-assembly-pipeline
This is a hybrid assembly snakemake pipeline using flye, medaka, and polypolish.  Intended for bacteria.  See dag.pdf for details.


\

In short, this pipeline will:

1) Run basic QC on your input reads, but will not process them (i.e., it will tell you if your reads are poor quality but will not do anything about it).
2) Run a flye 2.9 assembly
3) Polish with Medaka and Polypolish
4) Attempt to rotate molecules with Circlator fixstart (it will not run other circlator programs)
5) Do a final polish round with Polypolish (principally to make sure the contig ends were appropriately polished)
6) Run CheckM and Prokka on the assemblies
7) Pull out reads that were not mapped to the assemblies (I find this useful for troubleshooting when things don't go well).

Caution: be sure to check your assembly graph with something like Bandage.  The files fed to circlator will be rotated regardless of evidence of circularity.  If your contigs are not circular, the polypolish output should be fine for further analysis.




\
Inputs:  This pipeline assumes that you have reads formatted in the following way
Long reads: {strain}.fastq.gz 
Short reads: {strain}.R1.fastq.gz {strain}.R2.fastq.gz
These should be in the folders "long_reads" and "short_reads", respectively.  
Any number of sets of reads can be assembled at once, as long as they match this format.  In the future, I will update this pipeline so that the file extensions are modifiable in the config file.
