# snakemake --rulegraph | dot -Tpdf > dag.pdf
# snakemake --forceall --dag | dot -Tpdf > dag.pdf

conda create -n seqtk -c bioconda seqtk pigz -y
cd short_reads/SRR21131491.fastq/
seqtk seq SRR21131491.fastq -1 | pigz -c > OC8.R1.fastq.gz  
seqtk seq SRR21131491.fastq -2 | pigz -c > OC8.R2.fastq.gz
mv OC8* ..
cd ../..



  seqtk sample -s100 short_reads/OC8.R1.fastq.gz 50000 > short_reads/OC8-sub.R1.fastq
  seqtk sample -s100 short_reads/OC8.R1.fastq.gz  50000 > short_reads/OC8-sub.R2.fastq

seqtk sample long_reads/OC8.fastq.gz 20000  > long_reads/OC8-sub.fastq

mkdir test_read_source
mv long_reads/*.gz test_read_source
mv short_reads/*.gz test_read_source

gzip short_reads/*.fastq
gzip long_reads/*.fastq