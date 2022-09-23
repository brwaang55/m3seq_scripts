#!/bin/bash
out_path="/Genomics/grid/users/bw21/runlogs/PA01_transcriptome"
genome_path="/Genomics/adamsonlab/brucewang/genomes/ncbi-genomes-2021-05-11/GCF_000022345.1_ASM2234v1_genomic.fna"
transcriptome_path="/Genomics/adamsonlab/brucewang/genomes/ncbi-genomes-2021-05-11/GCF_000022345.1_ASM2234v1_genomic.gff"
cmd="gffread -F -w $out_path -g $genome_path $transcriptome_path"
echo $cmd


