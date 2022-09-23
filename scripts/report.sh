#!/bin/bash

sample_name="bw_scifi2"
sample_name="dash_2"

sample_file="/Genomics/grid/users/bw21/gitrepos/scBactoSeq/data/Barcodes/${sample_name}.csv"
config="/Genomics/grid/users/bw21/gitrepos/scBactoSeq/metadata/default.yaml"
out_folder="/Genomics/grid/users/bw21/Data/${sample_name}_map_out"
#input_bam_glob="/Genomics/grid/users/bw21/Data/${sample_name}_out/bw_scifi#{sample_name}.bam"
input_bam_glob="/Genomics/grid/users/bw21/Data/${sample_name}_out/dash_2#{sample_name}.bam"
r2_barcode_dir="/Genomics/grid/users/bw21/gitrepos/scBactoSeq/metadata/"
#bwa_version="0.7.16a-r1181"
echo $out_folder

# set the bash variable needed for the command-line
 scifi \
    report \
	-c ${config} \
	-o ${out_folder} \
        /Genomics/grid/users/bw21/gitrepos/scBactoSeq/data/Barcodes/sample_annotation.csv


