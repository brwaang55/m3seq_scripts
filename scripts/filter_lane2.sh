#!/bin/bash

sample_fname="bw_scifi7"
sample_name="bw_scifi7_2"
#sample_fname="dash_4"
#sample_name="dash_4"

sample_file="/Genomics/grid/users/bw21/gitrepos/scBactoSeq/data/Barcodes/${sample_fname}.csv"
config="/Genomics/grid/users/bw21/gitrepos/scBactoSeq/metadata/default.yaml"
out_folder="/Genomics/adamsonlab/brucewang/data/${sample_name}_map_out_unique3"
#out_folder="/Genomics/adamsonlab/brucewang/data/${sample_name}_map_out_unique3"
input_bam_glob="/Genomics/grid/users/bw21/Data/${sample_name}_map_out_unique3_unique3/bw_scifi#{sample_name}.bam"
#input_bam_glob="/Genomics/adamsonlab/brucewang/data/${sample_name}_map_out_unique3/bw_scifi7#{sample_name}.bam"
r2_barcode_dir="/Genomics/grid/users/bw21/gitrepos/scBactoSeq/metadata/"
#bwa_version="0.7.16a-r1181"
echo $out_folder

# set the bash variable needed for the command-line
# scifi7 \
scifi \
    filter \
	-c ${config} \
	-o ${out_folder} \
        --overwrite \
        /Genomics/grid/users/bw21/gitrepos/scBactoSeq/data/Barcodes/sample_annotation.csv


