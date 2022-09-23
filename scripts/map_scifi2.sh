#!/bin/bash

sample_fname="dash_4"
sample_name="dash_4"
sample_fname="bw_scifi2_1"
sample_name="bw_scifi2"
##sample_fname="bw_scifi4_1"
##sample_name="bw_scifi4_1"
#
sample_file="/Genomics/grid/users/bw21/gitrepos/scBactoSeq/data/Barcodes/${sample_fname}.csv"
config="/Genomics/grid/users/bw21/gitrepos/scBactoSeq/metadata/default.yaml"
out_folder="/Genomics/adamsonlab/brucewang/data/${sample_name}_map_out_unique"
#out_folder="/Genomics/adamsonlab/brucewang/data/${sample_name}_map_out"
input_bam_glob="/Genomics/grid/users/bw21/Data/${sample_name}_out2/${sample_name}#{sample_name}.bam"
r2_barcode_dir="/Genomics/grid/users/bw21/gitrepos/scBactoSeq/metadata/"
echo $ilumina_run_folder
echo $input_bam_glob
#bwa_version="0.7.16a-r1181"
echo $out_folder
rm -r ${out_folder}
mkdir  ${out_folder}
echo $out_folder

# set the bash variable needed for the command-line
 scifi \
    map \
	-c ${config} \
	-o ${out_folder} \
        --input-bam-glob "${input_bam_glob}" \
        /Genomics/grid/users/bw21/gitrepos/scBactoSeq/data/Barcodes/sample_annotation.csv


sample_fname="bw_scifi2_2"
sample_name="bw_scifi2_2"
#sample_fname="bw_scifi4_2"
#sample_name="bw_scifi4_2"

sample_file="/Genomics/grid/users/bw21/gitrepos/scBactoSeq/data/Barcodes/${sample_fname}.csv"
config="/Genomics/grid/users/bw21/gitrepos/scBactoSeq/metadata/default.yaml"
out_folder="/Genomics/adamsonlab/brucewang/data/${sample_name}_map_out_unique"
#out_folder="/Genomics/adamsonlab/brucewang/data/${sample_name}_map_out"
input_bam_glob="/Genomics/grid/users/bw21/Data/${sample_name}_out2/${sample_name}#{sample_name}.bam"
r2_barcode_dir="/Genomics/grid/users/bw21/gitrepos/scBactoSeq/metadata/"
echo $ilumina_run_folder
echo $input_bam_glob
#bwa_version="0.7.16a-r1181"
echo $out_folder
rm -r ${out_folder}
mkdir  ${out_folder}
echo $out_folder

# set the bash variable needed for the command-line
 scifi \
    map \
	-c ${config} \
	-o ${out_folder} \
        --input-bam-glob "${input_bam_glob}" \
        /Genomics/grid/users/bw21/gitrepos/scBactoSeq/data/Barcodes/sample_annotation.csv


