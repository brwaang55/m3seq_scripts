#!/bin/bash

sample_name='bw_scifi2'
sample_file="/Genomics/grid/users/bw21/gitrepos/scBactoSeq/data/Barcodes/${sample_name}.csv"
config="/Genomics/grid/users/bw21/gitrepos/scBactoSeq/metadata/default.yaml"
out_folder="/Genomics/grid/users/bw21/Data/${sample_name}_map_out"
input_bam_glob="/Genomics/grid/users/bw21/Data/bw_scifi_out/${sample_name}#{sample_name}.bam"
r2_barcode_dir="/Genomics/grid/users/bw21/gitrepos/scBactoSeq/metadata/"
echo $ilumina_run_folder
#bwa_version="0.7.16a-r1181"
echo $out_folder
expression="${out_folder}/expression"
metrics="${out_folder}/metrics"
reads="${out_folder}/reads"
feature="${out_folder}/feature"

for i in $(find  /Genomics/grid/users/bw21/Data/${sample_name}_map_out/${sample_name} -name "*.featureCounts.bam") 
do 
	file_name=${i##*/}
	echo "$file_name"
        new_fname="${i}.sorted.bam"
        samtools sort $i > $new_fname
	output=$(samtools depth $new_fname -b  ~/Data/genomes/MG1655_BS168/MG1655_BS1682.bed)
	output_path="${i}.coverage.txt"
        echo "running depth for ${i}"
	echo -e "$output" >> $output_path
        rm $new_fname
#	mito_file_name="/Genomics/ayroleslab2/bruce/mitochondria/bamfiles/${file_name}.mitochondria.bam"
#	original_bam_path="/Genomics/ayroleslab2/lamaya/mitochondria/bamfiles/${file_name}"
done

