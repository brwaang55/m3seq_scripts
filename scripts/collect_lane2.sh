#!/bin/bash

sample_name="bw_scifi7_2"
sample_fname="bw_scifi7"
#sample_name="dash_4"

sample_file="/Genomics/grid/users/bw21/gitrepos/scBactoSeq/data/Barcodes/${sample_name}.csv"
config="/Genomics/grid/users/bw21/gitrepos/scBactoSeq/metadata/default.yaml"
#out_folder="/Genomics/grid/users/bw21/Data/${sample_name}_map_out_unique3"
out_folder="/Genomics/adamsonlab/brucewang/data/${sample_name}_map_out_unique3"
r2_barcode_dir="/Genomics/grid/users/bw21/gitrepos/scBactoSeq/metadata/"
echo $ilumina_run_folder
#bwa_version="0.7.16a-r1181"
echo $out_folder
expression="${out_folder}/expression"
metrics="${out_folder}/metrics"
reads="${out_folder}/reads"
feature="${out_folder}/feature"
mapping="${out_folder}/mapping"
coverage="${out_folder}/coverage"
bam="${out_folder}/bam"
summary="${out_folder}/summary"
#rm -r ${expression}
#rm -r ${metrics}
#rm -r ${reads}
#rm -r ${reads}
#rm -r ${feature}
#rm -r ${coverage}
#rm -r ${bam}
#rm -r ${summary}

mkdir ${expression}
mkdir ${metrics}
mkdir ${mapping}
#rm -r ${reads}
mkdir ${reads}
mkdir ${feature}
mkdir ${coverage}
mkdir ${bam}
mkdir ${summary}
for i in $( find  /Genomics/adamsonlab/brucewang/data/${sample_name}_map_out_unique3/${sample_fname}/ -name "*expression.csv.gz")
do 
	mv ${i} ${expression}
done 
for i in $( find  /Genomics/adamsonlab/brucewang/data/${sample_name}_map_out_unique3/${sample_fname}/ -name "*metrics.csv.gz")
do 
	mv ${i} ${metrics}
done 
for i in $( find  /Genomics/adamsonlab/brucewang/data/${sample_name}_map_out_unique3/${sample_fname}/ -name "*STAR.Log.final.out")
do 
	mv ${i} ${mapping}
done 
for i in $( find  /Genomics/adamsonlab/brucewang/data/${sample_name}_map_out_unique3/${sample_fname}/  -name "*reads.csv.gz")
do 
	mv ${i} ${reads}
done 
for i in $( find  /Genomics/adamsonlab/brucewang/data/${sample_name}_map_out_unique3/${sample_fname}/ -name "*.featureCounts.quant_gene.tsv")
do 
	mv ${i} ${feature}
done 
for i in $( find  /Genomics/adamsonlab/brucewang/data/${sample_name}_map_out_unique3/${sample_fname}/ -name "*.coverage.txt")
do 
	mv ${i} ${coverage}
done 
#for i in $( find  /Genomics/adamsonlab/brucewang/data/${sample_name}_map_out_unique3/${sample_fname}/ -name "*.bam.featureCounts.bam")
#do 
#	mv ${i} ${bam}
#done 
for i in $( find  /Genomics/adamsonlab/brucewang/data/${sample_name}_map_out_unique3/${sample_fname}/ -name "*summarized*")
do
        mv ${i} ${summary}
done



#find  /Genomics/grid/users/bw21/Data/bw_scifi_map_out_unique3/bw_scifi1/ -name "*expression.csv.gz" | xargs mv "${out_folder}/expression"
#find  /Genomics/grid/users/bw21/Data/bw_scifi_map_out_unique3/bw_scifi1/ -name "*expression.csv.gz" | xargs mv "${out_folder}/metrics"

