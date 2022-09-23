#!/bin/bash

#folder="bw_scifi3_2"
#sample_a="bw_scifi3"
folder="bw_scifi4_2"
sample_a="bw_scifi4"


in_folder="/Genomics/adamsonlab/brucewang/data/${folder}_out/${folder}.bam"
out_folder="/Genomics/adamsonlab/brucewang/data/${folder}_out_split/"
out_file=="/Genomics/adamsonlab/brucewang/data/${folder}_out_split/${folder}.bam"
metrics=="${folder}_metrics.tsv"
sample_annotation="/Genomics/grid/users/bw21/gitrepos/scBactoSeq/data/Barcodes/${sample_a}_samples.tsv"
sample_annotation_csv="/Genomics/grid/users/bw21/gitrepos/scBactoSeq/data/Barcodes/${sample_a}_samples.csv"
echo $ilumina_run_folder
mkdir $out_folder
#bwa_version="0.7.16a-r1181

csvformat -T $sample_annotation_csv > $sample_annotation

# set the bash variable needed for the command-line
picard_jar="/Genomics/grid/users/bw21/bin/picard-2.19.2-CeMM-all.jar" 
echo $out_folder
head -10 ${sample_annotation}
# One step to fix unpaired reads (just in case)
#java -Xmx1000m -jar $picard_path \
#    FixMateInformation \
#           I=$input_bam
# Convert the bam into an unmapped fastq for realignment

java \
     -Xmx20G \
     -Djava.io.tmpdir=./tmp \
     -jar ${picard_jar} \
     IlluminaSamDemux \
     INPUT= ${in_folder} \
     OUTPUT_DIR= ${out_folder} \
     OUTPUT_PREFIX= ${folder} \
     LIBRARY_PARAMS= ${sample_annotation} \
     METRICS_FILE= ${metrics} \
     TMP_DIR=./tmp \
     COMPRESSION_LEVEL=9 \
     CREATE_MD5_FILE=true \
     OUTPUT_FORMAT=bam \
     BARCODE_TAG_NAME=BC \
     BARCODE_QUALITY_TAG_NAME=QT \
     MAX_MISMATCHES=3 \
     MIN_MISMATCH_DELTA=3 \
     MAX_NO_CALLS=3 \
     MINIMUM_BASE_QUALITY=-1 \
     VERBOSITY=INFO \
     QUIET=false \
     VALIDATION_STRINGENCY=STRICT \
     MAX_RECORDS_IN_RAM=500000 \
     CREATE_INDEX=false \
     GA4GH_CLIENT_SECRETS=client_secrets.json \
     USE_JDK_DEFLATER=false \
     USE_JDK_INFLATER=false \
     DEFLATER_THREADS=4 \
     MATCHING_THREADS=4 \
     READ_STRUCTURE= 8M13B5S8B37M63T \
     TAG_PER_MOLECULAR_INDEX=RX \
     TAG_PER_MOLECULAR_INDEX=r2
#     READ_STRUCTURE= 8M13B5S8B14S16M20T54S \
#     READ_STRUCTURE= 8M13B5S8B14S16M74T \
#     READ_STRUCTURE= 8M13B5S8B16M20T54S \
