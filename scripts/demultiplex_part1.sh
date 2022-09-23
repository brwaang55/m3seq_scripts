#!/bin/bash

folder=$1
folder="bw_scifi2_2"
folder="dash_1"
folder="dash_4"
folder="bw_scifi3_1"
folder="bw_scifi4_1"
folder="bw_scifi7_1"


#illumina_run_folder="/Genomics/grid/users/bw21/Data/210317_M04292_0465_000000000-G7HRY"
#illumina_run_folder="/Genomics/grid/users/bw21/Data/210324_M04292_0469_000000000-G7M2W"
#illumina_run_folder="/Genomics/grid/users/bw21/Data/000000000-G7MHJ"
#illumina_run_folder="/Genomics/adamsonlab/brucewang/data/Sequencing/000000000-G7WL8"
illumina_run_folder="/Genomics/grid/users/bw21/Data/210616_A00741_0197_AH73FWDRXY"
illumina_run_folder="/Genomics/adamsonlab/brucewang/data/210630_A00741_0204_AH73GLDRXY"
#illumina_run_folder="/Genomics/adamsonlab/brucewang/data/210718_M04292_0003_000000000-DBVG7"
illumina_run_folder="/Genomics/adamsonlab/brucewang/data/210723_A00741_0215_BHHKF2DRXY"
out_folder="/Genomics/adamsonlab/brucewang/data/${folder}_out/${folder}.bam"
mkdir "/Genomics/adamsonlab/brucewang/data/${folder}_out/"
echo $ilumina_run_folder
#bwa_version="0.7.16a-r1181"

# set the bash variable needed for the command-line
picard_jar="/Genomics/grid/users/bw21/bin/picard-2.19.2-CeMM-all.jar" 
echo $out_folder
# One step to fix unpaired reads (just in case)
#java -Xmx1000m -jar $picard_path \
#    FixMateInformation \
#           I=$input_bam
# Convert the bam into an unmapped fastq for realignment
java \
    -Xmx20G \
    -Djava.util.concurrent.ForkJoinPool.common.parallelism=2 \
    -Djava.io.tmpdir=./tmp \
    -jar ${picard_jar} \
    IlluminaBasecallsToMultiplexSam \
    RUN_DIR= ${illumina_run_folder} \
    LANE=1 \
    OUTPUT= ${out_folder} \
    SEQUENCING_CENTER=BSF \
    NUM_PROCESSORS=2 \
    APPLY_EAMSS_FILTER=false \
    INCLUDE_NON_PF_READS=true \
    TMP_DIR=tmp \
    CREATE_MD5_FILE=false \
    FORCE_GC=false \
    MAX_READS_IN_RAM_PER_TILE=9000000 \
    MINIMUM_QUALITY=2 \
    VERBOSITY=INFO \
    QUIET=false \
    VALIDATION_STRINGENCY=STRICT \
    CREATE_INDEX=false \
    GA4GH_CLIENT_SECRETS=client_secrets.json
