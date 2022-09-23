
/Genomics/grid/users/bw21/miniconda3/envs/cutadaptenv/bin/python3.7 -u -m scifi.scripts.summarizer \
    --r1-annot /Genomics/grid/users/bw21/gitrepos/scBactoSeq/data/Barcodes/bw_scifi_plate_anno.csv \
    --r1-attributes plate_well,treatment \
    --cell-barcodes r1 r2 \
    --only-summary \
    --no-save-intermediate \
    --min-umi-output 0 \
    --expected-cell-number 100000 \
    --no-output-header \
    --save-gene-expression \
    --correct-r1-barcodes \
    --species-mixture --sample-name NEB_30MIN_FIX_PENTA_H12 \
    /Genomics/grid/users/bw21/Data/bw_scifi2_map_out2/bw_scifi2/NEB_30MIN_FIX_PENTA_H12/NEB_30MIN_FIX_PENTA_H12.ALL.STAR.Aligned.out.bam.featureCounts.bam \
    /Genomics/grid/users/bw21/Data/bw_scifi2_map_out2/bw_scifi2/NEB_30MIN_FIX_PENTA_H12/NEB_30MIN_FIX_PENTA_H12.ALL


