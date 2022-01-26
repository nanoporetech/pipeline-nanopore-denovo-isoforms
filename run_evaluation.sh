#!/bin/bash

CONF="evaluation/eval_config.yml"
FASTQ="evaluation/data/SIRV_PCS109_phmm_fl.fq"
REF="evaluation/data/SIRV_150601a.fasta"
BAM="evaluation/SIRV_PCS109_phmm_fl.bam"
CLS="evaluation/pipeline-isONclust2_SIRV_E0/final_clusters/clusters.tsv"
SIZES="evaluation/pipeline-isONclust2_SIRV_E0/final_clusters/clusters_info.tsv"
TSV="evaluation/pipeline-isONclust2_SIRV_E0/clustering_quality.tsv"
PDF="evaluation/pipeline-isONclust2_SIRV_E0/clustering_quality.pdf"

#rm -fr evaluation/pipeline-isONclust2_SIRV_E0
snakemake -p -j 20 all --configfile "evaluation/eval_config.yml"
#snakemake --jobs 1000 --configfile $CONF --latency-wait 600 --drmaa-log-dir "sge_logs" --drmaa ' -P applications -V -cwd -l h_vmem=300G,mem_free=50G -pe mt 70' all

minimap2 -ax splice -t 2 $REF $FASTQ | samtools view -q 2 -F 2304 -b - | samtools sort - -o $BAM
samtools index $BAM

./scripts/compute_cluster_quality.py --sizes $SIZES --outfile $TSV --ont --clusters $CLS --classes $BAM --report $PDF
