#!/bin/bash

CONF="evaluation/eval_config.yml"
FASTQ="evaluation/data/SIRV_PCS109_phmm_fl.fq"
REF="evaluation/data/SIRV_150601a.fasta"
BAM="evaluation/SIRV_PCS109_phmm_fl.bam"
TSV="evaluation/pipeline-isONclust2_SIRV_E0/clustering_quality.tsv"
PDF="evaluation/pipeline-isONclust2_SIRV_E0/clustering_quality.pdf"

rm -fr evaluation/pipeline-isONclust2_SIRV_E0
snakemake -j 20 all --configfile "evaluation/eval_config.yml"

minimap2 -ax splice -t 10 $REF $FASTQ | samtools view -q 5 -F 2304 -b - | samtools sort -@ 10 - -o $BAM
samtools index $BAM

./scripts/compute_cluster_quality.py --outfile $TSV --ont --clusters evaluation/pipeline-isONclust2_SIRV_E0/final_clusters/clusters.tsv --classes $BAM --report $PDF
