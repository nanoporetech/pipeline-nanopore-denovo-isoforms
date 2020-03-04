#!/bin/bash

FASTQ_URL="http://ftp.sra.ebi.ac.uk/vol1/fastq/ERR358/005/ERR3588905/ERR3588905_1.fastq.gz"
REF_URL="http://ftp.ensembl.org/pub/release-99/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.28.dna.toplevel.fa.gz"
CONF="evaluation/eval_config_dmel.yml"
OUT_DIR="evaluation/pipeline-isONclust2_dmel/results"
DATA_DIR="evaluation/pipeline-isONclust2_dmel/data"
FASTQ="$DATA_DIR/ERR3588905_1.fastq"
REF="$DATA_DIR/Drosophila_melanogaster.BDGP6.28.dna.toplevel.fa.gz"
BAM="${OUT_DIR}/dmel_pcs109.bam"
TSV="${OUT_DIR}/clustering_quality_dmel.tsv"
PDF=$OUT_DIR/clustering_report.pdf
CLS=$OUT_DIR/final_clusters/clusters.tsv
SIZES=$OUT_DIR/final_clusters/clusters_info.tsv
CORES=30

rm -fr $OUT_DIR/results
mkdir -p $OUT_DIR/data

if [ ! -f $REF ];
then (cd $DATA_DIR/data; wget $REF_URL)
fi

if [ ! -f $FASTQ ];
then (cd $DATA_DIR/data; wget $FASTQ_URL; gzip -d ${FASTQ}.gz)
fi

snakemake -j $CORES all --configfile "evaluation/eval_config_dmel.yml"

minimap2 -ax splice -t $CORES $REF $FASTQ | samtools view -q 5 -F 2304 -b - | samtools sort -@ 10 - -o $BAM
samtools index $BAM

./scripts/compute_cluster_quality.py --sizes $SIZES --outfile $TSV --ont --clusters $CLS --classes $BAM --report $PDF
