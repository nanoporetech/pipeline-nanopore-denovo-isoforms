#!/bin/bash

snakemake clean_workdir --configfile "evaluation/eval_config.yml"
#snakemake -j 1 all --configfile "evaluation/eval_config.yml"
snakemake -j 10 merge_clustering --configfile "evaluation/eval_config.yml"
