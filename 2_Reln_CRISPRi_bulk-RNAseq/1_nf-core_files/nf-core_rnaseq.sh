#!/bin/bash

# load Anaconda to activate nf-core tools and nextflow env
# load Singularity

module load Anaconda3/2021.11
module load Singularity

conda activate $USER_DATA/nf-core_nextflow_env/

# output dir
mkdir -p ./results

# execute the pipeline with cheaha config
# params values in param.yml and process args in config


nextflow run nf-core/rnaseq \
    -r 3.10 \
    --input samplesheet.csv \
    --outdir ./results \
    -profile cheaha \
    -params-file params.yml \
    -bg \ # this is similar to screen for running in background