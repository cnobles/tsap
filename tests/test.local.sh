#!/usr/bin/env bash
set -ev

# Test script
PREFIX=${HOME}/miniconda3
export PATH=${PATH}:${PREFIX}/bin
source activate vivi
CORES=${1-1}

# Create test analysis directory
snakemake analysis/test --configfile configs/test.config.yml -np
snakemake analysis/test --configfile configs/test.config.yml

# Move test sequence files to analysis directory
cp tests/Data/*.fastq.gz analysis/test/input_data/

# Generate test DAG graph
snakemake --configfile configs/test.config.yml -np
snakemake --configfile configs/test.config.yml --dag | dot -Tsvg > analysis/test/reports/test.dag.svg
snakemake --configfile configs/test.config.yml --latency-wait 30 --cores ${CORES}
cat analysis/test/processData/test-1.key.csv
cat analysis/test/output/unique_aligns.test.csv
