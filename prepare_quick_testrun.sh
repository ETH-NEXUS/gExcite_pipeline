#!/usr/bin/env bash

mv testdata/results_and_fastqs.tar.gz .
tar -xf results_and_fastqs.tar.gz
touch results/pooled_samples/citeseq_count/PBMC_D1.tags.tsv
touch results/pooled_samples/citeseq_count/PBMC_D1.run_report.yaml
