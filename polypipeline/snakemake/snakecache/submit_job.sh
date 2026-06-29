#!/bin/bash
# PIPELINE_ID est passé en env par le Perl
JOB_ID=$(sbatch \
    --cpus-per-task=20 \
    --partition=defq \
    --parsable \
    --job-name="pp_${PIPELINE_ID}" \
    "$@")
echo "$JOB_ID"