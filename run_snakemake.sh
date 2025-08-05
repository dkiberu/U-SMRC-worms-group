#!/usr/bin/env bash
#SBATCH --account=project0039          # run under your project
#SBATCH --partition=nodes              # only the “nodes” partition
#SBATCH --job-name=snk_ctl             # name of the controller job
#SBATCH --cpus-per-task=1              # the controller itself uses 1 core
#SBATCH --mem=2G                       # and very little RAM
#SBATCH --time=7-00:00:00              # plenty of wall-time
#SBATCH --output=controller.out
#SBATCH --error=controller.err

module load apps/apptainer/1.3.4
module load apps/python3/3.12.7/gcc-8.5.0
mkdir -p logs

snakemake \
  --executor slurm \
  --jobs 64 \
  --jobname "schisto.{rule}.{jobid}" \
  --slurm-logdir logs \
  --retries 2 \
  --latency-wait 60 \
  --rerun-incomplete \
  --keep-going \
  --default-resources \
      mem_mb=24000 \
      threads=4 \
      runtime=7200 \
      slurm_account=project0039 \
      slurm_partition=nodes \
  \
  --use-apptainer
