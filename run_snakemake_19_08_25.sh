#!/usr/bin/env bash
#SBATCH --account=project0039
#SBATCH --partition=nodes,smp                 # try smp; it's healthy. If you must use nodes, change both 'smp' to 'nodes'
#SBATCH --job-name=snk_ctl
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=5-00:00:00
#SBATCH --output=controller.out
#SBATCH --error=controller.err
set -euo pipefail

# 1) Pin the cluster's Slurm client/libs on THIS node
module purge
module load slurm/23.11                 # <-- match your site's minor version
module load apps/apptainer/1.3.4
module load apps/python3/3.12.7/gcc-8.5.0

# 2) Keep Slurm clean of stray libs
unset LD_LIBRARY_PATH
export PATH=/usr/bin:/bin:/usr/sbin:/sbin:$PATH
hash -r

# 3) Ensure 'conda' is available in batch shells
if ! command -v conda >/dev/null 2>&1; then
  source "$HOME/miniconda3/etc/profile.d/conda.sh"
fi

# 4) Preflight: prove sbatch works from this node (same acct/partition you'll use)
echo "=== SLURM preflight ==="
echo "sbatch: $(which sbatch)"
sbatch --version || true
ldd "$(which sbatch)" | grep -i slurm || true
scontrol ping || true
echo "Submitting 1s test…"
sbatch -A project0039 -p smp -n1 -t1 --wrap "hostname"
echo "=== preflight OK ==="

mkdir -p logs

# 5) Run Snakemake FROM your Conda env, without activating the env
conda run --no-capture-output -n snk \
  snakemake -s snakefile_8_8_25 \
    --executor slurm \
    --jobs 40 \
    --jobname "schisto.{rule}.{jobid}" \
    --slurm-logdir logs \
    --retries 3 \
    --latency-wait 60 \
    --rerun-incomplete \
    --keep-going \
    --printshellcmds \
    --default-resources \
        mem_mb=24000 \
        threads=4 \
        runtime=7200 \
        slurm_account=project0039 \
        slurm_partition=smp,nodes \
    --use-apptainer
