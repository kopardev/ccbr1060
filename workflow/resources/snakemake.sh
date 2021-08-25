#! /bin/bash
# this file is snakemake.sh
module load snakemake samtools || exit 1

sbcmd="sbatch --cpus-per-task={threads} --mem={cluster.mem}"
sbcmd+=" --time={cluster.time} --partition={cluster.partition}"
sbcmd+=" --out={cluster.out} {cluster.extra}"

snakemake -pr --keep-going --local-cores $SLURM_CPUS_PER_TASK \
    --jobs 10 --cluster-config cluster.json --cluster "$sbcmd" \
    --latency-wait 120 all
