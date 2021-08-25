#! /bin/bash
# this file is snakemake.sh
module load snakemake samtools cutadapt STAR || exit 1

sbcmd="sbatch --cpus-per-task={threads} --mem={cluster.mem}"
sbcmd+=" --time={cluster.time} --partition={cluster.partition}"
#sbcmd+=" --out={cluster.out} {cluster.extra}"

snakemake -pr --keep-going --local-cores $SLURM_CPUS_PER_TASK \
    --jobs 56 --cluster-config /data/CCBR/projects/ccbr1060/Hg38_shRNA_hybrid/resources/cluster.json --cluster "$sbcmd" \
    --latency-wait 120 all
